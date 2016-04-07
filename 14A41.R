library(caret)
library(ggplot2)
library(reshape2)
library(plyr)
#source("ch1scoring_fc.R")
#source("ch2scoring_fc.R")
require(synapseClient)

data.dir = "../data/"
output.dir = "../doc/"

#train.file = "Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv"
#train.file = "Drug_synergy_data/train_merged.csv" 
#test.file = "Drug_synergy_data/ch1_test_monoTherapy.csv"
test.file = "syn4925540"
train.ch1 = "syn4925542"
lb.ch1 = "syn5699550"
lb.ch2 = "syn5699551"

train.feature = "syn5757887"
test.feature.ch1 = "syn5757962"
test.feature.ch2 = "syn5785376"
#test.ch2 = "syn4925546"

cell.info = "syn4925560"
methylation = "syn4925573"

scoring.script.ch1 = "syn5785644"
scoring.script.ch2 = "syn5785645"

main<-function() {
    set.seed(142341) 
    #run_mode = "only_train"
    run_mode = "only_predict"
    #run_mode = "train_and_predict"
    #challenge = "ch1a"
    challenge = "ch1b_test" 
    #challenge = "ch2" 
    synapseLogin()
    mod.list = list(rfFit=NULL, gbmFit=NULL, modFit=NULL)
    if(run_mode == "only_train" | run_mode =="train_and_predict") {
	mod.list = train.model(challenge)
    } 
    if(run_mode == "train_and_predict" | run_mode =="only_predict") {
	rfFit = mod.list$rfFit
	gbmFit = mod.list$gbmFit
	modFit = mod.list$modFit
	get.predictions(challenge, rfFit, gbmFit, modFit, rebuild=T) 
    }
}


get.synergy.data<-function(challenge, is.train) {
    if(is.train) {
	# Use all available data (train ch1 / leaderboards ch1 & ch2)
	#file.name = train.file
	file.name = getFileLocation(synGet(train.ch1))
	f = read.csv(file.name)
	file.name = getFileLocation(synGet(lb.ch1))
	f = rbind(f, read.csv(file.name))
	file.name = getFileLocation(synGet(lb.ch2))
	f = rbind(f, read.csv(file.name))
    } else{
	#file.name = test.file
	file.name = getFileLocation(synGet(test.file))
	f = read.csv(file.name) 
    }
    #f = read.csv(paste0(data.dir, file.name))
    f$cat = f[,"SYNERGY_SCORE"]
    f = f[f$QA==1,]
    return(f);
}


process.features<-function(f, challenge, is.train=T) {
    # Add the features to the data frame
    if(is.train) {
	#file.name = paste0(data.dir, "features_train.dat")
	file.name = getFileLocation(synGet(train.feature))
	e = read.table(file.name, header=T)
    } else if(challenge == "ch2") {
	#file.name = paste0(data.dir, "features_test_ch2.dat")
	file.name = getFileLocation(synGet(test.feature.ch2))
	e = read.table(file.name, header=T)
    } else {
	#file.name  = paste0(data.dir, "features_test_ch1.dat")
	file.name = getFileLocation(synGet(test.feature.ch1))
	e = read.table(file.name, header=T)
    }
    if(challenge != "ch2") { # need to be in features file created by python
	#file.name = paste0(data.dir, "cell_info.csv")
	file.name = getFileLocation(synGet(cell.info))
	d = read.csv(file.name) 
	# Put tissue info as a categorized feature
	d = data.frame(row.names = d$Sanger.Name, tissue = gsub("[ /()]", ".", d$Tissue..General.))
	f$.t = d[f$CELL_LINE, "tissue"]
	f = cbind(f, model.matrix(~.t + 0, f))
	f = f[,-which(colnames(f) == ".t")]
	#ggplot(f, aes(x=CELL_LINE, y=SYNERGY_SCORE)) + geom_boxplot(aes(color=tissue))
	#print(summary(f))
    }

    e$gexp.diff = abs(e$gexpA.amed - e$gexpB.amed) 
    e$gexp = e$gexpA.amed * e$gexpB.amed 
    e$met.diff = abs(e$metA.amed - e$metB.amed)
    e$met = e$metA.amed * e$metB.amed 
    e$mut = ifelse(is.na(abs(e$mutA-e$mutB)), min(e$mutA, e$mutB, na.rm=T), abs(e$mutA-e$mutB))
    e$cnv = abs(e$cnvA - e$cnvB)
    e$k.diff = abs(e$kA-e$kB)
    e$k.min = apply(e[,c("kA", "kB")], 1, min)
    e$k.max = apply(e[,c("kA", "kB")], 1, max)
    e$kegg.in = (e$kegg.inA + e$kegg.inB) 
    e$cosmic.in = (e$cosmic.inA + e$cosmic.inB) 

    # Assign cell line / drug median synergy as a feature 
    #if(is.train) { # it does not exists for test data - need to store it in feature file
    #	a = ddply(f, ~ CELL_LINE, summarize, syn.med = median(cat))
    #	b = ddply(f, ~ COMBINATION_ID, summarize, syn.med = median(cat))
    #	#d = ddply(f, ~tissue, summarise, syn.med = median(SYNERGY_SCORE))
    #}

    if(challenge == "ch2" & is.train == F) {
	f = e
	f$cat = NA
    } else {
	m = nrow(f)
	n = ncol(f) + ncol(e)
	f.mod = data.frame(matrix(nrow=m, ncol=n))
	for(i in 1:m) {
	    comb = as.character(f[i, "COMBINATION_ID"])
	    cell = as.character(f[i, "CELL_LINE"])
	    f.mod[i,] = c(f[i,], e[e$comb.id == comb & e$cell.line == cell, ])
	}
	colnames(f.mod) = c(colnames(f), colnames(e))
	# To correct integer converted comb names and cell lines
	f.mod$comb.id = factor(f$COMBINATION_ID)
	f.mod$cell.line = factor(f$CELL_LINE)
	f = f.mod
	# Drug response based features
	f$einf = abs(f$Einf_A - f$Einf_B)
	f$conc = abs(f$MAX_CONC_A - f$MAX_CONC_B)
	f$ic50 = abs(f$IC50_A - f$IC50_B)
	f$h = abs(f$H_A - f$H_B)
	#print(head(f)) 
    }

    # Keep the descriptive information in f.mod 
    f.mod = f[,c("comb.id", "cell.line", "cat")]
    #print(head(f.mod))
    # Remove insensitive cell lines (Einf)
    if(is.train) {
	d = f
	d$einf = (d$Einf_A + d$Einf_B) / 2
	a = ddply(d, ~ cell.line, summarize, syn.sd = sd(cat), syn.min = min(cat), syn.med = median(cat), einf.sd = sd(einf), einf.min = min(einf), einf.med = median(einf))
	# Cells with high min Einf has lower synergy
	#b = cor.test(a$syn.med, a$einf.min, use="complete")
	#print(sprintf("Cor b/w einf.min and syn.med: %f %f", b$estimate[[1]], b$p.value)) # -0.235
	cutoff = 40 # 20 (8 cell lines) # 10 (22 cell lines) 
	# 22RV1 CAL-120 HCC1143 HCC1428 HCC1937 KU-19-19 UACC-812 VCaP
	cell.lines = as.vector(a[a$einf.min>cutoff,"cell.line"])
	indices = f$cell.line %in% cell.lines
	if(length(indices) > 0) {
	    f = f[-indices,] 
	    f.mod = f.mod[-indices,]
	}
	# Remove combinations with missing info / low synergy 
	#indices = which(abs(f$cat) < 5) # - no improvement
	indices = which(is.na(f$cnv))
	if(length(indices) > 0) {
	    f = f[-indices,] 
	    f.mod = f.mod[-indices,] 
	}
    }
    print(sprintf("Number of instances: %d", nrow(f)))

    # Choose features to include
    #!
    features = c()
    indices = which(grepl("^\\.[cp]", colnames(f)))  # mz (m/g/e/z/c/p likely to worsen on ch1 train if only gbm is used)
    #!features = c(colnames(f)[indices], features)
    features = c("cnv", "sim.target", "sim.chemical", features) 
    features = c("cosmic.in", "kegg.in", features)
    features = c("guild.common", "guild.med", "guild.max", features)
    features = c("cell.med", "comb.med", features) 
    features = c(".pcensus", ".pkegg", features) 
    # Not included due to poor predictor performance
    #features = c("mut", features) 
    #features = c("dAB", features) 
    #features = c("k.diff", "k.min", "k.max", "dAB", features) 
    # Somewhat reliable sub
    #features = c("gexp", "met", "cnv", "sim.target", "sim.chemical", "cosmic.in", "kegg.in") 
    #features = c("guild.common", "guild.med", "guild.max", features)
    #features = c("cell.med", "comb.med", features) 
    #features = c(".pcensus", ".pkegg", features) 
    # Not used anymore
    #features = c("kegg.gex.med", "kegg.gex.max", "kegg.cnv.med", "kegg.cnv.max", "cosmic.gexp.med", "cosmic.gexp.max", "cosmic.cnv.med", "cosmic.cnv.max") 

    # Use all features
    if(challenge == "ch1a") {
	indices = which(grepl("^\\.[t]", colnames(f)))  
	features = c(colnames(f)[indices], features)
	indices = which(grepl("^\\.[ge]", colnames(f)))
	features = c(colnames(f)[indices], features) 
	features = c("einf", "h", "conc", "ic50", features) # likely toworsen
	features = c("gexp", "met", features) 
	features = c("gexp.diff", "met.diff", features) 
    # Use anything except monotherapy data
    } else if(challenge == "ch2") {
	indices = which(grepl("^\\.[ge]", colnames(f)))
	features = c(colnames(f)[indices], features)
	#f[,indices] = abs(f[,indices]) # use abs for e (and for g)
	features = c("gexp", "met", features)
	features = c("gexp.diff", "met.diff", features) 
    # Use anything except monotherapy, gexp and met data
    } else if(grepl("^ch1b", challenge)) {
	# Nothing special
	indices = which(grepl("^\\.[t]", colnames(f)))  
	features = c(colnames(f)[indices], features)
    } else {    
	stop("Unrecognized challenge!")
    }
    indices = which(colnames(f) %in% c(features)) 
    f = f[, indices]
    print(sprintf("Features: %s", paste(features, collapse=", ")))

    if(is.train == T) {
	# Check variance 
	nzv = nearZeroVar(f, saveMetrics= TRUE)
	print(nzv) 
	#print(rownames(nzv[nzv$zeroVar,])) 
	zv.idx = which(colnames(f) %in% rownames(nzv[nzv$zeroVar,]))
	print(sprintf("Zero variance variables: %s", paste(colnames(f)[zv.idx], collapse=", "))) 
	if (length(zv.idx) > 0) {
	    f = f[, -zv.idx]
	}
    }

    # Imputing and scaling
    f = predict(preProcess(f, method = c("center", "scale", "knnImpute"), k=5), f) # "BoxCox" 

    if(is.train == T) {
	# Check correlated features
	cor.mat = cor(f)
	cor.idx = findCorrelation(cor.mat, cutoff = 0.75) 
	print(sprintf("Correlated: %s", paste(colnames(f)[cor.idx], collapse=", "))) 
	if(length(cor.idx) > 0) {
	    f = f[,-cor.idx]
	}
    }

    f$cat = f.mod$cat
    print(summary(f))

    # Categorize synergy as 0/1 for challenge 2
    if(challenge == "ch2" & is.train == T) {
	f$cat = factor(ifelse(f$cat > 20, 1, 0))
	# Balance the data (0: 80%, 1: 20%) 
	indices.positive = which(f$cat == 1)
	indices = which(f$cat == 0)
	indices = sample(indices, length(indices.positive))
	indices = c(indices.positive, indices)
	f = f[indices,]
	f.mod = f.mod[indices,]
    }
    # Models have their built-in feature selection
    # Feature selection for classification tasks
    #model <- train(cat~., data=f, method="lvq", preProcess=NULL, trControl=trainControl(method = "cv"))
    #importance <- varImp(model, scale=FALSE)
    # Summarize importance
    #print(importance)
    return(list(f=f, f.mod=f.mod));
}


train.model<-function(challenge) {
    # Get synergy training data
    f = get.synergy.data(challenge, is.train=T)
    # Create expression and mutation based features
    val = process.features(f, challenge, is.train=T)
    f = val$f
    f.mod = val$f.mod
    f.mod$COMBINATION_ID = f.mod$comb.id
    f.mod$CELL_LINE = f.mod$cell.line

    inTrain = createDataPartition(y = f$cat, p = 0.7, list=F) 
    training = f[inTrain, ] 
    testing = f[-inTrain, ]

    # Build model(s)
    if(T) { 
	rfFit = NULL
	gbmFit = NULL
	modFit = NULL
	if(challenge != "ch2") {
	    ml.methods = c("rf", "gbm") 
	    ml.methods = c("gaussprRadial", ml.methods)
	} else {
	    # for ch2 
	    ml.methods = c("rf", "gbm") 
	    ml.methods = c("gaussprRadial", "glmnet", ml.methods)
	    ml.methods = c("LogitBoost", "ada", ml.methods) 
	}
	ctrl = trainControl(method = "cv", number = 10) 
	#ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
	for(ml.method in ml.methods) {
	    print(ml.method)
	    modFit = train(cat ~ ., data = training, method = ml.method, trControl = ctrl, verbose=F) 
	    pred = predict(modFit, testing)
	    print(predictors(modFit))
	    #print(modFit)
	    if(challenge == "ch2") { 
		a = confusionMatrix(pred, testing$cat)
	    } else {
		a = cor(pred, testing$cat) 
	    }
	    print(ml.method)
	    print(a)

	    # For scoring
	    # Output observed values for test set (for ch2 the original synergy file is used)
	    if(challenge != "ch2") {
		output.predictions(challenge, f.mod[-inTrain,], testing, suffix=".test.observed")
	    }
	    # Output predictions for testing set
	    testing$cat = pred
	    output.predictions(challenge, f.mod[-inTrain,], testing, suffix=".test")
	    # Get organizer's scores
	    check.scoring(challenge) 
	    
	    if(ml.method == "rf") {
		rfFit = modFit
	    } else if(ml.method == "gbm") {
		gbmFit = modFit
	    }
	}
	# gam to combine to predictors e.g., rf gbm
	pred.1 = predict(rfFit, testing)
	pred.2 = predict(gbmFit, testing)
	pred.comb = data.frame(pred.1, pred.2, cat=testing$cat)
	modFit = train(cat ~ ., data=pred.comb, method = "gam")
	pred = predict(modFit, testing)
    } else {
	# trainControl: boot for bootstrapping, cv for x-validation # repeatedcv for repeated 10-fold CV 
	#ctrl = trainControl(method = "cv", number = 10) 
	ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
	prep = NULL # c("center", "scale") # Already preprocessed above
	# Random forest
	rfFit = train(cat ~ ., data=training, method = "rf", preProcess = prep, trControl = ctrl) # using default for kNN, k=5
	# List the chosen features
	print(predictors(rfFit))
	pred = predict(rfFit, testing)
	# Without imputation it is possible that some are not predicted due to NA
	#d = extractPrediction(rfFit, testing)
	#print(head(d))
	if(challenge == "ch2") { 
	    a = confusionMatrix(pred, testing$cat)
	} else {
	    a = cor(pred, testing$cat) 
	}
	print(a)
	# Tree boost
	gbmFit = train(cat ~ ., data=training, method = "gbm", preProcess = prep, trControl = ctrl, verbose=F)
	print(predictors(gbmFit))
	pred = predict(gbmFit, testing)
	if(challenge == "ch2") { 
	    a = confusionMatrix(pred, testing$cat)
	} else {
	    a = cor(pred, testing$cat) 
	}
	print(a)
	# gam to combine to predictors e.g., rf gbm
	pred.1 = predict(rfFit, testing)
	pred.2 = predict(gbmFit, testing)
	pred.comb = data.frame(pred.1, pred.2, cat=testing$cat)
	modFit = train(cat ~ ., data=pred.comb, method = "gam")
	pred = predict(modFit, testing)
    }

    if(challenge == "ch2") { 
	a = confusionMatrix(pred, testing$cat)
    } else {
	a = cor(pred, testing$cat) 
    }
    print(a)

    dir.create(paste0(output.dir, challenge))
    write.table(testing$cat, paste0(output.dir, challenge, "/", "observed.dat"))
    write.table(pred, paste0(output.dir, challenge, "/", "prediction.dat"))

    # Plot prediction consistency
    p = qplot(pred, cat, data=testing, main = paste0("Prediction PCC: ", format(a, digits=2)))
    png(paste0(output.dir, challenge, "/", "training.png")) 
    print(p)
    dev.off()

    # Output observed values for test set (for ch2 the original synergy file is used)
    if(challenge != "ch2") {
	output.predictions(challenge, f.mod[-inTrain,], testing, suffix=".test.observed")
    }
    # Output predictions for testing set
    testing$cat = pred
    output.predictions(challenge, f.mod[-inTrain,], testing, suffix=".test")
    # Get organizer's scores
    check.scoring(challenge) 

    return(list(rfFit=rfFit, gbmFit=gbmFit, modFit=modFit));
}


get.predictions<-function(challenge, rfFit, gbmFit, modFit, rebuild=F) {
    # Get test data
    if(challenge == "ch2") {
	#file.name = paste0(data.dir, "features_test_ch2.dat")
	file.name = getFileLocation(synGet(test.feature.ch2))
	f = read.table(file.name, header=T)
    } else { 
	f = get.synergy.data(challenge, is.train=F)
    }
    # Create expression and mutation based features
    val = process.features(f, challenge, is.train=F)
    testing = val$f
    testing.mod = val$f.mod

    # Build models using all the training data
    if(rebuild) { 
	# Get training data
	f.training = get.synergy.data(challenge, is.train=T)
	val = process.features(f.training, challenge, is.train=T)
	training = val$f
	training.mod = val$f.mod
	#ctrl = trainControl(method = "cv")
	ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
	prep = NULL #c("center", "scale")
	if(challenge == "ch2") {
	    rfFit = train(cat ~ ., data=training, method = "rf", preProcess = prep, trControl = ctrl)
	    #rfFit = NULL
	    gbmFit = NULL
	    modFit = rfFit
	} else {
	    # Random forest
	    rfFit = train(cat ~ ., data=training, method = "rf", preProcess = prep, trControl = ctrl)
	    pred.1 = predict(rfFit, training)  
	    # Tree boost
	    gbmFit = train(cat ~ ., data=training, method = "gbm", preProcess = prep, trControl = ctrl, verbose=F)
	    pred.2 = predict(gbmFit, training)  
	    pred.comb = data.frame(pred.1, pred.2, cat=training$cat)
	    modFit = train(cat ~ ., data=pred.comb, method = "gam")
	}
    }

    # Make predictions using model(s)
    if(is.null(rfFit) | is.null(gbmFit)) {
	pred = predict(modFit, testing)
    } else {
	pred.1 = predict(rfFit, testing)  
	pred.2 = predict(gbmFit, testing) 
	pred.comb = data.frame(pred.1, pred.2)
	pred = predict(modFit, pred.comb)
    }

    # Get predictions for leaderboard data
    testing$cat = pred
    #print(head(testing$cat)) 

    # Output predictions
    output.predictions(challenge, f, testing)
}


output.predictions<-function(challenge, f, testing, suffix="") {
    # Get confidence scores for learderboard data (assign lower confidence to larger values)
    # Consider assigning scores based on the cell line senstivity (i.e. Einf) or performance on training set
    if(challenge == "ch2") { 
	# Output is categorical, if not converted 1s and 2s 
	testing$cat = ifelse(testing$cat == "1", 1, 0) 
	# Not very meaningful, prediction 0/1 # check whether there is 0 div error
	testing$conf = 1-as.numeric(testing$cat)/max(as.numeric(testing$cat))
    } else {
	testing$conf = 1-abs(testing$cat)/max(abs(testing$cat))
    }

    dir.create(paste0(output.dir, challenge))
    # Output predictions
    if(challenge == "ch2") {
	f$cat = testing$cat 
	f$conf = testing$conf
	d = acast(f, comb.id~cell.line, value.var="cat", fun.aggregate = function(x){round(mean(x, na.rm=T))}) #{max(x, na.rm=T)} # Max returns -Inf if all NA
	#print(summary(d))
	file.name = paste0(output.dir, challenge, "/", "synergy_matrix.csv", suffix)
	write.csv(d, file.name)
	d = acast(f, comb.id~cell.line, value.var="conf", fun.aggregate = function(x){round(mean(x, na.rm=T))}) 
	file.name = paste0(output.dir, challenge, "/", "confidence_matrix.csv", suffix)
	write.csv(d, file.name)
    } else {
	f$PREDICTION = testing$cat
	f$CONFIDENCE = testing$conf
	file.name = paste0(output.dir, challenge, "/", "prediction.csv", suffix)
	write.csv(f[,c("CELL_LINE", "COMBINATION_ID", "PREDICTION")], file.name, row.names=F)
	file.name = paste0(output.dir, challenge, "/", "combination_priority.csv", suffix)
	a = f[,c("COMBINATION_ID", "CONFIDENCE")]
	a = aggregate(CONFIDENCE~COMBINATION_ID, a, mean)
	write.csv(a, file.name, row.names=F)
    }
}


check.scoring<-function(challenge) {

    if(challenge == "ch2") {
	file.name = getFileLocation(synGet(scoring.script.ch2))
	source(file.name)
	#obs.file = paste0(output.dir, challenge, "/", "synergy_matrix.csv.test.observed")
	obs.file = paste0(data.dir, train.file)
	pred.file = paste0(output.dir, challenge, "/", "synergy_matrix.csv.test")
	conf.file = paste0(output.dir, challenge, "/", "confidence_matrix.csv.test")
	# conf.file="none"
	a = getGlobalScore_ch2(obs.file, pred.file)
	print(sprintf("Global score: %.3f", a))
	a = getPrecision_ch2(obs.file, pred.file, threshold=20)
	print(sprintf("F1 score: %.3f", a["F1"])) #print(a) prec sens  npv spec  auc  phi  BAC   F1 aupr
	a = getOneDimScore_ch2(obs.file, pred.file, rows=T) #confidence=conf.file
	print(sprintf("Row score (row): %.3f", a["mean"]))
	a = getOneDimScore_ch2(obs.file, pred.file, rows=F) #confidence=conf.file
	print(sprintf("Column score (row): %.3f", a["mean"]))
    } else {
	file.name = getFileLocation(synGet(scoring.script.ch1))
	source(file.name)
	obs.file = paste0(output.dir, challenge, "/", "prediction.csv.test.observed")
	pred.file = paste0(output.dir, challenge, "/", "prediction.csv.test")
	conf.file = paste0(output.dir, challenge, "/", "combination_priority.csv.test")
	# conf.file="none"
	a = getGlobalScore_ch1(obs.file, pred.file)
	print(sprintf("Primary score: %.3f (%.3f)", a["final"], a["tiebreak"]))
	print(sprintf("Partial correlation: %.3f", a["score"]))
	a = getDrugCombiScore_ch1(obs.file, pred.file, confidence=conf.file)
	print(sprintf("Secondary score: %.3f", a["mean"]))
    }
}


get.methylation.values.for.genes<-function() {
    # Gene id - symbol map # For the reverse map: org.Hs.egSYMBOL2EG 
    library(org.Hs.eg.db)
    x = org.Hs.egSYMBOL 
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    if(length(xx) > 0) {
	xx[gene.ids]
    }

    # From https://www.biostars.org/p/135446/
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
    genes = genes(txdb)
    #cpg = data.frame(vals=c("chr1:91190489-91192804","chr22:20674706-20675153"))
    #gr = makeGRangesFromDataFrame(cpg) # need strand info
    cpg = data.frame(chrname=c("chr1", "chr22"), start=c(91190489, 20674706), end=c(91192804, 20675153))
    gr = GRanges(cpg$chrname, IRanges(cpg$start, cpg$end))
    gene.ids = names(genes[nearest(gr, genes)])
    #gene.ids = genes[precede(gr, genes)]

    #file.name = paste0(data.dir, "methyl_ilse_m.csv")
    file.name = getFileLocation(synGet(methylation))
    d = read.csv(file.name)
    a = strsplit(as.character(d$X), ":")
    e.chr = unlist(lapply(a, function(x){x[1]}))
    table(e.chr)
    b = strsplit(unlist(lapply(a, function(x){x[2]})), "-")
    e.start = unlist(lapply(b, function(x){x[1]}))
    e.end = unlist(lapply(b, function(x){x[2]}))
    cpg = data.frame(chrname=e.chr, start=as.numeric(e.start), end=as.numeric(e.end))
    gr = GRanges(cpg$chrname, IRanges(cpg$start, cpg$end))
    gene.ids = names(genes[nearest(gr, genes)])
    gene.symbols = as.character(xx[gene.ids])
    e = d[,-1]
    selection.function<-function(asample){ 
       return(tapply(asample, factor(gene.symbols), function(x) { x[which.max(abs(x))] }))
    }
    d.merged = apply(e, 2, selection.function)
    write.csv(e, paste0(data.dir, "methylation.csv"))
}


main()

