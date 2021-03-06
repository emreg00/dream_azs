library(caret)
library(ggplot2)
library(reshape2)
library(plyr)
source("ch1scoring_fc.R")
source("ch2scoring_fc.R")

data.dir = "../data/"
output.dir = "../doc/"
#train.file = "Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv"
train.file = "Drug_synergy_data/train_merged.csv" 
test.file = "Drug_synergy_data/ch1_test_monoTherapy.csv"

main<-function() {
    set.seed(142341) 
    #set.seed(555444)
    #set.seed(999999)
    #!
    #run_mode = "only_train"
    run_mode = "only_predict"
    #run_mode = "train_and_predict"
    #challenge = "ch1a"
    #challenge = "ch1b" 
    challenge = "ch2" 
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


get.biomarkers<-function() {
    # Get predictions
    file.name = paste0(output.dir, "ch1b/prediction.csv")
    d = read.csv(file.name)
    # Get features
    file.name = paste0(data.dir, "features_test_ch1.dat")
    f = read.table(file.name, header=T)
    # Get combinations of interest
    file.name = paste0(data.dir, "biomarkers.csv")
    e = read.csv(file.name)
    comb.ids = e$comb.id
    # Put tissue info as a categorized feature
    file.name = paste0(data.dir, "cell_info.csv")
    e = read.csv(file.name) 
    e = data.frame(row.names = e$Sanger.Name, tissue = gsub("[ /()]", ".", e$Tissue..General.))
    f$.t = e[f$cell.line, "tissue"]
    f = cbind(f, model.matrix(~.t + 0, f))
    f = f[,-which(colnames(f) == ".t")]
    # Preprocess features
    f$cnv = abs(f$cnvA - f$cnvB)
    f$kegg.in = (f$kegg.inA + f$kegg.inB) 
    f$cosmic.in = (f$cosmic.inA + f$cosmic.inB) 
    features = c("comb.id","cell.line","cnv","cosmic.in","kegg.in","sim.chemical","sim.target","guild.common","guild.med","guild.max") 
    features = c(features,"cell.med","comb.med")
    features = c(features,".tbreast",".tlung",".turinary.tract")
    features = c(features,".tblood.lymph",".tgastrointestinal.tract..lower.",".tgastrointestinal.tract..upper.",".tmale.genital.system",".tsoft.tissue")
    f = f[,features]
    #d[d$COMBINATION_ID=="AKT.BCL2L1",]
    #f[f$comb.id == "AKT.BCL2L1" & f$cell.line=="MDA-MB-453",]
    e = d[d$COMBINATION_ID %in% comb.ids,]
    m = nrow(e)
    n = ncol(e) + ncol(f)
    e.mod = data.frame(matrix(nrow=m, ncol=n))
    for(i in 1:m) {
	comb = as.character(e[i, "COMBINATION_ID"])
	cell = as.character(e[i, "CELL_LINE"])
	e.mod[i,] = c(e[i,], f[f$comb.id == comb & f$cell.line == cell, ])
    }
    colnames(e.mod) = c(colnames(e), colnames(f))
    # To correct integer converted comb names and cell lines
    e.mod$comb.id = factor(e$COMBINATION_ID)
    e.mod$cell.line = factor(e$CELL_LINE)
    e = e.mod[,-c(1,2)]
    fit = aov(PREDICTION ~ cnv + guild.common + guild.med + sim.chemical + kegg.in + cosmic.in, data=e)
    summary(fit)
    e.sub = e[abs(e$PREDICTION)>=20,]
    features.sub = c("cnv","sim.chemical","cosmic.in","kegg.in","guild.med","guild.common")
    for(feature in features.sub) {
	a = cor.test(e.sub$PREDICTION, e.sub[[feature]], method="spearman")
	print(c(feature, a$estimate))
    }
}

get.synergy.data<-function(challenge, is.train) {
    if(is.train) {
	file.name = train.file
    } else{
	file.name = test.file
    }
    f = read.csv(paste0(data.dir, file.name))
    #f$syn = f[,"SYNERGY_SCORE"]
    #f$syn.einf = f[,"SYNERGY_SCORE"] / (1 + f[, "Einf_A"] + f[, "Einf_B"])
    #f$syn.ic = f[,"SYNERGY_SCORE"] / (1 + f[, "IC50_A"] + f[, "IC50_B"])
    #f$cut = cut(f[,"syn"], breaks=c(-10000, -100, -50, -10, 0, 10, 50, 100, 10000), labels=c("high-negative", "negative", "low-negative", "low-positive", "positive", "high-positive"))
    f$cat = f[,"SYNERGY_SCORE"]
    f = f[f$QA==1,]
    return(f);
}


process.features<-function(f, challenge, is.train=T) {
    # Add the features to the data frame
    if(is.train) {
	e = read.table(paste0(data.dir, "features_train.dat"), header=T)
    } else if(challenge == "ch2") {
	e = read.table(paste0(data.dir, "features_test_ch2.dat"), header=T)
    } else {
	e = read.table(paste0(data.dir, "features_test_ch1.dat"), header=T)
    }
    if(challenge != "ch2") { # need to be in features file created by python
	d = read.csv(paste0(data.dir, "cell_info.csv")) 
	# Put tissue info as a categorized feature
	d = data.frame(row.names = d$Sanger.Name, tissue = gsub("[ /()]", ".", d$Tissue..General.))
	f$.t = d[f$CELL_LINE, "tissue"]
	f = cbind(f, model.matrix(~.t + 0, f))
	f = f[,-which(colnames(f) == ".t")]
	#ggplot(f, aes(x=CELL_LINE, y=SYNERGY_SCORE)) + geom_boxplot(aes(color=tissue))
	#print(summary(f))
    }

    #e$gexp = ifelse(is.na(abs(e$gexpA-e$gexpB)), min(abs(e$gexpA), abs(e$gexpB), na.rm=T), abs(e$gexpA-e$gexpB))
    e$gexp.diff = abs(e$gexpA.amed - e$gexpB.amed) 
    e$gexp = e$gexpA.amed * e$gexpB.amed 
    e$met.diff = abs(e$metA.amed - e$metB.amed)
    e$met = e$metA.amed * e$metB.amed 
    e$mut = ifelse(is.na(abs(e$mutA-e$mutB)), min(e$mutA, e$mutB, na.rm=T), abs(e$mutA-e$mutB))
    #e[is.na(e$mut), "mut"] = 0 
    #e$cnv = ifelse(is.na(abs(e$cnvA-e$cnvB)), min(e$cnvA, e$cnvB, na.rm=T), abs(e$cnvA-e$cnvB))
    e$cnv = abs(e$cnvA - e$cnvB)
    e$k.diff = abs(e$kA-e$kB)
    e$k.min = apply(e[,c("kA", "kB")], 1, min)
    e$k.max = apply(e[,c("kA", "kB")], 1, max)
    e$kegg.in = (e$kegg.inA + e$kegg.inB) 
    e$cosmic.in = (e$cosmic.inA + e$cosmic.inB) 
    #e$kegg.cnv = ifelse(is.na(abs(e$kegg.cnvA-e$kegg.cnvB)), min(e$kegg.cnvA, e$kegg.cnvB, na.rm=T), abs(e$kegg.cnvA-e$kegg.cnvB))
    #e$cosmic.cnv = ifelse(is.na(abs(e$cosmic.cnvA-e$cosmic.cnvB)), min(e$cosmic.cnvA, e$cosmic.cnvB, na.rm=T), abs(e$cosmic.cnvA-e$cosmic.cnvB))

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
    #features = c("kegg.gex.med", "kegg.gex.max", "kegg.cnv.med", "kegg.cnv.max", "cosmic.gexp.med", "cosmic.gexp.max", "cosmic.cnv.med", "cosmic.cnv.max") # not used anymore
    #features = c("gexp", "cnv", "sim.chemical", "cosmic.in", features) # somewhat reliable sub
    #!
    features = c()
    indices = which(grepl("^\\.[cp]", colnames(f)))  # mz (m/g/e/z/c/p likely to worsen on ch1 train if only gbm is used)
    features = c(colnames(f)[indices], features)
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
    #f = predict(preProcess(f, method = c("center", "scale")), f) # no imputation
    #if(challenge == "ch1a" & is.train==T) {
	# Impute using category (synergy) information
	#f = predict(preProcess(f, method = c("center", "scale", "knnImpute"), k=5), f) # scales synergy as well
	#idx = which(colnames(f) == "cat")
	#d = f[,-idx]
	#d = f
	#d = predict(preProcess(d, method = c("center", "scale")), d) 
	#d$cat = f$cat
	#d = predict(preProcess(d, method = c("knnImpute"), k=5), d) 
	#f = d
    #} else {
	#idx = which(colnames(f) == "cat")
	#d = f[,-idx]
	#d = predict(preProcess(d, method = c("center", "scale", "knnImpute"), k=5), d) # "BoxCox" 
	#d$cat = f$cat #as.vector(scale(f$cat))
	#f = d
    #}
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
	# Simple decision tree for developmental/debugging purposes
	#modFit = train(cat ~ ., data = training, method='ctree', tuneLength=10,
	#	  trControl=trainControl(method='cv', number=10)) #, classProbs=F, summaryFunction=twoClassSummary))
	#modFit = myClassifier(training)
	rfFit = NULL
	gbmFit = NULL
	modFit = NULL
	#ml.methods = c("gbm")
	#ml.methods = c("ctree")
	#ml.methods = c("glmnet")
	ml.methods = c("rf", "gbm") 
	#ml.methods = c("gaussprRadial", "glmnet", ml.methods)
	#ml.methods = c("ada")
	#ml.methods = c("rf", "gbm") # 43 / 38
	#ml.methods = c("gaussprRadial") # 31
	#ml.methods = c("rlm", "bayesglm", "gaussprLinear") # poor
	#ml.methods = c("avNNet", "enet") # poor
	#ml.methods = c("lasso", "relaxo", "ridge") # poor
	# for ch2 
	#ml.methods = c("rf", "gbm") # 62.6 / 59 (sub global: 0.27 / 0.53)
	#ml.methods = c("gaussprRadial") # 54 (SG: 0.76)
	#ml.methods = c("ada") # 59.5 (SG: 0.46)
	#ml.methods = c("glm", "glmnet") # 48 / 48.5
	#ml.methods = c("LogitBoost", "glmnet") # 53.7 (SG: 0.02) / 48.6 (NA)
	#ml.methods = c("logreg") failed
	#ml.methods = c("gamboost") # failed
	#ctrl = trainControl(method = "cv", number = 10) 
	ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
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
	f = read.table(paste0(data.dir, "features_test_ch2.dat"), header=T)
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


# From http://stackoverflow.com/questions/29113106/write-custom-classifier-in-r-and-predict-function
# create a function that returns an object of class myClassifierClass
myClassifier = function(trainingData, ...) {
    model = structure(list(x = (trainingData[, "cnv"]/2 + trainingData[, "cosmic.in"] /2 + trainingData[, "guild.med"] + trainingData[, "guild.max"]) > 2, y = trainingData[, "cat"]), class = "myClassifierClass") 
    return(model)
}

# create a method for function print for class myClassifierClass
predict.myClassifierClass = function(modelObject) {
    return(rlogis(length(modelObject$y)))
} 


check.scoring<-function(challenge) {
    if(challenge == "ch2") {
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

    d = read.csv(paste0(data.dir, "methyl_ilse_m.csv"))
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

