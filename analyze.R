library(caret)
library(ggplot2)
library(reshape2)

data.dir = "../data/"
output.dir = "../doc/"

main<-function() {
    set.seed(142341)
    #set.seed(555444)
    #set.seed(999999)
    #challenge = "ch1a"
    challenge = "ch1b"
    #challenge = "ch2"
    rfFit = NULL
    gbmFit = NULL
    modFit = NULL
    mod.list = train.model(challenge)
    return() #!
    rfFit = mod.list$rfFit
    gbmFit = mod.list$gbmFit
    modFit = mod.list$modFit
    get.predictions(challenge, rfFit, gbmFit, modFit, rebuild=F) #rebuild=T
}


get.synergy.data<-function(file.name, challenge, is.train) {
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

    e$gexp = ifelse(is.na(abs(e$gexpA-e$gexpB)), min(abs(e$gexpA), abs(e$gexpB), na.rm=T), abs(e$gexpA-e$gexpB))
    e$mut = ifelse(is.na(abs(e$mutA-e$mutB)), min(e$mutA, e$mutB, na.rm=T), abs(e$mutA-e$mutB))
    #e[is.na(e$mut), "mut"] = 0 
    e$cnv = ifelse(is.na(abs(e$cnvA-e$cnvB)), min(e$cnvA, e$cnvB, na.rm=T), abs(e$cnvA-e$cnvB))
    #e[is.na(e$cnv),"cnv"] = 0
    #e$kegg.cnv = ifelse(is.na(abs(e$kegg.cnvA-e$kegg.cnvB)), min(e$kegg.cnvA, e$kegg.cnvB, na.rm=T), abs(e$kegg.cnvA-e$kegg.cnvB))
    #e$cosmic.cnv = ifelse(is.na(abs(e$cosmic.cnvA-e$cosmic.cnvB)), min(e$cosmic.cnvA, e$cosmic.cnvB, na.rm=T), abs(e$cosmic.cnvA-e$cosmic.cnvB))

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
	f = f.mod
	#print(head(f))
    }

    # Choose features to include
    #features = c("gexp", "mut", "cnv", "guild.med", "guild.max", "sim.target", "sim.chemical", "kegg.gex.med", "kegg.gex.max", "kegg.cnv.med", "kegg.cnv.max", "cosmic.gexp.med", "cosmic.gexp.max", "cosmic.cnv.med", "cosmic.cnv.max", "kegg.in", "cosmic.in") # (no in) 41.9 / (in) 39.9
    #features = c("gexp", "mut", "cnv", "guild.med", "guild.max", "sim.target", "sim.chemical", "kegg.in", "cosmic.in") # (no in) 41.3 / (in) 44.9
    #features = c("gexp", "mut", "cnv", "guild.med", "guild.max", "sim.target", "sim.chemical", "kegg.in", "cosmic.in") 
    #features = c("gexp", "mut") 
    #features = c(features, colnames(e)[3:(which(colnames(e) == "gexpA")-1)]) 
    #!
    features = c()
    indices = which(grepl("^\\.[gmczp]", colnames(e))) #!
    features = colnames(e)[indices] #!
    features = c("guild.common", "guild.med", "guild.max", features)
    features = c("gexp", "mut", "cnv", "sim.target", "sim.chemical", "kegg.in", "cosmic.in", features)
    #features = c("mut", "cnv", "kegg.mut", "kegg.cnv") 67.3 (61 without kegg, 51.8 w/o cnv, 12 w/o mut, cosmic also lowers)
    #features = c("mut", "cnv", "gexp", features) # 69.9 (20), 71.2 (30) 76.8 (40) w/ sensitivity filtering below, 74.4 w/o removing kegg.mut (correlated)
    #features = c("gexp", "cnv") # 28.5 w/ refined feature # ch2 selection
    indices = which(colnames(f) %in% c(features)) # , "cat"

    print(sprintf("Features: %s", paste(features, collapse=", ")))

    # Remove insensitive cell lines (Einf)
    if(is.train) {
	library(plyr)
	d = f
	d$einf = (d$Einf_A + d$Einf_B) / 2
	a = ddply(d, ~ CELL_LINE, summarize, syn.sd = sd(cat), syn.min = min(cat), syn.med = median(cat), einf.sd = sd(einf), einf.min = min(einf), einf.med = median(einf))
	cutoff = 40 # 20 (8 cell lines) # 10 (22 cell lines)
	# 22RV1 CAL-120 HCC1143 HCC1428 HCC1937 KU-19-19 UACC-812 VCaP
	cell.lines = as.vector(a[a$einf.min>cutoff,"CELL_LINE"])
	f = f[!f$CELL_LINE %in% cell.lines,] 
	# Cells with high min Einf has lower synergy
	b = cor.test(a$syn.med, a$einf.min, use="complete")
	print(sprintf("Cor b/w einf.min and syn.med: %f %f", b$estimate[[1]], b$p.value)) # -0.235
    }

    # Keep a copy to assign cat values afterwards
    f.org = f
    # Use all features
    if(challenge == "ch1a") {
	f = f[, c(4:11, indices)]
    # Use anything except monotherapy data
    } else if(challenge == "ch1b" | challenge == "ch2") {
	f = f[, indices]
    } else {    
	stop("Unrecognized challenge!")
    }

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
	#print(summary(f))
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
	cor.idx = findCorrelation(cor.mat, cutoff = .75)
	print(sprintf("Correlated: %s", paste(colnames(f)[cor.idx], collapse=", "))) 
	if(length(cor.idx) > 0) {
	    f = f[,-cor.idx]
	}
    }

    f$cat = f.org$cat
    print(summary(f))

    # Categorize synergy as 0/1 for challenge 2
    if(challenge == "ch2" & is.train == T) {
	f$cat = factor(ifelse(f$cat > 30, 1, 0))
	# Balance the data (0: 80%, 1: 20%) #!
	indices.positive = which(f$cat == 1)
	indices = which(f$cat == 0)
	indices = sample(indices, length(indices.positive))
	indices = c(indices.positive, indices)
	f = f[indices,]
    }
    # Models have their built-in feature selection
    # Feature selection for classification tasks
    #model <- train(cat~., data=f, method="lvq", preProcess=NULL, trControl=trainControl(method = "cv"))
    #importance <- varImp(model, scale=FALSE)
    # Summarize importance
    #print(importance)
    return(f);
}


train.model<-function(challenge) {
    # Get synergy training data
    f = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", challenge, is.train=T)
    # Create expression and mutation based features
    f = process.features(f, challenge, is.train=T)

    inTrain = createDataPartition(y = f$cat, p = 0.7, list=F) 
    training = f[inTrain, ] 
    testing = f[-inTrain, ]

    # Build model(s)
    # trainControl: boot for bootstrapping, cv for x-validation # repeatedcv for repeated 10-fold CV 
    ctrl = trainControl(method = "cv") 
    #ctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
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
    # Ridge (regression L2 regularization)
    #ridgeFit = train(cat ~ ., data=training, method = "ridge", preProcess = prep, trControl = ctrl)
    # Linear regression
    #glmFit = train(cat ~ ., data=training, method = "glm", preProcess = prep, trControl = ctrl)
    # Lasso (regression with L1 regularization)
    #lassoFit = train(cat ~ ., data=training, method = "lasso", preProcess = prep, trControl = ctrl)
    # gam to combine to predictors e.g., rf glm
    pred.1 = predict(rfFit, testing)
    pred.2 = predict(gbmFit, testing)
    pred.comb = data.frame(pred.1, pred.2, cat=testing$cat)
    modFit = train(cat ~ ., data=pred.comb, method = "gam")
    pred = predict(modFit, testing)
    if(challenge == "ch2") { 
	a = confusionMatrix(pred, testing$cat)
    } else {
	a = cor(pred, testing$cat) 
    }
    print(a)

    write.table(testing$cat, paste0(output.dir, "/", "observed.dat"))
    write.table(pred, paste0(output.dir, challenge, "/", "prediction.dat"))

    # Plot prediction consistency
    p = qplot(pred, cat, data=testing, main = paste0("Prediction PCC: ", format(a, digits=2)))
    png(paste0(output.dir, challenge, "/", "training.png")) 
    print(p)
    dev.off()

    return(list(rfFit=rfFit, gbmFit=gbmFit, modFit=modFit));
}


get.predictions<-function(challenge, rfFit, gbmFit, modFit, rebuild=F) {
    # Get test data
    if(challenge == "ch2") {
	# will not be taken into account (since does not contain all cell lines)
	#f = get.synergy.data("Drug_synergy_data/ch2_leaderBoard_monoTherapy.csv", challenge, is.train=F)
	f = read.table(paste0(data.dir, "features_test_ch2.dat"), header=T)
    } else { 
	f = get.synergy.data("Drug_synergy_data/ch1_leaderBoard_monoTherapy.csv", challenge, is.train=F)
    }
    # Create expression and mutation based features
    testing = process.features(f, challenge, is.train=F)

    # Build models using all the training data
    if(rebuild) { 
	# Get training data
	f.training = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", challenge, is.train=T)
	f.training = process.features(f.training, challenge, is.train=T)
	training = f.training
	ctrl = trainControl(method = "cv")
	prep = NULL #c("center", "scale")
	# Random forest
	rfFit = train(cat ~ ., data=training, method = "rf", preProcess = prep, trControl = ctrl)
	pred.1 = predict(rfFit, training)  
	# Tree boost
	gbmFit = train(cat ~ ., data=training, method = "gbm", preProcess = prep, trControl = ctrl, verbose=F)
	pred.2 = predict(gbmFit, training)  
	pred.comb = data.frame(pred.1, pred.2, cat=training$cat)
	modFit = train(cat ~ ., data=pred.comb, method = "gam")
    }

    # Make predictions using model(s)
    pred.1 = predict(rfFit, testing)  
    #pred.1 = predict(rfFit, testing, type = "prob") # prob not meaningful for RF
    pred.2 = predict(gbmFit, testing) 
    pred.comb = data.frame(pred.1, pred.2)
    pred = predict(modFit, pred.comb)

    # Get predictions for leaderboard data
    testing$cat = pred

    # Get confidence scores for learderboard data (assign lower confidence to larger values)
    #! Consider assigning scores based on the cell line senstivity (i.e. Einf) or performance on training set
    if(challenge == "ch2") { 
	#! Not very meaningful, prediction 0/1 # check whether there is 0 div error
	testing$conf = 1-as.numeric(testing$cat)/max(as.numeric(testing$cat))
    } else {
	testing$conf = 1-abs(testing$cat)/max(abs(testing$cat))
    }

    # Output predictions
    if(challenge == "ch2") {
	#! check why output is 1s and 2s
	f$cat = as.numeric(testing$cat) #ifelse(testing$cat > 10, 1, 0)
	f$conf = testing$conf
	library(reshape2)
	d = acast(f, comb.id~cell.line, value.var="cat")
	file.name = paste0(output.dir, challenge, "/", "synergy_matrix.csv")
	write.csv(d, file.name)
	d = acast(f, comb.id~cell.line, value.var="conf")
	file.name = paste0(output.dir, challenge, "/", "confidence_matrix.csv")
	write.csv(d, file.name)
    } else {
	f$PREDICTION = testing$cat
	f$CONFIDENCE = testing$conf
	file.name = paste0(output.dir, challenge, "/", "prediction.csv")
	write.csv(f[,c("CELL_LINE", "COMBINATION_ID", "PREDICTION")], file.name, row.names=F)
	file.name = paste0(output.dir, challenge, "/", "combination_priority.csv")
	a = f[,c("COMBINATION_ID", "CONFIDENCE")]
	a = aggregate(CONFIDENCE~COMBINATION_ID, a, mean)
	write.csv(a, file.name, row.names=F)
    }
}


main()

