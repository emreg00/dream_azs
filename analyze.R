library(caret)
library(ggplot2)
library(reshape2)
library(plyr)

data.dir = "../data/"
output.dir = "../doc/"

main<-function() {
    set.seed(142341)
    #set.seed(555444)
    #set.seed(999999)
    #challenge = "ch1a"
    #challenge = "ch1b" # not allowed to  use response / gexp data 
    challenge = "ch2"
    rfFit = NULL
    gbmFit = NULL
    modFit = NULL
    mod.list = train.model(challenge)
    #return() #!
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

    #e$gexp = ifelse(is.na(abs(e$gexpA-e$gexpB)), min(abs(e$gexpA), abs(e$gexpB), na.rm=T), abs(e$gexpA-e$gexpB))
    e$gexp = abs(e$gexpA + e$gexpB)
    #e$gexp.diff = abs(e$gexpA-e$gexpB)
    e$mut = ifelse(is.na(abs(e$mutA-e$mutB)), min(e$mutA, e$mutB, na.rm=T), abs(e$mutA-e$mutB))
    #e[is.na(e$mut), "mut"] = 0 
    #e$cnv = ifelse(is.na(abs(e$cnvA-e$cnvB)), min(e$cnvA, e$cnvB, na.rm=T), abs(e$cnvA-e$cnvB))
    e$cnv = abs(e$cnvA-e$cnvB)
    e$k.diff = abs(e$kA-e$kB)
    e$k.min = apply(e[,c("kA", "kB")], 1, min)
    e$k.max = apply(e[,c("kA", "kB")], 1, max)
    #e[is.na(e$cnv),"cnv"] = 0
    #e$kegg.cnv = ifelse(is.na(abs(e$kegg.cnvA-e$kegg.cnvB)), min(e$kegg.cnvA, e$kegg.cnvB, na.rm=T), abs(e$kegg.cnvA-e$kegg.cnvB))
    #e$cosmic.cnv = ifelse(is.na(abs(e$cosmic.cnvA-e$cosmic.cnvB)), min(e$cosmic.cnvA, e$cosmic.cnvB, na.rm=T), abs(e$cosmic.cnvA-e$cosmic.cnvB))

    # Assign cell line avg synergy as a feature - no improvement
    #if(is.train) {
    #	a = ddply(f, ~ CELL_LINE, summarize, syn.med = median(cat))
    #	b = ddply(f, ~ COMBINATION_ID, summarize, syn.med = median(cat))
    #}

    if(challenge == "ch2" & is.train == F) {
	f = e
	f$cat = NA
    } else {
	m = nrow(f)
	n = ncol(f) + ncol(e)
	#n = n + 2 # for checking avg synergy over combinations / cell lines
	f.mod = data.frame(matrix(nrow=m, ncol=n))
	for(i in 1:m) {
	    comb = as.character(f[i, "COMBINATION_ID"])
	    cell = as.character(f[i, "CELL_LINE"])
	    f.mod[i,] = c(f[i,], e[e$comb.id == comb & e$cell.line == cell, ])
	    # Assign cell line avg synergy as a feature
	    #f.mod[i,] = c(f.mod[i,], a[a$CELL_LINE == cell, "syn.med"], b[b$COMBINATION_ID == comb, "syn.med"]) 
	}
	colnames(f.mod) = c(colnames(f), colnames(e)) #, "cell.med", "comb.med")
	# To correct integer converted comb names and cell lines
	f.mod$comb.id = factor(f$COMBINATION_ID)
	f.mod$cell.line = factor(f$CELL_LINE)
	f = f.mod
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
	f = f[-indices,] 
	f.mod = f.mod[-indices,]
	# Remove combinations with missing info / low synergy 
	#indices = which(abs(f$cat) < 5) # - no improvement
	indices = which(is.na(f$cnv))
	#indices = which(is.na(f$cnv) | is.na(f$guild.med))
	#indices = which(is.na(f$cnv) | is.na(f$sim.chemical))
	f = f[-indices,] 
	f.mod = f.mod[-indices,] 
    }
    print(sprintf("Number of instances: %d", nrow(f)))

    # Choose features to include
    #features = c("gexp", "mut", "cnv", "guild.med", "guild.max", "sim.target", "sim.chemical", "kegg.gex.med", "kegg.gex.max", "kegg.cnv.med", "kegg.cnv.max", "cosmic.gexp.med", "cosmic.gexp.max", "cosmic.cnv.med", "cosmic.cnv.max", "kegg.in", "cosmic.in") # (no in) 41.9 / (in) 39.9
    #features = c("gexp", "mut", "cnv", "guild.med", "guild.max", "sim.target", "sim.chemical", "kegg.in", "cosmic.in") # (no in) 41.3 / (in) 44.9
    #features = c("gexp", "mut", "cnv", "guild.med", "guild.max", "sim.target", "sim.chemical", "kegg.in", "cosmic.in") 
    #features = c("gexp", "mut") 
    #features = c(features, colnames(e)[3:(which(colnames(e) == "gexpA")-1)]) 
    #!
    features = c()
    indices = which(grepl("^\\.[gcz]", colnames(e))) # gcz
    features = colnames(e)[indices]
    #features = c("cnv", "sim.chemical", "cosmic.in", features)
    features = c("gexp", "cnv", "sim.target", "sim.chemical", "cosmic.in", "kegg.in", features)
    #features = c("cnv", "sim.target", "sim.chemical", "cosmic.in", "kegg.in", features)
    features = c("guild.common", "guild.med", "guild.max", features)
    features = c("k.diff", "k.min", "k.max", "dAB", features)
    #features = c(".cEGFR", ".cBRAF", ".cTOP2A", ".cTOP2B", ".cALK", ".cATR", ".cTUBB1", ".cBCL2L1", ".cIGF1R", ".cMAP2K3", ".cCHECK1", ".cADAM17", ".cPIP5KIA", ".cFGFR1", features)
    #features = c("gexp", "mut", "cnv", "sim.target", "sim.chemical", "kegg.in", "cosmic.in", features)
    #features = c("cell.med",  "comb.med", features)
    indices = which(colnames(f) %in% c(features)) # , "cat"

    print(sprintf("Features: %s", paste(features, collapse=", ")))


    # Before keeping a copy to assign cat values afterwards
    #f.mod = f
    # Use all features
    if(challenge == "ch1a") {
	f = f[, c(4:11, indices)]
    # Use anything except monotherapy data
    } else if(grepl("^ch1b", challenge) | challenge == "ch2") {
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
    f = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", challenge, is.train=T)
    # Create expression and mutation based features
    val = process.features(f, challenge, is.train=T)
    f = val$f
    f.mod = val$f.mod

    inTrain = createDataPartition(y = f$cat, p = 0.7, list=F) 
    training = f[inTrain, ] 
    testing = f[-inTrain, ]

    # Build model(s)
    if(T) { #!
	# Simple decision tree for developmental/debugging purposes
	#modFit = train(cat ~ ., data = training, method='ctree', tuneLength=10,
	#	  trControl=trainControl(method='cv', number=10)) #, classProbs=F, summaryFunction=twoClassSummary))
	modFit = train(cat ~ ., data = training, method='glmnet', trControl=trainControl(method='cv', number=10)) 
	#modFit = myClassifier(training)
	rfFit = NULL
	gbmFit = NULL
	pred = predict(modFit, testing)
	#print(predictors(modFit))
	#print(modFit)
    } else {
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
    }

    if(challenge == "ch2") { 
	a = confusionMatrix(pred, testing$cat)
    } else {
	a = cor(pred, testing$cat) 
    }
    print(a)

    return(list(rfFit=rfFit, gbmFit=gbmFit, modFit=modFit)); #!

    dir.create(paste0(output.dir, challenge))
    write.table(testing$cat, paste0(output.dir, "/", "observed.dat"))
    write.table(pred, paste0(output.dir, challenge, "/", "prediction.dat"))

    # Plot prediction consistency
    p = qplot(pred, cat, data=testing, main = paste0("Prediction PCC: ", format(a, digits=2)))
    png(paste0(output.dir, challenge, "/", "training.png")) 
    print(p)
    dev.off()

    f.mod$COMBINATION_ID = f.mod$comb.id
    f.mod$CELL_LINE = f.mod$cell.line
    output.predictions(challenge, f.mod[-inTrain,], testing)

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
    val = process.features(f, challenge, is.train=F)
    testing = val$f
    testing.mod = val$f.mod

    # Build models using all the training data
    if(rebuild) { 
	# Get training data
	f.training = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", challenge, is.train=T)
	val = process.features(f.training, challenge, is.train=T)
	training = val$f
	training.mod = val$f.mod
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
    if(is.null(rfFit)) {
	pred = predict(modFit, testing)
    } else {
	pred.1 = predict(rfFit, testing)  
	#pred.1 = predict(rfFit, testing, type = "prob") # prob not meaningful for RF
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


output.predictions<-function(challenge, f, testing) {
    # Get confidence scores for learderboard data (assign lower confidence to larger values)
    #! Consider assigning scores based on the cell line senstivity (i.e. Einf) or performance on training set
    if(challenge == "ch2") { 
	# Not very meaningful, prediction 0/1 # check whether there is 0 div error
	testing$conf = 1-as.numeric(testing$cat)/max(as.numeric(testing$cat))
    } else {
	testing$conf = 1-abs(testing$cat)/max(abs(testing$cat))
    }

    dir.create(paste0(output.dir, challenge))
    # Output predictions
    if(challenge == "ch2") {
	f$cat = testing$cat 
	print(head(f$cat)) #!
	f$conf = testing$conf
	d = acast(f, comb.id~cell.line, value.var="cat")
	print(head(d)) #!
	#! check why output is 1s and 2s
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


main()

