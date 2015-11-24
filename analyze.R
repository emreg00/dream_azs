library(caret)
library(ggplot2)
library(reshape2)

data.dir = "../data/"
output.dir = "../doc/"

main<-function() {
    set.seed(142341)
    #set.seed(999999)
    #explore()
    challenge = "ch1a"
    #challenge = "ch1b"
    mod.list = train.model(challenge)
    rfFit = mod.list$rfFit
    gbmFit = mod.list$gbmFit
    modFit = mod.list$modFit
    get.predictions(challenge, rfFit, gbmFit, modFit)
}


get.synergy.data<-function(file.name, challenge) {
    f = read.csv(paste0(data.dir, file.name))
    f$syn = f[,"SYNERGY_SCORE"]
    f$syn.einf = f[,"SYNERGY_SCORE"] / (1 + f[, "Einf_A"] + f[, "Einf_B"])
    f$syn.ic = f[,"SYNERGY_SCORE"] / (1 + f[, "IC50_A"] + f[, "IC50_B"])
    #f$cut = cut(f[,"syn"], breaks=c(-10000, -100, -50, -10, 0, 10, 50, 100, 10000), labels=c("high-negative", "negative", "low-negative", "low-positive", "positive", "high-positive"))
    f = f[f$QA==1,]
    f$cat = f$syn 
    #! Not sure why syn.einf is not working
    #if(challenge == "ch1a") {
    #    f$cat = f$syn.einf # Challange 1A
    #} else if(challenge == "ch1b") {
    #    f$cat = f$syn # Challange 1B
    #}
    return(f);
}


get.expression.info<-function(file.name) {
    d = read.table(file.name, sep=",", row.names=1, quote='"', header=T, check.names=F)
    m = apply(d, 1, mean)
    s = apply(d, 1, sd)
    d.norm = abs((d-m) / s)
    return(d.norm);
}


get.mutation.info<-function(file.name) {
    e = read.csv(file.name)
    # Categorize mutations as 0: silent/unknown 1: non-synonimous 2: cancer
    e$mut = ifelse(e$Mutation.Description == "Substitution - coding silent" | e$Mutation.Description == "Unknown", 0, ifelse(e$FATHMM.prediction=="CANCER", 2, 1))
    # Remove .1 and _ENST suffices in Gene name for mutations
    e$gene = sapply(strsplit(as.character(e$Gene.name), "[._]"), function(x) { x[1] })
    return(e);
}


get.drug.info<-function(file.name) {
    g = read.csv(file.name)
    rownames(g) = g$ChallengeName
    return(g);
}


# add expression of genes (possibly with PCA or use lasso or pathway based filtering)
# add mutation info of genes (possibly with PCA or use lasso or pathway based filtering)
# add target - mutated gene info (i.e. distance in interactome / GUILD etc...)
create.features<-function(f, is.train=T) {
    # Get expression info
    d = get.expression.info(paste0(data.dir, "gex.csv"))
    # Get mutation info 
    e = get.mutation.info(paste0(data.dir, "mutations.csv"))
    # Get drug target info
    g = get.drug.info(paste0(data.dir, "Drug_info_release.csv"))
    # Get GUILD scores
    if(is.train) {
	h = read.table(paste0(data.dir, "guild_noexp.dat"), header=T) 
    } else {
	h = read.table(paste0(data.dir, "guild_test_noexp.dat"), header=T)
    }
    x = data.frame()

    for(i in 1:nrow(f)) {
	cell.line = as.character(f[i, "CELL_LINE"])
	compound.A = as.character(f[i, "COMPOUND_A"])
	compound.B = as.character(f[i, "COMPOUND_B"])
	comb.id = as.character(f[i, "COMBINATION_ID"])
	targets.A = as.character(g[compound.A, "Target.Official.Symbol."])
	targets.B = as.character(g[compound.B, "Target.Official.Symbol."])
	a = unlist(strsplit(targets.A, ","))
	b = unlist(strsplit(targets.B, ","))
	#print(c(cell.line, a))
	#print(d[a,cell.line])
	f[i, "gexp_A"] = mean(d[a,cell.line], rm.na=T)
	f[i, "gexp_B"] = mean(d[b,cell.line], rm.na=T)
	f[i, "mut_A"] = mean(e[e$gene %in% targets.A & e$cell_line_name == cell.line,"mut"], rm.na=T)
	f[i, "mut_B"] = mean(e[e$gene %in% targets.B & e$cell_line_name == cell.line,"mut"], rm.na=T)
	a = h[h$comb.id == comb.id & h$cell.line == cell.line, c("med", "mean", "sd", "max", "min")] # "max.a", "max.b", 
	if(nrow(a) == 0) {
	    a[1,] = 0
	}
	x = rbind(x, a)
    }   
    f[is.na(f$gexp_A),"gexp_A"] = 0
    f[is.na(f$gexp_B),"gexp_B"] = 0
    f[is.na(f$mut_A),"mut_A"] = 0
    f[is.na(f$mut_B),"mut_B"] = 0
    f = cbind(f, x)
    return(f);
}


train.model<-function(challenge) {
    # Get synergy training data
    f = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", challenge)
    # Create expression and mutation based features
    f = create.features(f, is.train=T)

    indices = which(colnames(f) %in% c("gexp_A", "gexp_B", "mut_A", "mut_B", "cat"))
    indices = c(indices, which(colnames(f) %in% c("med", "mean", "sd", "max", "min")))
    #indices = which(colnames(f) %in% c("med", "mean", "sd", "max", "min", "cat")) # only guild
    # Use all features
    if(challenge == "ch1a") {
	inTrain = createDataPartition(y = f$cat, p = 0.7, list=F) 
	training = f[inTrain, c(4:11, indices)]
	testing = f[-inTrain, c(4:11, indices)]
    # Use anything except monotherapy data
    } else if(challenge == "ch1b") {
	inTrain = createDataPartition(y = f$cat, p = 0.7, list=F) 
	training = f[inTrain, indices]
	testing = f[-inTrain, indices]
    } else {    
	stop("Unrecognized challenge!")
    }

    # Build model(s)
    # trainControl: boot for bootstrapping, cv for x-validation # repeatedcv for repeated 10-fold CV 
    ctrl = trainControl(method = "cv") #method = "repeatedcv", number = 10, repeats = 10)
    prep = c("center", "scale")
    # Random forest
    rfFit = train(cat ~ ., data=training, method = "rf", preProcess = prep, trControl = ctrl)
    pred = predict(rfFit, testing)
    a = cor(pred, testing$cat) 
    print(a)
    # Tree boost
    gbmFit = train(cat ~ ., data=training, method = "gbm", preProcess = prep, trControl = ctrl)
    pred = predict(gbmFit, testing)
    a = cor(pred, testing$cat) 
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
    a = cor(pred, testing$cat) 
    print(a)

    # Plot prediction consistency
    p = qplot(pred, cat, data=testing, main = paste0("Prediction PCC: ", format(a, digits=2)))
    png(paste0(output.dir, challenge, "/", "training.png")) 
    print(p)
    dev.off()

    return(list(rfFit=rfFit, gbmFit=gbmFit, modFit=modFit));
}


get.predictions<-function(challenge, rfFit, gbmFit, modFit, rebuild=F) {

    # Get test data
    f = get.synergy.data("Drug_synergy_data/ch1_leaderBoard_monoTherapy.csv", challenge)
    # Create expression and mutation based features
    f = create.features(f, is.train=F)

    indices = which(colnames(f) %in% c("gexp_A", "gexp_B", "mut_A", "mut_B", "cat"))
    indices = c(indices, which(colnames(f) %in% c("med", "mean", "sd", "max", "min")))
    # Use all features
    if(challenge == "ch1a") {
	#! fun.scale = function(x) { x$syn.einf * (1 + x$Einf_A + x$Einf_B) } #! do these need to be scaled / centered?
	fun.scale = function(x) { x$syn } 
	testing = f[, c(4:11, indices)]
    # Use anything except monotherapy data
    } else if(challenge == "ch1b") {
	fun.scale = function(x) { x$syn } 
	testing = f[, indices]
    } else {    
	stop("Unrecognized challenge!")
    }

    # Build models using all the training data
    if(rebuild) { 
	# Get training data
	f.training = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", challenge)
	if(challenge == "ch1a") {
	    training = f.training[, c(4:11, indices)]
	} else if(challenge == "ch1b") {
	    training = f.training[, indices]
	}
	ctrl = trainControl(method = "cv")
	prep = c("center", "scale")
	# Random forest
	rfFit = train(cat ~ ., data=training, method = "rf", preProcess = prep, trControl = ctrl)
	# Tree boost
	gbmFit = train(cat ~ ., data=training, method = "gbm", preProcess = prep, trControl = ctrl)
    }

    # Make predictions using model(s)
    pred.1 = predict(rfFit, testing)  
    #pred.1 = predict(rfFit, testing, type = "prob") # prob not meaningful for RF
    pred.2 = predict(gbmFit, testing) 
    pred.comb = data.frame(pred.1, pred.2)
    pred = predict(modFit, pred.comb)

    # Get predictions for leaderboard data (normalize with einf if using sny.einf)
    testing$syn.einf = pred
    testing$syn = pred
    testing$cat = fun.scale(testing)

    # Get confidence scores for learderboard data (assign lower confidence to values >= 10)
    testing$conf = 1-abs(testing$cat)/max(abs(testing$cat))
    f$PREDICTION = testing$cat
    f$CONFIDENCE = testing$conf

    # Output predictions
    file.name = paste0(output.dir, challenge, "/", "prediction.csv")
    write.csv(f[,c("CELL_LINE", "COMBINATION_ID", "PREDICTION")], file.name, row.names=F)
    file.name = paste0(output.dir, challenge, "/", "combination_priority.csv")
    a = f[,c("COMBINATION_ID", "CONFIDENCE")]
    a = aggregate(CONFIDENCE~COMBINATION_ID, a, mean)
    write.csv(a, file.name, row.names=F)
}


explore<-function() {

    # Get data
    f = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", "ch1a")
    #f = get.synergy.data("Drug_synergy_data/ch1_leaderBoard_monoTherapy.csv", "ch1a") # test
    # Create expression and mutation based features
    f = create.features(f, is.train=T)
    # Get expression and 
    d = get.expression.info(paste0(data.dir, "gex.csv"))
    e = get.mutation.info(paste0(data.dir, "mutations.csv"))

    indices = which(colnames(f) %in% c("gexp_A", "gexp_B", "mut_A", "mut_B", "cat"))
    indices = c(indices, which(colnames(f) %in% c("med", "mean", "sd", "max", "min")))
    training = f[,indices]

    # Check variance
    nzv = nearZeroVar(training, saveMetrics= TRUE)
    print(nzv) 

    # Check correlated features
    cor.mat = cor(training)
    cor.idx <- findCorrelation(cor.mat, cutoff = .75)
    training.filtered = training[,-cor.idx]
    
    # Imputing
    #training.processed <- predict(preProcess(training, method = c("center", "scale")), training)
    training.processed <- predict(preProcess(training, method = c("center", "scale", "knnImpute"), k=5), training.filtered) # "BoxCox"
    print(head(training.processed))

    return() 

    ggplot(data=e, aes(cell_line_name)) + geom_histogram(stat="bin") + coord_flip()

    table(e$cell_line_name, e$mut)
    ggplot(data=e, aes(cell_line_name)) + geom_histogram(stat="bin", group="mut") + coord_flip()

    ggplot(data=f, aes(COMPOUND_A, COMPOUND_B)) + geom_point(aes(color=factor(QA))) + theme(axis.text.x = element_text(angle=90))
    #g = f
    #g$cut = cut(g$SYNERGY_SCORE, breaks=c(1000, 50, 10, 5, -5, -10, -50, -1000))

    #h = melt(g[,c("H_A", "H_B", "SYNERGY_SCORE")], "SYNERGY_SCORE")
    ggplot(data=g, aes(Einf_A, Einf_B)) + geom_point(aes(color=factor(cut)))
    #artsy
    #ggplot(data=g, aes(Einf_A, Einf_B)) + geom_point(aes(color=factor(cut)), alpha=0.5) + guides(color=F) + theme(text=element_text(size=0), axis.ticks=element_line(size=NA), panel.background=element_rect(fill=NA), panel.border=element_rect(size=NA, fill=NA), panel.grid=element_line(size=NA), panel.grid.minor=element_line(size=NA))
    #f$ic.2 = f[,"SYNERGY_SCORE"] * (f[, "IC50_A"] * f[, "IC50_B"])
    #f$ic = f[,"SYNERGY_SCORE"] * (1 + f[, "IC50_A"] + f[, "IC50_B"])
    
    # Check synergy consistency in the same cell line (61 cases with multiple instances)
    comb.count = table(f$COMBINATION_ID, f$CELL_LINE)
    k = which(comb.count>1, arr.ind=T)
    y.1 = data.frame(syn.1=NA, syn.2=NA)
    y.2 = data.frame(syn.1=NA, syn.2=NA)
    y = data.frame(ic.1=NA, ic.2=NA)
    for(i in 1:nrow(k)) {
	comb = rownames(comb.count)[k[i, 1]]
	cell = gsub("\\.", "-", colnames(comb.count)[k[i, 2]])
	#print(c(comb, cell))
	x = f[f$COMBINATION_ID == comb & f$CELL_LINE==cell,]
	for(m in 1:(nrow(x)-1)) {
	    for(n in (m+1):nrow(x)) {
		y.1 = rbind(y.1, c(x[m, "syn"], x[n, "syn"]))
		#y.1 = rbind(y.1, c(x[m, "ic"], x[n, "ic"]))
		y.2 = rbind(y.2, c(x[m, "syn.einf"], x[n, "syn.einf"]))
		#y.2 = rbind(y.2, c(x[m, "ic.2"], x[n, "ic.2"]))
		y = rbind(y, c(x[m, "IC50_A"], x[n, "IC50_A"]))
		y = rbind(y, c(x[m, "IC50_B"], x[n, "IC50_B"]))
	    }
	}
    }

    a = cor(y, use="complete", method="spearman")
    print(a[1,2])
    #png(paste(output.dir, "synergy_consistency.png", sep="")) #, width=10, height=6)
    #par(mfrow=c(2,1))
    a = cor(y.1, use="complete", method="spearman")
    print(a[1,2])
    #plot(y.1, main = paste("Raw (q: ", format(a[1,2], digits=2), ")", sep=""))
    a = cor(y.2, use="complete", method="spearman")
    print(a[1,2])
    #plot(y.2, main = paste("Normalized (q: ", format(a[1,2], digits=2), ")", sep=""))
    #dev.off()

    # Check synergy consistency in similar cell lines   
    selection.function<-function(col.values, mapping) {
	#tapply(col.values, mapping, function(x) { x[which.max(abs(x))] })
	tapply(col.values, mapping, function(x) { mean(x) })
    }
    arr.cor = cor(d)
    arr.cor[upper.tri(arr.cor, diag=T)] = NA
    k = which(arr.cor>0.9, arr.ind=T)
    for(i in 1:nrow(k)) {
	cell.1 = gsub("\\.", "-", rownames(arr.cor)[k[i, 1]])
	cell.2 = gsub("\\.", "-", colnames(arr.cor)[k[i, 2]])
	#cell.1 = "MDA-MB-231"
	#cell.2 = "T-24"
	print(c(cell.1, cell.2))
	x = f[f$CELL_LINE==cell.1,]
	y = f[f$CELL_LINE==cell.2,]
	comb.common = intersect(x$COMBINATION_ID, y$COMBINATION_ID)
	if(length(comb.common) == 0) {
	    next();
	}
	a = x[x$COMBINATION_ID %in% comb.common, ]
	b = y[y$COMBINATION_ID %in% comb.common, ]
	a = a[order(a$COMBINATION_ID),]
	b = b[order(b$COMBINATION_ID),]
	#print(a)
	#print(b)
	# average multiple combinations 
	mapping = factor(a$COMBINATION_ID)
	#a = apply(a[4:12], 2, function(val) { selection.function(val, mapping)})
	a = data.frame(lapply(a[4:12], function(val) { selection.function(val, mapping)}))
	mapping = factor(b$COMBINATION_ID)
	#b = apply(b[4:12], 2, function(val) { selection.function(val, mapping)})
	b = data.frame(lapply(b[4:12], function(val) { selection.function(val, mapping)}))
	#print(a)
	#print(b)
	par(mfrow=c(1,2))
	vals.a = a[,"SYNERGY_SCORE"]
	vals.b = b[,"SYNERGY_SCORE"]
	plot(vals.a, vals.b)
	print(cor(vals.a, vals.b))#, method="spearman"))
	vals.a = a[,"SYNERGY_SCORE"] / (1+a[,"Einf_A"]+a[,"Einf_B"])
	vals.b = b[,"SYNERGY_SCORE"] / (1+b[,"Einf_A"]+b[,"Einf_B"])
	plot(vals.a, vals.b)
	print(cor(vals.a, vals.b))#, method="spearman"))
	#readline(prompt = "Pause. Press <Enter> to continue...")
    }
}


old<-function() {
    set.seed(3456)
    library(iris)
    adata = iris
    inTrain = createDataPartition(y = adata$Species, p = 0.7, list=F) #[[1]]
    training = adata[inTrain,]
    testing = adata[-inTrain,]
    # GLM w/ PCA preprocessing
    m = train(training, training$diagnosis, preProcess="pca", method="glm")
    a = preProcess(training, method="pca", thresh=0.8)
    trainPC = predict(a, training)
    b = preProcess(testing, method="pca")
    testPC = predict(b, testing)
    confusionMatrix(testing$diagnosis, predict(m3, testPC))

    methods = c("rf", "gbm", "glm")
    for(method in methods) {
    print(method)
    if(method == "rf") {
	# Random forest
	m = train(variable.formula, data = training, method = "rf", prox = T) 
	# getTree for individual trees
	# rfcv
    } else if(method == "gbm") {
	# Boost w/ trees
	m = train(variable.formula, data = training, method = "gbm", verbose = F) 
    } else if(method == "glm") {
	m = train(variable.formula, data = training, method ="glm", preProcess="pca")
    }
    pred = predict(m, testing)
    print(pred)
    print(names(pred))
    testing$pred = pred

    variable = "Species"
    testing$flag = pred == testing[[variable]]
    a = table(pred, testing[[variable]])
    print(a)

    a = confusionMatrix(pred, testing[[variable]])
    print(a)

    pred = predict(m, training)
    a = confusionMatrix(pred, training[[variable]])
    print(a)
    }
    #p = ggplot(data=testing, aes_string("pred", variable)) + geom_point() 
    #print(p)
    #p = qplot(Petal.Width, Petal.Length, colour=flag, data=testing)
    #print(p)
    return()
}

main()




