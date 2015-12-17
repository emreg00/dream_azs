library(caret)
library(ggplot2)
library(reshape2)

data.dir = "../data/"
output.dir = "../doc/"

main<-function() {
    set.seed(142341)
    explore()
    synergy.consistency()
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


create.features<-function(f, challenge, is.train) {
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
	f[i, "mut_A"] = mean(e[e$gene %in% a & e$cell_line_name == cell.line,"mut"], rm.na=T)
	f[i, "mut_B"] = mean(e[e$gene %in% b & e$cell_line_name == cell.line,"mut"], rm.na=T)
	#f[i, "mut_A"] = mean(c(f[i, "mut_A"], f[i, "mut_B"]), rm.na=T) #
	a = h[h$comb.id == comb.id & h$cell.line == cell.line, c("med", "mean", "sd", "max", "min")] # "max.a", "max.b", 
	if(nrow(a) == 0) {
	    #a[1,] = 0
	    a[1,] = NA
	}
	x = rbind(x, a)
    }
    # Impute NA values using kNN (during training)
    #f[is.na(f$gexp_A),"gexp_A"] = 0
    #f[is.na(f$gexp_B),"gexp_B"] = 0
    #f[is.na(f$mut_A),"mut_A"] = 0
    #f[is.na(f$mut_B),"mut_B"] = 0
    f = cbind(f, x)
    return(f);
}


temp.process.features<-function(f, challenge, is.train=T) {
    # Choose features to include
    #indices = which(colnames(f) %in% c("cat"))
    indices = which(colnames(f) %in% c("mut_A", "mut_B", "cat"))
    #f$mut = (f$mut_A + f$mut_B) / 2
    #indices = which(colnames(f) %in% c("mut", "cat"))
    #indices = which(colnames(f) %in% c("gexp_A", "gexp_B", "mut_A", "mut_B", "cat"))
    #indices = c(indices, which(colnames(f) %in% c("med", "sd", "min")))
    #indices = which(colnames(f) %in% c("med", "mean", "sd", "max", "min", "cat")) # only guild
    #f = f[, c(1:8, indices)]
    f = f[, indices]
    # Imputing
    f = predict(preProcess(f, method = c("knnImpute"), k=5), f) # "BoxCox"
    return(f);
}


explore<-function() {
    # Get data and expression and mutation based features
    f = get.synergy.data("Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv", "ch1a", is.train=T)
    #f = get.synergy.data("Drug_synergy_data/ch1_leaderBoard_monoTherapy.csv", "ch1a", is.train=F) # test

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

    # Feature selection
    ctrl <- rfeControl(functions=rfFuncs, method="cv", number=10)
    # run the RFE algorithm
    results <- rfe(training[,-which(colnames(training)=="cat")], training$cat, sizes=c(1:(ncol(training)-1)), rfeControl=ctrl)
    print(results)
    # list the chosen features
    predictors(results)
    plot(results, type=c("g", "o"))

    # For categorical data
    #ctrl <- trainControl(method="repeatedcv", number=10, repeats=3)
    #model <- train(cat~., data=training, method="lvq", preProcess="scale", trControl=ctrl)
    # estimate variable importance
    #importance <- varImp(model, scale=FALSE)
    #print(importance)
    #plot(importance)

    # Explore mutation data
    e = get.mutation.info(paste0(data.dir, "mutations.csv"))
    ggplot(data=e, aes(cell_line_name)) + geom_histogram(stat="bin") + coord_flip()

    table(e$cell_line_name, e$mut)
    ggplot(data=e, aes(cell_line_name)) + geom_histogram(stat="bin", group="mut") + coord_flip()
    return() 
}


check.synergy.consistency<-function() {
    ggplot(data=f, aes(COMPOUND_A, COMPOUND_B)) + geom_point(aes(color=factor(QA))) + theme(axis.text.x = element_text(angle=90))
    #g = f
    #g$cut = cut(g$SYNERGY_SCORE, breaks=c(1000, 50, 10, 5, -5, -10, -50, -1000))

    #h = melt(g[,c("H_A", "H_B", "SYNERGY_SCORE")], "SYNERGY_SCORE")
    ggplot(data=g, aes(Einf_A, Einf_B)) + geom_point(aes(color=factor(cut)))
    #artsy
    #ggplot(data=g, aes(Einf_A, Einf_B)) + geom_point(aes(color=factor(cut)), alpha=0.5) + guides(color=F) + theme(text=element_text(size=0), axis.ticks=element_line(size=NA), panel.background=element_rect(fill=NA), panel.border=element_rect(size=NA, fill=NA), panel.grid=element_line(size=NA), panel.grid.minor=element_line(size=NA))
    #f$ic = f[,"SYNERGY_SCORE"] * (1 + f[, "IC50_A"] + f[, "IC50_B"])
    f$syn = f[,"SYNERGY_SCORE"]
    f$syn.einf = f[,"SYNERGY_SCORE"] / (1 + f[, "Einf_A"] + f[, "Einf_B"])

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

    # Get expression and mutation data
    d = get.expression.info(paste0(data.dir, "gex.csv"))
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


