# To remove the prompt (does not work with Jupyter)
options(prompt=" ", continue=" ")

# Get parameters
initialize<-function() {
    parameters = data.frame(data.dir = "../data/", output.dir = "../doc/")
    #set.seed(142341)
    file.name = "Drug_synergy_data/ch1_train_combination_and_monoTherapy.csv"
    parameters$synergy.file = paste0(parameters$data.dir, file.name)
    parameters$results.file = paste0(parameters$output.dir, "results.dat")
    return(parameters);
}

# Data overview
overview<-function(parameters, summarize=F) {
    f = read.csv(parameters$synergy.file)
    if(summarize == T) {
	print(summary(f))
    }
    print(sprintf("Number of samples: %d", nrow(f))) 
    print(sprintf("Number of drug pairs: %d, number of cell lines: %d", length(levels(f$COMBINATION_ID)), length(levels(f$CELL_LINE))))
    print(sprintf("Features: %s", paste(colnames(f), collapse=", ")))
    return(f);
}

# Filter low quality / insensitive samples / correlated features
filter<-function(f, cutoff) {
    # Low quality
    f$cat = f[,"SYNERGY_SCORE"]
    n.qa = nrow(f[f$QA!=1,])
    f = f[f$QA==1,]

    # Insensitive
    require(plyr)
    #cutoff = 40 # 20 (8 cell lines) # 10 (22 cell lines)
    d = f 
    d$einf = (d$Einf_A + d$Einf_B) / 2
    a = ddply(d, ~ CELL_LINE, summarize, syn.sd = sd(cat), syn.min = min(cat), syn.med = median(cat), einf.sd = sd(einf), einf.min = min(einf), einf.med = median(einf))
    # 22RV1 CAL-120 HCC1143 HCC1428 HCC1937 KU-19-19 UACC-812 VCaP
    # Cells with high min Einf has lower synergy
    b = cor.test(a$syn.med, a$einf.min, use="complete")
    cell.lines = as.vector(a[a$einf.min>cutoff,"CELL_LINE"])
    n.einf = nrow(f[f$CELL_LINE %in% cell.lines,])
    f = f[!f$CELL_LINE %in% cell.lines,] 

    require(caret)
    # Find correlated features
    cor.mat = cor(f[,4:11])
    cor.idx = findCorrelation(cor.mat, cutoff = .75)
    print(c("Correlated features:", colnames(f)[cor.idx]))
    print(sprintf("Correlation between einf.min and syn.med: %f %f", b$estimate[[1]], b$p.value)) # -0.235
    print(sprintf("Insensitive cell lines: %s", paste(cell.lines, collapse=", ")))
    print(sprintf("Number of samples with QA < 1: %d, Einf > %d: %d", n.qa, cutoff, n.einf))
    return(f);
}


results<-function(parameters) {
    d = read.csv(parameters$results.file, sep="\t", header=T)
    d
}


old<-function() {
    # Baseline predictions (using monotherapy response data)
    testing = read.table(paste0(output.dir, "observed.dat"))
    pred = read.table(paste0(output.dir, "/ch1a_base/", "prediction.dat"))
    a = cor(pred, testing)
    # Correlation between predicted and observed synergy scores for 10-fold cross-validation
    print(sprintf("Correlation: %.2f", a[[1]]))

    # Baseline + Expression based predictions
    pred = read.table(paste0(output.dir, "/ch1a_gexp/", "prediction.dat"))
    a = cor(pred, testing)
    print(sprintf("Correlation: %.2f", a[[1]]))

    # Baseline + Mutation based predictions
    pred = read.table(paste0(output.dir, "/ch1a_mut/", "prediction.dat"))
    a = cor(pred, testing)
    print(sprintf("Correlation: %.2f", a[[1]]))

    # Baseline + GUILD based predictions
    pred = read.table(paste0(output.dir, "/ch1a_guild/", "prediction.dat"))
    a = cor(pred, testing)
    print(sprintf("Correlation: %.2f", a[[1]]))

    options(repr.plot.width=5, repr.plot.height=5)
    # Baseline + Mutation + GUILD based predictions
    pred = read.table(paste0(output.dir, "/ch1a_mut_guild/", "prediction.dat"))
    a = cor(pred, testing)
    print(sprintf("Correlation: %.2f", a[[1]]))
    p = qplot(pred$x, testing$x, main = paste0("Prediction PCC: ", format(a, digits=2)))
    print(p)

    # Mutation + GUILD based predictions
    pred = read.table(paste0(output.dir, "/ch1b_mut_guild/", "prediction.dat"))
    a = cor(pred, testing)
    print(sprintf("Correlation: %.2f", a[[1]]))
    p = qplot(pred$x, testing$x, main = paste0("Prediction PCC: ", format(a, digits=2)))
    print(p)
}

