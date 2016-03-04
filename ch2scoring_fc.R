library(ROCR)
# ------------------------------------------------------------------------------------
# Description: AZ-Sanger Challenge scoring functions
# Authors: Michael P Menden, Julio Saez-Rodriguez, Mike Mason, Thomas Yu
# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
# Get observation format for Subchallenge 2
# ------------------------------------------------------------------------------------
getObs_ch2 <- function(ls,threshold = F) {
  combi <- unique(ls$COMBINATION_ID)
  cell <- unique(ls$CELL_LINE)
  mat <- matrix(NA, nrow=length(combi), ncol=length(cell),
                dimnames=list(combi, cell))
  for (i in 1:nrow(ls)) 
    mat[as.character(ls$COMBINATION_ID[i]), 
        as.character(ls$CELL_LINE[i])] <- ls$SYNERGY_SCORE[i]
  
  if (is.numeric(threshold)) {
    mat <- signif(mat,3)
    mat[mat < threshold] = 0
    mat[mat >= threshold ] = 1
  }
  return(mat)
}

getObs_ch2_matrix <- function(obs) {
  if (all(row.names(obs) == c(1:nrow(obs)))) {
    row.names(obs) = obs[,1]
    obs = obs[,-1]
  }
  obs <- as.matrix(obs) 
  return(obs)
}

# ------------------------------------------------------------------------------------
# Get prediction format for Subchallenge 2
# ------------------------------------------------------------------------------------
getPred_ch2 <- function(pred) {
  if (all(row.names(pred) == c(1:nrow(pred)))) {
    row.names(pred) = pred[,1]
    pred = pred[,-1]
  }
  pred <- as.matrix(pred) 
  return(pred)
}

# ------------------------------------------------------------------------------------
# Get unsigned score from one dimensional ANOVA
# ------------------------------------------------------------------------------------
getNegLog10pVal_ch2 <- function(fit, obs) {
  s <- 0
  if (!is.na(fit$coefficients[2]) & sum(!is.na(obs)) > 2)
      s <- -log10(anova(fit)['pred','Pr(>F)'])
  return(s)
}
# -----------------------------------------------------------------------
# Scores from confusion Matrix
# -----------------------------------------------------------------------
getPrecision_ch2 <- function(obs, pred, threshold=30) {
  obs <- read.csv(obs)
  obs <- getObs_ch2(obs,threshold)
  
  pred <- read.csv(pred,stringsAsFactors=F,check.names = F)
  pred <- getPred_ch2(pred)  
  pred <- pred[match(row.names(obs),row.names(pred)),]
  pred <- pred[,match(colnames(obs),colnames(pred))]

  #Remove all NA's
  pred <- as.numeric(pred)[!is.na(obs)]
  obs <- as.numeric(obs)[!is.na(obs)]

  preds <- prediction(pred,obs)
  prec <- performance(preds,"prec") #precision (Acc + )
  sens <- performance(preds,"sens") #True positive rate (Sensitivity) (Cov +)
  npv <- performance(preds,"npv") #Negative predictive value (Acc - )
  spec <- performance(preds,"spec") #True negative rate(specificity) (Cov -)
  auc <- performance(preds,"auc") #Area under curve (AUC)
  phi <- performance(preds,"phi") #phi correlation coefficient, (matthews)
  aupr <- performance(preds, "prec", "rec") #Area under precision recall (AUPR)
  
  prec_val <- unlist(prec@y.values)[2]
  sens_val <- unlist(sens@y.values)[2]
  npv_val <- unlist(npv@y.values)[2]
  spec_val <- unlist(spec@y.values)[2]
  auc_val <- unlist(auc@y.values)
  phi_val <- unlist(phi@y.values)[2]
  BAC <- (sens_val + spec_val)/2
  F1 <- 2*preds@tp[[1]][2]/(2*preds@tp[[1]][2] + preds@fn[[1]][2] + preds@fp[[1]][2])
  aupr_val <- unlist(aupr@y.values)[2]
  
  #If predictions are 0, ROCR treats 0 as the positive value when it is actually the negative
  if (all(pred == 0)) {
    prec_val <- unlist(prec@y.values)[1]
    sens_val <- unlist(sens@y.values)[1]
    npv_val <- unlist(npv@y.values)[1]
    spec_val <- unlist(spec@y.values)[1]
    auc_val <- unlist(auc@y.values)
    phi_val <- unlist(phi@y.values)[1]
    BAC <- (sens_val + spec_val)/2
    F1 <- 2*preds@tp[[1]][1]/(2*preds@tp[[1]][1] + preds@fn[[1]][1] + preds@fp[[1]][1])
    aupr_val <- unlist(aupr@y.values)[1]
  }
  
  return(round(c(prec=prec_val,
                 sens = sens_val,
                 npv = npv_val,
                 spec=spec_val,
                 auc=auc_val,
                 phi=phi_val,
                 BAC=BAC, # Tie breaking Metric
                 F1=F1,                 
                 aupr=aupr_val),2))

}

# ------------------------------------------------------------------------------------
# Get the drug combinations score of Subchallenge 2
# ------------------------------------------------------------------------------------
getOneDimScore_ch2 <- function(obs, pred, confidence="none", topX=10, rows=T) {
  obs <- read.csv(obs)
  obs <- getObs_ch2(obs)

  pred <- read.csv(pred,stringsAsFactors=F,check.names = F)
  pred <- getPred_ch2(pred)

  # Remove potentially filtered cell lines from the training set
  obs = obs[,which(colnames(obs) %in% colnames(pred))]
  pred <- pred[match(row.names(obs),row.names(pred)),]
  pred <- pred[,match(colnames(obs),colnames(pred))]
  #pred = pred[,-which(is.na(colnames(pred)))] # removing NA column (potentially filtered cell lines from the training set)

  #if(rows) {
      # Remove potentiall filtered combinations
      #indices = which(apply(pred, 1, function(x){all(is.na(x))}))
      #pred = pred[-indices,]
      #indices = which(apply(obs, 1, function(x){all(is.na(x))}))
      #obs = obs[-indices,]
      # Remove rows without any matching non NA cell lines
      #indices = apply(pred+obs, 1, function(x){all(is.na(x))}) # Generates NAs
      #pred = pred[-indices,]
      #obs = obs[-indices,]
  #}

  n <- ncol(obs)
  if (rows)
    n <- nrow(obs)
  
  s <- c()
  for (i in 1:n) {
    sign <- 1
    if (rows) {
      if(all(is.na(pred[i,]+obs[i,]))) { next(); } # Skip if no overlappin non-NA values
      fit <- aov(obs[i,] ~ pred[i,])
      nlp <- getNegLog10pVal_ch2(fit,obs[i,])
      if(is.na(nlp)) { next(); } # Skip if the value is NA
        if (nlp!=0 & (mean(obs[i, pred[i,]==1], na.rm=T) < mean(obs[i, pred[i,]==0], na.rm=T)))
          sign <- -1
    } else {
      fit <- aov(obs[,i] ~ pred[,i]) 
      nlp <- getNegLog10pVal_ch2(fit,obs[,i])
      if(is.na(nlp)) { next(); } # Skip if the value is NA
        if (nlp!=0 & (mean(obs[pred[,i]==1, i], na.rm=T) < mean(obs[pred[,i]==0, i], na.rm=T)))
          sign <- -1
    }
    s <- c(s, sign * nlp)
  }
  
  if (!file.exists(confidence))
    return(round(c(mean=mean(s),
             ste=sd(s)),5))
  
  confidence <- read.csv(confidence,stringsAsFactors=F,check.names = F)
  confidence <- getPred_ch2(confidence)
  confidence <- confidence[match(row.names(obs),row.names(confidence)),]
  confidence <- confidence[,match(colnames(obs),colnames(confidence))]
  
  if (rows) {
    nVal <- round(topX * (nrow(confidence) / 100))
  } else {
    nVal <- round(topX * (ncol(confidence) / 100))
  }
  
  nStep <- 1000
  boot_score <- rep(0, nVal)
  for (i in 1:nStep) {
    if (rows) {
      avgConf <- sapply(1:nrow(confidence), function(x) mean(confidence[x, !is.na(obs[x,])]))
    } else {
      avgConf <- sapply(1:ncol(confidence), function(x) mean(confidence[!is.na(obs[,x]), x]))
    }
    idx <- order(avgConf, sample(length(avgConf)), decreasing = T)[1:nVal]
    boot_score <- boot_score + s[idx]
  }
  
  return(round(c(mean=mean(boot_score/nStep),
                 ste=sd(boot_score/nStep)),5))
}

# ------------------------------------------------------------------------------------
# Get the performance score of Subchallenge 2
# ------------------------------------------------------------------------------------
getGlobalScore_ch2 <- function(obs, pred) { 
  obs <- read.csv(obs)
  obs <- getObs_ch2(obs)
  
  pred <- read.csv(pred,stringsAsFactors=F,check.names = F)
  pred <- getPred_ch2(pred)
  pred <- pred[match(row.names(obs),row.names(pred)),]
  pred <- pred[,match(colnames(obs),colnames(pred))]

  # regress out combination bias
  cov <- rep(rownames(obs), ncol(obs))
  
  c0 <- rep(rownames(obs), ncol(obs))
  c1 <- as.vector(matrix(colnames(obs), ncol=ncol(obs), nrow=nrow(obs), byrow=T))
  
  obs <- as.vector(obs)
  pred <- as.vector(pred)
  
  c0 <- c0[!is.na(obs)]
  c1 <- c1[!is.na(obs)]
  pred <- pred[!is.na(obs)]
  obs <- obs[!is.na(obs)]
  
  if(all(pred==0) | all(pred==1))
    return(0)

  # run anove with combination label as covariate
  fit <- aov(obs ~ c0 + c1 + pred)
  pVal <- -log10(anova(fit)['pred','Pr(>F)'])
  
  sign <- 1

  if (sum(!is.na(obs[pred==1])) >0  && sum(!is.na(obs[pred==0]))>0)
    if (mean(obs[pred==1], na.rm=T) < mean(obs[pred==0], na.rm=T))
      sign <- -1
  
  return(round(sign * pVal,2)) # Final Metric
}
