library(dsda)
library(randomForest)
library(MASS)
library(hdqr)
library(PRROC)

stan = function(train, validation=NULL, test=NULL) {
  # standardize data set
  # Args:
  #   train:      original training set
  #   validation: original validation set
  #   test:       original test set
  # Returns:
  #   train:      standardized training set
  #   validation: standardized validation set
  #   test:       standardized test set
  train_m     = colMeans(train$X)
  std_train_x = t(apply(train$X, 1, function(x) x - train_m))
  train_sd    = apply(std_train_x, 2,
                function(x) sqrt(x %*% x / length(x)))
  train_sd[train_sd==0] = 1
  train$X = t(apply(std_train_x, 1, function(x) x / train_sd))
  if (!is.null(validation)) validation$X =
    scale(validation$X, center=train_m, scale=train_sd)
  if (!is.null(test)) test$X =
    scale(test$X, center=train_m, scale=train_sd)
  rm.att = function(x) {
    attributes(x) = attributes(x)[c(1,2)]
    x
  }
  train$X = rm.att(train$X)
  validation$X = rm.att(validation$X)
  test$X = rm.att(test$X)
  # returns:
  list(train=train, validation=validation, test=test)
}



hdqr_z_rep <- function(train, test, q, tau, lam2, l_list, rep=10, seed){
  x <- as.matrix(train$x)
  nvars <- dim(x)[2]
  n <- dim(x)[1]
  y <- train$y
  new_x <- as.matrix(test$x)
  ntau <- length(tau)
  coefmat <- array(NA, c((nvars+1), rep, ntau))
  for (r in 1:rep){
    set.seed(seed+r)
    u <- runif(length(y))
    z <- y + u
    for (i in 1:ntau){
      set.seed(seed+r)
      cv.fit <- cv_z(x, y, z, lambda = l_list, lam2=lam2, tau = tau[i], q=q)
      lam <- cv.fit$lambda.min
      lam_list <- seq(lam+2, lam, length.out=5)
      fit <- hdqr(x, z, lambda = lam_list, lam2=lam2, tau=tau[i])
      coefs <- coef(fit, s=cv.fit$lambda.min)
      coefmat[,r,i] <- coefs[,1]
    }
  }
  coef_avg <- apply(coefmat, c(1, 3), mean)
  pred <- apply(coef_avg, 2, function(x) as.matrix(as.matrix(cbind2(1, new_x)) %*% x)-1)
  aucvec <- apply(pred, 2, function(x) roc(as.factor(test$y), as.vector(x))$auc)
  ind <- which(aucvec==max(aucvec))[1]
  predclass <- ifelse(pred[, ind]>q,1,0)
  predclass <- factor(predclass, levels = levels(factor(test$y)))
  resmat <- confusionMatrix(as.factor(predclass),as.factor(test$y), positive="1")
  res <- resmat$byClass
  tab <- resmat$table
  #recall precision f1 auc
  prauc <- pr.curve(scores.class0 = pred[test$y==1, ind], 
      scores.class1 = pred[test$y==0, ind])[[2]]
  recall <- ifelse(is.na(res["Recall"]), 0, res["Recall"])
  gmean <- sqrt(recall * res["Specificity"])
  res_list <- c(recall, ifelse(is.na(res[5]), 0, res[5]),
    ifelse(is.na(res[7]), 0, res[7]), aucvec[ind], prauc, gmean)
  list(tau=tau, res=res_list, tau.min=tau[ind], q=q, table=tab)
}


logistic_class <- function(train, test, wts, lambda, lam2, seed){  
  set.seed(seed)
  cv.fit <- cv_log(train$x, train$y, lambda=lambda, lam2=lam2)
  lam <- cv.fit$lambda.min
  alp <- lam/(lam+lam2)
  # lam_list <- c(1, lam+lam2)
  fit <- glmnet(train$x, train$y, family="binomial", weights=wts, 
     lambda=lam+lam2, alpha=alp)
  pred <- stats::predict(fit, test$x, type="response", s=lam) 
  predclass <- stats::predict(fit, test$x, type="class", s=lam) 
  pred <- as.matrix(pred)
  res <- confusionMatrix(as.factor(predclass),as.factor(test$y), positive="1")
  log_recall <- ifelse(is.na(res$byClass["Recall"]), 0, res$byClass["Recall"])
  log_precision <- ifelse(is.na(res$byClass["Precision"]), 0, res$byClass["Precision"])
  log_F1 <- ifelse(is.na(res$byClass["F1"]), 0, res$byClass["F1"])
  log_auc <- roc(as.factor(test$y), as.vector(pred))$auc
  log_prauc <- pr.curve(scores.class0 = pred[test$y==1], 
      scores.class1 = pred[test$y==0])[[2]]
  log_gmean <- sqrt(log_recall * res$byClass["Specificity"])
  list(res=c(log_recall, log_precision, log_F1, log_auc, log_prauc, log_gmean),
    table = res$table)
}

rf_class <- function(train, test, wts, seed){
  # wts <- 100 / table(train$y)
  set.seed(seed)
  classifier <- randomForest(x=train$x, y=as.factor(train$y), weights = wts)
  pred <- stats::predict(classifier, test$x, type="prob") # predict
  predclass <- stats::predict(classifier, test$x, type="response")
  res <- confusionMatrix(as.factor(predclass),as.factor(test$y), positive="1")
  rf_recall <- ifelse(is.na(res$byClass["Recall"]), 0, res$byClass["Recall"])
  rf_precision <- ifelse(is.na(res$byClass["Precision"]), 0, res$byClass["Precision"])
  rf_F1 <- ifelse(is.na(res$byClass["F1"]), 0, res$byClass["F1"])
  rf_auc <- roc(as.factor(test$y), as.vector(pred[,2]))$auc
  rf_prauc <- pr.curve(scores.class1 = pred[test$y==0,2],
    scores.class0 = pred[test$y==1,2])[[2]]
  rf_gmean <- sqrt(rf_recall * res$byClass["Specificity"])
  list(res=c(rf_recall, rf_precision, rf_F1, rf_auc, rf_prauc, rf_gmean),
    table = res$table)
}


dsda_class <- function(train, test, lambda, lam2, seed){ 
  set.seed(seed)
  cv.fit <- cv_dsda(train$x, train$y, lambda=lambda, lam2=lam2)
  lam <- cv.fit$bestlambda
  alp <- lam/(lam+lam2)
  dsda.fit <- dsda.path(train$x, train$y, lambda=lam+lam2, alpha=alp)
  pred <- predict.dsda(dsda.fit, test$x, type="response") 
  pred <- as.matrix(pred)
  predclass <- predict.dsda(dsda.fit, test$x, type="class") 
  res <- confusionMatrix(as.factor(predclass),as.factor(test$y), positive="1")
  recall <- ifelse(is.na(res$byClass["Recall"]), 0, res$byClass["Recall"])
  precision <- ifelse(is.na(res$byClass["Precision"]), 0, res$byClass["Precision"])
  F1 <- ifelse(is.na(res$byClass["F1"]), 0, res$byClass["F1"])
  auc <- roc(as.factor(test$y), as.vector(pred))$auc
  prauc <- pr.curve(scores.class0 = pred[test$y==1], 
      scores.class1 = pred[test$y==0])[[2]]
  gmean <- sqrt(recall * res$byClass["Specificity"])
  list(res=c(recall, precision, F1, auc, prauc, gmean), table = res$table)
}


smote_class <- function(train, test, seed){
  k <- sum(train$y==1)
  K <- ifelse(k<5, max(0, k-1), 5)
  p <- dim(train$x)[2]+1
  set.seed(seed)
  smote_data <- smotefamily::SMOTE(data.frame(train$x), train$y, K=K)
  set.seed(seed)
  classifier <- randomForest(x=smote_data$data[,-p], y=as.factor(smote_data$data[,p]))
  pred <- stats::predict(classifier, test$x, type="prob") # predict
  predclass <- stats::predict(classifier, test$x, type="response")
  res <- confusionMatrix(as.factor(predclass),as.factor(test$y), positive="1")
  recall <- ifelse(is.na(res$byClass["Recall"]), 0, res$byClass["Recall"])
  precision <- ifelse(is.na(res$byClass["Precision"]), 0, res$byClass["Precision"])
  F1 <- ifelse(is.na(res$byClass["F1"]), 0, res$byClass["F1"])
  auc <- roc(as.factor(test$y), as.vector(pred[,2]))$auc
  prauc <- pr.curve(scores.class1 = pred[test$y==0,2],
    scores.class0 = pred[test$y==1,2])[[2]]
  gmean <- sqrt(recall * res$byClass["Specificity"])
  smote_res <- c(recall, precision, F1, auc, prauc, gmean)
  list(res=smote_res,
    table = res$table)
}

quantileDA_class <- function(train, test, seed){
  set.seed(seed)
  out <- quantileCV(train$x, train$y)
  cltest <- quantilecl(train$x, test$x, train$y, theta = out$theta.choice, 
    cl.test = test$y, skew.correct="Galton")$cl.test
  res <- confusionMatrix(as.factor(cltest), as.factor(test$y), positive="1")
  recall <- ifelse(is.na(res$byClass["Recall"]), 0, res$byClass["Recall"])
  precision <- ifelse(is.na(res$byClass["Precision"]), 0, res$byClass["Precision"])
  F1 <- ifelse(is.na(res$byClass["F1"]), 0, res$byClass["F1"])
  auc <- NA
  prauc <- NA
  gmean <- sqrt(recall * res$byClass["Specificity"])
  quan_res <- c(recall, precision, F1, auc, prauc, gmean)
  list(res=quan_res, table = res$table)  
}


DataGen = function(sampsiz, pp=3000, p0=5,
           mu=0.7, rho=0.7, blanc=0.9, sina=3) {
  # Generate the simulation examples used in Wang et al (2006)
  #
  # Args:
  #   sampsiz: sample size of the datasets (i.e., n)
  #   pp:      dimensions of the datasets (i.e., p)
  #   p0:      theoretically significant variables
  #   mu:      mean shift
  #   rho:     correlation used in covariance matrix
  #   blanc:   blance of +1 and -1 classes
  #   sina:    one of three sinarios used in the paper.
  # Returns:
  #   generated datasets

  # The means of two classes.
  mup = rep(mu, p0)
  mum = rep(-mu, p0)
  # The var-cov matrix of two classes.
  if(sina == 1){
    sigma = diag(p0)
  }
  if(sina == 2){
    sigma = diag(p0)
    for(i in 1:p0) for(j in 1:p0) 
      sigma[i,j] = ifelse(i == j, 1, rho)
  }
  if(sina == 3){
    sigma = diag(p0)
    for(i in 1:p0) for(j in 1:p0) 
      sigma[i,j] = rho ^ abs(i - j)
  }
  eo = eigen(sigma, symmetric=TRUE)
  sigma.sqrt = tcrossprod(eo$vec %*% (diag(sqrt(eo$val))), eo$vec)

  dat = NULL
  dat$X = matrix(NA, sampsiz, pp)
  dat$y = c(1,0)[factor(rbinom(sampsiz, 1, blanc))]

  for(i in seq.int(sampsiz)) {
    if (dat$y[i] == 1) {
      dat$X[i,] = c(mup + sigma.sqrt %*% rnorm(p0), rnorm(pp - p0))
    } else {
      dat$X[i,] = c(mum + sigma.sqrt %*% rnorm(p0), rnorm(pp - p0))
    }
  }
  # returns:
  dat
}

