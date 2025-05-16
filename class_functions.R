library(hdqr)


cv_z <- function(x, y, z, tau, lambda = NULL, nfolds=5L,
                       foldid, delta=.125, lam2=0.01, q=0.5, ...){
  ############################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ###Now fit the nfold models and store them
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = x.row))
  } else nfolds <- max(foldid)
  if (nfolds < 3L)
    stop("nfolds must be at least 3; nfolds = 5 recommended")
  lambda <- sort(lambda, decreasing=TRUE)
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- hdqr(x=x[!which, , drop=FALSE], y=z[!which],
      lambda=lambda, tau=tau, ...)
    if (outlist[[i]]$jerr != 0)
      stop(paste("Error occurs when fitting the", i, "th folder."))
  }
  cvstuff <- f1.path_z(outlist, x, y, tau, lambda, foldid, x.row, delta, lam2, q, ...)
  cvf1 <- cvstuff$cvf1
  ## wrap up output
  cvmax <- max(cvf1, na.rm=TRUE)
  idmin <- cvf1 >= cvmax
  lambda.max <- max(lambda[idmin], na.rm=TRUE)
  out <- list(lambda=lambda, cvf1=cvf1, lambda.min=lambda.max)
  obj <- c(out)
  # class(obj) <- "cv.kqr"
  obj
}

cv <- function(x, y, tau, lambda = NULL, nfolds=5L,
                       foldid, delta=.125, lam2=0.01, q=0.5, ...){
  ############################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ###Now fit the nfold models and store them
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = x.row))
  } else nfolds <- max(foldid)
  if (nfolds < 3L)
    stop("nfolds must be at least 3; nfolds = 5 recommended")
  lambda <- sort(lambda, decreasing=TRUE)
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- hdqr(x=x[!which, , drop=FALSE], y=y[!which],
      lambda=lambda, tau=tau, ...)
    if (outlist[[i]]$jerr != 0)
      stop(paste("Error occurs when fitting the", i, "th folder."))
  }
  cvstuff <- f1.path(outlist, x, y, tau, lambda, foldid, x.row, delta, lam2, q, ...)
  cvf1 <- cvstuff$cvf1
  ## wrap up output
  cvmax <- max(cvf1, na.rm=TRUE)
  idmin <- cvf1 >= cvmax
  lambda.max <- max(lambda[idmin], na.rm=TRUE)
  out <- list(lambda=lambda, cvf1=cvf1, lambda.min=lambda.max)
  obj <- c(out)
  # class(obj) <- "cv.kqr"
  obj
}

f1.path_z <- function(outlist, x, y, tau, lambda, foldid, x.row,
                       delta, lam2, q, ...){
  nfolds <- max(foldid)
  predmat <- matrix(NA, x.row, length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[whichfold, , drop = FALSE]) -1
    nlami <- length(fitobj$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  # predmat <- apply(predmat, 2, phi_x, q=q)
  cvf1 <- apply(predmat, 2, f1, y=y)
  out <- list(cvf1=cvf1)
  out
}


f1.path <- function(outlist, x, y, tau, lambda, foldid, x.row,
                       delta, lam2, q, ...){
  nfolds <- max(foldid)
  predmat <- matrix(NA, x.row, length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[whichfold, , drop = FALSE]) 
    nlami <- length(fitobj$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  # predmat <- apply(predmat, 2, phi_x, q=q)
  cvf1 <- apply(predmat, 2, f1, y=y)
  out <- list(cvf1=cvf1)
  out
}


f1 <- function(pred, y){
  # pred <- matrix(pred, length(pred), 1)
  # pred <- apply(pred, 1, phi_x, q=q)
  # print(as.factor(pred))
  # print(as.factor(y))
  # res <- confusionMatrix(as.factor(pred), as.factor(y), positive="1")
  # f1 <- ifelse(is.na(res$byClass[7]), 0, res$byClass[7])
  rocc <- roc(as.factor(y), as.vector(pred))$auc
  return(rocc)
}

phi_x <- function(z,qval){
  ifelse(z>=qval, 1, 0)
}

cv.folds<-function (n, folds = 10) 
{
    split(sample(1:n), rep(1:folds, length = n))
}

 cv_dsda <- function (x, y, K = 5, standardize=FALSE, lambda=lambda,lam2, thresh=1e-7,
  lambda.opt="min", alpha=1, ...)
 {
     n<-length(y)
     n1<-n-sum(y)
     n2<-sum(y)

     if(K>n)stop("The number of folds should be smaller than the sample size.")

     all.folds <- cv.folds(length(y), K)
     
     if(missing(lambda)){
        fit <- glmnet(x, y,  family="gaussian",alpha=alpha,standardize=standardize,thresh=thresh)
        lambda<-fit$lambda
      }

     nlambda<-length(lambda)
     predmat <- matrix(0, n, nlambda)
     for(l in 1:nlambda){
        lam <- lambda[l] 
        for (i in seq(K)) {
          omit <- all.folds[[i]]
          fit <- dsda.path(x[-omit, , drop = FALSE], y[-omit], lambda=lam+lam2,
            standardize=standardize,
            thresh=thresh, alpha=lam/(lam+lam2), ...)
         preds <- predict.dsda(fit,x[omit,,drop=FALSE], type="response")
         nlami <- length(fit$lambda)
         predmat[omit, l] <- preds
        }
      }
    
     cvf1 <- apply(predmat, 2, f1, y=y)
     # residmat[is.na(residmat)]<-min(n1/n,n2/n)
     # residmat<-matrix(residmat,nrow=nlambda)
     # cv <- apply(residmat, 1, mean)
     # cv.error <- sqrt(apply(residmat, 1, var)/K) cv.error = cv.error,
     if(lambda.opt=="min"){
        bestlambda<-max(lambda[which(cvf1==max(cvf1))])}
     else{
        bestlambda<-max(lambda[which(cvf1==max(cvf1))])}
     object <- list(lambda = lambda, cv = cv, bestlambda=bestlambda)
     invisible(object)
 }


cv_log <- function(x, y, lambda, nfolds=5L,
                       foldid, alpha=1, lam2, ...){
  ############################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ###Now fit the nfold models and store them
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = x.row))
  } else nfolds <- max(foldid)
  if (nfolds < 3L)
    stop("nfolds must be at least 3; nfolds = 5 recommended")
  lambda <- sort(lambda, decreasing=TRUE)
  outlist <- as.list(seq(nfolds))
  nlambda<-length(lambda)
  predmat <- matrix(0, x.row, nlambda)
  wts2 <- x.row / table(train$y) 
  # wts <- ifelse(train$y==0, wts2[1], wts2[2])
  for(l in 1:nlambda){
    lam <- lambda[l] 
    for (i in seq(nfolds)) {
      which <- foldid == i
      wts <- ifelse(y[!which] ==0, wts2[1], wts2[2])
      fit <- glmnet(x[!which, , drop=FALSE], y[!which], family="binomial", 
        weights=wts, alpha=lam/(lam+lam2), lambda=lam+lam2, ...)
      preds <- predict(fit, x[which,,drop=FALSE], type="response")
      predmat[which, l] <- preds
    }
  }
  predmat <- matrix(as.numeric(predmat), x.row, nlambda)
  cvf1 <- apply(predmat, 2, f1, y=y)
  lambda.max <- max(lambda[which(cvf1==max(cvf1))])
  out <- list(lambda=lambda, cvf1=cvf1, lambda.min=lambda.max)
  obj <- c(out)
  obj
}














