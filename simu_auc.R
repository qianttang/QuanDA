library(caret)
library(glmnet)
library(pROC)
library(Matrix)
library(randomForest)
source("~/research/ImbalanceData/hdqr/nips2025/functions_old.R")
source("~/research/ImbalanceData/hdqr/nips2025/class_functions.R")
library(splines)
library(PRROC)


l_list <- 10^(seq(1,-4, length.out=30))
n = 200
pp = 10000
rho = 2
sina = 5
blanc = 0.9
tau=seq(max(0, blanc-0.05), max(1.0, blanc+0.05), 0.01)
lam2 = 1e-02
auc_res <- matrix(NA, 50, 5)
for(i in 1:50){
  set.seed(i)
  train_data = DataGen(n, pp=pp, rho=rho, blanc=blanc, sina=sina)
  tmp1 = stan(train_data)
  tmp1 <- tmp1$train
  train = list(x=tmp1$X, y=tmp1$y)
  set.seed(2024+i)
  test_data = DataGen(n, pp=pp, rho=rho, blanc=blanc, sina=sina)
  tmp2 = stan(test_data)
  tmp2 <- tmp2$train
  test = list(x=tmp2$X, y=tmp2$y)

  ans_rep <- hdqr_z_rep(train, test, q=0.0, tau, lam2, l_list, rep=10, i)
  auc_res[i,1] <- ans_rep$res[4]

  wts2 <- length(train$y) / table(train$y)
  wts <- ifelse(train$y==0, wts2[1], wts2[2])
  res2 <- logistic_class(train, test, wts, l_list, lam2, i)
  auc_res[i,2] <- res2$res[4]

  res3 <- dsda_class(train, test, l_list, lam2=lam2, i)
  auc_res[i,3] <- res3$res[4]
  
  res5 <- rf_class(train, test, wts,i)
  auc_res[i,4] <- res5$res[4]

  res6 <-  tryCatch(smote_class(train, test, i)
  auc_res[i,5] <- res6$res[4]

}