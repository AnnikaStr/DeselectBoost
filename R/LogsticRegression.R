library("mboost")
library("mvtnorm")
library('plyr')
library('pROC')

fsim <- function(seed, n, p, pinf, n.test, mu, sigma.y, corr, tau){
  source('DeselectBoost.R')
  
  selectedVar.risk = vector('list',9)
  true.positive.risk = vector('list',9) 
  false.positive.risk = vector('list',9)
  AUC <- vector('list',9)
  Brier <- vector('list',9)
  
  beta = c(-5,-2.5,-1,1,2.5,5,rep(0,p-6))
  mean_X = rep(mu,length = p)
  sigma_X = toeplitz(sapply(seq(0,p-1), function(x) corr^x))
  
  set.seed(seed) 
  
  X = rmvnorm(n, mean = mean_X, sigma = sigma_X)
  dat.train = data.frame(X, y = as.factor(rbinom(n, 1,  exp( X %*% beta) / ( 1 + exp(X %*% beta)))))

  X.test = rmvnorm(n.test ,mean = mean_X, sigma = sigma_X)
  dat.test = data.frame(X.test, y = as.factor(rbinom(n.test, 1,  exp( X.test %*% beta) / ( 1 + exp(X.test %*% beta)))))
  
  glm1 = glmboost(y ~., data = dat.train, family = Binomial())
  cv25 <- cv(model.weights(glm1),type = 'kfold')
  cvr = cvrisk(glm1, folds = cv25, grid = 1:10000, papply = mclapply)
  stopIT = mstop(cvr)
  glm1[mstop(cvr)]
  
  nameVar <- names(dat.train)[1:p] 
  true.var <- nameVar[1:pinf]
  false.var <- nameVar[(pinf+1):p]
  
  selectedVar.risk[[3]] <- names(coef(glm1)[-1])
  true.positive.risk[[3]] <- length(which(true.var %in% names(coef(glm1))))
  false.positive.risk[[3]] <- length(which(false.var %in% names(coef(glm1))))
  
  AUC[[3]] = roc(dat.test$y, as.numeric(predict(glm1, newdata = dat.test,type = "response")))$auc
  Brier[[3]] = mean((as.numeric(predict(glm1, newdata = dat.test,type = "response"))-(as.numeric(dat.test$y)-1))^2)

  #################################
  #---- DeselectBoost
  i = 3
  for(k in tau){
    i = i+1
    glm1_after <- DeselectBoost(glm1, data = dat.train, fam = Binomial(), tau = k)
    selectedVar.risk[[i]] <- names(coef(glm1_after$model)[-1])
    true.positive.risk[[i]] <- length(which(true.var %in% selectedVar.risk[[i]]))
    false.positive.risk[[i]] <- length(which(false.var %in% selectedVar.risk[[i]]))
    
    AUC[[i]] = roc(dat.test$y, as.numeric(predict(glm1_after$model, newdata = dat.test,type = "response")))$auc
    Brier[[i]] = mean((as.numeric(predict(glm1_after$model, newdata = dat.test,type = "response"))-(as.numeric(dat.test$y)-1))^2)
  }
  
  #################################
  #---- oSE 
  se_cv = sd(cvr[,stopIT])/sqrt(dim(cv25)[2])
  error_mstop <- mean(cvr[,stopIT])
  mean_cv <- apply(cvr[,1:stopIT], 2, mean)
  opt_stop_oSE =  which.min(mean_cv > se_cv + error_mstop)-1
  glm1_oSE = glm1[opt_stop_oSE]
  
  selectedVar.risk[[1]] <- names(coef(glm1_oSE)[-1])
  true.positive.risk[[1]] <- length(which(true.var %in% names(coef(glm1_oSE))))
  false.positive.risk[[1]] <- length(which(false.var %in% names(coef(glm1_oSE))))
  
  AUC[[1]] = roc(dat.test$y, as.numeric(predict(glm1_oSE, newdata = dat.test,type = "response")))$auc
  Brier[[1]] = mean((as.numeric(predict(glm1_oSE, newdata = dat.test,type = "response"))-(as.numeric(dat.test$y)-1))^2)

  #################################
  #---- RobustC
  c_rC = 1.1
  error_mstop <- mean(cvr[,stopIT])
  mean_cv <- apply(cvr[,1:stopIT], 2, mean)
  opt_stop_rC = which.min(mean_cv > c_rC * error_mstop)-1
  glm1_robustC = glm1[opt_stop_rC]
  
  selectedVar.risk[[2]] <- names(coef(glm1_robustC)[-1])
  true.positive.risk[[2]] <- length(which(true.var %in% names(coef(glm1_robustC))))
  false.positive.risk[[2]] <- length(which(false.var %in% names(coef(glm1_robustC))))
  
  AUC[[2]] = roc(dat.test$y, as.numeric(predict(glm1_robustC, newdata = dat.test,type = "response")))$auc
  Brier[[2]] = mean((as.numeric(predict(glm1_robustC, newdata = dat.test,type = "response"))-(as.numeric(dat.test$y)-1))^2)
  
  return(list(Variables= selectedVar.risk, AUC = AUC, Brier = Brier, seed = seed, true.positive.risk = true.positive.risk, false.positive.risk = false.positive.risk, mstop = stopIT,opt_stop_oSE = opt_stop_oSE, opt_stop_rC = opt_stop_rC, beta = beta))
  
}  

p = 20
pinf = 6
n = 500
n.test = 1000
mu = 0
sigma.y = 1
corr = 0.8
tau = c(0.01,0.025,0.05,0.075,0.1,0.125)
numCores = 1

results = mclapply(1, fsim, mc.cores = numCores, n = n, p = p, pinf = pinf, n.test = n.test, mu = mu, sigma.y = sigma.y, corr =  corr, tau = tau, mc.set.seed = TRUE,  mc.preschedule = FALSE)

results$n = n 
results$n.test = n.test
results$p = p
results$sigma.y = sigma.y
results$mu = mu
results$corr = corr
results$tau = tau

save(results, file="DeselectionLogisticRegression.RData")




