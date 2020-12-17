### glmboost for Gaussian distribution with non-linear Base Learners ###
require("mboost")
require('plyr')
require('mvtnorm')

fsim <- function(seed, n, p, pinf, n.test, mu, sigma.y ,corr){
  source('DeselectBoost.R')
  
  selectedVar.risk = vector('list',9) 
  true.positive.risk = vector('list',9) 
  false.positive.risk = vector('list',9)
  MSEP <- vector('list',9)
  
  loss <- function(y,f){(y - f)^2}
  
  beta = c(1.5,1,-0.5,0.5,-1,-1.5,rep(0,p-6))
  mean_X = rep(mu,length = p)
  sigma_X = toeplitz(sapply(seq(0,p-1), function(x) corr^x))
  
  set.seed(seed) 
  
  X = rmvnorm(n, mean = mean_X, sigma = sigma_X)
  dat.train = data.frame(X, y = rnorm(n, beta[1]*sin(X[,1]) + beta[2]*(X[,2])^3 + beta[3] * (X[,3]^2) +  beta[4]*X[,4] +  beta[5]*X[,5] +  beta[6]*X[,6] , sigma.y))
  
  X.test = rmvnorm(n.test ,mean = mean_X, sigma = sigma_X)
  dat.test = data.frame(X.test, y = rnorm(n.test, beta[1]*sin(X.test[,1]) + beta[2]*(X.test[,2])^3 + beta[3] * (X.test[,3]^2) +  beta[4]*X.test[,4] +  beta[5]*X.test[,5] +  beta[6]*X.test[,6] , sigma.y))

  glm1 = gamboost(y ~ ., data = dat.train, control = boost_control(mstop = 500)) 
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
  MSEP[[3]] = mean(loss(dat.test$y,predict(glm1, newdata = dat.test,type = 'response')))
  
  #################################
  #---- DeselectBoost
  i = 3
  for(k in c(0.01,0.025,0.05,0.1,0.125)){
    i = i+1
    
    glm1_after <- DeselectBoost(glm1, data = dat.train, fam = Gaussian(), tau = k)
    
    selectedVar.risk[[i]] <- names(coef(glm1_after$model)[-1])
    true.positive.risk[[i]] <- length(which(true.var %in% selectedVar.risk[[i]]))
    false.positive.risk[[i]] <- length(which(false.var %in% selectedVar.risk[[i]]))
    
    MSEP[[i]] = mean(loss(dat.test$y, predict(glm1_after$model, newdata = dat.test, type = 'response')))
    
  }
  
  #################################
  #----- oSE 
  se_cv = sd(cvr[,stopIT])/sqrt(dim(cv25)[2])
  error_mstop <- mean(cvr[,stopIT])
  mean_cv <- apply(cvr[,1:stopIT], 2, mean)
  opt_stop_oSE =  which.min(mean_cv > se_cv + error_mstop)-1
  glm1_oSE = glm1[opt_stop_oSE]
  
  selectedVar.risk[[1]] <- names(coef(glm1_oSE)[-1])
  true.positive.risk[[1]] <- length(which(true.var %in% names(coef(glm1_oSE))))
  false.positive.risk[[1]] <- length(which(false.var %in% names(coef(glm1_oSE))))
  
  MSEP[[1]] <- mean(loss(dat.test$y, predict(glm1_oSE, newdata = dat.test,type = 'response')))
  
  ###################################
  #---- RobustC
  c_rC = 1.05
  error_mstop <- mean(cvr[,stopIT])
  mean_cv <- apply(cvr[,1:stopIT], 2, mean)
  opt_stop_rC = which.min(mean_cv > c_rC * error_mstop)-1
  glm1_robustC = glm1[opt_stop_rC]
  
  selectedVar.risk[[2]] <- names(coef(glm1_robustC)[-1])
  true.positive.risk[[2]] <- length(which(true.var %in% names(coef(glm1_robustC))))
  false.positive.risk[[2]] <- length(which(false.var %in% names(coef(glm1_robustC))))
  MSEP[[2]] <-  mean(loss(dat.test$y, predict(glm1_robustC,newdata=dat.test,type = "response")))
  
  return(list(Variables= selectedVar.risk,  MSEP = MSEP, seed = seed, true.positive.risk=true.positive.risk, false.positive.risk = false.positive.risk, mstop = stopIT, opt_stop_oSE = opt_stop_oSE, opt_stop_rC = opt_stop_rC, beta = beta, prozent = c(0.01,0.025,0.05,0.1,0.125)))
  
}  

n = 500
p = 20
pinf = 6
n.test = 1000
mu = 0
sigma.y = 1
corr = 0.2
numCores = 1

results = mclapply(1, fsim, mc.cores = numCores, n = n, p = p, pinf = pinf, n.test = n.test, mu = mu, sigma.y = sigma.y, corr =  corr, mc.set.seed = TRUE,  mc.preschedule = FALSE)

results$n = n 
results$n.test = n.test
results$p = p
results$sigma.y = sigma.y
results$mu = mu
results$corr = corr


save(results, file="DeselectionNonLinearRegression.RData")
