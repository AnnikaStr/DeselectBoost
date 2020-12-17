library('gamboostLSS')
library("mvtnorm")
library('plyr')
library("MASS")


fsim <- function(seed, p, n , n.test, corr, Mu, tau){ 
  source('DeselectBoost.R')
  
  selectedVar.mu.risk = vector('list',8) 
  selectedVar.sigma.risk = vector('list',8)
  true.positive.mu.risk = vector('list',8) 
  false.positive.mu.risk = vector('list',8)
  true.positive.sigma.risk = vector('list',8) 
  false.positive.sigma.risk = vector('list',8)
  Likelihood <- vector('list',8)

  loss <- function(y,mu,f){ -sum(dnorm( y,  mu, exp(f), log = TRUE)) }
  
  beta_mu = c(- 2, 1.25 , 1, rep(0,p-3)) 
  beta_sigma = c(0, 0, 0 , 0.5, -0.5, .5,rep(0,p-6))
  
  mu = rep(Mu,p)
  sigma = toeplitz(sapply(seq(0,p-1), function(x) corr^x))
  
  set.seed(seed)
  X = rmvnorm(n , mu, sigma)
  dat.train = data.frame(X, y = rnorm(n, mean = X %*% beta_mu, sd = exp(-1 +X %*% beta_sigma)))
  
  X.test = mvrnorm(n.test , mu, sigma)
  dat.test = data.frame(X.test, y = rnorm(n.test, mean = X.test %*% beta_mu, sd = exp(-1 + X.test %*% beta_sigma)))
  
  LSS = glmboostLSS(list(mu = y ~ .,sigma = y ~ .), data = dat.train, families = GaussianLSS(), method = 'noncyclic',control = boost_control(mstop = 2000,nu=0.005))
  cv25 <- cv(model.weights(LSS),type = 'kfold')
  cvr = cvrisk(LSS, folds = cv25, grid = 1:2000, papply = mclapply, trace = F)  
  stopIT = mstop(cvr)
  
  if(stopIT == 9000){
    cvr = cvrisk(LSS, folds = cv25, grid = 1:20000, papply = mclapply, trace = F) 
  }
  
  stopIT = mstop(cvr)
  LSS[mstop(cvr)] 
  LSS = glmboostLSS(list(mu = y ~ .,sigma = y ~ .), data = dat.train, families = GaussianLSS(), method = 'noncyclic',control = boost_control(mstop = stopIT,nu=0.005)) # so that combined_risk has the correct length
  
  
  nameVar <- names(dat.train)[1:p] 
  true.var.mu = nameVar[beta_mu!=0]
  false.var.mu = nameVar[beta_mu==0]
  true.var.sigma = nameVar[beta_sigma!=0]
  false.var.sigma = nameVar[beta_sigma==0]
  
  selectedVar.mu.risk[[3]] <- names(coef(LSS$mu)[-1])
  selectedVar.sigma.risk[[3]] <- names(coef(LSS$sigma)[-1])
  true.positive.mu.risk[[3]] <- length(which(true.var.mu %in% names(coef(LSS$mu))))
  false.positive.mu.risk[[3]] <- length(which(false.var.mu %in% names(coef(LSS$mu))))
  true.positive.sigma.risk[[3]] <- length(which(true.var.sigma %in% names(coef(LSS$sigma))))
  false.positive.sigma.risk[[3]] <- length(which(false.var.sigma %in% names(coef(LSS$sigma))))
  
  Likelihood[[3]] = loss(y = dat.test$y,mu = predict(LSS$mu,newdata=dat.test),f = predict(LSS$sigma, type = "response", newdata = dat.test))
  
  #----- oSE 
  se_cv = sd(cvr[,stopIT])/sqrt(dim(cv25)[2])
  error_mstop <- mean(cvr[,stopIT])
  mean_cv <- apply(cvr[,1:stopIT], 2, mean)
  opt_stop_oSE =  which.min(mean_cv > se_cv + error_mstop) - 1
  LSS_oSE = LSS[opt_stop_oSE]

  selectedVar.mu.risk[[1]] <- names(coef(LSS_oSE$mu)[-1])
  selectedVar.sigma.risk[[1]] <- names(coef(LSS_oSE$sigma)[-1])
  true.positive.mu.risk[[1]] <- length(which(true.var.mu %in% names(coef(LSS_oSE$mu))))
  false.positive.mu.risk[[1]] <- length(which(false.var.mu %in% names(coef(LSS_oSE$mu))))
  true.positive.sigma.risk[[1]] <- length(which(true.var.sigma %in% names(coef(LSS_oSE$sigma))))
  false.positive.sigma.risk[[1]] <- length(which(false.var.sigma %in% names(coef(LSS_oSE$sigma))))
  
  Likelihood[[1]] <- loss(y = dat.test$y,mu = predict(LSS_oSE$mu,newdata=dat.test),f = predict(LSS_oSE$sigma, type = "response", newdata = dat.test))
  
  #---- RobustC
  c_rC = 1.05
  error_mstop <- mean(cvr[,stopIT])
  mean_cv <- apply(cvr[,1:stopIT], 2, mean)
  opt_stop_rC = which.min(mean_cv > c_rC * error_mstop)-1
  LSS_robustC = LSS[opt_stop_rC]
  
  selectedVar.mu.risk[[2]] <- names(coef(LSS_robustC$mu)[-1])
  selectedVar.sigma.risk[[2]] <- names(coef(LSS_robustC$sigma)[-1])
  true.positive.mu.risk[[2]] <- length(which(true.var.mu %in% names(coef(LSS_robustC$mu))))
  false.positive.mu.risk[[2]] <- length(which(false.var.mu %in% names(coef(LSS_robustC$mu))))
  true.positive.sigma.risk[[2]] <- length(which(true.var.sigma %in% names(coef(LSS_robustC$sigma))))
  false.positive.sigma.risk[[2]] <- length(which(false.var.sigma %in% names(coef(LSS_robustC$sigma))))
  
  Likelihood[[2]] <-  loss(y = dat.test$y,mu = predict(LSS_robustC$mu,newdata=dat.test),f = predict(LSS_robustC$sigma, type = "response", newdata = dat.test))
  
  #---- DeselectBoost
  i = 3
  for(k in tau){
    i = i+1
    
    LSS_after <- DeselectBoost(LSS, fam = GaussianLSS(), data = dat.train, tau = k)
    
    selectedVar.mu.risk[[i]] <-  names(coef(LSS_after$model$mu)[-1])
    true.positive.mu.risk[[i]] <- length(which(true.var.mu %in% selectedVar.mu.risk[[i]]))
    false.positive.mu.risk[[i]] <- length(which(false.var.mu %in% selectedVar.mu.risk[[i]]))
    
    selectedVar.sigma.risk[[i]] <-  names(coef(LSS_after$model$sigma)[-1])
    true.positive.sigma.risk[[i]] <- length(which(true.var.sigma %in% selectedVar.sigma.risk[[i]]))
    false.positive.sigma.risk[[i]] <- length(which(false.var.sigma %in% selectedVar.sigma.risk[[i]]))
    
    Likelihood[[i]] <- loss(y = dat.test$y, mu = predict(LSS_after$model$mu,newdata=dat.test, type = 'response'), f = predict(LSS_after$model$sigma, type = "response", newdata = dat.test))
  }
  
  return(list(mu = selectedVar.mu.risk, sigma = selectedVar.sigma.risk, Likelihood = Likelihood, seed = seed, true.positive.mu.risk = true.positive.mu.risk, false.positive.mu.risk = false.positive.mu.risk,  true.positive.sigma.risk = true.positive.sigma.risk, false.positive.sigma.risk = false.positive.sigma.risk , mstop = stopIT, opt_stop_oSE = opt_stop_oSE, opt_stop_rC = opt_stop_rC, beta_mu = beta_mu, beta_sigma = beta_sigma))
}


p = 20
n = 500
n.test = 1000
corr = 0.2
Mu = 0
tau = c(0.01,0.025,0.05,0.1,0.125)
numCores = 1

results = mclapply(1, fsim, mc.cores = numCores, n = n, p = p, n.test = n.test, Mu = Mu, corr =  corr, tau = tau, mc.set.seed = TRUE,  mc.preschedule = FALSE)

results$n = n 
results$n.test = n.test
results$p = p
results$Mu = Mu
results$corr = corr
results$tau = tau

save(results, file="deselectionGamboostLSS_low_low.RData")



