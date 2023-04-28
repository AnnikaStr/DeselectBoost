DeselectBoost <- function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  require('plyr')
  tau = ifelse(is.null(tau), 0.01, tau)

  if(any(class(object) %in% 'mboost')){
    DeselectBoost_1(object, data = data, fam = fam, tau = tau, method = method[1])
  }else{
    switch(length(names(object))-1,{
     DeselectBoostLSS_2(object, data = data, fam = fam, tau = tau, method = method[1])},{
     DeselectBoostLSS_3(object, data = data, fam = fam, tau = tau, method = method[1])},{
     DeselectBoostLSS_4(object, data = data, fam = fam, tau = tau, method = method[1])  
    })
  }
}  


DeselectBoost_1 <- function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  require('plyr')
  
  if(is.null(data) && class(object$model.frame()) == 'list'){return(stop("Please enter the data."))
  } else if(!is.null(data)){
    data = data
  }else{
    data <- object$model.frame()
  }

  nameVar <- names(coef(object,which = ''))[-1]
  which.response <- which(sapply(1:(dim(data)[2]), function(x){identical(as.numeric(data[,x]), as.numeric(object$response))}))
  name.response <- colnames(data)[which.response]
  
  mstop <- object$mstop()
  RiskRed <- object$risk()
  totalRiskRed <- RiskRed[1] - RiskRed[mstop+1] 
  diffRiskRed = sapply(seq(1:(mstop)), function(k){RiskRed[k]-RiskRed[k+1]})
  
  if(any(class(object) %in% "glmboost")){
    select = selected(object) - 1
    diffRiskRed = diffRiskRed[selected(object)-1 != 0]
  }else{
    select = selected(object)
  }
  
  select = select[select != 0]
  Var = count(select)[[1]]
  Risk.Var <- lapply(1:length(Var),function(j){sum(diffRiskRed[which(count(select)[[1]][j] == select)])})
  
  n.parameter <- c(names(object$coef()))
  if('(Intercept)' %in% n.parameter) n.parameter <- n.parameter[-which(n.parameter == '(Intercept)')]
  
  
  Risk.order <- data.frame(Var,n.parameter, as.numeric(Risk.Var))
  Risk.order <- Risk.order[order(Risk.order$as.numeric.Risk.Var.),]
  Risk.order$CumRisk <- cumsum(Risk.order$as.numeric.Risk.Var.)
  colnames(Risk.order) <- c( 'Var', 'VarName', 'Risk', 'CumRisk')
  
  perc <- ifelse(is.null(tau), 0.01, tau) 
  percRiskRed <- totalRiskRed * perc
  if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
  }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(empty(RiskRedOver)){form2 = as.formula(paste(name.response, "~ 1"))
  }else{form2 <-as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName, collapse= "+")))
  if(!is.null(environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]])){
    dfbase = environment(environment(object[["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["df"]]
  }}
  #if(is.null(object$call$family)){ fam <-  Gaussian()
  #}else fam <- eval(parse(text = object$call$family))
  
  if(any(class(object) %in% "glmboost")){
    model_after = glmboost(form2, data = data, weights = model.weights(object), family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
  }else{
    model_after = gamboost(form2, data = data, weights = model.weights(object), family = fam, control = boost_control(mstop = mstop, nu = object$control$nu, risk = object$control$risk))
  }
  
  out <- model_after
  out$tau  = tau
  out$deselectmethod = method[1] 
  class(out) <- c(class(out))
  
  return(out)
}

DeselectBoostLSS_2 <- function(object, data = NULL , fam, tau = NULL, method = c('attributable','cumulative')){
  require('plyr')

  data = attr(object,'data')
  mstop <- ifelse(any(class(object)%in%'nc_mboostLSS'), mstop(object)[1],sum(mstop(object)))
  
  if(any(class(object) %in% "betaboost")){
    
    parameter = names(object)
    which.response <- which(sapply(1:dim(data)[2], function(x){identical(as.numeric(data[,x]), as.numeric(object$mu$response))}))
    name.response <- colnames(data)[which.response]
    select <- selected(object, parameter = names(object))
    
    RiskRed <- risk(object)
    totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
    diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
    
    Var = lapply(1:length(parameter), function(i){count(select[[i]])[2]})
    Risk.Var <- lapply(1:length(parameter),function(i){sapply(1:dim(Var[[i]])[1],function(j){sum(diffRiskRed[[i]][which(count(select[[i]])[[1]][j] == select[[i]])])})})
    
    w.parameter <- c(rep(parameter[1],length(count(select[[1]])[[1]])),rep(parameter[2],length(count(select[[2]])[[1]])))
    n.parameter <- c(names(object[[1]]$coef()),names(object[[2]]$coef()))
    
    Risk.order <- data.frame(w.parameter, n.parameter, unlist(Risk.Var))
    colnames(Risk.order) <- c('parameter', 'VarName', 'Risk')
    Risk.order <- Risk.order[order(Risk.order$Risk),]
    Risk.order$CumRisk <- cumsum(Risk.order$Risk)
    
    perc <- ifelse(is.null(tau), 0.01, tau) 
    percRiskRed <- (totalRiskRed[[1]]+totalRiskRed[[2]]) * perc
    
    if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
    }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
    
    if(any(RiskRedOver$parameter == parameter[1]) && any(RiskRedOver$parameter == parameter[2])){  
      form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
      form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
    }else{
      if(any(RiskRedOver$parameter == parameter[1])){ 
        form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
        form2 <- as.formula(name.response~1)
      }else{
        form1 <- as.formula(name.response~1)
        form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
      }
    }
    
    model_after <- betaboost(form1, form2, data=data, 
                             iterations = mstop , method = "noncycl", sl=c(object[[1]]$control$nu,object[[2]]$control$nu), weights = model.weights(object), form.type="gamboost")  
    
  }else{
    if (length(risk(object, merge = T)) > (mstop + 2)) return(stop("risk cannot contain more entries than mstop"))
    
    parameter = names(object)
    which.response <- which(sapply(1:dim(data)[2], function(x){identical(as.numeric(data[,x]), as.numeric(object$mu$response))}))
    name.response <- colnames(data)[which.response]
    
    select <- selected(object, parameter = names(object))
    select1 <- selected(object, parameter = names(object))
    
    if(any(class(object) %in% 'glmboostLSS')){
      select[[1]] <- select[[1]]-1  
      select[[2]] <- select[[2]]-1  
      select[[1]] <- select[[1]][select[[1]] !=0]
      select[[2]] <- select[[2]][select[[2]] !=0]
      
      RiskRed <- risk(object)
      totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
      diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
      diffRiskRed[[1]] <- diffRiskRed[[1]][select1[[1]]-1 != 0]
      diffRiskRed[[2]] <- diffRiskRed[[2]][select1[[2]]-1 != 0]
      
    }else{
      RiskRed <- risk(object)
      totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
      diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
      
    }
    
    Var = lapply(1:length(parameter), function(i){count(select[[i]])[2]})
    zero <- sapply(1:length(parameter),function(i){length(select[[i]])})
    Risk.Var <- lapply(which(zero != 0),function(i){sapply(1:dim(Var[[i]])[1],function(j){sum(diffRiskRed[[i]][which(count(select[[i]])[[1]][j] == select[[i]])])})})
    
    w.parameter <- c(rep(parameter[1],length(count(select[[1]])[[1]])),rep(parameter[2],length(count(select[[2]])[[1]])))
    n.parameter <- c(names(object[[1]]$coef()),names(object[[2]]$coef()))
    if('(Intercept)' %in% n.parameter) n.parameter <- n.parameter[-which(n.parameter == '(Intercept)')]
    
    Risk.order <- data.frame(w.parameter, n.parameter, unlist(Risk.Var))
    colnames(Risk.order) <- c('parameter', 'VarName', 'Risk')
    Risk.order <- Risk.order[order(Risk.order$Risk),]
    Risk.order$CumRisk <- cumsum(Risk.order$Risk)
    
    perc <- ifelse(is.null(tau), 0.01, tau) # 0.01 is default value
    percRiskRed <- (totalRiskRed[[1]]+totalRiskRed[[2]]) * perc
    if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
    }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
    
    if(is.na(parameter[2])) parameter[2]<-0
    if(is.na(parameter[1])) parameter[1]<-0
    
    if(any(RiskRedOver$parameter == parameter[1]) && any(RiskRedOver$parameter == parameter[2])){  
      form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+"))) 
      form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
    }else{
      if(any(RiskRedOver$parameter == parameter[1])){ 
        form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
        form2 <- as.formula(paste(name.response, '~1'))
      }else{
        form1 <-  as.formula(paste(name.response, '~1'))
        form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
      }
    }
    formula <- list(form1,form2)
    names(formula)<- names(object)
    dfbase <- environment(environment(environment(object[[1]][["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["ret"]][["model.frame"]])[["df"]]
    
    
    if(any(class(object) %in% 'nc_mboostLSS')){
      if(any(class(object) %in% "glmboostLSS")){
        model_after = glmboostLSS(formula, data = data, families = fam,  method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop =  mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu))))
      }else{
        model_after = gamboostLSS(formula,  data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop =  mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu))))
      }
    }else{
      if(any(class(object) %in% "glmboostLSS")){
        model_after = glmboostLSS(formula, data = data, families = fam,  method = 'cyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu))))
      }else{
        model_after = gamboostLSS(formula, data = data, families = fam, method = 'cyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu))))
      }
    }
  }
  out <- model_after
  deselect_para <-  list(tau = perc, deselectmethod = method[1])
  out <-append(x = out, values = deselect_para)
  
  class(out) <- c(class(out))
  
  return(out)
}


DeselectBoostLSS_3 <- function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
 
  data = attr(object,'data')
  mstop <- ifelse(any(class(object)%in%'nc_mboostLSS'), mstop(object)[1],sum(mstop(object)))
  
  if (length(risk(object, merge = T)) > (mstop + 3)) stop("risk cannot contain more entries than mstop")
  
  parameter = names(object)
  which.response <- which(sapply(1:dim(data)[2], function(x){identical(as.numeric(data[,x]), as.numeric(object$mu$response))}))
  name.response <- colnames(data)[which.response]
 
  select <- selected(object, parameter = names(object))
  select1 <- selected(object, parameter = names(object))
  
  
  if(any(class(object) %in% 'glmboostLSS')){
    select[[1]] <- select[[1]]-1
    select[[2]] <- select[[2]]-1  
    select[[3]] <- select[[3]]-1
    select[[1]] <- select[[1]][select[[1]] !=0]
    select[[2]] <- select[[2]][select[[2]] !=0]
    select[[3]] <- select[[3]][select[[3]] !=0]
    
    RiskRed <- risk(object)
    totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
    diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
    diffRiskRed[[1]] <- diffRiskRed[[1]][select1[[1]]-1 != 0]
    diffRiskRed[[2]] <- diffRiskRed[[2]][select1[[2]]-1 != 0]
    diffRiskRed[[3]] <- diffRiskRed[[3]][select1[[3]]-1 != 0]
    
  }else{
    RiskRed <- risk(object)
    totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
    diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
    
  }
  
  Var = lapply(1:length(parameter), function(i){count(select[[i]])[2]})
  
  # If a parameter is not selected
  zero <- sapply(1:length(parameter),function(i){length(select[[i]])})
  Risk.Var <- lapply(which(zero != 0),function(i){sapply(1:dim(Var[[i]])[1],function(j){sum(diffRiskRed[[i]][which(count(select[[i]])[[1]][j] == select[[i]])])})})
  
  w.parameter <- c(rep(parameter[1],length(count(select[[1]])[[1]])),rep(parameter[2],length(count(select[[2]])[[1]])),rep(parameter[3],length(count(select[[3]])[[1]])))
  n.parameter <- c(names(object[[1]]$coef()),names(object[[2]]$coef()),names(object[[3]]$coef()))
  if('(Intercept)' %in% n.parameter)  n.parameter<-n.parameter[-which(n.parameter == '(Intercept)')]
  
  Risk.order <- data.frame(w.parameter, n.parameter, unlist(Risk.Var))
  colnames(Risk.order) <- c('parameter', 'VarName', 'Risk')
  Risk.order <- Risk.order[order(Risk.order$Risk),]
  Risk.order$CumRisk <- cumsum(Risk.order$Risk)
  
  perc <- ifelse(is.null(tau), 0.01, tau) 
  percRiskRed <- (totalRiskRed[[1]]+totalRiskRed[[2]]) * perc
  if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
  }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(is.na(parameter[2])) parameter[2]<-0
  if(is.na(parameter[1])) parameter[1]<-0
  if(is.na(parameter[3])) parameter[3]<-0
  
  if(any(RiskRedOver$parameter == parameter[1])){
    form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+"))) 
  }else{  form1 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[2])){
    form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
  }else{  form2 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[3])){
    form3 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[3]], collapse= "+")))
  }else{  form3 <-  as.formula(paste(name.response, '~1'))}
  
  formula <- list(form1,form2,form3) 
  names(formula)<- names(object)
  dfbase <- environment(environment(environment(object[[1]][["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["ret"]][["model.frame"]])[["df"]]
  
  
  if(any(class(object) %in% 'nc_mboostLSS')){
    if(any(class(object) %in% "glmboostLSS")){
      model_after = glmboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu))))
    }else{
      model_after = gamboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu))))}
  }else{
    if(any(class(object) %in% "glmboostLSS")){
      model_after = glmboostLSS(formula, data = data, families = fam, method = 'cyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu))))
    }else{
      model_after = gamboostLSS(formula, data = data, families = fam, method = 'cyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu))))
    }
  }
  
  out <- model_after
  deselect_para <-  list(tau = perc, deselectmethod = method[1])
  out <-append(x = out, values = deselect_para)
  
  class(out) <- c(class(out))
  
  return(out)
}


DeselectBoostLSS_4 <- function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
  data = attr(object,'data')
  
  mstop <- ifelse(any(class(object)%in%'nc_mboostLSS'), mstop(object)[1],sum(mstop(object)))
  
  if (length(risk(object, merge = T)) > (mstop + 4)) stop("risk cannot contain more entries than mstop")
  
  parameter = names(object)
  which.response <- which(sapply(1:dim(data)[2], function(x){identical(as.numeric(data[,x]), as.numeric(object$mu$response))}))
  name.response <- colnames(data)[which.response]

  select <- selected(object, parameter = names(object))
  select1 <- selected(object, parameter = names(object))
  if(any(class(object) %in% 'glmboostLSS')){
    select[[1]] <- select[[1]]-1  
    select[[2]] <- select[[2]]-1 
    select[[3]] <- select[[3]]-1
    select[[4]] <- select[[4]]-1
    select[[1]] <- select[[1]][select[[1]] !=0]
    select[[2]] <- select[[2]][select[[2]] !=0]
    select[[3]] <- select[[3]][select[[3]] !=0]
    select[[4]] <- select[[4]][select[[4]] !=0]
    
    RiskRed <- risk(object)
    totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
    diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
    diffRiskRed[[1]] <- diffRiskRed[[1]][select1[[1]]-1 != 0]
    diffRiskRed[[2]] <- diffRiskRed[[2]][select1[[2]]-1 != 0]
    diffRiskRed[[3]] <- diffRiskRed[[3]][select1[[3]]-1 != 0]
    diffRiskRed[[4]] <- diffRiskRed[[4]][select1[[4]]-1 != 0]
    
  }else{
    RiskRed <- risk(object)
    totalRiskRed <- lapply(1:length(parameter),function(i){RiskRed[[i]][1] - RiskRed[[i]][length(RiskRed[[i]])]})
    diffRiskRed <- lapply(1:length(parameter),function(j){sapply(seq(1:(length(RiskRed[[j]])-1)),function(k){RiskRed[[j]][k]-RiskRed[[j]][k+1]})})
  }
  
  Var = lapply(1:length(parameter), function(i){count(select[[i]])[2]})
  
  # If a parameter is not selected
  zero<-sapply(1:length(parameter),function(i){length(select[[i]])})
  Risk.Var <- lapply(which(zero != 0),function(i){sapply(1:dim(Var[[i]])[1],function(j){sum(diffRiskRed[[i]][which(count(select[[i]])[[1]][j] == select[[i]])])})})
  
  w.parameter <- c(rep(parameter[1],length(count(select[[1]])[[1]])),rep(parameter[2],length(count(select[[2]])[[1]])),rep(parameter[3],length(count(select[[3]])[[1]])),rep(parameter[4],length(count(select[[4]])[[1]])))
  n.parameter <- c(names(object[[1]]$coef()),names(object[[2]]$coef()),names(object[[3]]$coef()),names(object[[3]]$coef()))
  if('(Intercept)' %in% n.parameter)  n.parameter<-n.parameter[-which(n.parameter == '(Intercept)')]
  
  Risk.order <- data.frame(w.parameter, n.parameter, unlist(Risk.Var))
  colnames(Risk.order) <- c('parameter', 'VarName', 'Risk')
  Risk.order <- Risk.order[order(Risk.order$Risk),]
  Risk.order$CumRisk <- cumsum(Risk.order$Risk)
  
  perc <- ifelse(is.null(tau), 0.01, tau) 
  percRiskRed <- (totalRiskRed[[1]]+totalRiskRed[[2]]) * perc
  if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
  }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(is.na(parameter[2])) parameter[2]<-0
  if(is.na(parameter[1])) parameter[1]<-0
  if(is.na(parameter[3])) parameter[3]<-0
  if(is.na(parameter[4])) parameter[4]<-0
  
  if(any(RiskRedOver$parameter == parameter[1])){
    form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
  }else{  form1 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[2])){
    form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
  }else{  form2 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[3])){
    form3 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[3]], collapse= "+"))) 
  }else{  form3 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[4])){
    form4 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[4]], collapse= "+"))) 
  }else{  form4 <-  as.formula(paste(name.response, '~1'))}
  
  formula <- list(form1,form2,form3,form4) 
  names(formula)<- names(object)
  dfbase <- environment(environment(environment(object[[1]][["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["ret"]][["model.frame"]])[["df"]]
  
  
  if(any(class(object) %in% "glmboostLSS")){
    
    model_after = glmboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu))))
    
  }else{
    model_after = gamboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop =  mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu))))
  }
  
  out <- model_after
  deselect_para <-  list(tau = perc, deselectmethod = method[1])
  out <-append(x = out, values = deselect_para)
  
  class(out) <- c(class(out))
  
  return(out)
}

