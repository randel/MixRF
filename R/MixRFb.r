#' Mixed Logistic Random Forest for Binary Data
#' 
#' @param Y The outcome variable.
#' @param x A formula string contains the predictors.
#' @param random A string in lme4 format indicates the random effect model.
#' @param data The data set as a data frame.
#' @param initialRandomEffects The initial values for random effects. 
#' @param ErrorTolerance The tolerance for log-likelihood.
#' @param ErrorTolerance0 The tolerance for eta (penalized quasi-likelihood, PQL).
#' @param MaxIterations The maximum iteration times for each run of PQL.
#' @param MaxIterations0 The maximum iteration times for PQL.
#' @param verbose The option to monitor each run of PQL or not.
#' 
#' @return A list contains the random forest, mixed model, and random effects. 
#' See the example below for the usage. A predict() function is also available below.

#' @export
#' @examples
#' 
#' # example data (http://stats.stackexchange.com/questions/70783/how-to-assess-the-fit-of-a-binomial-glmm-fitted-with-lme4-1-0)
# dat <- read.table("http://pastebin.com/raw.php?i=vRy66Bif")
# 
# library(party)
# library(lme4)
# 
# source('MixRFb.r')
# system.time(tmp <- MixRFb(Y=dat$true, x='factor(distance) + consequent + factor(direction) + factor(dist)', random='(1|V1)', 
#                           data=dat, initialRandomEffects=0, 
#                           ErrorTolerance=1, MaxIterations=200,
#                           ErrorTolerance0=0.3, MaxIterations0=15, verbose=T))
# 
# # tmp$forest
# # tmp$MixedModel
# # tmp$RandomEffects
# 
# # eta
# pred1 = predict.MixRF(tmp, dat, EstimateRandomEffects=TRUE)
# prob = 1/(1+exp(-pred1))
# res = (prob>.5)
# 
# # classification
# table(res,dat$true)



MixRFb <- function(Y, x, random, data, initialRandomEffects=0, 
                   ErrorTolerance=0.001, MaxIterations=200,
                   ErrorTolerance0=0.001, MaxIterations0=15, verbose=FALSE) {
  
  # Condition that indicates the loop has not converged or run out of iterations
  ContinueCondition0 <- TRUE
  iterations0 = 0
  
  # Get initial values
  
  mu = rep(mean(Y),length(Y))
  eta = log(mu/(1-mu))
  y = eta + (Y-mu)/(mu*(1-mu))
  weights = mu*(1-mu)
  
  AdjustedTarget <- y - initialRandomEffects
  
  f1 = as.formula(paste0('AdjustedTarget ~ ', x)) 
  f0 = as.formula(paste0('resi ~ -1 + ', random))
  
  # mimic randomForest's mtry
  ncol = length(strsplit(x,split="[+]")[[1]])
  mtry = if (!is.null(y) && !is.factor(y))
    max(floor(ncol/3), 1) else floor(sqrt(ncol))
  
  oldLogLik = oldEta = -Inf
  
  # PQL
  while(ContinueCondition0) {
    
    iterations0 <- iterations0 + 1
    
    iterations = 0
    ContinueCondition <- TRUE
    
    # random forest + lmer
    while(ContinueCondition) {
      
      iterations <- iterations + 1
      
      # random Forest
      data$AdjustedTarget = AdjustedTarget
      rf = cforest(f1, data=data, weights = weights, control = cforest_unbiased(mtry = mtry))
      
      # y - X*beta (out-of-bag prediction)
      pred = predict(rf, OOB = TRUE)
      resi = y - pred
      
      ## Estimate New Random Effects and Errors using lmer
      lmefit <- lmer(f0, data=data, weights=weights)
      
      # check convergence
      LogLik <- as.numeric(logLik(lmefit))
      
      ContinueCondition <- (abs(LogLik-oldLogLik)>ErrorTolerance & iterations < MaxIterations)
      oldLogLik <- LogLik
      
      # Extract (the only) random effects part (Zb) to make the new adjusted target
      AllEffects <- predict(lmefit)
      
      #  y-Zb
      AdjustedTarget <- y - AllEffects
      
      # monitor the change the of logLikelihood
      if(verbose) print(c(iterations0,iterations,LogLik))
    }
    
    eta = pred + AllEffects
    mu = 1/(1+exp(-eta))
    y = eta + (Y-mu)/(mu*(1-mu))
    AdjustedTarget <- y - AllEffects
    weights = as.vector(mu*(1-mu))
    
    print(c(iter = iterations0, maxEtaChange=max(abs(eta-oldEta))))
    
    ContinueCondition0 <- (max(abs(eta-oldEta))>ErrorTolerance0 & iterations0 < MaxIterations0)
    oldEta <- eta
  }
  
  result <- list(forest=rf, MixedModel=lmefit, RandomEffects=ranef(lmefit),
                 IterationsUsed=iterations0)
  
  return(result)
}



# predict the link transformed response (eta)

predict.MixRF <- function(object, newdata, EstimateRandomEffects=TRUE){
  
  forestPrediction <- predict(object$forest,newdata=newdata,OOB=T)
  
  # If not estimate random effects, just use the forest for prediction.
  if(!EstimateRandomEffects){
    return(forestPrediction)
  }
  
  RandomEffects <- predict(object$MixedModel, newdata=newdata, allow.new.levels=T)
  
  completePrediction = forestPrediction + RandomEffects
  
  return(completePrediction)
}
