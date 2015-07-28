#' Mixed Random Forest
#' 
#' @param Y The outcome variable.
#' @param x A data frame contains the predictors.
#' @param random A string in lme4 format indicates the random effect model.
#' @param data The data set as a data frame.
#' @param initialRandomEffects The initial values for random effects. 
#' @param ErrorTolerance The tolerance for log-likelihood.
#' @param MaxIterations The maximum iteration times.
#'
#' @return A list contains the random forest, mixed model, and random effects. 
#' See the example below for the usage. A predict() function is available for request.

#' @export
#' @examples
#' 
# source('MixRF.r')
# 
# library(lme4)
# library(randomForest)
# data(sleepstudy)
# 
# tmp = MixRF(Y=sleepstudy$Reaction, x=as.data.frame(sleepstudy$Days), random='(Days|Subject)',
#             data=sleepstudy, initialRandomEffects=0, ErrorTolerance=0.01, MaxIterations=100)
# 
# tmp$forest
# tmp$MixedModel
# tmp$RandomEffects


MixRF <- function(Y, x, random, data, initialRandomEffects=0, 
                  ErrorTolerance=0.001, MaxIterations=1000) {
  
  Target = Y
  
  # Condition that indicates the loop has not converged or run out of iterations
  ContinueCondition <- TRUE
  
  iterations <- 0
  
  # Get initial values
  AdjustedTarget <- Target - initialRandomEffects
  oldLogLik <- -Inf
  
  while(ContinueCondition){
    
    iterations <- iterations+1
    
    # randomForest 
    rf = randomForest(x,AdjustedTarget)
    
    # y - X*beta (out-of-bag prediction)
    resi = Target - rf$predicted
    
    ## Estimate New Random Effects and Errors using lmer
    f0 = as.formula(paste0('resi~',random))
    lmefit <- lmer(f0, data=data)
    
    # check convergence
    newLogLik <- as.numeric(logLik(lmefit))
    
    ContinueCondition <- (abs(newLogLik-oldLogLik)>ErrorTolerance & iterations < MaxIterations)
    oldLogLik <- newLogLik
    
    # Extract random effects to make the new adjusted target
    AllEffects <- predict(lmefit)-fixef(lmefit)
    
    #  y-Zb
    AdjustedTarget <- Target - AllEffects
  }
  
  result <- list(forest=rf, MixedModel=lmefit, RandomEffects=ranef(lmefit),
                 IterationsUsed=iterations)
  
  return(result)
}
