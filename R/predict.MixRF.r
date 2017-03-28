#' Prediction Function for MixRF
#'
#' @param object The fitted MixRF object.
#' @param newdata A data frame contains the predictors for prediction.
#' @param id The group variable in the new data.
#' @param EstimateRE To use the estimated random effects in the prediction or not. The default is TRUE.
#'
#' @return A matrix (now for balanced data) contains the predicted values.

#' @export
#' @examples
#'
#' library(MixRF)
#'
#' library(lme4)
#' library(randomForest)
#' data(sleepstudy)
#'
#' X = as.data.frame(sleepstudy$Days)
#' colnames(X) = 'Days'
#' tmp = MixRF(Y=sleepstudy$Reaction, X=X, random='(Days|Subject)',
#'             data=sleepstudy, initialRandomEffects=0, ErrorTolerance=0.01, MaxIterations=100)
#'
#' pred = predict.MixRF(object=tmp, newdata=sleepstudy, EstimateRE=TRUE)

    
predict.MixRF = function(object, newdata, EstimateRE=TRUE){

  forestPrediction <- predict(object$forest, newdata=newdata, OOB=T)

  # If not estimate random effects, just use the forest for prediction.
  if(!EstimateRE) return(forestPrediction)
  
  RandomEffects <- predict(object$MixedModel, newdata=newdata, allow.new.levels=T)

  completePrediction = forestPrediction + RandomEffects

  return(completePrediction)
}
