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
# source('MixRF.r')
#'
#' library(lme4)
#' library(randomForest)
#' data(sleepstudy)
#'
#' tmp = MixRF(Y=sleepstudy$Reaction, x=as.data.frame(sleepstudy$Days), random='(Days|Subject)',
#'             data=sleepstudy, initialRandomEffects=0, ErrorTolerance=0.01, MaxIterations=100)
#'
#' pred = predict.MixRF(object=tmp, newdata=sleepstudy, EstimateRE=TRUE)

predict.MixRF = function(object, newdata, id=NULL, EstimateRE=TRUE){

  forestPrediction = predict(object$forest, newdata)

  # If we aren't estimating random effects, we
  # just use the forest for prediction.
  if(!EstimateRE) return(forestPrediction)

  # Get the group identifiers if necessary
  if(is.null(id)) id = as.matrix(object$EffectModel$groups)

  n = length(unique(id))
  nT = length(forestPrediction)/n

  completePrediction = matrix(forestPrediction, n, nT)

  # Get the identities of the groups in the data
  uniqueID = unique(id)

  estRE = matrix(0, nrow=n, ncol=1)
  oRE = object$RandomEffects[[1]]
  idx.re = as.integer(rownames(oRE))
  estRE[idx.re,1] = as.vector(oRE[[1]])

  for(i in uniqueID){
    completePrediction[i,] = completePrediction[i,] + estRE[i,1]
  }

  return(completePrediction)
}
