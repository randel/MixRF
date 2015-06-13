#' Impute a large number of genes using the MixRF algorithm with parallel computing
#'
#' This function impute the expression of a large number of genes using the MixRF algorithm with parallel computing.
#'
#' @param Ydat An array of expression data of dimension sample-by-gene-by-tissue, nxpxT, where n is sample size.
#'   p is the number of genes, and T is the number of tissues.  Ydat[,1,] is a matrix of the first gene expression
#'   in T tissues for n individuals, nxT. Ydat[,,1] is a nxp matrix of the expression data of p genes in the first tissue.
#' @param eqtl.lis A list of eQTL names of length p. Each of the element in the list contains the name of the eQTLs for 
#'   the corresponding gene. The order of the list should correspond to the order of genes in Ydat.
#' @param snp.dat A matrix of genotype. Each row is a sample and each column corresponds to one SNP. The column names should match eqtl.lis.
#' @param cov A matrix of covariates. Each row is a sample and each column corresponds to one covariate. For example, age, gender.
#' @param iPC The option is used only when iPC=TRUE. When it is,  the imputed PCs (iPCs) for each tissue type will be constructed 
#'   based on the combined observed and imputed data on the selected genes. The iPCs will be adjusted as covariates 
#'   in the imputation. 
#' @param idx.selected.gene.iPC The option is used only when iPC=TRUE. When it is, one may select a subset of genes and impute 
#'   those first to construct iPCs.
#' @param parallel.size A numerical value specifying the number of CPUs/cores available for parallel computing.
#'
#' @return An nxpxT array of imputed and observed expression data. The observed values in Ydat are still kept and the missing values in Ydat are imputed.

#' @export
#' @examples
#' \dontrun{
#' data(sim)
#' 
#' idx.selected.gene.iPC = which(sapply(sim$eqtl.lis,length)>=1)
#' 
#' Yimp = MixRF.impute(sim$Ydat, sim$eqtl.lis, sim$snp.dat, sim$cov, iPC=TRUE, idx.selected.gene.iPC,
#'                     parallel.size=4)
#' }
#'
MixRF.impute <- function(Ydat, eqtl.lis, snp.dat, cov=NULL, iPC=TRUE, idx.selected.gene.iPC=NULL, 
                         parallel.size=1){
  
  MixRF.impute.single <- function(Yi, eqtl.lis.i, snp.dat, cov=NULL, iPC.cov=NULL){
    
    MixRF <- function(data, initialRandomEffects=0, 
                      ErrorTolerance=0.001, MaxIterations=1000) {
      
      Target = data$Y
      ID = data$ID
      
      # Condition that indicates the loop has not converged or run out of iterations
      ContinueCondition <- TRUE
      
      iterations <- 0
      
      # Get initial values
      AdjustedTarget <- Target - initialRandomEffects
      oldvar <- -Inf
      
      # y - X*beta
      resi = rep(NA,length(Target))
      
      if(ncol(data)==3) {
        x = as.data.frame(data[,-c(1:2)])
        colnames(x) = colnames(data)[3]
      } else x = data[,-c(1:2)]
      
      
      while(ContinueCondition){
        
        iterations <- iterations+1
        
        # randomForest 
        rf = randomForest(x,AdjustedTarget)
        
        # out-of-bag prediction
        resi = Target - rf$predicted
        
        ## Estimate New Random Effects and Errors using lmer
        lmefit <- lmer(resi~1+(1|ID), data=data)
        
        # check on convergence
        newvar <- attr(VarCorr(lmefit)$'ID','stddev')^2
        
        ContinueCondition <- (newvar-oldvar>ErrorTolerance & iterations < MaxIterations)
        oldvar <- newvar
        
        # Extract random effects to make the new adjusted target
        AllEffects <- predict(lmefit)-fixef(lmefit)
        #  y-Zb
        AdjustedTarget <- Target - AllEffects
      }
      
      result <- list(forest=rf, EffectModel=lmefit, RandomEffects=ranef(lmefit),
                     IterationsUsed=iterations, ErrorTolerance=ErrorTolerance)
      
      return(result)
    }
    
    
    
    predict.MixRF <- function(object, newdata, id=NULL, EstimateRandomEffects=TRUE){
      
      # Base predictions from the forest part
      if(ncol(newdata)==3) {
        x = as.data.frame(newdata[,-c(1:2)])
        colnames(x) = colnames(newdata)[3]
      } else x = newdata[,-c(1:2)]
      
      treePrediction <- predict(object$forest,x)
      
      # If we aren't estimating random effects, we 
      # just use the forest for prediction.
      if(!EstimateRandomEffects){
        return(treePrediction)
      }
      
      # Get the group identifiers if necessary
      if(is.null(id)){
        id <- as.matrix(object$EffectModel$groups)
      }
      
      n = length(unique(id))
      nT = length(treePrediction)/n
      
      completePrediction <- matrix(treePrediction, n, nT)
      
      # Get the identities of the groups in the data
      uniqueID <- unique(id)
      
      estRE <- matrix(0, nrow=n, ncol=1)
      oRE=object$RandomEffects[[1]]
      idx.re=as.integer(rownames(oRE))
      estRE[idx.re,1] <- as.vector(oRE[[1]])
      
      for(i in uniqueID){
        completePrediction[i,] = completePrediction[i,] + estRE[i,1]
      }
      
      return(completePrediction)
    }
    
    
    n=nrow(Yi)
    nT=ncol(Yi)
    
    # transform covariates array into a matrix
    covMat=NULL
    if (!is.null(iPC.cov)){
      np = ncol(iPC.cov)
      for (i in 1:nT){
        oo= matrix(0, nrow=nT*n, ncol=np)
        oo[n*(i-1)+1:n, ] <- iPC.cov[,,i]
        covMat = cbind(covMat,oo)
      }
    }
    
    y_mean = apply(Yi,2,mean,na.rm=T)
    
    Ymat = matrix(Yi-matrix(rep(y_mean,n),byrow=T,nrow=n),ncol=1)
    ID = rep(1:n, nT)
    
    # Marginal effects: regress each tissue and gene on each eQTL one by one
    ID2 = rep(1:nT, rep(n,nT))
    
    comCovMat <- NULL
    if (!is.null(cov)) {
      for (i in 1:nT) comCovMat=rbind(comCovMat, cov)
    }
    
    XTmat = comCovMat
    
    ## include all eQTLs, then adaptively weight them
    
    
    if(!is.null(eqtl.lis.i)) {
      eqtl.lis.i <- unlist(eqtl.lis.i)
      XTmat= cbind(XTmat, matrix(snp.dat[,eqtl.lis.i], ncol=length(eqtl.lis.i))[ID,])
    }
    
    ## adaptive.weight
    XTmat_marginal = XTmat
    
    for (i in 1:nT){
      for(j in 1:ncol(XTmat)){
        if(sum(!is.na(Ymat[ID2==i]))>0){
          try({lm_marginal = lm(Ymat[ID2==i] ~ XTmat[ID2==i,j],na.action=na.exclude)},silent=T)
          weight = coef(lm_marginal)[2]
          
          # coef may be NA since all X are 0s
          if(!is.na(weight)) XTmat_marginal[ID2==i,j] = XTmat[ID2==i,j] * weight
        }
      }
    }
    
    XX = cbind(covMat, XTmat_marginal)
    
    newdata = as.data.frame(cbind(ID, Ymat, XX))
    colnames(newdata) <- c('ID',"Y",paste0("X",1:ncol(XX)))
    xnam = colnames(newdata)[-c(1:2)]
    
    dat = newdata[!is.na(Ymat),]
    
    ### MixRF
    REEMresult <- MixRF(data=dat)
    fitted.y = predict.MixRF(REEMresult, newdata, id=newdata$ID, EstimateRandomEffects=TRUE)
    
    #   resid.y = Ymat - as.vector(fitted.y)
    
    Ypred= fitted.y + rep(y_mean,rep(n,nT))
    
    return(Ypred)
  }
  
  
  snp.dat = snp.dat[,unique(unlist(eqtl.lis))]
  
  #set up parallel backend
  cl <- makeCluster(parallel.size)
  registerDoParallel(cl)
  
  if (iPC==TRUE){
    g = length(idx.selected.gene.iPC)
    
    #loop
    ls <- foreach(i=idx.selected.gene.iPC) %dopar% {
      MixRF.impute.single(Yi=Ydat[,i,], eqtl.lis=eqtl.lis[[i]], snp.dat=snp.dat, cov=cov,
                          iPC.cov=NULL)
    }
    
    Yimp = Ydat[,idx.selected.gene.iPC,]
    for(i in 1:g){
      Yimp[,i,] = ls[[i]]
    }
    ### observed values
    Yimp[!is.na(Ydat[,idx.selected.gene.iPC,])] = Ydat[,idx.selected.gene.iPC,][!is.na(Ydat[,idx.selected.gene.iPC,])]
    
    # get PCs
    npc = 5
    nT=dim(Yimp)[3]
    n=dim(Yimp)[1]
    iPC.cov = array(NA,dim=c(n,npc,nT))
    
    for (ti in 1:nT) {
      gid=1:(dim(Yimp)[2])
      Yt=Yimp[, gid, ti]
      av.sampt= !is.na(Yt[,1])
      ev.ti=matrix(NA, n, npc)
      ss= svd(t(Yt[av.sampt,]))
      ev.ti[av.sampt,] = ss$v[,1:npc]
      
      iPC.cov[,,ti] = ev.ti
    }
    
  }
  
  g = ncol(Ydat)
  
  #loop
  ls <- foreach(i=1:g) %dopar% {
    MixRF.impute.single(Yi=Ydat[,i,], eqtl.lis=eqtl.lis[[i]], snp.dat, cov=cov, 
                        iPC.cov=iPC.cov)
  }
  
  Yimp = Ydat
  for(i in 1:g){
    Yimp[,i,] = ls[[i]]
  }
  ### observed values
  Yimp[!is.na(Ydat)] = Ydat[!is.na(Ydat)]
  
  stopCluster(cl)
  return(Yimp)
}
