#' Calculate cis- and trans-eQTLs
#' 
#' @param ncore The number of cores for parallel computing.
#' @param Ynew  An array of expression data of dimension sample-by-gene-by-tissue, nxpxT, where n is sample size.
#'   p is the number of genes, and T is the number of tissues.
#' @param ssnpDat The genotype data matrix (n by SNP size).
#' @param snp.info Input for MatrixEQTL, with col.names snpID,  chr,  pos.
#' @param gene.info Input for MatrixEQTL, with col.names geneID, chr,  lpos, rpos.
#' @param cov The covariates matrix for MatrixEQTL.
#'
#' @return A list contains the cis- and trans-eQTLs for each gene.

#' @export
#' @examples
#' 
#' library(parallel)
#' library(MatrixEQTL)
#' eqtl_list = get_eqtl(ncore=2, Ynew, ssnpDat, snp.info, gene.info, cov)


get_eqtl = function(ncore, Ynew, ssnpDat, snp.info, gene.info, cov) {
  
  library(parallel)
  cl = makeCluster(ncore)
  tmp = capture.output({clusterEvalQ(cl, {
    library(MatrixEQTL)
    source("eqtl.r")
  })})
  
  nT = dim(Ynew)[3]
  
  
  ###### Stouffer's cis-list
  tmp = capture.output({gen.eQTL.files(cl=cl, Ynew, ssnpDat, snp.info, gene.info, cov=cov, 
                                       cisDist=1e6, cis.p=1, trans.p=0, filenam='cis')})
  cis_stouffer = get_commonCis(cl, nT=nT)
  
  names(cis_stouffer) <- cisnam <- sapply(cis_stouffer,function(x) x[1,1])
  cis_list <- list()
  
  gnam <- as.vector(gene.info$geneID)
  colnames(Ynew) <- gnam
  
  nGene=dim(Ynew)[2]
  
  cis_list[[length(gnam)]] <- character(0)
  names(cis_list) <- gnam
  oolist <- stouffer2cisList(cis_stouffer, pcut=1e-6)
  cis_list[names(oolist)] <- oolist
  
  
  # trans-list
  eqtl.new.lis <- list()
  snam <- colnames(ssnpDat)
  for (j in 1:ceiling(length(snam)/5000)) {
    print(j)
    idx <- (j-1)*5000+1:5000
    idx <- idx[idx<=length(snam)]
    snamj <- snam[idx]
    snp.infoS <- snp.info[idx,]
    ssnpDatS <- ssnpDat[,idx]
    # already parallized in the tissue dimension
    tmp = capture.output({gen.eQTL.files(cl=cl, Ynew, ssnpDatS, snp.infoS, gene.info, cov=cov, cisDist=1e6, cis.p=0.05, trans.p=0.05, 
                                         filenam='trans')})
    
    ## trans-eQTLs, note that only focus on the genes
    eqtl.new = NULL
    eqtl.new[[length(gnam)]] <- character(0)
    
    try(eqtl.new <- get_transOver(geneID=gnam, filenam='trans', nT=nT), silent=T)
    eqtl.new.lis[[j]] <- eqtl.new
  }
  
  
  trans_list <- list()
  for (i in 1:nGene) trans_list[[i]] <- unique(unlist(sapply(eqtl.new.lis,function(x) x[[i]])))
  
  
  eqtl.comb <- combo(cis_list, trans_list)
  
  stopCluster(cl)
  
  return(eqtl.comb)
}


## generate files for MatrixEQTL

gen.eQTL.dat <- function(ti, Ynew, ssnpDat, snp.info, gene.info, cov,geneID, cisDist,cis.p, trans.p, filenam="temp") {
  
  av.samp = which(!is.na(Ynew[,1,ti]))
  
  snpsMat= t(ssnpDat[av.samp,])
  snps = SlicedData$new()
  snps$CreateFromMatrix(snpsMat)
  
  ID=colnames(snpsMat)
  exprMat=t(Ynew[av.samp, ,ti])
  colnames(exprMat) <- ID
  rownames(exprMat) <- geneID
  expr = SlicedData$new()
  expr$CreateFromMatrix(exprMat)
  
  cvrtMat = t(cov[av.samp,])
  colnames(cvrtMat) <- ID
  cvrt = SlicedData$new()
  cvrt$CreateFromMatrix(cvrtMat)
  
  options(MatrixEQTL.dont.preserve.gene.object = TRUE);
  
  
  ### Run Matrix eQTL
  me = Matrix_eQTL_main(
    snps = snps,
    gene = expr,
    cvrt = cvrt,
    output_file_name = paste0(filenam,ti,"_trans.txt"),
    pvOutputThreshold = trans.p,
    useModel = modelLINEAR,
    errorCovariance = numeric(),
    verbose = TRUE,
    output_file_name.cis = paste0(filenam,ti,"_cis.txt"),
    pvOutputThreshold.cis = cis.p,
    snpspos = snp.info,
    genepos = gene.info,
    cisDist = cisDist,      #### cis window = 1mb this is the default in MatrixEQTL
    pvalue.hist = FALSE,
    noFDRsaveMemory = TRUE);
}


# cl is for parallel computing, generated from the parallel package
gen.eQTL.files <- function(cl=NULL, Ynew, ssnpDat, snp.info, gene.info, cov=cov, cisDist=1e6, cis.p=0.01, trans.p=10^(-5), 
                           filenam="temp") {
  
  nT=dim(Ynew)[3]
  geneID=as.vector(gene.info$geneID)
  
  if (is.null(cl)){
    for (ti in 1:nT) gen.eQTL.dat(ti, Ynew=Ynew, ssnpDat=ssnpDat, snp.info=snp.info, gene.info=gene.info,
                                  cov=cov, geneID=geneID, cisDist=cisDist,cis.p=cis.p, trans.p=trans.p, filenam=filenam)
  } else {
    parLapply(cl, 1:nT, gen.eQTL.dat, Ynew=Ynew, ssnpDat=ssnpDat, snp.info=snp.info, gene.info=gene.info,
              cov=cov, geneID=geneID, cisDist=cisDist, cis.p=cis.p, trans.p=trans.p, filenam=filenam)
  }
}



## Stouffer's method to combine tissue-specific eQTLs

Stouffer.test <- function(z, w) { # z is a vector of Z-stat
  if (missing(w)) {
    w <- rep(1, length(z))/length(z)
  } else {
    if (length(w) != length(z))
      stop("Length of z and w must equal!")
  }
  
  Z  <- sum(w*z)/sqrt(sum(w^2))
  p.val <- 2*(1-pnorm(abs(Z)))
  return(p.val)
}



## calculate cis-eQTLs
# for 1 gene
get_cis1gene <- function(gene_i, XX, nover=1) {
  XXi = matrix(XX[XX[,2]==gene_i,], ncol=3)
  snp = names(which(table(XXi[,1])>=nover))
  pval = as.data.frame(matrix(NA,length(snp),3))
  pval[,1] = rep(gene_i,length(snp))
  pval[,2] = snp
  
  j = 0
  for(i in snp){
    j = j + 1
    oo=as.numeric(XXi[XXi[,1]==i,3])
    pval[j,3] = Stouffer.test(z=oo)
  }
  return(pval)
}


# parallel version

get_commonCis <- function(cl=NULL, nT=9, nover=1){
  
  filenam="cis"
  typ='cis'
  XX = NULL
  for (ti in 1:nT){
    xx = matrix(scan(paste0(filenam,ti,"_",typ,".txt"), what=""), byrow=T, ncol=5)[, c(1,2,4)]
    XX = rbind(XX,xx[-1,])
  }
  
  gene = unique(XX[,2])
  
  system.time(tmp <- parLapply(cl,gene,get_cis1gene,XX=XX, nover=nover))
  
  return(tmp)
}


# convert results for Stouffer's method to an eQTL list
stouffer2cisList <- function(eQTL_stouffer, pcut){
  elist = list()
  idx = NULL
  j = 0
  for(i in 1:length(eQTL_stouffer)){
    index = eQTL_stouffer[[i]][,3]<=pcut
    if(sum(index)>0) {
      idx = c(idx,eQTL_stouffer[[i]][1,1]) 
      j = j + 1
      elist[[j]] = eQTL_stouffer[[i]][index,2]
    }
  }
  names(elist) <- idx
  return(elist)    
}



## calculate trans-eQTLs

get_transOver <- function(geneID, filenam, nT=9, pcut=NULL){
  typ='trans'
  
  xx2lis <- list()
  for (ti in 1:nT){
    xx = matrix(scan(paste0(filenam,ti,"_",typ,".txt"), what=""), byrow=T, ncol=5)[, c(1,2,5)]  
    xx=xx[-1,]
    if (is.null(pcut)) {
      xx=xx[,-3]
    } else {
      idx=which(as.numeric(xx[,3])<=pcut)
      xx=xx[idx,-3]
    } 
    
    xx2 <- NULL
    for (i in 1:ceiling(nrow(xx)/10000)){
      idx <- (i-1)*10000+1:10000
      idx <- idx[idx<=nrow(xx)]
      xx2= c(xx2, apply(xx[idx,],1,paste0,collapse=";"))
    } 
    xx2lis[[ti]] <- xx2
  }  
  
  overlap <- NULL
  for (i in 1:nT) {
    # find common elements from nT-1 tissues
    overlap <- c(overlap, Reduce(intersect,xx2lis[-i]))
  }  
  overlap <- unique(overlap)
  overlap <- matrix(unlist(strsplit(overlap,split=";")),byrow=T,ncol=2)
  
  etlis <- list()
  etlis[[length(geneID)]] <- character(0)
  for (i in 1:nrow(overlap)){
    ii <- which(geneID==overlap[i,2])
    etlis[[ii]] <- c(etlis[[ii]], overlap[i,1])
  }
  names(etlis) <- geneID
  
  return(etlis)
}



# combine cis, trans eQTLs
combo <- function(eQTL.lis.cis, eQTL.lis.trans) {
  eQTL.lis.all <- list()
  for (i in 1:length(eQTL.lis.cis)){
    elis=list()
    uu <- unique(c(eQTL.lis.cis[[i]], eQTL.lis.trans[[i]]))
    uu <- uu[!is.na(uu)]
    elis <- uu
    eQTL.lis.all[[i]] <- elis
  }
  return(eQTL.lis.all)
}
