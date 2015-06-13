\name{sim}
\alias{sim}
\title{
 Simulated data list
}
\description{
This simulated data list is for demonstration. 
}

\value{
\item{Ydat}{An array of expression data of dimension sample-by-gene-by-tissue, nxpxT, where n is sample size.
p is the number of genes, and T is the number of tissues.  Ydat[,1,] is a matrix of the first gene expression
in T tissues for n individuals, nxT. Ydat[,,1] is a nxp matrix of the expression data of p genes in the first tissue.}

\item{eqtl.lis}{A list of eQTL names of length p. Each of the element in the list contains the name of the eQTLs for
the corresponding gene. The order of the list should correspond to the order of genes in Ydat.}

\item{snp.dat}{A matrix of genotype. Each row is a sample and each column corresponds to one SNP. The column names should match eqtl.lis.}

\item{cov}{A matrix of covariates. Each row is a sample and each column corresponds to one covariate. For example, age, gender.}
}

\seealso{
\code{\link{MixRF.impute}}
}
