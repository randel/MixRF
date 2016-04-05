% Please edit documentation in R/eqtl.r
\name{get_eqtl}
\alias{get_eqtl}
\title{Calculate cis- and trans-eQTLs}
\usage{
get_eqtl(ncore, Ynew, ssnpDat, snp.info, gene.info, cov)
}
\arguments{
\item{ncore}{The number of cores for parallel computing.}

\item{Ynew}{An array of expression data of dimension sample-by-gene-by-tissue, nxpxT, where n is sample size.
p is the number of genes, and T is the number of tissues.}

\item{ssnpDat}{The genotype data matrix (n by SNP size).}

\item{snp.info}{Input for MatrixEQTL, with col.names snpID,  chr,  pos.}

\item{gene.info}{Input for MatrixEQTL, with col.names geneID, chr,  lpos, rpos.}

\item{cov}{The covariates matrix for MatrixEQTL.}
}
\value{
A list contains the cis- and trans-eQTLs for each gene.
}
\description{
Calculate cis- and trans-eQTLs
}
\examples{
\dontrun{
# a fake example

# eqtl_list = get_eqtl(ncore=2, Ynew, ssnpDat, snp.info, gene.info, cov)
}
}
