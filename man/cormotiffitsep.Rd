\name{cormotiffitsep}
\alias{cormotiffitsep}
\title{Individual Study Motif Fit}
\description{This function fits a mixture modified t-distribution model to each study seperately. }

\usage{
cormotiffitsep(exprs,groupid,compid, tol=1e-3, max.iter=100)
}

\arguments{
\item{exprs}{ a matrix, the expression data after normalization that is on log2 scale, each row of the matrix corresponds to a gene
and each column of the matrix corresponds to a sample array. }
\item{groupid}{ the group label for each sample array, two arrays in the same study with same experinment condition\eqn{(e.g. control)} have the same groupid. }
\item{compid}{ the study design and comparison matrix, each row of the matrix corresponds to one study with the first column being the first experinment condition and the second column being the second experinment condition}
\item{tol}{ the relative tolerance level of error. }
\item{max.iter}{ maximun number of iterations. }
}

\value{
\item{p.post}{ the posterior probability for each gene to be differentially expressed. }
\item{motif.prior}{ fitted values of the probability for genes to be differentially expressed in each study, a \eqn{1*D} vector, where \eqn{D} is the number of studies}
\item{loglike}{ log-likelihood of the fitted model. }
}

\author{Hongkai Ji, Yingying Wei}

\references{Ji, H., Wei, Y.,(2011) Correlation Motif. Unpublished}

\examples{
data(simudata2)
n<-nrow(simudata2)
m<-ncol(simudata2)
#the expression data is from the second column to m
exprs.simu2<-as.matrix(simudata2[,2:m])

#prepare the group ID number for each sample array
data(simu2_groupid)

#prepare the design matrix for each group of samples
data(simu2_compgroup)

#fit seperate models to each study
motif.fitted.sep<-cormotiffitsep(exprs.simu2, simu2_groupid,simu2_compgroup)
}
