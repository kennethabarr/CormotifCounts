\name{cormotiffitfull}
\alias{cormotiffitfull}
\title{Full Model Motif Fit}
\description{This function fits the data to the model with all \eqn{2^D} possible 0-1 patterns, where \eqn{D} is the number of studies. }

\usage{
cormotiffitfull(exprs,groupid,compid, tol=1e-3, max.iter=100)
}

\arguments{
\item{exprs}{ a matrix, the expression data after normalization that is on log2 scale, each row of the matrix corresponds to a gene
and each column of the matrix corresponds to a sample array. }
\item{groupid}{ the group label for each sample array, two arrays in the same study with same experinment condition\eqn{(e.g. control)} have the same groupid. }
\item{compid}{ the study design and comparison matrix, each row of the matrix corresponds to one study with the first column being the first experinment condition and the second column being the second experinment condition}
\item{tol}{ the relative tolerance level of error. }
\item{max.iter}{ maximun number of iterations. }
}

\details{The difference between \eqn{cormotiffitfull} and \eqn{cormotif(...,K=2^D,...)} is that \eqn{cormotiffitfull} forces motif to be one of the those 0-1 patterns. For \eqn{cormotiffit}, the motif does not necessarily to be of either 1 or 0, such as (0,1,..,0). It could be (0.9,0.4,...,0.2).  
}


\value{
\item{p.post}{ the posterior probability for each gene to be differentially expressed. }
\item{motif.prior}{ fitted values of the probability distribution of the \eqn{2^D} 0-1 motifs.}
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

#fit 2^D 0-1 motifs to the data
motif.fitted.sep<-cormotiffitfull(exprs.simu2, simu2_groupid,simu2_compgroup)
}
\keyword{correlation motif fit}