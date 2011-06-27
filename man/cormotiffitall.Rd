\name{cormotiffitall}
\alias{cormotiffitall}
\title{All Studies Correlation Motif Fit}
\description{This function assumes that a gene is either differentially expressed in all studies or is not differentially expressed in any study. It gives the fitted values for the probability distribution of motif (0,0,...0) and motif (1,1,...,1), and the posterior probability for each gene to be differentially expressed in all studies.}

\usage{
cormotiffitall(exprs,groupid,compid, tol=1e-3, max.iter=100)
}

\arguments{
\item{exprs}{ a matrix, the expression data after normalization that is on log2 scale, each row of the matrix corresponds to a gene
and each column of the matrix corresponds to a sample array. }
\item{groupid}{ the group label for each sample array, two arrays in the same study with same experinment condition\eqn{(e.g. control)} have the same groupid. }
\item{compid}{ the study design and comparison matrix, each row of the matrix corresponds to one study with the first column being the first experinment condition and the second column being the second experinment condition}
\item{tol}{ the relative tolerance level of error. }
\item{max.iter}{ maximun number of iterations. }
}

\details{The difference between \eqn{cormotiffitall} and \eqn{cormotif(...,K=2,...)} is that \eqn{cormotiffitall} forces the motif to be one of the two patterns but \eqn{cormotiffit} allows motif patterns other than (0,...,0) and (1,..,1). 
}


\value{
\item{p.post}{ the posterior probability for each gene to be differentially expressed }
\item{motif.prior}{ fitted values of the probability distribution of motif (0,0,...0) and motif (1,1,...,1)  }
\item{loglike}{ log-likelihood of the fitted model }
}

\author{Hongkai Ji, Yingying Wei}
\references{Ji, H., Wei, Y.,(2011) Correlation Motif. Unpublished}
\examples{
data(simudata2)
n<-nrow(simudata2)
m<-ncol(simudata2)
#the expression data is from the second column to m
exprs.simu2<-as.matrix(simudata2[,2:m])

#prepare the group label for each sample array
data(simu2_groupid)

#prepare the design matrix for each group of samples
data(simu2_compgroup)

#fit the two motifs (0,0,...0) and (1,1,...,1) to the data
motif.fitted.all<-cormotiffitall(exprs.simu2, simu2_groupid,simu2_compgroup)
}
