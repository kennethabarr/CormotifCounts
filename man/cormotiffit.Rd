\name{cormotiffit}
\alias{cormotiffit}
\title{Correlation Motif Fit}
\description{This function fits the Correlation Motif model to multiple expression studies. It gives the fitted values for the probability distribution of each motif, the fitted values of the given correlation matrix and the posterior probability for each gene to be differentially expressed in each study.}

\usage{
cormotiffit(exprs,groupid,compid,K=1, tol=1e-3, max.iter=100, BIC=TRUE)
}

\arguments{
\item{exprs}{ a matrix, the expression data after normalization that is on log2 scale, each row of the matrix corresponds to a gene
and each column of the matrix corresponds to a sample array. }
\item{groupid}{ the group label for each sample array, two arrays in the same study with same experinment condition\eqn{(e.g. control)} have the same groupid. }
\item{compid}{ the study design and comparison matrix, each row of the matrix corresponds to one study with the first column being the first experinment condition and the second column being the second experinment condition.}
\item{K}{ a vector, each element specifing the number of motifs a model wants to fit. }
\item{tol}{ the relative tolerance level of error. }
\item{max.iter}{ maximun number of iterations. }
\item{BIC}{	default is BIC=TRUE, selecting the model with the lowest BIC value among all fitted models; if BIC=FALSE, selecting the model with the lowest AIC value among all fitted models.}
}

\details{For the i^{th} element of \eqn{K}, the function fits total number of \eqn{K[i]} motifs to the data. Each gene can belong to one of the \eqn{K[i]} possible motifs according to prior probability distribution, \eqn{motif.prior}. For genes in motif \eqn{j}, the probability that they are differentially expressed in study \eqn{d} is \eqn{motif.q(j,d)}. One should indicate the groupid and compid for each study clearly.
}

\value{
\item{bestmotif$p.post}{ the posterior probability for each gene to be differentially expressed in each study for the best fitted model}
\item{bestmotif$motif.prior}{ fitted values of the probability distribution of different motifs for the best fitted model }
\item{bestmotif$motif.q}{ fitted values of the correlation motif matrix for the best fitted model}
\item{bestmotif$loglike}{ log-likelihood of the best fitted model }
\item{bic}{	the BIC values of all fitted models}
\item{aic}{	the AIC values of all fitted models}
\item{loglike}{ log-likelihood of all fitted models }

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

#fit 2 correlation motifs to the data
motif.fitted<-cormotiffit(exprs.simu2, simu2_groupid,simu2_compgroup,K=2)
}
\keyword{correlation motif fit}