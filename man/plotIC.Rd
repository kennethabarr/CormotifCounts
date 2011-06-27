\name{plotIC}
\alias{plotIC}
\title{BIC and AIC plot}
\description{This function plots BIC and AIC values for all fitted motif models.}

\usage{
plotIC(fitted_cormotif)
}

\arguments{
\item{fitted_cormotif}{ The object obtained from cormotiffit. }
}

\details{The left graph is the BIC plot and the right graph is the AIC plot.  }

\author{Hongkai Ji, Yingying Wei}

\references{Ji, H., Wei, Y.,(2011) Correlation Motif. Unpublished}

\examples{
data(simudata2)
n<-nrow(simudata2)
m<-ncol(simudata2)
#the expression data is from the second column to \eqn{m}
exprs.simu2<-as.matrix(simudata2[,2:m])

#prepare the group ID number for each sample array
data(simu2_groupid)

#prepare the design matrix for each group of samples
data(simu2_compgroup)

#fit 2 correlation motif to the data
motif.fitted<-cormotiffit(exprs.simu2, simu2_groupid,simu2_compgroup,K=2)

plotIC(motif.fitted)
}

