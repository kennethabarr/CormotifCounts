\name{plotMotif}
\alias{plotMotif}
\title{Correlation Motif plot}
\description{This function plots the Correlation Motif patterns and the associated prior probability distributions.}

\usage{
plotMotif(fitted_cormotif,title="")
}

\arguments{
\item{fitted_cormotif}{The object obtained from cormotiffit. }
\item{title}{	The title for the graph. }
}

\details{Each row in both graphs corresponds to one motif pattern. The left graph shows the correlation motif pattern. The grey color scale of cell \eqn{(k,d)} indicates the probability that motif \eqn{k} is differentially expressed in study \eqn{d}. Each row of the bar chart corresponds
to the motif pattern in the same row of the left pattern graph. The length of
the bar in the bar chart shows the number of genes of the given pattern in the
dataset, which is equal to \eqn{motif.fitted\$bestmotif\$motif.prior} multiplying
the number of total genes.}

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

#fit 2 correlation motif to the data
motif.fitted<-cormotiffit(exprs.simu2, simu2_groupid,simu2_compgroup,K=2)

plotMotif(motif.fitted)
}

