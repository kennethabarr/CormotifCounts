\name{generank}
\alias{generank}
\title{Rank genes based on statistics}
\description{This function rank the genes according to the decreasing order of the given statistics.}

\usage{
generank(x)
}

\arguments{
\item{x}{	A \eqn{G*D} matrix of statistics, the number of rows is the number of genes and the number of columns is the number of studies. }
}

\details{The function returns a \eqn{G*D} matrix of index of  top ranked genes in each study according to the decreasing order of  statistics in that study. }

\author{Hongkai Ji, Yingying Wei}

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
#give the gene index list according to the decreasing order of 
#posterior probability for a gene to be differentially expressed in each study 
generank(motif.fitted$bestmotif$p.post)
}

