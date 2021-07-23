## Fit limma model using code as it is found in the original cormotif code. It has 
## only been modified to add names to the matrix of t values, as well as the
## limma fits

limmafit.default <- function(exprs,groupid,compid) {
	limmafits  <- list()
	compnum    <- nrow(compid)
	genenum    <- nrow(exprs)	
	limmat     <- matrix(0,genenum,compnum)
	limmas2    <- rep(0,compnum)
	limmadf    <- rep(0,compnum)
	limmav0    <- rep(0,compnum)
	limmag1num <- rep(0,compnum)
	limmag2num <- rep(0,compnum)

	rownames(limmat)  <- rownames(exprs)
	colnames(limmat)  <- rownames(compid)
	names(limmas2)    <- rownames(compid)
	names(limmadf)    <- rownames(compid)
	names(limmav0)    <- rownames(compid)
	names(limmag1num) <- rownames(compid)
	names(limmag2num) <- rownames(compid)
	
	for(i in 1:compnum) {
		selid1 <- which(groupid == compid[i,1]) 
		selid2 <- which(groupid == compid[i,2])
		eset   <- new("ExpressionSet", exprs=cbind(exprs[,selid1],exprs[,selid2]))
		g1num  <- length(selid1)
		g2num  <- length(selid2)
		designmat <- cbind(base=rep(1,(g1num+g2num)), delta=c(rep(0,g1num),rep(1,g2num)))
		fit <- lmFit(eset,designmat)
		fit <- eBayes(fit)
		limmat[,i] <- fit$t[,2]
		limmas2[i] <- fit$s2.prior
		limmadf[i] <- fit$df.prior
		limmav0[i] <- fit$var.prior[2]
		limmag1num[i] <- g1num
		limmag2num[i] <- g2num
		limmafits[[i]] <- fit

		# log odds
		# w<-sqrt(1+fit$var.prior[2]/(1/g1num+1/g2num))
		# log(0.99)+dt(fit$t[1,2],g1num+g2num-2+fit$df.prior,log=TRUE)-log(0.01)-dt(fit$t[1,2]/w, g1num+g2num-2+fit$df.prior, log=TRUE)+log(w)
	}
	names(limmafits) <- rownames(compid)
    limmacompnum<-nrow(compid)
	result<-list(t       = limmat, 
	             v0      = limmav0, 
	             df0     = limmadf, 
	             s20     = limmas2, 
	             g1num   = limmag1num, 
	             g2num   = limmag2num,
	             compnum = limmacompnum,
	             fits    = limmafits)
}

limmafit.counts <- 
  function (exprs, groupid, compid, norm.factor.method = "TMM", voom.normalize.method = "none") 
  {
	limmafits  <- list()
	compnum    <- nrow(compid)
	genenum    <- nrow(exprs)	
	limmat     <- matrix(NA,genenum,compnum)
	limmas2    <- rep(0,compnum)
	limmadf    <- rep(0,compnum)
	limmav0    <- rep(0,compnum)
	limmag1num <- rep(0,compnum)
	limmag2num <- rep(0,compnum)

	rownames(limmat)  <- rownames(exprs)
	colnames(limmat)  <- rownames(compid)
	names(limmas2)    <- rownames(compid)
	names(limmadf)    <- rownames(compid)
	names(limmav0)    <- rownames(compid)
	names(limmag1num) <- rownames(compid)
	names(limmag2num) <- rownames(compid)
	
    for (i in 1:compnum) {
      message(paste("Running limma for comparision",i,"/",compnum))
      selid1 <- which(groupid == compid[i, 1])
      selid2 <- which(groupid == compid[i, 2])
      # make a new count data frame
      counts <- cbind(exprs[, selid1], exprs[, selid2])
      
      # remove NAs 
      not.nas <- which(apply(counts, 1, function(x) !any(is.na(x))) == TRUE)
      
      # runn voom/limma
      d <- DGEList(counts[not.nas,])
      d <- calcNormFactors(d, method = norm.factor.method)
      g1num <- length(selid1)
      g2num <- length(selid2)
      designmat <- cbind(base = rep(1, (g1num + g2num)), delta = c(rep(0, 
                                                                       g1num), rep(1, g2num)))
      
      y <- voom(d, designmat, normalize.method = voom.normalize.method)
      fit <- lmFit(y, designmat)
      fit <- eBayes(fit)
       
      limmafits[[i]] <- fit
      limmat[not.nas, i] <- fit$t[, 2]
      limmas2[i] <- fit$s2.prior
      limmadf[i] <- fit$df.prior
      limmav0[i] <- fit$var.prior[2]
      limmag1num[i] <- g1num
      limmag2num[i] <- g2num
    }
    limmacompnum <- nrow(compid)
    names(limmafits) <- rownames(compid)
    result <- list(t       = limmat, 
                   v0      = limmav0, 
                   df0     = limmadf, 
                   s20     = limmas2,
                   g1num   = limmag1num, 
                   g2num   = limmag2num, 
                   compnum = limmacompnum,
                   fits    = limmafits)
  }
  
limmafit.list <- 
  function (fitlist, cmp.idx=2) 
  {
   compnum    <- length(fitlist)

   genes <- c()
   for (i in 1:compnum) genes <- unique(c(genes, rownames(fitlist[[i]])))

   genenum    <- length(genes)	
   limmat     <- matrix(NA,genenum,compnum)
   limmas2    <- rep(0,compnum)
   limmadf    <- rep(0,compnum)
   limmav0    <- rep(0,compnum)
   limmag1num <- rep(0,compnum)
   limmag2num <- rep(0,compnum)
   
   rownames(limmat)  <- genes
   colnames(limmat)  <- names(fitlist)
   names(limmas2)    <- names(fitlist)
   names(limmadf)    <- names(fitlist)
   names(limmav0)    <- names(fitlist)
   names(limmag1num) <- names(fitlist)
   names(limmag2num) <- names(fitlist)
   	
   for (i in 1:compnum) {
     this.t <- fitlist[[i]]$t[,cmp.idx]
     limmat[names(this.t),i] <- this.t
     
     limmas2[i]    <- fitlist[[i]]$s2.prior
     limmadf[i]    <- fitlist[[i]]$df.prior
     limmav0[i]    <- fitlist[[i]]$var.prior[cmp.idx]
     limmag1num[i] <- sum(fitlist[[i]]$design[,cmp.idx]==0)
     limmag2num[i] <- sum(fitlist[[i]]$design[,cmp.idx]==1)
   }
   
   limmacompnum <- compnum
   result <- list(t       = limmat, 
                  v0      = limmav0, 
                  df0     = limmadf, 
                  s20     = limmas2,
                  g1num   = limmag1num, 
                  g2num   = limmag2num, 
                  compnum = limmacompnum,
                  fits    = limmafits)
      
  }
  
## Rank genes based on statistics
generank<-function(x) {
	xcol<-ncol(x)
	xrow<-nrow(x)
	result<-matrix(0,xrow,xcol)
	z<-(1:1:xrow)
	for(i in 1:xcol) {
		y<-sort(x[,i],decreasing=TRUE,na.last=TRUE)
		result[,i]<-match(x[,i],y)
		result[,i]<-order(result[,i])
	}
	result
}

## Log-likelihood for moderated t under H0
modt.f0.loglike<-function(x,df) {
	a<-dt(x, df, log=TRUE)
	result<-as.vector(a)
      	flag<-which(is.na(result)==TRUE)
      	result[flag]<-0
      	result
}

## Log-likelihood for moderated t under H1
## param=c(df,g1num,g2num,v0)
modt.f1.loglike<-function(x,param) {
	df<-param[1]
	g1num<-param[2]
	g2num<-param[3]
	v0<-param[4]
	w<-sqrt(1+v0/(1/g1num+1/g2num))
	dt(x/w, df, log=TRUE)-log(w)
	a<-dt(x/w, df, log=TRUE)-log(w)
      	result<-as.vector(a)
      	flag<-which(is.na(result)==TRUE)
      	result[flag]<-0
      	result
}

## Correlation Motif Fit
cmfit<-function(x, type, K=1, tol=1e-3, max.iter=100) {
	## initialize
	xrow <- nrow(x)			
    xcol <- ncol(x)			
    loglike0 <- list()			
    loglike1 <- list()			
    p <- rep(1, K)/K			
    q <- matrix(runif(K * xcol), K, xcol)			
    q[1, ] <- rep(0.01, xcol)			
    for (i in 1:xcol) {			
        f0 <- type[[i]][[1]]			
        f0param <- type[[i]][[2]]			
        f1 <- type[[i]][[3]]			
        f1param <- type[[i]][[4]]			
        loglike0[[i]] <- f0(x[, i], f0param)			
        loglike1[[i]] <- f1(x[, i], f1param)			
    }			
    condlike <- list()			
    for (i in 1:xcol) {			
        condlike[[i]] <- matrix(0, xrow, K)			
    }			
    loglike.old <- -1e+10			
    for (i.iter in 1:max.iter) {			
        err <- tol + 1			
        clustlike <- matrix(0, xrow, K)			

        templike1 <- vector("numeric",xrow)
        templike2 <- vector("numeric",xrow)
        for (j in 1:K) {			
            for (i in 1:xcol) {			
                templike1 <- log(q[j, i]) + loglike1[[i]]			
                templike2 <- log(1 - q[j, i]) + loglike0[[i]]			
                
                tempmax <- Rfast::Pmax(templike1, templike2)
                
                templike1 <- exp(templike1 - tempmax)			
                templike2 <- exp(templike2 - tempmax)			
                
                tempsum <- templike1 + templike2
                clustlike[, j] <- clustlike[, j] + tempmax + 			
                  log(tempsum)			
                condlike[[i]][, j] <- templike1/tempsum			
            }			
            clustlike[, j] <- clustlike[, j] + log(p[j])			
        }			
		
        tempmax <- Rfast::rowMaxs(clustlike)
        for (j in 1:K) {			
            clustlike[, j] <- exp(clustlike[, j] - tempmax)			
        }
		
        tempsum <- Rfast::rowsums(clustlike)	
        for (j in 1:K) {			
            clustlike[, j] <- clustlike[, j]/tempsum			
        }			
	
        p.new <- (Rfast::colsums(clustlike) + 1)/(xrow + K)		
        q.new <- matrix(0, K, xcol)			
        for (j in 1:K) {			
            clustpsum <- sum(clustlike[, j])			
            for (i in 1:xcol) {			
                q.new[j, i] <- (sum(clustlike[, j] * condlike[[i]][, 			
                  j]) + 1)/(clustpsum + 2)			
            }			
        }			
        err.p <- max(abs(p.new - p)/p)			
        err.q <- max(abs(q.new - q)/q)			
        err <- max(err.p, err.q)			
        loglike.new <- (sum(tempmax + log(tempsum)) + sum(log(p.new)) + 			
            sum(log(q.new) + log(1 - q.new)))/xrow			
        p <- p.new			
        q <- q.new			
        loglike.old <- loglike.new			
        if (err < tol) {			
            break			
        }			
    }			
    clustlike <- matrix(0, xrow, K)			
    for (j in 1:K) {			
        for (i in 1:xcol) {			
            templike1 <- log(q[j, i]) + loglike1[[i]]			
            templike2 <- log(1 - q[j, i]) + loglike0[[i]]			
            
            tempmax <- Rfast::Pmax(templike1, templike2)
            			
            templike1 <- exp(templike1 - tempmax)	
            templike2 <- exp(templike2 - tempmax)	
            	
            tempsum <- templike1 + templike2			
            clustlike[, j] <- clustlike[, j] + tempmax + log(tempsum)			
            condlike[[i]][, j] <- templike1/tempsum			
        }			
        clustlike[, j] <- clustlike[, j] + log(p[j])			
    }			
		
    tempmax <- Rfast::rowMaxs(clustlike)	
    for (j in 1:K) {			
        clustlike[, j] <- exp(clustlike[, j] - tempmax)			
    }			
		
    tempsum <- Rfast::rowsums(clustlike)
    for (j in 1:K) {			
        clustlike[, j] <- clustlike[, j]/tempsum			
    }			
    p.post <- matrix(0, xrow, xcol)			
    for (j in 1:K) {			
        for (i in 1:xcol) {			
            p.post[, i] <- p.post[, i] + clustlike[, j] * condlike[[i]][, 			
                j]			
        }			
    }			
    loglike.old <- loglike.old - (sum(log(p)) + sum(log(q) + 			
        log(1 - q)))/xrow			
    loglike.old <- loglike.old * xrow			
    result <- list(p.post = p.post, motif.prior = p, motif.q = q, 			
        loglike = loglike.old)		
}

## Fit using (0,0,...,0) and (1,1,...,1)
cmfitall<-function(x, type, tol=1e-3, max.iter=100) {
	## initialize
	xrow<-nrow(x)
	xcol<-ncol(x)
	loglike0<-list()
	loglike1<-list()
	p<-0.01
	
	## compute loglikelihood
	L0<-matrix(0,xrow,1)
	L1<-matrix(0,xrow,1)
	for(i in 1:xcol) {
		f0<-type[[i]][[1]]
		f0param<-type[[i]][[2]]
		f1<-type[[i]][[3]]
		f1param<-type[[i]][[4]]
		loglike0[[i]]<-f0(x[,i],f0param)
		loglike1[[i]]<-f1(x[,i],f1param)
		L0<-L0+loglike0[[i]]
		L1<-L1+loglike1[[i]]
	}


	## EM algorithm to get MLE of p and q
	loglike.old <- -1e10
	for(i.iter in 1:max.iter) {
		if((i.iter%%50) == 0) {
			print(paste("We have run the first ", i.iter, " iterations",sep=""))
		}
		err<-tol+1

		## compute posterior cluster membership
		clustlike<-matrix(0,xrow,2)
		clustlike[,1]<-log(1-p)+L0
		clustlike[,2]<-log(p)+L1
		
		tempmax<-apply(clustlike,1,max)
		for(j in 1:2) {
			clustlike[,j]<-exp(clustlike[,j]-tempmax)
		}
		tempsum<-apply(clustlike,1,sum)
		
	
		## update motif occurrence rate
		for(j in 1:2) {
			clustlike[,j]<-clustlike[,j]/tempsum
		}
	
		p.new<-(sum(clustlike[,2])+1)/(xrow+2)
		
		## evaluate convergence
		err<-abs(p.new-p)/p
		

		## evaluate whether the log.likelihood increases
		loglike.new<-(sum(tempmax+log(tempsum))+log(p.new)+log(1-p.new))/xrow
			
		loglike.old<-loglike.new
		p<-p.new

		if(err<tol) {
			break;
		}
	}

	## compute posterior p
	clustlike<-matrix(0,xrow,2)
	clustlike[,1]<-log(1-p)+L0
	clustlike[,2]<-log(p)+L1
		
	tempmax<-apply(clustlike,1,max)
	for(j in 1:2) {
		clustlike[,j]<-exp(clustlike[,j]-tempmax)
	}
	tempsum<-apply(clustlike,1,sum)

	for(j in 1:2) {
		clustlike[,j]<-clustlike[,j]/tempsum
	}

	p.post<-matrix(0,xrow,xcol)
	for(i in 1:xcol) {
		p.post[,i]<-clustlike[,2]
	}
	
	## return

	#calculate back loglikelihood
	loglike.old<-loglike.old-(log(p)+log(1-p))/xrow
	loglike.old<-loglike.old*xrow
	result<-list(p.post=p.post, motif.prior=p, loglike=loglike.old)
}

## Fit each dataset separately
cmfitsep<-function(x, type, tol=1e-3, max.iter=100) {
	## initialize
	xrow<-nrow(x)
	xcol<-ncol(x)
	loglike0<-list()
	loglike1<-list()
	p<-0.01*rep(1,xcol)
	loglike.final<-rep(0,xcol)	

	## compute loglikelihood
	for(i in 1:xcol) {
		f0<-type[[i]][[1]]
		f0param<-type[[i]][[2]]
		f1<-type[[i]][[3]]
		f1param<-type[[i]][[4]]
		loglike0[[i]]<-f0(x[,i],f0param)
		loglike1[[i]]<-f1(x[,i],f1param)
	}

	p.post<-matrix(0,xrow,xcol)

	## EM algorithm to get MLE of p
	for(coli in 1:xcol) {
		loglike.old <- -1e10
		for(i.iter in 1:max.iter) {
			if((i.iter%%50) == 0) {
				print(paste("We have run the first ", i.iter, " iterations",sep=""))
			}
			err<-tol+1

			## compute posterior cluster membership
			clustlike<-matrix(0,xrow,2)
			clustlike[,1]<-log(1-p[coli])+loglike0[[coli]]
			clustlike[,2]<-log(p[coli])+loglike1[[coli]]
		
			tempmax<-apply(clustlike,1,max)
			for(j in 1:2) {
				clustlike[,j]<-exp(clustlike[,j]-tempmax)
			}
			tempsum<-apply(clustlike,1,sum)
		
			## evaluate whether the log.likelihood increases
			loglike.new<-sum(tempmax+log(tempsum))/xrow
			
			## update motif occurrence rate
			for(j in 1:2) {
				clustlike[,j]<-clustlike[,j]/tempsum
			}
	
			p.new<-(sum(clustlike[,2]))/(xrow)
		
			## evaluate convergence
			err<-abs(p.new-p[coli])/p[coli]
			loglike.old<-loglike.new
			p[coli]<-p.new

			if(err<tol) {
				break;
			}
		}

		## compute posterior p
		clustlike<-matrix(0,xrow,2)
		clustlike[,1]<-log(1-p[coli])+loglike0[[coli]]
		clustlike[,2]<-log(p[coli])+loglike1[[coli]]

		tempmax<-apply(clustlike,1,max)
		for(j in 1:2) {
			clustlike[,j]<-exp(clustlike[,j]-tempmax)
		}
		tempsum<-apply(clustlike,1,sum)

		for(j in 1:2) {
			clustlike[,j]<-clustlike[,j]/tempsum
		}

		p.post[,coli]<-clustlike[,2]
		loglike.final[coli]<-loglike.old
	}


	## return
	loglike.final<-loglike.final*xrow
	result<-list(p.post=p.post, motif.prior=p, loglike=loglike.final)
}

## Fit the full model
cmfitfull<-function(x, type, tol=1e-3, max.iter=100) {
	## initialize
	xrow<-nrow(x)
	xcol<-ncol(x)
	loglike0<-list()
	loglike1<-list()
	K<-2^xcol
	p<-rep(1,K)/K
	pattern<-rep(0,xcol)
	patid<-matrix(0,K,xcol)

	## compute loglikelihood
	for(i in 1:xcol) {
		f0<-type[[i]][[1]]
		f0param<-type[[i]][[2]]
		f1<-type[[i]][[3]]
		f1param<-type[[i]][[4]]
		loglike0[[i]]<-f0(x[,i],f0param)
		loglike1[[i]]<-f1(x[,i],f1param)
	}
	L<-matrix(0,xrow,K)
	for(i in 1:K)
	{	
		patid[i,]<-pattern
		for(j in 1:xcol) {
			if(pattern[j] < 0.5) {
				L[,i]<-L[,i]+loglike0[[j]]
			} else {
				L[,i]<-L[,i]+loglike1[[j]]
			}
		}

		if(i < K) {
			pattern[xcol]<-pattern[xcol]+1
			j<-xcol
			while(pattern[j] > 1) {
				pattern[j]<-0
				j<-j-1
				pattern[j]<-pattern[j]+1
			}
		}
	}
	
	## EM algorithm to get MLE of p and q
	loglike.old <- -1e10
	for(i.iter in 1:max.iter) {
		if((i.iter%%50) == 0) {
			print(paste("We have run the first ", i.iter, " iterations",sep=""))
		}
		err<-tol+1

		## compute posterior cluster membership
		clustlike<-matrix(0,xrow,K)
		for(j in 1:K) {
			clustlike[,j]<-log(p[j])+L[,j]
		}
		
		tempmax<-apply(clustlike,1,max)
		for(j in 1:K) {
			clustlike[,j]<-exp(clustlike[,j]-tempmax)
		}
		tempsum<-apply(clustlike,1,sum)
		

		## update motif occurrence rate
		for(j in 1:K) {
			clustlike[,j]<-clustlike[,j]/tempsum
		}

		p.new<-(apply(clustlike,2,sum)+1)/(xrow+K)
		
		## evaluate convergence
		err<-max(abs(p.new-p)/p)

		## evaluate whether the log.likelihood increases
		loglike.new<-(sum(tempmax+log(tempsum))+sum(log(p.new)))/xrow
		
		loglike.old<-loglike.new
		p<-p.new

		if(err<tol) {
			break;
		}
	}

	## compute posterior p
	clustlike<-matrix(0,xrow,K)
	for(j in 1:K) {
		clustlike[,j]<-log(p[j])+L[,j]
	}
		
	tempmax<-apply(clustlike,1,max)
	for(j in 1:K) {
		clustlike[,j]<-exp(clustlike[,j]-tempmax)
	}
	tempsum<-apply(clustlike,1,sum)

	for(j in 1:K) {
		clustlike[,j]<-clustlike[,j]/tempsum
	}

	p.post<-matrix(0,xrow,xcol)
	for(j in 1:K) {
		for(i in 1:xcol) {
			if(patid[j,i] > 0.5) {
				p.post[,i]<-p.post[,i]+clustlike[,j]
			}
		}
	}

	## return
	#calculate back loglikelihood
	loglike.old<-loglike.old-sum(log(p))/xrow
	loglike.old<-loglike.old*xrow
	result<-list(p.post=p.post, motif.prior=p, loglike=loglike.old)
}

generatetype<-function(limfitted)
{
	jtype<-list()
	df<-limfitted$g1num+limfitted$g2num-2+limfitted$df0
	for(j in 1:limfitted$compnum)
	{
  		jtype[[j]]<-list(f0=modt.f0.loglike, f0.param=df[j], f1=modt.f1.loglike, f1.param=c(df[j],limfitted$g1num[j],limfitted$g2num[j],limfitted$v0[j]))  
  	}
    	jtype
}

cormotiffit <- function(exprs, groupid=NULL, compid=NULL, K=1, tol=1e-3, 
                       max.iter=100, BIC=TRUE, norm.factor.method="TMM", 
                       voom.normalize.method = "none")
{
	# first I want to do some typechecking. Input can be either a normalized 
	# matrix, a count matrix, or a list of limma fits. Dispatch the correct
	# limmafit accordingly.
	limfitted <- list()
	if (is.matrix(exprs)) { # groupid and compid must be given
		if (is.null(groupid) | is.null(compid)) {
			stop("compid and groupid must be specified if exprs is not a list of limma fits")
		}
		if (all(X.counts == floor(X.counts), na.rm = TRUE)) {
			limfitted <- limmafit.counts(exprs,groupid,compid, norm.factor.method, voom.normalize.method)
		} else {
		    limfitted <- limmafit.default(exprs,groupid,compid)
		}
	} else if (is.list(exprs) {
		if (all(unlist(lapply(fitlist, function(i) attr(i, "class") == "MArrayLM")))) {
			limfitted <- limmafit.list(exprs)
		}
	} else {
		stop("input must be a matrix of gene expression or a list of limma fits")
	}
	
    jtype<-generatetype(limfitted)
	fitresult<-list()
	#for(i in 1:length(K))
	#	fitresult[[i]]<-cmfit(limfitted$t,type=jtype,K=K[i],max.iter=max.iter,tol=tol)
	bplapply(1:length(K), function(i) {
	  cmfit(limfitted$t,type=jtype,K=K[i],max.iter=max.iter,tol=tol)
	}) -> fitresult
	bic<-rep(0,length(K))
	aic<-rep(0,length(K))
	loglike<-rep(0,length(K))
	for(i in 1:length(K))
			loglike[i]<-fitresult[[i]]$loglike
	for(i in 1:length(K))
			bic[i]<--2*fitresult[[i]]$loglike+(K[i]-1+K[i]*limfitted$compnum)*log(dim(exprs)[1])
	for(i in 1:length(K))
			aic[i]<--2*fitresult[[i]]$loglike+2*(K[i]-1+K[i]*limfitted$compnum)
	if(BIC==TRUE)
	{
		bestflag=which(bic==min(bic))
	}
	else
	{
		bestflag=which(aic==min(aic))
	}
	result<-list(bestmotif=fitresult[[bestflag]],bic=cbind(K,bic),
			aic=cbind(K,aic),loglike=cbind(K,loglike))
	
} 

cormotiffitall<-function(exprs,groupid,compid, tol=1e-3, max.iter=100)
{
	limfitted<-limmafit(exprs,groupid,compid)
        jtype<-generatetype(limfitted)
	fitresult<-cmfitall(limfitted$t,type=jtype,tol=1e-3,max.iter=max.iter)
}

cormotiffitsep<-function(exprs,groupid,compid, tol=1e-3, max.iter=100)
{
	limfitted<-limmafit(exprs,groupid,compid)
        jtype<-generatetype(limfitted)
	fitresult<-cmfitsep(limfitted$t,type=jtype,tol=1e-3,max.iter=max.iter)
} 

cormotiffitfull<-function(exprs,groupid,compid, tol=1e-3, max.iter=100)
{
	limfitted<-limmafit(exprs,groupid,compid)
        jtype<-generatetype(limfitted)
	fitresult<-cmfitfull(limfitted$t,type=jtype,tol=1e-3,max.iter=max.iter)
} 

plotIC<-function(fitted_cormotif)
{
	oldpar<-par(mfrow=c(1,2))
	plot(fitted_cormotif$bic[,1], fitted_cormotif$bic[,2], type="b",xlab="Motif Number", ylab="BIC", main="BIC")
	plot(fitted_cormotif$aic[,1], fitted_cormotif$aic[,2], type="b",xlab="Motif Number", ylab="AIC", main="AIC")
}

plotMotif<-function(fitted_cormotif,title="")
{
	  layout(matrix(1:2,ncol=2))
          u<-1:dim(fitted_cormotif$bestmotif$motif.q)[2]
          v<-1:dim(fitted_cormotif$bestmotif$motif.q)[1]
          image(u,v,t(fitted_cormotif$bestmotif$motif.q),
          col=gray(seq(from=1,to=0,by=-0.1)),xlab="Study",yaxt = "n",
		ylab="Corr. Motifs",main=paste(title,"pattern",sep=" "))
	  axis(2,at=1:length(v))
          for(i in 1:(length(u)+1))
          {
                abline(v=(i-0.5))
          }
          for(i in 1:(length(v)+1)) 
          {
                abline(h=(i-0.5))
          }
	  Ng=10000
	  if(is.null(fitted_cormotif$bestmotif$p.post)!=TRUE)
		Ng=nrow(fitted_cormotif$bestmotif$p.post)
	  genecount=floor(fitted_cormotif$bestmotif$motif.p*Ng)
	  NK=nrow(fitted_cormotif$bestmotif$motif.q)
	  plot(0,0.7,pch=".",xlim=c(0,1.2),ylim=c(0.75,NK+0.25),
		frame.plot=FALSE,axes=FALSE,xlab="No. of genes",ylab="", main=paste(title,"frequency",sep=" "))
	  segments(0,0.7,fitted_cormotif$bestmotif$motif.p[1],0.7)
	  rect(0,1:NK-0.3,fitted_cormotif$bestmotif$motif.p,1:NK+0.3,
		col="dark grey")
	  mtext(1:NK,at=1:NK,side=2,cex=0.8)
	  text(fitted_cormotif$bestmotif$motif.p+0.15,1:NK,
	  labels=floor(fitted_cormotif$bestmotif$motif.p*Ng))
}
