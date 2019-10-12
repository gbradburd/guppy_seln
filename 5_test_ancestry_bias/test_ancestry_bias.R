################################################################
################################################################
# test for bias in the deviations from predicted ancestry 
#	in a set of candidate alleles involved in local adaptation 
#	to headwater habitat
################################################################
################################################################


################################
#	function library
################################

# calculate deviation from predicted allele frequency
get.deviation.to.anc <- function(P1,P2,f_hat,x){
	delta <- 0
	if(P1 < f_hat & f_hat < P2){
		delta <- f_hat - x
	} else if(P1 > f_hat & f_hat > P2){
		delta <- x - f_hat
	}
	return(delta)
}

# identify an allele with similar expected frequencies
#	within a specified tolerance
match.freq <- function(available.indices,mean.sim.freqs,x,tol){
	possibles <- which(abs(mean.sim.freqs[available.indices] - x) < tol)
	stopifnot(length(possibles) > 0)
	pick <- sample(possibles,1)
	return(available.indices[pick])
}

# identify a set of alleles with expected frequencies 
#	similar to that of a candidate set (within a specified tolerance)
#	by sampling without replacement
get.matched.freq.set <- function(sims,freqs2match,candidate.indices,tol){
	mean.sim.freqs <- rowMeans(sims)
	available.indices <- c(1:length(mean.sim.freqs))[-candidate.indices]
	freq.matches <- rep(NA,length(freqs2match))
	match.order <- sample(1:length(candidate.indices),length(candidate.indices),replace=FALSE)
	for(i in 1:length(match.order)){
		x <- freqs2match[match.order[i]]
		freq.matches[match.order[i]] <- match.freq(available.indices,mean.sim.freqs,x,tol)
		available.indices <- available.indices[-which(available.indices==freq.matches[match.order[i]])]
	}
	return(freq.matches)
}

# calculate the deviations between observed and predicted
#	frequencies for a frequency-matched "null" set of alleles
get.matched.freq.devs <- function(P1,P2,x,sims,freqs2match,candidate.indices,tol){
	matched.freq.set <- get.matched.freq.set(sims,freqs2match,candidate.indices,tol=tol)
	matched.freq.devs <- sapply(matched.freq.set,function(l){
							get.deviation.to.anc(P1=P1[l],
							  		 			 P2=P2[l],
							  		 			 f_hat=mean(sims[l,]),
							  		 			 x=x[l])
						 })
	return(
		list("indices" = matched.freq.set,
				"devs" = matched.freq.devs)
	)
}

# assigns a number of stars based on 
#	the significance of a test
pval.stars <- function(p){
	if(p > 0.1){
		stars <- 0
	}
	if(p > 0.05 & p < 0.1){
		stars <- 1
	}
	if(p < 0.05 & p > 0.001){
		stars <- 2
	}
	if(p < 0.001 & p > 0.0001){
		stars <- 3
	}
	if(p < 0.0001){
		stars <- 4
	}
	return(stars)
}

# writes text containing result of a paired t.test
pval.text <- function(a,b){
	pval <- t.test(a,b,paired=TRUE,alternative="greater")$p.value
	txt <- sprintf("%s (p = %s)",
				paste0(rep("*",pval.stars(pval)),collapse=""),
				  signif(pval,4))
	return(txt)
}

################################
# load data objects for both streams
################################

set.seed(123)

load("../3_reconcile_alleles/guppy.freqs.Robj")
load("../3_reconcile_alleles/ADMX.CA.freqs.Robj")
load("../3_reconcile_alleles/ADMX.TY.freqs.Robj")
load("../4_silico_sims/stream.sims.Robj")
load("../2_admixture/admix.props.Robj")
nLoci <- nrow(CA.sims)


################################
# identify candidate loci 
#	involved in local adaptation to 
#	headwater habitat
################################

# set a cutoff for how similar the frequencies 
#	of a particular allele should be in the 
#	pre-gene flow headwater populations
siml.cutoff <- 0.1

# set a cutoff for how dissiilar the frequencies 
#	of a particular allele should be between each  
#	pre-gene flow headwater populations and 
#	mainstem population
diff.cutoff <- 0.9

# identify candidate loci involved in local adaptation to 
#	headwater habitat
candidates <- which(abs(guppy.freqs[["NCA"]] - guppy.freqs[["NTY"]]) < siml.cutoff & 
				  abs(guppy.freqs[["NCA"]] - guppy.freqs[["SGS"]]) > diff.cutoff & 
	  			  abs(guppy.freqs[["NTY"]] - guppy.freqs[["SGS"]]) > diff.cutoff)


################################
# calculate ancestry-polarized 
#	deviations between observed and expected
#	frequencies for candidate loci
################################

devs2anc.ca <- sapply(candidates,function(l){
				get.deviation.to.anc(P1=ADMX.CA.freqs[l,1],
							  		 P2=ADMX.CA.freqs[l,2],
							  		 f_hat=mean(CA.sims[l,]),
							  		 x=guppy.freqs[["PCA"]][l])
				})


devs2anc.ty <- sapply(candidates,function(l){
				get.deviation.to.anc(P1=ADMX.TY.freqs[l,1],
							  		 P2=ADMX.TY.freqs[l,2],
							  		 f_hat=mean(TY.sims[l,]),
							  		 x=guppy.freqs[["PTY"]][l])
				})


################################
# calculate ancestry-polarized 
#	deviations between observed and expected
#	frequencies for frequency-matched "null" loci
################################

matched.devs.ca <- get.matched.freq.devs(P1=ADMX.CA.freqs[,1],
										 P2=ADMX.CA.freqs[,2],
										 x=guppy.freqs[["PCA"]],
										 sims=CA.sims,
										 freqs2match=guppy.freqs[["PCA"]][candidates],
										 candidate.indices=candidates,
										 tol=0.05)

matched.devs.ty <- get.matched.freq.devs(P1=ADMX.TY.freqs[,1],
										 P2=ADMX.TY.freqs[,2],
										 x=guppy.freqs[["PTY"]],
										 sims=TY.sims,
										 freqs2match=guppy.freqs[["PTY"]][candidates],
										 candidate.indices=candidates,
										 tol=0.05)


################################
# visualize some results
################################

# frequencies of alleles at candidate loci 
#	in headwater and mainstem populations 
pdf(file="../6_figs/candidate_loci_frequencies.pdf",width=6,height=6)
	plot(guppy.freqs[["NCA"]][candidates],
			guppy.freqs[["SGS"]][candidates],
				pch=19,col=adjustcolor(1,0.2),
		 xlab="frequency in both pre-gene flow headwater populations\n(frequencies are the same in each)",
		 ylab="frequency in mainstem population",
		 main=sprintf("Candidate loci (n=%s)",length(candidates)))
dev.off()

# example ancestry-polarized deviation 
#	from expected frequency
pdf(file="../6_figs/example.dev.pdf",width=8,height=4)
i <- 7
	plot(c(guppy.freqs[["NCA"]][candidates][i],
		   guppy.freqs[["SGS"]][candidates][i],
		   rowMeans(CA.sims[candidates,])[i],
		   guppy.freqs[["PCA"]][candidates][i]),
		 rep(1,4),
		 	col=c(4,2,"purple","purple"),
		 	pch=c(19,19,5,18),
			ylim=c(0,2),xlim=c(0,1),yaxt='n',
			ylab="",xlab="allele frequency",
			cex=c(rep(1.5,3),2))
		arrows(x0=rowMeans(CA.sims[candidates,])[i],
			   y0=1,
			   x1=guppy.freqs[["PCA"]][candidates][i]+0.013,
			   y1=1,
			   lwd=1.5,
			   length=0.07,
			   col="blue")
		par(xpd=TRUE)
		text(x=0.07,y=1.3,labels="pre-gene flow\nheadwater",col=4)
		text(x=0.62,y=0.6,labels="post-gene flow\nheadwater\n(observed)",col="purple")
		text(x=0.76,y=1.5,labels="post-gene flow\nheadwater\n(expected)",col="purple")
		text(x=0.92,y=0.78,labels="mainstem",col="red")
		legend(x="bottomleft",legend=c("deviation toward headwater","deviation toward mainstem"),lty=1,lwd=3,col=c(4,2))
dev.off()

# ancestry-polarized deviations for all 
# 	candidate loci in Caigual
pdf(file="../6_figs/candidate_devs.ca.pdf",width=12,height=6)
	par(mar=c(4,4.5,6,1))
	plot(1:length(candidates),ylim=c(0,1),type='n',xlab="candidate loci",ylab="frequency",cex.lab=1.5)
		points(1:length(candidates),guppy.freqs[["NCA"]][candidates],pch=20,col="blue")
		points(1:length(candidates),guppy.freqs[["SGS"]][candidates],pch=20,col="red")
		points(1:length(candidates),rowMeans(CA.sims[candidates,]),pch=5,col="purple")
		points(1:length(candidates),guppy.freqs[["PCA"]][candidates],pch=18,col="purple",cex=1.6)
		arrows(x0=1:length(candidates),
				 y0=rowMeans(CA.sims[candidates,]),
				 x1=1:length(candidates),
				 y1=guppy.freqs[["PCA"]][candidates],
				 lty=1,
				 lwd=1,
				 col=ifelse(devs2anc.ca>0,"blue","red"),
				 length=0.07)
	mtext(side=3,adj=0,padj=-0.1,
			text="Caigual: deviation in candidate loci\ntoward headwater ancestry",
			font=2,cex=1.5)
	par(xpd=TRUE)
		legend(x=103.9,y=1.35,
				pch=c(19,19,5,18,NA,NA),
				lty=c(NA,NA,NA,NA,1,1),
				col=c(4,2,"purple","purple",4,2),
				pt.cex=c(1,1,1,1.6,NA,NA),
				legend=c("frequency in pre-gene flow headwater",
						 "frequency in mainstem",
						 "expected frequency",
						 "frequency in post-gene flow headwater",
						 "deviation from mean toward headwater freq",
						 "deviation from mean toward mainstem freq"),
			  	cex=0.75)
dev.off()

# ancestry-polarized deviations for all 
# 	frequency-matched null loci in Caigual
pdf(file="../6_figs/null_devs.ca.pdf",width=12,height=6)
	par(mar=c(4,4.5,6,1))
	plot(1:length(matched.devs.ca$indices),ylim=c(0,1),type='n',xlab="null loci",ylab="frequency",cex.lab=1.5)
		points(1:length(matched.devs.ca$indices),guppy.freqs[["NCA"]][matched.devs.ca$indices],pch=20,col="blue")
		points(1:length(matched.devs.ca$indices),guppy.freqs[["SGS"]][matched.devs.ca$indices],pch=20,col="red")
		points(1:length(matched.devs.ca$indices),rowMeans(CA.sims[matched.devs.ca$indices,]),pch=5,col="purple")
		points(1:length(matched.devs.ca$indices),guppy.freqs[["PCA"]][matched.devs.ca$indices],pch=18,col="purple",cex=1.6)
		arrows(x0=1:length(matched.devs.ca$indices),
				 y0=rowMeans(CA.sims[matched.devs.ca$indices,]),
				 x1=1:length(matched.devs.ca$indices),
				 y1=guppy.freqs[["PCA"]][matched.devs.ca$indices],
				 lty=1,
				 lwd=1,
				 col=ifelse(matched.devs.ca$devs>0,"blue","red"),
				 length=0.07)
	mtext(side=3,adj=0,padj=-0.1,
			text="Caigual: deviation in null loci\ntoward headwater ancestry",
			font=2,cex=1.5)
	par(xpd=TRUE)
		legend(x=103.9,y=1.35,
				pch=c(19,19,5,18,NA,NA),
				lty=c(NA,NA,NA,NA,1,1),
				col=c(4,2,"purple","purple",4,2),
				pt.cex=c(1,1,1,1.6,NA,NA),
				legend=c("frequency in pre-gene flow headwater",
						 "frequency in mainstem",
						 "expected frequency",
						 "frequency in post-gene flow headwater",
						 "deviation from mean toward headwater freq",
						 "deviation from mean toward mainstem freq"),
			  	cex=0.75)
dev.off()

# ancestry-polarized deviations for all 
# 	candidate loci in Taylor
pdf(file="../6_figs/candidate_devs.ty.pdf",width=12,height=6)
	par(mar=c(4,4.5,6,1))
	plot(1:length(candidates),ylim=c(0,1),type='n',xlab="candidate loci",ylab="frequency",cex.lab=1.5)
		points(1:length(candidates),guppy.freqs[["NTY"]][candidates],pch=20,col="blue")
		points(1:length(candidates),guppy.freqs[["SGS"]][candidates],pch=20,col="red")
		points(1:length(candidates),rowMeans(TY.sims[candidates,]),pch=5,col="purple")
		points(1:length(candidates),guppy.freqs[["PTY"]][candidates],pch=18,col="purple",cex=1.6)
		arrows(x0=1:length(candidates),
				 y0=rowMeans(TY.sims[candidates,]),
				 x1=1:length(candidates),
				 y1=guppy.freqs[["PTY"]][candidates],
				 lty=1,
				 lwd=1,
				 col=ifelse(devs2anc.ty>0,"blue","red"),
				 length=0.07)
	mtext(side=3,adj=0,padj=-0.1,
			text="Taylor: deviation in candidate loci\ntoward headwater ancestry",
			font=2,cex=1.5)
	par(xpd=TRUE)
		legend(x=103.9,y=1.35,
				pch=c(19,19,5,18,NA,NA),
				lty=c(NA,NA,NA,NA,1,1),
				col=c(4,2,"purple","purple",4,2),
				pt.cex=c(1,1,1,1.6,NA,NA),
				legend=c("frequency in pre-gene flow headwater",
						 "frequency in mainstem",
						 "expected frequency",
						 "frequency in post-gene flow headwater",
						 "deviation from mean toward headwater freq",
						 "deviation from mean toward mainstem freq"),
			  	cex=0.75)
dev.off()

# ancestry-polarized deviations for all 
# 	frequency-matched null loci in Taylor
pdf(file="../6_figs/null_devs.ty.pdf",width=12,height=6)
	par(mar=c(4,4.5,6,1))
	plot(1:length(matched.devs.ty$indices),ylim=c(0,1),type='n',xlab="null loci",ylab="frequency",cex.lab=1.5)
		points(1:length(matched.devs.ty$indices),guppy.freqs[["NTY"]][matched.devs.ty$indices],pch=20,col="blue")
		points(1:length(matched.devs.ty$indices),guppy.freqs[["SGS"]][matched.devs.ty$indices],pch=20,col="red")
		points(1:length(matched.devs.ty$indices),rowMeans(TY.sims[matched.devs.ty$indices,]),pch=5,col="purple")
		points(1:length(matched.devs.ty$indices),guppy.freqs[["PTY"]][matched.devs.ty$indices],pch=18,col="purple",cex=1.6)
		arrows(x0=1:length(matched.devs.ty$indices),
				 y0=rowMeans(TY.sims[matched.devs.ty$indices,]),
				 x1=1:length(matched.devs.ty$indices),
				 y1=guppy.freqs[["PTY"]][matched.devs.ty$indices],
				 lty=1,
				 lwd=1,
				 col=ifelse(matched.devs.ty$devs>0,"blue","red"),
				 length=0.07)
	mtext(side=3,adj=0,padj=-0.1,
			text="Taylor: deviation in null loci\ntoward headwater ancestry",
			font=2,cex=1.5)
	par(xpd=TRUE)
		legend(x=103.9,y=1.35,
				pch=c(19,19,5,18,NA,NA),
				lty=c(NA,NA,NA,NA,1,1),
				col=c(4,2,"purple","purple",4,2),
				pt.cex=c(1,1,1,1.6,NA,NA),
				legend=c("frequency in pre-gene flow headwater",
						 "frequency in mainstem",
						 "expected frequency",
						 "frequency in post-gene flow headwater",
						 "deviation from mean toward headwater freq",
						 "deviation from mean toward mainstem freq"),
			  	cex=0.75)
dev.off()

# ancestry-polarized deviations for candidate loci 
#	compared to null set for both drainages
pdf(file="../6_figs/devs2anc.pdf",width=12,height=6)
	par(mfrow=c(1,2))
		hist(matched.devs.ca$devs,
				xlim=range(matched.devs.ca$devs,devs2anc.ca)+c(-0.03,0.03),
				xlab="",
				main="Caigual",
				ylab="",
				freq=FALSE,
				breaks=8,
				col=adjustcolor(1,0.4))
			abline(v=mean(matched.devs.ca$devs),col="white",lwd=4)
			abline(v=mean(matched.devs.ca$devs),col=1,lwd=3)
			hist(devs2anc.ca,col=adjustcolor(4,0.45),freq=FALSE,add=TRUE,breaks=8)
			abline(v=mean(devs2anc.ca),col="white",lwd=4)
			abline(v=mean(devs2anc.ca),col=4,lwd=3)
			text(x=0.11,y=7.9,labels=pval.text(devs2anc.ca,matched.devs.ca$devs))
		legend(x="topleft",pch=c(15,15,NA),lty=c(NA,NA,1),
				pt.cex=c(3,3,NA),
				lwd=c(NA,NA,4),
				col=c(adjustcolor(4,0.4),adjustcolor(1,0.4),1),
				legend=c("candidate loci","frequency-matched\nnon-candidates","mean"),
				cex=0.8,
				bty='n')
		mtext(text="counts",
				side=2,padj=-3.5,adj=0.5,cex=1.25)
		hist(matched.devs.ty$devs,
				xlim=range(matched.devs.ty$devs,devs2anc.ty)+c(-0.03,0.03),
				xlab="",
				main="Taylor",
				ylab="",
				freq=FALSE,
				col=adjustcolor(1,0.4))
			abline(v=mean(matched.devs.ty$devs),col="white",lwd=4)
			abline(v=mean(matched.devs.ty$devs),col=1,lwd=3)
			hist(devs2anc.ty,col=adjustcolor(4,0.45),freq=FALSE,add=TRUE)
			abline(v=mean(devs2anc.ty),col="white",lwd=4)
			abline(v=mean(devs2anc.ty),col=4,lwd=3)
			text(x=0.14,y=7.64,labels=pval.text(devs2anc.ty,matched.devs.ty$devs))
		mtext(text="mean deviation toward headwater frequencies",
				side=1,padj=4,adj=-5.75,cex=1.25)
dev.off()

# one panel figure that contains results of example
pdf(file="../6_figs/ancestry_deviation_Sfig.pdf",width=8.5,height=11)
layout(matrix(c(rep(1,8),rep(2,12),rep(3,12)),nrow=8,ncol=4,byrow=TRUE))
	################
	# Example Figure
	################
	par(mar=c(4.5,4.5,5.5,1))
	i <- 7
	plot(c(guppy.freqs[["NCA"]][candidates][i],
		   guppy.freqs[["SGS"]][candidates][i],
		   rowMeans(CA.sims[candidates,])[i],
		   guppy.freqs[["PCA"]][candidates][i]),
		 rep(1,4),
		 	col=c(4,2,"purple","purple"),
		 	pch=c(19,19,5,18),
			ylim=c(0,2),xlim=c(0,1),yaxt='n',
			ylab="",xlab="",
			cex=c(rep(1.5,3),2),
			cex.axis=1.5)
		arrows(x0=rowMeans(CA.sims[candidates,])[i],
			   y0=1,
			   x1=guppy.freqs[["PCA"]][candidates][i]+0.013,
			   y1=1,
			   lwd=1.5,
			   length=0.07,
			   col="blue")
		mtext(side=3,text="Deviation in candidate loci\ntoward headwater ancestry",font=2,cex=1.5,padj=-0.1)
	legend(x="bottomleft",legend=c("deviation toward headwater","deviation toward mainstem"),lty=1,lwd=3,col=c(4,2))
		par(xpd=TRUE)
		text(x=0.07,y=1.3,labels="pre-gene flow\nheadwater",col=4,cex=1.2)
		text(x=0.62,y=0.55,labels="post-gene flow\nheadwater\n(observed)",col="purple",cex=1.2)
		text(x=0.81,y=1.5,labels="post-gene flow\nheadwater\n(expected)",col="purple",cex=1.2)
		text(x=0.917,y=0.79,labels="mainstem",col="red",cex=1.2)
	mtext(side=3,adj=0,padj=-0.1,
			text="Example locus",
			font=2,cex=1.2)
	mtext(side=1,padj=2.5,cex=1.2,text="allele frequency")
	################
	# Caigual
	################
	par(mar=c(3.5,4.5,2,1))
	plot(1:length(candidates),ylim=c(0,1),type='n',xlab="",ylab="",cex.axis=1.5)
		points(1:length(candidates),guppy.freqs[["NCA"]][candidates],pch=20,col="blue")
		points(1:length(candidates),guppy.freqs[["SGS"]][candidates],pch=20,col="red")
		points(1:length(candidates),rowMeans(CA.sims[candidates,]),pch=5,col="purple")
		points(1:length(candidates),guppy.freqs[["PCA"]][candidates],pch=18,col="purple",cex=1.6)
		arrows(x0=1:length(candidates),
				 y0=rowMeans(CA.sims[candidates,]),
				 x1=1:length(candidates),
				 y1=guppy.freqs[["PCA"]][candidates],
				 lty=1,
				 lwd=1,
				 col=ifelse(devs2anc.ca>0,"blue","red"),
				 length=0.07)
	mtext(side=3,adj=0,padj=-0.1,
			text="Caigual",
			font=2,cex=1.2)
	mtext(side=2,adj=-0.3,padj=-2.4,
			text="frequency",
			font=1,cex=1.3)
	################
	# Taylor
	################
	par(mar=c(4.5,4.5,1,1))
	plot(1:length(candidates),ylim=c(0,1),type='n',xlab="",ylab="",cex.axis=1.5)
		points(1:length(candidates),guppy.freqs[["NTY"]][candidates],pch=20,col="blue")
		points(1:length(candidates),guppy.freqs[["SGS"]][candidates],pch=20,col="red")
		points(1:length(candidates),rowMeans(TY.sims[candidates,]),pch=5,col="purple")
		points(1:length(candidates),guppy.freqs[["PTY"]][candidates],pch=18,col="purple",cex=1.6)
		arrows(x0=1:length(candidates),
				 y0=rowMeans(TY.sims[candidates,]),
				 x1=1:length(candidates),
				 y1=guppy.freqs[["PTY"]][candidates],
				 lty=1,
				 lwd=1,
				 col=ifelse(devs2anc.ty>0,"blue","red"),
				 length=0.07)
	mtext(side=3,adj=0,padj=-0.1,
			text="Taylor",
			font=2,cex=1.2)
	mtext(side=1,text="candidate loci",padj=2.5,
			font=1,cex=1.2)
dev.off()