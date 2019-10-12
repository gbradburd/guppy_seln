################################################################
################################################################
#	run admixture simulations
################################################################
################################################################

################################
# function library
################################

# simulate a chromosome
make.a.chromo <- function(w,p){
	#recover()
	L <- nrow(p)
	l1 <- sample(1:L,w[1]*L)
	l2 <- c(1:L)[-l1]
	chromo <- rep(NA,L)
	chromo[l1] <- rbinom(length(l1),1,p[l1,1])
	chromo[l2] <- rbinom(length(l2),1,p[l2,2])
	return(chromo)
}

# simulate a diploid individual
make.an.ind <- function(w,p){
	geno <- make.a.chromo(w,p) + 
			make.a.chromo(w,p)
	return(geno)
}

# simulate a population sample of diploid individuals
make.a.stream <- function(w,p){
	stream <- t(apply(w,1,
				function(a){
					make.an.ind(a,p)
			  	}))
	return(stream)
}

################################
# load data objects for both streams
################################

load("../2_admixture/admix.props.Robj")
load("../3_reconcile_alleles/guppy.freqs.Robj")
load("../3_reconcile_alleles/ADMX.CA.freqs.Robj")
load("../3_reconcile_alleles/ADMX.TY.freqs.Robj")
n.samples <- list("NCA" = length(which(grepl("NCA",row.names(CA.w)))),
				  "NTY" = length(which(grepl("NTY",row.names(TY.w)))),
				  "SGS" = length(which(grepl("SGS",row.names(TY.w)))),
				  "PCA" = length(which(grepl("PCA",row.names(CA.w)))),
				  "PTY" = length(which(grepl("PTY",row.names(TY.w)))))


################################
# run 1000 simulations for each 
#	post-gene flow headwater population
################################

set.seed(123)

# for Caigual
CA.sims <- replicate(1000,
			colMeans(
				make.a.stream(
					CA.w[which(grepl("PCA",row.names(CA.w))),],
					ADMX.CA.freqs)
			)/2)

# for Taylor
TY.sims <- replicate(1000,
			colMeans(
				make.a.stream(
					TY.w[which(grepl("PTY",row.names(TY.w))),],
					ADMX.TY.freqs)
			)/2)

# save stream simulations
save(CA.sims,TY.sims,file="stream.sims.Robj")


################################
# visualize simulation output
################################

pdf(file="../6_figs/sim_freqs.pdf",width=10,height=5)
	par(mfrow=c(1,2))
	plot(guppy.freqs[["PCA"]],CA.sims[,1],
			pch=20,col=adjustcolor(1,0.1),
			xlim=c(0,1),ylim=c(0,1),
			xlab="sampled post_CA frequencies",
			ylab="simulated post_CA frequencies")
		abline(0,1,col=2)
	plot(guppy.freqs[["PTY"]],CA.sims[,1],
			pch=20,col=adjustcolor(1,0.1),
			xlim=c(0,1),ylim=c(0,1),
			xlab="sampled post_TY frequencies",
			ylab="simulated post_TY frequencies")
		abline(0,1,col=2)
		mtext(text="Simulation freqs match sample freqs",font=2,cex=1.3,side=3,adj=21,padj=-2)
dev.off()


# calculate empirical 95% quantiles for simulation output 
#	and compare to observed frequencies
sim.CIs.ca <- apply(CA.sims,1,function(x){quantile(x,c(0.025,0.975))})
sim.CIs.ty <- apply(TY.sims,1,function(x){quantile(x,c(0.025,0.975))})

pdf(file="../6_figs/sim_freqs_CIs.pdf",width=10,height=5)
	par(mfrow=c(1,2))
	plot(guppy.freqs[["PCA"]],type='n',
			xlim=c(0,1),ylim=c(0,1),
			xlab="sampled post_CA frequencies",
			ylab="simulated post_CA frequencies")
		segments(x0=guppy.freqs[["PCA"]],
				 y0=sim.CIs.ca[1,],
				 x1=guppy.freqs[["PCA"]],
				 y1=sim.CIs.ca[2,],
				 lwd=0.5)
		abline(0,1,col=2)
		legend(x="topleft",
				lwd=1,
				col=1,
				legend="95% quantile",
				cex=0.8)
	plot(guppy.freqs[["PTY"]],type='n',
			xlim=c(0,1),ylim=c(0,1),
			xlab="sampled post_TY frequencies",
			ylab="simulated post_TY frequencies")
		segments(x0=guppy.freqs[["PTY"]],
				 y0=sim.CIs.ty[1,],
				 x1=guppy.freqs[["PTY"]],
				 y1=sim.CIs.ty[2,],
				 lwd=0.5)
		abline(0,1,col=2)
dev.off()





