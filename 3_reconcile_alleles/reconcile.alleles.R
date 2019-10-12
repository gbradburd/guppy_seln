################################################################
################################################################
#	script for reconciling counted alleles between datasets and 
#		ADMIXTURE output
################################################################
################################################################

################################
# read in stream freqs
################################
CA.freqs <- read.table("../1_data/CA.frq.strat",stringsAsFactors=FALSE,header=TRUE)
TY.freqs <- read.table("../1_data/TY.frq.strat",stringsAsFactors=FALSE,header=TRUE)
CA.pop.freqs <- lapply(c("NCA","SGS","PCA"),
					function(x){
						CA.freqs$MAF[which(CA.freqs$CLST==x)]
					})
names(CA.pop.freqs) <- c("NCA","SGS","PCA")
TY.pop.freqs <- lapply(c("NTY","SGS","PTY"),
					function(x){
						TY.freqs$MAF[which(TY.freqs$CLST==x)]
					})
names(TY.pop.freqs) <- c("NTY","SGS","PTY")


################################
# read in ADMIXTURE freqs
################################

#ADMIXTURE apparently models the frequency of the major allele, 
#	and not the minor allele, 
#	so to get Minor AF just subtract from 1
ADMX.CA.freqs <- 1 - read.table("../2_admixture/CA_admx.2.P",stringsAsFactors=FALSE,header=FALSE)
ADMX.TY.freqs <- 1 - read.table("../2_admixture/TY_admx.2.P",stringsAsFactors=FALSE,header=FALSE)

################################
# plot for sanity check
################################
pdf(file="../6_figs/freq_check.pdf",width=12,height=6.5)
	layout(matrix(c(1,2,5,5,3,4,5,5),nrow=2,ncol=4,byrow=TRUE))
		plot(CA.pop.freqs$NCA,ADMX.CA.freqs[,1],xlab="pre_CA",ylab="CA_clst_1") ; abline(0,1,col=2)
		plot(CA.pop.freqs$SGS,ADMX.CA.freqs[,2],xlab="mainstem",ylab="CA_clst_2") ; abline(0,1,col=2)

		plot(TY.pop.freqs$NTY,ADMX.TY.freqs[,1],xlab="pre_TY",ylab="TY_clst_1") ; abline(0,1,col=2)
		plot(TY.pop.freqs$SGS,ADMX.TY.freqs[,2],xlab="mainstem",ylab="TY_clst_2") ; abline(0,1,col=2)

		plot(CA.pop.freqs$SGS,TY.pop.freqs$SGS,xlab="mainstem in CA dataset",ylab="mainstem in TY dataset") ; abline(0,1,col=2)
dev.off()

################################
# reconcile freq files
################################

# Notice that, because the specific makeup of 
#	each ADMIXTURE run is different, there are 
#	some major alleles in one analysis that are 
#	being modeled as minor alleles in the other
#	
# We need to make sure that we're being consistent 
#	in which allele we're tracking across the independent 
#	runs, so we need to reconcile these differences

# using CA as our reference
prob.loci <- which(CA.pop.freqs$SGS != TY.pop.freqs$SGS)

ADMX.TY.freqs[prob.loci,] <- 1 - ADMX.TY.freqs[prob.loci,]
TY.pop.freqs$SGS[prob.loci] <- 1 - TY.pop.freqs$SGS[prob.loci]
TY.pop.freqs$NTY[prob.loci] <- 1 - TY.pop.freqs$NTY[prob.loci]
TY.pop.freqs$PTY[prob.loci] <- 1 - TY.pop.freqs$PTY[prob.loci]

pdf(file="../6_figs/reconciled_freq_check.pdf",width=12,height=6.5)
	layout(matrix(c(1,2,5,5,3,4,5,5),nrow=2,ncol=4,byrow=TRUE))
		plot(CA.pop.freqs$NCA,ADMX.CA.freqs[,1],xlab="pre_CA",ylab="CA_clst_1") ; abline(0,1,col=2)
		plot(CA.pop.freqs$SGS,ADMX.CA.freqs[,2],xlab="mainstem",ylab="CA_clst_2") ; abline(0,1,col=2)

		plot(TY.pop.freqs$NTY,ADMX.TY.freqs[,1],xlab="pre_TY",ylab="TY_clst_1") ; abline(0,1,col=2)
		plot(TY.pop.freqs$SGS,ADMX.TY.freqs[,2],xlab="mainstem",ylab="TY_clst_2") ; abline(0,1,col=2)

		plot(CA.pop.freqs$SGS,TY.pop.freqs$SGS,xlab="mainstem in CA dataset",ylab="mainstem in TY dataset") ; abline(0,1,col=2)
dev.off()

# the alleles are now reconciled between datasets

################################
# plot summary of observed/estimated 
#	allele frequencies for different populations
################################

pdf(file="../6_figs/ADMIXTURE_freqs_comp.pdf",width=14,height=5)
	#quartz(width=14,height=5)
	layout(matrix(c(1,2,5,5,6,6,3,4,5,5,6,6),nrow=2,ncol=6,byrow=TRUE))
		par(mar=c(4.5,5,4.5,2.5),mgp=c(2,0.5,0))
		plot(CA.pop.freqs$NCA,ADMX.CA.freqs[,1],xaxt='n',yaxt='n',cex.lab=1.5,xlab="frequency in pre_CA",ylab="frequency in CA\nCluster 1",pch=20,col=adjustcolor(1,0.2),xlim=c(0,1),ylim=c(0,1)) ; abline(0,1,col=2,lwd=0.5,lty=2)
			axis(side=1,at=c(0,0.5,1)) ; axis(side=2,at=c(0,0.5,1))
		plot(CA.pop.freqs$SGS,ADMX.CA.freqs[,2],xaxt='n',yaxt='n',cex.lab=1.5,xlab="frequency in mainstem",ylab="frequency in CA\nCluster 2",pch=20,col=adjustcolor(1,0.2),xlim=c(0,1),ylim=c(0,1)) ; abline(0,1,col=2,lwd=0.5,lty=2)
			axis(side=1,at=c(0,0.5,1)) ; axis(side=2,at=c(0,0.5,1))
			mtext(side=3,text="Caigual: observed vs. estimated",font=2,adj=1.2,padj=-2,cex=1.2)
		plot(TY.pop.freqs$NTY,ADMX.TY.freqs[,1],xaxt='n',yaxt='n',cex.lab=1.5,xlab="frequency in pre_TY",ylab="frequency in TY\nCluster 1",pch=20,col=adjustcolor(1,0.2),xlim=c(0,1),ylim=c(0,1)) ; abline(0,1,col=2,lwd=0.5,lty=2)
			axis(side=1,at=c(0,0.5,1)) ; axis(side=2,at=c(0,0.5,1))
		plot(TY.pop.freqs$SGS,ADMX.TY.freqs[,2],xaxt='n',yaxt='n',cex.lab=1.5,xlab="frequency in mainstem",ylab="frequency in TY\nCluster 2",pch=20,col=adjustcolor(1,0.2),xlim=c(0,1),ylim=c(0,1)) ; abline(0,1,col=2,lwd=0.5,lty=2)
			axis(side=1,at=c(0,0.5,1)) ; axis(side=2,at=c(0,0.5,1))
			mtext(side=3,text="Taylor: observed vs. estimated",font=2,adj=1.2,padj=-2,cex=1.2)
		plot(ADMX.CA.freqs[,1],ADMX.TY.freqs[,1],cex.lab=1.5,xlab="frequency in Cluster 1 in CA dataset",ylab="frequency in Cluster 1 in TY dataset",pch=20,col=adjustcolor(1,0.2),xlim=c(0,1),ylim=c(0,1)) ; abline(0,1,col=2,lwd=0.5,lty=2)
		plot(ADMX.CA.freqs[,2],ADMX.TY.freqs[,2],cex.lab=1.5,xlab="frequency in Cluster 2 in CA dataset",ylab="frequency in Cluster 2 in TY dataset",pch=20,col=adjustcolor(1,0.2),xlim=c(0,1),ylim=c(0,1)) ; abline(0,1,col=2,lwd=0.5,lty=2)
			mtext(side=3,text="Comparing estimated cluster frequencies",font=2,adj=1.2,padj=-0.8,cex=2)
dev.off()

################################
# save reconciled freq files
################################
guppy.freqs <- list("SGS" = CA.pop.freqs$SGS,
					"NCA" = CA.pop.freqs$NCA,
					"PCA" = CA.pop.freqs$PCA,
					"NTY" = TY.pop.freqs$NTY,
					"PTY" = TY.pop.freqs$PTY)

save(guppy.freqs,file="guppy.freqs.Robj")

save(ADMX.CA.freqs,file="ADMX.CA.freqs.Robj")
save(ADMX.TY.freqs,file="ADMX.TY.freqs.Robj")
