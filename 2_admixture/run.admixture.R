################################################################
################################################################
# run ADMIXTURE on pre-contact stream datasets
################################################################
################################################################

# run a specified number of replicate ADMIXTURE analyses
run.admixture <- function(infile,n.reps,prfx){
	lapply(1:n.reps,
			function(i){
				seed <- sample(1:1e4,1)
				system(paste0("admixture -s ",seed," ",infile," ",2," > ",prfx,"_",i,".log.txt"))
				file.rename(paste0(prfx,".2.Q"),paste0(prfx,"_",i,".2.Q"))
				file.rename(paste0(prfx,".2.P"),paste0(prfx,"_",i,".2.P"))
		})
}

# get the log likelihoods associated with each replicate run
get.lnLs <- function(prfx){
	log.files <- list.files(pattern="log.txt")[which(grepl(prfx,list.files(pattern="log.txt")))]
	lnLs <- lapply(log.files,
				function(log){
					logtxt <- scan(log,what=character(),sep="\n")
					return(strsplit(logtxt[grepl("^Loglikelihood", logtxt)],": ")[[1]][2])
			})
	return(unlist(lnLs))
}

# for a set of runs, pick the one with the best log likelihood
#	and discard all others
post.process <- function(prfx){
	best <- which.max(get.lnLs(prfx))
	outfiles <- list.files(pattern=prfx)
	outfiles <- outfiles[!grepl(sprintf("_%s\\.",best),outfiles)]
	file.remove(outfiles)
	outfiles <- list.files(pattern=prfx)
	file.rename(outfiles[grepl(".P",outfiles)],sprintf("%s_admx.2.P",prfx))
	file.rename(outfiles[grepl(".Q",outfiles)],sprintf("%s_admx.2.Q",prfx))
}


################################
#	run 10 replicate ADMIXTURE analyses
################################

set.seed(123)

# for Caigual {NCA,PCA,SGS}
run.admixture("../1_data/CA_par_admx.bed",10,"CA_par_admx")
	post.process("CA")

# for Taylor {NTY,PTY,SGS}
run.admixture("../1_data/TY_par_admx.bed",10,"TY_par_admx")
	post.process("TY")


################################
#	read and save admixture proportion results
################################

CA.pa.inds <- as.matrix(read.table("../1_data/CA_par_admx.fam",stringsAsFactors=FALSE))[,1]
TY.pa.inds <- as.matrix(read.table("../1_data/TY_par_admx.fam",stringsAsFactors=FALSE))[,1]
CA.w <- as.matrix(read.table("CA_admx.2.Q",stringsAsFactors=FALSE))
	row.names(CA.w) <- CA.pa.inds
TY.w <- as.matrix(read.table("TY_admx.2.Q",stringsAsFactors=FALSE))
	row.names(TY.w) <- TY.pa.inds
save(CA.w,TY.w,file="admix.props.Robj")



################################
#	visualize results by stream
################################

delineate.streams <- function(w,names,plotNames){
#	recover()
	pop.breaks <- lapply(names,function(x){range(which(grepl(x,row.names(w))))})
	lapply(1:length(pop.breaks),
		function(i){
			x1 <- pop.breaks[[i]][1]-1
			x2 <- pop.breaks[[i]][2]
			abline(v=x1,col=1,lwd=5)
			abline(v=x1,col="goldenrod1",lwd=4)
			abline(v=x2,col=1,lwd=5)
			abline(v=x2,col="goldenrod1",lwd=4)
			par(xpd=NA)
			text(x=mean(c(x1,x2)),y=1.1,labels=plotNames[i],font=2,cex=1.5)
			par(xpd=TRUE)
		})
	return(invisible("delineated"))
}


# note - this section makes use of plotting functions from the 
#	R pkg conStruct

pdf(file="../6_figs/par_admx_admixture.pdf",width=10,height=8)
	par(mfrow=c(2,1))
	conStruct::make.structure.plot(CA.w)
		delineate.streams(CA.w,c("NCA","PCA","SGS"),c("pre_CA","post_CA","mainstem"))
	conStruct::make.structure.plot(TY.w)
		delineate.streams(TY.w,c("NTY","PTY","SGS"),c("pre_TY","post_TY","mainstem"))
dev.off()