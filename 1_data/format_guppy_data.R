################################################################
################################################################
#
#  script for formatting guppy data for ancestry analysis
#
################################################################
################################################################

################################
# population codes:
#	SGS - mainstem ("source" of immigrants to headwater translocation sites)
#	NCA - post gene flow headwater Caigual
#	NTY - pre gene flow headwater Taylor
#	PCA - post gene flow headwater Caigual
#	PTY - post gene flow headwater Taylor
################################

stream.names <- c("NCA","NTY","PCA","PTY","SGS")

################################
# export data in plink binary format using vcfTools and plink for each stream pop
################################
plink.bin.export <- function(prefix){
	call <- sprintf("vcftools --vcf nchr_%s_mapped.vcf --plink --out %s_mapped_counts",prefix,prefix)
		system(call)
	call <- sprintf("plink --noweb --file %s_mapped_counts --make-bed --out %s",prefix,prefix)
		system(call)
	file.remove(sprintf("%s_mapped_counts.map",prefix))
	file.remove(sprintf("%s_mapped_counts.ped",prefix))
	file.remove(sprintf("%s.nosex",prefix))
	return(invisible("exported"))
}

lapply(stream.names,
		function(x){
			plink.bin.export(x)
		})

################################
# create datasets for ADMIXTURE
#	consisting of: 
#		Caigual drainage: {NCA,SGS,PCA}
#		Taylor drainage: {NTY,SGS,PTY}
################################

# make plink files containing pre-contact parental populations and 
#	each of the two post-contact admixed populations

merge.list.ca <- c(c(paste0("SGS",c(".bed",".bim",".fam")),"\n"),
				c(paste0("PCA",c(".bed",".bim",".fam")),"\n"))
merge.list.ty <- c(c(paste0("SGS",c(".bed",".bim",".fam")),"\n"),
				c(paste0("PTY",c(".bed",".bim",".fam")),"\n"))
cat(merge.list.ca,file="merge.list.ca.txt",sep=" ")
cat(merge.list.ty,file="merge.list.ty.txt",sep=" ")

call <- "plink --noweb --bfile NCA --merge-list merge.list.ca.txt --make-bed --out CA_par_admx"
system(call)
call <- "plink --noweb --bfile NTY --merge-list merge.list.ty.txt --make-bed --out TY_par_admx"
system(call)

# export frequency data for each pop from each stream data subset
ca.famfile <- read.table("CA_par_admx.fam",stringsAsFactors=FALSE,header=FALSE)
ca.famfile[,1] <- unlist(lapply(strsplit(ca.famfile[,1],"-"),"[[",1))
write.table(ca.famfile,file="CA_par_admx.fam",quote=FALSE,col.names=FALSE,row.names=FALSE)

ca.famfile[,3] <- ca.famfile[,1]
write.table(ca.famfile[,1:3],file="ca.famfile.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

ty.famfile <- read.table("TY_par_admx.fam",stringsAsFactors=FALSE,header=FALSE)
ty.famfile[,1] <- unlist(lapply(strsplit(ty.famfile[,1],"-"),"[[",1))
write.table(ty.famfile,file="TY_par_admx.fam",quote=FALSE,col.names=FALSE,row.names=FALSE)

ty.famfile[,3] <- ty.famfile[,1]
write.table(ty.famfile[,1:3],file="ty.famfile.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


call <- "plink --noweb --bfile CA_par_admx --within ca.famfile.txt --freq --out CA"
system(call)
call <- "plink --noweb --bfile TY_par_admx --within ty.famfile.txt --freq --out TY"
system(call)


################################
# clean up intermediate files
################################

clean.up <- function(stream){
	burn.files <- paste0(stream,c(".bed",".bim",".fam",".log",".nosex"))
	file.remove(burn.files)
}

lapply(stream.names,clean.up)
file.remove(
	c("ca.famfile.txt",
	  "ty.famfile.txt",
	  "merge.list.ca.txt",
	  "merge.list.ty.txt",
	  "TY.nosex",
	  "TY_par_admx.nosex",
	  "CA.nosex",
	  "CA_par_admx.nosex")
)