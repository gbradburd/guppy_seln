################################################################
################################################################
#	master script that runs all steps of analysis
################################################################
################################################################

# requires:
#	R package 'conStruct'
#	vcftools
#	plink
#	ADMIXTURE
#
#	Note that vcftools, plink, and ADMIXTURE 
#	are called via the system() R command, 
#	so you may have to update R's PATH to include 
#	these software packages.

setwd("1_data")
	source("format_guppy_data.R")
setwd("..")

setwd("2_admixture")
	source("run.admixture.R")
setwd("..")

setwd("3_reconcile_alleles")
	source("reconcile.alleles.R")
setwd("..")

setwd("4_silico_sims")
	source("binomial.admx.sims.R")
setwd("..")

setwd("5_test_ancestry_bias")
	source("test_ancestry_bias.R")
setwd("..")

# all output figures will be generated in the "6_figs" directory
