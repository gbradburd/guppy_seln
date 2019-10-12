## README for guppy popgen selection analyses

This repo contains the code/data used to test for the signature 
of selection on maintaining locally adapted variants in the face 
of gene flow in Trinidadian guppies.

These analyses are associated with the manuscript:
Fitzpatrick et al. (2019) Genetic rescue without genomic swamping in wild populations.

 * [manuscript on bioRxiv](https://www.biorxiv.org/content/10.1101/701706v1)
 
## Using this repo

The data and individual R scripts used for the analyses in the paper 
are divided into the different steps in the pipeline.
You can examine each (and the output they produce) within 
the relevant subfolder, or you can run all steps at once using 
the `exe.master.R` script.

The current contents of each folder contain the output 
expected after successfully running the `exe.master.R` script.

## Contents

* 1_data
	* R script for formatting raw data
	* raw vcf data files for each stream
	* plink binary files for running ADMIXTURE
* 2_admixture
	* R script for running ADMIXTURE
	* estimated admixture proportions for each drainage
	* estimated allele frequencies in each ancestral population
* 3_reconcile_alleles
	* R script for ensuring the same alleles are being counted in each dataset
	* allele frequencies in each population
* 4_silico_sims
	* R script for simulating admixed headwater populations
	* R object with saved simulations
* 5_test_ancestry_bias
	* R script for testing for selection on putatively locally adapted alleles
* 6_figs
	* PDF output figures for different stages in the analysis


## Contact

Please direct all queries to bradburd (at) msu.edu, 
or post as issues on the git repo.