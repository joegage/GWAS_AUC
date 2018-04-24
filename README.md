# GWAS_AUC
This repository contains scripts to accompany Gage, de Leon, and Clayton 2018 (_cite_).

* **1_simulate_phenotypes.R**: Simulates phenotypes controlled by 10, 100, or 1,000 causative SNPs selected from the genotypic data described in Mazaheri et al. (2018; _paper doi_), which can be found at doi:10.5061/dryad.n0m260p.
* **simPheno.R**: the workhorse script for 1_simulate_phenotypes.R.  This is a more general phenotype simulation script.  Will generate simulated phenotypes given the desired heritability, number of causative loci, and population size (number of individuals).  The user can either feed in the allele frequencies and genotypes of real markers, or if those arguments are not included the script will simulate them as well.
* **2_prep_dirs_run_GAPIT.sh**: reading from a file with a list of trait names (one per line), this script creates directories for results, runs GWAS in GAPIT (Lipka et al., 2012; doi:10.1093/bioinformatics/bts444).
* **run_GAPIT.R**: workhorse script for 2_prep_dirs_run_GAPIT.sh.  Runs GWAS and writes out results.
* **3_run_pROC.R**: runs ROC_for_pROC.R
* **ROC_for_pROC.R**: creates ROC curves from GAPIT results, tests them against each other, generates averaged ROC curves for traits with the same parameters, and returns everything in formats conducive to plotting.
* **4_plot_ROC_tests.R**: Makes plots used in the paper.
* **5_correlation_figs.R**: Plots correlations between image-based and manual measurements of TL, SL, BN, and TW.
