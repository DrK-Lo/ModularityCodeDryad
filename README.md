README

This repository contains all the data and scripts needed to reproduce the
co-association network analysis for the paper: "Modularity of genes involved in local adaptation to climate despite physical linkage" published on bioarxiv https://www.biorxiv.org/content/early/2018/01/26/202481 and in review at Genome Biology.

Note that there are some large data files that could not be uploaded to GitHub, and are stored in the Dryad Respository (link to be provided upon acceptance of manuscript for publication). Note also that there to reproduce the results you will need a folder for "results", which is ignored by git.

### Steps to run empirical analyses:


1) open analysis/dryad.Rproj
	- This will set your working directory

2) Run AnalysisPEGA_3_dryad.Rmd
	- Contains all code to reproduce the co-association network analysis

3) Run AnalysisPEGA_4_recombin_dryad.Rmd
	- Creates gametic disequilibrium and recombination rate plots.
	- Note that this must be run after #2

Other scripts:

* galaxy2.R function for galaxy biplots

* myFunDist.R contains functions for making co-association networks from correlation data

* plot2Dcov.R another function for making galaxy biplots


### Descriptions of Data:

* large_files/Pine_Alpha_AveRho_WithSuperLogical_andAncestralFlip.RData
	* dataframe containing all association analyses (only a subset used for this study)
	
* large_files/Pine_Alpha_AveRho_WithSuperLogical_andAncestralFlip.metadata
	* metadata for corresponding file

* large_files/pine_super_outlier_LD_p4_all_genemeans.txt
	* results from Yeaman 2016 Science for the top candidate contigs based on enrichment for associations
	
* large_files/pine_super_outlier_LD_p4_all_genemeans.metadata
	* metadata for the corresponding file
	
* orthologs/parallelism_genes_pine_spruce_q05.txt
	* results from Yeaman 2016 Science for orthologs in pine and spruce that are
	* locally adapting to freezing
	
* orthologs/parallelism_genes_pine_spruce_q05.metadata
	* metadata for the corresponding file
	
* PineHeatmap.Rdata
	* contains correlation matrix for heatmap of environments
	
* PineHeatmap.metadata
	* metadata for the corresponding file
		
* recombinationrates/katie_contigs_recombination.txt
	* recombination rates for SNPs that mapped to our top candidate contigs
	
* recombinationrates/katie_contigs_metadata.txt
	* metadata for the corresponding file

* tair_pine_clusterall
	* GO annotations for top candidate contigs
	
* tair_pine_clusterall.metadata
	* metadata for the corresponding file


### Steps to run analyses of simulations:

* Change directory to the /simulations folder.

* Galaxy_Plots_on_AdaptreeEnvi_1R_aOverlayEnvi.Rmd 
	* runs the code to overlay the Adaptree 
	* environments onto the simulations. This takes a few minutes.

* Galaxy_Plots_on_AdaptreeEnvi_1R_bCalcCorrsV2.Rmd
	* runs the co-association network analysis on the simulations.
	* This will take a while.

NOTE: If you get an error message while installing packages when running the 
	Markdown documents, try installing them directly in the console first.
