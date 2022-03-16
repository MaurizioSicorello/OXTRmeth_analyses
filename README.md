# OXTRmeth_analyses
Open code for the paper entitled: 
"The DNA methylation landscape of the human oxytocin receptor gene (OXTR): Recommendations for future research"
by M.Sc. Svenja Müller, M. Sc. Maurizio Sicorello, Dr. Dirk Moser, M. Sc. Leonard Frach, 
B. Sc. Alicia Limberg, Dr. Anja M. Gumpp, Dr. Laura Ramo-Fernandez, Prof. Dr. Robert Kumsta, 
Prof. Dr. Iris -Tatjana Kolassa

Permission for openly sharing the data was not granted. Please contact the authors for access requests.

Folders:
-Figures (output figures of analyses)
-Functions (helper functions called from the analysis scripts)
-OXTR-shinyapp (shiny web application, available under: https://msicorello.shinyapps.io/oxtr-shinyapp/)
-renv (snapshot of the R package environment to facilitate reproducibility. See details below)
-Results (results from the main analyses)
----DMPs_and_DMRs (results from CpG-wise analyses)
----multimedBootResults_cat (bootstrapping results of multivariate mediation with categorical trauma IV)
----multimedBootResults_cont (bootstrapping results of multivariate mediation with continuous trauma IV)
----multivarMediation_loadings (CpG loadings for both multivaritate mediation models)
----OXTR_table S1_CpG numbering (supplementary table)
----plsPerm_all_CTQ/CTQcat/geneExpr (PLS permutations for the three outcomes)
----shinyPrepareOXTRstructure (results table for shinyapp)

Most analyses can be performed by opening the R-Project file "OXTRmeth_analyses" and then opening OXTRstructure.R.
Before running this script, make sure the package renv is installed using install.packages("renv"). Then, use
renv::restore() to restore our package versions. Packages should be installed automatically when you run the code.
Similarly, paths will be set in a cross-platform compatible manner, relative to the R-project file, 
using the package "here". Results were last checked with the following session info:

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22000)

The script "singleCpgs_and_DMRs.R" can be used to check the CpG-wise analyses.

The script "MediationAnalysis.m" contains the mediation analyses, performed in matlab (v2020a).