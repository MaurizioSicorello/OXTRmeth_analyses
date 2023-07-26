# The DNA methylation landscape of the human oxytocin receptor gene (OXTR): data-driven clusters and their relation to gene expression and childhood adversity
by Dr. Svenja Müller, M. Sc. Maurizio Sicorello, Dr. Dirk Moser, M. Sc. Leonard Frach,
B. Sc. Alicia Limberg, Dr. Anja M. Gumpp, Dr. Laura Ramo-Fernandez, Dr. Franziska Köhler-Dauner, Prof. Dr. Jörg M. Fegert, Prof. Dr. Christiane Waller, Prof. Dr. Robert Kumsta, Prof. Dr. Iris-Tatjana Kolassa

https://doi.org/10.1038/s41398-023-02548-6

## Open Code
Permission for openly sharing the data was not granted. Please contact the authors for access requests.

Folders: <br/>
-Figures (output figures of analyses) <br/>
-Functions (helper functions called from the analysis scripts) <br/>
-OXTR-shinyapp (shiny web application, available under: https://msicorello.shinyapps.io/oxtr-shinyapp/) <br/>
-renv (snapshot of the R package environment to facilitate reproducibility. See details below) <br/>
-Results (results from the main analyses) <br/>
----DMPs_and_DMRs (results from CpG-wise analyses) <br/>
----multimedBootResults_cat (bootstrapping results of multivariate mediation with categorical trauma IV) <br/>
----multimedBootResults_cont (bootstrapping results of multivariate mediation with continuous trauma IV) <br/>
----multivarMediation_loadings (CpG loadings for both multivaritate mediation models) <br/>
----OXTR_table S1_CpG numbering (supplementary table) <br/>
----plsPerm_all_CTQ/CTQcat/geneExpr (PLS permutations for the three outcomes) <br/>
----shinyPrepareOXTRstructure (results table for shinyapp) <br/>

Most analyses can be performed by opening the R-Project file "OXTRmeth_analyses" and then opening OXTRstructure.R.
Before running this script, make sure the package renv is installed using install.packages("renv"). Then, use
renv::restore() to restore our package versions. Packages should be installed automatically when you run the code.
Similarly, paths will be set in a cross-platform compatible manner, relative to the R-project file,
using the package "here". Results were last checked with the following session info: 

R version 4.1.2 (2021-11-01) <br/>
Platform: x86_64-w64-mingw32/x64 (64-bit) <br/>
Running under: Windows 10 x64 (build 22000) <br/>

The script "singleCpgs_and_DMRs.R" can be used to check the CpG-wise analyses.

The script "MediationAnalysis.m" contains the mediation analyses, performed in matlab (v2020a).
