

#################################
# Load and prepare packages/data

if (!require("pacman")) install.packages("pacman")
pacman::p_load(corrplot, mclust, party, psych, fpc, dbscan, stringr, ggplot2, EnvStats, stats, factoextra, reshape, caret, pls, plyr, here)


# load data and convert % to decimal
df = read.csv(here("..", "Data", "RUB_OXTR_Daten_26.4.csv"), header = T, sep = ";")

# subset mother CpGs and convert to decimal
df_CpG_m = df[,grepl("CpG_m", names(df))]
df_CpG_m = as.data.frame(apply(df_CpG_m,2, function(x){as.numeric(sub("%", "", x, fixed=TRUE))/100}))
df_outcomes = df[complete.cases(df_CpG_m), c("ctq_dich01_multi01", "ctq_sum", "Genexp_OXTR_mother", "t0_age_mo")]
df_CpG_m = df_CpG_m[complete.cases(df_CpG_m),]

# load identifier of gene sections
dfcodes = read.csv(here("..", "Data", "OXTR_segmentCodes.csv"), header = T, sep = ";")
dfcodes = dfcodes[dfcodes$CpG %in% str_replace(names(df_CpG_m), "_m_", "_"),]

# read Cpg glossar for chromosomal position
glossar <- read.csv(here("..", "Data", "CpG_glossar.csv"))

#################################
# qualitative description of OXTR

# plot methylation across CpGs to check for sufficient variance
boxplot(df_CpG_m, range = 1.5)

# count/plot number of outliers per CpG (3*IQR)
out_count = as.data.frame(matrix(nrow = ncol(df_CpG_m), ncol = 2, dimnames = list(NULL, c("CpG", "nOut"))))
out_count$CpG = str_replace(names(df_CpG_m), "_m_", "_")
for(i in 1:ncol(df_CpG_m)){out_count[i,2] = length(boxplot.stats(df_CpG_m[,i], coef = 3)$out)}
df_out = merge(out_count, dfcodes, by = "CpG", sort = F)
ggplot(data = df_out, aes(x = reorder(CpG, 1:nrow(df_out)), y = nOut, colour = segment2)) + 
  geom_point()

df_CpGwiseSDs <- as.numeric(read.csv2(here("..", "Data", "OXTRStabwSingleCpGs.csv"))[,3])^2
df_CpGwiseSDs <- c(df_CpGwiseSDs, rep(0.0228^2, times = ncol(df_CpG_m)-length(df_CpGwiseSDs)))
dfchi <- nrow(df_CpG_m)-1
CpGvariancePs <- numeric(ncol(df_CpG_m))

for(i in 1:ncol(df_CpG_m)){
  
  chi_square <- dfchi*var(df_CpG_m[,i])/df_CpGwiseSDs[i]
  CpGvariancePs[i] = pchisq(chi_square, df = dfchi, lower.tail = F)
  
}

pf(2.5, 109, 11, lower.tail = F)

sum(CpGvariancePs > 0.05)
sum(p.adjust(CpGvariancePs, method = "fdr") > 0.05)
dfvar_plot = data.frame(dfcodes, CpGvariancePs)
dfvar_plot$insuffVar = ifelse(CpGvariancePs > 0.05, 1, 0)

ggplot(data = dfvar_plot, aes(x = reorder(CpG, 1:nrow(df_out)), y = insuffVar, colour = segment2)) +
  geom_point()

# continue without CpGs with insufficient variance
# df_CpG_m <- df_CpG_m[, CpGvariancePs <= 0.05]
# dfcodes <- dfcodes[CpGvariancePs <= 0.05, ]
# 
# # filter glossar
# glossar <- glossar[glossar$CpG_Nr %in% dfcodes$CpG, ]

# matrix with coverage values
coverage <- t(glossar$coverage_perCent_CpGs)
coverage <- data.frame(coverage[rep(seq_len(nrow(coverage)), nrow(df_CpG_m)), ])

names(df_outcomes) <- c("group", "ctq", "genExpr", "age_mother")

# make factor out of CTQ variable (0 = control group, 1 = early adversity? oder andersherum?)
df_outcomes$group <- factor(df_outcomes$group, levels = 0:1, labels = c("CG", "EA"))



###################################

## ANALYSES - Identify differentially-methylated positions (DMPs) between groups using limma (other methods can be used)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("limma"))
  BiocManager::install("limma")

library(limma)

# matrix with methylation values
methylation <- t(df_CpG_m)

design <- model.matrix(~0 + group, data = df_outcomes)
colnames(design) <- c(levels(df_outcomes$group))

# fit the linear model for methylation values
(fit.m <- lmFit(methylation, design))

# create a contrast matrix for specific comparisons
(contMatrix.m = makeContrasts(CGvsEA = CG - EA,levels = design))

# fit the contrasts
fit.m2 <- contrasts.fit(fit.m, contMatrix.m)
(fit.m2 <- eBayes(fit.m2))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.m2))

DMPs <- topTable(fit.m2,  num = Inf)
head(DMPs, n = 10)

# plot the top 10 most significantly differentially methylated CpGs
png(filename = "Figures/DMPs_group_top10.png", height = 3000, width = 2500, res = 400)
par(mfrow = c(4,3))
sapply(rownames(DMPs)[1:10], function(cpg){minfi::plotCpg(methylation, cpg = cpg, pheno = df_outcomes$group, ylab = "beta values")})
dev.off()

#repeat with covariates (age))

# matrix with methylation values

design_age <- model.matrix(~1 + group + age_mother, data = df_outcomes)

# fit the linear model for methylation values
fit.age <- lmFit(methylation, design_age)
(fit.age2 <- eBayes(fit.age))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.age2))

DMPs.age <- topTable(fit.age2,  num = Inf, coef = "groupEA")
head(DMPs.age, n = 10)

DMPs.age[, c(1, 3)] <- -DMPs.age[, c(1, 3)]

##############################

## Differentially methylated regions (DMRs) using DMRcate - one of the most commonly used methods for DMR analysis, which comes with some options for sequencing data

options(connectionObserver = NULL)

if (!requireNamespace("DMRcate"))
  BiocManager::install("DMRcate")
if (!requireNamespace("bsseq"))
  BiocManager::install("bsseq")

library(DMRcate)
library(bsseq)


# exclude individuals with IDs B0031, B0132 and B0163 because of many missing values
#data_excl <- subset(, subset = data_mother$Code != "B0163" & data_mother$Code != "B0132" & data_mother$Code != "B0031")

# get chromosomal position by merging with glossar
m <- data.frame(methylation)
m$CpG_Nr <- row.names(m)
m$CpG_Nr <- gsub('m_', '', m$CpG_Nr)

methylation <- dplyr::left_join(m, glossar)
row.names(methylation) <- methylation$Chromosomal_Location
position <- row.names(methylation)

methylation <- as.matrix(methylation[, 1:110])

position <- stringr::str_split_fixed(position, ":", 2)

chr <- position[,1]
pos <- as.numeric(position[,2])
rownames(methylation) <- paste(chr, pos, sep = ":")

design_excl <- model.matrix(~0 + group, data = df_outcomes)

methdesign <- edgeR::modelMatrixMeth(design_excl)
methcont <- makeContrasts(CGvsEA = groupCG - groupEA, levels = methdesign)

coverage <- t(coverage[c(1:nrow(df_outcomes)), ])

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation, Cov = coverage, pos = pos, chr = chr, sampleNames = rownames(df_outcomes))

# test cpgs
myAnnotation_group <- sequencing.annotate(bseq_obj, methdesign = methdesign, all.cov = T, contrasts = T, cont.matrix = methcont, fdr = 0.05, coef = "CGvsEA")
str(myAnnotation_group)

# look for DMRs within 500bp
DMRs_group <- dmrcate(myAnnotation_group, lambda = 500)

# results table
(results.ranges <- extractRanges(DMRs_group))

# get relevant estimates 
results.stats <- data.frame(myAnnotation_group@ranges@elementMetadata@listData)
results.cpgs <- data.frame(myAnnotation_group@ranges@ranges@start)

# individual significant cpgs (after FDR correction)
(singleCpgs <- results.cpgs[results.stats$is.sig == T, ])

# get positions of cpgs within the DMR
dmrOx <- pos[pos >= results.ranges@ranges@start & pos <= results.ranges@ranges@start + results.ranges@ranges@width - 1]

# add chromosome to cpg position so set as name
dmrOx <- paste("chr3", dmrOx, sep = ".")

# get all cpgs within the ranges
m <- data.frame(t(methylation))

dmr <- m[names(m) %in% dmrOx]

# add variables group and total ctq score
dmr <- cbind(df_outcomes$group, df_outcomes$ctq, dmr)
test <- data.frame(dmr[, 2:11])

# get a correlation matrix of all cpgs within the ranges and ctq sum score
psych::corr.test(test)

# correlation table output

#corDMR <- corr.test(test)
#apaTables::apa.cor.table(corDMR, filename = "outputs/correlations_DMR_with_ctq.doc")

dmr$meanDMR <- rowMeans(dmr[3:11], na.rm = T)
names(dmr)[1:2] <- c("group", "ctq")

with(dmr, plot(meanDMR ~ ctq))

# group comparisons
with(dmr, plot(meanDMR ~ group))

with(dmr, car::leveneTest(meanDMR, group))
with(dmr, nortest::lillie.test(meanDMR))

with(dmr, t.test(meanDMR ~ group, var.equal = T))

# pval for mean DMR methylation is significant but potential problem of a mean value as correlations among (many) cpgs vary largely


## repeat DMR analysis controlling for age

design_excl_age <- model.matrix(~1 + group + age_mother, data = df_outcomes)
methdesign_age <- edgeR::modelMatrixMeth(design_excl_age)

bseq_obj <- BSseq(M = methylation, Cov = coverage, pos = pos, chr = chr, sampleNames = rownames(df_outcomes))

myAnnotation_group_age <- sequencing.annotate(bseq_obj, methdesign = methdesign_age, all.cov = T, contrasts = F, fdr = 0.10, coef = "groupEA")
str(myAnnotation_group_age)
DMRs_group_age <- dmrcate(myAnnotation_group_age, lambda = 500)
DMRs_group_age@meandiff <- -DMRs_group_age@meandiff 
DMRs_group_age@maxdiff <- -DMRs_group_age@maxdiff 

# no individual significant cpgs after adjusting for multiple testing if FDR = 0.05
# for FDR = 0.10 see results below (13 individual cpgs that have FDR < 0.05)

# results table
(results.ranges_age <- extractRanges(DMRs_group_age)) # 13 cpgs included in DMR 8809133-8809238 

results.stats <- data.frame(myAnnotation_group_age@ranges@elementMetadata@listData)
results.cpgs <- data.frame(myAnnotation_group_age@ranges@ranges@start)
 
(singleCpgs <- results.cpgs[results.stats$is.sig == T, ])
 
dmrOx_age <- pos[pos >= results.ranges@ranges@start & pos <= results.ranges_age@ranges@start + results.ranges_age@ranges@width - 1]
 
# add chromosome to cpg position so set as name
dmrOx_age <- paste("chr3", dmrOx_age, sep = ".")
 
# get all cpgs within the ranges
dmrOx_age <- m[names(m) %in% dmrOx_age]
 
# add variables group and total ctq score
dmr_age <- cbind(df_outcomes$group, df_outcomes$ctq, dmrOx_age)
 
dmr_age$meanDMR <- rowMeans(dmr_age[3:13], na.rm = T)
names(dmr_age)[1:2] <- c("group", "ctq")
 
with(dmr_age, t.test(meanDMR ~ group, var.equal = T))


## Sensitivity analysis - repeat DMR analysis using different base pair cut-offs defining a region

# look for DMRs within 500bp
DMRs_group250 <- dmrcate(myAnnotation_group, lambda = 250)
(results.ranges250 <- extractRanges(DMRs_group250)) # only 3 cpgs within 8809179-8809225   

DMRs_group750 <- dmrcate(myAnnotation_group, lambda = 750)
(results.ranges750 <- extractRanges(DMRs_group750)) # 16 cpgs within 8809055-8809160


##############################

## Same analyses for gene expression instead of group (here statistically as predictor, although conceptually gene expression is the outcome)

# single CpGs

# subset of mothers with gene expression data
df_expr <- df_outcomes[complete.cases(df_outcomes$genExpr), ]
df_CpG_m_expr <- df_CpG_m[rownames(df_CpG_m) %in% rownames(df_expr), ]

methylation2 <- t(df_CpG_m_expr)

design_expr <- model.matrix(~1 + genExpr, data = df_expr)

# fit the linear model for methylation values
fit.expr <- lmFit(methylation2, design_expr)
(fit.expr2 <- eBayes(fit.expr))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.expr2)) # no individually significant SNPs (low power with N = 69)

DMPs_geneExpr <- topTable(fit.expr2,  num = Inf)
head(DMPs_geneExpr, n = 20)


## DMRs ##

# get chromosomal position by merging with glossar
m <- data.frame(methylation2)
m$CpG_Nr <- row.names(m)
m$CpG_Nr <- gsub('m_', '', m$CpG_Nr)

methylation2 <- dplyr::left_join(m, glossar)
row.names(methylation2) <- methylation2$Chromosomal_Location
position <- row.names(methylation2)

methylation2 <- as.matrix(methylation2[, 1:69])

position <- stringr::str_split_fixed(position, ":", 2)

chr <- position[,1]
pos <- as.numeric(position[,2])
rownames(methylation2) <- paste(chr, pos, sep = ":")

design_excl <- model.matrix(~1 + genExpr, data = df_expr)

methdesign <- edgeR::modelMatrixMeth(design_excl)
coverage2 <- coverage[, c(1:nrow(df_expr))]

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation2, Cov = coverage2, pos = pos, chr = chr, sampleNames = rownames(df_expr))

# test cpgs
myAnnotation_geneExpr <- sequencing.annotate(bseq_obj, methdesign = methdesign, all.cov = T, contrasts = F, fdr = 0.05, coef = "genExpr")
str(myAnnotation_geneExpr)

# look for DMRs within 500bp
DMRs_geneExpr <- dmrcate(myAnnotation_geneExpr, lambda = 500) # no individual CpG that are significant, so DMRcate does not work

# only significant CpGs at FDR 0.30 (NOT RELIABLE AT ALL, RESULTS ARE NOT MEANINGFUL!!)

DMRs_geneExprFDR30 <- dmrcate(changeFDR(myAnnotation_geneExpr,FDR = 0.3), lambda = 500)

(results.ranges_geneExprFDR30 <- extractRanges(DMRs_geneExprFDR30)) # including 17 CpGs, but again - NOT MEANINGFUL 



##############################

# repeat main analyses (single CpGs) using ctq sum score instead of group variable

# matrix with methylation values
methylation <- t(df_CpG_m)

design_ctq <- model.matrix(~1 + ctq, data = df_outcomes)

# fit the linear model for methylation values
(fit.ctq <- lmFit(methylation, design_ctq))
(fit.ctq2 = eBayes(fit.ctq))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.ctq2))

DMPs.ctq <- topTable(fit.ctq2,  num = Inf, coef = "ctq")
head(DMPs.ctq, n = 10)

# results differ compared to the analysis of group differences 


#repeat with covariates (age))

# matrix with methylation values

design_age_ctq <- model.matrix(~1 + ctq + age_mother, data = df_outcomes)

# fit the linear model for methylation values
fit.age.ctq <- lmFit(methylation, design_age_ctq)
(fit.age2.ctq <- eBayes(fit.age.ctq))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.age2.ctq))

DMPs.age.ctq <- topTable(fit.age2.ctq,  num = Inf, coef = "ctq")
head(DMPs.age.ctq, n = 10)


# DMR analysis

m <- data.frame(methylation)
m$CpG_Nr <- row.names(m)
m$CpG_Nr <- gsub('m_', '', m$CpG_Nr)

methylation <- dplyr::left_join(m, glossar)
row.names(methylation) <- methylation$Chromosomal_Location
position <- row.names(methylation)

methylation <- as.matrix(methylation[, 1:110])

position <- stringr::str_split_fixed(position, ":", 2)

chr <- position[,1]
pos <- as.numeric(position[,2])
rownames(methylation) <- paste(chr, pos, sep = ":")

methdesign_ctq <- edgeR::modelMatrixMeth(design_ctq)

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation, Cov = coverage, pos = pos, chr = chr, sampleNames = rownames(df_outcomes))


# test cpgs
myAnnotation_ctq <- sequencing.annotate(bseq_obj, methdesign = methdesign_ctq, all.cov = T, contrasts = F, fdr = 0.05, coef = "ctq")
str(myAnnotation_ctq)

# look for DMRs within 500bp
DMRs_ctq <- dmrcate(myAnnotation_ctq, lambda = 500)

(results.ranges.ctq <- extractRanges(DMRs_ctq))

results.stats.ctq <- data.frame(myAnnotation_ctq@ranges@elementMetadata@listData)
results.cpgs.ctq <- data.frame(myAnnotation_ctq@ranges@ranges@start)

(singleCpgs_ctq <- results.cpgs.ctq[results.stats.ctq$is.sig == T, ])

dmrOxCtq <- pos[pos >= results.ranges.ctq@ranges@start & pos <= results.ranges.ctq@ranges@start + results.ranges.ctq@ranges@width - 1]

# add chromosome to cpg position so set as name
dmrOxCtq <- paste("chr3", dmrOxCtq, sep = ".")

# get all cpgs within the ranges
m <- data.frame(t(methylation))

dmrCtq <- m[names(m) %in% dmrOxCtq]


#############################

## check if results of single Cpgs (DMPs) are similar when using t.tests/linear models instead of limma


# group differences
cpgsLev <- lapply(df_CpG_m[,], function(x) car::leveneTest(x ~ df_outcomes$group)$`Pr(>F)`)
sort(unlist(cpgsLev)) # no equal variances for 7 cpgs

cpgsTtest <- lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group))

## get parameters for summary table

# differences in mean methylation
diff <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$estimate[2])) - 
  unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$estimate[1]))

##########


# differences in mean methylation
diff <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$estimate[1])) - 
  unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$estimate[2]))

# 95% confidence intervals
CI_lower <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$conf.int[1]))
CI_upper <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$conf.int[2]))

# t-statistic
t <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$statistic[1]))

# degrees of freedom
df <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$parameter))

# standard error
SE <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$stderr))

# p-value
p <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$p.value))


# create summary table
results.ttest <- data.frame(diff, CI_lower, CI_upper, t, df, SE, p)
row.names(results.ttest) <- names(df_CpG_m)

results.ttest <- results.ttest[order(results.ttest$p),]
results.ttest$FDR <- p.adjust(results.ttest$p, method = "BH", n = dim(results.ttest)[1])

# calculate pooled SD
SD <- SE  * sqrt(dim(df_outcomes)[1]) 

# calculate Cohen's D
d <- diff / SD

## calculate more robust Hedges' g

# empty vector
g <- rep(NA, times = dim(df_CpG_m)[2])

# loop to calculate g for every cpg site
for (i in 1:dim(df_CpG_m)[2]) {
  g[i] <- effsize::cohen.d(df_CpG_m[, i], 
                           df_outcomes$group, 
                           df_CpG_m[, i] ~ df_outcomes$group, 
                           pooled = T, 
                           paired = F, 
                           na.rm = F, 
                           hedges.correction = T)$estimate
}

# correlate d and g
cor(d, g) 

# almost perfect correlation 0.99993

# reverse sign
results.ttest$Hedges_g <- g

g_ci_lower <- rep(NA, times = dim(df_CpG_m)[2]) 

# loop to calculate g for every cpg site
for (i in 1:dim(df_CpG_m)[2]) {
  g_ci_lower[i] <- effsize::cohen.d(df_CpG_m[, i], 
                           df_outcomes$group, 
                           df_CpG_m[, i] ~ df_outcomes$group, 
                           pooled = T, 
                           paired = F, 
                           na.rm = F, 
                           hedges.correction = T)$conf.int[1]
}

g_ci_upper <- rep(NA, times = dim(df_CpG_m)[2]) 

# loop to calculate g for every cpg site
for (i in 1:dim(df_CpG_m)[2]) {
  g_ci_upper[i] <- effsize::cohen.d(df_CpG_m[, i], 
                                    df_outcomes$group, 
                                    df_CpG_m[, i] ~ df_outcomes$group, 
                                    pooled = T, 
                                    paired = F, 
                                    na.rm = F, 
                                    hedges.correction = T)$conf.int[2]
}

results.ttest$CI_g_lower <- g_ci_lower
results.ttest$CI_g_upper <- g_ci_upper


names(results.ttest)[7] <- "Pval"

summary(results.ttest$Hedges_g)
summary(results.ttest$diff)


# linear models using CTQ to predict DNAm


cpgsLm <- lapply(df_CpG_m[,], function(x) lm(df_outcomes$ctq ~ x))
results.lm.ctq <- t(data.frame(lapply(df_CpG_m[,], function(x) summary(lm(scale(x) ~ scale(df_outcomes$ctq)))$coefficients[2, ])))
results.lm.ctq <- data.frame(results.lm.ctq)

results.lm.ctq <- results.lm.ctq[order(results.lm.ctq$Pr...t..),]
results.lm.ctq$FDR <- p.adjust(results.lm.ctq$Pr...t.., method = "BH", n = dim(results.lm.ctq)[1])

names(results.lm.ctq)[c(2,4)] <- c("SE", "Pval")

results.lm.ctq$CI_lower <- with(results.lm.ctq, Estimate - 1.96 * SE)
results.lm.ctq$CI_upper <- with(results.lm.ctq, Estimate + 1.96 * SE)


results.lm.ctq <- results.lm.ctq[, c(1, 6:7, 2:5)]

summary(results.lm.ctq$Estimate)

# linear models predicting gene expression

cpgsLm <- lapply(df_CpG_m[,], function(x) lm(x ~ df_CpG_m))
results.lm <- t(data.frame(lapply(df_CpG_m_expr[,], function(x) summary(lm(scale(df_expr$genExpr) ~ scale(x)))$coefficients[2, ])))
results.lm <- data.frame(results.lm)

results.lm <- results.lm[order(results.lm$Pr...t..),]
results.lm$FDR <- p.adjust(results.lm$Pr...t.., method = "BH", n = dim(results.lm)[1])


results.lm$CI_lower <- with(results.lm, Estimate - 1.96 * SE)
results.lm$CI_upper <- with(results.lm, Estimate + 1.96 * SE)


results.lm <- results.lm[, c(1, 6:7, 2:5)]

##############################


# write outputs in excel sheets
require(openxlsx)
list_results <- list("DMPs_Group" = DMPs,
                     "DMPs_ctq" = DMPs.ctq, 
                     "DMPs_geneExpression" = DMPs_geneExpr[, 1:6],
                     "DMPs_Group_ageControlled" = DMPs.age,
                     "DMPs_CTQ_ageControlled" = DMPs.age.ctq,
                     "T-tests_Group" = results.ttest,
                     "LinearModels_CTQ" = results.lm.ctq,
                     "LinearModels_GeneExpr" = results.lm,
                     "DMR_group" = results.ranges,
                     "DMR_group250bp" = results.ranges250,
                     "DMR_group750bp" = results.ranges750,
                     "DMR_group_ageControlled" = results.ranges_age,
                     "DMR_ctq" = results.ranges.ctq)
write.xlsx(list_results, file = "Results/DMPs_and_DMRs.xlsx", colNames = T, rowNames = T)

##############################

## PLOTS

# single CpGs "manhattan plot"

pLog <- -log10(DMPs$P.Value)
cg_names <- row.names(DMPs)
oxtrPlot <- data.frame(cg_names, pLog)

dmpsG <- DMPs_geneExpr

dmpsG$pLog_genExpr <- -log10(dmpsG$P.Value)
dmpsG$cg_names <- row.names(dmpsG)
oxtrPlot <- merge(oxtrPlot, dmpsG)[c(1,2,9)]

names(oxtrPlot) <- c("CpG", "log10_P_Group", "log10_P_GeneExpr")

# merge with 'glossar' data frame to get different names for the cpgs
oxtrPlot$CpG <- gsub('m_', '', oxtrPlot$CpG)
oxtrPlot <- dplyr::left_join(oxtrPlot, glossar, by = c("CpG" = "CpG_Nr"))
oxtrPlot <- dplyr::left_join(oxtrPlot, dfcodes)

# split to get the bp position
position <- str_split_fixed(oxtrPlot$Chromosomal_Location, ":", 2)
oxtrPlot$pos <- as.numeric(position[,2])

oxtrPlot <- oxtrPlot[with(oxtrPlot, order(pos, decreasing = T)), ]

dmr_pos <- which(oxtrPlot$pos >= results.ranges@ranges@start & oxtrPlot$pos <= results.ranges@ranges@start + results.ranges@ranges@width - 1)

oxtrPlot$pos <- factor(oxtrPlot$pos, levels = unique(oxtrPlot$pos))

oxtrPlot$CpG <- factor(oxtrPlot$CpG, labels = 1:length(oxtrPlot$CpG))

png(filename = "Figures/OXTR_singleCpgs_Group_and_genExpr.png", width = 2000, height = 1500, type = "cairo", res = 300)
ggplot() + 
  geom_point(data = oxtrPlot[-dmr_pos, ], aes(oxtrPlot$pos[-dmr_pos], oxtrPlot$log10_P_Group[-dmr_pos], colour = factor(oxtrPlot$segment2[-dmr_pos], levels = unique(segment2))), shape = 2) +
  geom_point(data = oxtrPlot[dmr_pos, ], aes(oxtrPlot$pos[dmr_pos], oxtrPlot$log10_P_Group[dmr_pos], colour = factor(oxtrPlot$segment2[dmr_pos], levels = unique(segment2))), shape = 15) +
  geom_point(data = oxtrPlot, aes(pos, log10_P_GeneExpr, colour = factor(segment2, levels = unique(segment2))), shape = 1) +
  geom_hline(yintercept = (-log10(0.05)), linetype = 2) + # nominal p < 0.05
  ylim(0, 5) + labs(x = "OXTR CpG Site", y = "-log10 P-Value") + 
  theme_classic() +
  theme(legend.position = "right") +
  labs(colour = '') +
  scale_x_discrete(limits = rev, labels = NULL) +
  scale_color_brewer(palette = "Spectral")
dev.off()


######

# same plot but with pvalues from the ctq models instead of group variable

pLog <- -log10(DMPs.ctq$P.Value)
cg_names <- row.names(DMPs.ctq)
oxtrPlot <- data.frame(cg_names, pLog)

dmpsG <- DMPs_geneExpr

dmpsG$pLog_genExpr <- -log10(dmpsG$P.Value)
dmpsG$cg_names <- row.names(dmpsG)
oxtrPlot <- merge(oxtrPlot, dmpsG)[c(1,2,9)]

names(oxtrPlot) <- c("CpG", "log10_P_CTQ", "log10_P_GeneExpr")

# merge with 'glossar' data frame to get different names for the cpgs
oxtrPlot$CpG <- gsub('m_', '', oxtrPlot$CpG)
oxtrPlot <- dplyr::left_join(oxtrPlot, glossar, by = c("CpG" = "CpG_Nr"))
oxtrPlot <- dplyr::left_join(oxtrPlot, dfcodes)

# split to get the bp position
position <- str_split_fixed(oxtrPlot$Chromosomal_Location, ":", 2)
oxtrPlot$pos <- as.numeric(position[,2])

oxtrPlot <- oxtrPlot[with(oxtrPlot, order(pos, decreasing = T)), ]

dmr_pos <- which(oxtrPlot$pos >= results.ranges.ctq@ranges@start & oxtrPlot$pos <= results.ranges.ctq@ranges@start + results.ranges.ctq@ranges@width - 1)

oxtrPlot$pos <- factor(oxtrPlot$pos, levels = unique(oxtrPlot$pos))

oxtrPlot$CpG <- factor(oxtrPlot$CpG, labels = 1:length(oxtrPlot$CpG))

png(filename = "Figures/OXTR_singleCpgs_CTQ_and_genExpr.png", width = 2000, height = 1500, type = "cairo", res = 300)
ggplot() + 
  geom_point(data = oxtrPlot[-dmr_pos, ], aes(oxtrPlot$pos[-dmr_pos], oxtrPlot$log10_P_CTQ[-dmr_pos], colour = factor(oxtrPlot$segment2[-dmr_pos], levels = unique(segment2))), shape = 2) +
  geom_point(data = oxtrPlot[dmr_pos, ], aes(oxtrPlot$pos[dmr_pos], oxtrPlot$log10_P_CTQ[dmr_pos], colour = factor(oxtrPlot$segment2[dmr_pos], levels = unique(segment2))), shape = 15) +
  geom_point(data = oxtrPlot, aes(pos, log10_P_GeneExpr, colour = factor(segment2, levels = unique(segment2))), shape = 1) +
  geom_hline(yintercept = (-log10(0.05)), linetype = 2) + # nominal p < 0.05
  ylim(0, 7) + labs(x = "OXTR CpG Site", y = "-log10 P-Value") + 
  theme_classic() +
  theme(legend.position = "right") +
  labs(colour = '') +
  scale_x_discrete(limits = rev, labels = NULL) +
  scale_color_brewer(palette = "Spectral")
dev.off()

sessionInfo()

