

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
df_CpG_m <- df_CpG_m[, CpGvariancePs <= 0.05]
dfcodes <- dfcodes[CpGvariancePs <= 0.05, ]

# filter glossar
glossar <- glossar[glossar$CpG_Nr %in% dfcodes$CpG, ]

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
(contMatrix.m = makeContrasts(CGvsEA = EA - CG,levels = design))

# fit the contrasts
fit.m2 <- contrasts.fit(fit.m, contMatrix.m)
(fit.m2 <- eBayes(fit.m2))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.m2))

DMPs <- topTable(fit.m2,  num = Inf)
head(DMPs, n = 10)

# repeat with m values (recommended)
#methylation[methylation == 0] <- 0.000001

#methylation <- wateRmelon::Beta2M(methylation) 
# values of 0% are mapped to -Inf -> maybe allocate values close to zero more precisely?

# plot the top 10 most signifcantly differentially methylated CpGs
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

##############################

## Differentially methylated regions (DMRs) using DMRcate - one of the most commonly used methods for DMR analysis, which comes with some options for sequencing data

options(connectionObserver = NULL)

if (!requireNamespace("DMRcate"))
  BiocManager::install("DMRcate")
if (!requireNamespace("DSS"))
  BiocManager::install("DSS")
if (!requireNamespace("bsseq"))
  BiocManager::install("bsseq")

library(DMRcate)
library(DSS)
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
methcont <- makeContrasts(CGvsEA = groupEA - groupCG, levels = methdesign)

coverage <- t(coverage[c(1:nrow(df_outcomes)), ])

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation, Cov = coverage, pos = pos, chr = chr, sampleNames = rownames(df_outcomes))

# test cpgs
myAnnotation_group <- sequencing.annotate(bseq_obj, methdesign = methdesign, all.cov = T, contrasts = T, cont.matrix = methcont, fdr = 0.05, coef = "CGvsEA")
str(myAnnotation_group)

# look for DMRs within 500bp
DMRs_group <- dmrcate(myAnnotation_group, lambda = 500)

(results.ranges <- extractRanges(DMRs_group))

results.stats <- data.frame(myAnnotation_group@ranges@elementMetadata@listData)
results.cpgs <- data.frame(myAnnotation_group@ranges@ranges@start)

(singleCpgs <- results.cpgs[results.stats$is.sig == T, ])

dmrOx <- pos[pos >= results.ranges@ranges@start & pos <= results.ranges@ranges@start + results.ranges@ranges@width - 1]

# add chromosome to cpg position so set as name
dmrOx <- paste("chr3", dmrOx, sep = ".")

# get all cpgs within the ranges
m <- data.frame(t(methylation))

dmr <- m[names(m) %in% dmrOx]

# add variables group and total ctq score
dmr <- cbind(df_outcomes$group, df_outcomes$ctq, dmr)
test <- data.frame(dmr[, 2:15])

# get a correlation matrix of all cpgs within the ranges and ctq sum score
psych::corr.test(test)

# correlation table output

#corDMR <- corr.test(test)
#apaTables::apa.cor.table(corDMR, filename = "outputs/correlations_DMR_with_ctq.doc")

dmr$meanDMR <- rowMeans(dmr[3:15], na.rm = T)
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

myAnnotation_group_age <- sequencing.annotate(bseq_obj, methdesign = methdesign_age, all.cov = T, contrasts = F, fdr = 0.05, coef = "groupEA")
str(myAnnotation_group_age)
DMRs_group_age <- dmrcate(myAnnotation_group_age, lambda = 500)

(results.ranges_age <- extractRanges(DMRs_group_age))

# dmrcate suggests 3 individual cpgs that are significant when controlling for age (instead of 4), but the DMR still contains 13 cpgs 

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

dmr_age$meanDMR <- rowMeans(dmr_age[3:15], na.rm = T)
names(dmr_age)[1:2] <- c("group", "ctq")

with(dmr_age, t.test(meanDMR ~ group, var.equal = T))


## Sensitivity analysis - repeat DMR analysis using different base pair cut-offs defining a region

# look for DMRs within 500bp
DMRs_group250 <- dmrcate(myAnnotation_group, lambda = 250)
(results.ranges250 <- extractRanges(DMRs_group250)) # only 4 cpgs within 8809179-8809231   

DMRs_group750 <- dmrcate(myAnnotation_group, lambda = 750)
(results.ranges750 <- extractRanges(DMRs_group750)) # only 8 cpgs within 8809018-8809113


##############################

## Same analyses for gene expression instead of group (here statistically as predictor, although conceptually gene expression is the outcome)

# single CpGs

# subset of mothers with gene expression data
df_expr <- df_outcomes[complete.cases(df_outcomes$genExpr), ]
df_CpG_m_expr <- df_CpG_m[rownames(df_CpG_m) %in% rownames(df_expr), ]

methylation <- t(df_CpG_m_expr)

design_expr <- model.matrix(~1 + genExpr, data = df_expr)

# fit the linear model for methylation values
fit.expr <- lmFit(methylation, design_expr)
(fit.expr2 <- eBayes(fit.expr))

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit.expr2)) # no individually significant SNPs (low power with N = 69)

DMPs_geneExpr <- topTable(fit.expr2,  num = Inf)
head(DMPs_geneExpr, n = 20)


# DMRs

# get chromosomal position by merging with glossar
m <- data.frame(methylation)
m$CpG_Nr <- row.names(m)
m$CpG_Nr <- gsub('m_', '', m$CpG_Nr)

methylation <- dplyr::left_join(m, glossar)
row.names(methylation) <- methylation$Chromosomal_Location
position <- row.names(methylation)

methylation <- as.matrix(methylation[, 1:69])

position <- stringr::str_split_fixed(position, ":", 2)

chr <- position[,1]
pos <- as.numeric(position[,2])
rownames(methylation) <- paste(chr, pos, sep = ":")

design_excl <- model.matrix(~1 + genExpr, data = df_expr)

methdesign <- edgeR::modelMatrixMeth(design_excl)
coverage2 <- coverage[, c(1:nrow(df_expr))]

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation, Cov = coverage2, pos = pos, chr = chr, sampleNames = rownames(df_expr))

# test cpgs
myAnnotation_geneExpr <- sequencing.annotate(bseq_obj, methdesign = methdesign, all.cov = T, contrasts = F, fdr = 0.05, coef = "genExpr")
str(myAnnotation_geneExpr)

# look for DMRs within 500bp
DMRs_geneExpr <- dmrcate(myAnnotation_geneExpr, lambda = 500) # no individual CpG that are significant, so DMRcate does not work

# only significant CpGs at FDR 0.30 (NOT RELIABLE AT ALL, RESULTS ARE NOT MEANINGFUL!!)

DMRs_geneExprFDR30 <- dmrcate(changeFDR(myAnnotation_geneExpr,FDR = 0.3), lambda = 500)

(results.ranges_geneExprFDR30 <- extractRanges(DMRs_geneExprFDR30)) # including 40 (!) CpGs, but again - NOT MEANINGFUL 



#############################

## check if results of single Cpgs (DMPs) are similar when using t.tests/linear models instead of limma

# group differences
cpgsLev <- lapply(df_CpG_m[,], function(x) car::leveneTest(x ~ df_outcomes$group)$`Pr(>F)`)
sort(unlist(cpgsLev)) # no equal variances for 7 cpgs

cpgsTtest <- lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group))

t <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$statistic[1]))
df <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$statistic[1]))
SE <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$stderr))
p <- unlist(lapply(df_CpG_m[,], function(x) t.test(x ~ df_outcomes$group)$p.value))

results.ttest <- data.frame(t, df, SE, p)

results.ttest <- results.ttest[order(results.ttest$p),]
results.ttest$FDR <- p.adjust(results.ttest$p, method = "BH", n = dim(results.ttest)[1])
names(results.ttest)[4] <- "Pval"


# linear models predicting gene expression

cpgsLm <- lapply(df_CpG_m_expr[,], function(x) lm(df_expr$genExpr ~ x))
results.lm <- t(data.frame(lapply(df_CpG_m_expr[,], function(x) summary(lm(scale(df_expr$genExpr) ~ scale(x)))$coefficients[2, ])))
results.lm <- data.frame(results.lm)

results.lm <- results.lm[order(results.lm$Pr...t..),]
results.lm$FDR <- p.adjust(results.lm$Pr...t.., method = "BH", n = dim(results.lm)[1])

names(results.lm)[c(2,4)] <- c("SE", "Pval")

##############################


# write outputs in excel sheets
require(openxlsx)
list_results <- list("DMPs_Group" = DMPs, "DMPs_Group_ageControlled" = DMPs.age,
                     "DMPs_geneExpression" = DMPs_geneExpr[, 1:6],
                     "DMR_group" = results.ranges,
                     "DMR_group250bp" = results.ranges250,
                     "DMR_group750bp" = results.ranges750,
                     "DMR_group_ageControlled" = results.ranges_age,
                     "T-tests_Group" = results.ttest,
                     "LinearModels_GeneExpr" = results.lm)
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

png(filename = "Figures/OXTR_singleCpgs_Group_and_genExpr.png", width = 2000, height = 1500, type = "cairo", res = 300)
ggplot() + 
  geom_point(data = oxtrPlot[-dmr_pos, ], aes(oxtrPlot$pos[-dmr_pos], oxtrPlot$log10_P_Group[-dmr_pos], colour = factor(oxtrPlot$segment2, levels = unique(oxtrPlot$segment2))[-dmr_pos]), shape = 2) +
  geom_point(data = oxtrPlot[dmr_pos, ], aes(oxtrPlot$pos[dmr_pos], oxtrPlot$log10_P_Group[dmr_pos], colour = factor(oxtrPlot$segment2, levels = unique(oxtrPlot$segment2))[dmr_pos]), shape = 15) +
  geom_point(data = oxtrPlot, aes(pos, log10_P_GeneExpr, colour = factor(segment2, levels = unique(segment2))), shape = 1) +
  geom_hline(yintercept = (-log10(0.001)), linetype = 1) + # fdr 5%
  #geom_hline(yintercept = (-log10(0.007)), linetype = 2) + # fdr 10%
  geom_hline(yintercept = (-log10(0.05)), linetype = 2) + # nominal p < 0.05
  ylim(0, 4) + labs(x = "Position on chromosome 3", y = "-log10 P-Value") + 
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(colour = '') +
  scale_x_reverse()
dev.off()
