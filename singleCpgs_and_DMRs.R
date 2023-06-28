

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
colnames(coverage) <- glossar$Chromosomal_Location

coverage <- data.frame(coverage[rep(seq_len(nrow(coverage)), nrow(df_CpG_m)), ])

colnames(coverage) <- str_replace(colnames(coverage), "chr3.", "chr3:")

names(df_outcomes) <- c("group", "ctq", "genExpr", "age_mother")

# make factor out of CTQ variable (0 = control group, 1 = early adversity? oder andersherum?)
df_outcomes$group <- factor(df_outcomes$group, levels = 0:1, labels = c("CG", "EA"))



###################################

## ANALYSES - Identify differentially-methylated positions (DMPs) between groups using limma (other methods can be used)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("limma")) BiocManager::install("limma")

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

#options(connectionObserver = NULL)

if (!requireNamespace("DMRcate")) BiocManager::install("DMRcate")
if (!requireNamespace("bsseq")) BiocManager::install("bsseq")

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
row.names(methylation) <- paste(chr, pos, sep = ":")

design_excl <- model.matrix(~0 + group, data = df_outcomes)

methdesign <- edgeR::modelMatrixMeth(design_excl)
methcont <- makeContrasts(CGvsEA = groupCG - groupEA, levels = methdesign)

coverage <- t(coverage[c(1:nrow(df_outcomes)), ])
colnames(methylation) <- row.names(df_outcomes)
colnames(coverage) <- row.names(df_outcomes)


coverage <- coverage[c(str_order(row.names(coverage), decreasing = T)), ]

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation, Cov = coverage, pos = pos, chr = chr, sampleNames = row.names(df_outcomes))

# test cpgs
myAnnotation_group <- sequencing.annotate(bseq_obj, methdesign = methdesign, all.cov = T, contrasts = T, cont.matrix = methcont, fdr = 0.05, coef = "CGvsEA")
str(myAnnotation_group)

# look for DMRs within 500bp
DMRs_group <- dmrcate(myAnnotation_group, lambda = 500)

# results table
tools::R_user_dir("ExperimentHub", which="cache")
moveFiles<-function(package){
  olddir <- path.expand(rappdirs::user_cache_dir(appname=package))
  newdir <- tools::R_user_dir(package, which="cache")
  dir.create(path=newdir, recursive=TRUE)
  files <- list.files(olddir, full.names =TRUE)
  moveres <- vapply(files,
                    FUN=function(fl){
                      filename = basename(fl)
                      newname = file.path(newdir, filename)
                      file.rename(fl, newname)
                    },
                    FUN.VALUE = logical(1))
  if(all(moveres)) unlink(olddir, recursive=TRUE)
}

package="ExperimentHub"
moveFiles(package)


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
colnames(methylation2) <- row.names(df_expr)
colnames(coverage2) <- row.names(df_expr)

# create bseq object needed for sequencing.annotate() function
bseq_obj <- BSseq(M = methylation2, Cov = coverage2, pos = pos, chr = chr, sampleNames = rownames(df_expr))

# test cpgs
myAnnotation_geneExpr <- sequencing.annotate(bseq_obj, methdesign = methdesign, all.cov = T, contrasts = F, fdr = 0.05, coef = "genExpr")
str(myAnnotation_geneExpr)

# look for DMRs within 500bp
# DMRs_geneExpr <- dmrcate(myAnnotation_geneExpr, lambda = 500) # no individual CpG that are significant, so DMRcate does not work

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
colnames(methylation) <- row.names(df_outcomes)

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

results.ttest <- results.ttest[order(results.ttest$p),]
results.ttest$FDR <- p.adjust(results.ttest$p, method = "BH", n = dim(results.ttest)[1])

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

names <- rownames(results.lm.ctq)[1:10]

par(mfrow = c(5, 4))
for (i in 1:5) {
 plot(lm(scale(df_CpG_m[names[i]]) ~ scale(df_outcomes$ctq)))
}

par(mfrow = c(5, 4))
for (i in 6:10) {
  plot(lm(scale(df_CpG_m[names[i]]) ~ scale(df_outcomes$ctq)))
}


# linear models predicting gene expression

results.lm <- t(data.frame(lapply(df_CpG_m_expr[,], function(x) summary(lm(scale(df_expr$genExpr) ~ scale(x)))$coefficients[2, ])))
results.lm <- data.frame(results.lm)

results.lm <- results.lm[order(results.lm$Pr...t..),]
results.lm$FDR <- p.adjust(results.lm$Pr...t.., method = "BH", n = dim(results.lm)[1])

names(results.lm)[c(2,4)] <- c("SE", "Pval")


results.lm$CI_lower <- with(results.lm, Estimate - 1.96 * SE)
results.lm$CI_upper <- with(results.lm, Estimate + 1.96 * SE)


results.lm <- results.lm[, c(1, 6:7, 2:5)]


names <- rownames(results.lm)[1:10]

par(mfrow = c(5, 4))
for (i in 1:5) {
  plot(lm(scale(df_outcomes$genExpr) ~ scale(df_CpG_m[names[i]])))
}

par(mfrow = c(5, 4))
for (i in 6:10) {
  plot(lm(scale(df_outcomes$genExpr) ~ scale(df_CpG_m[names[i]])))
}

##############################

results.lm.ctq$CpG_Nr <- row.names(results.lm.ctq)
results.ttest$CpG_Nr <- row.names(results.ttest)
results.lm$CpG_Nr <- row.names(results.lm)

results.lm.ctq$CpG_Nr <- stringr::str_replace(results.lm.ctq$CpG_Nr, "_m_", "_")
results.ttest$CpG_Nr <- stringr::str_replace(results.ttest$CpG_Nr, "_m_", "_")
results.lm$CpG_Nr <- stringr::str_replace(results.lm$CpG_Nr, "_m_", "_")


results.lm.ctq2 <- dplyr::left_join(results.lm.ctq, glossar)
results.lm.ctq2$position <- stringr::str_split_fixed(results.lm.ctq2$Chromosomal_Location, ":", 2)[, 2]
results.lm.ctq2[, 1:5] <- round(results.lm.ctq2[, 1:5], 3)
results.lm.ctq2[, 6] <- round(results.lm.ctq2[, 6], 7)
results.lm.ctq2[, 7] <- round(results.lm.ctq2[, 7], 6)

results.ttest2 <- dplyr::left_join(results.ttest, glossar)
results.ttest2$position <- stringr::str_split_fixed(results.ttest2$Chromosomal_Location, ":", 2)[, 2]
results.ttest2[, c(1:6, 8:9)] <- round(results.ttest2[, c(1:6, 8:9)], 3)
results.ttest2[, 7] <- round(results.ttest2[, 7], 5)

results.lm2 <- dplyr::left_join(results.lm, glossar)
results.lm2$position <- stringr::str_split_fixed(results.lm2$Chromosomal_Location, ":", 2)[, 2]
results.lm2[, 1:7] <- round(results.lm2[, 1:7], 3)

require(openxlsx)
write.xlsx(results.lm.ctq2, file = "Results/CTQ.xlsx", colNames = T, rowNames = T)
write.xlsx(results.ttest2, file = "Results/Group.xlsx", colNames = T, rowNames = T)
write.xlsx(results.lm2, file = "Results/mRNAexpr.xlsx", colNames = T, rowNames = T)




#############################


# write outputs in excel sheets

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
dmr_pos_name <- oxtrPlot$pos[dmr_pos]

oxtrPlot$pos <- factor(oxtrPlot$pos, levels = unique(oxtrPlot$pos))

oxtrPlot$CpG <- factor(oxtrPlot$CpG, labels = 1:length(oxtrPlot$CpG))

library(tidyverse)
oxtrPlot2 <- oxtrPlot %>% gather("Association_with", "Pvalue", c(log10_P_Group, log10_P_GeneExpr))

####

oxtrPlot2$Association_with[pos %in% dmr_pos_name] <- "DMR"

oxtrPlot2$Association_with <- dplyr::recode(oxtrPlot2$Association_with,
                                            "log10_P_Group" = "Group (single CpGs)",
                                            "log10_P_GeneExpr" = "mRNA_expr (single CpGs)",
                                            "DMR" = "Group (DMR)")
oxtrPlot2$Association_with <- factor(oxtrPlot2$Association_with ,
                                     levels = c("Group (single CpGs)", "mRNA_expr (single CpGs)", "Group (DMR)"))


png(filename = "Figures/FigureS12.png", width = 6000, height = 4500, type = "cairo", res = 900)
ggplot() + 
  geom_point(data = oxtrPlot2, 
             aes(pos, 
                 Pvalue, 
                 colour = factor(segment2, levels = unique(segment2)),
                 shape = Association_with,
                 group = interaction(segment2,
                                     Association_with)
             )
  ) +
  geom_hline(yintercept = (-log10(0.05)), linetype = 2) + # nominal p < 0.05
  geom_hline(yintercept = (-log10(0.001)), linetype = 1) +
  ylim(0, 5) + labs(x = "OXTR CpG Site", y = "-log10 P-Value") + 
  theme_classic() +
  theme(legend.position = "right") +
  labs(colour = '', shape = '') +
  scale_x_discrete(labels = NULL) +
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
dmr_pos_name <- oxtrPlot$pos[dmr_pos]

oxtrPlot$pos <- factor(oxtrPlot$pos, levels = unique(oxtrPlot$pos))

oxtrPlot$CpG <- factor(oxtrPlot$CpG, labels = 1:length(oxtrPlot$CpG))

library(tidyverse)
oxtrPlot2 <- oxtrPlot %>% gather("Association_with", "Pvalue", c(log10_P_CTQ, log10_P_GeneExpr))

oxtrPlot2$Association_with[pos %in% dmr_pos_name] <- "DMR"

oxtrPlot2$Association_with <- dplyr::recode(oxtrPlot2$Association_with,
                                            "log10_P_CTQ" = "CTQ (single CpGs)",
                                            "log10_P_GeneExpr" = "mRNA_expr (single CpGs)",
                                            "DMR" = "CTQ (DMR)")
oxtrPlot2$Association_with <- factor(oxtrPlot2$Association_with ,
                                     levels = c("CTQ (single CpGs)", "mRNA_expr (single CpGs)", "CTQ (DMR)"))


png(filename = "Figures/Figure5.png", width = 6000, height = 4500, type = "cairo", res = 900)
ggplot() + 
  geom_point(data = oxtrPlot2, 
             aes(pos, 
                 Pvalue, 
                 colour = factor(segment2, levels = unique(segment2)),
                 shape = Association_with,
                 group = interaction(segment2,
                                     Association_with)
             )
  ) +
  geom_hline(yintercept = (-log10(0.05)), linetype = 2) + # nominal p < 0.05
  geom_hline(yintercept = (-log10(0.001)), linetype = 1) +
  ylim(0, 7) + labs(x = "OXTR CpG Site", y = "-log10 P-Value") + 
  theme_classic() +
  theme(legend.position = "right") +
  labs(colour = '', shape = '') +
  scale_x_discrete(labels = NULL) +
  scale_color_brewer(palette = "Spectral")
dev.off()




### power curves 

library(pwr)
library(tidyverse)
library(broom)


effect_sizes <- seq(0.05, 0.8, 0.05) 
sample_sizes <- seq(10, 500, 10) # going from 10 to 500 in increments of 10
sig_level <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

input_df_n <- crossing(effect_sizes, sample_sizes)
input_df <- crossing(effect_sizes, sig_level)


get_power <- function(df, n){
  power_result <- pwr.t.test(sig.level = df %>% pull(sig_level), 
                             d = df %>% pull(effect_sizes), n = n) %>% tidy()
  df <- df %>% mutate(power = power_result %>% pull(power))
  return(df)
}


power_curves <- get_power(input_df, n = 55)


power_curves <- power_curves %>% mutate(effect_sizes = as.factor(effect_sizes)) 


cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
          "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
          "#1B9E77", "#D95F02", "#7570B3", "#E7298A") 

png(filename = "Figures/power_ttest_fixed_n.png", width = 2000, height = 1500, type = "cairo", res = 250)
pwr_t <- ggplot(power_curves, 
       aes(x = sig_level,
           y = power, 
           colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0.001, 0.5)) + #this changes the limits of x
  scale_color_manual(values =  cols) +
  labs(col = "Cohen's d", x = "Significance level", y = "Power")
pwr_t
dev.off()

get_power_n <- function(df){
  power_result <- pwr.t.test(n = df %>% pull(sample_sizes), 
                             d = df %>% pull(effect_sizes)) %>% tidy()
  df <- df %>% mutate(power = power_result %>% pull(power))
  return(df)
}


power_curves2 <- get_power_n(input_df_n)


power_curves2 <- power_curves2 %>% mutate(effect_sizes = as.factor(effect_sizes)) 

png(filename = "Figures/power_ttest_fixed_alpha.png", width = 7500, height = 4500, type = "cairo", res = 750)
pwr_t2 <- ggplot(power_curves2, 
       aes(x = sample_sizes,
           y = power, 
           colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 500, 25), limits = c(10, 500)) + #this changes the limits of x
  scale_color_manual(values =  cols) +
  labs(col = "Cohen's d", x = "Sample size", y = "Power")
pwr_t2
dev.off()


get_power <- function(df, n){
  power_result <- pwr.r.test(sig.level = df %>% pull(sig_level), 
                             r = df %>% pull(effect_sizes), n = n) %>% tidy()
  df <- df %>% mutate(power = power_result %>% pull(power))
  return(df)
}


power_curves <- get_power(input_df, n = 110)


power_curves <- power_curves %>% mutate(effect_sizes = as.factor(effect_sizes)) 


png(filename = "Figures/power_cor_fixed_n.png", width = 2000, height = 1500, type = "cairo", res = 250)
pwr_r <- ggplot(power_curves, 
       aes(x = sig_level,
           y = power, 
           colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0.001, 0.5)) + #this changes the limits of x
  scale_color_manual(values =  cols) +
  labs(col = "Correlation", x = "Significance level", y = "Power")
pwr_r 
dev.off()


get_power_n <- function(df){
  power_result <- pwr.r.test(n = df %>% pull(sample_sizes), 
                             r = df %>% pull(effect_sizes)) %>% tidy()
  df <- df %>% mutate(power = power_result %>% pull(power))
  return(df)
}


power_curves2 <- get_power_n(input_df_n)


power_curves2 <- power_curves2 %>% mutate(effect_sizes = as.factor(effect_sizes)) 

png(filename = "Figures/power_cor_fixed_alpha.png", width = 2500, height = 1500, type = "cairo", res = 250)
pwr_r2 <- ggplot(power_curves2, 
       aes(x = sample_sizes,
           y = power, 
           colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 500, 25), limits = c(10, 500)) + #this changes the limits of x
  scale_color_manual(values =  cols) +
  labs(col = "Correlation", x = "Sample size", y = "Power")
pwr_r2 + theme_classic()
dev.off()


png(filename = "Figures/FigureS13.png", width = 7000, height = 4000, type = "cairo", res = 600)
cowplot::plot_grid(pwr_t, pwr_r, labels = c("a", "b"), ncol = 2, nrow = 1)
dev.off()

png(filename = "Figures/FigureS14.png", width = 9000, height = 4000, type = "cairo", res = 700)
cowplot::plot_grid(pwr_t2, pwr_r2, labels = c("a", "b"), ncol = 2, nrow = 1)
dev.off()



#### less levels of effect sizes

effect_sizes <- seq(0.1, 0.8, 0.1) 
sample_sizes <- seq(10, 500, 10) # going from 10 to 500 in increments of 10
sig_level <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)

input_df_n <- crossing(effect_sizes, sample_sizes)
input_df <- crossing(effect_sizes, sig_level)

power_curves <- get_power(input_df, n = 55)

power_curves <- power_curves %>% mutate(effect_sizes = as.factor(effect_sizes)) 


pwr_t <- ggplot(power_curves, 
                aes(x = sig_level,
                    y = power, 
                    colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0.001, 0.5)) + #this changes the limits of x
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Cohen's d", x = "Significance level", y = "Power")


power_curves2 <- get_power_n(input_df_n)

power_curves2 <- power_curves2 %>% mutate(effect_sizes = as.factor(effect_sizes)) 

pwr_t2 <- ggplot(power_curves2, 
                 aes(x = sample_sizes,
                     y = power, 
                     colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 500, 25), limits = c(10, 500)) + #this changes the limits of x
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Cohen's d", x = "Sample size", y = "Power")

power_curves <- get_power(input_df, n = 110)

power_curves <- power_curves %>% mutate(effect_sizes = as.factor(effect_sizes)) 


pwr_r <- ggplot(power_curves, 
                aes(x = sig_level,
                    y = power, 
                    colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0.001, 0.5)) + #this changes the limits of x
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Correlation", x = "Significance level", y = "Power")


power_curves2 <- get_power_n(input_df_n)

power_curves2 <- power_curves2 %>% mutate(effect_sizes = as.factor(effect_sizes)) 

pwr_r2 <- ggplot(power_curves2, 
                 aes(x = sample_sizes,
                     y = power, 
                     colour = effect_sizes)) + 
  geom_line()  + theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + #this line change the frequency of tick marks
  scale_x_continuous(breaks = seq(0, 500, 25), limits = c(10, 500)) + #this changes the limits of x
  scale_color_brewer(palette = "Spectral") +
  labs(col = "Correlation", x = "Sample size", y = "Power")

png(filename = "Figures/power_fixed_n_spectral.png", width = 3500, height = 2000, type = "cairo", res=300)
cowplot::plot_grid(pwr_t, pwr_r, labels = c("a", "b"), ncol = 2, nrow = 1)
dev.off()

png(filename = "Figures/power_fixed_alpha_spectral.png", width = 4500, height = 2000, type = "cairo", res=350)
cowplot::plot_grid(pwr_t2, pwr_r2, labels = c("a", "b"), ncol = 2, nrow = 1)
dev.off()


sessionInfo()

# R version 4.1.2 (2021-11-01)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 11.5.2
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
# 
# attached base packages:
#   [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] broom_0.7.12                pwr_1.3-0                   forcats_0.5.1               dplyr_1.0.8                
# [5] purrr_0.3.4                 readr_2.1.2                 tidyr_1.2.0                 tibble_3.1.6               
# [9] tidyverse_1.3.1             DMRcatedata_2.12.0          ExperimentHub_2.2.1         AnnotationHub_3.2.2        
# [13] BiocFileCache_2.2.1         dbplyr_2.1.1                bsseq_1.30.0                SummarizedExperiment_1.24.0
# [17] Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.61.0          GenomicRanges_1.46.1       
# [21] GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.3            BiocGenerics_0.40.0        
# [25] DMRcate_2.8.5               limma_3.50.1                here_1.0.1                  plyr_1.8.6                 
# [29] pls_2.8-0                   caret_6.0-90                lattice_0.20-45             reshape_0.8.8              
# [33] factoextra_1.0.7            EnvStats_2.6.0              ggplot2_3.3.5               stringr_1.4.0              
# [37] dbscan_1.1-10               fpc_2.2-9                   psych_2.1.9                 party_1.3-9                
# [41] strucchange_1.5-2           sandwich_3.0-1              zoo_1.8-9                   modeltools_0.2-23          
# [45] mvtnorm_1.1-3               mclust_5.4.9                corrplot_0.92               pacman_0.5.1               
# 
# loaded via a namespace (and not attached):
#   [1] Hmisc_4.6-0                                         class_7.3-20                                       
# [3] Rsamtools_2.10.0                                    foreach_1.5.2                                      
# [5] rprojroot_2.0.2                                     crayon_1.5.0                                       
# [7] MASS_7.3-55                                         rhdf5filters_1.6.0                                 
# [9] nlme_3.1-155                                        backports_1.4.1                                    
# [11] reprex_2.0.1                                        rlang_1.0.1                                        
# [13] XVector_0.34.0                                      readxl_1.3.1                                       
# [15] minfi_1.40.0                                        DSS_2.42.0                                         
# [17] filelock_1.0.2                                      BiocParallel_1.28.3                                
# [19] rjson_0.2.21                                        bit64_4.0.5                                        
# [21] glue_1.6.2                                          rngtools_1.5.2                                     
# [23] parallel_4.1.2                                      AnnotationDbi_1.56.2                               
# [25] haven_2.4.3                                         tidyselect_1.1.2                                   
# [27] XML_3.99-0.9                                        GenomicAlignments_1.30.0                           
# [29] xtable_1.8-4                                        magrittr_2.0.2                                     
# [31] cli_3.2.0                                           zlibbioc_1.40.0                                    
# [33] rstudioapi_0.13                                     doRNG_1.8.2                                        
# [35] rpart_4.1.16                                        effsize_0.8.1                                      
# [37] ensembldb_2.18.3                                    IlluminaHumanMethylationEPICanno.ilm10b4.hg19_0.6.0
# [39] shiny_1.7.1                                         xfun_0.30                                          
# [41] askpass_1.1                                         multtest_2.50.0                                    
# [43] cluster_2.1.2                                       KEGGREST_1.34.0                                    
# [45] interactiveDisplayBase_1.32.0                       ggrepel_0.9.1                                      
# [47] base64_2.0                                          biovizBase_1.42.0                                  
# [49] scrime_1.3.5                                        listenv_0.8.0                                      
# [51] Biostrings_2.62.0                                   png_0.1-7                                          
# [53] permute_0.9-7                                       future_1.24.0                                      
# [55] ipred_0.9-12                                        withr_2.4.3                                        
# [57] bitops_1.0-7                                        cellranger_1.1.0                                   
# [59] AnnotationFilter_1.18.0                             hardhat_0.2.0                                      
# [61] pROC_1.18.0                                         pillar_1.7.0                                       
# [63] bumphunter_1.36.0                                   cachem_1.0.6                                       
# [65] GenomicFeatures_1.46.5                              multcomp_1.4-18                                    
# [67] fs_1.5.2                                            flexmix_2.3-17                                     
# [69] kernlab_0.9-29                                      DelayedMatrixStats_1.16.0                          
# [71] vctrs_0.3.8                                         ellipsis_0.3.2                                     
# [73] generics_0.1.2                                      nortest_1.0-4                                      
# [75] lava_1.6.10                                         tools_4.1.2                                        
# [77] foreign_0.8-82                                      munsell_0.5.0                                      
# [79] DelayedArray_0.20.0                                 fastmap_1.1.0                                      
# [81] compiler_4.1.2                                      abind_1.4-5                                        
# [83] httpuv_1.6.5                                        rtracklayer_1.54.0                                 
# [85] beanplot_1.2                                        Gviz_1.38.3                                        
# [87] GenomeInfoDbData_1.2.7                              prodlim_2019.11.13                                 
# [89] gridExtra_2.3                                       edgeR_3.36.0                                       
# [91] utf8_1.2.2                                          later_1.3.0                                        
# [93] recipes_0.2.0                                       jsonlite_1.8.0                                     
# [95] scales_1.1.1                                        carData_3.0-5                                      
# [97] sparseMatrixStats_1.6.0                             genefilter_1.76.0                                  
# [99] lazyeval_0.2.2                                      promises_1.2.0.1                                   
# [101] car_3.0-12                                          latticeExtra_0.6-29                                
# [103] R.utils_2.11.0                                      checkmate_2.0.0                                    
# [105] nor1mix_1.3-0                                       cowplot_1.1.1                                      
# [107] statmod_1.4.36                                      siggenes_1.68.0                                    
# [109] dichromat_2.0-0                                     BSgenome_1.62.0                                    
# [111] HDF5Array_1.22.1                                    survival_3.3-0                                     
# [113] yaml_2.3.5                                          prabclus_2.3-2                                     
# [115] htmltools_0.5.2                                     memoise_2.0.1                                      
# [117] VariantAnnotation_1.40.0                            BiocIO_1.4.0                                       
# [119] locfit_1.5-9.4                                      quadprog_1.5-8                                     
# [121] digest_0.6.29                                       assertthat_0.2.1                                   
# [123] mime_0.12                                           rappdirs_0.3.3                                     
# [125] RSQLite_2.2.10                                      future.apply_1.8.1                                 
# [127] data.table_1.14.2                                   blob_1.2.2                                         
# [129] R.oo_1.24.0                                         preprocessCore_1.56.0                              
# [131] splines_4.1.2                                       Formula_1.2-4                                      
# [133] labeling_0.4.2                                      Rhdf5lib_1.16.0                                    
# [135] illuminaio_0.36.0                                   ProtGenerics_1.26.0                                
# [137] RCurl_1.98-1.6                                      hms_1.1.1                                          
# [139] modelr_0.1.8                                        rhdf5_2.38.0                                       
# [141] colorspace_2.0-3                                    base64enc_0.1-3                                    
# [143] BiocManager_1.30.16                                 mnormt_2.0.2                                       
# [145] tmvnsim_1.0-2                                       libcoin_1.0-9                                      
# [147] nnet_7.3-17                                         GEOquery_2.62.2                                    
# [149] Rcpp_1.0.8                                          coin_1.4-2                                         
# [151] fansi_1.0.2                                         tzdb_0.2.0                                         
# [153] parallelly_1.30.0                                   ModelMetrics_1.2.2.2                               
# [155] R6_2.5.1                                            lifecycle_1.0.1                                    
# [157] curl_4.3.2                                          robustbase_0.93-9                                  
# [159] Matrix_1.4-0                                        TH.data_1.1-0                                      
# [161] org.Hs.eg.db_3.14.0                                 RColorBrewer_1.1-2                                 
# [163] iterators_1.0.14                                    gower_1.0.0                                        
# [165] htmlwidgets_1.5.4                                   biomaRt_2.50.3                                     
# [167] missMethyl_1.28.0                                   rvest_1.0.2                                        
# [169] globals_0.14.0                                      openssl_2.0.0                                      
# [171] htmlTable_2.4.0                                     codetools_0.2-18                                   
# [173] IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0  lubridate_1.8.0                                    
# [175] gtools_3.9.2                                        prettyunits_1.1.1                                  
# [177] R.methodsS3_1.8.1                                   gtable_0.3.0                                       
# [179] DBI_1.1.2                                           httr_1.4.2                                         
# [181] stringi_1.7.6                                       progress_1.2.2                                     
# [183] reshape2_1.4.4                                      farver_2.1.0                                       
# [185] diptest_0.76-0                                      annotate_1.72.0                                    
# [187] timeDate_3043.102                                   xml2_1.3.3                                         
# [189] restfulr_0.0.13                                     BiocVersion_3.14.0                                 
# [191] DEoptimR_1.0-10                                     bit_4.0.4                                          
# [193] jpeg_0.1-9                                          pkgconfig_2.0.3                                    
# [195] knitr_1.37     
