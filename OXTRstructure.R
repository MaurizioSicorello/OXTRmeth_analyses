

# corrplots mit allen CpGs

# meeting:
# talk about shinyapp
# plan paper

# ideas for shinyapp:
# variance compared to sensitivity, percent outliers, skew, p-value DMR, bivariate correlations, Bayes faktoren, cluster assignment, 
# PLS loadings (?), both for kids and mothers





#################################
# Load and prepare packages/data

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(corrplot, mclust, party, psych, fpc, dbscan, stringr, ggplot2, EnvStats, stats, factoextra, reshape, caret, pls, plyr, here)

# load custom functions
source(here("Functions", "OXTRmeth_helperFuns.R"))

# load data and convert % to decimal
df = read.csv(here("Data", "RUB_OXTR_Daten_26.4.csv"), header = T, sep = ";")

# subset mother CpGs and convert to decimal
df_CpG_m = df[,grepl("CpG_m", names(df))]
df_CpG_m = as.data.frame(apply(df_CpG_m,2, function(x){as.numeric(sub("%", "", x, fixed=TRUE))/100}))
df_outcomes = df[complete.cases(df_CpG_m), c("ctq_sum", "Genexp_OXTR_mother", "ctq_dich01_multi01")]
df_CpG_m = df_CpG_m[complete.cases(df_CpG_m),]

# load identifier of gene sections
dfcodes = read.csv(here("Data", "OXTR_segmentCodes.csv"), header = T, sep = ";")
dfcodes = dfcodes[dfcodes$CpG %in% str_replace(names(df_CpG_m), "_m_", "_"),]


#################################
# qualitative description of OXTR

# plot methylation across CpGs to check for sufficient variance
boxplot(df_CpG_m, range = 1.5)


# methylation on Exon 1
median(as.matrix(df_CpG_m[,which(dfcodes$segment2 == "exon1")]))
iqr(as.matrix(df_CpG_m[,which(dfcodes$segment2 == "exon1")]))
describe(df_CpG_m)
mean(describe(df_CpG_m)[which(dfcodes$segment2 == "exon1"), "skew"])


# methylation on MT2
median(as.matrix(df_CpG_m[,which(dfcodes$segment2 == "MT2")]))
iqr(as.matrix(df_CpG_m[,which(dfcodes$segment2 == "MT2")]))
describe(df_CpG_m)
mean(describe(df_CpG_m)[which(dfcodes$segment2 == "MT2"), "skew"])


# methylation on enhancer
median(as.matrix(df_CpG_m[,which(dfcodes$segment2 == "enhancer")]))
iqr(as.matrix(df_CpG_m[,which(dfcodes$segment2 == "enhancer")]))
describe(df_CpG_m)
mean(describe(df_CpG_m)[which(dfcodes$segment2 == "enhancer"), "skew"])


# count/plot number of outliers per CpG (3*IQR)
out_count = as.data.frame(matrix(nrow = ncol(df_CpG_m), ncol = 2, dimnames = list(NULL, c("CpG", "nOut"))))
out_count$CpG = str_replace(names(df_CpG_m), "_m_", "_")
for(i in 1:ncol(df_CpG_m)){out_count[i,2] = length(boxplot.stats(df_CpG_m[,i], coef = 3)$out)}
df_out = merge(out_count, dfcodes, by = "CpG", sort = F)
ggplot(data = df_out, aes(x = reorder(CpG, 1:nrow(df_out)), y = nOut, colour = segment2)) + 
  geom_point()


# correlation plot 
CpG_corr = cor(df_CpG_m, use = "pairwise.complete.obs")
colnames(CpG_corr) <- NULL
rownames(CpG_corr) <- substring(rownames(CpG_corr), 7)

pdf(here("Figures", "corrPlot_mother.pdf"))
corrplot::corrplot(CpG_corr, tl.pos = "d", tl.cex = 0.25, method = "color", type = "upper", tl.col = "black")
dev.off()

mean(CpG_corr[132:182, 132:182])
mean(CpG_corr[53:64, 53:64])






#################################
# unsupervised datra-driven description of OXTR

# z score variables
df_CpG_m_z = scale(df_CpG_m)


#################
# cluster analysis (dbscan)

# transpose data
df_CpG_m_z_trans = t(df_CpG_m_z)

# determine eps. try k = 3, k = 5 (default), k = 8, k = 12
dbscan::kNNdistplot(df_CpG_m_z_trans, k =  3)
abline(h = 5, lty = 2) # for k = 3
dbscan::kNNdistplot(df_CpG_m_z_trans, k =  5)
abline(h = 5.8, lty = 2) # for k = 5
dbscan::kNNdistplot(df_CpG_m_z_trans, k =  8)
abline(h = 7.7, lty = 2) # for k = 8
dbscan::kNNdistplot(df_CpG_m_z_trans, k =  12)
abline(h = 8, lty = 2) # for k = 8

# apply clustering
res.fpc3 <- fpc::dbscan(df_CpG_m_z_trans, eps = 5, MinPts = 3)
print(res.fpc3)
res.fpc5 <- fpc::dbscan(df_CpG_m_z_trans, eps = 5.8, MinPts = 5)
print(res.fpc5)
res.fpc8 <- fpc::dbscan(df_CpG_m_z_trans, eps = 7.7, MinPts = 8)
print(res.fpc8)
res.fpc12 <- fpc::dbscan(df_CpG_m_z_trans, eps = 8, MinPts = 12)
print(res.fpc12)

# write.csv(data.frame(res.fpc5$cluster, names(df_CpG_m)), "clusterAssignment.csv") 

# plot
df_clustMem = data.frame(dfcodes$CpG, dfcodes$segment2, res.fpc3$cluster, res.fpc5$cluster, res.fpc8$cluster, res.fpc12$cluster)
df_clustMem_long = melt(df_clustMem, id.vars = c("dfcodes.CpG", "dfcodes.segment2"))

ggplot(data = df_clustMem_long, aes(x = reorder(dfcodes.CpG, rep((1:nrow(df_clustMem)), length.out = nrow(df_clustMem_long))), y = value, colour = dfcodes.segment2, shape = variable)) + 
  geom_jitter(width = 0, height = 0.3) + scale_y_continuous(breaks=c(0, 1, 2)) + scale_x_discrete(breaks=NULL) +
  xlab("CpG Sites") + ylab("Cluster Assignment") + theme_classic() + theme(legend.title = element_blank())
ggsave(here("Figures", "DBSCAN_clustMembership.pdf"), device = "pdf")
ggsave(here("Figures", "DBSCAN_clustMembership.png"), device = "png")

#################
# hierarchical clustering (hclust)

# dendrogram
dd <- dist(df_CpG_m_z_trans, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

fviz_dend(hc)
fviz_nbclust(df_CpG_m_z_trans, FUN = hcut, method = "wss")
fviz_nbclust(df, FUN = hcut, method = "silhouette")

pdf(here("Figures", "Dendrogram.pdf"))
fviz_dend(
  hc,
  k = 3,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.3
)
dev.off()


png(here("Figures", "Dendrogram_mothers.png"))
fviz_dend(
  hc,
  k = 3,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.3
)
dev.off()



clustAssign_mothers <- cutree(hc, k = 3)

# proportions of genetic sections represented in the three clusters
df_hclustAssign_mothers <- data.frame(dfcodes, clustAssign_mothers) 
round(table(df_hclustAssign_mothers[df_hclustAssign_mothers$clustAssign_mothers == 1, "segment2"])/sum(clustAssign_mothers == 1), 2)
round(table(df_hclustAssign_mothers[df_hclustAssign_mothers$clustAssign_mothers == 2, "segment2"])/sum(clustAssign_mothers == 2), 2)
round(table(df_hclustAssign_mothers[df_hclustAssign_mothers$clustAssign_mothers == 3, "segment2"])/sum(clustAssign_mothers == 3), 2)


df_hclustAssign_mothers$DBSCANclust <- df_clustMem$res.fpc5.cluster
df_hclustAssign_mothers$hclustRec <- ifelse(df_hclustAssign_mothers$clustAssign_mothers == 1, 
                                    0,
                                    ifelse(df_hclustAssign_mothers$clustAssign_mothers == 2, 
                                           1,
                                           2)
                                    )

round(sum(df_hclustAssign_mothers$DBSCANclust == df_hclustAssign_mothers$hclustRec)/nrow(df_hclustAssign_mothers), 2)



#################################
# supervised data-driven description of OXTR



# # mean sensitivity from Moser et al. (2020), Table S7, rightmost column
# meanSens = mean(c(0.67, 1.37, 1.39, 1.51, 1.24, 2.28))/100
# 
# # test which CpGs have variance higher than chance
# varTestFun <- function(x, referenceVar){
#   chi_square = ((length(x)-1)*var(x))/referenceVar
#   p = pchisq(chi_square, df = length(x)-1, lower.tail = F)
#   p
# }
# CpGvariancePs = apply(df_CpG_m, 2, varTestFun, meanSens^2)
# sum(CpGvariancePs > 0.05)
# sum(p.adjust(CpGvariancePs, method = "fdr") > 0.05)
# dfvar_plot = data.frame(dfcodes, CpGvariancePs)
# dfvar_plot$insuffVar = ifelse(CpGvariancePs > 0.05, 1, 0)
# 
# ggplot(data = dfvar_plot, aes(x = reorder(CpG, 1:nrow(df_out)), y = insuffVar, colour = segment2)) + 
#   geom_point()


df_CpGwiseVars <- as.numeric(read.csv2(here("Data", "OXTRStabwSingleCpGs.csv"))[,3])^2
df_CpGwiseVars <- c(df_CpGwiseVars, rep(0.0228^2, times = ncol(df_CpG_m)-length(df_CpGwiseVars)))
dfchi <- nrow(df_CpG_m)-1
CpGvariancePs <- numeric(ncol(df_CpG_m))

# descriptives for sensitivity
median(sqrt(df_CpGwiseVars))
iqr(sqrt(df_CpGwiseVars))
range(sqrt(df_CpGwiseVars))


for(i in 1:ncol(df_CpG_m)){
  
  chi_square <- dfchi*var(df_CpG_m[,i])/df_CpGwiseVars[i]
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



#################
# ML settings
maxcomp = 30
repeats = 2
folds = 5
perms = 1000

plsGrid <- expand.grid(ncomp = seq(1, maxcomp))

#################
# gene expression

# do nested k-fold PLS
PLSnested_Genexpr = PLSnestedCV(df_outcomes$Genexp_OXTR_mother, df_CpG_m, nrepeats = repeats, nfolds = folds, maxComps = maxcomp, setSeed = 1000)
PLSnested_Genexpr[[1]]

# create permutation samples [takes ~2h on my machine]
# plsPerm_all_geneExpr <- permutePLSnestedCV(df_outcomes$Genexp_OXTR_mother, df_CpG_m, nrepeats = repeats, nfolds = folds, nperms = perms)
# write.csv(plsPerm_all_geneExpr, here("Results", "plsPerm_all_geneExpr.csv"))

# p-value
df_geneExpr_perms <- read.csv(here("Results", "plsPerm_all_geneExpr.csv"))
p_geneExpr <- sum(df_geneExpr_perms[,2] >= PLSnested_Genexpr[[1]])/length(df_geneExpr_perms[,2])
p_geneExpr

# plot ncomp solutions for individual folds
df_GeneExpr_optComp <- melt(PLSnested_Genexpr[[3]])
df_GeneExpr_optComp$ncomp <- rep(c(1:maxcomp), times = repeats*folds)
df_GeneExpr_optComp_min <- ddply(df_GeneExpr_optComp, "variable", subset, value == min(value))
ggplot(data = df_GeneExpr_optComp, aes(y = value, x = ncomp, group = variable)) +
  geom_line(colour = "darkgrey") + geom_point(size = 1, colour = "darkgrey") +
  
  geom_point(data = df_GeneExpr_optComp_min, colour = "darkblue") + 
  
  theme_classic() + ylab("RMSE (cross-validated)") + xlab("number of components")
ggsave(here("figures", "numberComponents_geneExpr.png"), device = "png")


# create final model with optimal number of components
geneExpr_finalModel = plsr(DV ~ .,
                            data = PLSnested_Genexpr$dat,
                            ncomp = 1)
summary(geneExpr_finalModel)


# plot loadings
pls_geneExpr_loadings = geneExpr_finalModel$loadings[,1:ncol(geneExpr_finalModel$loadings)]
dfplot_geneExpr_loadings <- data.frame(dfcodes$CpG, pls_geneExpr_loadings, dfcodes$segment2)

ggplot(data = dfplot_geneExpr_loadings, aes(x = reorder(dfcodes.CpG, rep((1:nrow(dfplot_geneExpr_loadings)), length.out = nrow(dfplot_geneExpr_loadings))), y = pls_geneExpr_loadings, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_path(group = 1, size = 0.8) +
  
  theme_classic() + ylim(-max(abs(pls_geneExpr_loadings))*1.5, max(abs(pls_geneExpr_loadings))*1.5) 

# plot bivariate correlations
dfplot_geneExpr_cors <- dfplot_geneExpr_loadings
dfplot_geneExpr_cors$pls_geneExpr_loadings  <- cor(PLSnested_Genexpr$dat)[1, -1]

ggplot(data = dfplot_geneExpr_cors, aes(x = reorder(dfcodes.CpG, rep((1:nrow(dfplot_geneExpr_cors)), length.out = nrow(dfplot_geneExpr_cors))), y = pls_geneExpr_loadings, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_path(group = 1, size = 0.8) +
  
  theme_classic() + ylim(-max(abs(pls_geneExpr_loadings))*1.5, max(abs(pls_geneExpr_loadings))*1.5) 


#################
# childhood trauma [continuous]


# do nested k-fold PLS
PLSnested_CTQ = PLSnestedCV(df_outcomes$ctq_sum, df_CpG_m, nrepeats = repeats, nfolds = folds, maxComps = maxcomp, setSeed = 1000)
PLSnested_CTQ[[1]]


# create permutation samples [takes ~2h on my machine]
# plsPerm_all_CTQ <- permutePLSnestedCV(df_outcomes$ctq_sum, df_CpG_m, nrepeats = repeats, nfolds = folds, nperms = perms)
# write.csv(plsPerm_all_CTQ, here("Results", "plsPerm_all_CTQ.csv"))

# p-value
df_CTQ_perms <- read.csv(here("Results", "plsPerm_all_CTQ.csv"))
p_CTQ <- sum(df_CTQ_perms[,2] >= PLSnested_CTQ[[1]])/length(df_CTQ_perms[,2])
p_CTQ

# plot ncomp solutions for individual folds
df_CTQ_optComp <- melt(PLSnested_CTQ[[3]])
df_CTQ_optComp$ncomp <- rep(c(1:maxcomp), times = repeats*folds)
df_CTQ_optComp_min <- ddply(df_CTQ_optComp, "variable", subset, value == min(value))
ggplot(data = df_CTQ_optComp, aes(y = value, x = ncomp, group = variable)) +
  geom_line(colour = "darkgrey") + geom_point(size = 1, colour = "darkgrey") +
  
  geom_point(data = df_CTQ_optComp_min, colour = "darkblue") + 
  
  theme_classic() + ylab("RMSE (cross-validated)") + xlab("number of components")

ggsave(here("figures", "numberComponents_CTQ.png"), device = "png")

# create final model with optimal number of components
CTQ_finalModel = plsr(DV ~ .,
                           data = PLSnested_CTQ$dat,
                           ncomp = 1)
summary(CTQ_finalModel)


# plot loadings
pls_CTQ_loadings = CTQ_finalModel$loadings[,1:ncol(CTQ_finalModel$loadings)]
dfplot_CTQ_loadings <- data.frame(dfcodes$CpG, pls_CTQ_loadings, dfcodes$segment2)

ggplot(data = dfplot_CTQ_loadings, aes(x = reorder(dfcodes.CpG, rep((1:nrow(dfplot_CTQ_loadings)), length.out = nrow(dfplot_CTQ_loadings))), y = pls_CTQ_loadings, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_path(group = 1, size = 0.8) +
  
  theme_classic() + ylim(-max(abs(pls_CTQ_loadings))*1.5, max(abs(pls_CTQ_loadings))*1.5) 

# plot bivariate correlations
dfplot_CTQ_cors <- dfplot_CTQ_loadings
dfplot_CTQ_cors$pls_CTQ_loadings  <- cor(PLSnested_CTQ$dat)[1, -1]

ggplot(data = dfplot_CTQ_cors, aes(x = reorder(dfcodes.CpG, rep((1:nrow(dfplot_CTQ_cors)), length.out = nrow(dfplot_CTQ_cors))), y = pls_CTQ_loadings, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_path(group = 1, size = 0.8) +
  
  theme_classic() + ylim(-max(abs(pls_CTQ_loadings))*1.5, max(abs(pls_CTQ_loadings))*1.5) 



#################
# childhood trauma [categorical]

# do nested k-fold PLS
PLSnested_CTQcat = PLSnestedCV(df_outcomes$ctq_dich01_multi01, df_CpG_m, nrepeats = repeats, nfolds = folds, maxComps = maxcomp, setSeed = 1000, classification = T)
PLSnested_CTQcat[[1]]

outcome = df_outcomes$ctq_dich01_multi01
predictors = df_CpG_m

# create permutation samples [takes ~2h on my machine]
# plsPerm_all_CTQcat <- permutePLSnestedCV(df_outcomes$ctq_dich01_multi01, df_CpG_m, nrepeats = repeats, nfolds = folds, nperms = perms, classification = T)
# write.csv(plsPerm_all_CTQcat, here("Results", "plsPerm_all_CTQcat.csv"))

# p-value
df_CTQcat_perms <- read.csv(here("Results", "plsPerm_all_CTQcat.csv"))
p_CTQcat <- sum(df_CTQcat_perms[,2] >= PLSnested_CTQcat[[1]])/length(df_CTQcat_perms[,2])
p_CTQcat

# plot ncomp solutions for individual folds
df_CTQcat_optComp <- melt(PLSnested_CTQcat[[3]])
df_CTQcat_optComp$ncomp <- rep(c(1:maxcomp), times = repeats*folds)
df_CTQcat_optComp_min <- ddply(df_CTQcat_optComp, "variable", subset, value == max(value))
ggplot(data = df_CTQcat_optComp, aes(y = value, x = ncomp, group = variable)) +
  geom_line(colour = "darkgrey") + geom_point(size = 1, colour = "darkgrey") +
  
  geom_point(data = df_CTQcat_optComp_min, colour = "darkblue") + 
  
  theme_classic() + ylab("Accuracy (cross-validated)") + xlab("number of components")

ggsave(here("figures", "numberComponents_traumaCat.png"), device = "png")



# create final model with optimal number of components
CTQcat_finalModel = plsda(x = PLSnested_CTQcat$dat[,-1], 
                          y = PLSnested_CTQcat$dat[,1],
                      ncomp = 1)
summary(CTQcat_finalModel)


# plot loadings
pls_CTQcat_loadings = CTQcat_finalModel$loadings[,1:ncol(CTQcat_finalModel$loadings)]
dfplot_CTQcat_loadings <- data.frame(dfcodes$CpG, pls_CTQcat_loadings, dfcodes$segment2)
dfplot_CTQcat_loadings$dfcodes.CpG <- factor(dfplot_CTQcat_loadings$dfcodes.CpG, levels = unique(dfplot_CTQcat_loadings$dfcodes.CpG))
dfplot_CTQcat_loadings$dfcodes.segment2 <- factor(dfplot_CTQcat_loadings$dfcodes.segment2, levels = unique(dfplot_CTQcat_loadings$dfcodes.segment2))

ggplot(data = dfplot_CTQcat_loadings, aes(x = reorder(dfcodes.CpG, rep((1:nrow(dfplot_CTQcat_loadings)), length.out = nrow(dfplot_CTQcat_loadings))), y = pls_CTQcat_loadings, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_path(group = 1, size = 0.8) +
  
  theme_classic() + ylim(-max(abs(pls_CTQcat_loadings))*1.5, max(abs(pls_CTQcat_loadings))*1.5) + ylab("PLS Loadings") + xlab("OXTR CpG Site") +
  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_line(colour = "black"), legend.title = element_blank()) +
  scale_color_brewer(palette="Spectral")

ggsave(here("figures", "PLSloadings_catTrauma.pdf"), device = "pdf")
ggsave(here("figures", "PLSloadings_catTrauma.png"), device = "png")



#################
# combined analyses on geneexpr and CTQ

# combined plot of PLS loadings
names(dfplot_geneExpr_loadings)[2] <- "loadings"
names(dfplot_CTQ_loadings)[2] <- "loadings"
pls_all_loadings_long <- rbind(dfplot_geneExpr_loadings, dfplot_CTQ_loadings)
pls_all_loadings_long$dfcodes.CpG <- factor(pls_all_loadings_long$dfcodes.CpG, levels = unique(pls_all_loadings_long$dfcodes.CpG))
pls_all_loadings_long$dfcodes.segment2 <- factor(pls_all_loadings_long$dfcodes.segment2, levels = unique(pls_all_loadings_long$dfcodes.segment2))
pls_all_loadings_long$outcome <- rep(c("geneExpr", "CTQ"), each = nrow(pls_all_loadings_long)/2)

ggplot(data = pls_all_loadings_long, aes(x = dfcodes.CpG, y = loadings, group = outcome, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_line(size = 0.8) +
  geom_label(x = 25, y = 0.30, label = paste("CTQ:\n", "R² = ", round(PLSnested_CTQ[[1]], 2), "\n", "p = ", p_CTQ), colour = "black") +
  geom_label(x = 30, y = -0.25, label = paste("Gene expression:\n", "R² = ", round(PLSnested_Genexpr[[1]], 2), "\n", "p = ", p_geneExpr), colour = "black") +
  
  theme_classic() + ylim(-max(abs(pls_CTQ_loadings))*1.5, max(abs(pls_CTQ_loadings))*1.5) + ylab("PLS Loadings") + xlab("OXTR CpG Site") + 
  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_line(colour = "black"), legend.title = element_blank()) +
  scale_color_brewer(palette="Spectral")

ggsave("PLSloadings.pdf", device = "pdf")


# correlation between loadings
cor(dfplot_CTQ_loadings$loadings, dfplot_geneExpr_loadings$loadings)
plot(dfplot_CTQ_loadings$loadings, dfplot_geneExpr_loadings$loadings)
cor(dfplot_CTQ_loadings[dfplot_CTQ_loadings$dfcodes.segment2 != "exon3_transl", "loadings"], 
    dfplot_geneExpr_loadings[dfplot_geneExpr_loadings$dfcodes.segment2 != "exon3_transl", "loadings"])

# create df for mediation via PLS components
dfMediation <- data.frame(df_outcomes[!is.na(df_outcomes$Genexp_OXTR_mother), ], 
                          as.numeric(CTQ_finalModel$scores[!is.na(df_outcomes$Genexp_OXTR_mother)]), as.numeric(geneExpr_finalModel$scores))
names(dfMediation)[4:5] <- c("CTQ_CpG_scores", "GeneExpr_CpG_scores") 
dfMediation <- as.data.frame(scale(dfMediation))
cor(dfMediation)

# create df for multivariate mediation in matlab
dfMediation_Mat <- data.frame(dfMediation[,1:3], PLSnested_Genexpr$dat[,-1])
write.csv(dfMediation_Mat, file = "Data/MediationData.csv", row.names = F)




#################
# check pattern with unrelated variable

dfRandPlot <- data.frame(
  matrix(
    ncol = 3,
    nrow = 0
    )
)

for(i in 1:20){
  
  set.seed(i+1000)
  
  randOut <- rnorm(nrow(df_CpG_m))
  dfPLSage <- data.frame(randOut, df_CpG_m)
  
  testModel = plsr(randOut ~ .,
                   data = dfPLSage,
                   ncomp = 1)
  summary(testModel)
  
  pls_testModel_loadings = testModel$loadings[,1:ncol(testModel$loadings)]
  dfplot_testModel_loadings <- data.frame(dfcodes$CpG, pls_testModel_loadings, dfcodes$segment2)
  
  dfRandPlot <- rbind(dfRandPlot, dfplot_testModel_loadings)
  
}

dfRandPlot$iteration <- rep(1:20, each = nrow(dfplot_testModel_loadings))
dfRandPlot$dfcodes.segment2 <- factor(dfRandPlot$dfcodes.segment2, levels = unique(dfRandPlot$dfcodes.segment2))

ggplot(data = dfRandPlot, aes(x = reorder(dfcodes.CpG, rep((1:nrow(dfRandPlot)), length.out = nrow(dfRandPlot))), y = pls_testModel_loadings, group = iteration, colour = factor(dfcodes.segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_line(size = 0.8) +

  theme_classic() + ylim(-max(abs(pls_testModel_loadings))*1.5, max(abs(pls_testModel_loadings))*1.5) + ylab("PLS Loadings") + xlab("OXTR CpG Site") + 
  
  theme(axis.text.x=element_blank(), axis.ticks.x=element_line(colour = "black"), legend.title = element_blank()) +
  scale_color_brewer(palette="Spectral")

ggsave(here("figures", "PLSloadings_random.pdf"), device = "pdf")
ggsave(here("figures", "PLSloadings_random.png"), device = "png")





##################################
# plot multivariate mediation loadings

dfMedResults <- read.table(here("results", "multivarMediation_loadings.txt"), sep = ",", header = T)
dfMedResults_plotCat <- data.frame(dfMedResults[, c("Wcat", "Pcat")], dfcodes)

ggplot(data = dfMedResults_plotCat, aes(x = reorder(CpG, 1:nrow(dfMedResults_plotCat)), y = Wcat, colour = factor(segment2))) +
  geom_hline(yintercept=0, colour = "darkgrey") +
  geom_path(group = 1, size = 0.8) +
  
  theme_classic() + ylim(-0.4, 0.4) 

multiMedfdr <- p.adjust(dfMedResults$Pcat, "fdr")

multiMedSig <- ifelse(multiMedfdr < 0.05, round(multiMedfdr, 3), "")
data.frame(dfMedResults, multiMedSig)




##################################################################
# repeat unsupervised analyses for children

# subset mother CpGs and convert to decimal
df_CpG_k = df[,grepl("CpG_k", names(df))]
df_CpG_k = as.data.frame(apply(df_CpG_k,2, function(x){as.numeric(sub("%", "", x, fixed=TRUE))/100}))
df_CpG_k = df_CpG_k[complete.cases(df_CpG_k),]

# load identifier of gene sections
dfcodes = read.csv(here("Data", "OXTR_segmentCodes.csv"), header = T, sep = ";")
dfcodes = dfcodes[dfcodes$CpG %in% str_replace(names(df_CpG_k), "_k_", "_"),]


#################################
# qualitative description of OXTR

# plot methylation across CpGs to check for sufficient variance
boxplot(df_CpG_k, range = 1.5)

# count/plot number of outliers per CpG (3*IQR)
out_count = as.data.frame(matrix(nrow = ncol(df_CpG_k), ncol = 2, dimnames = list(NULL, c("CpG", "nOut"))))
out_count$CpG = str_replace(names(df_CpG_k), "_k_", "_")
for(i in 1:ncol(df_CpG_k)){out_count[i,2] = length(boxplot.stats(df_CpG_k[,i], coef = 3)$out)}
df_out = merge(out_count, dfcodes, by = "CpG", sort = F)
ggplot(data = df_out, aes(x = reorder(CpG, 1:nrow(df_out)), y = nOut, colour = segment2)) + 
  geom_point()


# correlation plot 
CpG_corr_k = cor(df_CpG_k, use = "pairwise.complete.obs")
colnames(CpG_corr_k) <- NULL
rownames(CpG_corr_k) <- substring(rownames(CpG_corr_k), 7)

pdf(here("Figures", "corrPlot_child.pdf"))
corrplot::corrplot(CpG_corr_k, tl.pos = "d", tl.cex = 0.25, method = "color", type = "upper", tl.col = "black")
dev.off()

# combined corrplot
pdf(here("Figures", "corrPlot_combined.pdf"))
par(mfrow=c(1,2))
corrplot::corrplot(CpG_corr, tl.pos = "d", tl.cex = 0.25, method = "color", type = "upper", tl.col = "black", title = "Mothers")
corrplot::corrplot(CpG_corr_k, tl.pos = "d", tl.cex = 0.25, method = "color", type = "upper", tl.col = "black", title = "Children")
dev.off()


mean(CpG_corr_k[132:182, 132:182])
mean(CpG_corr_k[53:64, 53:64])

# correspondance between mother and child correlation structure
exon3ind <- CpG_corr_k
exon3ind[132:182, 132:182] <- "exon 3 cluster"
exon3ind[53:64, 53:64] <- "intron 1 cluster"
exon3ind <- ifelse(exon3ind == "exon 3 cluster" | exon3ind == "intron 1 cluster", exon3ind, "other section")
CpG_corr_mk <- data.frame(CpG_corr[lower.tri(CpG_corr)], CpG_corr_k[lower.tri(CpG_corr_k)], exon3ind[lower.tri(exon3ind)])
names(CpG_corr_mk) <- c("mother", "child", "geneSection")

cor(CpG_corr_mk$mother, CpG_corr_mk$child)

ggplot(data = CpG_corr_mk, aes(y = child, x = mother, group = geneSection, colour = geneSection)) +
  geom_point() +
  geom_point(data = subset(CpG_corr_mk, geneSection == "intron 1 cluster"))

ggsave(here("Figures", "motherChildCorrespondance.pdf"), device = "pdf")

cor(CpG_corr_mk[CpG_corr_mk$geneSection == "exon 3 cluster", "mother"], CpG_corr_mk[CpG_corr_mk$geneSection == "exon 3 cluster", "child"])
cor(CpG_corr_mk[CpG_corr_mk$geneSection == "intron 1 cluster", "mother"], CpG_corr_mk[CpG_corr_mk$geneSection == "intron 1 cluster", "child"])
cor(CpG_corr_mk[CpG_corr_mk$geneSection == "other section", "mother"], CpG_corr_mk[CpG_corr_mk$geneSection == "other section", "child"])




df_CpGwiseSDs <- as.numeric(read.csv2(here("Data", "OXTRStabwSingleCpGs.csv"))[,3])^2
df_CpGwiseSDs <- c(df_CpGwiseSDs, rep(0.0228^2, times = ncol(df_CpG_k)-length(df_CpGwiseSDs)))
dfchi <- nrow(df_CpG_k)-1
CpGvariancePs <- numeric(ncol(df_CpG_k))

for(i in 1:ncol(df_CpG_k)){
  
  chi_square <- dfchi*var(df_CpG_k[,i])/df_CpGwiseSDs[i]
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
df_CpG_k <- df_CpG_k[, CpGvariancePs <= 0.05]
dfcodes <- dfcodes[CpGvariancePs <= 0.05, ]



#################################
# unsupervised datra-driven description of OXTR

# z score variables
df_CpG_k_z = scale(df_CpG_k)


#################
# factor analysis
KMO(df_CpG_k_z)
fa.parallel(df_CpG_k_z, fa="pc")
pca_all = principal(df_CpG_k, nfactors=7, rotate="oblimin")
pca_all # factor intercorrelations neglegible
pca_all = principal(df_CpG_k, nfactors=7, rotate="varimax")
pca_all$loadings


#################
# cluster analysis (dbscan)

# transpose data
df_CpG_k_z_trans = t(df_CpG_k_z)

# determine eps. try k = 3, k = 5 (default), k = 8, k = 12
dbscan::kNNdistplot(df_CpG_k_z_trans, k =  3)
abline(h = 5, lty = 2) # for k = 3
dbscan::kNNdistplot(df_CpG_k_z_trans, k =  5)
abline(h = 5.8, lty = 2) # for k = 5
dbscan::kNNdistplot(df_CpG_k_z_trans, k =  8)
abline(h = 7.7, lty = 2) # for k = 8
dbscan::kNNdistplot(df_CpG_k_z_trans, k =  12)
abline(h = 8, lty = 2) # for k = 8

# apply clustering
res.fpc3 <- fpc::dbscan(df_CpG_k_z_trans, eps = 5, MinPts = 3)
print(res.fpc3)
res.fpc5 <- fpc::dbscan(df_CpG_k_z_trans, eps = 5.8, MinPts = 5)
print(res.fpc5)
res.fpc8 <- fpc::dbscan(df_CpG_k_z_trans, eps = 7.7, MinPts = 8)
print(res.fpc8)
res.fpc12 <- fpc::dbscan(df_CpG_k_z_trans, eps = 8, MinPts = 12)
print(res.fpc12)

# write.csv(data.frame(res.fpc5$cluster, names(df_CpG_k)), "clusterAssignment.csv") 

# plot
df_clustMem = data.frame(dfcodes$CpG, dfcodes$segment2, res.fpc3$cluster, res.fpc5$cluster, res.fpc8$cluster, res.fpc12$cluster)
df_clustMem_long = melt(df_clustMem, id.vars = c("dfcodes.CpG", "dfcodes.segment2"))

ggplot(data = df_clustMem_long, aes(x = reorder(dfcodes.CpG, rep((1:nrow(df_clustMem)), length.out = nrow(df_clustMem_long))), y = value, colour = dfcodes.segment2, shape = variable)) + 
  geom_jitter(width = 0, height = 0.3) + scale_y_continuous(breaks=c(0, 1, 2)) + scale_x_discrete(breaks=NULL) +
  xlab("CpG Sites") + ylab("Cluster Assignment") + theme_classic() + theme(legend.title = element_blank())
ggsave(here("Figures", "DBSCAN_clustMembership_children.pdf"), device = "pdf")
ggsave(here("Figures", "DBSCAN_clustMembership_children.png"), device = "png")

#################
# hierarchical clustering (hclust)

# dendrogram
dd <- dist(df_CpG_k_z_trans, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

fviz_dend(hc)
fviz_nbclust(df_CpG_k_z_trans, FUN = hcut, method = "wss")
fviz_nbclust(df, FUN = hcut, method = "silhouette")

pdf(here("Figures", "Dendrogram_children.pdf"))
fviz_dend(
  hc,
  k = 3,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.3
)
dev.off()

png(here("Figures", "Dendrogram_children.png"))
fviz_dend(
  hc,
  k = 3,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.3
)
dev.off()


clustAssign_child <- cutree(hc, k = 3)

# proportions of genetic sections represented in the three clusters
df_hclustAssign_child <- data.frame(dfcodes, clustAssign_child) 
round(table(df_hclustAssign_child[df_hclustAssign_child$clustAssign_child == 1, "segment2"])/sum(clustAssign_child == 1), 2)
round(table(df_hclustAssign_child[df_hclustAssign_child$clustAssign_child == 2, "segment2"])/sum(clustAssign_child == 2), 2)
round(table(df_hclustAssign_child[df_hclustAssign_child$clustAssign_child == 3, "segment2"])/sum(clustAssign_child == 3), 2)


df_hclustAssign_child$DBSCANclust <- df_clustMem$res.fpc5.cluster
df_hclustAssign_child$hclustRec <- ifelse(df_hclustAssign_child$clustAssign_child == 1, 
                                    0,
                                    ifelse(df_hclustAssign_child$clustAssign_child == 2, 
                                           1,
                                           2)
)

round(sum(df_hclustAssign_child$DBSCANclust == df_hclustAssign_child$hclustRec)/nrow(df_hclustAssign_child), 2)


# mother-child agreement in clustering solutions
round(sum(df_hclustAssign_child$hclustRec == df_hclustAssign_mothers$hclustRec)/nrow(df_hclustAssign_child), 2)
round(sum(df_hclustAssign_child$DBSCANclust == df_hclustAssign_mothers$DBSCANclust)/nrow(df_hclustAssign_child), 2)

