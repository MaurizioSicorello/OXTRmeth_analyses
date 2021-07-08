# analysen wiederholen für MT2. ggf. plots für genexpression. Leons Segmentmethode. change code to show RMSE plots

# Benennung der CpGs mit Robert klären

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
df_outcomes = df[complete.cases(df_CpG_m), c("ctq_sum", "Genexp_OXTR_mother")]
df_CpG_m = df_CpG_m[complete.cases(df_CpG_m),]

# load identifier of gene sections
dfcodes = read.csv(here("Data", "OXTR_segmentCodes.csv"), header = T, sep = ";")
dfcodes = dfcodes[dfcodes$CpG %in% str_replace(names(df_CpG_m), "_m_", "_"),]


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


df_CpGwiseSDs <- as.numeric(read.csv2(here("Data", "OXTRStabwSingleCpGs.csv"))[,3])^2
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

# correlation plot 
CpG_corr = cor(df_CpG_m, use = "pairwise.complete.obs")
colnames(CpG_corr) <- NULL
rownames(CpG_corr) <- substring(rownames(CpG_corr), 7)

pdf(here("Figures", "corrPlot.pdf"))
corrplot::corrplot(CpG_corr, tl.pos = "d", tl.cex = 0.25, method = "color", type = "upper", tl.col = "black")
dev.off()


#################################
# unsupervised datra-driven description of OXTR

# z score variables
df_CpG_m_z = scale(df_CpG_m)


#################
# factor analysis
KMO(df_CpG_m_z)
fa.parallel(df_CpG_m_z, fa="pc")
pca_all = principal(df_CpG_m, nfactors=7, rotate="oblimin")
pca_all # factor intercorrelations neglegible
pca_all = principal(df_CpG_m, nfactors=7, rotate="varimax")
pca_all$loadings


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


clustAssign <- cutree(hc, k = 3)

df_hclustAssign <- data.frame(dfcodes, clustAssign) 
round(table(df_hclustAssign[df_hclustAssign$clustAssign == 1, "segment2"])/sum(clustAssign == 1), 2)
round(table(df_hclustAssign[df_hclustAssign$clustAssign == 2, "segment2"])/sum(clustAssign == 2), 2)
round(table(df_hclustAssign[df_hclustAssign$clustAssign == 3, "segment2"])/sum(clustAssign == 3), 2)


#################################
# supervised data-driven description of OXTR

#################
# ML settings
maxcomp = 30
repeats = 2
repeats2 = 10
folds = 5
perms = 1000

fitControlFinal <- trainControl(method = "repeatedcv", number = folds, repeats = 10)
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
# childhood trauma

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
dfMediation <- data.frame(df_outcomes[!is.na(df_outcomes$Genexp_OXTR_mother), ], as.numeric(CTQ_finalModel$scores[!is.na(df_outcomes$Genexp_OXTR_mother)]), as.numeric(geneExpr_finalModel$scores))
names(dfMediation)[3:4] <- c("CTQ_CpG_scores", "GeneExpr_CpG_scores") 
dfMediation <- as.data.frame(scale(dfMediation))
cor(dfMediation)

# create df for multivariate mediation in matlab
dfMediation_Mat <- data.frame(dfMediation[,1:2], PLSnested_Genexpr$dat[,-1])
write.csv(dfMediation_Mat, file = "MediationData.csv", row.names = F)


library("lavaan")

Medmodel=
  "
  #Regressions
  CTQ_CpG_scores ~ a*ctq_sum
  GeneExpr_CpG_scores ~ b*CTQ_CpG_scores + c*ctq_sum
  Genexp_OXTR_mother ~ d*GeneExpr_CpG_scores + e*CTQ_CpG_scores + f*ctq_sum
  
  #Defined Parameters:
  full_indirect := a*b*d
  first_indirect := a*e
  second_indirect := c*d
  c_prime := f
  c_fullEst := (a*b*d)+(a*e)+(c*d) + f
"

fit=sem(Medmodel,dfMediation)
summary(fit)


Medmodel=
  "
  #Regressions
  GeneExpr_CpG_scores ~ a*ctq_sum
  Genexp_OXTR_mother ~ b*GeneExpr_CpG_scores + c*ctq_sum
  
  #Defined Parameters:
  ie := a*b
  de := c
"

fit=sem(Medmodel,dfMediation)
summary(fit)










##################################
# PLS MT2

#################
# gene expression
df_MT2 <- df_CpG_m[, which(names(df_CpG_m) == "CpG_m_27"):which(names(df_CpG_m) == "CpG_m_53")]
df_predExpr_MT2 = data.frame(df_outcomes$Genexp_OXTR_mother, df_MT2)
names(df_predExpr_MT2)[1] <- "geneExpr"
df_predExpr_MT2[,1] <- scale(df_predExpr_MT2[,1], scale = T, center = T)
df_predExpr_MT2 <- df_predExpr_MT2[!is.na(df_predExpr_MT2$geneExpr), ]


GeneExpr_out <- as.data.frame(matrix(nrow = nfolds*repeats, ncol = 2))
GeneExpr_out_foldwise <- as.data.frame(matrix(nrow = maxcomp, ncol = repeats*nfolds))
plsGrid <- expand.grid(ncomp = seq(1, maxcomp))
count = 0
seed = 1000

for(i in 1:repeats){
  
  set.seed(seed+count)
  plsFolds <- createFolds(df_predExpr_MT2$geneExpr, k = nfolds)
  
  for(j in 1:nfolds){
    
    count = count+1
    
    trainSet <- df_predExpr_MT2[-plsFolds[[j]],]
    testSet <- df_predExpr_MT2[plsFolds[[j]],]
    
    set.seed(seed+count)
    plsOptComp <- train(geneExpr ~ ., 
                        data = trainSet, 
                        trControl = fitControl, 
                        method = "pls", 
                        metric = "RMSE",
                        preProcess = c('scale', 'center'), 
                        na.action = na.omit,
                        tuneGrid = plsGrid)
    
    Rsquared <- 1 - sum((testSet$geneExpr - predict(plsOptComp, newdata = testSet))^2)/(var(testSet$geneExpr)*(nrow(testSet)-1))
    
    GeneExpr_out[count, 1] <- max(Rsquared, 0)
    GeneExpr_out[count, 2] <- plsOptComp$bestTune$ncomp
    GeneExpr_out_foldwise[, count] <- plsOptComp$results$Rsquared
    
  }
}

mean(GeneExpr_out[,1])


##########
# permutation test

savePerm <- numeric(nperm)

t1 <- Sys.time()

for(h in 1:nperm){
  
  dfpredExpr_perm <- df_predExpr
  dfpredExpr_perm$geneExpr <- sample(dfpredExpr_perm$geneExpr)
  
  (df_predExpr$geneExpr)
  
  GeneExpr_out <- as.data.frame(matrix(nrow = nfolds*repeats, ncol = 2))
  GeneExpr_out_foldwise <- as.data.frame(matrix(nrow = maxcomp, ncol = repeats*nfolds))
  plsGrid <- expand.grid(ncomp = seq(1, maxcomp))
  count = 0
  seed = 1000
  
  for(i in 1:repeats){
    
    set.seed(seed+count)
    plsFolds <- createFolds(dfpredExpr_perm$geneExpr, k = nfolds)
    
    for(j in 1:nfolds){
      
      count = count+1
      
      trainSet <- dfpredExpr_perm[-plsFolds[[j]],]
      testSet <- dfpredExpr_perm[plsFolds[[j]],]
      
      set.seed(seed+count)
      plsOptComp <- train(geneExpr ~ ., 
                          data = trainSet, 
                          trControl = fitControl, 
                          method = "pls", 
                          metric = "Rsquared",
                          preProcess = c('scale', 'center'), 
                          na.action = na.omit,
                          tuneGrid = plsGrid)
      
      Rsquared <- 1 - sum((predict(plsOptComp, newdata = testSet) - testSet$geneExpr)^2)/(var(testSet$geneExpr)*nrow(testSet))
      GeneExpr_out[count, 1] <- max(Rsquared, 0)
      GeneExpr_out[count, 2] <- plsOptComp$bestTune$ncomp
      GeneExpr_out_foldwise[, count] <- plsOptComp$results$Rsquared
      
    }
  }
  savePerm[h] <- mean(GeneExpr_out[,1], na.rm = T)
}  

write.csv(savePerm, "plsPerm_all.csv")

t2 <- Sys.time()
t2-t1





t = plsr(geneExpr ~ ., data = trainSet, validation = "CV", segments = 10, ncomp = 7)
Rsquared <- 1 - sum((predict(t, ncomp = 7, newdata = testSet)[1:14] - testSet$geneExpr)^2)/(var(testSet$geneExpr)*nrow(testSet))


p = predict(t, ncomp = 7, newdata = testSet)[1:14]



plsOptComp
plot(plsOptComp)
plsOptComp$results










# sparse partial least squares
splsGrid <- expand.grid(ncomp = seq(1, length(names(df_dum))-1))

spls <- train(geneExpr ~. , 
              data = df_predExpr, 
              trControl = fitControl, 
              method = "spls", 
              preProcess = c('scale', 'center'), 
              na.action = na.omit)

spls$results[spls$results$ncomp == as.numeric(spls$bestTune), ]
plot(spls, metric = "Rsquared")










### profilanalyse




plot(x = names(df_CpG_m), y = ifelse(CpGvariancePs > 0.05, 1, 0))


t = apply(df_CpG_m, 2, var)
sum(t <= meanSens^2)

varTestFun(df_CpG_m[,1], meanSens^2)


varTest(df_CpG_m[,1], alternative = "greater", conf.level = 0.95, 
        sigma.squared = meanSens^2)

((nrow(df_CpG_m)-1)*var(df_CpG_m[,1]))/meanSens^2






# corrmatrix total
names(df_CpG_m) <- substr(names(df_CpG_m), 7, 12)
CpG_corr = cor(df_CpG_m, use = "pairwise.complete.obs")    
corrplot(CpG_corr, tl.pos = "td", tl.cex = 0.5, method = "color", type = "upper")


CpG_corr = cor(df_CpG_m[, which(names(df_CpG_m) == "143"):which(names(df_CpG_m) == "187")], use = "pairwise.complete.obs")    
corrplot.mixed(CpG_corr)








# methylation profiles: https://bioconductor.org/packages/release/bioc/vignettes/BPRMeth/inst/doc/BPRMeth_vignette.html

# Computational and Statistical Epigenomics

# mclust full
df_CpG_m_trans = t(df_CpG_m)
BIC <- mclustBIC(df_CpG_m_trans)
plot(BIC)
gmm <- Mclust(df_CpG_m_trans, G = 1:30, x = BIC)
summary(gmm)
plot(gmm, what = "BIC")

# mclust sufficient variance on CpGs

df_CpG_m_sens = df_CpG_m[,apply(df_CpG_m, 2, sd) > meanSens]

df_CpG_m_sens_trans = t(df_CpG_m_sens)
BIC <- mclustBIC(df_CpG_m_trans)
plot(BIC)
gmm <- Mclust(df_CpG_m_trans, G = 1:30, x = BIC)
summary(gmm)
plot(gmm, what = "BIC")



# PCA

exon3transl = df_CpG_m_sens[,which(names(df_CpG_m_sens) == "125"):which(names(df_CpG_m_sens) == "193")]

fa.parallel(exon3transl, fa="pc")
KMO(exon3transl)

corrplot(cor(df_CpG_m_sens[,which(names(df_CpG_m_sens) == "143"):which(names(df_CpG_m_sens) == "193")], use = "pairwise.complete.obs"))





# USE BORUTA METHOD FOR FEATURE SELECTION?!

# random forest
df_rf_m = data.frame(df$Genexp_OXTR_mother, df_CpG_m[which(names(df_CpG_m) == "E1"):which(names(df_CpG_m) == "E21")])
df_rf_m = df_rf_m[!is.na(df_rf_m$df.Genexp_OXTR_mother), ]

rf_model = cforest(df.Genexp_OXTR_mother ~ ., data = df_rf_m)
correl <- cor(predict(rf_model, OOB = T), df_rf_m$df.Genexp_OXTR_mother)
R_squared <- if(correl < 0){0}else{correl^2}



corrplot.mixed(cor(df_rf_m, use = "pairwise.complete.obs"))



#################################################################
######### CODE TESTING

# simulate data

N = 1000
M = 168
df_test <- matrix(rnorm(N*M,mean=0,sd=1), N, M)
weights <- runif(M, -0.5, 0.5)

DV <- df_test%*%weights
DV <- scale(DV) + rnorm(N, sd = sqrt(5))

df_test <- data.frame(DV, df_test)

summary(lm(DV~., data = df_test))

test_out <- as.data.frame(matrix(nrow = nfolds*repeats, ncol = 2))
test_out_foldwise <- as.data.frame(matrix(nrow = maxcomp, ncol = repeats*nfolds))
plsGrid <- expand.grid(ncomp = seq(1, maxcomp))
count = 0

for(i in 1:repeats){
  
  plsFolds <- createFolds(df_test$DV, k = nfolds)
  
  for(j in 1:nfolds){
    
    count = count+1
    
    trainSet <- df_test[-plsFolds[[j]],]
    testSet <- df_test[plsFolds[[j]],]
    
    set.seed(seed+count)
    plsOptComp <- train(DV ~ ., 
                        data = trainSet, 
                        trControl = fitControl, 
                        method = "pls", 
                        metric = "RMSE",
                        preProcess = c('scale', 'center'), 
                        na.action = na.omit,
                        tuneGrid = plsGrid)
    
    Rsquared <- 1 - sum((testSet$DV - predict(plsOptComp, newdata = testSet))^2)/(var(testSet$DV)*(nrow(testSet)-1))
    
    test_out[count, 1] <- max(Rsquared, 0)
    test_out[count, 2] <- plsOptComp$bestTune$ncomp
    test_out_foldwise[, count] <- plsOptComp$results$Rsquared
    
  }
}

mean(test_out[,1])

