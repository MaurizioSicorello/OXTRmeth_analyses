
# nested k-fold CV evaluation of PLS
PLSnestedCV <- function(outcome, predictors, nrepeats, nfolds, maxComps = 30, setSeed = 1000){
  
  fitControl <- trainControl(method = "repeatedcv",   
                             number = nfolds, 
                             repeats = nrepeats)
  
  # prepare input df
  df_pred = data.frame(outcome, predictors)
  df_pred = scale(df_pred)
  df_pred = as.data.frame(df_pred[!is.na(df_pred[,1]), ])
  names(df_pred)[1] = "DV"
  
  # prepare output df
  df_out <- as.data.frame(matrix(nrow = nrepeats*nfolds, ncol = 2))
  df_foldwise <- as.data.frame(matrix(nrow = maxComps, ncol = nrepeats*nfolds))
  plsGrid <- expand.grid(ncomp = seq(1, maxComps))
  count = 0
  
  # main loop
  for(i in 1:nrepeats){
    
    set.seed(setSeed+count)
    plsFolds <- createFolds(df_pred[,1], k = nfolds)
    
    for(j in 1:nfolds){
      
      count = count+1
      
      trainSet <- df_pred[-plsFolds[[j]],]
      testSet <- df_pred[plsFolds[[j]],]
      
      f <- as.formula("DV ~ .")
      
      set.seed(setSeed+count)
      plsOptComp <- train(f,
                          data = trainSet,
                          trControl = fitControl,
                          method = "pls",
                          metric = "RMSE",
                          preProcess = c('scale', 'center'),
                          na.action = na.omit,
                          tuneGrid = plsGrid)
      
      Rsquared <- 1 - sum((testSet[,1] - predict(plsOptComp, newdata = testSet))^2)/(var(testSet[,1])*(nrow(testSet)-1))
      
      df_out[count, 1] <- max(Rsquared, 0)
      df_out[count, 2] <- plsOptComp$bestTune$ncomp
      df_foldwise[, count] <- plsOptComp$results$RMSE
      
    }
  }
  return(list(Rsquared = mean(df_out[,1]), outerFolds = df_out, innerFolds = df_foldwise, dat = df_pred))
}



# permutation test of PLS
permutePLSnestedCV <- function(outcome, predictors, nrepeats, nfolds, nperms, maxComps = 30, setSeed = 1000){
  
  savePerm <- numeric(nperms)
  
  for(i in 1:nperms){
    
    outcomePerm <- sample(outcome)
    savePerm[i] <- PLSnestedCV(outcomePerm, predictors, nrepeats, nfolds, maxComps, setSeed = round(runif(1, 1000, 4000)))[[1]]
    
  }
  
  return(savePerm)
}
