random_forest <- function(train_inputs, train_labels, folds.no = 5, sampling.no = 10, trees.no= 30){
  if(missing(trees.no)){
    trees.no <- seq(10, 30, 10)
  }
  models <- list()
  R2 <- NULL
  for(t in trees.no){
    model <- caret::train(train_inputs,
                          train_labels,
                          method="rf",
                          ntree=t,
                          trControl=trainControl(method="repeatedcv",
                                                 number=folds.no,
                                                 repeats=sampling.no,
                                                 search="grid"))
    rr <- model$results
    R2 <- c(R2, rr[which(rr$mtry == model$bestTune[[1]]), "Rsquared"])
    models[[t]] <- model
  }
  return(models[[trees.no[which.max(R2)]]])
}
