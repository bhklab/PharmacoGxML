svm <- function(train_inputs, train_labels, folds.no = 5, sampling.no = 10, kernel=c("Linear", "Radial")){
  if(missing(kernel)){
    kernel <- c("Linear", "Radial")
  }
  models <- list()
  R2 <- NULL
  for(k in kernel){
    model <- caret::train(train_inputs,
                          train_labels,
                          method=paste0("svm", k),
                          trControl=trainControl(method="repeatedcv",
                                                 number=folds.no,
                                                 repeats=sampling.no,
                                                 search="grid"))
    rr <- model$results
    if(k == "Radial"){
      R2 <- c(R2, rr[which(rr$sigma == model$bestTune[[1]] &
                             rr$C == model$bestTune[[2]]), "Rsquared"])
    }else{
      R2 <- c(R2, rr[which(rr$C == model$bestTune[[1]]), "Rsquared"])
    }
    models[[k]] <- model
  }
  return(models[[kernel[which.max(R2)]]])
}
