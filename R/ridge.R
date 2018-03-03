ridge <- function(train_inputs, train_labels, folds.no){
  model <- cv.glmnet(train_inputs,
                     train_labels,
                     alpha=0,
                     nfolds=folds.no)
  return(model)
}
