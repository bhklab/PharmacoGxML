lasso <- function(train_inputs, train_labels, folds.no = 5){
  model <- cv.glmnet(train_inputs,
                     train_labels,
                     alpha=1,
                     nfolds=folds.no)
  return(model)
}
