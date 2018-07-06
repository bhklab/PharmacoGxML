ridge <- function(train_inputs, train_labels, folds.no = 5){
  model <- cv.glmnet(train_inputs,
                     train_labels,
                     alpha=0,
                     nfolds=folds.no)
  return(model)
}
