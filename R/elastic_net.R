elastic_net <- function(train_inputs, train_labels, folds.no, alpha){
  model <- cv.glmnet(train_inputs,
                     train_labels,
                     alpha=alpha,
                     nfolds=folds.no)
  return(model)
}
