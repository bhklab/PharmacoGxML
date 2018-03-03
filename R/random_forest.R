random_forest <- function(train_inputs, train_labels, folds.no, trees.no){
  model <- caret::train(train_inputs,
                        train_labels,
                        method="rf",
                        ntree=trees.no,
                        trControl=trainControl(method="cv",
                                               number=folds.no,
                                               search="grid"))
  return(model)
}
