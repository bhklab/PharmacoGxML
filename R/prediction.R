prediction <- function(model,
                       method=c("ridge", "lasso", "random_forest", "svm"),
                       test.set)
{
  predicted_labels <- NULL

  if(method == "ridge" | method == "lasso" | method == "elastic_net"){
    features <- rownames(coef(model))[2:nrow(coef(model))]
    if(length(intersect(colnames(test.set), features)) != length(features)){
      return(rep(NA, nrow(test.set)))
    }
    predicted_labels <- cbind(predicted_labels, predict(model, newx=test.set[,features], s="lambda.min"))
  }else{
    features <- colnames(model$trainingData)[1:ncol(model$trainingData) -1 ]
    predicted_labels <- cbind(predicted_labels, predict(model, newdata=test.set[,features]))
  }
  return(predicted_labels)
}
