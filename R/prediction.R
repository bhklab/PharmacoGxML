prediction <- function(models, validation.set, method=c("ridge", "lasso"))
{
  predicted_labels <- NULL
  switch(method, "ridge"={
    for(i in 1:length(models)){
      model <- models[[i]]
      features <- rownames(coef(model))[2:nrow(coef(model))]
      predicted_labels <- cbind(predicted_labels, predict(model, newx=validation.set[,features], s="lambda.min"))
    }
  }, "lasso"={
    model <- models[[i]]
    features <- rownames(coef(model))[2:nrow(coef(model))]
    predicted_labels <- cbind(predicted_labels, predict(model, newx=validation.set[,features], s="lambda.min"))
  })
  return(predicted_labels)
}
