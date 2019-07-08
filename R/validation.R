validation <- function(model,
                       validation.set,
                       validation.labels,
                       method=c("ridge", "lasso", "random_forest", "svm"),
                       assessment=c("corr", "CI", "rCI"))
{
  predicted_labels <- NULL

  if(method == "ridge" | method == "lasso" | method == "elastic_net"){
    features <- rownames(coef(model))[2:nrow(coef(model))]
    predicted_labels <- cbind(predicted_labels, predict(model, newx=validation.set[,features], s="lambda.min"))
  }else{
    features <- colnames(model$trainingData)[1:ncol(model$trainingData) -1 ]
    predicted_labels <- cbind(predicted_labels, predict(model, newdata=validation.set[,features]))
  }
  #predicted_labels <- apply(predicted_labels, MARGIN=1, mean)
  cells <- intersect(names(validation.labels), rownames(predicted_labels))
  validation.labels <- validation.labels[cells]
  predicted_labels <- predicted_labels[cells, 1]

  plot(validation.labels,
       predicted_labels,
       main=sprintf("Validation\n%s\nmethod:%s", "lapatinib", method),
       cex.main=1, ylab="Predictions", xlab="drug sensitivity", pch=20, col="gray40",
       xlim=c(0, 1), ylim=c(0, 1))
  fit <- lm(predicted_labels ~ validation.labels)
  slope <- fit$coefficients[[2]]
  intercept <- fit$coefficients[[1]]

  intercept <- round(intercept, 2)
  slope <- round(slope, 2)
  equation <- paste("y = ", slope, "x + ", intercept, sep = "")

  legend_label <- NULL
  if(method == "ridge" | method == "lasso" | method == "elastic_net"){
    abline(intercept, slope, lty=2)
    legend_label <- c(legend_label, equation)
  }
  line_no <- -3
  switch(assessment,
         "rCI"={
           validation_mci <- paired.concordance.index(predictions=predicted_labels,
                                                      observations=validation.labels,
                                                      delta.pred=0,
                                                      delta.obs=.2,
                                                      alternative="greater",
                                                      logic.operator="or")$cindex
           legend_label <- c(legend_label, sprintf("rCI=%s", round(validation_mci, digits=2)))
         })
  legend("topright",
         legend=paste(legend_label, sep="\n"),
         bty="n")
  return(predicted_labels)
}
