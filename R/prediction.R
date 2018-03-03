prediction <- function(models,
                       validation.set,
                       validation.labels,
                       method=c("ridge", "lasso"),
                       assessment=c("corr", "CI", "mCI"))
{
  predicted_labels <- NULL
  for(i in 1:length(models)){
    model <- models[[i]]
    features <- rownames(coef(model))[2:nrow(coef(model))]

    if(method == "ridge" | method == "lasso" | method == "elastic_net"){
      predicted_labels <- cbind(predicted_labels, predict(model, newx=validation.set[,features], s="lambda.min"))
    }else{
      predicted_labels <- cbind(predicted_labels, predict(model, newx=validation.set[,features]))
    }
  }
  predicted_labels <- apply(predicted_labels, MARGIN=1, mean)
  cells <- intersect(names(validation.labels), names(predicted_labels))
  validation.labels <- validation.labels[cells]
  predicted_labels <- predicted_labels[cells]

  plot(validation.labels,
       predicted_labels,
       main=sprintf("Validationc\n%s\nmethod:%s", "lapatinib", method),
       cex.main=1, ylab="Predictions", xlab="drug sensitivity", pch=20, col="gray40")
  fit <- lm(predicted_labels ~ validation.labels)
  slope <- fit$coefficients[[2]]
  intercept <- fit$coefficients[[1]]

  abline(intercept, slope, lty=2)
  intercept <- round(intercept, 2)
  slope <- round(slope, 2)
  equation <- paste("y = ", slope, "x + ", intercept, sep = "")

  mtext(equation, 3, line=-2)
  line_no <- -3
  switch(assessment,
         "mCI"={
           validation_mci <- paired.concordance.index(predictions=predicted_labels,
                                                      observations=validation.labels,
                                                      delta.pred=0,
                                                      delta.obs=.2,
                                                      alternative="greater",
                                                      logic.operator="or")$cindex
           mtext(sprintf("mCI=%s", validation_mci), 3, line=line_no)
           line_no <- line_no - 1
         })
  return(predicted_labels)
}
