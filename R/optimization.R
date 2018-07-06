library(mCI)
library(glmnet)
library(Hmisc)
optimization <- function(train,
                         labels,
                         folds.no=7,
                         sampling.no=1,
                         features.no=100,
                         method=c("ridge", "lasso", "random_forest", "svm"),
                         feature.selection=c("mRMR", "variance"),
                         assessment=c("corr", "CI", "mCI", "r_squared"),
                         result.path){
  performance <- list()
  models <- list()
  predictions_range <- NULL
  for (drug in 1:nrow(labels))
  {
    drug_name <- rownames(labels)[drug]
    performance[[drug_name]] <- list()
    all_predicted <- NULL
    all_valid_labels <- NULL
    output <- NULL
    x <- train
    y <- labels[drug, ]

    # Remove NA's
    toRemove <- which(is.na(y))
    if (length(toRemove != 0))
    {
      y <- y[-toRemove]
      x <- x[-toRemove, ]
    }

    # Cross validation
    # par(mfrow = c(2,5))

    for(s in 1:sampling.no){
      if(sampling.no == 1){
        order_of_labels <- 1:length(y)
      }else{
        order_of_labels <- sample(1:length(y))
      }
      # cuts <- Hmisc::cut2(order_of_labels, m=1, g=folds.no, onlycuts= TRUE)
      cuts <- Hmisc::cut2(order_of_labels, m=1, g=folds.no)

      for (i in seq_along(unique(cuts)))
      {
        # Split the data into a training set and validation set
        
        train_inputs <- x[-order_of_labels[as.numeric(cuts) == i], , drop=F]
        train_labels <- y[-order_of_labels[as.numeric(cuts) == i]]
        valid_inputs <- x[order_of_labels[as.numeric(cuts) == i], , drop=F]
        valid_labels <- y[order_of_labels[as.numeric(cuts) == i]]

        ##Feature selection
        features <- PharmacoGxML::featureSelection(train_inputs, train_labels, method=feature.selection, features.no=features.no)
        train_inputs <- train_inputs[, features, drop=F]
        valid_inputs <- valid_inputs[, features, drop=F]
        #######
        current_output <- NULL

        # Get the cross validated elastic net fit
        #alpha=1 is lasso, alph=0 is ridge
        cvfit <- do.call(method, c(train_inputs, train_labels))
        predicted_labels <- make_predictions(cvfit, newdata=valid_inputs)
        
        # print(sprintf("%s, %s, method: %s,   fold#: %s  sampling#: %s", drug, rownames(labels)[drug], method, i, s))

        if(length(predicted_labels) == length(valid_labels)){
          all_predicted <- c(all_predicted, predicted_labels)
          all_valid_labels <- c(all_valid_labels, valid_labels)
        }
        #plot(predicted_labels, valid_labels, main = paste("Fold", i))
      }
      fit <- lm(all_valid_labels ~ all_predicted)
      slope <- fit$coefficients[[2]]
      intercept <- fit$coefficients[[1]]

      plot(all_valid_labels,
           all_predicted,
           main=sprintf("%s\nmethod:%s", drug_name, method),
           cex.main=1, ylab="Predictions", xlab="drug sensitivity", pch=20, col="gray40")

      intercept <- round(intercept, 2)
      slope <- round(slope, 2)
      equation <- paste("y = ", slope, "x + ", intercept, sep = "")
      legend_label <- NULL
      if(method == "ridge" | method == "lasso" | method == "elastic_net"){
        abline(intercept, slope, lty=2)
        legend_label <- c(legend_label, equation)
      }
      if("corr" %in% assessment){
        corr <- round(cor(all_predicted, all_valid_labels, use="pairwise.complete.obs", method = "pearson"), digits=2)
        legend_label <- c(legend_label, sprintf("r=%s", corr))
        performance[[drug_name]][["corr"]] <- c(performance[[drug_name]][["corr"]], corr)
      }
      if("mCI" %in% assessment){
        mci <- round(mCI::paired.concordance.index(all_predicted, all_valid_labels, delta.pred=0, delta.obs=.2)$cindex, digits=2)
        legend_label <- c(legend_label, sprintf("mCI=%s", mci))
        performance[[drug_name]][["mCI"]] <- c(performance[[drug_name]][["mCI"]], mci)
      }
      if("CI" %in% assessment){
        ci <- round(mCI::paired.concordance.index(all_predicted, all_valid_labels, delta.pred=0, delta.obs=0)$cindex, digits=2)
        legend_label <- c(legend_label, sprintf("CI=%s", ci))
        performance[[drug_name]][["CI"]] <- c(performance[[drug_name]][["CI"]], ci)
      }

      if("r_squared" %in% assessment){
        r2 <- round((all_predicted - all_valid_labels) ^ 2, digits=2)
        legend_label <- c(legend_label, sprintf("R2=%s", r2))
        performance[[drug_name]][["r_squared"]] <- c(performance[[drug_name]][["r_squared"]], r2)
      }
      legend("topright",
             legend=paste(legend_label, sep="\n"),
             bty="n")
      predictions_range <- rbind(predictions_range, range(all_predicted, na.rm=T))
    }
    ##build the model for predicting future cases
    features <- featureSelection(x, y, method=feature.selection, features.no=features.no)
    train_set <- x[, features, drop=F]

    switch(method,
           "ridge"={
             model <- ridge(train_set, y)
           },
           "lasso"={
             model <- lasso(train_set, y)
           },
           "random_forest"={
             model <- random_forest(train_set, y)
           },
           "svm"={
             model <- svm(train_set, y)
           })
    models[[drug_name]] <- model

  }
  if(!missing(result.path)){
    save(predictions_range, file=sprintf("%s/prediction_ranges.RData", dir))
  }

  return(list("performance"=performance,  "model"=models))
}

