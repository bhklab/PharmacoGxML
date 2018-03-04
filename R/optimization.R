library(mCI)
library(glmnet)
library(Hmisc)
optimization <- function(train,
                         labels,
                         result.path,
                         folds.no=7,
                         sampling.no=1,
                         features.no=100,
                         method=c("ridge", "lasso"),
                         feature.selection=c("mRMR", "variance"),
                         assessment=c("corr", "CI", "mCI", "r_squared")){
  performance <- list()
  models <- list()
  if(!missing(result.path)){
    dir <- sprintf("%s/alpha_%s", result.path, gsub("[.]", "_", as.character(alpha)))
    if(!dir.exists(dir)){
      dir.create(dir)
    }
    pdf(sprintf("%s/glmnet_drugs_alpha_%s.pdf", result.path, alpha), width=12, height=120)
    #par(mfrow=c(30, 3))
  }
  predictions_range <- NULL
  for (drug in 1:nrow(labels))
  {
    drug_name <- rownames(labels)[drug]
    performance[[drug_name]] <- list()
    all_predicted <- NULL
    all_valid_labels <- NULL
    output <- NULL
    x <- train
    y <- as.matrix(labels[drug, ])

    # Remove NA's
    toRemove <- which(is.na(y))
    if (length(toRemove != 0))
    {
      y <- as.matrix(y[-toRemove, ])
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
      cuts <- Hmisc::cut2(order_of_labels, m=1, g=folds.no, onlycuts= TRUE)
      for (i in 1:(length(cuts)-1))
      {
        # Split the data into a training set and validation set
        start <- cuts[i]
        end <- cuts[i + 1] - 1
        if (i == (length(cuts) - 1))
        {
          end <- cuts[i + 1]
        }


        train_inputs <- x[-order_of_labels[start:end], , drop=F]
        train_labels <- y[-order_of_labels[start:end]]
        valid_inputs <- x[order_of_labels[start:end], , drop=F]
        valid_labels <- y[order_of_labels[start:end]]

        ##Feature selection
        features <- featureSelection(train_inputs, train_labels, method=feature.selection, features.no=features.no)
        train_inputs <- train_inputs[, features, drop=F]
        valid_inputs <- valid_inputs[, features, drop=F]
        #######
        current_output <- NULL

        # Get the cross validated elastic net fit
        #alpha=1 is lasso, alph=0 is ridge
        switch(method,
              "ridge"={
                cvfit <- ridge(train_inputs, train_labels, folds.no=5)
                predicted_labels <- predict(cvfit, newx=valid_inputs, s="lambda.min")
              },
              "lasso"={
                cvfit <- lasso(train_inputs, train_labels, folds.no=5)
                predicted_labels <- predict(cvfit, newx=valid_inputs, s="lambda.min")
              },
              "random_forest"={
                cvfit <- random_forest(train_inputs, train_labels, folds.no=5, sampling.no=10, trees.no=30)
                predicted_labels <- predict(cvfit, newx=valid_inputs)
              },
              "svm"={
                cvfit <- svm(train_inputs, train_labels, folds.no=5, sampling.no=10)
                predicted_labels <- predict(cvfit, newx=valid_inputs)
              })
        model[[drug_name]][[sprintf("sampling_%s",s)]][[sprintf("fold_%s",i)]] <- cvfit
        print(sprintf("%s, %s, method: %s,   fold#: %s  sampling#: %s", drug, rownames(labels)[drug], method, i, s))

        all_predicted <- c(all_predicted, predicted_labels)
        all_valid_labels <- c(all_valid_labels, valid_labels)
        #plot(predicted_labels, valid_labels, main = paste("Fold", i))
      }
      fit <- lm(all_valid_labels ~ all_predicted)
      slope <- fit$coefficients[[2]]
      intercept <- fit$coefficients[[1]]

      plot(all_valid_labels,
           all_predicted,
           main=sprintf("%s\nmethod:%s", drug_name, method),
           cex.main=1, ylab="Predictions", xlab="drug sensitivity", pch=20, col="gray40")

      abline(intercept, slope, lty=2)
      intercept <- round(intercept, 2)
      slope <- round(slope, 2)
      equation <- paste("y = ", slope, "x + ", intercept, sep = "")

      mtext(equation, 3, line=-2)
      line_no <- -3
      if("corr" %in% assessment){
        corr <- round(cor(all_predicted, all_valid_labels, use="pairwise.complete.obs", method = "pearson"), digits=2)
        mtext(sprintf("r=%s", corr), 3, line=line_no)
        performance[[drug_name]][["corr"]] <- c(performance[[drug_name]][["corr"]], corr)
        line_no <- line_no - 1
      }
      if("mCI" %in% assessment){
        mci <- round(mCI::paired.concordance.index(all_predicted, all_valid_labels, delta.pred=0, delta.obs=.2)$cindex, digits=2)
        mtext(sprintf("mCI=%s", mci), 3, line=line_no)
        performance[[drug_name]][["mCI"]] <- c(performance[[drug_name]][["mCI"]], mci)
        line_no <- line_no - 1
      }
      if("CI" %in% assessment){
        ci <- round(mCI::paired.concordance.index(all_predicted, all_valid_labels, delta.pred=0, delta.obs=0)$cindex, digits=2)
        mtext(sprintf("CI=%s", ci), 3, line=line_no)
        performance[[drug_name]][["CI"]] <- c(performance[[drug_name]][["CI"]], ci)
        line_no <- line_no - 1
      }

      if("r_squared" %in% assessment){
        r2 <- round((all_predicted - all_valid_labels) ^ 2, digits=2)
        mtext(sprintf("R2=%s", r2), 3, line=line_no)
        performance[[drug_name]][["r_squared"]] <- c(performance[[drug_name]][["r_squared"]], r2)
        line_no <- line_no - 1
      }
      predictions_range <- rbind(predictions_range, range(all_predicted, na.rm=T))
    }
    ##build the model for predicting future cases
    features <- featureSelection(x, y, method=feature.selection, features.no=features.no)
    train_set <- x[, features, drop=F]

    switch(method,
           "ridge"={
             model <- ridge(train_set, y, folds.no=5)
           },
           "lasso"={
             model <- lasso(train_set, y, folds.no=5)
           },
           "random_forest"={
             model <- random_forest(train_set, y, folds.no=5, sampling.no=10, trees.no=30)
           },
           "svm"={
             model <- svm(train_set, y, folds.no=5, sampling.no=10)
           })
    models[[drug_name]] <- model

  }
  if(!missing(result.path)){
    dev.off()
  }
  if(!missing(result.path)){
    save(predictions_range, file=sprintf("%s/prediction_ranges.RData", dir))
  }

  return(list("performance"=performance,  "model"=models))
}

