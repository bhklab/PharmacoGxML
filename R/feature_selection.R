featureSelection <- function(x, y, method, features.no){
  switch(method,
         "mRMR"={
           #var_features <- apply(x, MARGIN=2, var)
           #features <- names(sort(var_features, decreasing=T))[1:(features.no*10)]
           f_data <- mRMR.data(data=as.data.frame(cbind(x[,features], y), stringAsFactor=FALSE))
           system.time(
           features <- mRMR.ensemble(data=f_data,
                                     target_indices=ncol(x) + 1,
                                     feature_count=features.no,
                                     solution_count=1)
           )
           features <- features@feature_names[unlist(features@filters)]
         },
         "variance"={
           var_features <- apply(x, MARGIN=2, var)
           features <- names(sort(var_features, decreasing=T))[1:features.no]
         })
  return(features)
}
