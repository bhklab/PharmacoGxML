featureSelection <- function(x, y, method, features.no, shrink=T, x_raw=NULL){
  if(is.null(x_raw)){
    x_raw <- x
  }
  switch(method,
         "mRMR"={
           if(shrink){
             var_features <- apply(x_raw, MARGIN=2, sd, na.rm=T)
             mad_features <- apply(x_raw, MARGIN=2, mad, na.rm=T)
             
             x <- x[, which(var_features > quantile(var_features, .75, na.rm=T) & 
                            mad_features > quantile(mad_features, .75, na.rm=T)), drop=FALSE]
           }
           
           f_data <- mRMR.data(data=as.data.frame(cbind(x, y), stringAsFactor=FALSE))
           features <- mRMR.ensemble(data=f_data,
                                     target_indices=ncol(x) + 1,
                                     feature_count=features.no,
                                     solution_count=1)
           features <- features@feature_names[unlist(features@filters)]
         },
         "variance"={
           var_features <- apply(x_raw, MARGIN=2, var)
           features <- names(sort(var_features, decreasing=T))[1:features.no]
         },
         "mad"={
           var_features <- apply(x_raw, MARGIN=2, mad)
           features <- names(sort(var_features, decreasing=T))[1:features.no]
         })
  return(features)
}
