make_predictions <- function(model, newdata){
	UseMethod("make_predictions", model)
}

make_predictions.cv.glmnet <- function(model, newdata){
	return(predict(model, newx=newdata, s="lambda.min"))
}

make_predictions.train <- function(model, newdata){
	return(predict(model, newdata=newdata, s="lambda.min"))
}

