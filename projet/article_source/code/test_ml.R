library(e1071)
library(randomForest)
load(file = "act7009.rdata")
set.seed(1337)
validation_percentage = 0.15

assign_validation_split <- function(dataset, validation_percentage){
    nb_validation = floor(nrow(dataset)*validation_percentage)
    nb_train = nrow(dataset) - nb_validation
    validation_items = c(rep(1,nb_validation), rep(0, nb_train))
    validation_items = sample(validation_items, length(validation_items))
    dataset$sub_dataset = validation_items
        
    return(dataset)
}

mse <- function(actual, expected){
    mse = sqrt(sum((actual-expected)^2))

    return(mse)           
}

dataset_to_use = assign_validation_split(dataset, 0.15)

train_dataset = dataset_to_use[dataset_to_use$sub_dataset==0, ]
train_x = train_dataset[,1:(ncol(train_dataset)-4)]
train_y = train_dataset[,ncol(train_dataset)-1]
test_dataset = dataset_to_use[dataset_to_use$sub_dataset==1, ]
test_x = test_dataset[,1:(ncol(train_dataset)-4)]
test_y = test_dataset[,ncol(train_dataset)-1]

article_model = colSums(train_x*train_y)/nrow(train_dataset)
article_predicted = as.matrix(test_x) %*% article_model
article_error = article_predicted - test_y

svm_model = svm(x = train_x,y = train_y, kernel= "polynomial")
svm_predicted = predict(svm_model, test_x)
svm_error = svm_predicted-test_y

rf_model = randomForest(x = train_x,y = train_y, mtry = 10, nodesize = 10, ntree = 30)
rf_predicted = predict(rf_model, test_x)
rf_error = rf_predicted - test_y 

linear_regression_model_intersect_0 = lm(formula = train_y ~  . +0, data = cbind(train_x, train_y) )
linear_model_prediction_intersection_0 = predict(linear_regression_model_intersect_0, test_x)
lm_error_intersection_0 = linear_model_prediction_intersection_0 - test_y

linear_regression_model = lm(formula = train_y ~  ., data = cbind(train_x, train_y) )
linear_model_prediction = predict(linear_regression_model, test_x)
lm_error = linear_model_prediction - test_y

par(mfrow=c(1,5))

plot(article_error, ylim=c(-150,150), pch=20)
plot(svm_error, ylim=c(-150,150), pch=20)
plot(rf_error, ylim=c(-150,150), pch=20)
plot(lm_error_intersection_0, ylim=c(-150,150), pch=20)
plot(lm_error, ylim=c(-150,150), pch=20)

par(mfrow=c(1,1))

mse(test_y, article_predicted)
mse(test_y, svm_predicted)
mse(test_y, rf_predicted)
mse(test_y, linear_model_prediction)



