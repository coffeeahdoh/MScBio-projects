# Title: CodingExample_3
# Author: MichaelAngelo Confiado
# Last update: 2020-06-22
# Detail: Coding example from a course project to perform a quality check on a previous analyst's findings
#         of miRNA expression in different cancer types. Additionally, the analyses were extended to 
#         include other machine learning methods. Original data set is NDA protected, prompts from 
#         professor omitted, and only some analyses were kept."





########## previous analyst's code, with QC comments ##########
setwd("~/Documents/R projects/Project-2_AIinBio")

dat<-read.delim( file="miRNAs200TypesofCancerData.txt", 
                 sep= "\t",
                 stringsAsFactors = F) # calling the matrix and assigning to object "dat"
dat[1:10,1:10] # a glimpse of "dat"
dim(dat) # dimensions of dat

source("PCAforWideMatrices.R") # obtain f(x)s from PCAforWideMatrices.R

for(i in 3:ncol(dat)){
  dat[is.na(dat[,i]), i] <- mean(dat[,3], na.rm = TRUE)
} 

# forloop assigned means of only col 3 to ALL NAs of data set
# original analyst does not include column 2 (var2) 

test <- PCAforWideMatrices(dat[,3:ncol(dat)])

# perform PCA for dataset, saved to vector "test"
# analyst does not include column 2 in this as well

group <- dat$Type 

# create group vector to group data by type of cancer

library(ggplot2) 

GGPlotPCAWide(test, Groups = group, IndsNames = FALSE) #plot the results

### The previous analyst does not include column 2 throughout the analyses, which
### will produce inaccurate results. Also, the mean of col 3 was used
### to replace ALL NAs in the data set, affecting any further calculations. 
### In addition, the analyst did not scale the raw data which is
### an important step to PCA as the resulting plot is sensitive to the differing
### magnitudes of each variable.
### 
### Corrections are supplied below.





########## proposed updates ##########

# analyst did not represent all columns (3:ncol(dat))
# included col 2 in proceeding codes for PCA
for(i in 2:ncol(dat)){  
  dat[is.na(dat[,i]),i] <- mean(dat[,i], na.rm = TRUE)
}

# analyst did not scale data, PCA susceptible to non-uniform variance
# scaling data
dat_scaled <- scale(dat[,2:ncol(dat)], center = T, scale = T)
dat_scaled <- data.frame("Type" = dat[,1], dat_scaled)
summary(dat_scaled)

# PCA plot, inc. col 2
# renamed "test" to mirna_pca
mirna_pca <- PCAforWideMatrices(dat_scaled[,2:ncol(dat_scaled)])

# groups for PCA plot, renamed to mirna_pca_grps
mirna_pca_grps <- dat_scaled$Type

# new PCA plot, data scaled, inc. ALL cols 
GGPlotPCAWide(mirna_pca, Groups = mirna_pca_grps, IndsNames = FALSE)

### Compared to the previous analysts' plot, the new plot shows that the cancer
### types cluster more tightly. Extending out of the cluster center, there appears
### to be more mixture of components away from the center. In the previous analyst's
### plot, Mesenchymal cancer types appear to express associated RNAs differently than the other
### cancer types. In the new PCA analysis, projections suggest that all the cancer types are
### similar among the cancer types with few outliers. However, the PCA plot only 
### retains about 21% of the variance, indicating underlying factors may be influencing 
### miRNA expression profiles, and correlation based on only the two highest PCs
### does not fully encompass the data.





########## Naive Bayes performed by previous analyst ##########

library (e1071)

set.seed(456)
ind <- sample(2, nrow(dat),
              replace = TRUE,
              prob = c(0.55, 0.45))
training <- dat[ind==1,]
testing <- dat[ind==2,]
dim(training)
dim(testing)

training2<-data.frame(training)
training2[1:5,1:5]

str(training2)[1]
classf<-naiveBayes(x = training2[,-1], y = as.factor(training2[,1]),laplace = 1 )

testing2<-data.frame(testing)
testing2[1:5,1:5]

NBmodelpredict<-predict(classf, testing2[,-1], type = "class")
NBmodelpredict
NBtab<-table(testing[,1],NBmodelpredict);NBtab
sum(diag(NBtab))/sum(NBtab)





########## Attempt to improve previous NB, by transforming to categorical data ##########

# Categorizing training data set, where value < mean = low, value > mean = high
# Naive Bayes may better classify with categorized data

#categorized training set
train_vals <- training[,-1]
train_cats <- NULL
for(i in 1:ncol(train_vals)){ 
  train_cats[i] <- list(ifelse(train_vals[i] < mean(train_vals[,i]),"low","high"))
  train_cats <- as.data.frame(train_cats)
}
colnames(train_cats) <- colnames(training[,-1])
train_cats <- cbind(Type = training[,1], train_cats); dim(train_cats); dim(training)

#categorized testing set
test_vals <- testing[,-1]
test_cats<- NULL
for(i in 1:ncol(test_vals)){ 
  test_cats[i]<- list(ifelse(test_vals[i] < mean(test_vals[,i]),"low","high"))
  test_cats <- as.data.frame(test_cats)
}
colnames(test_cats) <- colnames(testing[,-1])
test_cats <- cbind(Type = testing[,1], test_cats); dim(test_cats); dim(testing)

# Naive Bayes classifier model, w/ laplace = 1
NB_model_cat <- naiveBayes(x = train_cats[,-1], y = as.factor(train_cats[,1]), laplace = 1 )

# Predictions and accuracy using testing object as new input data
NB_model_cat_pred <- predict(NB_model_cat, test_cats[,-1], type = "class")
NB_model_cat_pred
NB_model_cat_tab <- table(test_cats[,1], NB_model_cat_pred); NB_model_cat_tab
sum(diag(NB_model_cat_tab))/sum(NB_model_cat_tab)

# Kappa comparisons
#install.packages("vcd")
library("vcd")
Kappa(NBtab)
Kappa(NB_model_cat_tab)

### The previous analyst's code for Naive Bayes classifier correctly accounts for all data,
### however, it has poor agreement with their model's prediction compared to actual 
### observations (Kappa = 0.1976, unweighted; Kappa = 0.1862, weighted; accuracy = 39.07%).
### Using their model would not be any more efficient than diagnosing cancer based on chance.
### An improvement of the model is done by categorizing the data set. In this case,
### assigning "low" to values below the mean per variable and "high" to values above the
### mean per variable. This results in Kappa = .3909, unweighted; Kappa = 0.4790, weighted 
### ; and an accuracy of ~53% which represents a fair to moderate model for predicting
### cancer type based on Naive Bayes. Using Naive Bayes to classify cancer types
### based on miRNA data may not be efficient as categorizing data eliminates minute
### differences.





########## Neural Network performed by previous analyst ##########

# load neuralnet
library("neuralnet")

# min-max normalization
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
dat_norm<-data.frame(sapply(dat[,-1], normalize)); summary(dat_norm)
dat_norm<-cbind("Type" = dat[,1], dat_norm)

# Creation of training and testing sets based on normalized data
train_nn <- dat_norm[ind==1,]  
test_nn <- dat_norm[ind==2,]  

# ensuring training and testing data sets are the same
# as those of the previous analyst's 
dim(train_nn); rownames(training) %in% rownames(train_nn)
dim(test_nn); rownames(testing) %in% rownames(test_nn)

# Artifical Neural Network model, 5 neurons in 1 hidden layer
# -- per suggestion of previous analyst
set.seed(123)
nn_fivone_model <- neuralnet(as.character(Type)~.,
                data = train_nn,
                hidden = 5,
                err.fct = "ce",
                linear.output = F)
plot(nn_fivone_model)

# Prediction -- training model
nn_fivone_pred <- predict(nn_fivone_model, train_nn[,-1])
head(nn_fivone_pred)
nn_fivone_pred <- apply(nn_fivone_pred, MARGIN = 1, which.max); nn_fivone_pred
nn_fivone_pred <- ifelse(nn_fivone_pred == 1, "Classical", 
                        ifelse(nn_fivone_pred == 2, "Mesenchymal", 
                              ifelse(nn_fivone_pred == 3, "Neural", 
                                    "Proneural")))
nn_fivone_tab <- table(nn_fivone_pred, train_nn$Type); nn_fivone_tab

# Accuracy & Error Rate -- training model
sum(diag(nn_fivone_tab))/sum(nn_fivone_tab)
1 - sum(diag(nn_fivone_tab ))/sum(nn_fivone_tab)

# Prediction -- testing
nn_fivone_pred_test <- predict(nn_fivone_model, test_nn[,-1])
nn_fivone_pred_test <- apply(nn_fivone_pred_test, MARGIN = 1, which.max)
nn_fivone_pred_test <- ifelse(nn_fivone_pred_test == 1, "Classical", 
                 ifelse(nn_fivone_pred_test == 2, "Mesenchymal", 
                        ifelse(nn_fivone_pred_test == 3, "Neural", 
                               "Proneural")))
nn_fivone_tab_test<- table(nn_fivone_pred_test, test_nn$Type)
nn_fivone_tab_test

# Accuracy & Error Rate 
sum(diag(nn_fivone_tab_test))/sum(nn_fivone_tab_test)
1 - sum(diag(nn_fivone_tab_test))/sum(nn_fivone_tab_test)

# Kappa
Kappa(nn_fivone_tab_test)

### Compared to NB + categorization, an ANN w/ 5 neurons in one hidden layer
### was less accurate less agreeable with expected (actual) data
### (Kappa = 0.3320, unweighted; Kappa = 0.4248, weighted; Accuracy = 0.507). 
### Although it has fair agreement. This particular model will not be more
### efficient than the previous model. 





########## Attempt to improve previous Neural Network ##########

# Artifical Neural Network model, 10 neurons/hidden layer, hidden layers =  2
# --New analysis
set.seed(123)
nn_tentwo_model <- neuralnet(as.character(Type)~.,
                data = train_nn,
                hidden = c(10,10),
                err.fct = "ce",
                linear.output = F)
plot(nn_tentwo_model)

# Prediction -- training model
nn_tentwo_pred <- predict(nn_tentwo_model, train_nn[,-1])
head(nn_tentwo_pred)
nn_tentwo_pred <- apply(nn_tentwo_pred, MARGIN = 1, which.max); nn_tentwo_pred
nn_tentwo_pred <- ifelse(nn_tentwo_pred == 1, "Classical", 
               ifelse(nn_tentwo_pred == 2, "Mesenchymal", 
                      ifelse(nn_tentwo_pred == 3, "Neural", 
                             "Proneural")))
nn_tentwo_tab <- table(nn_tentwo_pred, train_nn$Type); nn_tentwo_tab 

# Accuracy & Error Rate -- training model
sum(diag(nn_tentwo_tab))/sum(nn_tentwo_tab)
1 - sum(diag(nn_tentwo_tab))/sum(nn_tentwo_tab)

# Prediction -- testing
nn_tentwo_pred_test <- predict(nn_tentwo_model, test_nn[,-1])
nn_tentwo_pred_test <- apply(nn_tentwo_pred_test, MARGIN = 1, which.max)
nn_tentwo_pred_test <- ifelse(nn_tentwo_pred_test == 1, "Classical", 
                 ifelse(nn_tentwo_pred_test== 2, "Mesenchymal", 
                        ifelse(nn_tentwo_pred_test == 3, "Neural", 
                               "Proneural")))
nn_tentwo_tab_test <- table(nn_tentwo_pred_test, test_nn$Type)
nn_tentwo_tab_test

# Accuracy & Error Rate -- testing
sum(diag(nn_tentwo_tab_test))/sum(nn_tentwo_tab_test)
1 - sum(diag(nn_tentwo_tab_test))/sum(nn_tentwo_tab_test)

# Kappa
Kappa(tab8a)

### An ANN with 10 neurons per hidden layer with 2 hidden layers is more efficient
### in classifying cancer types based on miRNA expression than NB + categorization
### (Kappa = 0.4090, unweighted; Kappa = 0.4536, weighted; Accuracy =  0.5674419).
### This ANN model has moderate agreement with actual results of the testing 
### set data. 





########## Decision Tree classification attempt ##########

library("C50")
library("caret")

# Creation of control object
ctrl <- caret::trainControl(method = "cv", number = 10,
                            selectionFunction = "best")

# Decision Tree performance evaluation
set.seed(123)
dt_model <- caret::train(as.factor(Type) ~ ., 
                      data = training, method = "C5.0",
                      metric = "Kappa",
                      trControl = ctrl)
dt_model
plot(dt_model)

### Best parameters --  model = tree  winnow = FALSE   
###                     trials = 20   
###                     accuracy = 0.6244850 
###                     kappa = 0.4839259

# Predictions -- testing, not used to train model
dt_pred<- predict(dt_model, testing[,-1])
dt_tab<- table(tree_pred,testing$Type); dt_tab

# Accuracy & Error Rate -- testing, not used to train model
sum(diag(dt_tab))/sum(dt_tab) # 0.5767442
1-sum(diag(dt_tab))/sum(dt_tab)

# Kappa
Kappa(dt_tab)

plot(dt_model)





########## Decision Tree performance, ROC curves and AUC per cancer type ##########

library(ROCR)

# Accessing levels of Types
lvls <- levels(as.factor(dat$Type))

# Base plot for ROC curve
plot(x = NA, y= NA,
     xlim=c(0,1), ylim=c(0,1),
     bty = "n", 
     main = "ROC Curves for Cancer Types",
     ylab = "True positive rate",
     xlab = "False positive rate")

# for loop -- creating binary classifications for each cancer type (one v. all) 
for (type.id in 1:4) {
  facts = as.factor(training$Type == lvls[type.id])
  
  set.seed(123)
  dt_model_roc <- C5.0(x = training[, -1], y = facts, iterations= 20)
  dt_pred_roc <- predict(dt_model_roc, testing[,-1], type='prob')
  
  score <- dt_pred_roc[, 'TRUE']
  actual_class <- testing$Type == lvls[type.id]
  
  score_pred <- prediction(score, actual_class)
  c50perf <- performance(score_pred, "tpr", "fpr")
  
  roc.x = unlist(c50perf@x.values)
  roc.y = unlist(c50perf@y.values)
  lines(roc.y ~ roc.x, col=type.id+1, lwd=2)
  
  c50auc <- performance(score_pred, "auc")
  
  c50auc <- unlist(slot(c50auc, "y.values"))
  aucs[type.id] <- c50auc
}
abline(a = 0, b = 1, lwd = 2, lty = 3)
legend(0.7, 0.5, legend=lvls,
       col=c(2, 3, 4, 5), lty=1, cex=0.8, lwd = 2)

AUC.Classical <- aucs[1]; AUC.Classical 
AUC.Mesenchymal <- aucs[2]; AUC.Mesenchymal 
AUC.Neural <- aucs[3]; AUC.Neural 
AUC.Proneural <- aucs[4]; AUC.Proneural 

### An initial impression of the curves suggest that the decision tree model 
### cannot discrimination neural cancer types from other types. But has fair performance
### in discriminating Classical, Mesenchymal, and Proneural types from others.





########## Attempt to improve decision tree performance with Bagging ##########

# Accessing Bagging function from ipred
library(ipred)

# Bagging model
set.seed(123)
ctrl_bag <- trainControl(method = "cv", number = 10)
dt_bag_model <- train(as.factor(Type) ~ ., data = training, method = "treebag", trControl = ctrl)

# Predictions -- testing
dt_bag_pred <- predict(dt_bag_model,testing)
dt_bag_tab <- table (dt_bag_pred, testing$Type); dt_bag_tab

# Accuracy of the bagging model
sum(diag(dt_bag_tab))/sum(dt_bag_tab) 

# Kappa
Kappa(dt_bag_tab)

### Arguments used in train() were method = "treebag", method name for bagging in the 
### ipred package, and trControl = ctrl to call the 10-fold cross-validation resampling.
### This results in Kappa values 0.4195 (unweighted) and 0.4540 (weighted) with an 
### accuracy of 0.5767442. These results are similar to the decision tree analysis
### which is reasonable since bagging is an aggregating process that results in 
### the mean of all sub-models, and specifically for this bagging method, it is 
### based on the decision tree. 





########## Attempt to improve decision tree performance with Boosting ##########

# Accessing adaboost from adabag
library ("adabag")

# data prep -- previous analyst's training and testing data sets 
training[2:ncol(training)] <- sapply(training[2:ncol(training)],as.numeric)
training$Type <- as.factor(training$Type)
str(training)
testing[2:ncol(testing)] <- sapply(testing[2:ncol(testing)],as.numeric)
testing$Type <- as.factor(testing$Type)


# AdaBoost model
set.seed(123)
dt_boost_model <- caret::train(as.factor(Type) ~ ., 
                        data = training, 
                        method = "AdaBag",
                        metric = "Kappa", 
                        trControl = ctrl_bag)

# Predictions & Accuracy
dt_boost_pred<- predict(dt_boost_model, testing[,-1]); dt_boost_pred
dt_boost_tab<- table(dt_boost_pred, testing$Type); dt_boost_pred
sum(diag(dt_boost_tab))/sum(dt_boost_tab)

# Kappa
Kappa(dt_boost_tab)

### AdaBoost results in Kappas 0.4565 (unweighted) and 0.4738 (weighted) with an 
### accuracy of ~60.5%. This boosting model appears better than both Bagging and 
### the original Decision Tree models. Decision trees are applied in series for this 
### model in addition to strengthening (boosting) weak classifiers per cycle. This improves
### the performance based on the previous decision tree as each prior tree contains
### a "boosted" weak classifier.





########## Random Forest classification attempt ##########

# install.packages("randomForest")
library("randomForest")

# Random Forest performance, using same training and testing sets as the ones
# created earlier to maintain consistency.
set.seed(123)
rf_model <- caret::train(as.factor(Type) ~ ., 
                        data = training, 
                        method = "rf",
                        metric = "Kappa", 
                        trControl = ctrl)
rf_model

# Predict, RF
rf_pred <- predict(rf_model, testing)
rf_tab<- table(rf_pred, testing$Type); rf_tab

# Accuracy and Error Rate, RF
sum(diag(rf_tab))/sum(rf_tab)
1-sum(diag(rf_tab))/sum(rf_tab)

# Kappa
Kappa(rf_tab)

### The random forest model provided the best outcome for accuracy and kappa values.
### Accuracy is about 63.3% and kappa values are 0.4951 (unweighted) and 0.5311
### (weighted). Even though kappa values are only moderately agreeable with expected
### outcomes, this model offers the most efficient model out of the models trialed in
### terms of classifying cancer types of this data set.





########## Attempt to improve Random Forest ##########

ctrl_cv20 <- trainControl(method = "cv", number = 20)

grid <- expand.grid(.mtry = c(101, 110, 115, 120, 125, 130, 150, 200))

set.seed(123)
rf_model_imp <- caret::train(as.factor(Type) ~ ., 
                        data = training, 
                        method = "rf",
                        metric = "Kappa", 
                        trControl = ctrl_cv20,
                        tuneGrid = grid)
rf_model_imp

rf_pred_imp <- predict(rf_model_imp, testing)
rf_tab_imp <- table(rf_pred_imp, testing$Type); rf_tab_imp

sum(diag(rf_tab_imp))/sum(rf_tab_imp)
1-sum(diag(rf_tab_imp))/sum(rf_tab_imp)





########## Comparison of Models ##########

### 
## Model                                                    Accuracy w/ new data (testing)  Kappa(unweighted) Kappa(weighted)
### Naive Bayes (previous analyst)                            46.0%                          0.2937            0.2651
### Naive Bayes (new)                                         54.0%                          0.3909            0.4790
### Artificial Neural Network, 5 neurons, 1 hidden layer      50.7%                          0.3320            0.4248
### Artificial Neural Network, 10 neurons / layer, 2 layers   56.7%                          0.4090            0.4536
### Decision Tree, 20 iterations                              57.7%                          0.4195            0.4928
### Bagging                                                   57.7%                          0.4195            0.4540
### AdaBoost                                                  60.5%                          0.4565            0.4738
### Random Forest                                             63.3%                          0.4951            0.5311
###
###
### Among the models tested Random Forest offered the most efficient classification
### model. While all models were efficient in training their respective algorithms,
### The resulting testing sets were the main indicators of performance.
### With 63.3% accuracy in classifying new data (testing set) and moderate agreement
### with the known types of cancer for the testing data set random forest outperformed 
### prior models. 