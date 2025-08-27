setwd("F:\\Green University\\Personal\\MSc in RUET\\MSc Thesis\\Thesis Work\\Journal\\Validation")
dataset <- read.csv("Zscores_and_stages_of_influential_genes.csv")
dataset <- dataset[2:11]

library(caTools)
set.seed(123)
split <- sample.split(dataset$Stage, SplitRatio = 0.7)
training <- subset(dataset, split == TRUE)
testing <- subset(dataset, split == FALSE)

training[-1] <- scale(training[-1])
testing[-1] <- scale(testing[-1])

#---------------------------------------------------------
# Apply SVM Algorithm
#---------------------------------------------------------
library(e1071)  

classifier_svm <- svm(formula = factor(Stage) ~ .,
                      data = training,
                      scale = TRUE,
                      type = 'C-classification',
                      kernel = 'radial',
                      cost = 1,
                      gamma = 0.001,
                      probability = TRUE)  # Enable probability estimation

y_pred_svm <- predict(classifier_svm, newdata = testing[-1], probability = TRUE)
svm_prob <- attr(y_pred_svm, "probabilities")[,2]  # Extract probability for the positive class

#---------------------------------------------------------
# Apply KNN Algorithm
#---------------------------------------------------------
library(class)

classifier_knn <- knn(train = training[-1],
                      test = testing[-1],
                      cl = training$Stage,
                      k = 7,
                      prob = TRUE)  # Enable probability estimation

knn_prob <- attr(classifier_knn, "prob")  # Extract probability of the winning class
#knn_prob <- ifelse(classifier_knn == levels(testing$Stage)[1], 1 - knn_prob, knn_prob)  # Adjust probability based on class

#---------------------------------------------------------
# Apply Naïve Bayes Algorithm
#---------------------------------------------------------
library(caret)

classifier_nb <- naiveBayes(Stage ~ ., data = training)
nb_prob <- predict(classifier_nb, newdata = testing, type = "raw")[,2]  # Extract probability for positive class

#---------------------------------------------------------
# Apply Random Forest Algorithm
#---------------------------------------------------------
library(randomForest)

classifier_rf <- randomForest(formula = factor(Stage) ~ .,
                              data = training,
                              scale = TRUE,
                              type = 'C-classification')

rf_prob <- predict(classifier_rf, newdata = testing, type = "prob")[,2]  # Extract probability for positive class

#---------------------------------------------------------
# ROC Curve for All Models
#---------------------------------------------------------
library(ROCR)

actual <- as.numeric(as.factor(testing$Stage)) - 1  # Convert Stage to binary (0,1)

# SVM
pred_svm <- prediction(svm_prob, actual)
perf_svm <- performance(pred_svm, "tpr", "fpr")

# KNN
pred_knn <- prediction(knn_prob, actual)
perf_knn <- performance(pred_knn, "tpr", "fpr")

# Naïve Bayes
pred_nb <- prediction(nb_prob, actual)
perf_nb <- performance(pred_nb, "tpr", "fpr")

# Random Forest
pred_rf <- prediction(rf_prob, actual)
perf_rf <- performance(pred_rf, "tpr", "fpr")

# Plot all ROC curves
plot(perf_svm, col = "blue", lwd = 2, main = "ROC Curves for Classification Models")
plot(perf_knn, col = "green", lwd = 2, add = TRUE)
plot(perf_nb, col = "red", lwd = 2, add = TRUE)
plot(perf_rf, col = "purple", lwd = 2, add = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray")  # Diagonal line

legend("bottomright",
       legend = c("SVM", "KNN", "Naïve Bayes", "Random Forest"),
       col = c("blue", "green", "red", "purple"),
       lwd = 2)