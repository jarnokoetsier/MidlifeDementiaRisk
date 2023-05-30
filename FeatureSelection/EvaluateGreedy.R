# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(caret)
library(pROC)
library(randomForest)

# Load data
load("~/selectedFeatures.RData") # Selected features by Greedy approach
load("~/Data/X_CAIDE1.RData")
load("~/Data/Y_CAIDE1.RData")
load("~/Data/X_test.RData")
load("~/Data/Y_test.RData")

# Training data:
Y_train <- Y_CAIDE1

# Convert to M values
X_train_M <- log2(X_CAIDE1[selectedFeatures,]/(1-X_CAIDE1[selectedFeatures,]))
X_test_M <- log2(X_test[selectedFeatures,]/(1-X_test[selectedFeatures,]))

# Scale data
X_train_s <- t((X_train_M -rowMeans(X_train_M))/apply(X_train_M,1,sd))
X_test_s <- t((X_test_M -rowMeans(X_train_M))/apply(X_train_M,1,sd))

# Linear model + RFE
perf_EN <- rep(NA, length(selectedFeatures))

# Recursive Feature elimination
for (f in 1:length(selectedFeatures)){
  # Start with all features
  features_temp <- selectedFeatures[1:f]
  
  # Make formula
  Y <- Y_train$CAIDE
  X_fit <- cbind.data.frame(Y, as.data.frame(X_train_s[,features_temp]))
  colnames(X_fit) <- c("Y", features_temp)
  
  # Fit logistic regression model
  model <- lm(Y ~ ., data = X_fit)
  
  # Predict
  predX <- as.data.frame(X_test_s[,features_temp])
  colnames(predX) <- features_temp
  pred <- predict(model, predX)
  
  # Get AUC
  perf_EN[f] <- RMSE(pred = pred, obs = Y_test$CAIDE)
}


# Random Forest + RFE
perf_RF <- rep(NA, length(selectedFeatures))

# Recursive Feature elimination
for (f in 1:length(selectedFeatures)){
  # Start with all features
  features_temp <- selectedFeatures[1:f]
  
  # Make formula
  Y <- Y_train$CAIDE
  X_fit <- cbind.data.frame(Y, as.data.frame(t(X_train_M)[,features_temp]))
  colnames(X_fit) <- c("Y", features_temp)
  
  # Fit logistic regression model
  set.seed(123)
  model <- randomForest(Y ~ ., data = X_fit)
  
  # Predict
  predX <- as.data.frame(t(X_test_M)[,features_temp])
  colnames(predX) <- features_temp
  pred <- predict(model, predX)
  
  # Get AUC
  perf_RF[f] <- RMSE(pred = pred, obs = Y_test$CAIDE)
}

# Prepare data for plotting
plotData <- data.frame(nFeatures = c(1:length(perf_EN),
                                     1:length(perf_RF)),
                       RMSE = c(perf_EN, perf_RF),
                       Method = c(rep("Linear Regression", length(perf_EN)),
                                  rep("Random Forest", length(perf_RF))))

# Construct plot
p <- ggplot(plotData) +
  geom_point(aes(x = nFeatures, y = RMSE, color = Method)) +
  geom_line(aes(x = nFeatures, y = RMSE, color = Method), alpha = 0.5) +
  geom_hline(yintercept = 1.916394, linetype = "dashed", color = "#A50F15")+
  geom_hline(yintercept = 1.966543, linetype = "dashed", color = "#CC4C02")+
  geom_hline(yintercept = 1.902888, linetype = "dashed", color = "#6A51A3")+
  scale_color_brewer(palette = "Dark2") +
  xlab("# Features") +
  ggtitle("CAIDE1", subtitle = "Performance in Test Set") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.subtitle = element_text(hjust = 0.5,
                                     size = 10,
                                     face = "italic"))


# Save plot
ggsave(p, file = "GreedyResults.png", width = 8, height = 6)



