# Load packages
library(caret)
library(glmnet)
library(foreach)
library(prospectr)
library(doParallel)
library(tidyverse)
library(ggrepel)
library(ggpubr)

# Set data and output directory
dataDir <- "E:/Thesis/EXTEND/"
outputDir <- dataDir

# Load data
load(paste0(dataDir, "Methylation/normDataMatrix.RData"))
load(paste0(dataDir, "Phenotype/CAIDE.RData"))
load(paste0(dataDir, "Phenotype/LIBRA.RData"))
load(paste0(dataDir, "Phenotype/metaData_ageFil.RData"))

#*****************************************************************************#
# Split samples with CAIDE1 score
#*****************************************************************************#

# Prepare data
all_X <- as.matrix(t(normDataMatrix))[CAIDE$ID,]
all_Y <- CAIDE
all(rownames(all_X) == all_Y$ID)

# Split data into training and test set
nTrain_male <- round(0.8*nrow(all_X[all_Y$Sex == "Male",]))
selectedSamples_male <- prospectr::kenStone(
  X = all_X[all_Y$Sex == "Male",], 
  k = nTrain_male,
  pc = 20,
  .center = TRUE,
  .scale = TRUE
)

nTrain_female <- round(0.8*nrow(all_X[all_Y$Sex == "Female",]))
selectedSamples_female <- prospectr::kenStone(
  X = all_X[all_Y$Sex == "Female",], 
  k = nTrain_female,
  pc = 20,
  .center = TRUE,
  .scale = TRUE
)


# Make training and test data
X_train <- t(all_X[c(selectedSamples_male$model, selectedSamples_female$model),])
Y_train <- all_Y[c(selectedSamples_male$model, selectedSamples_female$model),]

X_test <- t(all_X[c(selectedSamples_male$test, selectedSamples_female$test),])
Y_test <- all_Y[c(selectedSamples_male$test, selectedSamples_female$test),]


#*****************************************************************************#
# Check age distribution
#*****************************************************************************#

# Age
ggplot() +
  geom_histogram(data = Y_train, aes(x = Age, fill = "Train"), bins = 30, alpha = 0.5) +
  geom_histogram(data = Y_test, aes(x = Age, fill = "Test"), bins = 30, alpha = 0.5) +
  ylab("Count") +
  xlab("Age") +
  labs(fill = NULL) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom")

# Sex
table(Y_train$Sex)/sum(table(Y_train$Sex))
table(Y_test$Sex)/sum(table(Y_test$Sex))

# Cell type composition
Cell_train <- Y_train[,9:14]
rownames(Cell_train) <- Y_train$IID
Cell_train_gather <- gather(Cell_train)
Cell_train_gather$ID <- rep(rownames(Cell_train), ncol(Cell_train))

Cell_test <- Y_test[,9:14]
rownames(Cell_test) <- Y_test$IID
Cell_test_gather <- gather(Cell_test)
Cell_test_gather$ID <- rep(rownames(Cell_test), ncol(Cell_test))

Cell_total <- rbind.data.frame(Cell_train_gather, Cell_test_gather)
Cell_total$Data <- c(rep("Train", nrow(Cell_train_gather)),
                     rep("Test", nrow(Cell_test_gather)))


ggplot(Cell_total) +
  geom_point(aes(x = key, y = value, color = Data), 
             position = position_jitterdodge(), alpha = 0.2) +
  geom_violin(aes(x = key, y = value, fill = Data),
              alpha = 0.5, draw_quantiles = 0.5, color = "black", size = 0.5) +
  xlab("Cell Type") +
  ylab("Relative abundance") +
  labs(fill = NULL) +
  guides(color = "none") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1")


# CAIDE score
ggplot() +
  geom_histogram(data = Y_train, aes(x = SmokingScore, fill = "Train"), bins = 30, alpha = 0.5) +
  geom_histogram(data = Y_test, aes(x = SmokingScore, fill = "Test"), bins = 30, alpha = 0.5) +
  ylab("Count") +
  xlab("CAIDE1 Score") +
  labs(fill = NULL) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom")
