# Clear workspace and console
rm(list = ls())
cat("\014") 

# Install packages
install.packages("prospectr")
install.packages("ggridges")

# Load packages
library(prospectr)
library(tidyverse)

# Load data
load("CAIDE.Rdata")
load("EPILIBRA.Rdata")
load("metaData_ageFil.RData")
load("methSet_allNorm_fil.RData")
load("cellType.RData")
load("gt_results.RData")

#*****************************************************************************#
# 1.1. Apply Kennard-Stone algorithm to split the data
#*****************************************************************************#

# 1. Get samples with all scores available

# All samples (age filtered)
all_samples <- unique(dat$ID)

# CAIDE1 samples
CAIDE1_samples <- unique(CAIDE$ID)
all(CAIDE1_samples %in% all_samples)

# CAIDE2 samples
CAIDE2 <- CAIDE[CAIDE$ID %in% unique(names(gt_rs7412)),]
CAIDE2_samples <- unique(names(gt_rs7412)[names(gt_rs7412) %in% CAIDE1_samples])
all(CAIDE2_samples %in% all_samples)

# LIBRA samples
LIBRA_samples <- EPILIBRA$ID
all(LIBRA_samples %in% all_samples)

# Get intersection
intersect_samples <- intersect(intersect(CAIDE1_samples, CAIDE2_samples),LIBRA_samples)

# Get Basename instead of ID
intersect_samples <- dat$Basename[dat$ID %in% intersect_samples]


# Prepare data
all_X <- as.matrix(t(methSet_allNorm_fil))[intersect_samples,]
all_Y <- dat[dat$Basename %in% intersect_samples, 
             c("ID", "Basename", "Position", "Plate", "Age", "Sex")]
all_Y <- inner_join(all_Y, cellType, by = c("Basename" = "ID"))
all(rownames(all_X) == all_Y$Basename)

# Split data into training and test set
nTrain_male <- round(0.8*nrow(all_X[all_Y$Sex == 1,]))
selectedSamples_male <- prospectr::kenStone(
  X = all_X[all_Y$Sex == 1,], 
  k = nTrain_male,
  pc = 20,
  .center = TRUE,
  .scale = TRUE
)

nTrain_female <- round(0.8*nrow(all_X[all_Y$Sex == 2,]))
selectedSamples_female <- prospectr::kenStone(
  X = all_X[all_Y$Sex == 2,], 
  k = nTrain_female,
  pc = 20,
  .center = TRUE,
  .scale = TRUE
)

female_Y <- all_Y[all_Y$Sex == 2,]
male_Y <- all_Y[all_Y$Sex == 1,]
TestTrain <- list(Train = c(male_Y$Basename[selectedSamples_male$model], female_Y$Basename[selectedSamples_female$model]),
                  Test = c(male_Y$Basename[selectedSamples_male$test], female_Y$Basename[selectedSamples_female$test]))

# Save outcome
save(TestTrain, file = "TestTrain.RData")

# Make training and test data
load("TestTrain.RData")


#*****************************************************************************#
# 1.2. Make all data combinations
#*****************************************************************************#

# All
all_X <- as.matrix(t(methSet_allNorm_fil))
all_Y <- dat[,c("ID", "Basename", "Position", "Plate", "Age", "Sex")]
all_Y <- inner_join(all_Y, cellType, by = c("Basename" = "ID"))
all(rownames(all_X) == all_Y$Basename)

# CAIDE1
length(intersect(CAIDE$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE1 <- t(all_X[setdiff(CAIDE$Basename, TestTrain$Test),])
Y_CAIDE1 <- all_Y[all_Y$Basename %in% setdiff(CAIDE$Basename, TestTrain$Test),]
all(colnames(X_CAIDE1) == Y_CAIDE1$Basename)

# CAIDE2
length(intersect(CAIDE2$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE2 <- t(all_X[setdiff(CAIDE2$Basename, TestTrain$Test),])
Y_CAIDE2 <- all_Y[all_Y$Basename %in% setdiff(CAIDE2$Basename, TestTrain$Test),]
all(colnames(X_CAIDE2) == Y_CAIDE2$Basename)

# LIBRA
length(intersect(EPILIBRA$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_LIBRA <- t(all_X[setdiff(EPILIBRA$Basename, TestTrain$Test),])
Y_LIBRA <- all_Y[all_Y$Basename %in% setdiff(EPILIBRA$Basename, TestTrain$Test),]
all(colnames(X_LIBRA) == Y_LIBRA$Basename)

# Non-test
length(intersect(all_Y$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_nonTest <- t(all_X[setdiff(all_Y$Basename, TestTrain$Test),])
Y_nonTest <- all_Y[all_Y$Basename %in% setdiff(all_Y$Basename, TestTrain$Test),]
all(colnames(X_nonTest) == Y_nonTest$Basename)

# Test
X_test <- t(all_X[TestTrain$Test,])
Y_test <- all_Y[all_Y$Basename %in% TestTrain$Test,]



#*****************************************************************************#
# 1.2. Check distribution of age, sex, cell type, and CAIDE1 score
#*****************************************************************************#

#=============================================================================#
# Which samples belong where
#=============================================================================#

# Make data frame for plotting
plotOverlap <- data.frame(ID = c(unique(all_Y$ID),
                                 unique(c(Y_CAIDE1$ID, Y_test$ID)),
                                 unique(c(Y_CAIDE2$ID,Y_test$ID)),
                                 unique(c(Y_LIBRA$ID,Y_test$ID))),
                          Source = c(rep("All",length(unique(all_Y$ID))),
                                     rep("CAIDE1", length(c(Y_CAIDE1$ID, Y_test$ID))),
                                     rep("CAIDE2",length(unique(c(Y_CAIDE2$ID,Y_test$ID)))),
                                     rep("LIBRA",length(unique(c(Y_LIBRA$ID,Y_test$ID)))))
) 
plotOverlap$Test <- ifelse(plotOverlap$ID %in% Y_test$ID, "Test", "Training")

# Order of samples in the plot
orderID <- names(rev(sort(table(plotOverlap$ID))))
plotOverlap$ID <- factor(plotOverlap$ID, levels = orderID)
plotOverlap$Source <- factor(plotOverlap$Source,
                             c("All", "CAIDE1", "CAIDE2","LIBRA"))

# Color in plot
plotOverlap$Color <- paste(plotOverlap$Test, plotOverlap$Source, sep = " ")
plotOverlap$Color[plotOverlap$Test == "Test"] <- "Test"

# Make plot
overlap_plot <- ggplot() +
  geom_tile(data = plotOverlap, aes(x = ID, y = Source, fill = Color)) +
  facet_grid(cols = vars(Test), scales = "free", space = "free") +
  xlab("Samples") + 
  ylab(NULL) +
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = "none",
        axis.text.x = element_blank())

# Save plot
ggsave(overlap_plot, file = "OverlapPlot.png", width = 8, height = 6)

#=============================================================================#
# Age
#=============================================================================#
plotData <- rbind.data.frame(Y_test, Y_nonTest, Y_CAIDE1, Y_CAIDE2, Y_LIBRA)
plotData$Set <- c(rep("Test", nrow(Y_test)),
                  rep("Training All", nrow(Y_nonTest)),
                  rep("Training CAIDE1", nrow(Y_CAIDE1)),
                  rep("Training CAIDE2", nrow(Y_CAIDE2)),
                  rep("Training LIBRA", nrow(Y_LIBRA)))

# Make plot
plot_age <- ggplot(plotData, aes(x = Age, y = Set, fill = Set)) + 
  ggridges::geom_density_ridges(stat = "binline", bins = 30, 
                                scale = 0.95) +
  ylab(NULL) +
  xlab("Age") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.position = "none")

# Save plot
ggsave(plot_age, file = "Age_TrainingAndTest.png", width = 8, height = 6)

# Two-sample Kolmogorov-Smirnov test
ks.test(Y_test$Age, Y_nonTest$Age)
ks.test(Y_test$Age, Y_CAIDE1$Age)
ks.test(Y_test$Age, Y_CAIDE2$Age)
ks.test(Y_test$Age, Y_LIBRA$Age)

#=============================================================================#
# Sex
#=============================================================================#

# Format data
perc <- table(paste0(plotData$Sex, plotData$Set))/rep(table(plotData$Set),2)
percData <- data.frame(Percentage = perc,
                       Sex = c(rep("Male",5), rep("Female",5)),
                       Set = rep(unique(plotData$Set),2))

# Make plot
plot_sex <- ggplot(percData, aes(x = Set, fill = Sex)) +
  geom_bar(aes(y = Percentage.Freq), color = "black", stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  ylab("Sex Proportion") +
  xlab(NULL) +
  theme_classic()

# Save plot
ggsave(plot_sex, file = "Sex_TrainingAndTest.png", width = 8, height = 6)


#=============================================================================#
# Cell type composition
#=============================================================================#

plotCell <- plotData[,7:12]
colnames(plotCell) <- c("CD8 T-cells", "CD4 T-cells", "NK cells", "B-cells", 
                        "Monocytes","Neutrophils")
plotCell <- gather(plotCell)
plotCell$Set <- rep(plotData$Set,6)
plotCell$ID <- rep(plotData$ID,6)
plotCell$Test <- rep(plotData$Test,6)

plot_cell <- ggplot(plotCell, aes(x = key, y = value, fill = Set, color = Set)) +
  geom_point(position = position_jitterdodge(), alpha = 0.2) +
  geom_violin(alpha = 0.5, draw_quantiles = 0.5, color = "black", size = 0.5) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  ylab("Cell Type Proportion") +
  xlab(NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(plot_cell, file = "Cell_TrainingAndTest.png", width = 10, height = 6)



#=============================================================================#
# CAIDE1 score
#=============================================================================#

# Make plot
plot_caide <- ggplot() +
  geom_histogram(data = CAIDE[CAIDE$ID %in% Y_CAIDE1$ID,], aes(x = CAIDE, fill = "Train"), bins = 10, alpha = 0.5) +
  geom_histogram(data = CAIDE[CAIDE$ID %in% Y_test$ID,], aes(x = CAIDE, fill = "Test"), bins = 10, alpha = 0.5) +
  ylab("Count") +
  xlab("CAIDE1 Score") +
  labs(fill = NULL) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom")

ks.test(CAIDE$CAIDE[CAIDE$ID %in% Y_CAIDE1$ID], CAIDE$CAIDE[CAIDE$ID %in% Y_test$ID])

# Save plot
ggsave(plot_caide, file = "CAIDE1_TrainingAndTest.png", width = 8, height = 6)

#=============================================================================#
# LIBRA score
#=============================================================================#

# Make plot
plot_libra <- ggplot() +
  geom_histogram(data = EPILIBRA[EPILIBRA$ID %in% Y_LIBRA$ID,], aes(x = LIBRA, fill = "Train"), bins = 10, alpha = 0.5) +
  geom_histogram(data = EPILIBRA[EPILIBRA$ID %in% Y_test$ID,], aes(x = LIBRA, fill = "Test"), bins = 10, alpha = 0.5) +
  ylab("Count") +
  xlab("LIBRA Score") +
  labs(fill = NULL) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom")

ks.test(EPILIBRA$LIBRA[EPILIBRA$ID %in% Y_LIBRA$ID], EPILIBRA$LIBRA[EPILIBRA$ID %in% Y_test$ID])

# Save plot
ggsave(plot_libra, file = "LIBRA_TrainingAndTest.png", width = 8, height = 6)


#=============================================================================#
# PCA
#=============================================================================#

# Unit scale the training Data
trainingData_scaled <- t((X_nonTest - rowMeans(X_nonTest))/(apply(X_nonTest,1,sd)))

# Make PCA model
pcaList_train <-  prcomp(trainingData_scaled,        
                         retx = TRUE,
                         center = TRUE,
                         scale = TRUE,
                         rank. = 6)

# Get the PCA scores of the training data
scores_train <- as.data.frame(pcaList_train$x)
scores_train$ID <- rownames(scores_train)

# Calculate the explained variance of the PCs
explVar <- round(((pcaList_train$sdev^2)/sum(pcaList_train$sdev^2))*100,2)

# Scale test data (using the standard deviation and mean of the training data)
testData_scaled <- t((X_test - rowMeans(X_nonTest))/(apply(X_nonTest,1,sd)))

# Calculate the scores of the test data
scores_test <- as.data.frame(as.matrix(testData_scaled) %*% as.matrix(pcaList_train$rotation))
scores_test$ID <- rownames(scores_test)

# Combine the scores of test and training data in a single data frame
scores_all <- rbind.data.frame(scores_train, scores_test)
scores_all$Train <- c(rep("Training All", nrow(scores_train)), rep("Test", nrow(scores_test)))

# Plot the scores of the training data and the project scores of the test data 
PCA_TrainTest_12 <- ggplot() +
  geom_point(data = scores_all, aes(x = PC1, y = PC2, shape = Train, color = Train), 
             size = 2, alpha = 0.8) +
  scale_shape_manual(values = c(15,17,0,2)) +
  xlab(paste0("PC1 (", explVar[1], "%)")) +
  ylab(paste0("PC2 (", explVar[2], "%)")) +
  labs(title = NULL) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16),
        plot.caption = element_text(hjust = 0.5,
                                    size = 10,
                                    face = "italic")) +
  scale_color_brewer(palette = "Set1")

# Save plot
ggsave(PCA_TrainTest_12, file = "PC1vs2_TrainingAndTest.png", width = 8, height = 6)


# Residuals
# Split data into reconstruction and residuals based on the first ten PC's
loadings <- pcaList_train$rotation

reconstruction_train <- as.matrix(scores_train[,1:2]) %*% t(loadings[,1:2])
residuals_train <- trainingData_scaled - reconstruction_train

reconstruction_test <- as.matrix(scores_test[,1:2]) %*% t(loadings[,1:2])
residuals_test <- testData_scaled - reconstruction_test

# Calculate the orthogonal distances
ortDist_train <- sqrt(rowSums(residuals_train^2))
ortDist_test <- sqrt(rowSums(residuals_test^2))

plotOrt <- data.frame(Distance = c(ortDist_train,ortDist_test),
                      Set = c(rep("Training All", length(ortDist_train)),
                              rep("Test", length(ortDist_test))))

# Make histogram of orthogonal distances
resPlot <- ggplot() +
  geom_histogram(data = plotOrt[plotOrt$Set == "Training All",], 
                 aes(x = Distance, fill = Set), alpha = 0.5, color = "white") +
  geom_histogram(data = plotOrt[plotOrt$Set == "Training All",], 
                 aes(x = Distance, fill = Set), alpha = 0.5) +
  geom_histogram(data = plotOrt[plotOrt$Set == "Test",],
                 aes(x = Distance, fill = Set), alpha = 0.5) +
  xlab("Residual Distance") +
  ylab("Count") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank())

# Save plot
ggsave(resPlot, file = "residualDistances.png", width = 8, height = 6)



