
# Install packages
install.packages("prospectr")

# Load data
load("CAIDE.Rdata")
load("methSet_allNorm_fil.RData")
load("cellType.RData")

#*****************************************************************************#
# Split samples with CAIDE1 score
#*****************************************************************************#

# Prepare data
all_X <- as.matrix(t(methSet_allNorm_fil))[CAIDE$Basename,]
all_Y <- CAIDE
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


# Make training and test data
X_train <- t(all_X[c(selectedSamples_male$model, selectedSamples_female$model),])
Y_train <- all_Y[c(selectedSamples_male$model, selectedSamples_female$model),]

X_test <- t(all_X[c(selectedSamples_male$test, selectedSamples_female$test),])
Y_test <- all_Y[c(selectedSamples_male$test, selectedSamples_female$test),]


#*****************************************************************************#
# Check age distribution
#*****************************************************************************#

# Age
plot_age <- ggplot() +
  geom_histogram(data = Y_train, aes(x = Age, fill = "Train"), bins = 30, alpha = 0.5) +
  geom_histogram(data = Y_test, aes(x = Age, fill = "Test"), bins = 30, alpha = 0.5) +
  ylab("Count") +
  xlab("Age") +
  labs(fill = NULL) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(plot_age, file = "Age_TrainingAndTest.png", width = 8, height = 6)


# Sex
table(Y_train$Sex)/sum(table(Y_train$Sex))
table(Y_test$Sex)/sum(table(Y_test$Sex))

# Cell type composition
Cell_train <- Y_train[,17:22]
rownames(Cell_train) <- Y_train$IID
Cell_train_gather <- gather(Cell_train)
Cell_train_gather$ID <- rep(rownames(Cell_train), ncol(Cell_train))

Cell_test <- Y_test[,17:22]
rownames(Cell_test) <- Y_test$IID
Cell_test_gather <- gather(Cell_test)
Cell_test_gather$ID <- rep(rownames(Cell_test), ncol(Cell_test))

Cell_total <- rbind.data.frame(Cell_train_gather, Cell_test_gather)
Cell_total$Data <- c(rep("Train", nrow(Cell_train_gather)),
                     rep("Test", nrow(Cell_test_gather)))


plot_cell <- ggplot(Cell_total) +
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

ggsave(plot_cell, file = "CellType_TrainingAndTest.png", width = 8, height = 6)



# CAIDE score
plot_caide <- ggplot() +
  geom_histogram(data = Y_train, aes(x = CAIDE, fill = "Train"), bins = 10, alpha = 0.5) +
  geom_histogram(data = Y_test, aes(x = CAIDE, fill = "Test"), bins = 10, alpha = 0.5) +
  ylab("Count") +
  xlab("CAIDE1 Score") +
  labs(fill = NULL) +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(plot_caide, file = "CAIDE1_TrainingAndTest.png", width = 8, height = 6)

#******************************************************************************#
# 1.2. Visualize test and training set with PCA
#******************************************************************************#

# Unit scale the training Data
trainingData_scaled <- t((X_train - rowMeans(X_train))/(apply(X_train,1,sd)))

# Make PCA model
pcaList_train <-  prcomp(trainingData_scaled,        
                         retx = TRUE,
                         center = FALSE,
                         scale = FALSE)

# Get the PCA scores of the training data
scores_train <- as.data.frame(pcaList$x)
scores_train$ID <- rownames(scores_train)

# Calculate the explained variance of the PCs
explVar <- round(((pcaList$sdev^2)/sum(pcaList$sdev^2))*100,2)

# Scale test data (using the standard deviation and mean of the training data)
testData_scaled <- t((X_test - rowMeans(X_train))/(apply(X_train,1,sd)))

# Calculate the scores of the test data
scores_test <- as.data.frame(as.matrix(testData_scaled) %*% as.matrix(pcaList$rotation))
scores_test$ID <- rownames(scores_test)

# Combine the scores of test and training data in a single data frame
scores_all <- rbind.data.frame(scores_train, scores_test)
scores_all$Train <- c(rep("Training", nrow(scores_train)), rep("Test", nrow(scores_test)))

# Plot the scores of the training data and the project scores of the test data 
PCA_TrainTest_56 <- ggplot() +
  geom_point(data = scores_all, aes(x = PC5, y = PC6, shape = Train, color = Train), 
             size = 2, alpha = 0.8) +
  scale_shape_manual(values = c(15,17,0,2)) +
  xlab(paste0("PC5 (", explVar[5], "%)")) +
  ylab(paste0("PC6 (", explVar[6], "%)")) +
  labs(title = NULL, 
       caption = "NOTE: The PCA model is constructed using the training data only. The test data is projected.") +
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
  #scale_color_manual(values = rev(RColorBrewer::brewer.pal("Dark2", n = 3)[c(1,2)]))

# Save plot
ggsave(PCA_TrainTest_56, file = "PC5vs6_TrainingAndTest.png", width = 8, height = 6)
