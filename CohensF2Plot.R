# Variables:
# finalModel:   ElasticNet (glmnet) model object
# X_CAIDE1:     methylation data:
#                 -rows: probes/CpGs/features
#                 -columns: samples
# Y_CAIDE1:     meta data: contains values for CAIDE1 factors (column 14-20)

# Load packages
library(tidyverse)
library(patchwork)

# Extract features from model
modelFeatures <- names(coef(finalModel)[-1,1])[abs(coef(finalModel)[-1,1])>0]

# Extract coefficients
modelCoefs <- data.frame(CpG = modelFeatures,
                         coefValue = coef(finalModel)[modelFeatures,1])

# Data matrix with model features only
X_CAIDE1_modelFeatures <- X_CAIDE1[modelFeatures,]

# Check whether samples are in correct order
all(colnames(X_CAIDE1_modelFeatures) == Y_CAIDE1$Basename)

# Make formula for regression model
formula <- paste0("cbind(",paste(modelFeatures, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(Y_CAIDE1[,14:20]), collapse = " + ")))

# Mean center the data
X_coefs_scaled <- (X_CAIDE1_top - rowMeans(X_CAIDE1_top))

# Combine data matrix with independent variables (CAIDE1 factors)
dataMatrix <- cbind(t(X_coefs_scaled),Y_CAIDE1[,14:20])

# Fit model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix))

# Get fitted and residual values
fittedValues <- fitted(model)
residualValues <- residuals(model)

# Calculate R-squared
sse <- colSums(residualValues^2)
ssr <- colSums(fittedValues^2)
Rsquared = ssr/(ssr + sse)

# Get global effect size
globalEffect <- (Rsquared)/(1-Rsquared)

# Calculate Cohen's f2 statistic for each CAIDE1 factor
cohenF <- list()
factors <- colnames(Y_CAIDE1[,14:20])
for (i in 1:length(factors)){
  
  # Formula without factor
  formula <- paste0("cbind(",paste(topFeatures, collapse = ", "),") ~ ", 
                    paste0("0 + ", paste(factors[-i], collapse = " + ")))
  
  # Fit model
  model_i <- lm(as.formula(formula), data = as.data.frame(dataMatrix))
  
  # Get fitted and residual values
  fittedValues <- fitted(model_i)
  residualValues <- residuals(model_i)
  
  # Calculate Rsquared
  sse <- colSums(residualValues^2)
  ssr <- colSums(fittedValues^2) # Data is mean centered, so no need to subtract the mean to get the sum of squares regression
  Rsquared_i <- ssr/(ssr + sse)
  
  # Calculate cohen's f2 statistic (local effect size)
  cohenF[[i]] <- ((Rsquared - Rsquared_i)/(1-Rsquared))
}

# Combine results into data frame
factorNames <- c("Age", "Sex", "Edu", "BP", "BMI", "Chol", "Physical")
plotDF <- data.frame(cohenF = c(unlist(cohenF), globalEffect),
                     Effect = rep(c(factorNames, "Global"), each = nrow(X_coefs_scaled)),
                     CpG = rep(topFeatures,length(factorNames) +1))

# Reorder CAIDE1 factors for plot
plotDF$Effect <- factor(plotDF$Effect, levels = c(factorNames, "Global"))

# Combine with coefficient values in final model
plotDF <- inner_join(plotDF, topCoefs, by = c("CpG" = "CpG"))

# Reorder the CpG sites based on their regression coefficient in final model
plotDF$CpG <- factor(plotDF$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))


# Make main part of the plot (cohen's f2 values)
main <- ggplot(plotDF) +
  geom_bar(aes(y = cohenF, x = CpG, fill = Effect), stat="identity", alpha = 1) +
  facet_grid(rows = vars(Effect)) +
  xlab("Probes") +
  ylab(expression("Cohen's " ~ f^2)) +
  scale_fill_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill= "grey", linewidth = 1, 
                                        linetype="solid"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 13))

# make top part of the plot (regression coefficients) 
topCoefs$CpG <- factor(topCoefs$CpG, levels = unique(arrange(plotDF, coefValue)$CpG))
top <- ggplot(topCoefs) +
  geom_bar(aes(x = CpG, y = coefValue, fill = coefValue), stat = "identity", color = "black") +
  ylab("Coefficients\nFinal Model") +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5), oob = scales::squish) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# Combine main and top plots
p <- top / main +
  plot_layout(heights = c(1,7))

# Save plot
ggsave(p,file = "cohenF_ElasticNetModel_Cor_CAIDE1_all.png", width = 8, height = 8)

