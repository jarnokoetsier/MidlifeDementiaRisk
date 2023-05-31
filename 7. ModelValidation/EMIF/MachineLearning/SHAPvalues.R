
#######################################################################################

# SHAP values: MRSs

#######################################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library("Numero")
library("DALEX") # SHAP values
library(fmsb) # Radar Chart

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")


###############################################################################

# Get SHAP values

###############################################################################

# Load MCI model
load("EMIF/Fit_EMIF_MCI_RF.RData")

# Make explained object
explainer <- DALEX::explain(fit, data = X_train, y = NULL)

# Calculate SHAP values
for (i in 1:nrow(X_test)){
  shap <- predict_parts(explainer = explainer, 
                        new_observation = X_test[i,], 
                        type = "shap",
                        B = 50)
  
  # get mean contribution (i.e., SHAP value)
  shap_fil <- shap %>%
    group_by(variable_name) %>%
    reframe(
      Contribution = mean(contribution))
  
  if (i == 1){
    output <- shap_fil
  } else{
    output <- inner_join(output, shap_fil, by = c("variable_name" = "variable_name"))
  }
  
}
output <- as.data.frame(output)
rownames(output) <- output$variable_name
output <- output[,-1]
colnames(output) <- rownames(X_test)

# Save SHAP values 
save(output, file = "EMIF/output_shap_RF.RData")


###############################################################################

# Compare SHAP values between different models

###############################################################################

#==============================================================================#
# MCI models from EMIF
#==============================================================================#

# Models
models <- c("EN", "sPLS", "RF")

plotDF <- NULL
for (i in 1:length(models)){
  
  # Load SHAP values
  load(paste0("EMIF/output_shap_", models[i], ".RData"))
  
  # Scale SHAP values per individual
  output_scaled <- t(t(output)/colSums(abs(output)))
  
  # Get average SHAP values
  temp <- data.frame(Variable = rownames(output_scaled),
                        value = rowMeans(abs(output_scaled)))
  
  if (i > 1){
    plotDF <- inner_join(plotDF, temp, by = c("Variable" = "Variable"))
  } else{
    plotDF <- temp
  }
}

# Set row names
rownames(plotDF) <- c("Age", "Alcohol Intake", "BMI", "Depression", "Type II Diabetes",
                      "Dietary Intake", "Education", "HDL Chol.", "Heart Disease",
                      "Physical Act.", "Sex", "Smoking", "Syst. Blood Pressure", "Total Chol.")
plotDF <- plotDF[,-1]

# Set column names
colnames(plotDF) <- c("ElasticNet", "sPLS", "Random Forest")

# Prepare data for plotting
data <- as.data.frame(t(plotDF))
data <- rbind(rep(0.5,5) , rep(0,5) , data)

# Color vector
colors_border=rev(c("#EF3B2C","#CB181D", "#99000D") )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart( data  , axistype=0 , 
            #custom polygon
            pcol=colors_border  , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.4,0.1), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

#==============================================================================#
# Epi-CAIDE1 and Epi-LIBRA 
#==============================================================================#

# Models
models <- c("CAIDE1", "LIBRA")

plotDF <- NULL
for (i in 1:length(models)){
  
  # Load SHAP values
  load(paste0("EMIF/output_shap_", models[i], ".RData"))
  
  # Scale SHAP values per individual
  output_scaled <- t(t(output)/colSums(abs(output)))
  
  # Get average SHAP values
  temp <- data.frame(Variable = rownames(output_scaled),
                     value = rowMeans(abs(output_scaled)))
  
  if (i > 1){
    plotDF <- inner_join(plotDF, temp, by = c("Variable" = "Variable"))
  } else{
    plotDF <- temp
  }
}

# Set row names
rownames(plotDF) <- c("Age", "Alcohol Intake", "BMI", "Depression", "Type II Diabetes",
                      "Dietary Intake", "Education", "HDL Chol.", "Heart Disease",
                      "Physical Act.", "Sex", "Smoking", "Syst. Blood Pressure", "Total Chol.")
plotDF <- plotDF[,-1]

# Set column names
colnames(plotDF) <- c("CAIDE1", "LIBRA")

# Prepate data for plotting
data <- as.data.frame(t(plotDF))
data <- rbind(rep(0.5,5) , rep(0,5) , data)

# Color vector
colors_border=rev(c("#2171B5","#084594"))

# plot with default options:
radarchart( data  , axistype=0 , 
            #custom polygon
            pcol=colors_border  , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8,
            #custom labels
            vlcex=0.8 
)



#######################################################################################

# SHAP values: MRSs + CSF biomarkers

#######################################################################################

# Load packages
library(tidyverse)
library(caret)
library(glmnet)
library(spls)
library(ranger)
library("Numero")
library("DALEX") # SHAP values
library(fmsb) # Radar Chart

# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load data
load("EMIF/X_train_EMIF.RData")
load("EMIF/Y_train_EMIF.RData")
load("EMIF/X_test_EMIF.RData")
load("EMIF/Y_test_EMIF.RData")

# Prepare data
load("EMIF/metaData_fil.RData")
rownames(metaData_fil) <- metaData_fil$X
CSFbio <- metaData_fil[,c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "Age")]
colnames(CSFbio) <- c("Ptau_ASSAY_Zscore", "Ttau_ASSAY_Zscore", "AB_Zscore", "ChrAge")
samples <- rownames(CSFbio)[(!is.na(CSFbio$Ptau_ASSAY_Zscore)) & 
                              (!is.na(CSFbio$AB_Zscore)) &
                              (!is.na(CSFbio$Ttau_ASSAY_Zscore))]

# Add CSF info
Y_train <- Y_train[intersect(samples, rownames(Y_train)),]
Y_test <- Y_test[intersect(samples, rownames(Y_test)),]

X_train <- cbind.data.frame(X_train[rownames(Y_train),], CSFbio[rownames(Y_train),])
X_test <- cbind.data.frame(X_test[rownames(Y_test),], CSFbio[rownames(Y_test),])


###############################################################################

# Get SHAP values

###############################################################################

# Load MCI model
load("EMIF/Fit_EMIF_CI_sPLS_CSFbio.RData")

# Make explained object
explainer <- DALEX::explain(fit, data = X_train, y = NULL)

# Calculate SHAP values
for (i in 1:nrow(X_test)){
  shap <- predict_parts(explainer = explainer, 
                        new_observation = X_test[i,], 
                        type = "shap",
                        B = 50)
  
  # get mean contribution (i.e., SHAP value)
  shap_fil <- shap %>%
    group_by(variable_name) %>%
    reframe(
      Contribution = mean(contribution))
  
  if (i == 1){
    output <- shap_fil
  } else{
    output <- inner_join(output, shap_fil, by = c("variable_name" = "variable_name"))
  }
  
}
output <- as.data.frame(output)
rownames(output) <- output$variable_name
output <- output[,-1]
colnames(output) <- rownames(X_test)

# Save SHAP values 
save(output, file = "EMIF/output_shap_sPLS_CSFbio.RData")


###############################################################################

# Compare SHAP values between different models

###############################################################################

#==============================================================================#
# MCI models from EMIF
#==============================================================================#

# Models
models <- c("EN", "sPLS", "RF")

plotDF <- NULL
for (i in 1:length(models)){
  
  # Load SHAP values
  load(paste0("EMIF/output_shap_", models[i], "_CSFbio.RData"))
  
  # Scale SHAP values per individual
  output_scaled <- t(t(output)/colSums(abs(output)))
  
  # Get average SHAP values
  temp <- data.frame(Variable = rownames(output_scaled),
                     value = rowMeans(abs(output_scaled)))
  
  if (i > 1){
    plotDF <- inner_join(plotDF, temp, by = c("Variable" = "Variable"))
  } else{
    plotDF <- temp
  }
}
plotDF <- plotDF[plotDF$Variable != "ChrAge",]

# Set row names
rownames(plotDF) <- c("Age", "Alcohol Intake", "BMI", "Depression", "Type II Diabetes",
                      "Dietary Intake", "Education", "HDL Chol.", "Heart Disease",
                      "Physical Act.", "Sex", "Smoking", "Syst. Blood Pressure", "Total Chol.")
plotDF <- plotDF[,-1]

# Set column names
colnames(plotDF) <- c("ElasticNet", "sPLS", "Random Forest")

# Prepare data for plotting
data <- as.data.frame(t(plotDF))
data <- rbind(rep(0.5,5) , rep(0,5) , data)

# Color vector
colors_border=rev(c("#EF3B2C","#CB181D", "#99000D") )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart( data  , axistype=0 , 
            #custom polygon
            pcol=colors_border  , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,0.4,0.1), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)



