# Clear workspace and console
rm(list = ls())
cat("\014") 

load("CAIDE1_Cor/X_CAIDE1_CorCV.RData")
load("Data/Y_CAIDE1.RData")
load("~/selectedFeatures_v3.RData")

# Convert to M values
X_CAIDE1_M <- log2(X_CAIDE1_CorCV/(1-X_CAIDE1_CorCV))

# Scale data
X_train = t((X_CAIDE1_M -rowMeans(X_CAIDE1_M))/apply(X_CAIDE1_M,1,sd))

# Phenotypes
Y_train = Y_CAIDE1$CAIDE

# Data Matrix
dataMatrix <- cbind.data.frame(X_train,Y_train)

# number of features
nFeatures <- 500

#selectedFeatures <- NULL
nonSelectedFeatures <- setdiff(colnames(X_train), selectedFeatures)
for (i in 1:nFeatures){
  coefficients <- rep(NA,length(nonSelectedFeatures))
  for (f in 1:length(nonSelectedFeatures)){
    # Get new feature to be tested
    newFeature <- nonSelectedFeatures[f]
    
    # Add new feature to already selected features
    formula <- paste0("Y_train ~ ", paste(c(selectedFeatures,newFeature), collapse = " + "))
    
    # fit linear model
    model <- lm(as.formula(formula), data = dataMatrix)
    
    # Get absolute coefficients
    coefficients[f] <- abs(coef(model)[newFeature])
  }
  
  # Select feature with highest absolute coefficient
  selectedFeatures <- c(selectedFeatures, nonSelectedFeatures[which.max(coefficients)])
  
  # Remove selected feature from nonSelectedFeatures object
  nonSelectedFeatures <- nonSelectedFeatures[-which.max(coefficients)]
  save(selectedFeatures, file = "selectedFeatures_v3.RData")
}


