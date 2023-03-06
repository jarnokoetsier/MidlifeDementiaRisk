# Load data
load("~/ADNI/ADNI_metadata_classes_640.Rdata")
load("ADNI/X_ADNI.RData")
load("ADNI/lumi_dpval_ADNI.RData")
load("ADNI/MetaData_ADNI.RData")

# Keep midlife samples only:
keepSamples<- MetaData_baseline$Basename[MetaData_baseline$Age <= 75]

# Keep probes with low detection pvalue
keepProbes <- c(removeSamples,rownames(lumi_dpval)[rowSums(lumi_dpval>0.1) == 0])

X_ADNI_M <- log2(X_ADNI/(1-X_ADNI))
X_ADNI_fil <- X_ADNI_M[keepProbes,keepSamples]


predictedScore <- matrix(NA, nrow = ncol(X_ADNI_fil), ncol = length(factors))
colnames(predictedScore) <- factors
rownames(predictedScore) <- colnames(X_ADNI_fil)
for (f in 1:length(factors)){
  f1 <- factors[f]
  model <- allModels[[f1]]
  
  # Get features needed for model fitting
  features <- colnames(model$trainingData)[-10001]
  
  # Get missing features
  missingFeatures <- features[!(features %in% rownames(X_ADNI_fil))]
  
  # Set coefficient of missing features to zero
  addFeatures <- matrix(0,nrow = length(missingFeatures), ncol = ncol(X_ADNI_fil))
  rownames(addFeatures) <- missingFeatures
  colnames(addFeatures) <- colnames(X_ADNI_fil)
  X_ADNI_all <- cbind(X_ADNI_fil, addFeatures)
  
  # Make predictions
  predictedScore[,f] <- predict(model, t(X_ADNI_all))
  
}
predictedAge <- agep(X_ADNI_fil, method='all')

predictedScore$Age <- predictedAge[,1]