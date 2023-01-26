# Greedy sparse regression
library(matlib)
Ytrain
Xtrain <- t(dataMatrix)
n <- nrows(dataMatrix)
Iexact <- NULL

I <- 1:n
I <- setdiff(I, Iexact)

# While the number of non-zero features are below sparsity level..
j = 0
while(length(Iexact) < sparsity){
  j = j + 1
  # L02: regression coefficients: set to zero
  L02 <- rep(0,length(I))
  
  # For each feature that is not yet selected....
  for (i in 1:length(I)){
    
    # Get all non-zero coefficients and the "new" feature
    Iloop <- c(Iexact, I[i])
    
    # Calculate regression coefficient of "new" feature
    L02[i] <- t(Y_train)%*%X_train[,Iloop]%*%inv(t(X_train[,Iloop])%*%X_train[,Iloop])%*%(t(X_train[,Iloop])%*%Y_train)
  }
  
  # select feature with maximum regression coefficient
  Imax <- which.max(L02)
  
  # Add selected feature to coefficient list
  Iexact <- c(Iexact, I[Imax])
  
  if (j %% 10 == 0){
    # train model
    X_CV <- X_train[,Iexact]
    Y_CV <- Y_train
  }
  
  # Remove feature from not-yet-selected features
  I <- I[-Imax]
}
x <- rep(0,n)
x[idx] <- inv(t(X_train[,Iexact])%*%X_train[,Iexact])%*%t(X_train[,Iexact])%*%Y_train
