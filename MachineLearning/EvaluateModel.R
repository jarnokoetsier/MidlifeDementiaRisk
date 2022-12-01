
load("E:/Thesis/EXTEND/Phenotypes/CAIDE.RData")
# Multiple regression
features <- rownames(dataMatrix_S)[sample(1:nrow(dataMatrix_S),200)]
formula <- paste0("cbind(",paste(features, collapse = ", "),") ~ ", 
                  paste0("0 + ", paste(colnames(CAIDE[,10:16]), collapse = " + ")))

# Make data matrix
dataMatrix <- dataMatrix_S[features, CAIDE$Basename]

# Scale data matrix
dataMatrix_scaled <- t((dataMatrix - rowMeans(dataMatrix))/(apply(dataMatrix,1,sd)))

# Get factprs
factors <- CAIDE[,10:16]
#factors_scaled <- t(t(factors)/(apply(factors,2,max)))


# Add factors
dataMatrix_scaled <- cbind(dataMatrix_scaled,factors)

# Fit model
model <- lm(as.formula(formula), data = as.data.frame(dataMatrix_scaled))
coeff <- coef(model)

# scale coefficients
coeff_scaled <- (coeff - rowMeans(coeff))/(apply(coeff,1,sd))
coeff_scaled <- as.data.frame(t(coeff_scaled))


ggplot(coeff_scaled) +
  geom_point(aes(x = age_c, y = Sex_c, color = BMI_c))


library(ggradar) # install from github
library(scales)


test <- coeff_scaled %>% 
  mutate_at(vars(-age_c), rescale) %>%
  mutate_at(vars(age_c), rescale)

test <- cbind.data.frame(factor(rownames(test)),test)
ggradar(test[1:3,]) +
  theme(legend.position = "none",)



fittedValues <- fitted(model)
residualValues <- residuals(model)

ssr <- sum(colSums(residualValues^2))
sse <- sum(colSums(fittedValues^2))