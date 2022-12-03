
library(tidyverse)
#*****************************************************************************#
# Make all data combinations
#*****************************************************************************#
# Make training and test data
load("CAIDE.Rdata")
load("CAIDE2.Rdata")
load("EPILIBRA.Rdata")
load("metaData_ageFil.RData")
load("methSet_allNorm_fil.RData")
load("cellType.RData")
load("gt_results.RData")
load("TestTrain.RData")

# All
all_X <- as.matrix(t(methSet_allNorm_fil))
all_Y <- dat[,c("ID", "Basename", "Position", "Plate", "Age", "Sex")]
all_Y <- inner_join(all_Y, cellType, by = c("Basename" = "ID"))
all(rownames(all_X) == all_Y$Basename)

# CAIDE1
length(intersect(CAIDE$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_CAIDE1 <- t(all_X[setdiff(CAIDE$Basename, TestTrain$Test),])
Y_CAIDE1 <- all_Y[all_Y$Basename %in% setdiff(CAIDE$Basename, TestTrain$Test),]
Y_CAIDE1 <- inner_join(Y_CAIDE1, CAIDE[,c(1,10:17)], by = c("ID" = "ID"))
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
Y_LIBRA <- inner_join(Y_LIBRA, EPILIBRA[,c(1,10:21)], by = c("ID" = "ID"))
all(colnames(X_LIBRA) == Y_LIBRA$Basename)

# Non-test
length(intersect(all_Y$Basename, TestTrain$Test)) == length(TestTrain$Test)
X_nonTest <- t(all_X[setdiff(all_Y$Basename, TestTrain$Test),])
Y_nonTest <- all_Y[all_Y$Basename %in% setdiff(all_Y$Basename, TestTrain$Test),]
all(colnames(X_nonTest) == Y_nonTest$Basename)

# Test
X_test <- t(all_X[TestTrain$Test,])
Y_test <- all_Y[all_Y$Basename %in% TestTrain$Test,]
Y_test <- inner_join(Y_test, EPILIBRA[,c(1,10:21)], by = c("ID" = "ID"))
Y_test <- inner_join(Y_test, CAIDE[,c(1,10:17)], by = c("ID" = "ID"))

X_test <- X_test[,Y_test$Basename]
all(colnames(X_test) == Y_test$Basename)

#*****************************************************************************#
# Variance selection
#*****************************************************************************#

cpg_var <- apply(X_nonTest, 1, var)
cpg_selected_var <- names(tail(sort(cpg_var), 10000))
rm(cpg_var)
#all(rownames(X_nonTest) %in% probe_annotation$ID)

X_nonTest_var <- X_nonTest[cpg_selected_var, ]
X_CAIDE1_var <- X_CAIDE1[cpg_selected_var, ]
X_LIBRA_var <- X_LIBRA[cpg_selected_var, ]
X_test_var <- X_test[cpg_selected_var, ]

save(X_nonTest_var, file = "X_nonTest_var.RData")
save(X_CAIDE1_var, file = "X_CAIDE1_var.RData")
save(X_LIBRA_var, file = "X_LIBRA_var.RData")
save(X_test_var, file = "X_test_var.RData")

rm(X_nonTest_var)
rm(X_CAIDE1_var)
rm(X_LIBRA_var)
rm(X_test_var)

save(Y_nonTest, file = "Y_nonTest.RData")
save(Y_CAIDE1, file = "Y_CAIDE1.RData")
save(Y_LIBRA, file = "Y_LIBRA.RData")
save(Y_test, file = "Y_test.RData")
#*****************************************************************************#
# S-score selection
#*****************************************************************************#

calculate_S <- function(x){
  S = abs(mean(x)- 0.5)/var(x)
}

cpg_S <- apply(X_nonTest, 1, calculate_S)
cpg_selected_S <- names(tail(sort(cpg_S),10000))
rm(cpg_S)

X_nonTest_S <- X_nonTest[cpg_selected_S, ]
X_CAIDE1_S <- X_CAIDE1[cpg_selected_S, ]
X_LIBRA_S <- X_LIBRA[cpg_selected_S, ]
X_test_S <- X_test[cpg_selected_S, ]

save(X_nonTest_S, file = "X_nonTest_S.RData")
save(X_CAIDE1_S, file = "X_CAIDE1_S.RData")
save(X_LIBRA_S, file = "X_LIBRA_S.RData")
save(X_test_S, file = "X_test_S.RData")

rm(X_nonTest_S)
rm(X_CAIDE1_S)
rm(X_LIBRA_S)
rm(X_test_S)