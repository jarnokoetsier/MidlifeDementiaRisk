library(tidyverse)
setwd("E:/Thesis/EXTEND/Genotypes")
load("df_Result_PGS.RData")

setwd("E:/Thesis/EXTEND/Phenotypes")
load("CAIDE.Rdata")

CAIDE_fil <- CAIDE[CAIDE$ID %in% rownames(df_Result_PGS),]
rownames(CAIDE_fil) <- CAIDE_fil$ID

PGS_fil <- df_Result_PGS[CAIDE_fil$ID,]

all(rownames(PGS_fil) == CAIDE_fil$ID)


factors <- CAIDE_fil[,10:17]

corMatrix <- matrix(NA, nrow = ncol(PGS_fil), ncol = ncol(factors))
for (i in 1:ncol(factors)){
  corMatrix[,i] <- apply(PGS_fil,2,function(x){cor(x,factors[,i], method = "spearman")})
}

colnames(corMatrix) <- colnames(factors)
rownames(corMatrix) <- colnames(PGS_fil)



# Format data for plotting
plotCor <- gather(as.data.frame(corMatrix))
plotCor$key1 <- rep(rownames(corMatrix), ncol(corMatrix))

# Perform clustering to get sample order
model_rows <- hclust(dist(1-abs(corMatrix)), "ward.D2")
model_cols <- hclust(dist(t(1-abs(corMatrix))), "ward.D2")
order_rows <- model_rows$labels[model_rows$order]
order_cols <- model_cols$labels[model_cols$order]

plotCor$key <- factor(plotCor$key, levels = order_cols)
plotCor$key1 <- factor(plotCor$key1, levels = order_rows)

# Highest correlation for each factor
maxCor <- data.frame(
  Factor = names(apply(corMatrix,2,function(x){which.max(abs(x))})),
  MaxPGS = rownames(corMatrix)[apply(corMatrix,2,function(x){which.max(abs(x))})])

plotCor <- inner_join(plotCor, maxCor, by = c("key" =  "Factor"))
plotCor$MaxPGS[plotCor$MaxPGS != plotCor$key1] <- "No"
plotCor$MaxPGS[plotCor$MaxPGS == plotCor$key1] <- "Yes"

# Main correlation plot
main <- ggplot(plotCor) +
  geom_tile(aes(x = key, y = key1, fill = value, color = MaxPGS), width = 0.8, height = 0.9) +
  scale_fill_gradient2(low = "#000072", mid = "white", high = "red", midpoint = 0,
                       limits = c(-0.5,0.5)) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill = paste0(str_to_title("spearman"),"\nCorrelation"))+
  scale_color_manual(values = c("white","black")) +
  guides(color = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_classic()

# Dendrogram
dendroPlot <- ggplot() +
  geom_tile(data = plotCor, aes(x = as.numeric(key), y = 1), fill = "white", alpha = 0) +
  geom_dendro(model,xlim = c(1,length(unique(plotCor$key)))) +
  theme_void() +
  theme(legend.position = "none") 

# Combine plots
p <- dendroPlot + main +
  plot_layout(nrow = 2, ncol = 1,
              heights = c(1,4))