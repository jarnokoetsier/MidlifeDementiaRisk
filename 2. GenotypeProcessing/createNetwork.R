# Clear workspace and console
rm(list = ls())
cat("\014") 

# Load packages
library(tidyverse)
library(patchwork)
library(data.table)
library(corrr)
library(ggdendroplot)
library(igraph)
library(ggraph)
library(plotly)

# Load PGSs
load("E:/Thesis/EXTEND/df_list.RData")
PGSresult <- df_list$bayesr
PGS_names <- c("AD (w/o APOE)", "AD (stage I)", "AD (stage I-II)", "Adiponectin", "Alchol consumption",
               "ASD", "BD", "BL", "BMI", "BW", "BW (fetal genetics)", "BW (maternal genetics)",
               "CAD", "Childhood obsesity", "Chronotype", "CKD (European)", "CKD (trans-ethnic)",
               "DC2", "EA (Lee, 2014)", "EA (Okbay, 2022)","Extreme BMI", "Extreme height",
               "Extreme HWR", "AD (family history)", "HbA1c", "HC", "HDL", "Infant Head Circumference",
               "LDL", "LST","MDD", "MVPA","Obesity (class I)", "Obsesity (Class II)", "Obesity (Class III)",
               "Overweigth", "PD", "Pubertal Growth", "RA", "SBP (automated reading)", "SBP (manual reading)",
               "SHR", "Sleep duration", "T2D", "TC",
               "TG", "WC", "WHR")

colnames(PGSresult) <- PGS_names

# Calculate correlations
corrDF <- as.data.frame(correlate(PGSresult, diagonal = 0, method = "spearman"))
rownames(corrDF) <- corrDF$term
corrDF <- corrDF[,-1]

# Prepare data for plotting
corrDF1 <- gather(corrDF, key = from, value = weight)
corrDF1$to<- rep(rownames(corrDF), ncol(corrDF))
corrDF1$type <- ifelse(corrDF1$weight < 0, "Negative", "Positive")
corrDF1 <- corrDF1[abs(corrDF1$weight) >= 0.1,]

# Set edges of network
edges <- corrDF1[,c("from", "to", "weight", "type")]

# Set nodes of network
nodes <- data.frame(Name = colnames(corrDF))

# Create network
network <- graph_from_data_frame(edges, directed = FALSE, nodes)

# remove isolated nodes
Isolated <- which(degree(network) == 0)
network_fil <- delete.vertices(network, Isolated)

# Plot network
p <- ggraph(network_fil) + 
  geom_edge_link(aes(color = weight)) + 
  geom_node_point() +
  theme_graph() +
  geom_node_text(aes(label = name), repel=TRUE) +
  scale_edge_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(edge_color = "Spearman\ncorrelation")

# Save plot
ggsave(p, file = "E:/PRS-app/Network.png", width = 10, height = 8)


# Create interactive network using plotly:
L <- layout_nicely(network_fil)
vs <- V(network_fil)
rownames(L) <- names(vs)
Xn <- L[,1]
Yn <- L[,2]

p <- plot_ly(x = ~Xn, y = ~Yn, size = degree(network_fil), color = I("black"),
             mode = "markers", text = vs$name, hoverinfo = "text")

es <- as.data.frame(get.edgelist(network_fil))

Nv <- length(vs)
Ne <- length(es[1]$V1)

edge_shapes <- list()
for(i in 1:Ne) {
  v0 <- es[i,]$V1
  v1 <- es[i,]$V2
  
  edge_shape = list(
    type = "line",
    line = list(color = "#030303", width = 0.3),
    x0 = Xn[v0],
    y0 = Yn[v0],
    x1 = Xn[v1],
    y1 = Yn[v1]
  )
  
  edge_shapes[[i]] <- edge_shape
}

axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)

fig <- layout(
  p,
  title = 'Disease Network',
  shapes = edge_shapes,
  xaxis = axis,
  yaxis = axis
)

fig