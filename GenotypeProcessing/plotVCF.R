# Load packages
library(tidyverse)
library(ggpubr)
library(grid)

# Set working directory
setwd("E:/Thesis/EXTEND/Genotypes/Plots")

# Load data
load("gt_results.RData")

# Put results into data frame

# Get phased genotypes
all(names(gt_rs429358) == names(gt_rs7412))

GT_rs429358 <- unlist(lapply(str_split(gt_rs429358, ":"), `[[`, 1))[-1]
GT_rs429358 <- str_replace_all(GT_rs429358, "0", "T")
GT_rs429358 <- str_replace_all(GT_rs429358, "1", "C")

GT_rs7412 <- unlist(lapply(str_split(gt_rs7412, ":"), `[[`, 1))[-1]
GT_rs7412 <- str_replace_all(GT_rs7412, "0", "C")
GT_rs7412 <- str_replace_all(GT_rs7412, "1", "T")


# Get haplotypes
hap1 <- paste0(unlist(lapply(str_split(GT_rs429358, "|"), `[[`, 2)),
               unlist(lapply(str_split(GT_rs7412, "|"), `[[`, 2)))

hap2 <- paste0(unlist(lapply(str_split(GT_rs429358, "|"), `[[`, 4)),
               unlist(lapply(str_split(GT_rs7412, "|"), `[[`, 4)))

e1 <- str_detect(hap1, "CT") + str_detect(hap2, "CT")
e2 <- str_detect(hap1, "TT") + str_detect(hap2, "TT")
e3 <- str_detect(hap1, "TC") + str_detect(hap2, "TC")
e4 <- str_detect(hap1, "CC") + str_detect(hap2, "CC")


P_rs429358 <- unlist(lapply(str_split(gt_rs429358, ":"), `[[`, 4))[-1]
P_rs429358_00 <- as.numeric(unlist(lapply(str_split(P_rs429358, ","), `[[`, 1)))
P_rs429358_01 <- as.numeric(unlist(lapply(str_split(P_rs429358, ","), `[[`, 2)))
P_rs429358_11 <- as.numeric(unlist(lapply(str_split(P_rs429358, ","), `[[`, 3)))

P_rs7412 <- unlist(lapply(str_split(gt_rs7412, ":"), `[[`, 4))[-1]
P_rs7412_00 <- as.numeric(unlist(lapply(str_split(P_rs7412, ","), `[[`, 1)))
P_rs7412_01 <- as.numeric(unlist(lapply(str_split(P_rs7412, ","), `[[`, 2)))
P_rs7412_11 <- as.numeric(unlist(lapply(str_split(P_rs7412, ","), `[[`, 3)))

resultsDF <- data.frame(
  SampleID = rep(names(gt_rs7412)[-1],3),
  GT_rs429358 = rep(GT_rs429358,3),
  P_rs429358 = c(P_rs429358_00,P_rs429358_01,P_rs429358_11),
  GT_P_rs429358 = c(rep("0|0",length(GT_rs429358)), 
                    rep("0|1",length(GT_rs429358)), 
                    rep("1|1",length(GT_rs429358))),
  GT_rs7412 = rep(GT_rs7412,3),
  P_rs7412 = c(P_rs7412_00,P_rs7412_01,P_rs7412_11),
  GT_P_rs7412 = c(rep("0|0",length(GT_rs7412)), 
                    rep("0|1",length(GT_rs7412)), 
                    rep("1|1",length(GT_rs7412)))
)



resultsDF <- data.frame(
  SampleID = rep(names(gt_rs7412)[-1],6),
  GT = c(rep(GT_rs429358,3),rep(GT_rs7412,3)),
  P = c(P_rs429358_00,P_rs429358_01,P_rs429358_11,
        P_rs7412_00,P_rs7412_01,P_rs7412_11),
  GT_P = c(rep("T/T",length(GT_rs429358)), 
           rep("T/C",length(GT_rs429358)), 
           rep("C/C",length(GT_rs429358)),
           rep("C/C",length(GT_rs7412)), 
           rep("C/T",length(GT_rs7412)), 
           rep("T/T",length(GT_rs7412))),
  SNP = c(rep("rs429358",3*length(GT_rs429358)),
          rep("rs7412",3*length(GT_rs7412))),
  e1 = rep(e1,6),
  e2 = rep(e2, 6),
  e3 = rep(e3,6),
  e4 = rep(e4,6)
)



main1 <- ggplot(resultsDF[resultsDF$SNP == "rs429358",]) +
  geom_tile(aes(x = GT_P, y = fct_reorder(SampleID,e4), fill =  P)) +
  scale_fill_viridis_c() +
  labs(fill = "Posterior\nProbability") +
  ggtitle("rs429358") +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

main2 <- ggplot(resultsDF[resultsDF$SNP == "rs7412",]) +
  geom_tile(aes(x = GT_P, y = fct_reorder(SampleID,e4), fill =  P)) +
  scale_fill_viridis_c() +
  labs(fill = "Posterior\nProbability") +
  ggtitle("rs7412") +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

main3 <- ggplot(resultsDF) +
  geom_tile(aes(x = "\u03b51",y = fct_reorder(SampleID,e4), fill = as.factor(e1))) +
  geom_tile(aes(x = "\u03b52",y = fct_reorder(SampleID,e4), fill = as.factor(e2))) +
  geom_tile(aes(x = "\u03b53",y = fct_reorder(SampleID,e4), fill = as.factor(e3))) +
  geom_tile(aes(x = "\u03b54",y = fct_reorder(SampleID,e4), fill = as.factor(e4))) +
  ylab("Samples") +
  ggtitle("APOE status") +
  scale_fill_brewer(palette = "Reds") +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

p <- ggarrange(main3, main1, main2,
          ncol = 3, 
          nrow = 1,
          heights = c(1,1),
          widths = c(5,5,5),
          align = "v",
          common.legend = FALSE)

ggsave(p, file = "APOEstatus.png", height = 8, width = 12)


# Get legends

legendPlot <- ggplot(resultsDF[resultsDF$SNP == "rs429358",]) +
  geom_tile(aes(x = GT_P, y = fct_reorder(SampleID,e4), fill =  P)) +
  scale_fill_viridis_c() +
  labs(fill = "Posterior\nProbability") +
  ggtitle("rs429358") +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))

legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)


legendPlot <-ggplot(resultsDF) +
  geom_tile(aes(x = "\u03b51",y = fct_reorder(SampleID,e4), fill = as.factor(e1))) +
  geom_tile(aes(x = "\u03b52",y = fct_reorder(SampleID,e4), fill = as.factor(e2))) +
  geom_tile(aes(x = "\u03b53",y = fct_reorder(SampleID,e4), fill = as.factor(e3))) +
  geom_tile(aes(x = "\u03b54",y = fct_reorder(SampleID,e4), fill = as.factor(e4))) +
  ylab("Samples") +
  ggtitle("APOE status") +
  labs(fill = "APOE\nstatus") +
  scale_fill_brewer(palette = "Reds", labels = c(0,1,2)) +
  theme_void() +
  theme(axis.text.x = element_text(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5,
                                  face = "bold",
                                  size = 16))
legend <- cowplot::get_legend(legendPlot)
grid.newpage()
grid.draw(legend)

            
            
            
            
          