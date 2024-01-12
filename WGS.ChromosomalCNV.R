library(gplots)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
#------------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
Sample <- args[1]
Path <- getwd()
#-----------------------------------------------------------------------------------#
Data <- read.table(paste0(Path, '/05.SV/', sprintf('%s', Sample), '.Chromosome.CNV.txt'),
                   sep='\t',
                   header=T)
#-----------------------------------------------------------------------------------#
Chromosome <- unique(Data$Chr)
max_values <- c()
median_values <- c()
#-----------------------------------------------------------------------------------#
for (chr in Chromosome) {
  Sub <- subset(Data, Chr == chr)
  max_values <- c(max_values, max(Sub$Order))
  median_values <- c(median_values, median(Sub$Order))
}
max_values <- max_values[1:length(max_values)-1]
#-----------------------------------------------------------------------------------#
labels <- Chromosome
data <- data.frame(x = median_values, y = 1)
#-----------------------------------------------------------------------------------#
CNV <- ggplot(Data, aes(x=Order, y=TPM)) +
  geom_point(size=0.5) +
  ggtitle('Chromosomal copy numbers\n') + 
  xlab('') +
  ylab(Log[10]~NormalizedCNV) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_hline(yintercept=0,
             linetype='dashed',
            color='red',
            alpha=0.5) + 
  geom_vline(xintercept = max_values, 
             linetype = "solid", 
             color = "gray", 
             alpha=0.5) +
  annotate(geom = "text", x = median_values, y = -6, label = unique(Data$Chr), size = 4.5) +
  annotate(geom = "text", x = 2975.5, y = 6, label = 'test', size = 4) +
  coord_cartesian(ylim = c(-5, 5), expand = T, clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(2,2,5,2), "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=15, face='bold'),
        plot.title = element_text(size=20, hjust=0.5, face = 'bold'),
        axis.title = element_text(size = 10),
        axis.text=element_text(color="black"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
#-----------------------------------------------------------------------------------#
ggsave(paste0(Path, '/05.SV/', sprintf('%s', Sample), '.Chromosome.pdf'),
       plot=CNV,
       width=15)
#-----------------------------------------------------------------------------------#