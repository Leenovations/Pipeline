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
Batch <- read.table(paste0(Path, '/', Sample, '.batch.config'), sep='\t', col.names='Option')
Type <- strsplit(Batch[1,1], '=')[[1]][2]

if (Type=='WGBS'){
  Data <- read.table(paste0(Path, '/04.ChromosomeCNV/', sprintf('%s', Sample), '.Chromosome.CNV.txt'),
                     sep='\t',
                     header=T)
} else if (Type=='WGS'){
  Data <- read.table(paste0(Path, '/05.SV', '/00.ChromosomeCNV/', sprintf('%s', Sample), '.Chromosome.CNV.txt'),
                       sep='\t',
                       header=T)
}
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
max_values <- max_values[1 : length(max_values)-1] # nolint: seq_linter.
#-----------------------------------------------------------------------------------#
labels <- Chromosome
data <- data.frame(x = median_values, y = 1)
#-----------------------------------------------------------------------------------#
model <- lm(Norm ~ poly(Order, 2, raw = TRUE), data = Data)
#-----------------------------------------------------------------------------------#
CNV <- ggplot(Data, aes(x=Order, y=Norm)) +
  geom_point(size=0.5) +
  ggtitle('Chromosomal copy numbers\n', subtitle = sprintf('%s', Sample)) + 
  xlab('') +
  ylab('Normalized Depth') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0), labels=c(0)) +
  geom_hline(yintercept=0,
             linetype='dashed',
            color='red',
            alpha=0.5) + 
  geom_vline(xintercept = max_values, 
             linetype = "solid", 
             color = "gray", 
             alpha=0.5) +
  annotate(geom = "text", x = median_values, y = -0.001, label = unique(Data$Chr), size = 4.5) +
  coord_cartesian(ylim = c(-0.0003, 0.005), expand = T, clip = "off") +
  theme_bw() +
  theme(plot.margin = unit(c(2,2,5,2), "lines"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=15, face='bold'),
        plot.title = element_text(size=20, hjust=0.5, face = 'bold'),
        plot.subtitle = element_text(size=10, hjust=0.5),
        axis.title = element_text(size = 10),
        axis.text=element_text(color="black"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank())
#-----------------------------------------------------------------------------------#
if (Type=='WGBS'){
  ggsave(paste0(Path, '/04.ChromosomeCNV/', sprintf('%s', Sample), '.ChromosomeCNV.pdf'),
         plot=CNV,
         width=15)
  ggsave(paste0(Path, '/04.ChromosomeCNV/', sprintf('%s', Sample), '.ChromosomeCNV.png'),
         plot=CNV,
         width=15)
} else if (Type=='WGS'){
  ggsave(paste0(Path, '/05.SV/00.ChromosomeCNV/', sprintf('%s', Sample), '.ChromosomeCNV.pdf'),
         plot=CNV,
         width=15)
  ggsave(paste0(Path, '/05.SV/00.ChromosomeCNV/', sprintf('%s', Sample), '.ChromosomeCNV.png'),
         plot=CNV,
         width=15)
}