library(ggplot2)
#------------------------------------------------------------------------------#
Data <- read.table('/media/node03-HDD01/03.WGS/00.RawData/SRR10354232/05.SV/01.GeneCNV/SRR10354232.Norm.CNV.txt',
                   sep='\t',
                   header=T)
Range <- read.table('/media/src/hg19/08.bed/whole.exome.exon.bed',
                    header = FALSE, 
                    col.names=c('Chr', 'Start', 'End', 'Gene', 'range'))

Gene <- unique(Data$Gene)
Gene <- 'ISG15'
for (gene in Gene){
  Norm <- subset(Data, Gene==gene)
  Start_data <- Norm[, c(1,2,4,5)]
  colnames(Start_data) <- c('Chr', 'Start', 'Gene', colnames(Norm)[5])
  End_data <- Norm[, c(1,3,4,5)]
  colnames(End_data) <- c('Chr', 'Start', 'Gene', colnames(Norm)[5])
  New <- rbind(Start_data, End_data)
  New$Start <- as.numeric(New$Start)
  New <- New[order(as.numeric(New$Start)), ]
  Range_sub <- subset(Range, Gene==gene)
  Xmin <- Range_sub$Start 
  Xmax <- Range_sub$End 
  
  CNV <- ggplot(New, aes(x=Start, y=SRR10354232)) +
    geom_line(color='purple') +
    ggtitle(paste0(gene, '\n')) +
    xlab('') +
    ylab(Log[10]~NormalizedCNV) +
    scale_x_discrete(expand = c(0, 0)) +
    geom_hline(yintercept=0,
               linetype='dashed',
               color='red',
               alpha=0.5) +
    geom_vline(xintercept=Range_sub$End,
               linetype='solid',
               color='gray',
               alpha=0.5) +
    theme_bw() +
    theme(plot.margin = unit(c(2,2,10,2), "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(size=15, face='bold'),
          plot.title = element_text(size=20, hjust=0.5, face = c('bold.italic')),
          axis.title = element_text(size = 10),
          axis.text=element_text(color="black"),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank()) +
    annotate(geom = "rect", xmin=Xmin, xmax=Xmax, ymin=-6 , ymax=-6.5, color="black", fill='lavender') +
    coord_cartesian(ylim = c(-5, 5), expand = T, clip = "off")
  
  ggsave(paste0('/media/node03-HDD01/03.WGS/00.RawData/SRR10354232/05.SV/01.GeneCNV/',gene, '.pdf'),
         plot=CNV, width=15)
}