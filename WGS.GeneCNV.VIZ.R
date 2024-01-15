library(ggplot2)
library(cowplot)
library(ggbiplot)
library(gridExtra)
library(patchwork)
#------------------------------------------------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
Sample <- args[1]
Target_gene <- args[2]
Path <- getwd()
Target_path <- paste0(Path, '/05.SV/01.GeneCNV/')
#------------------------------------------------------------------------------#
Coverage <- read.table(paste0(Target_path, sprintf('%s', Sample), '.', sprintf('%s', Target_gene), '.bedcov'),
                       sep='\t',
                       header=FALSE,
                       col.names=c('Chr', 'Start', 'End', 'Gene', 'Exon', 'Strand', 'Coverage'))

Normalized <- read.table(paste0(Target_path, sprintf('%s', Sample), '.', sprintf('%s', Target_gene), '.Norm.CNV.txt'),
                         sep='\t',
                         header=TRUE)
colnames(Normalized) <- c('Chr', 'Start', 'End', 'Gene', 'Exon', 'Strand', 'Norm')

Range <- read.table('/media/src/hg19/04.cnv/NCBI.RefSeq.Selected.Exon.Chr.X.bed',
                    header = FALSE, 
                    col.names=c('Chr', 'Start', 'End', 'Gene', 'Exon', 'Strand'))
#------------------------------------------------------------------------------#
Coverage$Start <- as.character(Coverage$Start)
Normalized$Start <- as.character(Normalized$Start)

Range_sub <- subset(Range, Gene==Target_gene)
Exon_graph <- Range_sub
Exon_graph$median <- apply(Exon_graph[, c("Start", "End")], 1, median)
Range_sub$Start <- as.character(Range_sub$Start)

Exon_graph$Start <- Exon_graph$Start + 5
Exon_graph$End <- Exon_graph$End - 5
Exon_graph$Start <- as.character(Exon_graph$Start)
Exon_graph$End <- as.character(Exon_graph$End)
Exon_graph$median <- as.character(Exon_graph$median)
Exon_location <- Exon_graph$median

if (Exon_graph$Strand[1] == '+') {
  Number <- 1:nrow(Exon_graph)
  Exon_graph$Exon_label <- Number
} else {
  Number <- nrow(Exon_graph):1
  Exon_graph$Exon_label <- Number
}

Xmin <- Exon_graph$Start
Xmax <- Exon_graph$End

COV <- ggplot(Coverage, aes(x=Start, y=Coverage)) +
  geom_line(color='purple', group=1) +
  ggtitle(paste0('\n', Target_gene, '\n')) +
  xlab('') +
  ylab('Coverage') +
  scale_x_discrete(expand = c(0, 0)) +
  geom_vline(xintercept=Range_sub$Start,
             linetype='solid',
             color='gray',
             alpha=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=13),
        plot.title = element_text(size=20, hjust=0.5, face = c('bold.italic')),
        axis.title = element_text(size = 10),
        axis.text=element_text(color="black"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank())

CNV <- ggplot(Normalized, aes(x=Start, y=Norm)) +
  geom_line(color='purple', group=1) +
  xlab('') +
  ylab(Log[2]~NormalizedCNV) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_hline(yintercept=0,
             linetype='dashed',
             color='red',
             alpha=0.5) +
  geom_vline(xintercept=Range_sub$Start,
             linetype='solid',
             color='gray',
             alpha=0.5) +
  theme_bw() +
  theme(plot.margin = unit(c(10,10,30,10), 'mm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=13, face='bold'),
        plot.title = element_text(size=20, hjust=0.5, face = c('bold.italic')),
        axis.title = element_text(size = 10),
        axis.text=element_text(color="black"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank()) +
  annotate(geom = "rect", xmin=Xmin, xmax=Xmax, ymin=-6 , ymax=-7, color="black", fill='lavender') +
  annotate(geom = "text", x = Exon_location, y = -6.5, label = Exon_graph$Exon_label, size = 3) +
  coord_cartesian(ylim = c(-5, 5), expand = T, clip = "off")

Plot <- COV + CNV + plot_layout(ncol = 1, heights = c(4,4))
ggsave(paste0(Target_path, Sample, '.' ,Target_gene, '.cnv.pdf'),
       width = 10,
       plot=Plot)