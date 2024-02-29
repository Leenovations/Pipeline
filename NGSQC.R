library(ggplot2)
library(cowplot)
library(gridExtra)
library(patchwork)
library(ggprism)
#---------------------------------------------------------------#
Sample <- read.table('SampleSheet.txt', sep='\t', header=F)
Name <- Sample[1, 1]
R1 <- Sample[1, 2]
R2 <- Sample[1, 3]
#---------------------------------------------------------------#
Data <- read.table(sprintf('02.BamQC/%s.coverage.txt', Name),
                   sep='\t',
                   header=F)[1:24,]
#---------------------------------------------------------------#
Order <- c('X', 'Y')
for (i in 1:22) {
  Order <- append(Order, i)
}
Order <- c(Order[3:length(Order)], 'X', 'Y')
#---------------------------------------------------------------#
Depth <- ggplot(Data, aes(x=V1, y=V7)) +
  geom_bar(stat="identity", fill='lightcyan3') +
  xlab('') +
  ylab('Mean Depth') +
  scale_x_discrete(limits=Order) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(size=10, face='bold'))
#---------------------------------------------------------------#
Coverage <- ggplot(Data, aes(x=V1, y=V6)) +
  geom_bar(stat="identity", fill='lightgoldenrod2') +
  xlab('') +
  ylab('Coverage') +
  ylim(0,100) +
  scale_x_discrete(limits=Order) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(size=10, face='bold'))
#---------------------------------------------------------------#
#---------------------------------------------------------------#
Data <- read.table(sprintf('02.BamQC/%s.DepthofCoverage.txt', Name),
                   sep='\t',
                   header=F,
                   col.names = c('Depth', 'Count'))
#---------------------------------------------------------------#
DepthOfCoverage <- ggplot(Data, aes(x=Depth, y=Count)) + 
  geom_bar(stat="identity", fill='lavenderblush2', color='black') +
  ggtitle('Depth of Coverage (x)') +
  theme_classic() +
  scale_y_continuous(guide = "prism_offset") +
  scale_x_continuous(guide = "prism_offset", 
                     limits = c(0, 30), 
                     breaks = c(0, 5, 10, 15, 20 ,25, 30),
                     labels = c(0, 5, 10, 15, 20 ,25, 30)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(hjust = 0.5, size=10, face='bold'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.left = element_blank(),
        axis.line.x.bottom = element_line(lineend = 'round'))
#---------------------------------------------------------------#
QCplot <- DepthOfCoverage + Depth + Coverage + plot_layout(ncol = 1)
#---------------------------------------------------------------#
ggsave('03.Plots/QCplot.png',
        height=5,
        width=5,
        plot=QCplot)