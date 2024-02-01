argvs <- commandArgs( trailingOnly=TRUE );

fn <- as.character(argvs[1]) ;
# res/fstWC-land_cult.tbl
# fn <- 'res/fstWC-land_cult.tbl'
a1 <- read.table(fn, header=T, stringsAsFactors=F)
# CHROM   POS     SVID    SVLEN   WEIR_AND_COCKERHAM_FST
# CLV01_Chr01     34777   v3023_723       -723    0.0945242
# CLV01_Chr01     43095   v3025_1447      1447    0.00114114
a1$CHROM <- gsub('CLV01_', '', a1$CHROM)
cutoff <- quantile(a1$WEIR_AND_COCKERHAM_FST, probs=c(0.99), na.rm=T)[1]
pdf(file=paste0(fn, "-histo.pdf"))
hist(a1$WEIR_AND_COCKERHAM_FST, xlim=c(0,1), breaks=seq(0,1,0.01), xlab= "Fst", ylab= "Frequency", main="")
hist(a1$WEIR_AND_COCKERHAM_FST, xlim=c(0,1), breaks=seq(0,1,0.01), xlab= "Fst", ylab= "Frequency", main=paste0(fn, " cutoff:", round(cutoff, digits=4)))
dev.off()

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(showtext))
showtext_auto()

pdf(file=paste0(fn, "-alongChr.pdf"), height=8.5, width= 11)
# options(bitmapType="cairo")
p <- ggplot(data = a1, aes(x= POS/1e6, y= WEIR_AND_COCKERHAM_FST, colour= factor(CHROM))) +
  ylim(0,1) + labs(x= 'Chromosome', y= 'Fst') +
  geom_point(stat = 'identity', size=1) +
  facet_grid(~CHROM, scales= 'free_x', space= 'free_x', switch= 'x') +
  theme_classic()
print(p)
print(p + geom_hline(yintercept= cutoff, linetype= 'dashed', col= 'red'))
dev.off()

