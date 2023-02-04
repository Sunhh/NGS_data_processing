# Set up parameters.
fn <- 'kmer_cnt_chr/21QDX551/21QDX551_Chr01.tbl'
argvs <- commandArgs( trailingOnly=TRUE );
if ( is.na(argvs[1]) ) {
  message("###########################################################");
  message("Rscript this.R  <kmer_cnt/Chr01.tbl>  <output_prefix>");
  message("######");
  message("  <Chr01.tbl>: is the output of count_kmer_distr.pl\n");
  message("###########################################################");
  q();
}
fn        <- as.character(argvs[1])
opref     <- as.character(argvs[2])

library(ggplot2)
library(dplyr)
library(reshape2)
options(bitmapType="cairo")

# Set up functions.
plot1 <- function(dat, scale_txt='mb') {
  if (scale_txt == 'mb') {
    dvd = 1e6
    seqBy = 0.25
    seqRound = 2
  } else if (scale_txt == 'kb') {
    dvd = 1e3
    seqBy = 10
    seqRound = 0
  } else {
    message("[Err] unknown scale_txt [", scale_txt, "]\n")
    quit()
  }
  dat$position = dat$position/dvd
  minx = min(dat$position)
  maxx = max(dat$position)
  brkx = round(seq(minx, maxx, by= seqBy), digits= seqRound)
  ggplot(dat, aes(x= position, y=value, fill= variable)) +
    geom_area(size= 0.5, colour= 'black', alpha= 1) +
    xlab(paste0('window start (/', scale_txt, ')')) +
    scale_x_continuous(expand=c(0,0), breaks= brkx) +
    scale_y_continuous(expand=c(0,0))
}

# Load data
a1 <- read.table(fn, header=T, stringsAsFactors=F)

# Plot the whole chromosome.
b1 <- round(a1[,-c(1,2,3,4)] / a1[,4] * 100, digits= 1)
b1 <- cbind(position=a1[,c(2)], b1)
c1 <- melt(b1, id= "position")
png(filename= paste0(opref, '.whole.png'), width= 480*10, type= 'cairo-png')
plot1(dat= c1, scale_txt = 'mb')
dev.off()

# Plot first 1 mb
b2 <- b1[b1$position <= 1e6+1, ]
c2 <- melt(b2, id= "position")
png(filename= paste0(opref, '.1st_1mb.png'), width=480*10, type= 'cairo-png')
plot1(dat= c2, scale_txt= 'kb')
dev.off()

# Plot last 1 mb
b3 <- b1[b1$position >= max(b1$position)-1e6+1, ]
c3 <- melt(b3, id= "position")
png(filename= paste0(opref, '.last_1mb.png'), width=480*10, type= 'cairo-png')
plot1(dat= c3, scale_txt= 'kb')
dev.off()

