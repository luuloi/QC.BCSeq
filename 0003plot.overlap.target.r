library(ggplot2)

# get inpath and outpath
options(echo=TRUE) 	# if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
pin <- args[1]

# pin <- "/share/ClusterScratch/phuluu/captureSeq/15092015/tmp/bedtools/overlap_target/count.tsv"
df <- read.table(pin)

# plot
ggplot(data=df, aes(x=V2, y=V1)) +
geom_point(size = 4) + # Shape depends on cond
xlab("Times of overlapping with capture regions") + ylab("Number of primary regions")
ggsave(filename=gsub(".tsv", ".png", pin))

# plot
df$V3 <- 100*df$V1/sum(df$V1) 
ggplot(data=df, aes(x=V2, y=V3)) +
geom_point(size = 4) + # Shape depends on cond
xlab("Times of overlapping with capture regions") + ylab("Number of primary regions (%)")
ggsave(filename=gsub(".tsv", "_percent.png", pin))
