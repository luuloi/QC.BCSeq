library(ggplot2)
library(gridExtra)

# get inpath and outpath
options(echo=TRUE) 	# if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
pin <- args[1]   	# duplicate histogram

# pin <- "/share/ClusterScratch/phuluu/captureSeq/15092015/tmp/bedtools/CpG_target/count.tsv"
df <- read.table(pin)
tb <- data.frame(t(quantile(df$V5)))
row.names(tb) <- c("Number.of.CpGs")
# Histogram overlaid with kernel density curve
p1 <- ggplot(df, aes(x=V5)) + geom_histogram(aes(y=..density..), binwidth=.5, colour="black", fill="white") +
geom_density(alpha=.2, fill="#FF6666") + # Overlay with transparent density plot
xlab("Number of CpGs in each target regions") + ylab("Frequency")
ggsave(filename=gsub(".tsv", ".CpG.png", pin))

# Plot both hist and table
# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tb <- tableGrob(tb, theme=tt)
# Plot chart and table into one object
png(gsub(".tsv", ".CpG.table.png", pin), res=300, width=1840, height=1840)
grid.arrange(p1, tb, nrow=2, as.table=TRUE, heights=c(3,1))
dev.off()
