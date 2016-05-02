library(ggplot2)
library(gridExtra)

# get inpath and outpath
options(echo=TRUE) 	# if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
pin <- args[1]   	# duplicate histogram
nCpG <- args[2]
nCpG_capture <- args[3]
oi <- args[4]
oi_capture <- args[5]

# pin <- "/share/ClusterScratch/phuluu/captureSeq/15092015/tmp/bedtools/overlap_target/"
# nCpG <- 774294
# nCpG_capture <- 809747
# oi <- oi_capture <- 0

# List all the name of saturate count files
fns <- grep("length.", dir(pin), value=TRUE)
# Read all the saturate count files
tmp <- sapply(fns, function(x){
	sample <- gsub("length.", "", x)
	col <- list()
	col[[sample]] <- read.table(paste0(pin, x), stringsAsFactors=FALSE)[,1]
	return(col)
	}, USE.NAMES=FALSE)
# form table
tb <- sapply(tmp, function(x){
	return(c(length(x), sum(x)))
	})
row.names(tb) <- c("nTarget", "Width of targets(bp)")
tb <- data.frame(tb)
tb["nCpG", "primary"] <- nCpG
tb["nCpG", "capture"] <- nCpG_capture
tb["Internal Overlap", "primary"] <- oi
tb["Internal Overlap", "capture"] <- oi_capture
for (i in 1:dim(tb)[1]){ for (j in 1:dim(tb)[2]){tb[i,j] <- prettyNum(tb[i,j],big.mark=",",scientific=FALSE)}}
tb <- data.frame(t(tb))
# Plot both hist and table
# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tb <- tableGrob(tb, theme=tt)
# Plot chart and table into one object
png(paste0(pin, "metrics_tbl.png"), res=300, width=1840, height=800)
grid.arrange(tb, nrow=1, as.table=TRUE, heights=1)
dev.off()

# convert to long form
width <- do.call(rbind, lapply(names(tmp), function(x){
	return(data.frame(Target=rep(x, length(tmp[[x]])), width=tmp[[x]]))
	}))
# Density plots with semi-transparent fill
p1 <- ggplot(width, aes(x=width, fill=Target)) + geom_density(alpha=.3) + xlab("Width of targets") + ylab("Frequency")
ggsave(filename=paste0(pin, "width.dist.png"))

tb.width <- sapply(names(tmp), function(x){
		return(quantile(tmp[[x]]))
		})
tb.width <- data.frame(t(tb.width))

# Plot both hist and table
# Set theme to allow for plotmath expressions
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tb.width <- tableGrob(tb.width, theme=tt)
# Plot chart and table into one object
png(paste0(pin, "width.dist_tbl.png"), res=300, width=1840, height=1840)
grid.arrange(p1, tb.width, nrow=2, as.table=TRUE, heights=c(3,1))
dev.off()

