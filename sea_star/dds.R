source("http://bioconductor.org/workflows/R")
workflowInstall("rnaseqGene")

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

browseVignettes("DESeq2")

install.packages("ggplot2")

library("DESeq2")
library("ggplot2")

countsTable <- read.delim('countsdata_trim.txt', header = TRUE, stringsAsFactors = TRUE, row.names = 1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header = TRUE, stringsAsFactors = TRUE, row.names = 1)
head(conds)
colData <- as.data.frame(conds)
head(colData)


dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~day + location + health)

dim(dds)

dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)

colSums(counts(dds))
hist(colSums(counts(dds)), breaks = 80, xlim=c(0,max(colSums(counts(dds)))))

colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S"))

dds <- DESeq(dds)

summary(res)

d <- plotCounts(dds, gene="TRINITY_DN46124_c1_g2_TRINITY_DN46124_c1_g2_i6_g.21322_m.21322", intgroup=(c("location","health","day")), returnData=TRUE)
d
p <- ggplot(d, aes(x= health, y=count, shape = location, colour = day, fill=location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100,4000)) + scale_x_discrete(limits=c("H","S")) + scale_shape_manual(values=c(21,24))
p

p <- ggplot(d, aes(x= health, y=count, shape = location, fill=location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size = 4, alpha=0.9) + scale_y_log10(breaks=c(25,100,4000)) + scale_x_discrete(limits=c("H","S")) + scale_shape_manual(values=c(21,24))
p

ggsave("dot_plot-TRINITY_DN46124_c1_g2_TRINITY_DN46124_c1_g2_i6_g.21322_m.21322.png", p, width=8, height=4, dpi=300)

