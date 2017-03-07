library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)
colDataINT<-subset(colData, colData$location=="int")
colDataSUB<-subset(colData, colData$location=="sub")
                   
countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)
                   
countDataINT<-countData[, which(colnames(countData) %in% row.names(colDataINT))]
countDataSUB<-countData[, -which(colnames(countData) %in% row.names(colDataSUB))]
dim(countDataINT)
dim(countDataSUB)

#Model 1 - FULL
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ location + health)

dim(dds)
#[1] 26550    65

dds <- dds[ rowSums(counts(dds)) > 100, ]
dim(dds)
#[1] 13334    65

colData(dds)$health <- factor(colData(dds)$health, levels=c("H","S")) #sets that "healthy is the reference

dds <- DESeq(dds) 

res <- results(dds)
res <- res[order(res$padj),]
head(res)

summary(res)

#out of 13318 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 70, 0.53% 
#LFC < 0 (down)   : 18, 0.14% 
#outliers [1]     : 424, 3.2% 
#low counts [2]   : 10968, 82% 
#(mean count < 46)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#Model 2 INTERTIDAL

ddsINT <- DESeqDataSetFromMatrix(countData = countDataINT, colData = colDataINT, design = ~ health)

dim(ddsINT)
#[1] 26550    39

ddsINT <- ddsINT[ rowSums(counts(ddsINT)) > 100, ]
dim(ddsINT)
#[1] 12082    39 

colData(ddsINT)$health <- factor(colData(ddsINT)$health, levels=c("H","S")) #sets that "healthy is the reference

ddsINT <- DESeq(ddsINT) 

resINT <- results(ddsINT)
resINT <- resINT[order(resINT$padj),]
head(resINT)

summary(resINT)
#out of 12062 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 82, 0.68% 
#LFC < 0 (down)   : 8, 0.066% 
#outliers [1]     : 0, 0% 
#low counts [2]   : 11010, 91% 
#(mean count < 77)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#Model 3 SUBTIDAL

ddsSUB <- DESeqDataSetFromMatrix(countData = countDataSUB, colData = colDataSUB, design = ~ health)

dim(ddsSUB)
#[1] 26550    26

ddsSUB <- ddsSUB[ rowSums(counts(ddsSUB)) > 100, ]
dim(ddsSUB)
#[1] 11507    26 # at > 100; we loose many fewer genes

colData(ddsSUB)$health <- factor(colData(ddsSUB)$health, levels=c("H","S")) #sets that "healthy is the reference

ddsSUB <- DESeq(ddsSUB) 

resSUB <- results(ddsSUB)
resSUB <- resSUB[order(resSUB$padj),]
head(resSUB)

summary(resSUB)
#out of 11499 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)     : 11, 0.096% 
#LFC < 0 (down)   : 196, 1.7% 
#outliers [1]     : 772, 6.7% 
#low counts [2]   : 1793, 16% 
#(mean count < 6)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results

#Plot counts 

d <- plotCounts(dds, returnData=TRUE)
d

norm.counts <- counts(dds, normalized=TRUE)
dim(norm.counts)

############## PCA plots
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("health","location"))

vsdINT <- varianceStabilizingTransformation(ddsINT, blind=FALSE)
plotPCA(vsdINT, intgroup=c("health"))

vsdSUB <- varianceStabilizingTransformation(ddsSUB, blind=FALSE)
plotPCA(vsdSUB, intgroup=c("health"))

