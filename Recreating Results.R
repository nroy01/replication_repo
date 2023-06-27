#read in data as a matrix

cnts <- read.csv('raw_counts_Cooper_data_final.csv')

rownames(cnts) <- cnts$X
cnts$X <- NULL

cnts <- as.matrix(cnts)

column_data <- read.csv('Cooper_Data_Meta.csv')

rownames(column_data) <- column_data$Sample

column_data <- column_data[,c("Tissue","Treatment","Stage")]

column_data$Tissue <- factor(column_data$Tissue)
column_data$Treatment <- factor(column_data$Treatment)
column_data$Stage <- factor(column_data$Stage)
#column_data$Treatment <- NULL
#column_data$Stage <- NULL
#column_data$Sample <- NULL


head(cnts,3)


column_data

all(rownames(column_data) %in% colnames(cnts))

all(rownames(column_data) == colnames(cnts))


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cnts,
                              colData = column_data,
                              design = ~ Tissue)
dds


keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds$Tissue <- relevel(dds$Tissue, ref = "Benign_Tissue")

dds <- DESeq(dds)
res <- results(dds)
res

summary(res)

resultsNames(dds)

'''
resLFC <- lfcShrink(dds, coef="Tissue_Serous_EOC_vs_Benign_Tissue", type="apeglm")
resLFC


resOrdered <- res[order(res$pvalue),]

sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

# (unevaluated code chunk)
library("IHW")
resIHW <- results(dds, filterFun=ihw)
summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult


plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

class(res05$padj)

#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]
'''
library("AnnotationDbi")
library("org.Hs.eg.db")

final_data <- read.csv('Final_Data.csv')

final_data

final_data$symbol <- mapIds(org.Hs.eg.db, keys = final_data$X, column = 'SYMBOL', keytype = 'ENSEMBL')



#rownames(final_data)


#write.csv(res, "/Users/neel/Desktop/R Studio Files\\ Final_Data.csv", row.names=TRUE)

install.packages("devtools")
library(devtools)
install_github("stephens999/ashr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")



render("Replicating Results.Rmd", "html_document")
