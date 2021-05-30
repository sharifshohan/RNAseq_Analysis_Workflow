#########RNAseq Analysis Workflow

#At first you need to download the data

#This is where you will ge the data for analysis
#https://github.com/sharifshohan/RNAseq_Analysis_Workflow/blob/main/count_matrix.txt

#This is a mouse expression dataset with control and treatment

# At first I am going to load the data
countdata <- read.csv("count_matrix.txt", header = TRUE, row.names = 1, sep = "\t")
View(countdata)
#or you can view with head command
head(countdata)

##Install and load libraries
#As I have installed all the required libraries before. I am just going to 
#load the libraries here

library(DESeq2)

#Now let me create expermental labels (two conditions)
coldata <- DataFrame(condition = factor(c("ctrl","ctrl","ctrl","treat",
                                          "treat","treat")))
coldata

#Now lets create a DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata,
                              design = ~condition)
dds

##Now lets run DEseq
dds <- DESeq(dds)

#Lets get differentially expressed genes
res <- results(dds)

#Lets look at the table
head(res)

#We have got gene names, log2Foldchange and adjusted p vlaues which are the things
#We will require going forward

#Lets order according to adjusted p values

resOrdered <- res[order(res$padj),]
#Let's see the top of oredered data
head(resOrdered)

#Now for the downstream analysis 
#Lets get differentially expressed genes with specific fold change and padj

sig <- resOrdered[!is.na(resOrdered$padj) &
                        resOrdered$padj < 0.10 &
                        abs(resOrdered$log2FoldChange) >= 1,]
#Lets see whats in sig
head(sig)

###How to create a heatmap from this data
##At first we are going to select genes

selected <- rownames(sig)
selected

#Now lets load libraries for heatmap. There are several ways to do this.
#I am going to show one here

library(RColorBrewer)
library(gplots)

#Lets set color first
hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
#Lets plot the heatmap now. It has some inbuilt functions we are going to use

heatmap.2(log2(counts(dds, normalized = TRUE)[rownames(dds) %in% selected,]),
          col = hmcol, scale = "row",
          Rowv = TRUE, Colv = FALSE,
          dendrogram = "row",
          trace = "none",
          margin = c(4,6), cexRow = 0.5, cexCol = 1, keysize = 1 )

#So we have got a very nice heatmap for our data


#Now we are going to do Gene ontology of the selected genes
#Universe

universe <- rownames(resOrdered)

#Load mouse annotation and ID library
library(org.Mm.eg.db)

#Now convert the gene names to EntrezID

genemap <- select(org.Mm.eg.db, selected, "ENTREZID","SYMBOL")
univmap <- select(org.Mm.eg.db, universe, "ENTREZID","SYMBOL")

#Now lets load GO scoring package

library(GOstats)

#Lets set up the analysis

param<- new("GOHyperGParams", geneIds= genemap, universeGeneIds= univmap,
            annotation = "org.Mm.eg.db", ontology = "BP", pvalueCutoff=0.01,
            conditional=FALSE, testDirection = "over")

#you can ignore the warnings

#Run analysis

hyp <- hyperGTest(param)
#this step takes some time. So do not worry
#Check the summary of the analysis

summary(hyp)

#### thank you for watching. 
