library(edgeR)
setwd("d:/RNA-Seq/trimmed-fastq/featureCounts")
allCounts = NULL
#Loop through sample names
for (sample in c('CD24Mid-1', 'CD24Neg-2', 'CD24Mid-3', 'CD24Neg-4', 'CD24Mid-5', 'CD24Neg-6', 'CD24Mid-7', 'CD24Neg-8')) {
  #Define filename from sample name
  filename = paste(sample, '.counts2.txt', sep='')
  #Read feature counts result
  counts = read.table(filename, header=TRUE)
  #Remove all other columns except geneid and counts
  counts = counts[, -c(2:6), drop=FALSE]
  #Rename columns
  colnames(counts) = c('Geneid', sample)
  #Merge across samples
  if (is.null(allCounts)) {
    allCounts = counts #No need to merge for 1st sample
  } else {
    allCounts = merge(allCounts, counts, by="Geneid")
  }
}

#Set row names to Gene IDS
row.names(allCounts) = allCounts[,1]
#Drop redundant Gene ID columns
allCounts = allCounts[,-1]

#Read into DGEList for edgeR
counts = DGEList(counts = allCounts)
group <- factor(c(1, 2, 1, 2, 1, 2, 1, 2))
design <- model.matrix(~group)
counts = DGEList(counts = allCounts, group = group)
counts
design

#Filter out lowly expressed genes
keep <- filterByExpr(counts)
counts <- counts[keep, , keep.lib.sizes=FALSE]
summary(keep)

counts$samples

#Normalise library size
counts <- calcNormFactors(counts)
counts$samples
plotMD(cpm(counts, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
plotMD(cpm(counts, log=TRUE), column=2)
abline(h=0, col="red", lty=2, lwd=2)
plotMD(cpm(counts, log=TRUE), column=3)
abline(h=0, col="red", lty=2, lwd=2)

#Data exploration
colours <- rep(c("blue", "red"), 2)
plotMDS(counts, col=colours[group], method="bcv")
plotMDS(counts, col=colours[group], main = "MDS plot of gene counts")
#Export allCounts to excel
write.csv(allCounts, "d:/RNA-Seq/trimmed-fastq/featureCounts/allCounts.csv")

#Estimating common dispersion
counts = estimateDisp(counts, design)
plotBCV(counts, main = "BCV plot of gene counts")
plotMDS(counts, col=colours[group], main = "MDS plot of gene counts")
#Exact test
et <- exactTest(counts)
topresults <- topTags(et)

#Set p-value threshold
pvalue = 0.05
#Set logFC threshold
logfc = 0
#For each gene set 0 for no change, 1 for up and -1 for down
etAdjusted = decideTestsDGE(et, adjust.method="fdr", p.value=pvalue, lfc=logfc)
#Convert to vector for adding to dataframe
etAdjusted = as.vector(etAdjusted)
#Retrieve all gene information in table format
allGenes = et$table
#Add column indicating whether significant up(1), down(-1) or no change (0)
allGenes$change = etAdjusted
#Gene list of up and down
up <- row.names(allGenes[allGenes$change == 1,])
down <- row.names(allGenes[allGenes$change == -1,])

#Export gene list
write.csv(up, "d:/RNA-Seq/trimmed-fastq/featureCounts/up-counts.csv")
write.csv(down, "d:/RNA-Seq/trimmed-fastq/featureCounts/down-counts.csv")
write.csv(allGenes, "d:/RNA-Seq/trimmed-fastq/featureCounts/allGenes.csv")


summary(etAdjusted)
summary(allGenes)
summary(allGenes$change)

#Goanalysis
genechange <- c(up, down)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GO.db")
library("GO.db")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
#Converting ensembl gene id to entrez gene id
convert2EntrezID <- function(IDs, orgAnn, ID_type="ensembl_gene_id")
{
  GOgenome = sub(".db","",orgAnn)
  if (ID_type == "ensembl_gene_id")
  {
    orgAnn <- get(paste(GOgenome,"ENSEMBL2EG", sep=""))
  }
  else if (ID_type == "gene_symbol")
  {
    orgAnn <- get(paste(GOgenome,"SYMBOL2EG", sep=""))
  }
  else if (ID_type == "refseq_id")
  {
    orgAnn <- get(paste(GOgenome,"REFSEQ2EG", sep=""))
  }
  else
  {
    stop("Currently only the following type of IDs are supported: ensembl_gene_id, refseq_id and gene_symbol!")
  }
  if(!is(orgAnn, "AnnDbBimap"))
  {
    stop("orgAnn is not a valid annotation dataset! For example, orgs.Hs.eg.db package for human
			 and the org.Mm.eg.db package for mouse.")
  }
  xx <- as.list(orgAnn)
  #IDs = unique(IDs[!is.na(IDs)])
  x = do.call(rbind, lapply(IDs,function(x1)
  {
    r= xx[names(xx)==x1]
    if (length(r) >0)
    {
      r[[1]][1]
    } else {
      c('')
    }
  }))
  # unique(x)
  x
}
genenames = row.names(et$table)
genes_entrez = convert2EntrezID(IDs=genenames, orgAnn="org.Mm.eg.db",
                                ID_type="ensembl_gene_id")

length(genes_entrez)
length(genenames)
go = goana(et, geneid = genes_entrez, FDR = 0.05, species = "Mm")
topGO(go, sort = "up")
topGO(go, sort = "down")

#KEGG pathway analysis
keg = kegga(et, geneid = genes_entrez, FDR = 0.05, species = "Mm")
topKEGG(keg, sort = "up")
topKEGG(keg, sort = "down")

plotMD(et, main = "MD plot of DGE between CD24Neg samples versus CD24Mid samples", legend = "bottomright", ylab = "LogFC")

#Exporting GO and KEGG analyses
write.csv(go, "d:/RNA-Seq/trimmed-fastq/featureCounts/GO_analysis.csv")
write.csv(keg, "d:/RNA-Seq/trimmed-fastq/featureCounts/KEGG_analysis.csv")
MytopGO <- topGO(go, sort='up')
MybottomGO <- topGO(go, sort = "down")
write.csv(MytopGO, "d:/RNA-Seq/trimmed-fastq/featureCounts/topGO_analysis.csv")
write.csv(MybottomGO, "d:/RNA-Seq/trimmed-fastq/featureCounts/topGOdown_analysis.csv")
MytopKEGG <- topKEGG(keg, sort = "up")
MybottomKEGG <- topKEGG(keg, sort = "down")
write.csv(MytopKEGG, "d:/RNA-Seq/trimmed-fastq/featureCounts/topKEGG_analysis.csv")
write.csv(MybottomKEGG, "d:/RNA-Seq/trimmed-fastq/featureCounts/topKEGGdown_analysis.csv")

#Top significantly up/down DE genes
library(edgeR)
topTags(et)
topTags(et, sort='p.value')
topTags(et, n=20, sort='p.value')

#Top highest logFC genes
topTags(et, n=10, sort='logFC')
all.res <- as.data.frame(topTags(et, n=Inf))
up.res <- all.res[all.res$logFC > 0,]
top10Up = head(up.res, 10)
down.res <- all.res[all.res$logFC < 0,]
top10Down = head(down.res, 10)
top10Up
top10Down
write.csv(topTags(et, n=20, sort='p.value'), "d:/RNA-Seq/trimmed-fastq/featureCounts/significant_DGE.csv")
write.csv(top10Up, "d:/RNA-Seq/trimmed-fastq/featureCounts/top10Up_logFC.csv")
write.csv(top10Down, "d:/RNA-Seq/trimmed-fastq/featureCounts/top10Down_logFC.csv")

install.packages("corrplot")
library("corrplot")
corrplot(cor(allCounts), method = 'circle', order = 'hclust', cl.lim=c(0,1))
install.packages("pheatmap")
library("pheatmap")
pheatmap(cor(allCounts), main = "Heatmap of normalised gene counts", display_number = TRUE)

plot(allCounts [,1:2], xlab="CD24Mid-1 gene counts", ylab="CD24Neg-2 gene counts")

MytopGO
#Obtain the gene list
gene_list <- genes(topGO)
head(gene_list)
