rm(list = ls())
setwd("D:/project/AML")
counts<-read.table("rawcounts_ensembl.txt")
head(counts)
library('biomaRt')
library("curl")
library(ggrepel)
library(dplyr)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))#####小鼠：mmusculus_gene_ensembl
my_ensembl_gene_id<-row.names(counts)
options(timeout = 4000000)
#my_ensembl_gene_id<-row.names(raw_count_filt)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),filters = 'ensembl_gene_id', values= my_ensembl_gene_id, mart = mart)
#readcount<-read.csv(file="raw_count_filt.csv",header = TRUE)
head(mms_symbols)
head(counts)
counts$gene<-rownames(counts)

colnames(counts)<-c("ensembl_gene_id","MOLM13-T-1","MOLM13-T-2","MOLM13-T-3","MOLM13-U-1","MOLM13-U-2","MOLM13-U-3",
                    "MV411-T-1","MV411-T-2","MV411-T-3","MV411-U-1","MV411-U-2","MV411-U-3",
                    "U937-T-1","U937-T-2","U937-T-3","U937-U-1","U937-U-2","U937-U-3" )
merge<-merge(counts,mms_symbols,by="ensembl_gene_id")

head(merge)
merge<-merge[,-19:-20]
merge<-merge[!duplicated(merge$external_gene_name), ]
rownames(merge)<-merge$external_gene_name
str(merge)
write.table(merge,"rawcounts_genesymbol.txt")

######find DEG###########
merge<-read.table("rawcounts_genesymbol.txt",header = T)
head(merge)
library(tidyverse)
library(DESeq2)
cell <- factor(c(rep("MOLM13",6),rep("MV411",6),rep("U937",6)))
condition <- factor(c(rep("T",3),rep("U",3),rep("T",3),rep("U",3),rep("T",3),rep("U",3)))
colData <- data.frame(row.names=colnames(merge),condition=condition,cell=cell)
colData
colData$condition <- relevel(colData$condition, ref="U")

countData <- merge[apply(merge, 1, sum) > 1 , ] 

dds<-DESeqDataSetFromMatrix(countData,colData, formula(~cell+condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds,normalized=T) 
head(normalized_counts)
write.table(normalized_counts,"normalized_counts.txt")
dds <- DESeq(dds)
res <- results(dds)
head(res)
resOrdered <- res[order(res$padj),]
resultsNames(dds)

######MOLM13####
head(countData)
MOLM13<-countData[,1:6]
condition <- factor(c(rep("T",3),rep("U",3)))
colData <- data.frame(row.names=colnames(MOLM13),condition=condition)
colData
colData$condition <- relevel(colData$condition, ref="U")
dds<-DESeqDataSetFromMatrix(MOLM13,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds,normalized=T) 
head(normalized_counts)
dds <- DESeq(dds)
MOLM13res <- results(dds)
resOrdered <- MOLM13res[order(MOLM13res$padj),]
resultsNames(dds)
sum(MOLM13res$padj < 0.05, na.rm=TRUE)  #有多少padj小于0.1的 
head(MOLM13res)
diff_gene <-subset(MOLM13res, padj < 0.01 & abs(log2FoldChange)> 1) 
diff<-rownames(diff_gene)



MOLM13res <-subset(MOLM13res, padj < 0.01)
MOLM13res<-MOLM13res[order(MOLM13res$log2FoldChange),]
down<-rownames(MOLM13res[1:25,])

MOLM13res<-MOLM13res[order(MOLM13res$log2FoldChange,decreasing=T),]

up<-rownames(MOLM13res[1:25,])
diff<-c(up,down)
diff
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
count<-count[diff,]
#ECcount<-na.omit(ECcount)
pdf("MOLM13-diff-logFC1-heatmap.pdf",width=8,height=20)
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
down
########upgene############
MOLM13diff_gene <-subset(MOLM13res, padj < 0.01 & log2FoldChange > 2) 
str(MOLM13diff_gene)
MOLM13<-rownames(MOLM13diff_gene)
MOLM13:973genes;
library(clusterProfiler)
library(org.Hs.eg.db)

pdf("/public/home/nieyg/project/AML/GOandKEGG/MOLM13-upgene-GO.pdf")
gene.df <- bitr(MOLM13, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MOLM13-upgene-GO-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MOLM13-upgene-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MOLM13-upgene-GO-CC.csv")
######KEGG##########
pdf("/public/home/nieyg/project/AML/GOandKEGG/MOLM13-upgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MOLM13-KEGG.csv")



########downgene############
MOLM13diff_gene <-subset(MOLM13res, padj < 0.01 & log2FoldChange < -2) 
str(MOLM13diff_gene)
MOLM13<-rownames(MOLM13diff_gene)
MOLM13
library(clusterProfiler)
library(org.Hs.eg.db)

pdf("/public/home/nieyg/project/AML/GOandKEGG/MOLM13-downgene-GO.pdf")
gene.df <- bitr(AF4, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MOLM13-downgene-GO-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)

write.csv(ego,"MOLM13-downgene-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MOLM13-downgene-GO-CC.csv")
######KEGG##########
pdf("/public/home/nieyg/project/AML/GOandKEGG/MOLM13-downgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MOLM13-downgene-KEGG.csv")


####火山图#####
pdf("/public/home/nieyg/project/AML/GOandKEGG/MOLM13-volcano.pdf")

#padj < 0.01 & log2FoldChange > 1
head(MOLM13res)
data<-as.data.frame(MOLM13res)
head(data)
data<-na.omit(data)
data$change = ifelse(data$pvalue < 0.01 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'Up','Down'),
                     'Stable')
#data$change<-as.factor(data$change)

p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-5, 5)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="MOLM13 DEG")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p
data$label=ifelse(data$padj< 0.001 & abs(data$log2FoldChange) >= 3,rownames(data),"")
p+geom_text_repel(data = data, aes(x =data$log2FoldChange, 
                                   y = -log10(data$pvalue), 
                                   label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE)


#######MV411##########
head(countData)
MV411<-countData[,7:12]
condition <- factor(c(rep("T",3),rep("U",3)))
colData <- data.frame(row.names=colnames(MV411),condition=condition)
colData
colData$condition <- relevel(colData$condition, ref="U")
dds<-DESeqDataSetFromMatrix(MV411,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds,normalized=T) 
head(normalized_counts)
dds <- DESeq(dds)
MV411res <- results(dds)
resOrdered <- MV411res[order(MV411res$padj),]
resultsNames(dds)
sum(MV411res$padj < 0.05, na.rm=TRUE)  #有多少padj小于0.1的 
MV411diff_gene <-subset(MV411res, padj < 0.01 & abs(log2FoldChange) >1) 
o<-rownames(MV411diff_gene)
m<-intersect(diff,o)
head(MV411res)
MV411res<-MV411res[order(MV411res$log2FoldChange),]
down<-rownames(MV411res[1:25,])
MV411res<-MV411res[order(MV411res$log2FoldChange,decreasing=T),]
up<-rownames(MV411res[1:25,])
diff<-c(up,down)
MV411<-rownames(MV411diff_gene)
MV411

targetcount<-normalized_counts[which(rownames(normalized_counts)%in%o),]
count=t(scale(t(targetcount),scale = T,center = T))
#ECcount<-na.omit(ECcount)
pdf("MV411-diff-logFC1-heatmap.pdf",width=8,height=20)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
dev.off()
##########upgene###########

pdf("overlap-downgene-GO-BP.pdf")
dev.off()
gene.df <- bitr(m, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.3,
                qvalueCutoff = 0.3,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"overlap-GO-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.3,
                qvalueCutoff = 0.3,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"overlap-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.3,
                qvalueCutoff = 0.3,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"overlap-GO-CC.csv")
######KEGG##########
pdf("overlap-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.3,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.3)
ego
barplot(ego, showCategory=20)
write.table(ego,"overlap-KEGG.csv")
dev.off()
########downgene############
MV411diff_gene <-subset(MV411res, padj < 0.01 & log2FoldChange < -1) 
str(MV411diff_gene)
MV411<-rownames(MV411diff_gene)
MV411
pdf("/public/home/nieyg/project/AML/GOandKEGG/MV411-downgene-GO.pdf")
gene.df <- bitr(AF4, fromType = "SYMBOL",
        toType = c("ENSEMBL", "ENTREZID"),
        OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MV411-downgene-GO-BP.csv")


ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "MF", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MV411-downgene-GO-MF.csv")

ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "CC", ###BP,MF,CC
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)

ego
barplot(ego, showCategory=20)
write.csv(ego,"MV411-downgene-GO-CC.csv")
######KEGG##########
pdf("/public/home/nieyg/project/AML/GOandKEGG/MV411-downgene-KEGG.pdf")

ego <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.05)
ego
barplot(ego, showCategory=20)
write.table(ego,"MV411-downgene-KEGG.csv")

pdf("/public/home/nieyg/project/AML/GOandKEGG/MV411-volcano.pdf")

head(MV411res)
data<-as.data.frame(MV411res)
head(data)
data<-na.omit(data)
data$change = ifelse(data$pvalue < 0.01 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'Up','Down'),
                     'Stable')
p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-5, 5)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="MV411 DEG")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p

#####U937#######
head(countData)
U937<-countData[,13:18]
condition <- factor(c(rep("T",3),rep("U",3)))
colData <- data.frame(row.names=colnames(U937),condition=condition)
colData
colData$condition <- relevel(colData$condition, ref="U")
dds<-DESeqDataSetFromMatrix(U937,colData, formula(~condition)) 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds,normalized=T) 
head(normalized_counts)
dds <- DESeq(dds)
U937res <- results(dds)
resOrdered <- U937res[order(U937res$padj),]
resultsNames(dds)
sum(U937res$padj < 0.05, na.rm=TRUE)  #有多少padj小于0.05的 
U937diff_gene <-subset(U937res, padj < 0.05 & log2FoldChange < -1) 
str(U937diff_gene)
U937<-rownames(U937diff_gene)
U937
data<-as.data.frame(U937res)
head(data)
data<-na.omit(data)
data$change = ifelse(data$pvalue < 0.05 & abs(data$log2FoldChange) >= 1, 
                     ifelse(data$log2FoldChange> 1 ,'Up','Down'),
                     'Stable')
pdf("/public/home/nieyg/project/AML/GOandKEGG/U937-volcano.pdf")

p <- ggplot(data = data, 
            aes(x = data$log2FoldChange, 
                y = -log10(data$pvalue), 
                colour=change,
                label = rownames(data))) +
  geom_point(alpha=1, size=1) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-5, 5)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.01),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="U937 DEG")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p


Moverlap<-intersect(MV411,MOLM13)
Moverlap
last<-setdiff(Moverlap,U937)
last
last<-setdiff(U937,Moverlap)
last

mart<-read.csv("Human_genemart_export.csv")

########heatmap###########
head(normalized_counts)
library(pheatmap)
overlap<-normalized_counts[which(rownames(normalized_counts)%in%Moverlap),]
count=t(scale(t(overlap),scale = T,center = T))
#ECcount<-na.omit(ECcount)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)



MLL<-c("HOXA9","HOXA10","MEIS1","PBX3","MEF2C","BCL2","CDK6","CDK9","NFKB1","MYC")
MLL<-normalized_counts[which(rownames(normalized_counts)%in%MLL),]
MLL
count=t(scale(t(MLL),scale = T,center = T))
#ECcount<-na.omit(ECcount)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)



U937res <-subset(U937res, padj < 0.05)
diff<-rownames(U937res)
targetcount<-normalized_counts[which(rownames(normalized_counts)%in%diff),]
count=t(scale(t(targetcount),scale = T,center = T))
head(count)
count<-count[diff,]
pdf("U937-diff-top25-heatmap.pdf",width=8,height=20)
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)
dev.off()
down



AF9<-c("USP32","NUDT3","LRRC41","RCHY1","RSPH9","FAM20B","PACSIN2","VPS18","UHRF1BP","RNF220","SFT2D3","STT3A","SYNGAP1","MXI1","SMYD3","SAFB2","SRP9","MTCH1","TEX2","ASCC2","TPM3","SOBP","CEP85","PIK3AP1","IL17RA","ZC3H11A","ATP6V1C","AKR1A1","B3GNT2","ARHGEF1","TP73","ATG7","CCDC25","ITGA11","AP1B1","CCDC12","INTS7","CXXC5","COG2","FUCA2","BRD4","MRFAP1","PTPRE","POC1A","OGG1","RPA3","UBTF","FNDC3B","SLC6A6","AK4	","FAF1","MRPL9","FAR2","E2F7","ANKRD18","GGA1","UBE2L3","IKBKE","COX7A2","B3GALNT","INTS9","TEAD3","TXNDC12","TRIB1","HAUS6","BCAN","TARBP1","LSM1","CTBP2","DERL1","SYNGR2","TRIM11","SLC38A2","RBM14-R","TMEM53","PGM2","SKIV2L","PRDX3","CD300A","CS	","FBXO5","NELFB","PSMA1","TAP2","PRCC","PPM1F","SDE2","DESI2","CHSY1","ZNF395","NAT9","CDCA5","UBAC1","MEF2D","TM2D2","GLO1","SERPINB","MLLT1","NEMP2","ARRB1","TACC1","TMEM104","TXLNA","MUTYH","CITED2","USP7","TOP1MT","SRSF3","FOXP2","FHOD1","MACROD1","ST3GAL1","ZNF250","DPH7","CALU","NFE2L1","POMK","RPL36AL","DENND1A","NEK6","TRIM35","NCK2","RARA","AARS2","SLC45A3","SLC22A1","COQ2","IPO9","RUNX1","REEP4","ORC1","FDFT1","QSOX1","ST3GAL4","DCPS","ARMC5","C1orf13","BATF3","KLF10","PMVK","HPS6","TNPO1","FSCN1","MTFR2","FAAP100","ARF1","SLC38A1","SH3GL1","PVT1","SRGN","SMG5","PRKCA","UBAP2","PARP1","MOB3A","SLC20A2","FANCE","PNPLA3","EZR	","GPR68","KCNAB2","SETD7","MYEOV","FOXP4-A","CALM2","HEATR1","DUSP6","GFI1","AGTRAP","RUSC1-A","EHBP1","GRK6","TPRN","VEGFA","PLAGL2","IRF8","RPTOR","SRRT","NOP9","UQCC2","RAD54L","CCDC71","CHST11","ACAT2","SUPT16H","EHD1","QSOX2","PYCR2","PHGDH","FAM89A","UCK2","WDR5","PRKDC","GTF3A","CNIH4","AMD1","MIR155H","SASH1","TRA2B","FLAD1","MRPL3","HRH2","ERI3","TAGLN2","ISG20L2","ZBTB24","NOL6","SCARB1","TOMM5","TFB2M","DHX9","ZDHHC14","ATP1A1","PXT1","DNAAF2","KCNQ5","EIF4EBP","SFXN2","LGALS12","SNHG16","PUF60","AKAP1","UCHL5","SQLE","EXOSC4","EXOSC3","PPIF","NR1D1","ENO1","WDR77","EIF4G1","WT1","PFDN2","KIF21B","XPO5","RCL1","OAF","MAK16","FADS2","SIGMAR1","PLK3","HBEGF","PES1","KANK1","BMP8B","CCL3","MIR17HG","NT5DC2","PDF","SLC16A1","MEX3A","KCNK5","SLC39A1","MYC","SLC29A")
AF4<-c("CD1C","LINC01222","KAT6A","RPL23AP53","AP3M1","AP3D1","NBN","DROSHA","PHIP","ZDHHC20","HMGN1","TMEM242","EFR3A","RAB2A","GEMIN8P4","ASH2L","KLHL26","GIGYF2","ORC2","MRFAP1","CEP85","CYSTM1","RNLS","FNTA","NDFIP2-AS1","TMBIM6","MIR181A1HG","NFATC1","IQGAP1",
       "IFT80","SLC25A44","PAK2","SMG1P3","NFIX","C5orf51","MAPK6","TUBGCP3","SDF2","JOSD1","PURB","VAPA	","CARS2","USP4	","ENTPD4","GDAP2","PRPF39","WWP1","DDX39B","SNHG21","THAP1","CTU1","LINC00861","EYA3","KCNE1","SEPHS2","ZC3H13","ABCC4","LRRC8D","CDCA7","WNK1","SLC1A3","SAFB2","C16orf87","METTL17","COPS4","AP4E1","TTC33","EP400","UIMC1","LRIF1","BAALC-AS2","CA3-AS1",
       "UBTF","TMEM68","SAMD12","INO80","CTHRC1","SETD7","DIDO1","ST13","TPP2","ZNF622","YWHAZ","SH2D4A","RNU5F","CAP1","ESD","LDLRAD","UPF1","STARD3N","PAIP1","ATE1","NCBP2","UBC","AP2S1","RHOBTB","PARN","TRIM35","ECSIT","MIER2","CRADD","IL17RA","MSRA","TBC1D7","DOCK9","TMED5","CDC37","ERG","CBLL1","USP22","PRPF38B","ABHD10","TNFRSF10","DLEU1","TAF2","FLT3LG","GLYR1","NR6A1","PRPF4B","CCNT2-AS","DCXR","C9orf64","SLTM","RRP8","HNRNPUL","CELF2","BAX","SIRT6","PLCB4","COPS5","PROSER1","GNA13","SENP5","RXRB","MFHAS1","TRAPPC1","PET117","C19orf4","PDS5B","TM9SF3","PAM16","RCAN1","UBP1","USP12","PCBP1-AS","KPTN","LINC0100","MCM9","CYP19A1","ZBTB7A","OPTN","PLEKHF1","PPP6R1","PDCD2L","TMEM238","LSM14A","ADIPOR2","CHD1","BTBD1","TMEM126","CAPRIN1","SLC35F2","MED10",
       "IL31RA","SAP18","CAAP1","MTF1","MLLT1","LINC0053","LYST","HMGCS1","TCEA1","DENND4A","FAM27C","AZIN1","ANKRD33","LATS2","MICU2","GTF3A","FBXO33","TMEM249","R3HCC1","GLS","ZC3H4","BRD4","MRPS31","CCT8","PRICKLE","IL1RAP","SRSF4","CTNNB1","ZNF507","ANKH","NUDT5","RANBP3","SLC43A1","FSD2","PAF1","KPNA3","PIGF","TPM4","FAH","FAM157A","MYO10","DDX50","DNAJC15","TTC26","GSR","KLHL8","TOP1","LILRA5","CAB39","TMED8","RUNX1","PTDSS1","PSMA6","LINC0021","ZNF595","GTF2E2","UBALD2","BIN3","CCSAP","CCDC102","TOX","PSMB7","DTNBP1","ANKRD40","AKAP8","DDX11L1","TACC1","BCL2","EMC2","FLT1","SUB1","MIR155H","BCLAF1","GTF2A1","UBE2D3","ICE1","ARF1","POLG","QSER1","RSL24D1","CLIC4","TM2D2","ARPP19","KLLN","HEATR3","TMCO1","RBM14-RB","FAM207A","MORF4L1","THUMPD3","HIVEP3","ZC3H8","MTUS2","RAB29","EIF3M","SPG21","CNOT7","RFXAP","PPP2CB","NR3C1","USP7","CHAMP1","PNRC2","CBX5","CCDC25","KBTBD6","ANKRD27","NSUN2","CLNS1A","CUL4A","SLC2A1","DHX29","DNAJC7","DCSTAMP","KLHL2","SLC30A5","WDR75","HRH2","HNRNPLL","PRC1","PPP6R3","PPP1CC","NOP16","INTS8","ATP6V1C","PSMG3","KCTD3","KBTBD2","SLC38A2","MRPS28","ADAM9","SGTA","WRNIP1","NUP54","TGDS","PSPC1","TNPO1","MTMR6","CHD7","ZBTB2","XPO4","URI1","CPSF6","PYGM","ATG3","SCLT1","MYBL1","PPP2R2A","TMEM209","PLAUR","ERGIC2","CDC27","VPS35","OTULIN","CHERP","SETMAR","DDX23","CRLS1","SLC25A1","SNX5","NFKBID","FRG1","ACSL1","NLRP12","RAP2A","DNAJC8","ATF3","HOMER1","ANKLE2","EBAG9","NEMP2","STARD13","REEP4","SFPQ","CCDC168","XPO6","EXOSC5","SH3RF3","RHEB","NMD3","G3BP2","RPF1","CAGE1","CCDC117","SUCLA2","CEBPG","MRPL36","BRCA2","SUGT1","TSPAN3","SMARCA4","PINX1","SH3GL1","MTDH","CNIH4","POLR3A","PCID2","STRAP","NDUFAF6","SACS","XRCC5","PPP1R8","PSIP1","CENPJ","VDAC3","SOD1","LDLR","YBX1","SLC38A1","ALG5","FASTKD3","DLEU2","KIFC1","SAR1A","MRPL15","EIF4EBP","SYCE2","NT5C3A","U2AF2","GATAD2A","ANP32B","NADK2","CHRAC1","MCUR1","EIF4G2","PPAT","PDE12","TRMT11","TMEM64","ARMC1","TUBB1","TMX1","SLC36A4","PPP1R14","ZNF706","NCAPH2","PPIF","KHSRP","C1orf10","ENO1","KIR3DX1","HNRNPD","CEBPZ","LONRF1","KIF21B","NUCKS1","ACAT1","UBE2V2","LPCAT1","SCRIB","ZC3H15","OXCT1","DPM1","PKM","SHQ1","ANGPT1","TXNDC9","SNTB1","SNHG16","NUP37","TSEN2","RAD54B","LAP3","LONP1","PIH1D2","SLC16A1","IPO5","EGR1","MARS2","TCTE3","DUT","LETM1","SLC7A1","SUPT16H","PRPF19","TRMT10C","MPO","TEX10","CBWD1","RGCC","PSMA3","RFTN1","CASC11","MTHFD2","PSMD1","PRKDC","TFB2M","SRSF7","SLC25A1","COQ2","SLC1A4","DDX21","DGCR8","NUDT15","NAA50","ELAVL1","SLC3A2","GSPT1","PLK3","SAE1","DNMT1","RAD21","NEXN","TRA2B","IL7R","DHX9","CMC2","ASS1","RBL1","TFAM","NUFIP1","TGFBR1","AMD1","CDCA7L","KIF14","TEX30","RRP12","POMP","HNRNPA3","ATP1A1","NLE1","PTGES3","DHFR","EIF2S2","NUDCD1","CDC42EP","AP3M2","NOP58","FARSB","RAD1","HSPA9","IVNS1AB","CBWD5","ERI1","SNHG4","LZTS1","BUB1B","WDR43","KCTD15","WT1","USP31","LMNB1","ANKRD18D","UCHL5","HNRNPAB","WDR77","TFDP1","DUSP6","TUBA1B","HAUS8","FABP5","FBXO5","CDC25C","PNO1","HSPD1","LDHA","STMN1","CKS2","FAM72B","HSPH1","NDC1","TRIB3","PAICS","WASF1","OSM","FAM72D","MIR17HG","DIAPH3","UHRF1","HSPA8","ESPL1","OCSTAMP","PBK")
m<-intersect(AF9,AF4)
m
