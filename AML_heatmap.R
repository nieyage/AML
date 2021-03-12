####test####
rm(list = ls())
setwd("D:/project/AML")
library(pheatmap)
counts<-read.table("normalized_counts.txt")
head(counts)
MOLM13<-counts[,1:6]
MV411<-counts[,7:12]
U937<-counts[,13:18]

human_AML<-read.csv("human-AML.csv")
head(human_AML)
AF9<-human_AML$MLL.AF9.MOLM13..target.genes
AF4<-human_AML$MLL.AF4.MV4.11..target.genes
str(AF4)
#######MOLM13 target gene heatmap #############
MOL<-MOLM13[which(rownames(MOLM13)%in%AF9),]
head(MOL)
count=t(scale(t(MOL),scale = T,center = T))
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
#######MV411 target gene heatmap #############
MV<-MV411[which(rownames(MV411)%in%pro),]
head(MV)
count=t(scale(t(MV),scale = T,center = T))
count<-na.omit(count)
pheatmap(count,cluster_cols = T,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         #cellwidth = 10, cellheight = 10,
         show_rownames=F,show_colnames=T)
########cell cycle###########
G1<-c("ANAPC2","CCND1","CCNE1","CDC34","CDK4","CDK6","CDKN1B","CDKN3","CUL1","CUL2","CUL3","SKP2")
S<-c("ABL1","MCM2","MCM3","MCM4","MCM5","PCNA","RPA3","SUMO1","UBA1")
G2<-c("ANAPC2","ANAPC4","DIRAS3","BCCIP","BIRC5","CCNB1","CCNG1","CCNH","CCNT1","CCNT2","CDK5R1","CDK5RAP1","CDK7","CDKN3","CKS1B","CKS2","DDX11","DNM2","GTF2H1","GTSE1","HERC5","KPNA2","MNAT1","SERTAD1")
M<-c("CCNB2","CCNF","CDK1","CDC16","CDC20","MRE11A","RAD51")
checkpoint<-c("ATM","ATR","BRCA1","BRCA2","CCNG2","CDC2","CDC34","CDK2","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN3","CHEK1","CHEK2","CUL1","CUL2","CUL3","GADD45A","HUS1","KNTC1","MAD2L1","MAD2L2","NBN","RAD1","RAD17","RAD9A","RB1","RBBP8","TP53")
ne<-c("ATM","BAX","BRCA1","CDKN2B","RBL1","RBL2","TP53")

proliferation<-c("KRAS","FGF3","ERBB2","FGF4","IL3","CDC20","FLT3","CSF2","JUP","CDK4","CDKN3","CUL3","RPA3","SALL4","PRL-3","KIT","FLT3","PML",
                 "IL9","MYC","E2F1","PDGFB","EGFR","CSF1R","ROS1","KRAS","SRC","YES","FES","ABL1","ZHX2","MOS","CRK","FOS","THRA","REL","MYB","CCND3",
                 "CCNE2","CCNF","CYR61","MECOM")
inhibit_proli<-c("RB1","APC","CDKN2A","NF2","VHL","NTRK1","E2F1","HRAS","MAP2K1","ZBTB16","CEBPA","CDK6","CDKN1B","SOCS1","IDH1 ","WT1","TP53","PML","AREG","CCNG2","CDKN1A","CDKN1C","GADD45A")
pro_apoptosis<-c("CTSB","TNFSF10","CASP10","CASP9","RIPK1","BAD","CASP8","APAF1","TRAF1","ERN1","EIF2AK3","CASP6","BNIP3","BNIP3L","FOXO3A","DAPK2","GADD45G","HINT2","Mien1","P2rx7","Ifng","SIAH1","CASP1","MYC","TP53","FOS","REL","CREB1","WT1")
inhibit_apoptosis<-c("KRAS","BID","IRAK1","MALT1","NAIP","BIRC3","BIRC2","CSF1R","BIRC5","E2F1","ADRB2","SCIN","SIX1","Syce3","TSC22D3","VEGFA","MECOM")
m
MV<-counts[which(rownames(counts)%in% m),]
head(MV)
#count<-MV
count=t(scale(t(MV[,1:12]),scale = T,center = T))
count<-na.omit(count)
library(pheatmap)
pheatmap(count,cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         cellwidth = 10, cellheight = 10,
         show_rownames=T,show_colnames=T)
head(MOLM13)
test<-U937
head(test)
test<-na.omit(test)
test[,7]<-"na"
test<-test[,c(7,1:6)]
write.table(test,"U937_GSEA.txt",sep = "\t")



#########
