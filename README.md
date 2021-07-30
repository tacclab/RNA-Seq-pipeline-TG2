# RNA-Seq-pipeline-TG2
##### Pipeline adapted from [1]

----------------------------
This is a pipeline for RNA-Seq  data we used in order analyze TG2-lncRNA expression

# QUALITY CONTROL, PREPROCESSING, AND ALIGNEMNT
Fastqc, Tophat2, samtools, HTSeq, R (DESeq2) 

## A.	FastQC, Tophat2, samtools, HTSeq
### FastQC tool provides quality control checks metrics on row reads data. It takes all the fastq files coming from RNA-sequencing (in our case 12 fastq were given, 6 for the forward strand, 6 for the reverse strand, for a total of 6 samples, 3 treated and 3 controls).
$ fastqc Fastq/*.fastq -o FastQC
### It removes sequences with low quality. This step should be repeated for each 
sample.
$ awk -v s=10 -v e=0 ‘ {if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; } ‘ Fastq/DLB_S26_R1_001.fastq > Fastq/Trimmed_DLB_S26_R1_001.fastq
$ awk -v s=10 -v e=0 ‘ {if (NR%2 == 0) print substr($0, s+1, length($0)-s-e); else print $0; } ‘ Fastq/DLB_S26_R2_001.fastq > Fastq/Trimmed_DLB_S26_R2_001.fastq

### Tophat tool is used to align row read to a reference genome
$ tophat2 -r 200 -p 8 -o name/Thophat_out -G 
~/Indexes/.../UCSC/.../Genes/genes.gtf
~/Indexes/.../UCSC/.../Sequence/Bowtie2Index/genome

### This step should be repeated for each sample.
Fastq/Trimmed_DLB_S26_R1_001.fastq Fastq/Trimmed_DLB_S26_R2_001.fastq

### Bam files have to be sorted with samtools in order to be mapped on a reference genome:
$ samtools sort DLB_S26_accepted_hits.bam -o /mnt/Chiara4TB/BAMs_sorted/DLB_S26_sorted.bam -O BAM -n -@ 4 &

### HTSeq and samtools are used to quantify the number of reads mapped on each gene and to assemble gene expression:
$ samtools view DLB_S26/Tophat_Out/DLB_S26_accepted_hits.sorted.bam | python -m HTSeq.scripts.count -q -s no - ~/Bowtie2Index/Homo_sapiens_nuovo/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2015-07-17-14-32-32/Genes/genes.gtf > DLB_S26/DLB_S26.count.txt

## B.	R (DESeq2)
### Run DESeq2 in RStudio environment  to do the differential gene expression analysis.
### Load necessary libraries
> library(gplots)
> library(DESeq2)
> library(pheatmap)

### Create necessary data object
> setwd("/Users/Jimi/Desktop/analysis")
> sample.names <-sort(paste(c("MT","WT"),rep(1:3,each=2),sep="")) > file.names <- paste(sample.names,".count.txt",sep="") 
> condizioni <- factor(c(rep("MT",3),rep("WT",3))
> sampleTable <- data.frame (sampleName = sample.names, fileName = file.names, condition = factor(c(rep("MT",3),rep("WT",3))))

### HTSeq count data
> ddsHTSeq <- DESeqDataSetFromHTSeqCount (sampleTable = sampleTable, directory=".", design =~ condition)

### Differential gene analysis
> ddsHTSeq<-ddsHTSeq[rowSums(counts(ddsHTSeq))>10,]
> dds<-DESeq(ddsHTSeq)

### Quality check on the samples
> rld<-rlogTransformation(dds,blind=FALSE)

### Plot
> pdf("plotPCA_real.pdf") [Figure.1]
> plotPCA(rld,intgroup="condition",ntop=nrow(counts(ddsHTSeq)))
> dev.off()

### Plot correlation heatmap
> cU<-cor(as.matrix(assay(rld)))
> cols<-c("dodgerblue3","firebrick3")[condizioni]
> pdf("heatmap2_real.pdf")[Figure.2]
> heatmap.2(cU, symm = TRUE, col = colorRampPalette(c("white","blue")) (100), labCol = colnames(cU), labRow = colnames(cU), distfun = function(c) as.dist(1-c), trace = "none", Colv = TRUE, cexRow = 0.9, cexCol = 0.9, key = F, font = 2, RowSideColors = cols, ColSideColors=cols)
> dev.off()

### Differential gene anlysis results
> res <- results(dds,contrast=c("condition","MT","WT"))
> grp.mean <- sapply(levels(dds$condition),function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition==lvl]))
> norm.counts<-counts(dds,normalized=TRUE)
> all<-data.frame(res,grp.mean,norm.counts)
> write.table(all,file="DESeq2_all_rm.txt",sep="\t")

### Figures on significantly differentiated genes
> pdf("plotMA_real.pdf")[Figure.3]
> plotMA(res,ylim=c(-5,5),alpha=0.01)
> dev.off()
> topGene <- rownames(res) [res$padj <= sort(res$padj)[5]  & !is.na(res$padj)]
>with(res[topGene, ], {
points(baseMean, log2FoldChange, col="dodgerblue", cex=1.5, lwd=2)
text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
> sig.dat<-assay(rld)[res$padj < 0.01 & !is.na(res$padj),]
> annC <- data.frame(condition = condizioni)
> rownames(annC)<-colnames(sig.dat)
> pdf('rplot_real.pdf') [Figure.4]
> pheatmap(sig.dat, scale = "row", fontsize_row = 9, annotation_col = annC)
> dev.off()

### Filtered results based on FoldChange (>=1.2 or >=0.5 & <=0.9)
> sig.dat$fc < ((sig.dat$MT1 + sig.dat$MT2 + sig.dat$MT3)/3)/((sig.dat$WT1 + sig.dat$WT2 + sig.dat$WT3)/3)
> sig.dat2 <- as.data.frame(sig.dat)
> sig.dat2$fc <- ((sig.dat2$MT1 + sig.dat2$MT2 + sig.dat2$MT3)/3)/((sig.dat2$WT1 + sig.dat2$WT2 + sig.dat2$WT3)/3)
> sig.dat3 <- sig.dat2[sig.dat2$fc >= 1.2 | (sig.dat2$fc >= 0.5 & sig.dat2$fc <= 0.9),]

### Remove the Fc column 
> sig.dat4 <- sig.dat3[,c(1:6)]
> pdf("rplot_filtered_real.pdf")[Figure.5]
> pheatmap(sig.dat4, scale="row", fontsize_row = 3, annotation_col = annC)
> dev.off()


# REFERENCE
1. Yalamanchili HK, Wan YW, Liu Z. Data Analysis Pipeline for RNA-seq Experiments: From Differential Expression to Cryptic Splicing. Curr Protoc Bioinformatics. 2017;59:11.15.1-11.15.21. Published 2017 Sep 13. doi:10.1002/cpbi.33


