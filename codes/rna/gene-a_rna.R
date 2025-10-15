#write by ~~~ at 2023.02.13
#01、搭环境"myenrichplot",
{
  for (i in c("data.table","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/BED/special/rRNA.bed"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#03、跑命令
{
  ##去接头、质检
  {
    shell <- paste0(path,"src/run_a.sh")
    cat("#!/bin/bash\n", file = shell)
    prex  <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*1.fq.gz$*"),".1.fq.gz"), "[",1)
    lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
    lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
    cmd_01 <- paste0("trim_galore --phred33 --fastqc --illumina --paired ",lay1," ",lay2," -o ",paste0(path,"trim\n"))
    cmd_02 <- paste0("fastqc --extract -t 30 -o ",path,"QC ",c(paste0(path,"trim/",prex,".R1_val_1.fq.gz\n"),paste0(path,"trim/",prex,".R2_val_2.fq.gz\n")))
    for (i in cmd_01) {cat(i, file = shell, append = T)}
    for (i in cmd_02) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
  }
  ##回帖基因组
  {
    shell <- paste0(path,"src/run_b.sh")
    cat("#!/bin/bash\n", file = shell)
    prex  <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*1.fq.gz$*"),".1.fq.gz"), "[",1)[c(1,2,5,6)]
    lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
    lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
    cmd_01 <- paste0("#trim_galore --phred33 --fastqc --illumina ",
                     "--paired ",lay1," ",lay2," -o ",path,"trim","\n")
    lay1 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
    lay2 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
    cmd_02 <- paste0("STAR --runThreadN 30 --readFilesCommand zcat --outSAMtype ",
                     "BAM SortedByCoordinate --outMultimapperOrder Random --genomeDir ",inx1," ",
                     "--outFileNamePrefix ",path,"bam/",prex," ",
                     "--outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 ",
                     "--readFilesIn ",lay1," ",lay2,"\n")
    cmd_03 <- paste0("mv ",path,"bam/",prex,"Aligned.sortedByCoord.out.bam ",path,
                     "bam/",prex,".bam","\n")
    cmd_04 <- paste0("samtools index ",path,"bam/",prex,".bam","\n")
    cmd_05 <- paste0("TEtranscripts -t ",
                     paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                     paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                     " --GTF ",inx3," --TE ",inx4," --sortByPos --mode uniq ",
                     "--project uniqpr --outdir ",path,"count/","\n")
    for (i in cmd_01) {cat(i, file = shell, append = T)}
    for (i in cmd_02) {cat(i, file = shell, append = T)}
    for (i in cmd_03) {cat(i, file = shell, append = T)}
    for (i in cmd_04) {cat(i, file = shell, append = T)}
    for (i in cmd_05) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
  }
  ##bamCoverage
  {
    shell <- paste0(path,"src/run_bw.sh")
    cat("#!/bin/bash\n", file = shell)
    prex  <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*R1.fq.gz$*"),".R1.fq.gz"), "[",1)
    lay1 <- paste0(path,"fq/",prex,".R1.fq.gz")
    lay2 <- paste0(path,"fq/",prex,".R2.fq.gz")
    cmd_31 <- paste0("bamCoverage  --normalizeUsing CPM --numberOfProcessors 20 --minMappingQuality 1 --samFlagExclude 256 --outFileFormat bigwig --binSize 20 -b ",
                     paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
    for (i in cmd_31) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
  }
  ##统计reads
  {
    shell <- paste0(path,"src/run_c.sh")
    cat("#!/bin/bash\n", file = shell)
    prex <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*R1.fq.gz$*"),".R1.fq.gz"), "[",1)
    cmd_31 <- paste0("#TEtranscripts -t ",
                     paste(paste0(path,"bam/",prex[3:6],".bam"),collapse = " ")," -c ",
                     paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " "),
                     " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                     "--project ccpr --outdir ",paste0(path,"count/TEtranscripts_geneName"),"\n")
    cmd_32 <- paste0("#TEcount -b ",paste0(path,"bam/",prex,".bam")," --GTF ",inx3,
                     " --TE ",inx4," --sortByPos --mode multi --outdir ",
                     paste0(path,"count/TEcount")," --project ",
                     sapply(str_split(prex,pattern = "_"), "[",1),"\n")
    cmd_33 <- paste0("#featureCounts -g 'transcript_id' -F 'GTF' -T 30 -a ", inx4,
                     " -o ",paste0(path,"count/featureCounts_sites/","ccpr.sites.count.txt "),
                     paste0(path,"bam/",prex,".bam", collapse = " "),"\n")
    cmd_34 <- paste0("#TElocal -b ",
                     paste0(path,"bam/",prex,".bam"),
                     " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                     paste0(path, "count/TElocal/","TEsites.",prex),"\n")
    cmd_35 <- paste0("paste ",paste0(path, "count/TElocal/* |"),
                     "awk -F '\\t' -v p='\\t' '{print $1p$2p$4p$6p$8p$10p$12}' |grep -v ENSMU* ",
                     ">",paste0(path, "count/TElocal/pr.TE.sites.txt"),"\n")
    for (i in cmd_31) {cat(i, file = shell, append = T)}
    for (i in cmd_32) {cat(i, file = shell, append = T)}
    for (i in cmd_33) {cat(i, file = shell, append = T)}
    for (i in cmd_34) {cat(i, file = shell, append = T)}
    for (i in cmd_35) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
  }
  ##组装转录组 & 统计reads
  {
    ##stringTie version 1.3.5
    shell <- paste0(path,"src/run_stringTie.sh")
    cat("#!/bin/bash\n", file = shell)
    prex  <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*R1.fq.gz$*"),".R1.fq.gz"), "[",1)
    cmd_31 <- paste0("#stringtie ", paste0(path,"bam/",prex,".bam")," ",
                     "-p 30 -G ", inx3, " -o ", paste0(path,"stringTie/",prex,".gtf"),"\n")
    cmd_32 <- paste0("#stringtie --merge -o ",paste0(path,"stringTie/v_KO74.merge.gtf")," -G ",inx3," ",
                     paste(paste0(path,"stringTie/",prex[1:4],".gtf"), collapse = " "),"\n")
    cmd_33 <- paste0("#stringtie --merge -o ",paste0(path,"stringTie/v_KO92.merge.gtf")," -G ",inx3," ",
                     paste(paste0(path,"stringTie/",prex[c(1,2,5,6)],".gtf"), collapse = " "),"\n")
    cmd_34 <- paste0("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/TEtranscripts -t ",
                     c(paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                       paste(paste0(path,"bam/",prex[5:6],".bam"),collapse = " ")),
                     " -c ", paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " "),
                     " --GTF ",paste0(path,"stringTie/v_KO",c(74,92),".merge.gtf"),
                     " --TE ",inx4," --sortByPos --mode multi ",
                     "--project v_",c(74,92),"merge"," --outdir ", paste0(path,"stringTie/count"),"\n")
    for (i in cmd_31) {cat(i, file = shell, append = T)}
    for (i in cmd_32) {cat(i, file = shell, append = T)}
    for (i in cmd_33) {cat(i, file = shell, append = T)}
    for (i in cmd_34) {cat(i, file = shell, append = T)}
    print(paste0("nohup bash ",shell, " >>",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
  }
}
#04、差异Gene、TE、TE sites from DESeq2
{
  #reads from TEtranscripts
  data <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/TEtranscripts_geneName/ccpr.cntTable",
                     header = T, check.names = F, stringsAsFactors = F)
  gene <- data[grep(data[,1], pattern = ":", invert = T),]
  tete <- data[grep(data[,1], pattern = ":", invert = F),]
  rc_a <- gene[,c(1,4,2,7,6)]
  rc_b <- gene[,c(1,5,3,7,6)]
  rc_c <- tete[,c(1,4,2,7,6)]
  rc_d <- tete[,c(1,5,3,7,6)]
  #reads from featureCounts
  data <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/featureCounts_sites/ccpr.sites.count.txt",
                     header = T, check.names = F, stringsAsFactors = F, skip = 1)
  data <- data[,-c(2,3,4,5,6)]
  rc_e <- data[,c(1,4,5,2,3)]
  rc_f <- data[,c(1,6,7,2,3)]
  #reads from featureCounts
  data <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/TElocal/pr.TE.sites.txt",
                     header = T, check.names = F, stringsAsFactors = F)
  rc_e <- data[,c(1,4,5,2,3)]
  rc_f <- data[,c(1,6,7,2,3)]
  #自定义的差异基因分析的函数
  shy_diff <- function(mydata, Class){
    counts <- mydata
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    namess <- sapply(str_split(basename(colnames(counts)),pattern = ".bam"), "[", 1)
    colnames(counts) <- namess
    ##DESeq2
    yyyCOL <- data.frame(row.names = colnames(counts), Class)
    ddsdds <- DESeqDataSetFromMatrix(countData = counts, colData = yyyCOL, design = ~Class)
    ddsdds <- ddsdds[rowMeans(counts(ddsdds)) > 5,]
    ddsdds <- DESeq(ddsdds)
    resres <<- results(ddsdds, contrast = c("Class","KO","WT"))
    #return(resres)
  }
  #临时：
  shy_diff(mydata <- data[,c(1,4,2,7,6)],Class = factor(c("KO","KO","WT","WT")))
  tmp1 <- as.data.frame(resres)
  tmp1[grep(x = rownames(tmp1), pattern = "RLTR1B"),]
  ##NO.1
  shy_diff(mydata = rc_a, Class = factor(c("KO","KO","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","pr.KO74.gene.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","pr.KO74.gene.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","pr.KO74.gene.down.csv"))
  ##NO.2
  shy_diff(mydata = rc_b, Class = factor(c("KO","KO","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","pr.KO92.gene.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","pr.KO92.gene.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","pr.KO92.gene.down.csv"))
  ##NO.3
  shy_diff(mydata = rc_c, Class = factor(c("KO","KO","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","pr.KO74.fullte.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","pr.KO74.fullte.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","pr.KO74.fullte.down.csv"))
  ##NO.4
  shy_diff(mydata = rc_d, Class = factor(c("KO","KO","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","pr.KO92.fullte.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","pr.KO92.fullte.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","pr.KO92.fullte.down.csv"))
  ##NO.5
  shy_diff(mydata = rc_e, Class = factor(c("KO","KO","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","pr.KO74.teSites.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","pr.KO74.teSites.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","pr.KO74.teSites.down.csv"))
  ##NO.6
  shy_diff(mydata = rc_f, Class = factor(c("KO","KO","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","pr.KO92.teSites.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","pr.KO92.teSites.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","pr.KO92.teSites.down.csv"))
}
#05、画图：火山图 for Genes
{
  shy_volcano_gene <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold[resdata$log2FoldChange > log2(1.5) &  resdata$padj <0.05] = "Upregulated genes"
    resdata$threshold[resdata$log2FoldChange < -log2(1.5) & resdata$padj <0.05] = "Downregulated genes"
    resdata$threshold[!(resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05) & !(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05)] = "Other genes"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated genes","Downregulated genes","Other genes"))
    p_1 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","gray"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.2,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    ggsave(plot = p_1, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_gene(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.volcano.svg")
  shy_volcano_gene(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.volcano.pdf")
}
#06、画图：火山图 for TEs
{
  shy_volcano_te <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold[resdata$log2FoldChange > log2(1.5) &  resdata$padj <0.05] = "Upregulated TEs"
    resdata$threshold[resdata$log2FoldChange < -log2(1.5) & resdata$padj <0.05] = "Downregulated TEs"
    resdata$threshold[!(resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05) & !(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05)] = "Other TEs"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated TEs","Downregulated TEs","Other TEs"))
    resdata[which(resdata$threshold=="Other TEs"),"X"] <- NA
    labb <- resdata
    if (unname(table(resdata$threshold)[1]) + unname(table(resdata$threshold)[2]) >20){
      labb <- labb[which(labb$padj<10^(-50)),]
    }
    p_2 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
           labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
           geom_point(size=0.8) +
           geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
           geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
           geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +#ylim(0,max(-log10(resdata$padj)[is.finite(-log10(resdata$padj))])+20) +
           theme_bw() + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
           scale_color_manual(values = c("red","blue","black"),
                              labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                         paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                         paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
           theme(legend.position = c(0.8,0.8), 
                 legend.background = element_rect(fill = NA),
                 legend.spacing.x = unit(0,"mm"),
                 legend.spacing.y = unit(0,"mm"),
                 legend.key = element_blank(),
                 legend.text = element_text(size = 13),
                 axis.title = element_text(size = 13)) +
           geom_text_repel(data = labb,
                           aes(x=log2FoldChange, y=-log10(padj), label=X), hjust=-0.125, size=4, colour="black", family="serif") +
           guides(color = guide_legend(override.aes = list(size = 3)))
    p_2
    ggsave(plot = p_2, saveFile, units = "cm", width = 18, height = 16)
}
shy_volcano_te(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.fullte.csv",
               saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/pr.KO74.fullte.volcano.pdf")
shy_volcano_te(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.fullte.csv",
               saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/pr.KO92.fullte.volcano.pdf")
}
#07、画图：富集分析
list.files(paste0(path,"DESeq2"), pattern = "*gene.csv$", full.names = T)
##条形图
{
  shy_fuji_bar <- function(resFile, change, savePath){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    if (change=="upup"){gexx <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]}
    if (change=="down"){gexx <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]}
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(gexx), keytype="SYMBOL", column="ENTREZID"))
    PATH_a <- enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kegg <<- PATH_a@result[1:15,]
    PATH_b <- ReactomePA::enrichPathway(gene = geyy, organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    reac <<- PATH_b@result[1:15,]
    GO   <- enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    term <<- GO@result[1:15,]
    ##画图
    kegg$Description <- factor(kegg$Description, levels = rev(rbind(kegg$Description)))
    reac$Description <- factor(reac$Description, levels = rev(rbind(reac$Description)))
    term$Description <- factor(term$Description, levels = rev(rbind(term$Description)))
    p_1 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
      geom_bar(stat = "identity", fill="royalblue", width = 0.8) +
      theme_classic() +
      ylab(NULL) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    p_2 <- ggplot(reac, aes(x=-log10(pvalue), y=Description))+
      geom_bar(stat = "identity", fill="royalblue", width = 0.8) +
      theme_classic() +
      ylab(NULL) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    p_3 <- ggplot(term, aes(x=-log10(pvalue), y=Description))+
      geom_bar(stat = "identity", fill="royalblue", width = 0.8) +
      theme_classic() +
      ylab(NULL) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    ggsave(plot = p_1, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".KEGG.",str_split(basename(resFile),".csv")[[1]][1],"r.pdf"))
    ggsave(plot = p_2, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".ReactomePA.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
    ggsave(plot = p_3, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
  }
  shy_fuji_bar(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv",
               savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/",
               change  = "upup")
  shy_fuji_bar(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv",
               savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/",
               change  = "down")
}
##气泡图
{
  shy_fuji_bubble <- function(resFile, change, savePath){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    if (change=="upup"){gexx <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]}
    if (change=="down"){gexx <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]}
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(gexx), keytype="SYMBOL", column="ENTREZID"))
    PATH_a <- enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kegg <<- PATH_a@result[1:15,]
    PATH_b <- ReactomePA::enrichPathway(gene = geyy, organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    reac <<- PATH_b@result[1:15,]
    GO   <- enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    term <<- GO@result[1:15,]
    ##画图
    kegg$Description <- factor(kegg$Description, levels = rev(rbind(kegg$Description)))
    reac$Description <- factor(reac$Description, levels = rev(rbind(reac$Description)))
    term$Description <- factor(term$Description, levels = rev(rbind(term$Description)))
    p_1 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      ylab(NULL) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    p_2 <- ggplot(reac, aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      ylab(NULL) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    p_3 <- ggplot(term, aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      ylab(NULL) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    ggsave(plot = p_1, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".KEGG.",str_split(basename(resFile),".csv")[[1]][1],".bubble.pdf"))
    ggsave(plot = p_2, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".ReactomePA.",str_split(basename(resFile),".csv")[[1]][1],".bubble.pdf"))
    ggsave(plot = p_3, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".bubble.pdf"))
  }
  shy_fuji_bubble(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv",
                  savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/",
                  change  = "upup")
  shy_fuji_bubble(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv",
                  savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/",
                  change  = "down")
}
##Temp
{
  shy_fuji_temp <- function(resFile, savePath){
    kegg <- as.data.frame(read.table(resFile, header = T, stringsAsFactors = F, sep = "\t"))
    ##画图
    kegg$Description <- factor(kegg$Description, levels = rev(rbind(kegg$Description)))
    p_3 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      ylab(NULL) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    ggsave(plot = p_3, 
           units = "cm", width = 20, height = 12,
           filename = paste0(savePath, str_split(basename(resFile),".txt")[[1]][1],".bubble.pdf"))
  }
  shy_fuji_temp(resFile="/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/Enrichment/MM_pull_downKEGG.txt",
                savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/Enrichment/")
}
###排版图表时的修改 [KEGG]
{
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/"
  resFile <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv"
  #resFile <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/shangshandalaohu.geneOnly.csv"
  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  for (symb in seq(2)) {
    ovgene <- unlist(list(mtacup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    kegg <- kegg@result
    kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
    assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.01),])
  }
  #
  kta1$Description <- factor(kta1$Description, levels = rev(kta1$Description))
  p_05 <- ggplot(kta1[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(8)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
    labs(y="",title = "this is title--upup") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  p_05
  #
  kta2$Description <- factor(kta2$Description, levels = rev(kta2$Description))
  p_06 <- ggplot(kta2[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(10)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    labs(y="",title = "this is title--down") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
          axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
  p_06
  plot <- p_05 +p_06 +plot_layout(nrow = 2)
  ggsave(plot = plot, 
         units = "cm", width = 22, height = 15,
         filename = paste0(path,"PLOT/Enrichment/","pr.kegg.pdf"))
}
##排版图表时的修改  [GOBP、GOMF]
{
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/"
  resFile <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv"
  #resFile <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/shangshandalaohu.geneOnly.csv"
  resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
  resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
  mtacup <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]
  kddown <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]
  for (symb in seq(2)) {
    ovgene <- unlist(list(mtacup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", 
                                      ont = "MF",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
    gogo <- gogo@result
    assign(paste0("gta",symb), gogo[which(gogo$pvalue<0.01),])
  }
  #
  gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
  p_05 <- ggplot(gta1[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(13)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 12), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
    labs(y="",title = "this is title--upup") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  p_05
  #
  gta2$Description <- factor(gta2$Description, levels = rev(gta2$Description))
  p_06 <- ggplot(gta2[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(10)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
    labs(y="",title = "this is title--down") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
          axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
  p_06
  plot <- p_05 +p_06 +plot_layout(nrow = 2)
  ggsave(plot = plot, 
         units = "cm", width = 22, height = 15,
         filename = paste0(path,"PLOT/Enrichment/","pr.gomf.pdf"))
}
#08、画图：表达量热图
{
shy_pheatmap <- function(readsFile, gene_resFile, tete_resFile, savePath, mark){
  print("reads counted by TEtranscripts")
  data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = "_"), "[", 1)
  gene <- data[grep(rownames(data), pattern = ":", invert = T),c(3,1,4,2,6,5)]
  gene <- gene[rowSums(gene) > 5,]
  gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
  tete <- data[grep(rownames(data), pattern = ":", invert = F),c(3,1,4,2,6,5)]
  tete <- tete[rowSums(tete) > 5,]
  tete <- as.data.frame(apply(tete, 2, function(x){x/sum(x)*1000000}))
  ##读取DESeq2结果
  resG <- as.data.frame(read.csv(gene_resFile, header = T, stringsAsFactors = F))
  resG <- na.omit(resG)
  sigG <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
  sigG_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
  sigG_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
  resG <- as.data.frame(read.csv(tete_resFile, header = T, stringsAsFactors = F))
  resG <- na.omit(resG)
  sigT <- resG[which(abs(resG$log2FoldChange)>log2(1.5) & resG$padj<0.05), "X"]
  sigT_upup <- resG[which(resG$log2FoldChange>log2(1.5) & resG$padj<0.05), "X"]
  sigT_down <- resG[which(resG$log2FoldChange < -log2(1.5) & resG$padj<0.05), "X"]
  ##开始画图
  pdf(paste0(savePath, str_split(basename(readsFile), pattern = "\\.")[[1]][1],".by_",mark, "_pheatmap.pdf"), 
      width = 8, height = 10)
  ##gene表达量热图
  pheatmap(mat = gene,
           main="the heatmap of all Genes",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(gene), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##siginificant gene表达量热图
  pheatmap(mat = gene[which(rownames(gene) %in% sigG),],
           main="the heatmap of significant Genes",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(gene), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##Upregulated gene表达量热图
  pheatmap(mat = gene[which(rownames(gene) %in% sigG_upup),],
           main="the heatmap of Upregulated Genes",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(gene), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##Downregulated gene表达量热图
  pheatmap(mat = gene[which(rownames(gene) %in% sigG_down),],
           main="the heatmap of Downregulated Genes",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(gene), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##tete表达量热图
  pheatmap(mat = tete,
           main="the heatmap of all TEs",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(tete), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##siginificant tete表达量热图
  pheatmap(mat = tete[which(rownames(tete) %in% sigT),],
           main="the heatmap of significant TEs",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(tete), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##Upregulated tete表达量热图
  pheatmap(mat = tete[which(rownames(tete) %in% sigT_upup),],
           main="the heatmap of Upregulated TEs",
           scale = "row", cellwidth = 60, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(tete), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  ##Downregulated tete表达量热图
  pheatmap(mat = tete[which(rownames(tete) %in% sigT_down),],
           main="the heatmap of Downregulated TEs",
           scale = "row", cellwidth = 60, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(tete), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  
  dev.off()
}
shy_pheatmap(readsFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/TEtranscripts_geneName/ccpr.cntTable",
             gene_resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv",
             tete_resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.fullte.csv",
             savePath  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/",
             mark = "74")
shy_pheatmap(readsFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/TEtranscripts_geneName/ccpr.cntTable",
             gene_resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv",
             tete_resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.fullte.csv",
             savePath  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/",
             mark = "92")

}
#09、画图：质谱数据 [DIY]
{
  ##读质谱数据（X轴）
  trea <- read.delim("/disk5/aaSHY/1shyJOB/4protein/data/LUask/treat/480_FPEP210002742-1A_Proteins-pr-flag.txt",
                     sep = "\t", header = T, stringsAsFactors = F)[,c(1,2,7)]
  trea <- tidyr::separate_rows(data = trea, "gene_name", sep = ";")
  ctrl <- read.delim("/disk5/aaSHY/1shyJOB/4protein/data/LUask/ctrl/480_FPEP210002741-1A_Proteins-OV-flag.txt",
                     sep = "\t", header = T, stringsAsFactors = F)[,c(1,2,7)]
  ctrl <- tidyr::separate_rows(data = ctrl, "gene_name", sep = ";")
  group_t <- trea[which(trea$gene_name %in% intersect(trea$gene_name, ctrl$gene_name)),]
  group_c <- ctrl[which(ctrl$gene_name %in% intersect(trea$gene_name, ctrl$gene_name)),]
  group_t <- stats::aggregate(group_t$X..Peptides, by=list(type=group_t$gene_name), sum)
  group_c <- stats::aggregate(group_c$X..Peptides, by=list(type=group_c$gene_name), sum)
  group_m <- base::merge(x = group_t, y = group_c, by = "type")
  treas <- trea %>% dplyr::group_by(gene_name) %>% summarise(X..Peptides = mean(X..Peptides))
  onlyTre <- dplyr::setdiff(treas$gene_name, group_m$type)
  temp <- data.frame(type = onlyTre, 
                     x.x = as.integer(treas[which(treas$gene_name %in% onlyTre),]$X..Peptides),
                     x.y = 0)
  group_m <- rbind(group_m, temp)
  group_m$ratio <- log2(group_m$x.x/group_m$x.y)
  ##读表达量数据（Y轴）
  data <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/TEtranscripts_geneName/ccpr.cntTable",
                     header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = "_"), "[", 1)
  gene <- data[grep(rownames(data), pattern = ":", invert = T),c(6,5)]#取WT
  gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
  gene$expr <- log10((gene[,1]+gene[,2])/2+1)
  gene$type <- rownames(gene)
  mydata <- base::merge(x=group_m, y=gene, by="type")
  mydata$diff[(mydata$x.x/mydata$x.y) >=2] = "significant"
  mydata$diff[(mydata$x.x/mydata$x.y) < 2] = "un-significant"
  mydata[1:2,]
  ##读差异分析的数据（红点的显著基因）
  #diff_g <- read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv")
  #significant <- diff_g[which(diff_g$log2FoldChange > log2(1.5) & diff_g$padj < 0.05), "X"]
  #mydata$diff[mydata$type %in% significant] = "Up_regulated"
  #mydata$diff[!(mydata$type %in% significant)] = "Non_significant"
  #mydata[which(mydata$type=="ff"),]
  ##开始画图
  okdata <- mydata[,-c(2,3,5,6)]
  a1 <- setdiff(okdata[is.infinite(okdata$ratio),"type"], "pprr8")
  okdata <- okdata[which(!(okdata$type %in% a1)),]
  okdata$ratio[is.infinite(okdata$ratio)] <- 5
  mark = c("pprr8","ff")
  p_1 <- ggplot(data = okdata) +
    geom_point(aes(x=ratio, y=expr, color=diff)) +
    geom_text_repel(data = okdata[which(okdata$type %in% mark),], aes(x=ratio,y=expr, label=type)) +
    xlab(label = "TreatGroup / ControlGroup") +
    ylab(label = "Log10 (CPM of WT group)") +
    scale_color_manual(values = c("indianred","gray")) +
    scale_x_continuous(breaks = c(-4,-3,-2,-1,0,1,2,3,4,5),
                       labels = c(-4,-3,-2,-1,0,1,2,3,4,"Inf (x / 0)")) +
    theme_classic() +
    labs(title = "TITLE: plot of Protein MS ... ...") +
    theme(legend.text  = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          legend.background = element_rect(fill = "gray95"),
          legend.position = c(0.85,0.95),
          axis.text    = element_text(family = "serif", size = 12),
          axis.title   = element_text(family = "serif", size = 12))
  p_1
  ggsave(plot = p_1, filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/ProteinMS/MS.Inf.v2.pdf", 
         units = "cm", width = 16, height = 12)
  ##筛出目标数据看一下
  okdata[which(okdata$ratio > 0 & okdata$ratio < 2 & okdata$expr>3),"type"]
  okdata[which(okdata$ratio > 0 & okdata$ratio < 2 & okdata$expr > 1 & okdata$expr < 2),"type"]
  okdata[which(okdata$diff=="Up_regulated" & okdata$expr > 1 & okdata$expr < 2 & okdata$ratio > -1),"type"]
  write.csv(x = okdata, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/Protein/peptides.vs.74.csv")
  ##Unique展示：只展示treat组
  unia <- trea[which(!(trea$Accession %in% ctrl$Accession)),]
  unib <- stats::aggregate(unia$X..Peptides, by=list(type=unia$gene_name), sum)
  unic <- merge(x=unib, y=gene, by="type")
  unic$diff[unic$type %in% significant] = "Up_regulated"
  unic$diff[!(unic$type %in% significant)] = "Non_significant"
  colnames(unic)[2] <- "num_peptides"
  ggplot(data = unic) +
    geom_point(aes(x=log2(num_peptides), y=expr, color=diff)) +
    geom_text_repel(data = unic[which(unic$expr < 2 & unic$expr > 1 & log2(unic$num_peptides)>0),], aes(x=log2(num_peptides),y=expr, label=type))
}
#10、画图：富集分析 [GSEA]
{
shy_gsea <- function(resFile, savePath){
  aa <- read.csv(file = resFile,
                 header = T, stringsAsFactors = F)[,c(1,3)]
  ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = T))
  ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ae <- ad$rank; names(ae) <- ad$ENTREZID
  rm(aa,ab,ac,ad)
  af_gogo <- suppressWarnings(clusterProfiler::gseGO(geneList = ae, ont = "BP", seed = T, maxGSSize = 1000, verbose = T, keyType = "ENTREZID", pvalueCutoff = 1, OrgDb = "org.Mm.eg.db", by = "fgsea"))
  testxxx <- af_gogo@result
  testxxx[,10] <- apply(af_gogo@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
  write.csv(testxxx, row.names = F, quote = F,
            file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".result.GOBP.csv"))
  af_kegg <- suppressWarnings(clusterProfiler::gseKEGG(geneList = ae, maxGSSize = 1000, organism = "mmu", pvalueCutoff = 1, use_internal_data = T, seed = T, by = "fgsea"))
  testxxx <- af_kegg@result
  testxxx[,10] <- apply(af_kegg@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
  write.csv(testxxx, row.names = F, quote = F,
            file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".result.KEGG.csv"))
  ##2-cell Genes
  aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt")
  aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
  ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
  af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab))
  rm(aa,ab)
  #af_reac <- suppressWarnings(ReactomePA::gsePathway(geneList = ae, maxGSSize = 1000, eps = NA, organism = "mouse", pvalueCutoff = 1, verbose = T, seed = T, by = "fgsea"))
  ##开始画图
  #dev.new()
  p_a <- myenrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for 2-cell Genes", pvalue_table = T)
  p_1 <- myenrichplot::gseaplot2(x = af_gogo, 1:5, title = "GO_BP_GSEA: top 5", rel_heights = c(1.5, 0.3, 0.6), color = ggsci::pal_lancet()(5), pvalue_table = T)
  p_2 <- myenrichplot::ridgeplot(x = af_gogo, showCategory = 10, fill = "pvalue", label_format = 100)
  p_3 <- myenrichplot::cnetplot(setReadable(af_gogo, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID"), foldChange = ae)
  p_4 <- myenrichplot::gseaplot2(x = af_kegg, 1:5, title = "KEGG__GSEA: top 5", rel_heights = c(1.5, 0.3, 0.6), color = ggsci::pal_lancet()(5), pvalue_table = T)
  p_5 <- myenrichplot::ridgeplot(x = af_kegg, showCategory = 10, fill = "pvalue", label_format = 100)
  p_6 <- myenrichplot::cnetplot(setReadable(af_kegg, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID"), foldChange = ae)
  pdf(file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".GSEA.pdf"), width = 12, height = 10)
  print(p_a)
  print(p_1)
  print(p_2)
  print(p_3)
  print(p_4)
  print(p_5)
  print(p_6)
  dev.off()
  path=paste0(savePath,"pathview/",str_split(basename(resFile), "\\.")[[1]][2])
  dir.create(path, recursive = T); setwd(path)
  for (i in seq(10)){
    suppressWarnings(pathview::pathview(gene.data = ae, pathway.id = af_kegg@result$ID[i], same.layer = F,
                                        pdf.size = c(12,12), kegg.native = F, species = "mmu", map.symbol = T, gene.annotpkg = "org.Mm.eg.db"))
  }
}
shy_gsea(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv",
         savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/GSEA/")
shy_gsea(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv",
         savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/GSEA/")
}
#11、画图：PCA [DIY]
{
#----------------------------------PCA-batch----------------------------------
#hisat2 + featureCounts + DESeq2 + ggplot2
path = "/ChIP_seq_1/sht/dot1L/shy_dot1ReviewRNAseq/"
pprr <- read.table("/RNA_seq_1/zx/XD/Count/gencode.vM21.primary_assembly.annotation.gtf_cut",header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)
dot1 <- read.table("/RNA_seq_2/sht/N1906264_GZZ_80-354364506/shyTEST/count/Dot1lRNAseq_KOvsWT.count", header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)
tblc <- read.table(paste0(path, "TBLCs/count/gene.count.txt"), header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)[,-c(1:5)]
tclc <- read.table(paste0(path, "TCLCs/count/gene.count.txt"), header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)[,-c(1:5)]
ercl <- read.table("/RNA_seq_1/zx/TOMATO/Count/gencode.vM21.primary_assembly.annotation.gtf_cut",
                   header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)#[,c(2,3)]
vtax <- cbind(pprr,dot1,tblc,tclc,ercl); colnames(vta)  <- c(paste0("i_wt",seq(2)),paste0("i_ko_a",seq(2)),paste0("i_ko_b",seq(2)),paste0("ko",seq(3)),paste0("wt",seq(3)),paste0("tb",seq(3)),paste0("tc",seq(4)),paste0("erc",seq(2)))
tabl <- data.frame(names = colnames(vtax),
                   condition = sapply(str_split(colnames(vta),"[0-9]"), "[", 1), 
                   batch = c(rep("A",2), rep("B",2), rep("C",2), rep("D",3), rep("E",3), rep("F",3), rep("G",4), rep("H",2)))
countmatrix <- as.matrix(vtax)
dds1 <- DESeqDataSetFromMatrix(countmatrix, colData = tabl, design = ~condition)[rowSums(counts(dds1))>5]
##https://cloud.tencent.com/developer/article/1625223
dds2 <- estimateSizeFactors(dds1)
rawx <- SummarizedExperiment(counts(dds2, normalized = F), colData=colData(dds2))
norx <- SummarizedExperiment(counts(dds2, normalized = T), colData=colData(dds2))
rldf <- rlog(dds2)
vsdf <- vst(dds2)
##去除批次效应
#View(mat)
#mat <- assay(rldf)
#mat <- limma::removeBatchEffect(mat, rldf$batch)
#assay(rldf) <- mat
pd_1 <- plotPCA(DESeqTransform(rawx), intgroup = c("condition"), returnData = T)
pd_2 <- plotPCA(DESeqTransform(norx), intgroup = c("condition"), returnData = T)
pd_3 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
pd_4 <- plotPCA(vsdf, intgroup = c("condition"), returnData = T)
shy_linshi <- function(df, i){
  mark = c("raw","normlization","rlog", "vst")[i]
  perc <- round(100 * attr(df, "percentVar"))
  df$condition <- factor(as.character(df$condition), 
                           levels = c("i_ko_a","i_ko_b","i_wt","ko","wt","erc","tb","tc"),
                           labels = c("KO74","KO92","mESC-pr","mESC-KO", "mESC-WT", "2CLC","TBLC","TLSC"))
  ppp <<- ggplot(df, aes(PC1, PC2, color=condition)) +
          geom_point(size=5) +
          xlab(paste0("PC1: ",perc[1],"% variance")) +
          ylab(paste0("PC2: ",perc[2],"% variance")) +
          scale_color_nejm() +
          theme_classic() +
          labs(color="Cell type",title = paste0("TITLE: plot of PCA ... ... ... ... ... ...",mark)) +
          theme(legend.text  = element_text(family = "serif", size = 12),
                legend.title = element_text(family = "serif", size = 12),
                legend.background = element_rect(fill = "gray95"),
                legend.position = c(0.8,0.75),
                axis.text    = element_text(family = "serif", size = 12),
                axis.title   = element_text(family = "serif", size = 12))
  print(ppp)
}
pdf(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/PCA/pr.PCA.pdf",width = 10, height = 8)
shy_linshi(df = pd_1, i = 1)
shy_linshi(df = pd_2, i = 2)
shy_linshi(df = pd_3, i = 3)
shy_linshi(df = pd_4, i = 4)
dev.off()
rm(pprr,dot1,tblc,tclc,ercl,dds1,dds2,ppp,pd_1,pd_2,pd_3,pd_4,rld1,rldf,vsdf,rawx,norx,tabl,vtax)
}
#12、画图：Corelation
{
shell <- paste0(path,"src/run_Correlation.sh")
cat("#!/bin/bash\n", file = shell)
aa <- list.files(paste0(path,"bamCoverage"), pattern = "*bw$", full.names = T)
ab <- list.files("/ChIP_seq_1/sht/dot1L/shy_dot1ReviewRNAseq/corelation/bw", pattern = "*bw$", full.names = T)
ac <- list.files("/RNA_seq_1/Pim3_KO/bamCoverage", pattern = "*CPM.bw$", full.names = T)
ad <- list.files("/RNA_seq_1/Ssrp1/bamCoverage", pattern = "*CPM.bw$", full.names = T)
cmd_a1 <- paste0("multiBigwigSummary bins -b ",paste0(c(aa,ab,ac,ad), collapse = " ")," -o ",
                 paste0(path,"PLOT/Correlation/","pr.Correlation.npz","\n"))
cmd_a2 <- paste0("plotCorrelation -in ",
                 paste0(path,"PLOT/Correlation/","pr.Correlation.npz "),
                 "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers --colorMap bwr -o ",
                 paste0(path,"PLOT/Correlation/","pr.Correlation.pdf "),
                 "--outFileCorMatrix ",paste0(path,"PLOT/Correlation/","pr.Correlation.inx.txt\n"))
for (m in cmd_a1){cat(m, file = shell, append = T)}
for (m in cmd_a2){cat(m, file = shell, append = T)}
rm(aa,ab,ac,ad,cmd_a1,cmd_a2)
print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
#13、画图：UCSC Track
{
  #1, upload bigwig file to cyverse, paste the links
  #2, make UCSC custom track
  #3, tuning parameters
  ####https://genome.ucsc.edu/s/hysun088/pr
}
#14、画图：Up TEs loci number [TElocal]
{
shy_lociNum <- function(resFile, saveFile){
  ress <- read.csv(file = resFile, header = T, row.names = 1, stringsAsFactors = F)
  tenm <- as.data.frame(table(sapply(str_split(string = rownames(ress),pattern = ":|_dup"), "[",1)))
  rela <- paste(sapply(str_split(string = rownames(ress),pattern = ":"), "[",2),
                sapply(str_split(string = rownames(ress),pattern = ":"), "[",3), sep = "__")
  tenm$clas <- sapply(str_split(as.character(as.data.frame(table(rela))[,1]),"__"),"[",2)
  tenm <- tenm[order(tenm$Freq,decreasing = T),]
  forP <- tenm[1:10,]
  forP$Var1 = factor(forP$Var1, levels = as.character(forP$Var1))
  p_1 <- ggplot(forP, aes(Var1, Freq, fill = clas)) + 
         geom_bar(stat = "identity",width = 0.6) +
         theme_classic() +
         theme(axis.title = element_text(family = "serif", color = "black", size = 14)) +
         theme(axis.text = element_text(size = 12, family = "serif", color = "black")) +
         theme(legend.title = element_blank(), 
               title = element_text(family = "serif", color = "black", size = 14)) +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) +
         theme(legend.position = c(0.7,0.7)) +
         ylab("Number of Loci") + 
         xlab(NULL) +
         labs(title = expression(Upregulated~TEs~"in"~italic(pr)^{"-/-"}~ESCs)) +
         scale_fill_manual(breaks = unique(forP$clas),values = ggsci::pal_aaas(palette = "default")(length(unique(forP$clas))))
  ggsave(plot = p_1, filename = saveFile, units = "cm", width = 12, height = 12)
}
shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.teSites.upup.csv",
            saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/TElociNum/pr.KO74.teSites.upup.pdf")
shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.teSites.upup.csv",
            saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/TElociNum/pr.KO92.teSites.upup.pdf")
shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.teSites.down.csv",
            saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/TElociNum/pr.KO74.teSites.down.pdf")
shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.teSites.down.csv",
            saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/TElociNum/pr.KO92.teSites.down.pdf")
}
#15、画图：MM fusion genes expression
{
##DESeq2 diff Gene
shy_diff_stringTie <- function(countFile, Class){
  mydata <- read.table(countFile, header = T, stringsAsFactors = F)
  counts <- na.omit(mydata[grep(x=mydata[,1],pattern = ":", invert=T),])
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  ##DESeq2
  yyyCOL <- data.frame(row.names = colnames(counts), Class)
  ddsdds <- DESeqDataSetFromMatrix(countData = counts, colData = yyyCOL, design = ~Class)
  ddsdds <- ddsdds[rowMeans(counts(ddsdds)) > 5,]
  ddsdds <- DESeq(ddsdds)
  resres <- results(ddsdds, contrast = c("Class","KO","WT"))
  resres <<- na.omit(as.data.frame(resres))
}
shy_diff_stringTie(countFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/count/v_74merge.cntTable", 
                   Class = c("KO","KO","WT","WT"))
##the Fusion Gene & PLOT
gtff <- import(inx4)
shy_fusion_stringTie <- function(gtfFile, whichTEs, mark, savePath){
  gf00 <- rtracklayer::import(gtfFile)
  #gf00 <- as.data.frame(gf00[which(gf00$type=="transcript")])
  #dfaa <- gf00 %>% 
  #        group_by(gene_id) %>% 
  #        summarise(seqnames=unique(seqnames), start=min(start), end=max(end))
  #gfaa <- makeGRangesFromDataFrame(dfaa, keep.extra.columns = T)
  gfaa <- gf00[which(gf00$type=="transcript")]
  gfbb <- gtff[grep(gtff$gene_id, pattern = paste0(whichTEs, ".*"))]
  oovv <- suppressWarnings(GenomicRanges::findOverlapslaps(query = gfaa, subject = gfbb))
  thid <- gfaa[unique(queryHits(oovv))]$gene_id
  resres$colr[rownames(resres) %in% thid] = "fusioned"
  resres$colr[!(rownames(resres) %in% thid)] = "oth"
  resA <- resres[which(resres$colr=="fusioned"),]
  resB <- resres[which(resres$colr=="oth"),]
  #开始作图
  p_1 <- ggplot() +
         theme_bw() +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border     = element_rect(colour = "black", size = 0.8),
               axis.line        = element_line(),
               axis.text        = element_text(family = "serif")) +
         xlab(label = expression(Log2~(Fold~Change~of~italic(pr)^{"-/-"}~ESCs))) +
         ylab(label = expression(Log10~(adjusted~italic(P)~plain(value)))) +
         geom_point(data=resB,aes(x=log2FoldChange,y=-log10(padj),colour="Other Genes"),size=1.2) +
         scale_color_manual(values = c("red","gray")) +
         geom_point(data=resA,aes(x=log2FoldChange,y=-log10(padj),colour=paste0(whichTEs,"-fushion Genes")),size=1.2) +
         theme(legend.position = c(0.2,0.9), 
               legend.background = element_rect(fill = NA),
               legend.spacing.x = unit(0,"mm"),
               legend.spacing.y = unit(0,"mm"),
               legend.key = element_blank(),
               legend.title = element_blank(),
               legend.text  = element_text(size = 13),
               axis.title   = element_text(size = 13)
         ) +
         geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed",colour="royalblue") +
         geom_vline(aes(xintercept=-log2(1.5)),   linetype="dashed",colour="royalblue") +
         geom_vline(aes(xintercept= log2(1.5)),   linetype="dashed",colour="royalblue") +
         guides(color = guide_legend(override.aes = list(size = 3),order=1))
  ggsave(plot = p_1, 
         filename = paste0(path, "PLOT/fusionTEs/",mark,"_fusion_",whichTEs,".pdf"), 
         width = 15, height = 15, units = "cm")
}
shy_fusion_stringTie(gtfFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf",
                     whichTEs = "MT2",#MT2#MM
                     mark = "KO74",
                     savePath = paste0(path, "PLOT/fusionTEs/"))
shy_fusion_stringTie(gtfFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO92.merge.gtf",
                     whichTEs = "MM",#MT2#MM
                     mark = "KO92",
                     savePath = paste0(path, "PLOT/fusionTEs/"))
}
#16、画图：repSites of selected TEs
{
shy_TEsites_plot <- function(resFile, saveFile){
  gf00 <- read.csv(file = resFile,
                   header = T, stringsAsFactors = F, row.names = 1)
  gf00 <- na.omit(gf00)
  gfaa <- gf00[grep(rownames(gf00), pattern = ".*MM-int.*"),]
  gfbb <- gf00[grep(rownames(gf00), pattern = ".*MT2.*"),]
  gfcc <- gf00[grep(rownames(gf00), pattern = ".*MM-int.*|.*MT2.*", invert = T),]
  #开始作图
  p_1 <- ggplot() +
         theme_bw() +
         theme(panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "black", size = 0.8),
               axis.line = element_line(), axis.text = element_text(family = "serif")) +
         xlab(label = expression(Log2~(Fold~Change~of~italic(pr)^{"-/-"}~ESCs))) +
         ylab(label = expression(Log10~(adjusted~italic(P)~plain(value)))) +
         scale_color_manual(breaks = c("MM","MT2","Others"),values = c("red","purple","gray")) +
         geom_point(data=gfcc,aes(x=log2FoldChange,y=-log10(padj),colour="Others"),size=1.2) +
         geom_point(data=gfbb,aes(x=log2FoldChange,y=-log10(padj),colour="MT2"),size=1.2) +
         geom_point(data=gfaa,aes(x=log2FoldChange,y=-log10(padj),colour="MM"),size=1.2) +
         theme(legend.position = c(0.2,0.9), 
               legend.background = element_rect(fill = NA),
               legend.spacing.x = unit(0,"mm"), legend.spacing.y = unit(0,"mm"),
               legend.key = element_blank(), legend.title = element_blank(),
               legend.text  = element_text(size = 13), axis.title   = element_text(size = 13)) +
         geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed",colour="royalblue") +
         geom_vline(aes(xintercept=-log2(1.5)),   linetype="dashed",colour="royalblue") +
         geom_vline(aes(xintercept= log2(1.5)),   linetype="dashed",colour="royalblue") +
         guides(color = guide_legend(override.aes = list(size = 3),order=1))
  ggsave(plot = p_1, filename = saveFile, width = 15, height = 15, units = "cm")
}
shy_TEsites_plot(resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.teSites.csv",
                 saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/TEsitesPlot/KO74.MM_mt2.pdf")
shy_TEsites_plot(resFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.teSites.csv",
                 saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/PLOT/TEsitesPlot/KO92.MM_mt2.pdf")
}
##
##
##
#17、看看pr RNA-seq有多少reads落在rDNA
{
  inx6 = "/Reference/aaSHY/zOther/Rn45s/INDEX/hisat2/Rn45s"
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202306/"
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/"
  prex <- unique(sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*fq.gz$*"),"_"), "[",1))
  cmd_a <- paste0("hisat2 -q --very-sensitive -p 20 -k 5 -x ",inx6,
                  " -1 ",paste0(path,"trim/",prex,"*.R1_val_1.fq.gz"),
                  " -2 ",paste0(path,"trim/",prex,"*.R2_val_2.fq.gz"),
                  " --no-mixed --no-discordant",
                  " |samtools view -F 4 -b |samtools sort > ",
                  paste0(path,"bam/rn45s_",prex,".bam"),"\n")
  cmd_b <- paste0("samtools index ",paste0(path,"bam/rn45s_",prex,".bam"),"\n")
  cmd_c <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 --samFlagExclude 256 --outFileFormat bigwig --binSize 20 -b ",
                  paste0(path,"bam/rn45s_",prex,".bam -o "),
                  paste0(path,"bamCoverage/rn45s_",prex,".bcv.bw\n"))
  shell <- paste0("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/run_rn45s.sh")
  cat("#!/bin/bash\n", file = shell)
  for (i in cmd_a) {cat(i, append = T, file = shell)}
  for (i in cmd_b) {cat(i, append = T, file = shell)}
  for (i in cmd_c) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_rn45s.log"," 2>&1 &"))
}
{
  cmd_7 <- paste0("computeMatrix reference-point --referencePoint center -p 30 -a ",
                  1000," -b ",1000," --missingDataAsZero --skipZeros ",
                  "-R ", inx2, " -S ",
                  paste0(path,"bamCoverage/", prex,".bw", collapse = " ")," -o ",
                  paste0(path,"profile/",basename(inx2),".mat.gz"),"\n")
  cmd_8 <- paste0("plotProfile --plotHeight 8 --plotWidth 8 -m ",
                  paste0(path,"profile/",basename(inx2),".mat.gz")," -out ",
                  paste0(path,"profile/",basename(inx2),".pdf "),
                  "--perGroup","\n")
}
##
##
##
#18、看看pr RNA-seq有多少MM-int loci是上调的
{
  loci <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.teSites.upup.csv", header = T, stringsAsFactors = F)
  loci <- loci[grep(x = loci$X, pattern = "MM-int"),]
  upup <- sapply(str_split(loci$X, pattern = ":"),"[",1)
  tefs <- rtracklayer::import(con = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf", format = "gtf")
  tefs <- tefs[grep(x = tefs$gene_id, pattern = "MM-int")]
  fina <- as.data.frame(tefs[which(tefs$transcript_id %in% upup)])
  fina$source <- paste0("chr",fina$seqnames, ":", fina$start, "-", fina$end)
  write.csv(x = fina[,c(6,4,5,11)], file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.MMintSites.upup.csv", quote = F, row.names = F, col.names = F)
}
##
##
##
#19、对pr KO后激活的2细胞基因排个序
{
  pprr <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv", sep = ",", header = T, row.names = 1, stringsAsFactors = F)
  twoc <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed", sep = "\t", header = F, stringsAsFactors = F)
  wr01 <- pprr[which(rownames(pprr) %in% twoc$V4),]
  wr01 <- wr01[order(wr01$log2FoldChange, decreasing = T),]
  wr01$geneName <- rownames(wr01)
  openxlsx::write.xlsx(x = wr01, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.2C.gene.xlsx", asTable = T)
}
#######################################################################
#22、可变剪接
{
  bam1 <- list.files(paste0(path,"bam"), pattern = "^KO-92.*bam$",full.names = T)
  bam2 <- list.files(paste0(path,"bam"), pattern = "^E14.*bam$",full.names = T)
  cat(file = paste0(path,"rMATS/KO92/opt.txt"), sep = ",", bam1)
  cat(file = paste0(path,"rMATS/KO92/ctr.txt"), sep = ",", bam2)
  paste0("rmats.py --nthread 10 -t paired --readLength 140 ",
         "--gtf ",inx3," --b1 ",path,"rMATS/KO92/opt.txt ","--b2 ",path,
         "rMATS/KO92/ctr.txt"," --od ",path,"rMATS/KO92/ --tmp ",path,"rMATS/KO92/tmp","\n")
}
#######################################################################
#23、python GSEA
{
  temp <- read.csv(paste0(path, "aaprKO74/DESeq2/prKO.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, file = paste0(path, "aaprKO74/DESeq2/KO74.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "aaprKO74/DESeq2/KO74.res.rnk"), " -p 10 -o ", 
                   paste0(path, "aaprKO74/PLOT/pythonGSEA/"))#conda activate base gseapy 0.10.8
}
################################################################################
#24、MM-融合转录本、1165基因的重合, 得到159个基因, 这些基因在ff KD和pr KD的热图
{
  ge01 <- import(inx3)
  ge01 <- ge01[which(ge01$type == "gene")]
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  myta <- ge01[which(ge01$gene_name %in% myta$name)]
  fusg <- read.table(file = "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed", sep = "\t")
  fusg <- makeGRangesFromDataFrame(df = fusg, seqnames.field = "V1", start.field = "V2",
                                   end.field = "V3", strand.field = "V6")
  ov01 <- findOverlaps(query = ge01, subject = fusg, ignore.strand = T)
  tmp1 <- ge01[unique(queryHits(ov01))]
  #length(unique(tmp1$gene_name))：1490个MM融合基因
  ov02 <- findOverlaps(query = myta, subject = fusg, ignore.strand = T)
  tmp2 <- as.data.frame(myta[unique(queryHits(ov02))])
  wr02 <- tmp2[,c(1,2,3,12,4,5)]
  wr02$width <- "."
  write.table(x = wr02, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/tagMM.MMFus.159.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
  ne02 <- unique(wr02$gene_name)
}
{
  pprr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/count/prKO.cpm.onlyGene.xlsx")
  pprr <- pprr[which(pprr$theName %in% ne02),]
  rownames(pprr) <- pprr$theName; pprr <- pprr[,c(3,4,1,2)]
  ffs <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/202306/count/ff.cpm.onlyGene.xlsx")
  ffs <- ffs[which(ffs$theName %in% ne02),]
  rownames(ffs) <- ffs$theName; ffs <- ffs[,c(3,4,1,2)]; colnames(ffs)[1:2] <- c("ffWT-1","ffWT-2")
  cnt1 <- cbind(pprr, ffs)
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nthe 159 genes -- pr-ff",
           scale = "row", cellwidth = 30,
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}

