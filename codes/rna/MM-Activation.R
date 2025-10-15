#write by ~~~ at 2024.03.05
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/STARv279a"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/STARv279a"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie/count","bamCoverage","trim","DESeq2","dumpROOM","bam","count","PLOT")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*1.fq.gz$*"),"_R1.fq.gz"), "[",1)
  lay1 <- paste0(path,"fq/",prex,"_R1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_R2.fq.gz")
  cmd_01 <- paste0("#trim_galore --phred33 -j 8 ",
                   "-o ", path, "trim --paired ",lay1," ",lay2,"\n")
  lay1 <- paste0(path,"trim/",prex,"_R1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_R2_val_2.fq.gz")
  cmd_03 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_04 <- paste0("#mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_05 <- paste0("#samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_06 <- paste0("#bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_07 <- paste0("#TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:3],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[4:6],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ccMM- --outdir ",paste0(path,"count"),"\n")
  cmd_08 <- paste0("#TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_09 <- paste0("#stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_10 <- paste0("#stringtie --merge -o ",
                   path,"stringtie/MM-.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_11 <- paste0("#cat ", paste0(path,"stringtie/MM-.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/MM-.merge.R.gtf"),"\n")
  cmd_12 <- paste0("#TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:3],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[4:6],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/MM-.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project MM-Merge ",
                   "--outdir ",path,"stringtie/count","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_03) {cat(i, file = shell, append = T)}
  for (i in cmd_04) {cat(i, file = shell, append = T)}
  for (i in cmd_05) {cat(i, file = shell, append = T)}
  for (i in cmd_06) {cat(i, file = shell, append = T)}
  for (i in cmd_07) {cat(i, file = shell, append = T)}
  for (i in cmd_08) {cat(i, file = shell, append = T)}
  for (i in cmd_09) {cat(i, file = shell, append = T)}
  for (i in cmd_10) {cat(i, file = shell, append = T)}
  for (i in cmd_11) {cat(i, file = shell, append = T)}
  for (i in cmd_12) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
#######################################################################
#04、DESeq2：Gene、TE、TE sites
{
  shy_diff <- function(mark, geneFile, siteFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = VS)
    write.csv(geneAA, paste0(path,"DESeq2/",mark, ".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = VS)
    write.csv(teOnly, paste0(path,"DESeq2/",mark, ".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.down.csv"))
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
    colLen <- ceiling(length(colnames(tmp2))/2)
    teGene <- tmp2
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.down.csv"))
  }
  shy_diff(mark = "MM-", 
           geneFile = paste0(path, "count/ccMM-.cntTable"),
           siteFile = paste0(path, "count/ccMM-Site.cntTable"),
           Class = factor(c("OE","OE","OE","WT","WT","WT")), VS = c("Class","OE","WT"))
}
#######################################################################
#05、CPM Excel
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("MMKD-1","MMKD-2","MM-KD-1","MM-KD-2","ctrl-1","ctrl-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cpm.allGeneTE.xlsx")
  #
  cnt2 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/ccMM-.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt2) <- c("MM-Ac-1","MM-Ac-2","MM-Ac-3","ctrlAc-1","ctrlAc-2","ctrlAc-3")
  tmp1 <- cnt2[grep(x = rownames(cnt2), pattern=":", invert = F),]#TE
  tmp2 <- cnt2[grep(x = rownames(cnt2), pattern=":", invert = T),]#Gene
  tmp3 <- cnt2#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/MM-Ac.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/MM-Ac.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/MM-Ac.cpm.allGeneTE.xlsx")
}
#######################################################################
#06、Volcano for Genes
{
  shy_volcano_gene <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold[resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05] = "Upregulated genes"
    resdata$threshold[resdata$log2FoldChange < -log2(1.5) & resdata$padj <0.05] = "Downregulated genes"
    resdata$threshold[!(resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05) & !(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05)] = "Other genes"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated genes","Downregulated genes","Other genes"))
    p_1 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-6, 9), breaks = seq(-6, 9, 2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
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
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/MM-.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/MM-.gene.volcano.pdf"))
}
#######################################################################
#07、Volcano for TEs
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
      labb <- labb[which(labb$padj<10^(-5)),]#or -10
    }
    p_2 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-4,4)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +#ylim(0,max(-log10(resdata$padj)[is.finite(-log10(resdata$padj))])+20) +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.2,0.8), 
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
  shy_volcano_te(resFile  = paste0(path, "DESeq2/MM-.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/MM-.te.GeneBackGround.volcano.pdf"))
}
#######################################################################
#08、PCA
{
  #----------------------------------PCA-batch----------------------------------
  cnt1 <- read.table(paste0(path, "count/ccMM-.cntTable"),
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
  colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
  #
  base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
  myta <- cbind(cnt1, forPCA)
  #
  condition = sapply(str_split(colnames(myta),"-"), "[", 1)
  btch = c(rep("ours",length(colnames(cnt1))), condition[c(length(colnames(cnt1))+1):length(myta)])
  clas <- data.frame(names = colnames(myta), condition = condition, batch = btch)
  dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(myta, colData = clas, design = ~condition))#[rowSums(counts(vtax))>5]
  dds2 <- BiocGenerics::estimateSizeFactors(dds1)
  rldf <- DESeq2::rlog(dds2)
  pd_1 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
  rldf <- DESeq2::vst(dds2)
  pd_2 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)
  #去除批次效应
  #SummarizedExperiment::assay(rldf) <- limma::removeBatchEffect(assay(rldf), batch=rldf$batch)
  #https://cloud.tencent.com/developer/article/1625223
  #豆包说“通常情况下，数据集小于30个样品时可以用rlog，数据集大于30个样品时用vst
  #vsdf <- vst(dds2)
  #rawx <- SummarizedExperiment(counts(dds2, normalized = F), colData=colData(dds2))
  #norx <- SummarizedExperiment(counts(dds2, normalized = T), colData=colData(dds2))
  #pd_3 <- plotPCA(DESeqTransform(rawx), intgroup = c("condition"), returnData = T)
  #pd_4 <- plotPCA(DESeqTransform(norx), intgroup = c("condition"), returnData = T)
  shy_linshi <- function(df, i){
    mark <- c("rlog", "vst")[i]
    perc <- round(100 * attr(df, "percentVar"))
    df$condition <- factor(as.character(df$condition), levels = as.character(df$condition),
                           labels = as.character(df$condition))
    ppp <<- ggplot(df, aes(PC1, PC2, color=condition)) +
      geom_point(size=5) +
      xlab(paste0("PC1: ",perc[1],"% variance")) +
      ylab(paste0("PC2: ",perc[2],"% variance")) +
      scale_color_d3() +
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
  pdf(file = paste0(path,"PLOT/PCA.pdf"),width = 10, height = 8)
  shy_linshi(df = pd_1, i = 1)
  shy_linshi(df = pd_2, i = 2)
  dev.off()
}
#######################################################################
#09、GO BP, KEGG, ReactomePA
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
    kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
    reac$Description <- gsub(x = reac$Description,fixed = T, pattern = "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)", replacement = "Regulation of IGF transport and uptake by IGFBPs")
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
           filename = paste0(savePath, change, ".KEGG.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
    ggsave(plot = p_2, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".ReactomePA.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
    ggsave(plot = p_3, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
  }
  shy_fuji_bar(resFile  = paste0(path, "DESeq2/MM-.geneOnly.csv"),
               savePath = paste0(path, "PLOT/Enrichment/"),
               change   = "upup")
  shy_fuji_bar(resFile  = paste0(path, "DESeq2/MM-.geneOnly.csv"),
               savePath = paste0(path, "PLOT/Enrichment/"),
               change   = "down")
}
#######################################################################
#10、pheatmap
{
  shy_pheatmap <- function(readsFile, gene_resFile, tete_resFile, savePath, mark){
    print("reads counted by TEtranscripts")
    #
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = ".bam"), "[", 1)
    gene <- data[grep(rownames(data), pattern = ":", invert = T),]
    gene <- gene[rowSums(gene) > 5,]
    gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
    tete <- data[grep(rownames(data), pattern = ":", invert = F),]
    tete <- tete[rowSums(tete) > 5,]
    tete <- as.data.frame(apply(tete, 2, function(x){x/sum(x)*1000000}))
    #
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
    pdf(paste0(savePath, str_split(basename(readsFile), pattern = "\\.")[[1]][1],".by_",mark, "_pheatmap.v2.pdf"), 
        width = 8, height = 10)
    ##gene表达量热图
    pheatmap(mat = gene,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##siginificant gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG),],
             main="the heatmap of significant Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Upregulated gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG_upup),],
             main="the heatmap of Upregulated Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Downregulated gene表达量热图
    pheatmap(mat = gene[which(rownames(gene) %in% sigG_down),],
             main="the heatmap of Downregulated Genes",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(gene), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##tete表达量热图
    pheatmap(mat = tete,
             main="the heatmap of all TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##siginificant tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT),],
             main="the heatmap of significant TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Upregulated tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT_upup),],
             main="the heatmap of Upregulated TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    ##Downregulated tete表达量热图
    pheatmap(mat = tete[which(rownames(tete) %in% sigT_down),],
             main="the heatmap of Downregulated TEs",
             scale = "row", cellwidth = 60, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
             labels_col = colnames(tete), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
    dev.off()
  }
  shy_pheatmap(mark = "heat", readsFile = paste0(path, "count/ccMM-.cntTable"),
               gene_resFile = paste0(path, "DESeq2/MM-.geneOnly.csv"),
               tete_resFile = paste0(path, "DESeq2/MM-.te.GeneBackGround.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
#######################################################################
#11、GSEA
{
  shy_gsea <- function(resFile, savePath){
    aa <- read.csv(file = resFile, header = T, stringsAsFactors = F)[,c(1,3)]
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
    af_kegg <- suppressWarnings(clusterProfiler::gseKEGG(geneList = ae, maxGSSize = 1000, organism = "mmu", pvalueCutoff = 1, use_internal_data = F, seed = T, by = "fgsea"))
    temp <- af_kegg@result$Description
    temp <- sapply(str_split(temp, pattern = " - Mus"), "[", 1)
    af_kegg@result$Description <- temp
    testxxx <- af_kegg@result
    testxxx[,10] <- apply(af_kegg@result, 1, function(x){chartr(old = ", ", new = "; ", x[10])})
    write.csv(testxxx, row.names = F, quote = F,
              file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".result.KEGG.csv"))
    ##2-cell Genes
    aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt")
    aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
    ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
    af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab, pvalueCutoff = 1))
    rm(aa,ab)
    #af_reac <- suppressWarnings(ReactomePA::gsePathway(geneList = ae, maxGSSize = 1000, eps = NA, organism = "mouse", pvalueCutoff = 1, verbose = T, seed = T, by = "fgsea"))
    ##开始画图
    #dev.new()
    p_a <- enrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for 2-cell Genes", pvalue_table = T)
    p_1 <- enrichplot::gseaplot2(x = af_gogo, 1:5, title = "GO_BP_GSEA: top 5", rel_heights = c(1.5, 0.3, 0.6), color = ggsci::pal_lancet()(5), pvalue_table = T)
    p_2 <- enrichplot::ridgeplot(x = af_gogo, showCategory = 50, fill = "pvalue", label_format = 100)
    p_3 <- enrichplot::cnetplot(setReadable(af_gogo, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID"), color.params = list(foldChange = ae))
    p_4 <- enrichplot::gseaplot2(x = af_kegg, 1:5, title = "KEGG__GSEA: top 5", rel_heights = c(1.5, 0.3, 0.6), color = ggsci::pal_lancet()(5), pvalue_table = T)
    p_5 <- enrichplot::ridgeplot(x = af_kegg, showCategory = 50, fill = "pvalue", label_format = 100)
    p_6 <- enrichplot::cnetplot(setReadable(af_kegg, OrgDb = "org.Mm.eg.db", keyType = "ENTREZID"), color.params = list(foldChange = ae))
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
  shy_gsea(resFile = paste0(path, "DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv"), 
           savePath = paste0(path, "PLOT/EnrichmentGSEA/"))
  {
    temp <- read.csv(paste0(path, "DESeq2/MM-.geneOnly.csv"), header = T)
    temp <- temp[, c("X", "log2FoldChange")]
    temp <- temp[order(temp$log2FoldChange, decreasing = T),]
    write.table(x = temp, file = paste0(path, "DESeq2/MM-.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
    cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                     paste0(path, "DESeq2/MM-.res.rnk"), " -p 10 -o ", 
                     paste0(path, "PLOT/EnrichmentGSEA/"))#conda activate base gseapy 0.10.8
  }
}
#######################################################################
#12、Corelation
{
  shell <- paste0(path,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  bw01 <- list.files(paste0(path,"bamCompare"), pattern = "*bw$", full.names = T)
  bw02 <- list.files(path = "/ChIP_seq_1/aaSHY/likeTwoC",
                     recursive = T, pattern = ".bw$", full.names = T)
  bw02 <- bw02[grep(bw02, pattern="bamCompare")]
  cmd_01 <- paste0("multiBigwigSummary bins ", 
                   "-p 20 -b ",paste0(c(bw01,bw02), collapse = " "),
                   " --smartLabels",
                   " -o ",paste0(path,"PLOT/Correlation/gata3.cor.npz","\n"))
  cmd_02 <- paste0("plotCorrelation -in ",path,
                   "PLOT/Correlation/gata3.cor.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", path,"PLOT/Correlation/gata3.cor.pdf ",
                   "--outFileCorMatrix ",path,"PLOT/Correlation/gata3.cor.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
#######################################################################
#13、UCSC Track
{
  prex <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*1.fq.gz$*"),"_R1.fq.gz"), "[",1)
  chro <- c(seq(19),"X","Y")
  for (i in prex) {
    print(i)
    tmp1 <- rtracklayer::import(con = paste0(path,"stringtie/",i,".gtf"))
    tmp2 <- as.data.frame(tmp1[which(tmp1@seqnames %in% chro)])
    tmp2$seqnames <- paste0("chr",tmp2$seqnames)
    tmp3 <- makeGRangesFromDataFrame(df = tmp2, keep.extra.columns = T)
    rtracklayer::export(object = tmp3, con = paste0(path,"stringtie/ucscTrack/",i,".gtf"), format = "gtf")
  }
}
#######################################################################
#14、MA plot
{
  shy_teMA <- function(difffile, savefile){
    rf01 <- read.csv(file = difffile, row.names = 1)
    rf01$subFam <- sapply(str_split(rownames(rf01),pattern = ":"),"[",1)
    rf01$family <- sapply(str_split(rownames(rf01),pattern = ":"),"[",2)
    rf01$classs <- sapply(str_split(rownames(rf01),pattern = ":"),"[",3)
    rf01$colo <- "other"
    for (i in c("ERV1","ERVK","ERVL","Gypsy")){
      rf01$colo[which(rf01$padj<0.05 & abs(rf01$log2FoldChange)>log2(1.5) & rf01$family==i)] <- i
    }
    rf02 <- rf01
    rf02$colo <- factor(rf02$colo,levels=c("ERV1","ERVK","ERVL","Gypsy","other"))
    myCo <- list(ERV1="sienna1",ERVK="steelblue3",ERVL="palevioletred2",Gypsy="green",Others="grey80")
    myCo <- unlist(myCo)
    #
    p_01 <- ggplot() +
      geom_point(data = rf02, aes(x=log2(baseMean+1), y=log2FoldChange, color=colo, shape=NULL),size=2.5) +
      guides(color = guide_legend(title = "repFamliy")) +
      labs(x = "log2(Expression in baseMean)", y = "log2Foldchange") +
      geom_hline(aes(yintercept =  log2(1.5)), linetype="dashed", colour="grey") +
      geom_hline(aes(yintercept = -log2(1.5)), linetype="dashed", colour="grey") +
      geom_text_repel(data = rf02[which(rf02$colo!="other"),], 
                      aes(x=log2(baseMean+1), y=log2FoldChange, label=subFam), 
                      hjust=-0.125, size=4, colour="black", family="serif") +
      scale_color_manual(values = myCo)+
      scale_y_continuous(limits = c(floor(min(rf02$log2FoldChange)),ceiling(max(rf02$log2FoldChange))),
                         breaks = seq(floor(min(rf02$log2FoldChange)),ceiling(max(rf02$log2FoldChange)))) +
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(),
            axis.text = element_text(family = "serif"),
            panel.border = element_rect(colour = "black", linewidth = 0.8))+
      theme(legend.position.inside = c(0.91,0.25),legend.background = element_blank(),
            legend.text = element_text(family = "serif"))
    ggsave(plot = p_01,
           units = "cm", height = 12, width = 13.5577, filename = savefile)
  }
  shy_teMA(difffile = paste0(path,"DESeq2/rmSampAc2Ctrl1/MM-.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/MM-.te.GeneBackGround.MA.r3.pdf"))
}
#######################################################################
#15、MA plot调整为CPM点图 [以及后面去掉了control 1和activation 2样本]
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("MMKD-1","MMKD-2","MM-KD-1","MM-KD-2","ctrl-1","ctrl-2")
  cnt1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]
  cnt1 <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
  cnt2 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/ccMM-.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt2) <- c("MM-Ac-1","MM-Ac-2","MM-Ac-3","ctrlAc-1","ctrlAc-2","ctrlAc-3")
  cnt2 <- cnt2[grep(x = rownames(cnt2), pattern=":", invert = F),]
  cnt2 <- as.data.frame(apply(cnt2, 2, function(x){x/sum(x)*1000000}))
  cnts <- cbind(cnt1[,c(1,2,5,6)], cnt2[,c(1,3,5,6)])
  myta <- data.frame(MMKD = log2(rowMeans(cnts[,c(1,2)]) +1),
                     MMNC = log2(rowMeans(cnts[,c(3,4)]) +1),
                     MM-Acti = log2(rowMeans(cnts[,c(5,6)]) +1),
                     MM-NCaa = log2(rowMeans(cnts[,c(7,8)]) +1))
  #
  {
    myte <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyTE.MM.csv", row.names = 1)#不带基因背景
    myte$mark <- NA
    myte$mark[which(myte$log2FoldChange > 0 & myte$pvalue < 0.05)] <- "upup"
    myte$mark[which(myte$log2FoldChange < 0 & myte$pvalue < 0.05)] <- "down"
    te01 <- rownames(myte[which(myte$mark=="upup"),])
    te02 <- rownames(myte[which(myte$mark=="down"),])
    myte <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.teOnly.csv", row.names = 1)
    myte$mark <- NA
    myte$mark[which(myte$log2FoldChange > 0 & myte$pvalue < 0.05)] <- "upup"
    myte$mark[which(myte$log2FoldChange < 0 & myte$pvalue < 0.05)] <- "down"
    te03 <- rownames(myte[which(myte$mark=="upup"),])
    te04 <- rownames(myte[which(myte$mark=="down"),])
    #
    myta$MMMark <- "other"
    myta$MMMark[which(rownames(myta) %in% te01)] <- "upup"
    myta$MMMark[which(rownames(myta) %in% te02)] <- "down"
    myta$MM-Mark <- "other"
    myta$MM-Mark[which(rownames(myta) %in% te03)] <- "upup"
    myta$MM-Mark[which(rownames(myta) %in% te04)] <- "down"
  }
  labe <- myta[grep(x = rownames(myta), pattern = "MM-int|MM-_Mm"),]
  myta$MMMark <- factor(x = myta$MMMark, levels = c("upup","down","other"))
  p_07 <- ggplot() +
    geom_point(data = myta[which(myta$MMMark =="other"),], aes(x = MMNC, y = MMKD), color = "gray") +
    geom_point(data = myta[which(myta$MMMark !="other"),], aes(x = MMNC, y = MMKD, color = MMMark)) +
    scale_color_manual(values = c("red","blue")) +
    geom_point(data = labe, aes(x = MMNC, y = MMKD),size=2,shape=2,color = "blue") +
    geom_text_repel(data = labe, aes(x = MMNC, y = MMKD, label=rownames(labe))) +
    labs(x="Log2(CPM+1)(NC)",y="Log2(CPM+1)(MM KD)",title = "TE: MM KD vs Control") +
    guides(color = guide_legend(title = "Condition")) +
    theme_few() +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          plot.title = element_text(size = 16,family = "serif",color = "black"),
          legend.position = c(0.8,0.2),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black")) +
    geom_abline(slope = 1,intercept = 0,lty="dashed")
  ggsave(plot = p_07, 
         units = "cm", width = 12, height = 12,
         filename = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/PLOT/MMKD.cpm.point.pdf")
  #
  heta <- myta
  heta$MM-Mark <- factor(x = heta$MM-Mark, levels = c("upup","down","other"))
  p_08 <- ggplot() +
    geom_point(data = heta[which(heta$MM-Mark =="other"),], aes(x = MM-NCaa, y = MM-Acti), color = "gray") +
    geom_point(data = heta[which(heta$MM-Mark !="other"),], aes(x = MM-NCaa, y = MM-Acti, color = MM-Mark)) +
    scale_color_manual(values = c("red","blue")) +
    geom_point(data = labe, aes(x = MM-NCaa, y = MM-Acti),size=2,shape=2,color = "red") +
    geom_text_repel(data = labe, aes(x = MM-NCaa, y = MM-Acti, label=rownames(labe))) +
    labs(x="Log2(CPM+1)(NC)",y="Log2(CPM+1)(MM- Activation)",title = "TE: MM- Activation vs Control") +
    guides(color = guide_legend(title = "Condition")) +
    theme_few() +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          plot.title = element_text(size = 16,family = "serif",color = "black"),
          legend.position = c(0.8,0.2),
          axis.title = element_text(size = 12,family = "serif",color = "black"),
          axis.text  = element_text(size = 12,family = "serif",color = "black")) +
    geom_abline(slope = 1,intercept = 0,lty="dashed")
  ggsave(plot = p_08, 
         units = "cm", width = 12, height = 12,
         filename = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/PLOT/MM-Ac.cpm.point.pdf")
}
#######################################################################
#16、MA plot TE sites
{
  shy_repSite <- function(resdFile, saveFile){
    gf00 <- read.csv(file = resdFile, header = T, stringsAsFactors = F, row.names = 1)
    gf00 <- na.omit(gf00)
    gfaa <- gf00[grep(rownames(gf00), pattern = ".*MM-int.*"),]
    gfbb <- gf00[grep(rownames(gf00), pattern = ".*MM-_Mm.*"),]
    gfcc <- gf00[grep(rownames(gf00), pattern = ".*MM-int.*|.*MM-_Mm.*", invert = T),]
    #开始作图
    p_1 <- ggplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 0.8),
            axis.line = element_line(), axis.text = element_text(family = "serif")) +
      ylab(label = expression(Log2Fold~Change)) +
      xlab(label = expression(Log2~(baseMean + 1))) +
      scale_color_manual(breaks = c("MM-int","MM-_Mm","Others"),values = c("red","purple","gray")) +
      geom_point(data=gfcc,aes(x=log2(baseMean+1),y=log2FoldChange,colour="Others"),size=1.2) +
      geom_point(data=gfbb,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MM-_Mm"),size=1.2) +
      geom_point(data=gfaa,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MM-int"),size=1.2) +
      theme(legend.position = c(0.7,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"), legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(), legend.title = element_blank(),
            legend.text  = element_text(size = 13), axis.title   = element_text(size = 13)) +
      geom_hline(aes(yintercept=-log2(1.5)),   linetype="dashed",colour="royalblue") +
      geom_hline(aes(yintercept= log2(1.5)),   linetype="dashed",colour="royalblue") +
      guides(color = guide_legend(override.aes = list(size = 3),order=1))
    ggsave(plot = p_1, filename = saveFile, width = 15, height = 15, units = "cm")
  }
  shy_repSite(resdFile = paste0(path,"DESeq2/MM-.te.site.GeneBackGround.csv"),
              saveFile = paste0(path,"PLOT/MM-.te.site.GeneBackGround.MM_MM-.pdf"))
}
#######################################################################
#17、Up TEs loci number [TElocal]
{
  shy_lociNum <- function(resdFile, saveFile){
    ress <- read.csv(file = resdFile, header = T, row.names = 1, stringsAsFactors = F)
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
      labs(title = expression(deRegulated~TEs~"in"~italic(MM-)^{"+/+"}~ESCs)) +
      scale_fill_manual(breaks = unique(forP$clas),values = ggsci::pal_aaas(palette = "default")(length(unique(forP$clas))))
    ggsave(plot = p_1, filename = saveFile, units = "cm", width = 12, height = 12)
  }
  shy_lociNum(resdFile = paste0(path, "DESeq2/MM-.te.site.GeneBackGround.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/MM-.lociNum.upup.pdf"))
  shy_lociNum(resdFile = paste0(path, "DESeq2/MM-.te.site.GeneBackGround.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/MM-.lociNum.down.pdf"))
}
#######################################################################
#18、MM fusion genes expression
{
  shy_fusion_stringTie <- function(countFile, Class, gtfFile, whichTE, contrast){
    tmp1 <- read.table(countFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = contrast)
    resres <- na.omit(as.data.frame(geneAA))
    #
    gtff <- rtracklayer::import(inx4)
    gf00 <- rtracklayer::import(gtfFile)
    gfaa <- gf00[which(gf00$type=="transcript")]
    gfbb <- gtff[grep(gtff$gene_id, pattern = paste0(whichTE, ".*"))]
    oovv <- suppressWarnings(GenomicRanges::findOverlaps(query = gfaa, subject = gfbb))
    thid <- gfaa[unique(queryHits(oovv))]$gene_name
    resres$colr[rownames(resres) %in% thid] = "fusioned"
    resres$colr[!(rownames(resres) %in% thid)] = "oth"
    resA <- resres[which(resres$colr=="fusioned"),]
    resB <- resres[which(resres$colr=="oth"),]
    #
    p_02 <- ggplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 0.8),
            axis.line = element_line(), axis.text = element_text(family = "serif")) +
      xlab(label = expression(Log2~(Fold~Change~of~italic(MM-)^{"+/+"}~ESCs))) +
      ylab(label = expression(-Log10~(adjusted~italic(P)~plain(value)))) +
      geom_point(data=resB,aes(x=log2FoldChange,y=-log10(padj),colour="Other Genes"),size=1.2) +
      scale_color_manual(values = c("red","gray")) +
      scale_x_continuous(limits = c(-6, 10)) +
      geom_point(data=resA,aes(x=log2FoldChange,y=-log10(padj),colour=paste0(whichTE,"-fushion Genes")),size=1.2) +
      theme(legend.position = c(0.3,0.9), 
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
    ggsave(plot = p_02, 
           width = 15, height = 15, units = "cm",
           filename = paste0(path, "PLOT/chimeric/MM-.chimeric.",whichTE,".pdf"))
  }
  shy_fusion_stringTie(countFile = paste0(path, "stringtie/count/MM-Merge.cntTable"), 
                       Class = c("OE","OE","OE","WT","WT","WT"),
                       contrast = c("Class","OE","WT"),
                       gtfFile = paste0(path, "stringtie/MM-.merge.R.gtf"),
                       whichTE = "MM-int")#MM-_Mm#MM-int)
}
#######################################################################
#19、把MM- activation上调、下调的基因/TE和MM KD下调、上调的基因/TE重合下
#先计算
#做个PCA(所有Gene), 去除离群的样本 [保留activation 1,3、control 2,3]
{
  cnt2 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/ccMM-.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt2) <- c("MM-Ac-1","MM-Ac-2","MM-Ac-3","ctrlAc-1","ctrlAc-2","ctrlAc-3")
  dta3 <- cnt2[grep(rownames(cnt2),pattern=":",invert = T),]
  tabl <- data.frame(names = colnames(dta3),
                     condition = sapply(str_split(colnames(dta3),"-"), "[", 1), 
                     batch = c(rep("A",3), rep("C",3)))
  countmatrix <- as.matrix(dta3)
  dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(countmatrix, colData = tabl, design = ~condition))#[rowSums(counts(vtax))>5]
  ##https://cloud.tencent.com/developer/article/1625223
  dds2 <- BiocGenerics::estimateSizeFactors(dds1)
  #rawx <- SummarizedExperiment(counts(dds2, normalized = F), colData=colData(dds2))
  #norx <- SummarizedExperiment(counts(dds2, normalized = T), colData=colData(dds2))
  rldf <- rlog(dds2)#vsdf <- vst(dds2), 豆包说“通常情况下，数据集小于30个样品时可以用rlog，数据集大于30个样品时用vst。”
  plotPCA(rldf, intgroup = c("condition"), returnData = F)
  plotPCA(rldf, intgroup = c("condition"), returnData = T)
}
#基于组合样本，做差异分析
{
  shy_diff <- function(mark, geneFile, siteFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    tmp1 <- tmp1[,c(1,3,5,6)]
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = VS)
    write.csv(geneAA, paste0(path,"DESeq2/rmSampAc2Ctrl1/",mark, ".geneOnly.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = VS)
    write.csv(teOnly, paste0(path,"DESeq2/rmSampAc2Ctrl1/",mark, ".teOnly.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    geneTE <- teGene[grep(rownames(teGene), pattern = ":", invert = T),]
    write.csv(geneTE, paste0(path,"DESeq2/rmSampAc2Ctrl1/",mark, ".gene.teBackGround.csv"))
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/rmSampAc2Ctrl1/",mark, ".te.GeneBackGround.csv"))
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
    tmp2 <- tmp2[,c(1,3,5,6)]
    colLen <- ceiling(length(colnames(tmp2))/2)
    teGene <- tmp2
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/rmSampAc2Ctrl1/",mark, ".te.site.GeneBackGround.csv"))
  }
  shy_diff(mark = "MM-", 
           geneFile = paste0(path, "count/ccMM-.cntTable"),
           siteFile = paste0(path, "count/ccMM-Site.cntTable"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
}
#根据差异结果，展示维恩关系
{
  te01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyTE.MM.csv", row.names = 1)#不带基因背景
  te01 <- te01[grep(rownames(te01), pattern=":", invert = F),]
  te01$mark <- NA
  te01$mark[which(te01$log2FoldChange > 0 & te01$pvalue < 0.05)] <- "upup"
  te01$mark[which(te01$log2FoldChange < 0 & te01$pvalue < 0.05)] <- "down"
  #pvalue
  te02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.teOnly.csv", row.names = 1)
  te02$mark <- NA
  te02$mark[which(te02$log2FoldChange > 0 & te02$pvalue < 0.05)] <- "upup"
  te02$mark[which(te02$log2FoldChange < 0 & te02$pvalue < 0.05)] <- "down"
  #
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv", row.names = 1)#不带TE背景
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(te01[which(te01$mark=="upup"),])
  tmp1[["b"]] = rownames(te02[which(te02$mark=="down"),])
  tmp2 <- list()
  tmp2[["a"]] = rownames(ge01[which(ge01$mark=="upup"),])
  tmp2[["b"]] = rownames(ge02[which(ge02$mark=="down"),])
  #
  tmp3 <- list()
  tmp3[["a"]] = rownames(te01[which(te01$mark=="down"),])
  tmp3[["b"]] = rownames(te02[which(te02$mark=="upup"),])
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    main = "TE: upup KD_MM -- down OE_MM-",
                                    sub  = "upup: FC is 0 & pvalue<.05; down: FC is 0 & pvalue<.05", 
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  p_02 <- VennDiagram::venn.diagram(x = tmp2, 
                                    main = "Gene: upup KD_MM -- down OE_MM-",
                                    sub  = "upup: FC is .58 & padj<.05; down: FC is .58 & padj<.05", 
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  p_03 <- VennDiagram::venn.diagram(x = tmp3, 
                                    main = "TE: down KD_MM -- upup OE_MM-",
                                    sub  = "down: FC is 0 & pvalue<.05; upup: FC is 0 & pvalue<.05", 
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  p_04 <- VennDiagram::venn.diagram(x = tmp4,
                                    main = "Gene: down KD_MM -- upup OE_MM-",
                                    sub  = "down: FC is .58 & padj<.05; upup: FC is .58 & padj<.05",
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v1.r2.pdf"))
  grid.draw(p_01); dev.off()
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v2.r2.pdf"))
  grid.draw(p_02); dev.off()
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v3.r2.pdf"))
  grid.draw(p_03); dev.off()
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v4.r2.pdf"))
  grid.draw(p_04); dev.off()
  #
  base::intersect(tmp1[["a"]],tmp1[["b"]])
  base::intersect(tmp3[["a"]],tmp3[["b"]])
  base::intersect(tmp2[["a"]],tmp2[["b"]])
  base::intersect(tmp4[["a"]],tmp4[["b"]])
  print(length(base::intersect(tmp1[["a"]],tmp1[["b"]])))
  print(length(base::intersect(tmp2[["a"]],tmp2[["b"]])))
  print(length(base::intersect(tmp3[["a"]],tmp3[["b"]])))
  print(length(base::intersect(tmp4[["a"]],tmp4[["b"]])))
}
#交集的225个基因 [MM KD下调 -- MM- activation上调]: 和2-cell gene的Venn
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv", row.names = 1)#不带TE背景
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #
  twoc <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed")
  #
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  tmp4[["c"]] = twoc$V4
  #225 genes
  ne01 <- base::intersect(tmp4[["a"]],tmp4[["b"]])
  #
  p_04 <- VennDiagram::venn.diagram(x = tmp4,
                                    main = "KD_MM down -- OE_MM- upup -- 2Cell",
                                    sub  = "FC is 1.5 & padj<.05; 2-cell makers",
                                    category.names = c("KD_MM", "OE_MM-", "2-cell"),
                                    fill=c("#FD763F","#23BAC5","#EECA40"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.twoc.Venn.v1.pdf"))
  grid.draw(p_04)
  dev.off()
}
#交集的225个基因 [MM KD下调 -- MM- activation上调]: 和MM tags gene (rank >3)的Venn
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv", row.names = 1)#不带TE背景
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  #
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  tmp4[["c"]] = myta$name
  #225 genes
  ne01 <- base::intersect(tmp4[["a"]],tmp4[["b"]])
  #
  p_04 <- VennDiagram::venn.diagram(x = tmp4,
                                    main = "KD_MM down -- OE_MM- upup -- tags",
                                    sub  = "FC is 1.5 & padj<.05; MM-tags genes",
                                    category.names = c("KD_MM", "OE_MM-", "MM-tags rank>3"),
                                    fill=c("#FD763F","#23BAC5","#EECA40"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.tags.Venn.v1.pdf"))
  grid.draw(p_04)
  dev.off()
}
#交集的225个基因 [MM KD上调 -- MM- activation下调]: 和MM tags gene (rank >3)的Venn
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv", row.names = 1)#不带TE背景
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  #
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="upup"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="down"),])
  tmp4[["c"]] = myta$name
  #225 genes
  ne01 <- base::intersect(tmp4[["a"]],tmp4[["b"]])
  ne02 <- base::intersect(tmp4[["a"]],tmp4[["c"]])
  ne03 <- base::intersect(tmp4[["b"]],tmp4[["c"]])
  ovov <- base::intersect(ne01,tmp4[["c"]])
  setdiff(ne02, ne01)
  setdiff(ne03,ne01)
  #
  p_05 <- VennDiagram::venn.diagram(x = tmp4,
                                    main = "KD_MM upup -- OE_MM- down -- tags",
                                    sub  = "FC is 1.5 & padj<.05; MM-tags genes",
                                    category.names = c("KD_MM", "OE_MM-", "MM-tags rank>3"),
                                    fill=c("#FD763F","#23BAC5","#EECA40"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.tags.Venn.v1.pdf"))
  grid.draw(p_05)
  dev.off()
}
#交集的225个基因 [MM KD下调 -- MM- activation上调]: 和MM fusion gene的Venn
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv", row.names = 1)#不带TE背景
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #MM基因: 基于2-cell的组装转录本 [第一种方法]
  {
    gf01 <- import(inx3);gf01 <- gf01[which(gf01$type=="gene")] 
    tf01 <- import(inx4);tf01 <- tf01[which(tf01$gene_id %in% c("MM-int","MM-_Mm"))]
    tf02 <- import("/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/stringtie/star/two.gtf")
    tf02 <- tf02[which(tf02$type == "transcript")]
    ov01 <- suppressWarnings(findOverlaps(query = tf02, subject = tf01, ignore.strand = T))
    tf02 <- tf02[unique(queryHits(ov01))]
    ov02 <- suppressWarnings(findOverlaps(query = gf01, subject = tf02, ignore.strand = T))
    ne02 <- gf01[unique(queryHits(ov02))]$gene_name
  }
  #MM基因: 基于2-cell的组装转录本 [第二种方法]
  {
    #gf01 <- import(inx3);gf01 <- gf01[which(gf01$type=="gene")] 
    #tf01 <- import(inx4);tf01 <- tf01[which(tf01$gene_id %in% c("MM-int","MM-_Mm"))]
    #tf02 <- import("/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/stringtie/star/two.gtf")
    #
    #tf03 <- tf02[which(tf02$type == "transcript")]
    #tf04 <- tf02[which(tf02$type == "exon")]
    #ov01 <- suppressWarnings(findOverlaps(query = tf04, subject = tf01, ignore.strand = T))
    #tf04 <- tf04[unique(queryHits(ov01))]
    #tf03 <- tf03[which(tf03$transcript_id %in% unique(tf04$transcript_id))]
    #ov02 <- suppressWarnings(findOverlaps(query = gf01, subject = tf03, ignore.strand = T))
    #ne02 <- gf01[unique(queryHits(ov02))]$gene_name
  }
  #
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  tmp4[["c"]] = ne02
  #225 genes
  ne01 <- base::intersect(tmp4[["a"]],tmp4[["b"]])
  #
  p_04 <- VennDiagram::venn.diagram(x = tmp4,
                                    main = "KD_MM down -- OE_MM- upup -- MMFus",
                                    sub  = "FC is 1.5 & padj<.05; MM Fusion gene",
                                    category.names = c("KD_MM", "OE_MM-", "MM-Fusion"),
                                    fill=c("#FD763F","#23BAC5","#EECA40"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.MMFus.Venn.v1.pdf"))
  grid.draw(p_04)
  dev.off()
  intersect(intersect(tmp4[["a"]],tmp4[["b"]]),tmp4[["c"]])
}
#交集的225个基因 [MM KD下调 -- MM- activation上调]: 做富集
{
  ne01 <- base::intersect(tmp4[["a"]],tmp4[["b"]])
  geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ne01), keytype="SYMBOL", column="ENTREZID"))
  #
  kegg <- enrichKEGG(gene = geyy,
                     organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  kegg <- kegg@result
  kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
  #
  gobp <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db,
                                    keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
  gobp <- gobp@result
  #
  gomf <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db,
                                    keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
  gomf <- gomf@result
  #
  #
  kegg$Description <- factor(kegg$Description, levels = rev(kegg$Description))
  p_01 <- ggplot(kegg[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(8)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    labs(y="",title = "KEGG\n225 genes: MM- ac up & MM KD down") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  p_01
  #
  gobp$Description <- factor(gobp$Description, levels = rev(gobp$Description))
  p_02 <- ggplot(gobp[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(8)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    labs(y="",title = "GO BP\n225 genes: MM- ac up & MM KD down") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  p_02
  gomf$Description <- factor(gomf$Description, levels = rev(gomf$Description))
  p_03 <- ggplot(gomf[1:10,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(8)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
    labs(y="",title = "GO MF\n225 genes: MM- ac up & MM KD down") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  p_03
  plot <- p_01 +p_02 +p_03 +plot_layout(nrow = 3)
  ggsave(plot = plot, 
         units = "cm", width = 24, height = 24,
         filename = paste0(path,"PLOT/MM-ac.MMkd.kegg.gobp.gpmf.pdf"))
}
#根据维恩关系，画出上述交集的21个TE的热图 (其中1个是repeat, 所以去掉)
{
  myte <- base::intersect(tmp3[["a"]],tmp3[["b"]])
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("MMKD-1","MMKD-2","MM-KD-1","MM-KD-2","ctrl-1","ctrl-2")
  cnt1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]
  cnt1 <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
  cnt2 <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/count/ccMM-.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt2) <- c("MM-Ac-1","MM-Ac-2","MM-Ac-3","ctrlAc-1","ctrlAc-2","ctrlAc-3")
  cnt2 <- cnt2[grep(x = rownames(cnt2), pattern=":", invert = F),]
  cnt2 <- as.data.frame(apply(cnt2, 2, function(x){x/sum(x)*1000000}))
  cnts <- cbind(cnt1[,c(1,2,5,6)], cnt2[,c(1,3,5,6)])
  dta1 <- cnts[which(rownames(cnts) %in% myte),]
  dta2 <- dta1#dta2[,1:2] <- as.data.frame(apply(dta1[,1:2],2,function(x){-(1/(x/rowMeans(dta1[,c(3,4)])))}))
  dta2[,1:4] <- as.data.frame(apply(dta1[,1:4],2,function(x){x/rowMeans(dta1[,c(3,4)])}))
  dta2[,5:8] <- as.data.frame(apply(dta1[,5:8],2,function(x){x/rowMeans(dta1[,c(7,8)])}))
  pdf("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/PLOT/venn/pheatmap.21TE.r2.pdf", width = 8, height = 10)
  dta3 <- dta2[,c(3,4,1,2)]
  dta3$family <- sapply(str_split(rownames(dta3), pattern = ":"),"[",2)
  rank <- c("ERVL","ERVL-MaLR","ERV1","ERVK","L1")
  myta <- data.frame()
  for (i in rank) {
    tst1 <- dta3[which(dta3$family == i),]
    myta <- rbind(myta, tst1)
  }
  rownames(myta) <- sapply(str_split(rownames(myta), pattern = ":"),"[",1)
  annotation_row <- data.frame(family = myta$family)
  rownames(annotation_row) <- rownames(myta)
  annotation_col <- list(class = c("ERVL" = "#F7A6AC", "ERVL-MaLR"="#F7B7D2",
                                   "ERV1" = "#EEC186", "ERVK" = "#EEF0A7",
                                   "L1" = "#B2DBB9"))
  myta <- myta[c(1,4,2,3,5:20),-5]
  annotation_row$family <- factor(annotation_row$family, levels = unique(annotation_row$family))
  pheatmap(mat = myta,
           main="heatmap of 21 venn TEs (MM KD)",
           scale = "row", cellwidth = 20,
           annotation_row = annotation_row, annotation_names_row = F,
           annotation_colors = annotation_col,
           show_rownames = T, show_colnames = T, cluster_row = F, cluster_cols = F,
           labels_col = colnames(myta), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  
  #
  dta4 <- dta2[,c(7,8,5,6)]
  dta4$family <- sapply(str_split(rownames(dta4), pattern = ":"),"[",2)
  rank <- c("ERVL","ERVL-MaLR","ERV1","ERVK","L1")
  myta <- data.frame()
  for (i in rank) {
    tst1 <- dta4[which(dta4$family == i),]
    myta <- rbind(myta, tst1)
  }
  rownames(myta) <- sapply(str_split(rownames(myta), pattern = ":"),"[",1)
  annotation_row <- data.frame(family = myta$family)
  rownames(annotation_row) <- rownames(myta)
  annotation_col <- list(class = c("ERVL" = "#F7A6AC", "ERVL-MaLR"="#F7B7D2",
                                   "ERV1" = "#EEC186", "ERVK" = "#EEF0A7",
                                   "L1" = "#B2DBB9"))
  myta <- myta[c(1,4,2,3,5:20),-5]
  annotation_row$family <- factor(annotation_row$family, levels = unique(annotation_row$family))
  pheatmap(mat = myta,
           main="heatmap of 21 venn TEs (MM- Act)",
           scale = "row", cellwidth = 20,
           annotation_row = annotation_row, annotation_names_row = F,
           annotation_colors = annotation_col,
           show_rownames = T, show_colnames = T, cluster_row = F, cluster_cols = F,
           labels_col = colnames(myta), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
#根据维恩关系，画出上述交集的21个TE的热图 [在早期胚胎发育 -- 不分早中晚]
{
  myte <- base::intersect(tmp3[["a"]],tmp3[["b"]])
  base::load(file = "/Reference/aaSHY/DatabaseFiles/exprPlotHeatmap.RData")
  exprHeat <- exprHeat[,colnames(exprHeat)[grep(x = colnames(exprHeat), pattern = "G04_")]]
  colnames(exprHeat) <- gsub(x = colnames(exprHeat), pattern = "G04_", replacement = "")
  dta3 <- exprHeat[which(rownames(exprHeat) %in% myte),]
  #
  dta3$family <- sapply(str_split(rownames(dta3), pattern = ":"),"[",2)
  rank <- c("ERVL","ERVL-MaLR","ERV1","ERVK","L1")
  myta <- data.frame()
  for (i in rank) {
    tst1 <- dta3[which(dta3$family == i),]
    myta <- rbind(myta, tst1)
  }
  rownames(myta) <- sapply(str_split(rownames(myta), pattern = ":"),"[",1)
  annotation_row <- data.frame(family = myta$family)
  rownames(annotation_row) <- rownames(myta)
  annotation_col <- list(class = c("ERVL" = "#F7A6AC", "ERVL-MaLR"="#F7B7D2",
                                   "ERV1" = "#EEC186", "ERVK" = "#EEF0A7",
                                   "L1" = "#B2DBB9"))
  myta <- myta[c(1,4,2,3,5:20),-7]
  annotation_row$family <- factor(annotation_row$family, levels = unique(annotation_row$family))
  pdf("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/PLOT/venn/pheatmap.21TE.embryogenesis.pdf", width = 8, height = 10)
  pheatmap(mat = myta,
           main="heatmap of 21 venn TEs (MM- Act)",
           scale = "row", cellwidth = 20,
           annotation_row = annotation_row, annotation_names_row = F,
           annotation_colors = annotation_col,
           show_rownames = T, show_colnames = T, cluster_row = F, cluster_cols = F,
           labels_col = colnames(myta), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
#根据维恩关系，画出交集基因的KEGG结果 [失败]
#KEGG结果不太理想，所以尝试带上TE背景 [失败]
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/geneTE.MM.csv", row.names = 1)#不带TE背景
  ge01 <- ge01[grep(rownames(ge01), pattern=":", invert = T),]
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.gene.teBackGround.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #
  tmp2 <- list()
  tmp2[["a"]] = rownames(ge01[which(ge01$mark=="upup"),])
  tmp2[["b"]] = rownames(ge02[which(ge02$mark=="down"),])
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  kdupup <- base::intersect(tmp2[["a"]],tmp2[["b"]])#111
  kddown <- base::intersect(tmp4[["a"]],tmp4[["b"]])#217
  for (symb in seq(2)) {
    ovgene <- unlist(list(kdupup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    kegg <- enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    kegg <- kegg@result
    kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
    assign(paste0("kta",symb), kegg[which(kegg$pvalue<0.01),])
  }
  kta1
  kta2
}
#单独画GOBP：MM KD后下调的1058 genes, MM- activate后上调的950 genes [采用]
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/"
  kddown <- tmp4[["b"]]
  mtacup <- tmp4[["a"]]
  for (symb in seq(2)) {
    ovgene <- unlist(list(mtacup,kddown)[symb])
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(ovgene), keytype="SYMBOL", column="ENTREZID"))
    gobp <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db,
                                      keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
    gobp <- gobp@result
    assign(paste0("gta",symb), gobp[which(gobp$pvalue<0.01),])
  }
  #
  gta1$Description <- factor(gta1$Description, levels = rev(gta1$Description))
  p_05 <- ggplot(gta1[1:10,])+
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
                labs(y="",title = "950 genes: MM- ac up") +
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
          geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
          geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
          geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
          labs(y="",title = "1058 genes: MM kd down") +
          theme(text = element_text(size = 12,family = "serif",color = "black"),
                panel.border = element_rect(linewidth = 1),
                plot.title = element_text(size = 12,family = "serif",color = "#5861AC",face = "bold"),
                axis.title = element_text(size = 12,family = "serif",color = "#5861AC"),
                axis.text  = element_text(size = 12,family = "serif",color = "#5861AC"))
  p_06
  plot <- p_05 +p_06 +plot_layout(nrow = 2)
  ggsave(plot = plot, 
         units = "cm", width = 24, height = 15,
         filename = paste0(path,"PLOT/venn/","MM-ac.MMkd.gobp.pdf"))
}
#针对MM KD、MM- activate做GSEA
{
  inx6 <- "/Reference/aaSHY/GSEAgmt/TwoGenes.gmt"
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.csv")
  ge01 <- ge01[order(ge01$log2FoldChange, decreasing = T),]
  ge01 <- ge01[,c(1,3)]
  write.table(x = ge01,
              file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.rnk",
              quote = F, col.names = F, row.names = F,sep = "\t")
  paste0("gseapy prerank -r ",path,"/DESeq2/rmSampAc2Ctrl1/MM-.geneOnly.rnk ",
         "-g ",inx6," -p 10 -o ",path,"PLOT/EnrichmentGSEA/rmSampAc2Ctrl1")
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv")
  ge02 <- ge02[order(ge02$log2FoldChange, decreasing = T),]
  ge02 <- ge02[,c(1,3)]
  write.table(x = ge02,
              file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/MM.geneOnly.rnk",
              quote = F, col.names = F, row.names = F,sep = "\t")
  paste0("gseapy prerank -r ",
         "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/MM.geneOnly.rnk ",
         "-g ",inx6," -p 10 -o ",
         "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/PLOT/enrichGSEA")
}
#######################################################################
#20、把MM- activation上调/下调的基因/TE和MM KD下调/上调的基因/TE重合下 [KEGG]
{
  {
    resf <- c("/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/MM-.geneOnly.csv",
              "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/gene_MM.csv")
    dta1 <- data.frame()
    for (diff in resf) {
      for (changes in c("upup","down")){
        resdata <- as.data.frame(read.csv(diff, header = T, stringsAsFactors = F))
        resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
        if (changes=="upup"){gexx <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]}
        if (changes=="down"){gexx <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]}
        geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(gexx), keytype="SYMBOL", column="ENTREZID"))
        kegg <- enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
        kegg <- kegg@result
        kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
        kegg <- kegg[which(kegg$pvalue<0.01),]
        temp <- data.frame(chang = changes,
                           terms = kegg$Description,
                           marks = str_split(diff, pattern = "/")[[1]][7])
        dta1 <- rbind(dta1, temp)
      }
    }
  }
  #
  dta1$terms <- gsub(dta1$terms, pattern=",", replacement="_")
  kd01 <- dta1[which(dta1$marks=="MMKD"),]
  ac01 <- dta1[which(dta1$marks=="MM-Activation"),]
  tmp1 <- list()
  tmp1[["a"]] = kd01$terms[which(kd01$chang=="upup")]
  tmp1[["b"]] = ac01$terms[which(ac01$chang=="down")]
  tmp2 <- list()
  tmp2[["a"]] = kd01$terms[which(kd01$chang=="down")]
  tmp2[["b"]] = ac01$terms[which(ac01$chang=="upup")]
  p_01 <- VennDiagram::venn.diagram(x = tmp1, 
                                    main = "Gene: upup KD_MM -- down OE_MM-",
                                    sub  = "FC is 0 & padj<.05; enrich pvalue is .01", 
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FFFFCC","#CCFFFF"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1.5, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  p_02 <- VennDiagram::venn.diagram(x = tmp2, 
                                    main = "Gene: down KD_MM -- upup OE_MM-",
                                    sub  = "FC is 0 & padj<.05; enrich pvalue is .01",
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FFFFCC","#CCFFFF"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 1.5, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v1.pdf"))
  grid.draw(p_01); dev.off()
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v2.pdf"))
  grid.draw(p_02); dev.off()
  #
  base::intersect(tmp1[["a"]],tmp1[["b"]])
  base::intersect(tmp2[["a"]],tmp2[["b"]])
  write.csv(x = dta1, 
            file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/dumpROOM/enrich.kd.ov.ac.csv",
            row.names = F, quote = F)
}
#######################################################################
#21、
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/gene_MM.csv", row.names = 1)
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/DESeq2/MM-.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  #
  tmp2 <- list()
  tmp2[["a"]] = rownames(ge01[which(ge01$mark=="upup"),])
  tmp2[["b"]] = rownames(ge02[which(ge02$mark=="down"),])
  #
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  #
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/MM-Activation/"
  p_02 <- VennDiagram::venn.diagram(x = tmp2, 
                                    main = "Gene: upup KD_MM -- down OE_MM-",
                                    sub  = "upup: FC is .58 & padj<.05; down: FC is .58 & padj<.05", 
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)
  #
  p_04 <- VennDiagram::venn.diagram(x = tmp4,
                                    main = "Gene: down KD_MM -- upup OE_MM-",
                                    sub  = "down: FC is .58 & padj<.05; upup: FC is .58 & padj<.05",
                                    category.names = c("KD_MM", "OE_MM-"), fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v2.pdf"))
  grid.draw(p_02); dev.off()
  pdf(file = paste0(path, "PLOT/MM-.MM.Venn.v4.pdf"))
  grid.draw(p_04); dev.off()
}
#######################################################################
#22、可变剪接
{
  bam1 <- list.files(paste0(path,"bam"), pattern = "ac.*bam$",full.names = T)
  bam2 <- list.files(paste0(path,"bam"), pattern = "co.*bam$",full.names = T)
  cat(file = paste0(path,"rMATS/opt.txt"), sep = ",", bam1)
  cat(file = paste0(path,"rMATS/ctr.txt"), sep = ",", bam2)
  paste0("rmats.py --nthread 10 -t paired --readLength 140 ",
         "--gtf ",inx3," --b1 ",path,"rMATS/opt.txt ","--b2 ",path,
         "rMATS/ctr.txt"," --od ",path,"rMATS --tmp ",path,"rMATS/tmp","\n")
}
################################################################################
#23、MM- ac在45S rDNA上的富集
{
  shell <- paste0(path,"src/run_suppl.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- basename(list.files(path = paste0(path, "fq"), pattern = "gz", recursive = T))
  prex <- unique(gsub(pattern = "_R.*.fq.gz", replacement = "", fixed = F, x = prex))
  lay1 <- paste0(path,"trim/",prex,"_R1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_R2_val_2.fq.gz")
  for (kk in c("MM","rRNA")) {
    if (kk == "MM") {
      inxM <- inxD
    }
    else {
      inxM <- inxC
    }
    cmd_04 <- paste0("STAR --genomeDir ",inxM, " --genomeLoad LoadAndExit","\n")
    cmd_05 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                     "--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 ",
                     "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                     "--outFileNamePrefix ",paste0(path,"consensus/",kk,"/",prex)," ",
                     "--outSAMprimaryFlag AllBestScore ",
                     "--genomeDir ",inxM," --readFilesIn ",lay1," ",lay2,"\n")
    cmd_06 <- paste0("mv ",path,"consensus/",kk,"/",prex,"Aligned.sortedByCoord.out.bam ",path,
                     "consensus/",kk,"/",prex,".bam","\n")
    cmd_07 <- paste0("samtools index ",paste0(path,"consensus/",kk,"/",prex,".bam"),"\n")
    cmd_08 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                     "consensus/",kk,"/",prex,".bam -o ",path,
                     "consensus/",kk,"/",prex,".bw","\n")
    for (i in cmd_04) {cat(i, file = shell, append = T)}
    for (i in cmd_05) {cat(i, file = shell, append = T)}
    for (i in cmd_06) {cat(i, file = shell, append = T)}
    for (i in cmd_07) {cat(i, file = shell, append = T)}
    for (i in cmd_08) {cat(i, file = shell, append = T)}
  }
  rm(kk)
  print(paste0("nohup bash ",shell," >",path,"Log/run_suppl.log"," 2>&1 &"))
}

