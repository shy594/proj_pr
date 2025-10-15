#write by ~~~ at 2024.04.11
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
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
  cmd_01 <- paste0("trim_galore --phred33 -j 8 ",
                   "-o ", path, "trim --paired ",lay1," ",lay2,"\n")
  lay1 <- paste0(path,"trim/",prex,"_R1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_R2_val_2.fq.gz")
  cmd_03 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_04 <- paste0("mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_05 <- paste0("samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_07 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project gfa --outdir ",paste0(path,"count"),"\n")
  cmd_08 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_09 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_10 <- paste0("stringtie --merge -o ",
                   path,"stringtie/gfa.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_11 <- paste0("cat ", paste0(path,"stringtie/gfa.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/gfa.merge.R.gtf"),"\n")
  cmd_12 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/gfa.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project gfa ",
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
  #
  #two-pass mode是否有影响
  cmd_13 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"dumpROOM/",prex," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_14 <- paste0("mv ",
                   paste0(path,"dumpROOM/",prex,"Aligned.sortedByCoord.out.bam "),
                   paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_15 <- paste0("samtools index ",paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_16 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"dumpROOM/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"dumpROOM/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project gfa --outdir ",paste0(path,"dumpROOM"),"\n")
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  for (i in cmd_16) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
#######################################################################
#0A、带和不带2-pass mode, 对差异倍数的影响
{
  tmp1 <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/count/gfa_sigdiff_gene_TE.txt"
  tmp1 <- read.table(tmp1)
  tmp2 <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/dumpROOM/gfa_sigdiff_gene_TE.txt"
  tmp2 <- read.table(tmp2)
  tmp3 <- cbind(tmp1[,], tmp2[match(rownames(tmp1),rownames(tmp2)),])
  tmp3 <- tmp3[,c(2,8)]
  tmp3$mark <- tmp3$log2FoldChange/tmp3$log2FoldChange.1
  tmp3 <- tmp3[order(tmp3$mark, decreasing = T),]
  tmp3[1:10,]
  ggplot(tmp3) +geom_point(aes(x = log2FoldChange, y = log2FoldChange.1)) +geom_abline(intercept = 0)
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
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
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
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.teOnly.down.csv"))
  }
  shy_diff(mark = "gfaOE", 
           geneFile = paste0(path, "count/gfa.cntTable"),
           siteFile = paste0(path, "count/gfa.site.cntTable"),
           Class = factor(c("OE","OE","WT","WT")), VS = c("Class","OE","WT"))
}
#######################################################################
#05、CPM Excel
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/count/gfa.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("gfaOE-1","gfaOE-2","pcag-1","pcag-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/count/gfaOE.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/count/gfaOE.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/count/gfaOE.cpm.allGeneTE.xlsx")
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
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/gfaOE.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/gfaOE.geneOnly.volcano.pdf"))
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
  shy_volcano_te(resFile  = paste0(path, "DESeq2/gfaOE.te.GeneBackGround.csv"),
                 saveFile = paste0(path, "PLOT/gfaOE.te.GeneBackGround.volcano.pdf"))
  shy_volcano_te(resFile  = paste0(path, "DESeq2/gfaOE.teOnly.csv"),
                 saveFile = paste0(path, "PLOT/gfaOE.teOnly.volcano.pdf"))
}
#######################################################################
#08、MA
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
  shy_teMA(difffile = paste0(path,"DESeq2/gfaOE.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/gfaOE.te.GeneBackGround.MA.pdf"))
}
#######################################################################
#09、MA plot TE sites
{
  shy_repSite <- function(resdFile, saveFile){
    gf00 <- read.csv(file = resdFile, header = T, stringsAsFactors = F, row.names = 1)
    gf00 <- na.omit(gf00)
    gfaa <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*"),]
    gfbb <- gf00[grep(rownames(gf00), pattern = ".*MT2_Mm.*"),]
    gfcc <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*|.*MT2_Mm.*", invert = T),]
    #开始作图
    p_1 <- ggplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 0.8),
            axis.line = element_line(), axis.text = element_text(family = "serif")) +
      ylab(label = expression(Log2Fold~Change)) +
      xlab(label = expression(Log2~(baseMean + 1))) +
      scale_color_manual(breaks = c("MERVL-int","MT2_Mm","Others"),values = c("red","purple","gray")) +
      geom_point(data=gfcc,aes(x=log2(baseMean+1),y=log2FoldChange,colour="Others"),size=1.2) +
      geom_point(data=gfbb,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MT2_Mm"),size=1.2) +
      geom_point(data=gfaa,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MERVL-int"),size=1.2) +
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
  shy_repSite(resdFile = paste0(path,"DESeq2/gfaOE.te.site.teOnly.csv"),
              saveFile = paste0(path,"PLOT/gfaOE.te.site.teOnly.mervl.pdf"))
}
#######################################################################
#10、PCA
{
  #----------------------------------PCA-batch----------------------------------
  cnt1 <- read.table(paste0(path, "count/gfa.cntTable"),
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  cnt1 <- cnt1[grep(rownames(cnt1), pattern = ":", invert = T),]
  colnames(cnt1) <- sapply(str_split(basename(colnames(cnt1)), ".bam"),"[",1)
  #
  base::load(file = "/Reference/aaSHY/DatabaseFiles/forPCA.RData")
  myta <- cbind(cnt1, forPCA)
  #
  condition = sapply(str_split(colnames(myta),"-"), "[", 1)
  btch = c(rep("ours",length(colnames(cnt1))), condition[c(length(colnames(cnt1))+1):length(myta)])
  clas <- data.frame(row.names = colnames(myta), condition = condition, batch = btch)
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
#11、pheatmap
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
    ######pheatmap(mat = tete[which(rownames(tete) %in% sigT_down),],
    ######         main="the heatmap of Downregulated TEs",
    ######         scale = "row", cellwidth = 60, 
    ######         show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
    ######         labels_col = colnames(tete), angle_col = "45",
    ######         breaks = seq(-2,2,by=0.01), border = F,
    ######         color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
    ######                    colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
    dev.off()
  }
  shy_pheatmap(mark = "heat", readsFile = paste0(path, "count/gfa.cntTable"),
               gene_resFile = paste0(path, "DESeq2/gfaOE.geneOnly.csv"),
               tete_resFile = paste0(path, "DESeq2/gfaOE.teOnly.csv"),
               savePath = paste0(path, "PLOT/pheatmap/"))
}
#######################################################################
#12、Up TEs loci number [TElocal]
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
      theme(legend.position = c(0.7,0.8)) +
      ylab("Number of Loci") + 
      xlab(NULL) +
      labs(title = expression(deRegulated~TEs~"in"~italic(gfa)^{"+/+"}~ESCs)) +
      scale_fill_manual(breaks = unique(forP$clas),values = ggsci::pal_aaas(palette = "default")(length(unique(forP$clas))))
    ggsave(plot = p_1, filename = saveFile, units = "cm", width = 12, height = 12)
  }
  shy_lociNum(resdFile = paste0(path, "DESeq2/gfaOE.te.site.teOnly.upup.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/gfa.lociNum.teOnly.upup.pdf"))
  shy_lociNum(resdFile = paste0(path, "DESeq2/gfaOE.te.site.teOnly.down.csv"),
              saveFile = paste0(path, "PLOT/TElociNum/gfa.lociNum.teOnly.down.pdf"))
}
#######################################################################
#13、MERVL fusion genes expression
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
      xlab(label = expression(Log2~(Fold~Change~of~italic(MT2)^{"+/+"}~ESCs))) +
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
           filename = paste0(path, "PLOT/chimeric/gfa.chimeric.",whichTE,".pdf"))
  }
  shy_fusion_stringTie(countFile = paste0(path, "stringtie/count/gfa.cntTable"), 
                       Class = c("OE","OE","WT","WT"),
                       contrast = c("Class","OE","WT"),
                       gtfFile = paste0(path, "stringtie/gfa.merge.R.gtf"),
                       whichTE = "MT2_Mm")#MT2_Mm#MERVL-int)
}
#######################################################################
#14、enrich KEGG
{
  resFile <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/DESeq2/gfaOE.geneOnly.csv"
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
  p_05 <- ggplot(kta1[1:7,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(8)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
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
         filename = paste0(path,"PLOT/enrich/","gfa.kegg.pdf"))
}
#######################################################################
#15、enrich GO BP, GO MF
{
  resFile <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202404/gfaOE/DESeq2/gfaOE.geneOnly.csv"
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
    #geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
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
  p_06 <- ggplot(gta2[1:9,])+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#5861AC", alpha=1.0) +
    theme_few() +
    scale_x_continuous(breaks = seq(10)) +
    geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    #geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
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
         filename = paste0(path,"PLOT/enrich/","gfa.gomf.pdf"))
}
#######################################################################
#16、GSEA
{
  shy_gsea <- function(resFile, savePath){
    aa <- read.csv(file = resFile, header = T, stringsAsFactors = F)[,c(1,3)]
    ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL",
                                                 toType = "ENTREZID",
                                                 OrgDb = "org.Mm.eg.db", drop = T))
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
  shy_gsea(resFile = paste0(path, "DESeq2/gfaOE.geneOnly.csv"), 
           savePath = paste0(path, "PLOT/enrichGSEA/"))
}
#######################################################################
#17、python GSEA
{
  temp <- read.csv(paste0(path, "DESeq2/gfaOE.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, file = paste0(path, "DESeq2/gfa.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/gfa.res.rnk"), " -p 10 -o ", 
                   paste0(path, "PLOT/enrichGSEA/"))#conda activate base gseapy 0.10.8
}
#######################################################################
#18、UCSC Track
#https://genome.ucsc.edu/s/hysun088/gfaOE
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
#19、Correlation
#保留了参数: #"bamCompare -p 10 --scaleFactorsMethod SES --sampleLength 100000 "
{
  shell <- paste0(path,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  bw01 <- list.files(paste0(path,"bamCompare"), pattern = "*bw$", full.names = T)
  bw02 <- list.files(path = "/ChIP_seq_1/aaSHY/likeTwoC",
                     recursive = T, pattern = ".bw$", full.names = T)
  bw02 <- bw02[grep(bw02, pattern="bamCompare")]
  bw03 <- list.files(path = "/ChIP_seq_2/aaSHY/pr/ask/MERVL/mt2Activation/bamCompare", full.names = T)
  cmd_01 <- paste0("multiBigwigSummary bins ", 
                   "-p 20 -b ",paste0(c(bw01,bw02,bw03), collapse = " "),
                   " --smartLabels",
                   " -o ",paste0(path,"PLOT/Correlation/gfa.cor.npz","\n"))
  cmd_02 <- paste0("plotCorrelation -in ",path,
                   "PLOT/Correlation/gfa.cor.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", path,"PLOT/Correlation/gfa.cor.pdf ",
                   "--outFileCorMatrix ",path,"PLOT/Correlation/gfa.cor.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
#######################################################################
#20、哪些是滋养层基因
{
  base::load(file = "/Reference/aaSHY/DatabaseFiles/GO_DATA.RData")
  #GO ID和GO Name
  tmp2 <- data.frame(goID = base::names(GO_DATA$PATHID2NAME),goName = base::unname(GO_DATA$PATHID2NAME))
  tmp2 <- tmp2[which(tmp2$goID %in% tmp1$goID[which(tmp1$goONT=="BP")]),]
  #GO ID和GO Genes Name
  tmp3 <- GO_DATA$PATHID2EXTID
  #
  myid <- tmp2[grep(x = tmp2$goName, pattern = "trophecto", ignore.case = T),]
  set1 <- unique(unname(unlist(tmp3[myid$goID])))
  #吕老师说有些不是滋养层基因
  #
  data <- read.table(paste0(path,"count/gfa.cntTable"), 
                     header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = ".bam"), "[", 1)
  gene <- data[grep(rownames(data), pattern = ":", invert = T),]
  gene <- gene[rowSums(gene) > 5,]
  gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
  gene <- gene[which(rownames(gene) %in% set1),]
  pheatmap(mat = gene,
           main="the heatmap of the Genes",
           scale = "row", cellwidth = 60, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = T,
           labels_col = colnames(gene), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
#######################################################################
#21、对gfa OE后激活的2细胞基因排个序
{
  gata <- read.csv(file = paste0(path,"DESeq2/gfaOE.geneOnly.csv"),
                   sep = ",", header = T, row.names = 1, stringsAsFactors = F)
  twoc <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed", sep = "\t", header = F, stringsAsFactors = F)
  wr01 <- gata[which(rownames(gata) %in% twoc$V4),]
  wr01 <- wr01[order(wr01$log2FoldChange, decreasing = T),]
  wr01$geneName <- rownames(wr01)
  openxlsx::write.xlsx(x = wr01, 
                       file = paste0(path,"DESeq2/gfaOE.geneOnly.2C.gene.xlsx"), asTable = T)
}
#######################################################################
#22、可变剪接
{
  bam1 <- list.files(paste0(path,"bam"), pattern = "ga.*bam$",full.names = T)
  bam2 <- list.files(paste0(path,"bam"), pattern = "pc.*bam$",full.names = T)
  cat(file = paste0(path,"rMATS/opt.txt"), sep = ",", bam1)
  cat(file = paste0(path,"rMATS/ctr.txt"), sep = ",", bam2)
  paste0("rmats.py --nthread 10 -t paired --readLength 140 ",
         "--gtf ",inx3," --b1 ",path,"rMATS/opt.txt ","--b2 ",path,
         "rMATS/ctr.txt"," --od ",path,"rMATS --tmp ",path,"rMATS/tmp","\n")
}




