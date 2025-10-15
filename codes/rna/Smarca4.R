#write by Sun Haiayng at 2024.05.17
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/publicData/Smarca4/"
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
  prex <- c("KD-1","KD-2","WT-1","WT-2")
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   list.files(path = paste0(path, "rawdata"), pattern = "sra", recursive = T, full.names = T),
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  cmd_03 <- paste0("trim_galore --phred33 -j 4 --clip_R1 4 ",
                   "-o ", path, "trim ", path, "fq/", prex, ".fq.gz", "\n")
  cmd_04 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1, " --readFilesIn ", path,
                   "trim/",prex,"_trimmed.fq.gz","\n")
  cmd_05 <- paste0("mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_06 <- paste0("samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_07 <- paste0("bamCoverage --normalizeUsing CPM -p 10 ",
                   "--minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   path,"bam/",prex,".bam"," -o ",path,"bamCoverage/",prex,".bw","\n")
  cmd_08 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project smarca4 --outdir ",paste0(path,"count"),"\n")
  cmd_09 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prex,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/",prex),"\n")
  cmd_10 <- paste0("stringtie ", path,"bam/",prex,".bam"," ",
                   "-p 20 -G ", inx3, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_11 <- paste0("stringtie --merge -o ",
                   path,"stringtie/smarca4.merge.gtf -G ",inx3," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_12 <- paste0("cat ", paste0(path,"stringtie/smarca4.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/smarca4.merge.R.gtf"),"\n")
  cmd_13 <- paste0("TEtranscripts ",
                   "-t ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[3:4],".bam"),collapse = " ")," ",
                   "--GTF ",path,"stringtie/smarca4.merge.R.gtf ",
                   "--TE ", inx4," --sortByPos --mode multi --project smarca4 ",
                   "--outdir ",path,"stringtie/count", "\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
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
  for (i in cmd_13) {cat(i, file = shell, append = T)}
  #
  #two-pass mode是否有影响
  cmd_14 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",path,"dumpROOM/",prex," ",
                   "--outSAMprimaryFlag AllBestScore --twopassMode Basic ",
                   "--genomeDir ",inx1," --readFilesIn ",path,
                   "trim/",prex,"_trimmed.fq.gz","\n")
  cmd_15 <- paste0("mv ",
                   paste0(path,"dumpROOM/",prex,"Aligned.sortedByCoord.out.bam "),
                   paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_16 <- paste0("samtools index ",paste0(path,"dumpROOM/",prex,".bam"),"\n")
  cmd_17 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"dumpROOM/",prex[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"dumpROOM/",prex[3:4],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ssrp1 --outdir ",paste0(path,"dumpROOM"),"\n")
  for (i in cmd_14) {cat(i, file = shell, append = T)}
  for (i in cmd_15) {cat(i, file = shell, append = T)}
  for (i in cmd_16) {cat(i, file = shell, append = T)}
  for (i in cmd_17) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log"," 2>&1 &"))
}
################################################################################
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
  shy_diff(mark = "smarca4-KD", 
           geneFile = paste0(path, "count/smarca4.cntTable"),
           siteFile = paste0(path, "count/smarca4.site.cntTable"),
           Class = factor(c("KD","KD","WT","WT")), VS = c("Class","KD","WT"))
}
################################################################################
#05、CPM Excel
{
  cnt1 <- read.table(paste0(path,"count/","smarca4",".cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("KD-1","KD-2","WT-1","WT-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path,"count/","smarca4-KD",".cpm.onlyTE.xlsx"))
  openxlsx::write.xlsx(x = tmp2, file = paste0(path,"count/","smarca4-KD",".cpm.onlyGene.xlsx"))
  openxlsx::write.xlsx(x = tmp3, file = paste0(path,"count/","smarca4-KD",".cpm.allGeneTE.xlsx"))
}
################################################################################
#06、MA
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
  shy_teMA(difffile = paste0(path,"DESeq2/smarca4-KD.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/smarca4-KD.te.GeneBackGround.MA.pdf"))
  shy_teMA(difffile = paste0(path,"DESeq2/smarca4-KD.teOnly.csv"),
           savefile = paste0(path,"PLOT/smarca4-KD.teOnly.MA.pdf"))
}
################################################################################
#07、部分基因的表达柱状图
{
  cnt1 <- read.table(paste0(path,"count/","smarca4",".cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("KD-1","KD-2","WT-1","WT-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  #
  #先画基因only Gene
  {
    dta1 <- tmp2[which(rownames(tmp2) %in% c("Duxf3","AW822073","Gm4981","Gm19459","Tcstv3")),]
    dta1$KD <- rowMeans(dta1[,1:2])/rowMeans(dta1[,3:4])
    dta1$WT <- rowMeans(dta1[,3:4])/rowMeans(dta1[,3:4])
    dta1$KDerr <- apply(dta1[,1:2], 1, function(x){sd(x)/sqrt(2)})
    dta1$WTerr <- apply(dta1[,3:4], 1, function(x){sd(x)/sqrt(2)})
    dta1 <- dta1[,-c(1:4)]
    dta1$name <- rownames(dta1)
    myta <- data.frame()
    for (i in dta1$name) {
      ssss <- data.frame(name = i,
                         class = c("KD","WT"),
                         cpm = as.numeric(dta1[which(dta1$name==i),1:2]),
                         err = sd(as.numeric(dta1[which(dta1$name==i),1:2]))/sqrt(2))
      myta <- rbind(myta, ssss)
    }
    myta[1,]
    myta$class <- factor(myta$class, levels = c("WT","KD"))
    ggplot(data = myta) +
      geom_bar(aes(x = name, y = cpm, fill = class), linewidth =1, color = "white", stat = "identity", position = "dodge") +
      theme_classic() +
      scale_fill_manual(values = c("gray30","steelblue4")) +
      geom_errorbar(aes(x = name, ymin = cpm, ymax = cpm+err, color = class),
                    stat = "identity", position = "dodge", width=0.2)
  }
  #
  #再画MM-int
  dta1 <- tmp3[grep(x = rownames(tmp3), pattern = "MM-int|MT2_Mm"),]
  dta1$KD <- rowMeans(dta1[,1:2])/rowMeans(dta1[,3:4])
  dta1$WT <- rowMeans(dta1[,3:4])/rowMeans(dta1[,3:4])
  dta1$KDerr <- apply(dta1[,1:2], 1, function(x){sd(x)/sqrt(2)})
  dta1$WTerr <- apply(dta1[,3:4], 1, function(x){sd(x)/sqrt(2)})
  dta1 <- dta1[,-c(1:4)]
  dta1$name <- rownames(dta1)
  myta <- data.frame()
  for (i in dta1$name) {
    ssss <- data.frame(name = i,
                       class = c("KD","WT"),
                       cpm = as.numeric(dta1[which(dta1$name==i),1:2]),
                       err = sd(as.numeric(dta1[which(dta1$name==i),1:2]))/sqrt(2))
    myta <- rbind(myta, ssss)
  }
  myta[1,]
  myta$class <- factor(myta$class, levels = c("WT","KD"))
  p_01 <- ggplot(data = myta) +
    geom_bar(aes(x = name, y = cpm, fill = class), linewidth =1, color = "white", stat = "identity", position = "dodge") +
    theme_classic() +
    scale_fill_manual(values = c("gray30","steelblue4")) +
    theme(legend.position = c(0.9,0.9), 
          legend.title = element_text(family = "serif", size = 15),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = NA),
          legend.spacing.x = unit(0,"mm"),
          legend.spacing.y = unit(0,"mm"),legend.key = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(family = "serif", size = 16),
          axis.text = element_text(family = "serif", size = 12)) +
    guides(fill = guide_legend(override.aes = list(size = 3)))
  ggsave(plot = p_01, 
         "/ChIP_seq_2/aaSHY/pr/publicData/Smarca4/PLOT/expr_MM.pdf", units = "cm", width = 18, height = 16)
}





