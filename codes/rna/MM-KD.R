#write by ~~~ at 2024.01.10
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","ggridges","ggrepel","KEGG.db","scales","stringr",
              "Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork",
              "ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager",
              "dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler",
              "topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/STARv279a"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/STARv279a"
}
#03、跑命令
{
  for (i in c("src","fq","Log","stringtie","DESeq2","bamCoverage","trim","dumpROOM","bam","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_b.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- basename(list.files(path = paste0(path, "fq"), pattern = "gz", recursive = T))
  prex <- unique(gsub(pattern = "_.*.fq.gz", replacement = "", fixed = F, x = prex))
  lay1 <- paste0(path,"fq/",prex,"_R1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_R2.fq.gz")
  cmd_03 <- paste0("trim_galore --phred33 -j 8 ",
                   "-o ", path, "trim --paired ", lay1, " ", lay2, "\n")
  lay1 <- paste0(path,"trim/",prex,"_R1_val_1.fq.gz")
  lay2 <- paste0(path,"trim/",prex,"_R2_val_2.fq.gz")
  cmd_04 <- paste0("STAR --genomeDir ",inx1, " --genomeLoad LoadAndExit","\n")
  cmd_05 <- paste0("STAR --runThreadN 10 --readFilesCommand zcat ",
                   "--genomeLoad LoadAndKeep --limitBAMsortRAM 30000000000 ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore ",
                   "--genomeDir ",inx1," --readFilesIn ",lay1," ",lay2,"\n")
  cmd_06 <- paste0("mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_07 <- paste0("samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_08 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 --samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".bam")," -o ",paste0(path,"bamCoverage/",prex,".bw"),"\n")
  cmd_09 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prex[3:6],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " "),
                   " --GTF ",inx2," --TE ",inx3," --sortByPos --mode multi ",
                   "--project ccMMKD --outdir ",paste0(path,"count"),"\n")
  cmd_10 <- paste0("TElocal -b ",path,"bam/",prex,".bam ",
                   "--GTF ",inx2," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ",path, "count/",prex,"\n")
  cmd_11 <- paste0("stringtie ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx2, " -o ", paste0(path,"stringtie/",prex,".gtf"),"\n")
  cmd_12 <- paste0("stringtie --merge -o ",
                   paste0(path,"stringtie/MMKD.merge.gtf")," -G ",inx2," ",
                   paste(paste0(path,"stringtie/",prex,".gtf"), collapse = " "),"\n")
  cmd_13 <- paste0("cat ", paste0(path,"stringtie/MMKD.merge.gtf "), 
                   "|sed 's/gene_name/gene_abcd/g; s/\\tgene_id/\\tgene_name/g' >",
                   paste0(path,"stringtie/MMKD.merge.R.gtf"),"\n")
  cmd_14 <- paste0("TEtranscripts ",
                   "-t ", paste(paste0(path,"bam/",prex[3:6],".bam"),collapse = " ")," ",
                   "-c ",paste(paste0(path,"bam/",prex[1:2],".bam"),collapse = " ")," ",
                   "--GTF ",paste0(path,"stringtie/MMKD.merge.R.gtf "),
                   "--TE ", inx3," --sortByPos --mode multi --project MMMerge ",
                   "--outdir ", paste0(path,"stringtie/count"),"\n")
  #
  #
  #
  cmd_aa <- paste0("TElocal -b ",path,"bam/",prex,".bam ",
                   "--GTF ",inx2," --TE ",inx4," --sortByPos --mode uniq ",
                   "--project ",path, "countUniq/",prex,"\n")
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_10) {cat(i, append = T, file = shell)}
  for (i in cmd_11) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  for (i in cmd_aa) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#04、差异分析：Gene、TE、TE sites (DESeq2)
{
  data <- read.table(paste0(path, "count/MMKD.cntTable"), header = T, check.names = F, stringsAsFactors = F)
  gene <- data[grep(data[,1], pattern = ":", invert = T),]
  te01 <- data[grep(data[,1], pattern = ":", invert = F),]
  gete <- data
  site <- read.table(paste0(path, "count/MMKD.site.cntTable"), header = T, check.names = F, stringsAsFactors = F)
  site <- data[grep(data[,1], pattern = ":", invert = F),]
  shy_diff <- function(mydata, marker, Class){
    counts <- mydata
    rownames(counts) <- counts[,1]
    counts <- counts[,-1]
    colnames(counts) <- sapply(str_split(basename(colnames(counts)), pattern = ".bam"), "[", 1)
    for (i in c("MM","mt2")) {
      sm01 <- which(c("MM","mt2")==i)*2
      sm01 <- c(sm01-1, sm01, 5,6)
      inputs <- counts[,sm01]
      colLen <- ceiling(length(colnames(inputs))/2)
      yyyCOL <- data.frame(row.names = colnames(inputs), Class, stringsAsFactors = T)
      ddsdds <- DESeq2::DESeqDataSetFromMatrix(countData = inputs, colData = yyyCOL, design = ~Class)
      #ddsdds <- ddsdds[rowSums(BiocGenerics::counts(ddsdds) > 2) >= colLen, ]
      ddsdds <- ddsdds[rowMeans(BiocGenerics::counts(ddsdds)) > 5,]
      ddsdds <- DESeq2::DESeq(ddsdds)
      resres <- DESeq2::results(ddsdds, contrast = c("Class","KD","WT"))
      write.csv(resres, paste0(path,"DESeq2/",marker,".",i,".csv"))
      upup <- subset(resres, log2FoldChange >  log2(1.5) & padj<0.05)
      down <- subset(resres, log2FoldChange < -log2(1.5) & padj<0.05)
      write.csv(upup, paste0(path,"DESeq2/",marker,".",i,".upup.csv"))
      write.csv(down, paste0(path,"DESeq2/",marker,".",i,".down.csv"))
    }
  }
  shy_diff(mydata = gene, Class = factor(c("KD","KD","WT","WT")),marker = "onlyGene")
  shy_diff(mydata = te01, Class = factor(c("KD","KD","WT","WT")),marker = "onlyTE")
  shy_diff(mydata = gete, Class = factor(c("KD","KD","WT","WT")),marker = "geneTE")
  shy_diff(mydata = site, Class = factor(c("KD","KD","WT","WT")),marker = "site.onlyTE")
}
################################################################################
#05、计算CPM
{
  data <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/count/ccMMKD.cntTable",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  colnames(data) <- sapply(str_split(basename(colnames(data)), ".bam"),"[",1)
  #放一起做标准化
  myta <- as.data.frame(apply(data, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- myta[,c(length(colnames(myta)),1:c(length(colnames(myta))-1))]
  openxlsx::write.xlsx(x = myta, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/count/MM-KD.CPM.xlsx")
  #gene单独标准化
  myta <- data[grep(rownames(data), pattern = ":", invert = T),]
  myta <- as.data.frame(apply(myta, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- myta[,c(length(colnames(myta)),1:c(length(colnames(myta))-1))]
  openxlsx::write.xlsx(x = myta, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/count/MM-KD.gene.CPM.xlsx")
  #tete单独标准化
  myta <- data[grep(rownames(data), pattern = ":", invert = F),]
  myta <- as.data.frame(apply(myta, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- myta[,c(length(colnames(myta)),1:c(length(colnames(myta))-1))]
  openxlsx::write.xlsx(x = myta, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/count/MM-KD.tete.CPM.xlsx")
}
################################################################################
#06、说MT2的一个样本和Control搞错了，画出PCA
{
  #----------------------------------PCA-batch----------------------------------
  data <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/count/ccMMKD.cntTable",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  colnames(data) <- sapply(str_split(basename(colnames(data)), ".bam"),"[",1)
  data <- data[grep(rownames(data), pattern = ":", invert = T),]
  tabl <- data.frame(names = colnames(data),
                     condition = sapply(str_split(colnames(data),"_"), "[", 1), 
                     batch = c(rep("A",2), rep("B",2), rep("C",2)))
  countmatrix <- as.matrix(data)
  dds1 <- suppressWarnings(DESeq2::DESeqDataSetFromMatrix(countmatrix, colData = tabl, design = ~condition))#[rowSums(counts(vtax))>5]
  ##https://cloud.tencent.com/developer/article/1625223
  dds2 <- BiocGenerics::estimateSizeFactors(dds1)
  #rawx <- SummarizedExperiment(counts(dds2, normalized = F), colData=colData(dds2))
  #norx <- SummarizedExperiment(counts(dds2, normalized = T), colData=colData(dds2))
  rldf <- rlog(dds2)#vsdf <- vst(dds2), 豆包说“通常情况下，数据集小于30个样品时可以用rlog，数据集大于30个样品时用vst。”
  ##去除批次效应
  #View(mat)
  #mat <- assay(rldf)
  #mat <- limma::removeBatchEffect(mat, rldf$batch)
  #assay(rldf) <- mat
  #
  #pd_1 <- plotPCA(DESeqTransform(rawx), intgroup = c("condition"), returnData = T)
  #pd_2 <- plotPCA(DESeqTransform(norx), intgroup = c("condition"), returnData = T)
  pd_3 <- plotPCA(rldf, intgroup = c("condition"), returnData = T)#pd_4 <- plotPCA(vsdf, intgroup = c("condition"), returnData = T)
  shy_linshi <- function(df, i){
    mark = c("raw","normlization","rlog")[i]
    perc <- round(100 * attr(df, "percentVar"))
    df$condition <- factor(as.character(df$condition), 
                           levels = as.character(df$condition),
                           labels = as.character(df$condition))
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
            legend.position = c(0.5,0.25),
            axis.text    = element_text(family = "serif", size = 12),
            axis.title   = element_text(family = "serif", size = 12))
    print(ppp)
  }
  pdf(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202306/PLOT/PCA/PCA.pdf",width = 10, height = 8)
  #shy_linshi(df = pd_1, i = 1)
  #shy_linshi(df = pd_2, i = 2)
  shy_linshi(df = pd_3, i = 3)#shy_linshi(df = pd_4, i = 4)
  dev.off()
}
#热图
{
  path="/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/"
  data <- read.table(paste0(path, "count/ccMMKD.cntTable"), 
                     header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(data) <- c("MMKD-1","MMKD-2","mt2KD-1","mt2KD-2","ctrl-1","ctrl-2"); #colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = "_"), "[", 1)
  gene <- data[grep(rownames(data), pattern = ":", invert = T),]
  gene <- gene[rowSums(gene) > 5,]
  gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
  #
  #pdf(paste0(path,"heatMap/deRegulateTE.pdf"), width = 12, height = 10)
  pheatmap(mat = gene,
           main="\n\nthe heatmap of all genes",
           scale = "row", cellwidth = 60, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = T,
           labels_col = colnames(gene), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  #dev.off()
}
################################################################################
#05、追加早期胚胎发育阶段TE的表达数据
{
  ff01 <- read.csv(file = paste0(path, "DESeq2/tete_MM.csv"))
  ff02 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/tete.cpmAndSD.xlsx")
  ff03 <- ff02[,c("geneName",colnames(ff02)[-1][c(6,12,8,9,11,10,5,4,7,1,3,2)])]
  ff01[,colnames(ff02)[-1][c(6,12,8,9,11,10,5,4,7,1,3,2)]] <- NA
  for (i in seq_along(ff01$X)) {
    ff01[i,8:19] <- ff03[which(ff03$geneName==ff01$X[i]),-1]
  }
  openxlsx::write.xlsx(x = ff01, file = paste0(path, "DESeq2/tete_MMKD_andStages.xlsx"))
}
################################################################################
#06、可视化部分TE的表达量变化
{
  #01、用2014 Deng的数据
  cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
  prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
  colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
  cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = F),]#invert = F, 那么输出TE的CPM
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  df03 <- data.frame()
  for (i in c("ETnERV3-int", "RLTR13B2", "RLTR19D", "L1M8", "L1MEb", "MM-int")) {
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i, ":")),]
    df03 <- rbind(df03, temp)
  }
  df03$stage <- factor(x = df03$stage, levels = prex)
  p_01 <- ggplot(data = df03) +
          geom_line(mapping = aes(x = stage, y = cpmExp, group = geneName), color = "steelblue", linewidth = 1) +
          geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 1.5) +
          geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                        width = 0.15, color = "steelblue", linewidth = 0.8) +
          labs(x = "", y = "mean CPM", title = "2014 Deng Early Mid Late") +
          facet_wrap(.~ geneName, scales = "free_y", nrow = 4) +
          theme_classic() +
          theme(plot.margin = margin(1,1,1,1,unit = "cm"),
                plot.title = element_text(size = 16, color = "black", face = "bold"),
                axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
                strip.text.x = element_text(size = 12, color = "red", face = "bold"),
                strip.background = element_blank(),
                panel.spacing.y = unit(0.8, "cm"))
  p_01
  #02、用pacbio corresponding illumina的数据
  stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  colnames(stag) <- sapply(str_split(colnames(stag), ".bam."), "[", 2)
  stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))
  stag <- stag[grep(x = rownames(stag), pattern = ":", invert = F), -7]
  my01 <- data.frame()
  for (i in c("ETnERV3-int", "RLTR13B2", "RLTR19D", "L1M8", "L1MEb", "MM-int")) {
    temp <- stag[grep(rownames(stag), pattern = paste0("^", i, ":")),]
    my01 <- rbind(my01, temp)
  }
  my01$name <- rownames(my01)
  my01 <- tidyr::gather(my01, key = "stage", value = "expres", -"name")
  my01$stage <- factor(x = my01$stage, levels = unique(my01$stage))
  p_02 <- ggplot(data = my01) + 
          geom_line(mapping = aes(x = stage, y = expres, group = name)) +
          geom_point(mapping = aes(x = stage, y = expres)) +
          labs(x="", y="CPM Expression",
               title = "PacBio corresponding Illumina data") +
          facet_wrap(.~ name, scales = "free_y") +
          theme(plot.margin = margin(1,1,1,1,unit = "cm"),
                plot.title = element_text(size = 16, color = "black", face = "bold"),
                axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
                strip.text.x = element_text(size = 12, color = "red", face = "bold"))
  p_02
  ##
  ##
  ##
  #03、把MM-int和RLTR1B-int画在一张子图里
  #####用2014 Deng的数据
  cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
  prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
  colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
  cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = F),]#invert = F, 那么输出TE的CPM
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  df03 <- data.frame()
  for (i in c("RLTR1B-int","MM-int")) {
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i, ":")),]
    df03 <- rbind(df03, temp)
  }
  df03 <- df03[which(!(df03$stage %in% c("twoEarly","twoMid","twoLate"))),]
  df03$stage <- factor(x = df03$stage, levels = prex)
  p_03 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = log2(cpmExp), group = geneName, color = geneName), linewidth = 1) +
    scale_color_aaas() +
    geom_point(mapping = aes(x = stage, y = log2(cpmExp)), color = "black", size = 1.5) +
    geom_errorbar(mapping = aes(x = stage, ymin = log2(cpmExp-cpmSD), ymax = log2(cpmExp+cpmSD)), width = 0.15, color = "black", linewidth = 0.8) +
    labs(x = "", y = "log2(mean CPM)", title = "2014 Deng") +
    theme_classic() +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          axis.text.y = element_text(family = "serif", size = 12, colour = "black"),
          axis.title = element_text(family = "serif", size = 12, colour = "black"),
          legend.text = element_text(family = "serif", size = 12, colour = "black"),
          legend.title = element_text(family = "serif", size = 12, colour = "black"),
          panel.spacing.y = unit(0.8, "cm"))
  p_03
  ggsave(plot = p_03, 
         filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/dumpROOM/MMRLTR-cpm-data1.pdf",
         units = "cm", height = 15, width = 20)
  #####用pacbio corresponding illumina的数据
  stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  colnames(stag) <- sapply(str_split(colnames(stag), ".bam."), "[", 2)
  stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))
  stag <- stag[grep(x = rownames(stag), pattern = ":", invert = F), -7]
  my01 <- data.frame()
  for (i in c("RLTR1B-int","MM-int")) {
    temp <- stag[grep(rownames(stag), pattern = paste0("^", i, ":")),]
    my01 <- rbind(my01, temp)
  }
  my01$name <- rownames(my01)
  my01 <- tidyr::gather(my01, key = "stage", value = "expres", -"name")
  my01$stage <- factor(x = my01$stage, levels = unique(my01$stage))
  p_04 <- ggplot(data = my01) + 
    geom_line(mapping = aes(x = stage, y = log2(expres), group = name, color=name), linewidth = 1) +
    scale_color_aaas() +
    geom_point(mapping = aes(x = stage, y = log2(expres)), size=1.5) +
    labs(x="", y="log2(CPM Expression)",
         title = "PacBio corresponding Illumina data") +
    theme_classic() +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          axis.text.y = element_text(family = "serif", size = 12, colour = "black"),
          axis.title = element_text(family = "serif", size = 12, colour = "black"),
          legend.text = element_text(family = "serif", size = 12, colour = "black"),
          legend.title = element_text(family = "serif", size = 12, colour = "black"),
          panel.spacing.y = unit(0.8, "cm"))
  p_04
  ggsave(plot = p_04, 
         filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/dumpROOM/MMRLTR-cpm-data2.pdf",
         units = "cm", height = 15, width = 20)
  ##
  ##
  ##
  my02 <- my01[which(my01$name=="RLTR1B-int:ERV1:LTR"),]
  my02 <- my01[which(my01$name=="MM-int:ERVL:LTR"),]
  p_05 <- ggplot(data = my02) + 
    geom_line(mapping = aes(x = stage, y = expres, group = name, color=name), linewidth = 1) +
    scale_color_aaas() +
    geom_point(mapping = aes(x = stage, y = expres), size=1.5) +
    labs(x="", y="log2(CPM Expression)",
         title = "PacBio corresponding Illumina data") +
    theme_classic() +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          legend.position = c(0.3, 0.8),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          axis.text.y = element_text(family = "serif", size = 12, colour = "black"),
          axis.title = element_text(family = "serif", size = 12, colour = "black"),
          legend.text = element_text(family = "serif", size = 12, colour = "black"),
          legend.title = element_text(family = "serif", size = 12, colour = "black"),
          panel.spacing.y = unit(0.8, "cm"))
  p_05
  ggsave(plot = p_05,
         filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/dumpROOM/RLTR-cpm-data2.pdf",
         units = "cm", height = 15, width = 16)
}
################################################################################
#07、MM KD上、下调的TE和pr上、下调的TE取交集
{
  ff01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.fullte.csv", row.names = 1)
  ff01$mark <- NA
  ff01$mark[which(ff01$log2FoldChange > log2(1.5) & ff01$pvalue < 0.05)] <- "upup"
  ff01$mark[which(ff01$log2FoldChange < -log2(1.5) & ff01$pvalue < 0.05)] <- "down"
  #
  ff02 <- read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/tete_MM.csv", row.names = 1)
  ff02$mark <- NA
  ff02$mark[which(ff02$log2FoldChange > log2(1.5) & ff02$pvalue < 0.05)] <- "upup"
  ff02$mark[which(ff02$log2FoldChange < -log2(1.5) & ff02$pvalue < 0.05)] <- "down"
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(ff01[which(ff01$mark=="upup"),])
  tmp1[["b"]] = rownames(ff02[which(ff02$mark=="upup"),])
  tmp2 <- list()
  tmp2[["a"]] = rownames(ff01[which(ff01$mark=="down"),])
  tmp2[["b"]] = rownames(ff02[which(ff02$mark=="down"),])
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "TE Overlap",
                                    sub  = "pr KO74 upup & MM KD upup", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr", "MM KD"), fill=c("#FFFFCC","#CCFFFF"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "Plot/MMGeneVenn.pdf"))
  grid.draw(p_01)
  dev.off()
  base::intersect(tmp1[["a"]],tmp1[["b"]])
  base::intersect(tmp2[["a"]],tmp2[["b"]])
}
################################################################################
#08、MM KD下调的基因和pr上调的基因取交集
{
  ff01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO92.gene.csv", row.names = 1)
  ff01$mark <- NA
  ff01$mark[which(ff01$log2FoldChange > log2(1.5) & ff01$pvalue < 0.05)] <- "upup"
  ff01$mark[which(ff01$log2FoldChange < -log2(1.5) & ff01$pvalue < 0.05)] <- "down"
  #
  ff02 <- read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/gene_MM.csv", row.names = 1)
  ff02$mark <- NA
  ff02$mark[which(ff02$log2FoldChange > log2(1.5) & ff02$pvalue < 0.05)] <- "upup"
  ff02$mark[which(ff02$log2FoldChange < -log2(1.5) & ff02$pvalue < 0.05)] <- "down"
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(ff01[which(ff01$mark=="upup"),])
  tmp1[["b"]] = rownames(ff02[which(ff02$mark=="down"),])
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "Gene Overlap",
                                    sub  = "pr KO92 upup & MM KD down", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr", "MM KD"), fill=c("#FFFFCC","#CCFFFF"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "Plot/MMGeneVenn.pdf"))
  grid.draw(p_01)
  dev.off()
  base::intersect(tmp1[["a"]],tmp1[["b"]])
  base::intersect(tmp2[["a"]],tmp2[["b"]])
}
################################################################################
#09、火山图
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
      scale_color_manual(values = c("red","blue","black"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      scale_x_continuous(limits = c(-8,8), breaks = seq(-8,8,by=1)) +
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
  shy_volcano_gene(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/gene_MM.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT/MMKD.gene.volcano.pdf")
}
################################################################################
#10、TE MA plot
{
  resdata = na.omit(read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/tete_MM.csv", header = TRUE, stringsAsFactors = F,row.names=1))
  rownames(resdata) <- sapply(str_split(rownames(resdata), pattern = ":"),"[", 1)
  resdata["MM-int",]
  get_MA_plot<-function(path,relation){
    annot<-read.delim(relation,header=T,stringsAsFactors=F)
    rownames(annot)<-annot$repName
    terc<-as.data.frame(resdata)
    terc<-terc[rownames(terc)%in%annot$repName,]
    annot<-annot[rownames(terc),]
    terc<-cbind(terc,annot)
    terc_s<-terc[which((terc$pvalue<0.05&(abs(terc$log2FoldChange)>log2(1.5)))&(terc$repClass=="LTR"|terc$repClass=="LINE"|terc$repClass=="SINE"|terc$repClass=="DNA")),]
    terc_o<-terc[which(!((terc$pvalue<0.05&(abs(terc$log2FoldChange)>log2(1.5)))&(terc$repClass=="LTR"|terc$repClass=="LINE"|terc$repClass=="SINE"|terc$repClass=="DNA"))),]
    xm=0
    for (i in terc_s$repFamily) {
      xm=xm+1
      if (i %in% c("ERV1","ERVK","ERVL","Gypsy")){
        terc_s$P[xm]<-i}
      else {
        terc_s$P[xm]<-"Others"
      }
    }
    terc_o$P<-rep("Others",dim(terc_o)[1])
    terc_a<-rbind(terc_o,terc_s)
    write.csv(terc_a, "as.csv")
    library(ggplot2)
    library(ggrepel)
    idx<-terc_a$repName%in%terc_s$repName
    terc_a$repName[!idx]<-NA
    terc_a$P<-factor(terc_a$P,levels=c("ERV1","ERVK","ERVL","Gypsy","Others"))
    color_s<-list(ERV1="sienna1",ERVK="steelblue3",ERVL="palevioletred2",
                  Gypsy="green",Others="grey80")
    color_s<-unlist(color_s)
    
    ggplot(terc_a,aes(x=log2(baseMean+1),y=terc_a$log2FoldChange,color=P, 
                      shape=NULL))+
      geom_point(size=2.5)+
      geom_text_repel(aes(label=repName),hjust=-0.125,size=4,colour="black",family="serif")+
      geom_hline(aes(yintercept= log2(1.5)),linetype="dashed",colour="grey")+
      geom_hline(aes(yintercept= -log2(1.5)),linetype="dashed",colour="grey")+
      scale_color_manual(values = color_s)+
      guides(color=guide_legend(title = "repFamliy"))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.line = element_line(),axis.text = element_text(family = "serif"),
            panel.border = element_rect(colour = "black", linewidth = 0.8))+
      theme(legend.position = c(0.91,0.2),legend.background = element_blank(),
            legend.text = element_text(family = "serif")) +
      xlab("log2(Expression in baseMean)")+ylab("log2Foldchange")
    ggsave(filename = paste0(path,"/fullTEsMAplot_padj.pdf"), height = 4.16, width = 4.7)
  }
  get_MA_plot(
    path="/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT",
    relation="/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/relation_fullTEs.txt"
  )
}
################################################################################
#11、附加上调TE在早中晚的表达图，不卡Fold Change
#图基于2014Deng、2020PacBioMateNGS的数据
{
  myta = na.omit(read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/tete_MM.csv", header = TRUE, stringsAsFactors = F,row.names=1))
  myta <- myta[which(myta$pvalue < 0.05 & myta$log2FoldChange >0),]
  #
  base::load(file = "/Reference/aaSHY/DatabaseFiles/exprPlotBar.RData")
  pd01 <- df01[which(df01$geneName %in% rownames(myta)),]
  pd01$stage <- factor(x = pd01$stage, levels = unique(pd01$stage))
  p_01 <- ggplot(data = pd01) +
    geom_bar(mapping = aes(x = stage, y = cpmExp), stat = "identity", fill = "steelblue", width = 0.8) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group=1)) +
    geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                  width = 0.2, color = "black", linewidth = 0.8) +
    labs(x = "", y = "CPM", title = "2014Deng") +
    facet_wrap(.~ geneName, scales = "free_y", nrow = 4) +
    theme_classic() +
    theme(plot.margin = margin(0,0,0,0,unit = "cm"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  ggsave(plot = p_01, 
         filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT/2014Deng_on_8UpTE.pdf", 
         units = "cm", height = 22, width = 18)
  #p_01
  ##
  pd04 <- df04[which(df04$geneName %in% rownames(myta)),]
  pd04$sample <- factor(x = pd04$sample, levels = unique(pd04$sample))
  p_04 <- ggplot(data = pd04) +
    geom_bar(mapping = aes(x = sample, y = cpm), stat = "identity", fill = "steelblue", width = 0.8) +
    geom_line(mapping = aes(x = sample, y = cpm, group=1)) +
    labs(x = "", y = "CPM", title = "2020PacBioMateNGS") +
    facet_wrap(.~ geneName, scales = "free_y", nrow = 4) +
    theme_classic() +
    theme(plot.margin = margin(0,0,0,0,unit = "cm"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  #p_04
  ggsave(plot = p_04, 
         filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT/2020PacBioMateNGS_on_8UpTE.pdf", 
         units = "cm", height = 22, width = 18)
}
################################################################################
#11、pvalue < 0.05 & log2FoldChange >0条件下共有8个TE；画下这几个ERV上H3K9me3的富集
{
  myta = na.omit(read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/tete_MM.csv", header = TRUE, stringsAsFactors = F,row.names=1))
  myta <- myta[which(myta$pvalue < 0.05 & myta$log2FoldChange >0),]
  te01 <- gsub(x = rownames(myta), pattern = ":", replacement = "::")
  bed2 <- paste0("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/", te01, ".bed", collapse = " ")
  #
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/"
  shell <- paste0(path,"src/run_onMMKD.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", bed2, " -S ",
                   paste0(path, "/bw/bcp/IP-1_Input.bw ",path, "/bw/bcp/IP-3_Input.bw "),
                   "-o ",path,"PLOT/MMKD_upTE.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red --perGroup -m ",path,
                   "PLOT/MMKD_upTE.mat.gz ",
                   "--samplesLabel ", "IP-1_Input IP-3_Input ",
                   "-out ",path,"PLOT/MMKD_upTE.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_onMMKD.log"," 2>&1 &"))
}
################################################################################
#12、热图pheatmap
{
  path="/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/"
  data <- read.table("/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cntTable", 
                     header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(data) <- c("MMKD-1","MMKD-2","mt2KD-1","mt2KD-2","ctrl-1","ctrl-2"); #colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = "_"), "[", 1)
  tete <- data[grep(rownames(data), pattern = ":", invert = F),c(1,2,5,6)]
  tete <- tete[rowSums(tete) > 5,]
  tete <- as.data.frame(apply(tete, 2, function(x){x/sum(x)*1000000}))
  myta <- na.omit(read.csv("/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/tete_MM.csv", header = TRUE, stringsAsFactors = F,row.names=1))
  upup <- rownames(myta[which(myta$pvalue < 0.05 & myta$log2FoldChange >0),])
  upup <- tete[upup,]
  down <- rownames(myta[which(myta$pvalue < 0.05 & myta$log2FoldChange <0),])
  down <- tete[down,]
  pdf(paste0(path,"heatMap/deRegulateTE.pdf"), width = 12, height = 10)
  pheatmap(mat = upup,
           main="\n\nthe heatmap of upup TEs",
           scale = "row", cellwidth = 60, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(upup), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  pheatmap(mat = down,
           main="\n\nthe heatmap of down TEs",
           scale = "row", cellwidth = 60, 
           show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(down), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
#13、MM KD on 45s rRNA
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
################################################################################
#14、GSEA
{
  aa <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/gene_MM.csv", header = T, stringsAsFactors = F)[,c(1,3)]
  ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = T))
  ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ae <- ad$rank; names(ae) <- ad$ENTREZID
  aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt")
  aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
  ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
  af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab, pvalueCutoff = 1))
  #dev.new()
  pdf(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT/GSEA_on_2C_marker.pdf", width = 15, height = 10)
  enrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for 2-cell Genes", pvalue_table = T)
  dev.off()
}
#
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/"
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
  }
  shy_gsea(resFile = paste0(path, "DESeq2/onlyGene.MM.csv"), 
           savePath = paste0(path, "PLOT/enrichGSEA/"))
}
################################################################################
#15、GO BP
#气泡图
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
           units = "cm", width = 30, height = 16,
           filename = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".bubble.pdf"))
  }
  shy_fuji_bubble(resFile  = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/gene_MM.csv",
                  savePath = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/PLOT/",
                  change  = "upup")
  shy_fuji_bubble(resFile  = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/gene_MM.csv",
                  savePath = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/PLOT/",
                  change  = "down")
}
#柱状图
{
  shy_fuji_bar <- function(resFile, change, savePath){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- resdata[grep(resdata[,1], pattern = ":", invert = T),]
    if (change=="upup"){gexx <- resdata[which(resdata$log2FoldChange > log2(1.5)  &  resdata$padj < 0.05),"X"]}
    if (change=="down"){gexx <- resdata[which(resdata$log2FoldChange < -log2(1.5) &  resdata$padj < 0.05),"X"]}
    geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(gexx), keytype="SYMBOL", column="ENTREZID"))
    ptha <- enrichKEGG(gene = geyy, organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    ptha <- DOSE::setReadable(ptha, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
    if (change=="upup"){kegg <- ptha@result[1:7,][-2,]}# for upup
    if (change=="down"){kegg <- ptha@result[1:6,]}# for down
    kegg$Description <- gsub(x = kegg$Description, pattern = " - Mus musculus (house mouse)", fixed = T, replacement = "")
    wr01 <- ptha@result
    write.csv(x = wr01, file = paste0(savePath, change, ".kegg.",str_split(basename(resFile),".csv")[[1]][1],".csv"))
    #
    #pthb <- ReactomePA::enrichPathway(gene = geyy, organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    #reac <- pthb@result[1:15,]
    #wr02 <- reac@result
    #write.csv(x = wr02, file = paste0(savePath, change, ".reac.",str_split(basename(resFile),".csv")[[1]][1],".csv"))
    #
    #gogo <- enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
    #term <- gogo@result[1:5,]
    #wr03 <- gogo@result
    #write.csv(x = wr03, file = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".csv"))
    ##画图
    kegg$Description <- factor(kegg$Description, levels = rev(rbind(kegg$Description)))
    #reac$Description <- factor(reac$Description, levels = rev(rbind(reac$Description)))
    #term$Description <- factor(term$Description, levels = rev(rbind(term$Description)))
    p_1 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
      geom_bar(aes(), fill="steelblue", stat = "identity") +#geom_point(aes(color=-log10(pvalue),size=Count)) +
      theme_classic() +
      ylab(NULL) +
      #scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif",color = "black"),
            axis.text  = element_text(size = 12,family = "serif",color = "black"))
    #p_2 <- ggplot(reac, aes(x=-log10(pvalue), y=Description))+
    #  geom_bar(aes(size=Count), fill="steelblue", stat = "identity") +
    #  theme_classic() +
    #  ylab(NULL) +
    #  #scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
    #  theme(text = element_text(size = 12,family = "serif",color = "black"),
    #        axis.title = element_text(size = 12,family = "serif",color = "black"),
    #        axis.text  = element_text(size = 12,family = "serif",color = "black"))
    #p_3 <- ggplot(term, aes(x=-log10(pvalue), y=Description))+
    #  geom_bar(aes(size=Count),fill="steelblue", stat = "identity") +
    #  theme_classic() +
    #  ylab(NULL) +
    #  #scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
    #  theme(text = element_text(size = 12,family = "serif",color = "black"),
    #        axis.title = element_text(size = 12,family = "serif",color = "black"),
    #        axis.text  = element_text(size = 12,family = "serif",color = "black"))
    ggsave(plot = p_1, 
           units = "cm", width = 30, height = 8,
           filename = paste0(savePath, change, ".KEGG.",str_split(basename(resFile),".csv")[[1]][1],".bar.pdf"))
    #ggsave(plot = p_2, 
    #       units = "cm", width = 20, height = 16,
    #       filename = paste0(savePath, change, ".ReactomePA.",str_split(basename(resFile),".csv")[[1]][1],".bar.pdf"))
    #ggsave(plot = p_3, 
    #       units = "cm", width = 30, height = 8,
    #       filename = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".bar.pdf"))
  }
  shy_fuji_bar(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/gene_MM.csv",
               savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT/",
               change  = "upup")
  shy_fuji_bar(resFile  = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/gene_MM.csv",
               savePath = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/PLOT/",
               change  = "down")
}
################################################################################
#22、可变剪接
{
  path = "/ChIP_seq_1/aaSHY/likeTwoC/Zscan4c/"
  bam1 <- list.files(paste0(path,"bam"), pattern = "OE.*bam$",full.names = T)
  bam2 <- list.files(paste0(path,"bam"), pattern = "WT.*bam$",full.names = T)
  cat(file = paste0(path,"rMATS/opt.txt"), sep = ",", bam1)
  cat(file = paste0(path,"rMATS/ctr.txt"), sep = ",", bam2)
  paste0("rmats.py --nthread 10 -t paired --readLength 140 ",
         "--gtf ",inx3," --b1 ",path,"rMATS/opt.txt ","--b2 ",path,
         "rMATS/ctr.txt ",
         "--od ",path,"rMATS/ --tmp ",path,
         "rMATS/tmp","\n")
}
################################################################################
#23、MM KD前后的Dux的reads数
{
  #reads数量
  cnt1 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.cntTable", 
                     sep = "\t", header = T, row.names = 1)
  colnames(cnt1) <- sapply(str_split(colnames(cnt1), pattern = "\\."),"[",8)
  cnt1 <- cnt1[,c(5,6,1,2)]
  cnt1[which(rownames(cnt1) %in% c("Duxf3","AW822073","Gm4981","Gm19459","Gm10807")),]
  #gene CPM
  cnt2 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.onlyGene.xlsx")
  cnt2 <- cnt2[which(cnt2$theName %in% c("Duxf3","AW822073","Gm4981","Gm19459","Gm10807")),c(7,5,6,1,2)]
  cnt2 <- tidyr::gather(data = cnt2, key = "class", value = "cpm", -"theName")
  ggplot(data = cnt2) +
    geom_bar(aes(x = theName, y = cpm, fill = class), stat = "identity", color="white", position = "dodge") +
    theme_classic() +
    scale_fill_manual(values = c("gray","gray","black","black")) +
    labs(x = "", y = "CPM") +
    theme(axis.title = element_text(size = 12, family = "serif"),
          axis.text  = element_text(size = 12, family = "serif"),
          strip.text = element_blank()) +
    facet_wrap(facets = ~theName, scales = "free")
}
################################################################################
#24、MT2 Activation前后的Dux的reads数
{
  #reads数量
  cnt1 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/count/ccMT2.cntTable",
                     sep = "\t", header = T, row.names = 1)
  colnames(cnt1) <- c("AC-1","AC-2","AC-3","WT-1","WT-2","WT-3")
  cnt1 <- cnt1[,c(5,6,1,3)]
  cnt1[which(rownames(cnt1) %in% c("Duxf3","AW822073","Gm4981","Gm19459","Gm10807")),]
  #gene CPM
  cnt2 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/count/mt2Ac.onlyGene.xlsx")
  cnt2 <- cnt2[which(cnt2$theName %in% c("Duxf3","AW822073","Gm4981","Gm19459","Gm10807")),c(7,5,6,1,3)]
  cnt2 <- tidyr::gather(data = cnt2, key = "class", value = "cpm", -"theName")
  ggplot(data = cnt2) +
    geom_bar(aes(x = theName, y = cpm, fill = class), stat = "identity", color="white", position = "dodge") +
    theme_classic() +
    scale_fill_manual(values = c("gray","gray","black","black")) +
    labs(x = "", y = "CPM") +
    theme(axis.title = element_text(size = 12, family = "serif"),
          axis.text  = element_text(size = 12, family = "serif"),
          strip.text = element_blank()) +
    facet_wrap(facets = ~theName, scales = "free")
}
################################################################################
#25、MM KD、MT2 AC RNA-seq在MM融合基因上的表达 [基于reads count]
#再次定义MM融合基因, 为了生成GTF文件
{
  gf01 <- import(inx2); gf01 <- gf01[which(gf01$type=="gene")]
  tf01 <- import(inx3); tf01 <- tf01[which(tf01$gene_id %in% c("MM-int","MT2_Mm"))]
  tf02 <- import("/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/stringtie/star/two.gtf")
  #
  tf03 <- tf02[which(tf02$type == "transcript")]
  tf04 <- tf02[which(tf02$type == "exon")]
  ov01 <- suppressWarnings(findOverlaps(query = tf04, subject = tf01, ignore.strand = T))
  tf04 <- tf04[unique(queryHits(ov01))]
  tf03 <- tf03[which(tf03$transcript_id %in% unique(tf04$transcript_id))]
  #
  ne01 <- as.data.frame(tf02[which(tf02$transcript_id %in% tf03$transcript_id)])
  ne02 <- as.data.frame(tf03)
  ne02$type <- "gene"
  df01 <- rbind(ne01,ne02)[,1:11]
  wr01 <- dplyr::arrange(.data = df01, gene_id, transcript_id, 
                         factor(x = df01$type, levels = c("gene","transcript","exon")))
  wr01$gene_name <- wr01$transcript_id
  wr01 <- makeGRangesFromDataFrame(df = wr01, keep.extra.columns = T)
  rtracklayer::export(object = wr01, format = "gtf",
                      con = "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.gtf")
}
{
  shell <- paste0(path,"src/run_countFus.sh")
  cat("#!/bin/bash\n", file = shell)
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/"
  cmd_01 <- paste0("TEtranscripts -t ",
                   path,"bam/MM_1.bam ",
                   path,"bam/MM_2.bam -c ",
                   path,"bam/Control_1.bam ",
                   path,"bam/Control_2.bam ",
                   "--GTF /Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.gtf ",
                   "--TE ",inx3," --sortByPos --mode multi ",
                   "--project MMKD-Fus --outdir ",paste0(path,"count"),"\n")
  bath <- "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/"
  cmd_02 <- paste0("TEtranscripts -t ",
                   bath,"bam/activation-1.bam ",
                   bath,"bam/activation-3.bam -c ",
                   bath,"bam/control-2.bam ",
                   bath,"bam/control-3.bam ",
                   "--GTF /Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.gtf ",
                   "--TE ",inx3," --sortByPos --mode multi ",
                   "--project mt2AC-Fus --outdir ",paste0(bath,"count"),"\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
#boxplot [MT2 Activation]
{
  cnt1 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/count/mt2AC-Fus.cntTable",
                     sep = "\t", header = T, row.names = 1)
  colnames(cnt1) <- c("AC-1","AC-2","WT-1","WT-2")
  cnt1 <- cnt1[grep(x = rownames(cnt1), pattern = ":", invert = T),]
  myta <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- tidyr::gather(data = myta, key = "class", value = "cpm", -"geneName")
  ggplot(data = myta) +
    geom_boxplot(aes(x = class, y = log2(cpm +0.1), fill = class)) +
    theme_classic() +
    scale_fill_manual(values = c("grey","grey","black","black")) +
    labs(x = "", y = "log2 (CPM +0.1)")
}
#boxplot [MM KD]
{
  cnt1 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD-Fus.cntTable",
                     sep = "\t", header = T, row.names = 1)
  colnames(cnt1) <- c("KD-1","KD-2","WT-1","WT-2")
  cnt1 <- cnt1[grep(x = rownames(cnt1), pattern = ":", invert = T),]
  myta <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- tidyr::gather(data = myta, key = "class", value = "cpm", -"geneName")
  ggplot(data = myta) +
    geom_boxplot(aes(x = class, y = log2(cpm +0.1), fill = class)) +
    theme_classic() +
    scale_fill_manual(values = c("grey","grey","black","black")) +
    labs(x = "", y = "log2 (CPM +0.1)")
}
################################################################################
#26、MM KD、MT2 AC RNA-seq在MM融合基因上的表达 [基于bamCoverage]
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/"
  bw01 <- list.files(path = path, pattern = ".bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bamCoverage")][c(1:4,7,9,11,12)]
  {
    bed1 <- bed1 <- "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/KD-AC-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/KD-AC-",basename(inx0),".mat.gz ",
                       "--samplesLabel KD-WT-1 KD-WT-2 KD-1 KD-2 AC-1 AC-2 AC-WT-1 AC-WT-2 ",
                       "-out ",path,"zempROOM/KD-AC-",basename(inx0),".r1.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MM/"
  bw01 <- list.files(path = path, pattern = ".bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bamCoverage")][c(1:4,7,9,11,12)]
  {
    bed1 <- bed1 <- "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (bwbw in bw01) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                       "--missingDataAsZero -R ", beds, " -S ",bwbw," ",
                       "-o ",path,"zempROOM/KD-AC-", basename(bwbw),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/KD-AC-",basename(bwbw),".mat.gz ",
                       "-out ",path,"zempROOM/KD-AC-",basename(bwbw),".r3.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#27、MM KD、MT2 AC RNA-seq中1165个基因的表达热图 [设WT CPM为1]
{
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  ne01 <- myta$name
  mervKD <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.onlyGene.xlsx")
  mervKD <- mervKD[which(mervKD$theName %in% ne01),]
  rownames(mervKD) <- mervKD$theName; mervKD <- mervKD[,c(1,2,5,6)]
  mervAC <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/count/mt2Ac.onlyGene.xlsx")
  mervAC <- mervAC[which(mervAC$theName %in% ne01),]
  rownames(mervAC) <- mervAC$theName; mervAC <- mervAC[,c(1,3,5,6)]
  cnt1 <- cbind(mervKD, mervAC)
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  cnt2 <- cnt1
  cnt2 <- cnt2 +0.001
  cnt2[1,]
  cnt3 <- as.data.frame(apply(cnt2[,c(1:4)], 2, function(x){x/rowMeans(cnt2[,c(3,4)])}))
  cnt4 <- as.data.frame(apply(cnt2[,c(5:8)], 2, function(x){x/rowMeans(cnt2[,c(7,8)])}))
  dfa1 <- cbind(cnt3,cnt4)
  dfa1 <- log2(dfa1 +0.001)
  dfa1 <- dfa1[which(rowSums(is.na(dfa1)) <1),]
  #wr01 <- dfa1
  #wr01$geneName <- rownames(wr01)
  #write.csv(x = wr01, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-ff/normCPM-1165.csv", quote = F)
  dfa1[1,]
  pheatmap(mat = dfa1,
           main="\nthe heatmap of the venn 1165 Genes",
           scale = "none", cellwidth = 30,
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(dfa1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#28、MM-融合转录本、1165基因的重合, 得到159个基因, 这些基因在MM disturb的热图
{
  ge01 <- import(inx2)
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
  mervKD <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/count/MMKD.onlyGene.xlsx")
  mervKD <- mervKD[which(mervKD$theName %in% ne02),]
  rownames(mervKD) <- mervKD$theName; mervKD <- mervKD[,c(1,2,5,6)]
  mervAC <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/count/mt2Ac.onlyGene.xlsx")
  mervAC <- mervAC[which(mervAC$theName %in% ne02),]
  rownames(mervAC) <- mervAC$theName; mervAC <- mervAC[,c(1,3,5,6)]
  cnt1 <- cbind(mervKD, mervAC)[,c(3,4,1,2,7,8,5,6)]
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nthe 159 genes -- MM disturb",
           scale = "row", cellwidth = 30,
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#
#
#
#
################################################################################
#29、MM KD后MM 位点表达的boxplot
{
  path <- "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/"
  prex <- c("KD-1","KD-2","WT-1","WT-2")
  cnt1 <- read.table(paste0(path,"countUniq/MM.site.cntTable"), header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- prex
  colLen <- length(colnames(cnt1))
  cnt1 <- cnt1[rowSums(cnt1 == 0) < colLen, ]
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  openxlsx::write.xlsx(x = tmp1, file = paste0(path, "countUniq/MM.cpm.onlyTE.site.xlsx"))
  #
  #
  #
  mervs <- openxlsx::read.xlsx(paste0(path,"countUniq/MM.cpm.onlyTE.site.xlsx"))
  MM <- mervs[grep("MM-int",mervs$theName),]
  MM <- MM[rowSums(MM[,1:4]==0) < 4,]
  MM <- tidyr::pivot_longer(data = MM, cols = -theName, names_to = "group", values_to = "cpms")
  vaue <- wilcox.test(x = MM$cpms[grep("WT",MM$group)],y = MM$cpms[grep("KD",MM$group)], alternative = "greater")
  print(vaue$p.value)
  #
  MM$marks <- sapply(str_split(MM$group,"-"),"[",1)
  MM$marks <- factor(MM$marks, levels = rev(unique(MM$marks)))
  p_001 <- ggplot(data = MM,aes(x = marks, y = log2(cpms+0.1))) +
    geom_boxplot(aes(fill = marks)) +
    labs(title = paste0("MM-int site Expr"),
         y = "Log2 (CPM+0.1)", x = "", fill="Group") +
    geom_signif(comparisons = list(c("WT", "KD")),
                map_signif_level=T, family = "serif",
                textsize=6, test = wilcox.test, step_increase = 0) +
    stat_summary(geom = "point", fun = "median",shape = 22,size = 2,
                 position = position_dodge(width = 1)) +
    theme_classic() +
    scale_fill_manual(values = c("#3FA0C0","#F0A19A"),
                      labels = c("WT", "MM-KD")) +
    theme(plot.margin = margin(unit = "cm",1,1,1,1),
          legend.text = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          axis.text.x = element_text(angle = 45,hjust = 1,family = "serif",size = 12),
          axis.text.y = element_text(family = "serif", size = 12)) +
    scale_y_continuous(limits = c(-6,12)) +
    #scale_y_continuous(limits = c(0,10),breaks = c(0,1,2,3,4,5)) +
    scale_x_discrete(labels = c("WT", "MM-KD"))
  print(p_001)
  ggsave(plot = p_001,
         filename = paste0(path,"countUniq/","MM-KD.uniq.MM-int.boxplot.pdf"),
         units = "cm", width = 20, height = 15)
}









