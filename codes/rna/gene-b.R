#write by ~~~ at 2023.06.05
#01、搭环境
{
  for (i in c("data.table","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
  #"myenrichplot",
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#03、跑命令
{
  shell <- paste0(path,"src/run_b.sh")
  cat("#!/bin/bash\n", file = shell)
  prefix  <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*_1.fq.gz$*"),"_1.fq.gz"), "[",1)
  one <- paste0(path,"fq/",prefix,"_1.fq.gz")
  two <- paste0(path,"fq/",prefix,"_2.fq.gz")
  cmd_31 <- paste0("trim_galore --phred33 --fastqc -a 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' ",
                   "-a2 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' --clip_R1 10 --clip_R2 10 ",
                   "--three_prime_clip_R1 3 --three_prime_clip_R2 3 ",
                   "--paired ",one," ",two," -o ",paste0(path,"trim\n"))
  one <- paste0(path,"trim/",prefix,"_1_val_1.fq.gz")
  two <- paste0(path,"trim/",prefix,"_2_val_2.fq.gz")
  cmd_32 <- paste0("STAR --runThreadN 30 --readFilesCommand zcat ",
                   "--outMultimapperOrder Random --genomeDir ",inx1," ",
                   "--outFileNamePrefix ",path,"bam/",prefix," ",
                   "--outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 --readFilesIn ",
                   one," ",two,"\n")
  cmd_33 <- paste0("samtools sort -@ 20 ",paste0(path,"bam/",prefix,"Aligned.out.sam")," -o ",paste0(path,"bam/",prefix,".bam"),"\n")
  cmd_34 <- paste0("samtools index ",paste0(path,"bam/",prefix,".bam"),"\n")
  cmd_35 <- paste0("rm ",paste0(path,"bam/",prefix,"Aligned.out.sam"),"\n")
  cmd_36 <- paste0("bamCoverage  --normalizeUsing CPM --numberOfProcessors 20 --minMappingQuality 1 --samFlagExclude 256  --outFileFormat bigwig --binSize 20 -b ",
                   paste0(path,"bam/",prefix,".bam")," -o ",paste0(path,"bamCoverage/",prefix,".bw"),"\n")
  cmd_37 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prefix[1:4],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prefix[5:6],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",inx4," --sortByPos --mode multi ",
                   "--project ccFbl --outdir ",paste0(path,"count"),"\n")
  cmd_38 <- paste0("TElocal -b ",
                   paste0(path,"bam/",prefix,".bam"),
                   " --GTF ",inx3," --TE ",inx5," --sortByPos --project ",
                   paste0(path, "count/TElocal/",prefix),"\n")
  cmd_41 <- paste0("stringtie ", paste0(path,"bam/",prefix,".bam")," ",
                   "-p 30 -G ", inx3, " -o ", paste0(path,"stringTie/",prefix,".gtf"),"\n")
  cmd_42 <- paste0("stringtie --merge -o ",paste0(path,"stringTie/v_KDFbl.merge.gtf")," -G ",inx3," ",
                   paste(paste0(path,"stringTie/",prefix[c(1,2,5,6)],".gtf"), collapse = " "),"\n")
  cmd_43 <- paste0("stringtie --merge -o ",paste0(path,"stringTie/v_KDPm8.merge.gtf")," -G ",inx3," ",
                   paste(paste0(path,"stringTie/",prefix[c(3,4,5,6)],".gtf"), collapse = " "),"\n")
  cmd_44 <- paste0("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/stringTie/TEtranscripts -t ",
                   c(paste(paste0(path,"bam/",prefix[1:2],".bam"),collapse = " "),
                     paste(paste0(path,"bam/",prefix[3:4],".bam"),collapse = " ")),
                   " -c ", paste(paste0(path,"bam/",prefix[5:6],".bam"),collapse = " "),
                   " --GTF ",paste0(path,"stringTie/v_KD",c("Fbl","Pm8"),".merge.gtf"),
                   " --TE ",inx4," --sortByPos --mode multi",
                   " --project v_",c("Fbl","Pm8"),"merge"," --outdir ", paste0(path,"stringTie/count"),"\n")
  for (i in cmd_31) {cat(i, file = shell, append = T)}
  for (i in cmd_32) {cat(i, file = shell, append = T)}
  for (i in cmd_33) {cat(i, file = shell, append = T)}
  for (i in cmd_34) {cat(i, file = shell, append = T)}
  for (i in cmd_35) {cat(i, file = shell, append = T)}
  for (i in cmd_36) {cat(i, file = shell, append = T)}
  for (i in cmd_37) {cat(i, file = shell, append = T)}
  for (i in cmd_38) {cat(i, file = shell, append = T)}
  for (i in cmd_41) {cat(i, file = shell, append = T)}
  for (i in cmd_42) {cat(i, file = shell, append = T)}
  for (i in cmd_43) {cat(i, file = shell, append = T)}
  for (i in cmd_44) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#04、优化了的差异分析流程
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
  shy_diff(mark = "fbl", 
           geneFile = paste0(path, "count/ccFbl.cntTable"),
           siteFile = paste0(path, "count/ccFbl.site.cntTable"),
           Class = factor(c("KD","KD","WT","WT")), VS = c("Class","KD","WT"))
}
################################################################################
#05、CPM Excel
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("fblKD-1","fblKD-2","WT-1","WT-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/fbl.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/fbl.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/fbl.cpm.allGeneTE.xlsx")
}
################################################################################
#06、TE CPM point [带基因背景]
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("fblKD-1","fblKD-2","ctrl-1","ctrl-2")
  cnt1 <- as.data.frame(apply(cnt1, 2, function(x){x/sum(x)*1000000}))
  cnt1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]
  #
  cnt2 <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202104/aaPrmt1KO74/count/Prmt1.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt2) <- c("prmt1KO-1","prmt1KO-2","WT-1","WT-2")
  cnt2 <- as.data.frame(apply(cnt2, 2, function(x){x/sum(x)*1000000}))
  cnt2 <- cnt2[grep(x = rownames(cnt2), pattern=":", invert = F),]
  cnts <- cbind(cnt1, cnt2)
  myta <- data.frame(ccFblKD = log2(rowMeans(cnts[,c(1,2)]) +1),
                     ccFblNC = log2(rowMeans(cnts[,c(3,4)]) +1),
                     prmt1KO = log2(rowMeans(cnts[,c(5,6)]) +1),
                     prmt1WT = log2(rowMeans(cnts[,c(7,8)]) +1))
  #
  {
    myte <- read.csv(file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/fbl.te.GeneBackGround.csv", row.names = 1)#不带基因背景
    myte$mark <- NA
    myte$mark[which(myte$log2FoldChange >  log2(1.5) & myte$padj < 0.05)] <- "upup"
    myte$mark[which(myte$log2FoldChange < -log2(1.5) & myte$padj < 0.05)] <- "down"
    te01 <- rownames(myte[which(myte$mark=="upup"),])
    te02 <- rownames(myte[which(myte$mark=="down"),])
    myte <- read.csv(file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202104/aaPrmt1KO74/DESeq2/Prmt1KO.te.GeneBackGround.csv", row.names = 1)
    myte$mark <- NA
    myte$mark[which(myte$log2FoldChange >  log2(1.5) & myte$padj < 0.05)] <- "upup"
    myte$mark[which(myte$log2FoldChange < -log2(1.5) & myte$padj < 0.05)] <- "down"
    te03 <- rownames(myte[which(myte$mark=="upup"),])
    te04 <- rownames(myte[which(myte$mark=="down"),])
    #
    myta$fblMark <- "other"
    myta$fblMark[which(rownames(myta) %in% te01)] <- "upup"
    myta$fblMark[which(rownames(myta) %in% te02)] <- "down"
    myta$prmt1Mark <- "other"
    myta$prmt1Mark[which(rownames(myta) %in% te03)] <- "upup"
    myta$prmt1Mark[which(rownames(myta) %in% te04)] <- "down"
  }
  labe <- myta[grep(x = rownames(myta), pattern = "MERVL-int|MT2_Mm"),]
  myta$fblMark <- factor(x = myta$fblMark, levels = c("upup","down","other"))
  p_07 <- ggplot() +
    geom_point(data = myta[which(myta$fblMark =="other"),], aes(x = ccFblNC, y = ccFblKD), color = "gray") +
    geom_point(data = myta[which(myta$fblMark !="other"),], aes(x = ccFblNC, y = ccFblKD, color = fblMark)) +
    scale_color_manual(values = c("red","blue")) +
    geom_point(data = labe, aes(x = ccFblNC, y = ccFblKD),size=2,shape=17,color = "red",fill="red") +
    geom_text_repel(data = labe, aes(x = ccFblNC, y = ccFblKD, label=rownames(labe))) +
    labs(x="Log2(CPM+1)(NC)",y="Log2(CPM+1)(Fbl KD)",title = "TE: Fbl KD vs Control (gene backGround)") +
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
         filename = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/fblKD.cpm.point.pdf")
  #
  heta <- myta
  heta$prmt1Mark <- factor(x = heta$prmt1Mark, levels = c("upup","down","other"))
  p_08 <- ggplot() +
    geom_point(data = heta[which(heta$prmt1Mark =="other"),], aes(x = prmt1WT, y = prmt1KO), color = "gray") +
    geom_point(data = heta[which(heta$prmt1Mark !="other"),], aes(x = prmt1WT, y = prmt1KO, color = prmt1Mark)) +
    scale_color_manual(values = c("red","blue")) +
    geom_point(data = labe, aes(x = prmt1WT, y = prmt1KO),size=2,shape=17,color = "red",fill="red") +
    geom_text_repel(data = labe, aes(x = prmt1WT, y = prmt1KO, label=rownames(labe))) +
    labs(x="Log2(CPM+1)(NC)",y="Log2(CPM+1)(Prmt1 KO)",title = "TE: Prmt1 KO vs Control (gene backGround)") +
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
         filename = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202104/aaPrmt1KO74/PLOT/prmt1KO.cpm.point.pdf")
}
################################################################################
#07、TE MA plot
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
  shy_teMA(difffile = paste0(path,"DESeq2/fbl.te.GeneBackGround.csv"),
           savefile = paste0(path,"PLOT/fbl.te.GeneBackGround.MA.pdf"))
}
################################################################################
#
#
#
#
#
#
################################################################################
#08、差异Gene、TE、TE sites from DESeq2
{
  #reads from TEtranscripts
  data <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable",
                     header = T, check.names = F, stringsAsFactors = F)
  #colnames(data) <- c("gene",sapply(str_split(basename(colnames(data)), ".bam"),"[",1)[2:7])
  gene <- data[grep(data[,1], pattern = ":", invert = T),]
  tete <- data[grep(data[,1], pattern = ":", invert = F),]
  rc_a <- gene[,c(1,4,3,7,6)]
  rc_b <- gene[,c(1,5,2,7,6)]
  rc_c <- tete[,c(1,4,3,7,6)]
  rc_d <- tete[,c(1,5,2,7,6)]
  #reads from TElocal
  data <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/TElocal/FblPrmt8.txt",
                     header = T, check.names = F, stringsAsFactors = F)
  #colnames(data) <- c("gene",sapply(str_split(basename(colnames(data)), ".bam"),"[",1)[2:7])
  rc_e <- data[,c(1,4,5,6,7)]
  rc_f <- data[,c(1,2,3,6,7)]
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
    resres <<- results(ddsdds, contrast = c("Class","KD","WT"))
    #return(resres)
  }
  ##NO.1
  shy_diff(mydata = rc_a, Class = factor(c("KD","KD","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","Prmt8.gene.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","Prmt8.gene.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","Prmt8.gene.down.csv"))
  ##NO.2
  shy_diff(mydata = rc_b, Class = factor(c("KD","KD","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","Fbl.gene.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","Fbl.gene.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","Fbl.gene.down.csv"))
  ##NO.3
  shy_diff(mydata = rc_c, Class = factor(c("KD","KD","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","Prmt8.fullte.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","Prmt8.fullte.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","Prmt8.fullte.down.csv"))
  ##NO.4
  shy_diff(mydata = rc_d, Class = factor(c("KD","KD","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","Fbl.fullte.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","Fbl.fullte.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","Fbl.fullte.down.csv"))
  ##NO.5
  shy_diff(mydata = rc_e, Class = factor(c("KD","KD","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","Prmt8.teSites.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","Prmt8.teSites.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","Prmt8.teSites.down.csv"))
  ##NO.6
  shy_diff(mydata = rc_f, Class = factor(c("KD","KD","WT","WT")))
  write.csv(resres,paste0(path,"DESeq2/","Fbl.teSites.csv"))
  upup <- subset(resres, log2FoldChange > log2(1.5) & padj<0.05)
  down <- subset(resres, log2FoldChange < -log2(1.5)& padj<0.05)
  write.csv(upup,paste0(path,"DESeq2/","Fbl.teSites.upup.csv"))
  write.csv(down,paste0(path,"DESeq2/","Fbl.teSites.down.csv"))
}
################################################################################
#09、计算Fbl CPM
{
  data <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable",
                     header = T, check.names = F, stringsAsFactors = F, row.names = 1)
  colnames(data) <- sapply(str_split(basename(colnames(data)), ".bam"),"[",1)
  data <- data[,c(4,1,6,5)]
  #放一起做标准化
  myta <- as.data.frame(apply(data, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- myta[,c(5,1:4)]
  openxlsx::write.xlsx(x = myta, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/Fbl-KD.CPM.xlsx")
  #gene单独标准化
  myta <- data[grep(rownames(data), pattern = ":", invert = T),]
  myta <- as.data.frame(apply(myta, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- myta[,c(5,1:4)]
  openxlsx::write.xlsx(x = myta, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/Fbl-KD.gene.CPM.xlsx")
  #tete单独标准化
  myta <- data[grep(rownames(data), pattern = ":", invert = F),]
  myta <- as.data.frame(apply(myta, 2, function(x){x/sum(x)*1000000}))
  myta$geneName <- rownames(myta)
  myta <- myta[,c(5,1:4)]
  openxlsx::write.xlsx(x = myta, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/Fbl-KD.tete.CPM.xlsx")
}
################################################################################
#10、画图：火山图 for Genes
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
  shy_volcano_gene(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.volcano.pdf")
  shy_volcano_gene(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/fbl.geneOnly.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/Volcano/fbl.geneOnly.volcano.pdf")
}
################################################################################
#11、画图：火山图 for TEs
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
      labb <- labb[which(labb$padj<10^(-20)),]#or -10
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
  shy_volcano_te(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.fullte.csv",
                 saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/Prmt8.fullte.volcano.pdf")
  shy_volcano_te(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.fullte.csv",
                 saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/Fbl.fullte.volcano.pdf")
}
################################################################################
#12、画图：富集分析
{
  list.files(paste0(path,"DESeq2"), pattern = "*gene.csv$", full.names = T)
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
           filename = paste0(savePath, change, ".KEGG.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
    ggsave(plot = p_2, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".ReactomePA.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
    ggsave(plot = p_3, 
           units = "cm", width = 20, height = 16,
           filename = paste0(savePath, change, ".GO.",str_split(basename(resFile),".csv")[[1]][1],".pdf"))
  }
  ##条形图
  shy_fuji_bar(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv",
               savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",
               change  = "upup")
  shy_fuji_bar(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv",
               savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",
               change  = "down")
  ##气泡图
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
  shy_fuji_bubble(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.gene.csv",
                  savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",
                  change  = "upup")
  shy_fuji_bubble(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.gene.csv",
                  savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",
                  change  = "down")
  ##Temp
  shy_fuji_temp <- function(resFile, savePath){
    kegg <- as.data.frame(read.table(resFile, header = T, stringsAsFactors = F, sep = "\t"))
    ##画图
    kegg$Description <- factor(kegg$Description, levels = rev(rbind(kegg$Description)))
    p_3 <- ggplot(kegg, aes(x=-log10(pvalue), y=Description))+
      geom_point(aes(color=-log10(pvalue), size=Count)) +
      theme_classic() +
      ylab(NULL) +
      scale_color_gradient2(low = "white", high = "red3", midpoint = 1) +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            axis.title = element_text(size = 12,family = "serif", color = "black"),
            axis.text  = element_text(size = 12,family = "serif", color = "black"))
    ggsave(plot = p_3, 
           units = "cm", width = 20, height = 12,
           filename = paste0(savePath, str_split(basename(resFile),".txt")[[1]][1],".bubble.pdf"))
  }
  shy_fuji_temp(resFile="/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/Enrichment/MERVL_pull_downKEGG.txt",
                savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/Enrichment/")
  
}
################################################################################
#13、画图：表达量热图
{
  shy_pheatmap <- function(readsFile, gene_resFile, tete_resFile, savePath, mark){
    print("reads counted by TEtranscripts")
    data <- read.table(readsFile, header = T, row.names = 1, check.names = F, stringsAsFactors = F)
    colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = "_"), "[", 1)
    gene <- data[grep(rownames(data), pattern = ":", invert = T),c(3,2,4,1,6,5)]
    gene <- gene[rowSums(gene) > 5,]
    gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
    tete <- data[grep(rownames(data), pattern = ":", invert = F),c(3,2,4,1,6,5)]
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
  shy_pheatmap(readsFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable",
               gene_resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv",
               tete_resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.fullte.csv",
               savePath  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",
               mark = "heat")
  shy_pheatmap(readsFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable",
               gene_resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.gene.csv",
               tete_resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.fullte.csv",
               savePath  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",
               mark = "heat")
  
}
################################################################################
#14、画图：质谱数据 [DIY]
{
  ##读质谱数据（X轴）
  trea <- read.delim("/disk5/aaSHY/1shyJOB/4protein/data/LUask/treat/480_FPEP210002742-1A_Proteins-Prmt1-flag.txt",
                     sep = "\t", header = T, stringsAsFactors = F)[,c(1,2,7)]
  trea <- tidyr::separate_rows(data = trea, "gene_name", sep = ";")
  ctrl <- read.delim("/disk5/aaSHY/1shyJOB/4protein/data/LUask/ctrl/480_FPEP210002741-1A_Proteins-OV-flag.txt",
                     sep = "\t", header = T, stringsAsFactors = F)[,c(1,2,7)]
  ctrl <- tidyr::separate_rows(data = ctrl, "gene_name", sep = ";")
  group_t <- trea[which(trea$gene_name %in% intersect(trea$gene_name, ctrl$gene_name)),]
  group_c <- ctrl[which(ctrl$gene_name %in% intersect(trea$gene_name, ctrl$gene_name)),]
  group_t <- stats::aggregate(group_t$X..Peptides, by=list(type=group_t$gene_name), sum)
  group_c <- stats::aggregate(group_c$X..Peptides, by=list(type=group_c$gene_name), sum)
  group_m  <- merge(x = group_t, y = group_c, by = "type")
  #View(group_m)
  group_m$ratio <- log2(group_m$x.x/group_m$x.y)#group_m$ratio <- log2(group_m$x.x/group_m$x.y)
  ##读表达量数据（Y轴）
  data <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable",
                     header = T, row.names = 1, check.names = F, stringsAsFactors = F)
  colnames(data) <- sapply(str_split(basename(colnames(data)), pattern = "_"), "[", 1)
  gene <- data[grep(rownames(data), pattern = ":", invert = T),c(6,5)]
  gene <- as.data.frame(apply(gene, 2, function(x){x/sum(x)*1000000}))
  gene$expr <- log10((gene[,1]+gene[,2])/2+1)
  gene$type <- rownames(gene)
  mydata <- merge(x=group_m, y=gene, by="type")
  mydata[1:2,]
  ##读差异分析的数据（红点的显著基因）
  diff_g <- read.csv("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv")
  significant <- diff_g[which(diff_g$log2FoldChange > log2(1.5) & diff_g$padj < 0.05), "X"]
  mydata$diff[mydata$type %in% significant] = "Up_regulated"
  mydata$diff[!(mydata$type %in% significant)] = "Non_significant"
  ##开始画图
  okdata <- mydata[,-c(2,3,5,6)]
  ggplot(data = okdata) +
    geom_point(aes(x=ratio, y=expr, color=diff)) +
    geom_text_repel(data = okdata[which(okdata$type=="Angel2"),], aes(x=ratio,y=expr, label=type))
  ##筛出目标数据
  okdata[which(okdata$ratio > 0 & okdata$ratio < 2 & okdata$expr>3),"type"]
  okdata[which(okdata$ratio > 0 & okdata$ratio < 2 & okdata$expr > 1 & okdata$expr < 2),"type"]
  okdata[which(okdata$diff=="Up_regulated" & okdata$expr > 1 & okdata$expr < 2 & okdata$ratio > -1),"type"]
  write.csv(x = okdata, file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/ProteinMS/peptides.vs.74.csv")
  ##Unique展示
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
################################################################################
#15、画图：富集分析 [GSEA]
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
    aa <- clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/2cellgene_GMT.gmt")
    aa <- suppressWarnings(clusterProfiler::bitr(aa$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
    ab <- data.frame(ont = rep("act_2cell_gene", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
    af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 1000, TERM2GENE = ab, pvalueCutoff = 1))
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
  shy_gsea(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv",
           savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/EnrichmentGSEA/")
  shy_gsea(resFile  = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.gene.csv",
           savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/EnrichmentGSEA/")
}
################################################################################
#16、画图：PCA [DIY]
{
  #----------------------------------PCA-batch----------------------------------
  #hisat2 + featureCounts + DESeq2 + ggplot2
  path = "/ChIP_seq_2/aaSHY/Dot1l/shy_dot1ReviewRNAseq/"
  prmt <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/ccFbl.cntTable",header=T, check.names = F, stringsAsFactors = F, row.names = 1)
  colnames(prmt) <- sapply(str_split(basename(colnames(prmt)), ".bam"),"[",1)
  prmt <- prmt[,c(6,5,3,2,4,1)]
  prmt <- prmt[grep(rownames(prmt), pattern = ":", invert = T),]
  dot1 <- read.table("/RNA_seq_2/sht/N1906264_GZZ_80-354364506/shyTEST/count/Dot1lRNAseq_KOvsWT.count", header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)
  dot1 <- dot1[which(rownames(dot1) %in% rownames(prmt)),]
  prmt <- prmt[which(rownames(prmt) %in% rownames(dot1)),]
  
  tblc <- read.table(paste0(path, "TBLCs/count/gene.count.txt"), header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)[,-c(1:5)]
  tblc <- tblc[which(rownames(tblc) %in% rownames(dot1)),]
  tclc <- read.table(paste0(path, "TCLCs/count/gene.count.txt"), header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)[,-c(1:5)]
  tclc <- tclc[which(rownames(tclc) %in% rownames(dot1)),]
  ercl <- read.table("/RNA_seq_1/zx/TOMATO/Count/gencode.vM21.primary_assembly.annotation.gtf_cut",
                     header=F, check.names = F, stringsAsFactors = F, row.names = 1, skip = 2)#[,c(2,3)]
  ercl <- ercl[which(rownames(ercl) %in% rownames(dot1)),]
  vtax <- cbind(prmt,dot1,tblc,tclc,ercl)
  colnames(vtax)  <- c(paste0("i_wt",seq(2)),paste0("i_ko_a",seq(2)),paste0("i_ko_b",seq(2)),paste0("ko",seq(3)),paste0("wt",seq(3)),paste0("tb",seq(3)),paste0("tc",seq(4)),paste0("erc",seq(2)))
  tabl <- data.frame(names = colnames(vtax),
                     condition = sapply(str_split(colnames(vtax),"[0-9]"), "[", 1), 
                     batch = c(rep("A",2), rep("B",2), rep("C",2), rep("D",3), rep("E",3), rep("F",3), rep("G",4), rep("H",2)))
  countmatrix <- as.matrix(vtax)
  dds1 <- DESeqDataSetFromMatrix(countmatrix, colData = tabl, design = ~condition)#[rowSums(counts(vtax))>5]
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
                           labels = c("Prmt8","Fbl","shCtrl","Dot1KO", "Dot1WT", "Tomato","TBLC","TLSC"))
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
  pdf(file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/PCA/PCA.pdf",width = 10, height = 8)
  shy_linshi(df = pd_1, i = 1)
  shy_linshi(df = pd_2, i = 2)
  shy_linshi(df = pd_3, i = 3)
  shy_linshi(df = pd_4, i = 4)
  dev.off()
  rm(prmt,dot1,tblc,tclc,ercl,dds1,dds2,ppp,pd_1,pd_2,pd_3,pd_4,rld1,rldf,vsdf,rawx,norx,tabl,vtax)
}
################################################################################
#17、画图：Corelation
{
  shell <- paste0(path,"src/run_Correlation.sh")
  cat("#!/bin/bash\n", file = shell)
  aa <- list.files(paste0(path,"bamCoverage"), pattern = "*bw$", full.names = T)
  ab <- list.files("/ChIP_seq_2/aaSHY/Dot1l/shy_dot1ReviewRNAseq/corelation/bw", pattern = "*bw$", full.names = T)
  ac <- list.files("/RNA_seq_1/Pim3_KD/bamCoverage", pattern = "*CPM.bw$", full.names = T)
  ad <- list.files("/RNA_seq_1/Ssrp1/bamCoverage", pattern = "*CPM.bw$", full.names = T)
  cmd_a1 <- paste0("multiBigwigSummary bins -b ",paste0(c(aa,ab,ac,ad), collapse = " ")," -o ",
                   paste0(path,"PLOT/Correlation/","Prmt1.Correlation.npz","\n"))
  cmd_a2 <- paste0("plotCorrelation -in ",
                   paste0(path,"PLOT/Correlation/","Prmt1.Correlation.npz "),
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers --colorMap bwr -o ",
                   paste0(path,"PLOT/Correlation/","Prmt1.Correlation.pdf "),
                   "--outFileCorMatrix ",paste0(path,"PLOT/Correlation/","Prmt1.Correlation.index.txt\n"))
  for (m in cmd_a1){cat(m, file = shell, append = T)}
  for (m in cmd_a2){cat(m, file = shell, append = T)}
  rm(aa,ab,ac,ad,cmd_a1,cmd_a2)
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#18、画图：UCSC Track
{
  #1, upload bigwig file to cyverse, paste the links
  #2, make UCSC custom track
  #3, tuning parameters
  ####https://genome.ucsc.edu/s/hysun088/Prmt1
}
################################################################################
#19、画图：Up TEs loci number [TElocal]
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
      labs(title = expression(Upregulated~TEs~"in"~italic(Fbl)^{"-/-"}~ESCs)) +
      scale_fill_manual(breaks = unique(forP$clas),values = ggsci::pal_aaas(palette = "default")(length(unique(forP$clas))))
    ggsave(plot = p_1, filename = saveFile, units = "cm", width = 12, height = 12)
  }
  shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.teSites.upup.csv",
              saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TElociNum/Prmt8.teSites.upup.pdf")
  shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.teSites.upup.csv",
              saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TElociNum/Fbl.teSites.upup.pdf")
  shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.teSites.down.csv",
              saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TElociNum/Prmt8.teSites.down.pdf")
  shy_lociNum(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.teSites.down.csv",
              saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TElociNum/Fbl.teSites.down.pdf")
}
################################################################################
#20、画图：MERVL fusion genes expression
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
    resres <- results(ddsdds, contrast = c("Class","KD","WT"))
    resres <<- na.omit(as.data.frame(resres))
  }
  shy_diff_stringTie(countFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/stringTie/count/v_Pm8merge.cntTable", Class = c("KD","KD","WT","WT"))
  ##the Fusion Gene & PLOT
  gtff <- rtracklayer::import(inx4)
  shy_fusion_stringTie <- function(gtfFile, whichTEs, mark, savePath){
    gf00 <- rtracklayer::import(gtfFile)
    #gf00 <- as.data.frame(gf00[which(gf00$type=="transcript")])
    #dfaa <- gf00 %>% 
    #        group_by(gene_id) %>% 
    #        summarise(seqnames=unique(seqnames), start=min(start), end=max(end))
    #gfaa <- makeGRangesFromDataFrame(dfaa, keep.extra.columns = T)
    gfaa <- gf00[which(gf00$type=="transcript")]
    gfbb <- gtff[grep(gtff$gene_id, pattern = paste0(whichTEs, ".*"))]
    oovv <- suppressWarnings(findOverlaps(query = gfaa, subject = gfbb))
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
            panel.border     = element_rect(colour = "black", linewidth = 0.8),
            axis.line        = element_line(),
            axis.text        = element_text(family = "serif")) +
      xlab(label = expression(Log2~(Fold~Change~of~italic(Prmt1)^{"-/-"}~ESCs))) +
      ylab(label = expression(Log10~(adjusted~italic(P)~plain(value)))) +
      geom_point(data=resB,aes(x=log2FoldChange,y=-log10(padj),colour="Other Genes"),size=1.2) +
      scale_color_manual(values = c("red","gray")) +
      geom_point(data=resA,aes(x=log2FoldChange,y=-log10(padj),colour=paste0(whichTEs,"-fushion Genes")),size=1.2) +
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
    ggsave(plot = p_1, 
           filename = paste0(path, "PLOT/fusionTEs/",mark,"_fusion_",whichTEs,".pdf"), 
           width = 15, height = 15, units = "cm")
  }
  shy_fusion_stringTie(gtfFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/stringTie/v_KDPm8.merge.gtf",
                       whichTEs = "MT2_Mm",#MT2#MERVL
                       mark = "Prtm8",
                       savePath = paste0(path, "PLOT/fusionTEs/"))
  shy_fusion_stringTie(gtfFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/stringTie/v_KDFbl.merge.gtf",
                       whichTEs = "MT2_Mm",#MT2#MERVL
                       mark = "Fbl",
                       savePath = paste0(path, "PLOT/fusionTEs/"))
}
################################################################################
#21、画图：repSites of selected TEs
{
  shy_TEsites_plot <- function(resFile, saveFile){
    gf00 <- read.csv(file = resFile,header = T, stringsAsFactors = F, row.names = 1)
    gf00 <- na.omit(gf00)
    gfaa <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*"),]
    gfbb <- gf00[grep(rownames(gf00), pattern = ".*MT2.*"),]
    gfcc <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*|.*MT2.*", invert = T),]
    #开始作图
    p_1 <- ggplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = linewidth(colour = "black", size = 0.8),
            axis.line = element_line(), axis.text = element_text(family = "serif")) +
      xlab(label = expression(Log2~(Fold~Change~of~italic(Prmt1)^{"-/-"}~ESCs))) +
      ylab(label = expression(Log10~(adjusted~italic(P)~plain(value)))) +
      scale_color_manual(breaks = c("MERVL","MT2","Others"),values = c("red","purple","gray")) +
      geom_point(data=gfcc,aes(x=log2FoldChange,y=-log10(padj),colour="Others"),size=1.2) +
      geom_point(data=gfbb,aes(x=log2FoldChange,y=-log10(padj),colour="MT2"),size=1.2) +
      geom_point(data=gfaa,aes(x=log2FoldChange,y=-log10(padj),colour="MERVL"),size=1.2) +
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
  shy_TEsites_plot(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.teSites.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TEsitesPlot/Prmt8.mervl_mt2.pdf")
  shy_TEsites_plot(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.teSites.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TEsitesPlot/Fbl.mervl_mt2.pdf")
}
################################################################################
#22、不带基因背景的TE的MA plot
{
  geff <- rtracklayer::import(con = inx3)
  geff <- geff[which(geff$type=="gene")]
  teff <- rtracklayer::import(con = inx4)
  mkov <- suppressWarnings(GenomicRanges::findOverlaps(query = geff, subject = teff, ignore.strand = T))
  teng <- teff[setdiff(seq(length(teff)),unique(subjectHits(mkov)))]
  rtracklayer::export(object = teng, 
                      con = "/Reference/aaSHY/appTools/TEtoolkit/noCHR/mm10_rmsk_TE_without_chr_noGene.gtf",
                      format = "gtf")
  cmd_37 <- paste0("featureCounts -a ",
                   "/Reference/aaSHY/appTools/TEtoolkit/noCHR/mm10_rmsk_TE_without_chr_noGene.gtf",
                   " -p -d 0 -D 1000 -T 20 -o ",
                   "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/featureCounts/Prmt8Fbl.reads.cntTable ",
                   paste(paste0(path,"bam/",prefix,".bam"),collapse = " "),"\n")
  data <- read.table("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/featureCounts/Prmt8Fbl.reads.txt",
                     header = T, check.names = F, stringsAsFactors = F, skip = 1)
  #colnames(data) <- c("gene",sapply(str_split(basename(colnames(data)), ".bam"),"[",1)[2:7])
  rc_e <- data[,c(1,4,5,6,7)]
  rc_f <- data[,c(1,2,3,6,7)]
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
    resres <<- results(ddsdds, contrast = c("Class","KD","WT"))
    #return(resres)
  }
  get_MA_plot<-function(relation){
    annot<-read.delim(relation,header=T,stringsAsFactors=F)
    rownames(annot)<-annot$repName
    terc<- as.data.frame(resres)
    terc<-terc[rownames(terc)%in%annot$repName,]
    annot<-annot[rownames(terc),]
    terc<-cbind(terc,annot)
    terc_s<-terc[which((terc$padj<0.05&(abs(terc$log2FoldChange)>log2(1.5)))&(terc$repClass=="LTR"|terc$repClass=="LINE"|terc$repClass=="SINE"|terc$repClass=="DNA")),]
    terc_o<-terc[which(!((terc$padj<0.05&(abs(terc$log2FoldChange)>log2(1.5)))&(terc$repClass=="LTR"|terc$repClass=="LINE"|terc$repClass=="SINE"|terc$repClass=="DNA"))),]
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
      theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                       axis.line = element_line(),axis.text = element_text(family = "serif"),
                       panel.border = element_rect(colour = "black", size = 0.8))+
      theme(legend.position = c(0.91,0.2),legend.background = element_blank(),
            legend.text = element_text(family = "serif")) +
      xlab("log2(Expression in baseMean)")+ylab("log2Foldchange")
    #ggsave("fullTEsMAplot_padj.pdf", height = 4.16, width = 4.7)
  }
  ##
  shy_diff(mydata = rc_e, Class = factor(c("KD","KD","WT","WT")))
  get_MA_plot(relation="/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/relation_fullTEs.txt")
  shy_diff(mydata = rc_f, Class = factor(c("KD","KD","WT","WT")))
  get_MA_plot(relation="/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/relation_fullTEs.txt")
  resres[which(rownames(resres)%in%c("MERVL-int","MT2_Mm")),]
}
################################################################################
#23、TE sites MA plot
{
  shy_TEsites_plot <- function(resFile, saveFile){
    gf00 <- read.csv(file = resFile, header = T, stringsAsFactors = F, row.names = 1)
    gf00 <- na.omit(gf00)
    gf00 <- gf00[which(sapply(stringr::str_split(rownames(gf00), pattern = ":"), "[", 1) %in% teng),]
    gfaa <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*"),]
    gfbb <- gf00[grep(rownames(gf00), pattern = ".*MT2.*"),]
    gfcc <- gf00[grep(rownames(gf00), pattern = ".*MERVL-int.*|.*MT2.*", invert = T),]
    #开始作图
    p_1 <- ggplot() +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 0.8),
            axis.line = element_line(), axis.text = element_text(family = "serif")) +
      ylab(label = expression(Log2Fold~Change)) +
      xlab(label = expression(Log2~(baseMean + 1))) +
      scale_color_manual(breaks = c("MERVL","MT2","Others"),values = c("red","purple","gray")) +
      geom_point(data=gfcc,aes(x=log2(baseMean+1),y=log2FoldChange,colour="Others"),size=1.2) +
      geom_point(data=gfbb,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MT2"),size=1.2) +
      geom_point(data=gfaa,aes(x=log2(baseMean+1),y=log2FoldChange,colour="MERVL"),size=1.2) +
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
  shy_TEsites_plot(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.teSites.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TEsitesPlot/Prmt8.noGene.mervl_mt2.pdf")
  shy_TEsites_plot(resFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Fbl.teSites.csv",
                   saveFile = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/TEsitesPlot/Fbl.noGene.mervl_mt2.pdf")
}
################################################################################
#24、MERVL-int locis changed in Prmt8 KD or Fbl KD (fold change)
{
  fbll <- read.csv(file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/v1/Fbl.teSites.csv", header = T, stringsAsFactors = F, row.names = 1)
  fbll <- na.omit(fbll)
  fbll <- fbll[grep(rownames(fbll), pattern = "MERVL-int"),]
  prmt <- read.csv(file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/v1/Prmt8.teSites.csv", header = T, stringsAsFactors = F, row.names = 1)
  prmt <- na.omit(prmt)
  prmt <- prmt[grep(rownames(prmt), pattern = "MERVL-int"),]
  fbll <- fbll[which(rownames(fbll) %in% unique(c(rownames(prmt),rownames(fbll)))),]
  prmt <- prmt[which(rownames(prmt) %in% unique(c(rownames(prmt),rownames(fbll)))),]
  myta <- data.frame()
  for (i in unique(c(rownames(prmt),rownames(fbll)))){
    temp <- data.frame(a1 = i,
                       a2 = prmt[i,2],
                       a3 = fbll[i,2])
    myta <- rbind(myta, temp)
  }
  myta[is.na(myta)] <- 0
  ggplot(data = myta) +geom_point(aes(x=a2,y=a3), color="steelblue") +
    labs(x="log2Foldchange of Prmt8 KD vs shCtrl", y="log2Foldchange of Fbl KD vs shCtrl", 
         title = "627 MERVL-int locis changed in Prmt8 KD or Fbl KD") +
    scale_y_continuous(limits = c(-4,6)) +
    scale_x_continuous(limits = c(-4,6)) +
    geom_hline(aes(yintercept =  log2(1.5)),linetype="dashed",colour="grey") +
    geom_hline(aes(yintercept = -log2(1.5)),linetype="dashed",colour="grey")+
    geom_vline(aes(xintercept =  log2(1.5)),linetype="dashed",colour="grey")+
    geom_vline(aes(xintercept = -log2(1.5)),linetype="dashed",colour="grey")+
    theme_classic()
  #部分基因在早期胚胎发育阶段的CPM表达值 [暂时找不到CPM表格]
  {
    teng <- sapply(str_split(as.character(myta[which(myta$a2>log2(1.5) & myta$a3>log2(1.5)),1]), ":"),"[",1)
    upte <- teff[which(teff$transcript_id %in% teng)]
    upte <- GenomicRanges::flank(x = upte, width = 10000, start = T)
    upte <- GenomicRanges::flank(x = upte, width = 10000, start = F)
    mtov <- suppressWarnings(GenomicRanges::findOverlaps(query = upte, subject = geff, ignore.strand=T))
    twoc <- base::intersect(geff[unique(subjectHits(mtov))]$gene_name,as.data.frame(clusterProfiler::read.gmt("/Reference/aaSHY/GSEAgmt/TwoGenes.gmt"))$gene)
    fjin <- geff[unique(subjectHits(mtov))]$gene_name
    stag <- read.csv("/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/count/star/cpmfile")
    temp <- stag[which(stag$X %in% fjin),]
    temp <- tidyr::gather(data = temp, key = "name", value = "cpm", -"X")
    ggplot(data = temp) +geom_boxplot(aes(x=name, y=cpm), fill="steelblue4") +ylim(0,50) +theme_classic() +labs(y="CPM")
    ggplot(data = temp[which(temp$name %in% c("X2cell","X4cell")),]) +geom_line(aes(x=name, y=cpm, group=X), color="steelblue4") +theme_classic() +labs(y="CPM")
    stag[which(stag$X=="Arap2"),]
    }
}
################################################################################
#25、Prmt8转录组在Prmt1上调基因的GSEA富集
{
  Prmt8Res <- "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/DESeq2/Prmt8.gene.csv"
  Prmt1Res <- "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202104/DESeq2/Prmt1.KO74.gene.csv"
  savePath = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/EnrichmentGSEA/"
  aa <- read.csv(file = Prmt8Res, header = T, stringsAsFactors = F)[,c(1,3)]
  ab <- suppressWarnings(clusterProfiler::bitr(geneID = aa$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = T))
  ac <- dplyr::distinct(.data = ab, SYMBOL, .keep_all = T)
  ac$rank <- apply(ac, 1, function(x){aa[which(aa$X==x[1]),2]})
  ad <- ac[order(ac$rank, decreasing = T),c(2,3)]
  ae <- ad$rank; names(ae) <- ad$ENTREZID
  rm(aa,ab,ac,ad)
  #
  aa <- read.csv(file = Prmt1Res, header = T, stringsAsFactors = F)
  aa <- aa[which(aa$log2FoldChange > log2(1.5) & aa$padj < 0.05),]
  aa <- suppressWarnings(clusterProfiler::bitr(aa$X, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
  ab <- data.frame(ont = rep("Prmt8UpGenes", length(unique(aa$ENTREZID))), gen = unique(aa$ENTREZID))
  af_gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = ae, maxGSSize = 3000, TERM2GENE = ab, pvalueCutoff = 1))
  rm(aa,ab)
  #
  p_a <- enrichplot::gseaplot2(x = af_gsea, geneSetID = 1, title = "GSEA for Prmt8UpGenes", pvalue_table = T)
  pdf(file = paste0(savePath, str_split(basename(resFile),".csv")[[1]][1], ".GSEA.pdf"), width = 12, height = 10)
  print(p_a)
  dev.off()
}
################################################################################
#26、不带基因背景 (先获取不带基因背景的TE的GTF，再使用TEtranscripts统计reads数量)
#########################方法不对, 舍弃#########################################
#########################方法不对, 舍弃#########################################
#########################方法不对, 舍弃#########################################
{
  gf01 <- rtracklayer::import(con = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf")
  gf01 <- gf01[which(gf01$type=="exon")]
  tf01 <- rtracklayer::import(con = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf")
  ov01 <- GenomicRanges::findOverlaps(query = gf01, subject = tf01, ignore.strand = T)
  tf02 <- tf01[setdiff(x = seq(length(tf01$type)), y = unique(subjectHits(ov01)))]
  rtracklayer::export(object = tf02, con = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE_noGene.gtf")
  #
  prefix  <- sapply(str_split(list.files((paste0(path,"fq/")),pattern = "*_1.fq.gz$*"),"_1.fq.gz"), "[",1)
  cmd_37 <- paste0("TEtranscripts -t ",
                   paste(paste0(path,"bam/",prefix[1:2],".bam"),collapse = " ")," -c ",
                   paste(paste0(path,"bam/",prefix[5:6],".bam"),collapse = " "),
                   " --GTF ",inx3," --TE ",
                   "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE_noGene.gtf ",
                   "--sortByPos --mode multi ",
                   "--project Fbl_noGene --outdir ",paste0(path,"count"),"\n")
  #
  red1 <- read.table(file = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/count/Fbl_noGene.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
  red1 <- red1[,c(1,4,5,6)]
  colnames(red1) <- c("KD-1","KD-2","Ctrl-1","Ctrl-2")
  #red2 <- as.data.frame(apply(red1, 2, function(x){x/sum(x)*1000000}))
  #red2[grep(rownames(red2), pattern = "MERVL-int"),]
  #
  Class <- base::factor(c(rep("KD",2),rep("Ctrl",2)))
  red2 <- red1#[grep(rownames(red1), pattern=":", invert = F),]
  inputs <- red2
  inputs <- inputs[base::max(inputs)>1, ]
  yyyCOL <- data.frame(row.names = colnames(inputs), Class, stringsAsFactors = T)
  ddsdds <- DESeq2::DESeqDataSetFromMatrix(countData = inputs, colData = yyyCOL, design = ~Class)
  ddsdds$Class = relevel(ddsdds$Class,ref="Ctrl")
  ddsdds <- DESeq2::DESeq(ddsdds)
  resres <- DESeq2::results(ddsdds, contrast = c("Class","KD","Ctrl"))
  abUpup <- as.data.frame(subset(resres, log2FoldChange > log2(1.5) & padj<0.05))
  resres[grep(rownames(resres), pattern="MERVL-int"),]
}
################################################################################
#27、画下Fbl在早期胚胎发育阶段的表达
{
  #01、用2014 Deng的数据
  cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
  prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
  colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
  cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = T),]#invert = F, 那么输出TE的CPM
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  df03 <- data.frame()
  for (i in c("Fbl")) {
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i, "$")),]
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
  ggsave(plot = p_01, 
         filename = "/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/Fpl-Embryogenesis-CPM.pdf",
         units = "cm", height = 15, width = 15)
}
{
  #01、用2014 Deng的数据
  cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
  cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
  prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
  prex <- prex[c(11,10,3,8,5,12,1,7,4)]
  colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
  cn01 <- cn01[,colnames(cn01)[grep("two-|blastoMid-|blastoLate-",colnames(cn01),invert = T)]]
  #
  #
  #
  cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
  cn02$geneName <- rownames(cn02)
  cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = T),]#invert = F, 那么输出TE的CPM
  df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
  df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
  df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
  #
  #
  #
  demo <- "Prmt1"
  df03 <- data.frame()
  for (i in demo) {
    temp <- df02[grep(df02$geneName, pattern = paste0("^", i, "$")),]
    df03 <- rbind(df03, temp)
  }
  df03$stage <- factor(x = df03$stage, levels = prex)
  p_01 <- ggplot(data = df03) +
    geom_line(mapping = aes(x = stage, y = cpmExp, group = geneName), color = "steelblue", linewidth = 1) +
    geom_point(mapping = aes(x = stage, y = cpmExp), color = "steelblue", size = 1.5) +
    geom_errorbar(mapping = aes(x = stage, ymin = cpmExp-cpmSD, ymax = cpmExp+cpmSD), 
                  width = 0.15, color = "steelblue", linewidth = 0.8) +
    labs(x = "", y = "mean CPM", title = paste0(demo," (","GSE45719",")")) +
    #facet_wrap(.~ geneName, scales = "free_y", nrow = 4) +
    scale_x_discrete(labels = c("Oocyte","Zygote","2-cell Early","2-cell Mid",
                                "2-cell Late","4-cell","8-cell","16-cell",
                                "Blastocyst Early")) +
    theme_classic() +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(size = 16, color = "black", face = "bold"),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          strip.text.x = element_text(size = 12, color = "red", face = "bold"),
          strip.background = element_blank(),
          panel.spacing.y = unit(0.8, "cm"))
  p_01
  ggsave(plot = p_01, 
         filename = paste0("/ChIP_seq_2/aaSHY/Prmt1/RNAseq/202306/PLOT/",demo,
                           "-Embryogenesis-CPM.pdf"),
         units = "cm", height = 15, width = 15)
}
################################################################################
#
#28、Fbl-KD KEGG、GO-bp、GO-mf
{
  {
    res1 <- as.data.frame(read.csv(paste0(path,"DESeq2/fbl.geneOnly.csv"), 
                                   header = T, stringsAsFactors = F))
    upup <- res1[which(res1$log2FoldChange > log2(1.5)  &  res1$padj < 0.05),"X"]
    down <- res1[which(res1$log2FoldChange < -log2(1.5) &  res1$padj < 0.05),"X"]
    for (symb in seq(2)) {
      gen1 <- unlist(list(upup,down)[symb])
      geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(gen1), keytype="SYMBOL", column="ENTREZID"))
      #
      kegg <- enrichKEGG(gene = geyy,
                         organism = "mmu", keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      kegg <- kegg@result
      kegg$Description <- sapply(str_split(kegg$Description, pattern = " - Mus"), "[", 1)
      assign(paste0("keg",symb), kegg)
      #
      gobp <- clusterProfiler::enrichGO(gene = geyy, ont = "BP", readable = T,
                                        keyType = "ENTREZID", OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05,qvalueCutoff = 0.05)
      gobp <- gobp@result
      assign(paste0("gbp",symb), gobp)
      #
      gomf <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, readable = T,
                                        keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05,qvalueCutoff = 0.05)
      gomf <- gomf@result
      assign(paste0("gmf",symb), gomf)
    }
  }
  #upup
  {
    keg1$Description <- factor(keg1$Description, levels = rev(keg1$Description))
    p_01 <- ggplot(keg1[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      theme_few() +
      scale_x_continuous(breaks = seq(12)) +
      geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      labs(y="",title = "KEGG\ngenes: Fbl KD deRegulated") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_01
    #
    gbp1$Description <- factor(gbp1$Description, levels = rev(gbp1$Description))
    p_02 <- ggplot(gbp1[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      theme_few() +
      scale_x_continuous(breaks = seq(1,28, by=2)) +
      geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 15), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 17), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 19), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 21), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 23), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 25), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 27), color = "white", linewidth  = 0.8) +
      labs(y="",title = "GO BP\ngenes: Fbl KD deRegulated") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_02
    gmf1$Description <- factor(gmf1$Description, levels = rev(gmf1$Description))
    p_03 <- ggplot(gmf1[1:10,])+
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
      geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 12), color = "white", linewidth  = 0.8) +
      labs(y="",title = "GO MF\ngenes: Fbl KD deRegulated") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_03
    plot <- p_01 +p_02 +p_03 +plot_layout(nrow = 3)
    ggsave(plot = plot, 
           units = "cm", width = 24, height = 24,
           filename = paste0(path,"PLOT/Enrichment/fbl.kegg.gobp.gpmf.upup.pdf"))
  }
  #down
  {
    keg2$Description <- factor(keg2$Description, levels = rev(keg2$Description))
    p_04 <- ggplot(keg2[2:11,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      theme_few() +
      scale_x_continuous(breaks = seq(12)) +
      geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      labs(y="",title = "KEGG\ngenes: Fbl KD deRegulated") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_04
    #
    gbp2$Description <- factor(gbp2$Description, levels = rev(gbp2$Description))
    p_05 <- ggplot(gbp2[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      theme_few() +
      scale_x_continuous(breaks = seq(1,24, by=2)) +
      geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 15), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 17), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 19), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 21), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 23), color = "white", linewidth  = 0.8) +
      labs(y="",title = "GO BP\ngenes: Fbl KD deRegulated") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_05
    gmf2$Description <- factor(gmf2$Description, levels = rev(gmf2$Description))
    p_06 <- ggplot(gmf2[1:10,])+
      geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
      theme_few() +
      scale_x_continuous(breaks = seq(20)) +
      geom_vline(aes(xintercept = 1), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 2), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 3), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 4), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 6), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 7), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 8), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 9), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 10), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 11), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 12), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 13), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 14), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 15), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 16), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 17), color = "white", linewidth  = 0.8) +
      geom_vline(aes(xintercept = 18), color = "white", linewidth  = 0.8) +
      labs(y="",title = "GO MF\ngenes: Fbl KD deRegulated") +
      theme(text = element_text(size = 12,family = "serif",color = "black"),
            panel.border = element_rect(linewidth = 1),
            plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
            axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
            axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
    p_06
    plot <- p_04 +p_05 +p_06 +plot_layout(nrow = 3)
    ggsave(plot = plot, 
           units = "cm", width = 24, height = 24,
           filename = paste0(path,"PLOT/Enrichment/fbl.kegg.gobp.gpmf.down.pdf"))
  }
}
################################################################################
#29、python GSEA
{
  temp <- read.csv(paste0(path, "DESeq2/fbl.geneOnly.csv"), header = T)
  temp <- temp[, c("X", "log2FoldChange")]
  temp <- temp[order(temp$log2FoldChange, decreasing = T),]
  write.table(x = temp, file = paste0(path, "DESeq2/fbl.res.rnk"), quote = F, sep = "\t", col.names = F, row.names = F)
  cmd_12 <- paste0("gseapy prerank -g ","/Reference/aaSHY/GSEAgmt/TwoGenes.gmt -r ",
                   paste0(path, "DESeq2/fbl.res.rnk"), " -p 10 -o ", 
                   paste0(path, "PLOT/EnrichmentGSEA/python/"))#conda activate base gseapy 0.10.8
}
################################################################################
#30、co-regulated GSEA









