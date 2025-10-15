#Write by ~~~ on 2023.09.04
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","ggforce","bedtoolsr","refGenome","yyplot","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inx5 = "/ChIP_seq_1/aaSHY/rnadna/RADICL/GenBank/rRNA.fa"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inx9 = "/Reference/aaSHY/BED/special/spliceF.txt"
  prex <- c("n1","n2","n3")
  bath <- paste0(path, "specialOut/")
}
##
##
##
#03、对于和MM RNA互作的DNA tags，统计每个TE上有多少数量，然后normalize
{
  bath <- paste0(path, "specialOut/")
  prex <- c("n1","n2","n3")
  tf01 <- rtracklayer::import(con = inx3)
  name <- unique(paste0(tf01$gene_id,"::",tf01$family_id,"::",tf01$class_id))
  gsiz <- read.table(file = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa.fai", sep = "\t", stringsAsFactors = F)
  gsiz <- sum(gsiz$V2)
  myta <- data.frame(teName = name, gsiz = gsiz)
  myta[,c("tsiz","n1","n2","n3")] <- NA
  myta[1,]
  for (i in myta$teName) {
    tf02 <- tf01[which(tf01$gene_id==str_split(i, pattern = "::")[[1]][1])]
    myta$tsiz[which(myta$teName==i)] <- sum(tf02@ranges@width)
    for (n in prex) {
      ff01 <- read.table(file = paste0(bath, n, ".MM.dtag.TE.txt"), sep = "\t", header = T, stringsAsFactors = F)
      if (i %in% ff01$Var1) {
        myta[which(myta$teName==i),n] <- ff01[which(ff01$Var1==i),"Freq"]
      }
      else {
        myta[which(myta$teName==i),n] <- 0
      }
    }
  }
  myta$nomz <- rowMeans(myta[,c("n1","n2","n3")])[1]/(myta$tsiz/myta$gsiz)
  tmp1 <- myta[order(myta$nomz, decreasing = T),]
  write.csv(x = tmp1, file = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/zemp/nomalizedNums.fa2.csv", quote = F, row.names = F)
}
##
##
##
#04、MM KD后的上调和下调的ERV中有哪些有MM RNA结合（总结合reads（n1+n2+n3）>100）?
{
  tmp1 <- read.csv(file = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/nomalizedNums.fa2.csv", stringsAsFactors = F)
  tmp1$repFamily <- sapply(str_split(tmp1$teName, pattern = "::"),"[",2)
  tmp2 <- na.omit(read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/tete_MM.csv", header = TRUE, stringsAsFactors = F,row.names=1))
  upup <- gsub(x = rownames(tmp2[which(tmp2$log2FoldChange >log2(1.5) & tmp2$pvalue <0.05),]), pattern = ":", replacement = "::")
  down <- gsub(x = rownames(tmp2[which(tmp2$log2FoldChange < -log2(1.5) & tmp2$pvalue <0.05),]), pattern = ":", replacement = "::")
  tmp1$deregulate <- NA
  tmp1$deregulate[which(tmp1$teName %in% upup)] <- "upupRegulated"
  tmp1$deregulate[which(tmp1$teName %in% down)] <- "downRegulated"
  tmp1$sumReads <- rowSums(x = tmp1[,c(4:6)])
  write.csv(x = tmp1, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/MMKD/DESeq2/tete_MMKD_RNA-DNA.csv", row.names = F, quote = F)
  #tmp1[grep(tmp1$teName, pattern="RLTR1B"),];tmp2[grep(rownames(tmp2), pattern="RLTR1B"),]
  #wr01 <- tmp1[!is.na(tmp1$deregulate),]
  #wr01 <- wr01[grep(wr01$repFamily, pattern="ERV"),]
  #wr01 <- wr01[order(wr01$sumReads, decreasing = T),]
  #View(wr01)
}
##
##
##
#05、有MERLV RNA结合的RLTR1B, RLTR4_MM-int, RLTR1B-int坐标可以给我下吗
{
  tf01 <- rtracklayer::import(con = inx3)
  tf02 <- tf01[grep(x = tf01$gene_id, pattern = "RLTR1B|RLTR4_MM-int|RLTR1B-int", ignore.case = F)]
  myta <- data.frame()
  for (n in c("n1","n2","n3")) {
    ff01 <- read.table(file = paste0(bath, n, ".MM.DNA.bed"), sep = "\t", header = F, stringsAsFactors = F)
    colnames(ff01) <- c("chr","start","end","name","score","strand")
    gf01 <- GenomicRanges::makeGRangesFromDataFrame(df = ff01, keep.extra.columns = T)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = tf02, subject = gf01, ignore.strand = T))
    df01 <- as.data.frame(tf02[queryHits(ov01)])
    temp <- data.frame(position = paste0(df01$seqnames,":",df01$start,"-",df01$end),
                       dupID = df01$transcript_id, sample = n,
                       teName = paste0(df01$gene_id, "::", df01$family_id,"::",df01$class_id))
    myta <- rbind(myta, temp)
  }
  write.csv(x = myta, file = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/sites.RLTR1B-MM.csv", quote = F, row.names = F)
}
##
##
##
#06、看MM RNA互作的DNA（+-3kb）上是否有H3K9me3富集
for (n in prex) {
  ff01 <- read.table(file = paste0(bath, n, ".MM.DNA.bed"), sep = "\t", header = F, stringsAsFactors = F)
  ff01$V2 <- ff01$V2-5000
  ff01$V3 <- ff01$V3+5000
  write.table(x = ff01, file = paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/", n,".expand5000bp.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
} #已完成，弃用
{
  shell <- paste0(path,"run_intersEnrich_v1.sh")
  cat("#!/bin/bash\n", file = shell)
  bath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bed1 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n1.noChr.expand3000bp.bed")
  bed2 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n2.noChr.expand3000bp.bed")
  bed3 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n3.noChr.expand3000bp.bed")
  cmd_09 <- paste0("computeMatrix reference-point -p 20 -a 3000 -b 3000 --referencePoint center ",
                   "--missingDataAsZero -R ", paste0(bed1," ",bed2), " -S ",
                   paste0(bath, c("WT"),"/bw/bcp/IP-1_Input.bw", collapse = " ")," ",
                   "-o ",path,"intersEnrich/fa2.n1n2.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --perGroup --colorList ",
                   "white,red -min -0.8 -max 0.8 -m ", path, "intersEnrich/fa2.n1n2.mat.gz ",
                   "--samplesLabel ", paste0(c("WT"), "IP-1_Input", collapse = " ")," ",
                   "-out ",path,"intersEnrich/fa2.n1n2.h3k9me3.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_intersEnrich_v1.log"," 2>&1 &"))
}
##
##
##
#看MM RNA互作的DNA（+-5kb）上是否有H3K9me2富集，数据来自GSE131014，和GSE177058
for (n in prex) {
  ff01 <- read.table(file = paste0(bath, n, ".MM.DNA.bed"), sep = "\t", header = F, stringsAsFactors = F)
  ff01$V2 <- ff01$V2-5000
  ff01$V3 <- ff01$V3+5000
  write.table(x = ff01, file = paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/", n,".expand5000bp.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
} #已完成，弃用
{
  shell <- paste0(path,"run_intersEnrich.sh")
  cat("#!/bin/bash\n", file = shell)
  bath = "/ChIP_seq_2/aaSHY/pr/publicData/H3K9me2/"
  bed1 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n1.expand5000bp.bed")
  bed2 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n2.expand5000bp.bed")
  bed3 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n3.expand5000bp.bed")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 3000 -b 3000 -m 10000 ",
                   "--missingDataAsZero -p 20 -R ", paste0(bed1," ",bed2," ",bed3), " -S ",
                   paste0(bath, c("condition1","condition2"),"/bw/bcp/IP-1_Input.bw", collapse = " ")," ",
                   paste0(bath, c("condition1","condition2"),"/bw/bcp/IP-2_Input.bw", collapse = " ")," ",
                   "-o ",path,"intersEnrich/fa2.n1n2n3.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "intersEnrich/fa2.n1n2n3.mat.gz ",
                   "--samplesLabel ", paste0(c("condition1","condition2"), "IP-1_Input", collapse = " ")," ",
                   paste0(c("condition1","condition2"), "IP-2_Input", collapse = " ")," ",
                   "-out ",path,"intersEnrich/fa2.n1n2n3.h3k9me2.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_intersEnrich.log"," 2>&1 &"))
}
##
##
##
#看MM RNA互作的DNA（+-3kb）上是否有H3K9me3, Setdb1*2, Kap1, Hdac1, G9a, H2Bub1的富集
for (n in prex) {
  ff01 <- read.table(file = paste0(bath, n, ".MM.DNA.bed"), sep = "\t", header = F, stringsAsFactors = F)
  #ff01$V1 <- paste0("chr",ff01$V1)
  ff01$V2 <- as.integer(ff01$V2-3000)
  ff01$V3 <- as.integer(ff01$V3+3000)
  write.table(x = ff01, 
              file = paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/", n,".noChr.expand3000bp.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
} #已完成，弃用
{
  shell <- paste0(path,"run_intersEnrich_v2.sh")
  cat("#!/bin/bash\n", file = shell)
  bed1 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n1.noChr.expand3000bp.bed")
  bed2 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n2.noChr.expand3000bp.bed")
  bed3 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n3.noChr.expand3000bp.bed")
  #tmp1 <- read.table(file = bed3, sep = "\t", header = F, stringsAsFactors = F)
  #min(tmp1$V3-tmp1$V2)
  #class(tmp1$V2);class(tmp1$V3)
  bw01 <- c("/ChIP_seq_2/aaSHY/Dot1l/shy_dot1ChIP/review/bw/Dot1lIP-Input.bw",
            "/ChIP_seq_2/aaSHY/Npm1/ChIP-seq/bw/bcp/Npm1IP-Input.bw",
            "/ChIP_seq_2/Ssrp1_ChIP_bp/bamCompare/D419_Ssrp1_mm10.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP007651_Kap1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP136128_Hdac1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP021942_G9a_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP068136_H2Bub1_pseudocount_1_log2.bw",
            "/ChIP_seq_2/aaSHY/Dot1l/shy_H3K79me3_dot1lReview/bw/H3K79me3vs.bw")
  for (i in bw01) {
    m=gsub(x = basename(i), pattern = "_pseudocount_._log2.bw", replacement = "")
    print(m)
    cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 3000 -b 3000 ",
                     "--missingDataAsZero -R ", paste0(bed1," ",bed2," ",bed3), " -S ",i," ",
                     "-o ",path,"intersEnrich/fa2.n1n2n3.",m,".mat.gz","\n")
    cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "intersEnrich/fa2.n1n2n3.",m,".mat.gz ",
                     "--samplesLabel ", m," ",
                     "-out ",path,"intersEnrich/fa2.n1n2n3.",m,".pdf","\n")
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_12) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_intersEnrich_v2.log"," 2>&1 &"))
}




