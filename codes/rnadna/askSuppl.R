#Write by ~~~ on 2024.04.12
#跑RNA-DNA的unique模式
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","ggforce","bedtoolsr","refGenome","yyplot","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inx5 = "/ChIP_seq_1/aaSHY/rnadna/RADICL/GenBank/rRNA.fa"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inx9 = "/Reference/aaSHY/BED/special/spliceF.txt"
}
###########################################################################
#03、跑命令
{
  shel <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shel)
  zemp <- list.files(path = paste0(path, "rawdata"), full.names = T)
  prex <- list.files(path = paste0(path, "rawdata"))
  cmd_01 <- paste0("#parallel-fastq-dump -s ",zemp," ",#--gzip 
                   "-t 10 --split-3 -O ", path,"fq/","\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path,"fq/*.fastq","\n")
  pair <- paste0(path,"zemp/pair.txt"); cat("", file = pair)
  for (m in paste0(path,"fq/",prex,".fq\n")) {cat(m, file =pair, append = T)}
  cmd_03 <- paste0("#cat ",pair," |head -n 2 >",path,"zemp/ff01.txt","\n")
  cmd_04 <- paste0("#cat ",pair," |head -n 4 |tail -n 2 >",path,"zemp/ff02.txt","\n")
  cmd_05 <- paste0("#cat ",pair," |tail -n 2 >",path,"zemp/ff03.txt","\n")
  cmd_06 <- paste0("#fastuniq -i ", path, "zemp/", c("ff01", "ff02", "ff03"), 
                   ".txt", " -t q -o ", path, "fq/n",c(1,2,3),"DNA.s1.fq", " -p ", path,
                   "fq/n",c(1,2,3), "RNA.s1.fq", "\n")
  cmd_07 <- paste0("#rRNAdust ", inx5, " ", path, "fq/", prex, ".s1.fq", " -o ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_08 <- paste0("#bwa aln -t 10 -f ",path,"bam/", prex, ".s2.sai ", inx4, " ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_09 <- paste0("#bwa samse ", inx4, " ", path, "bam/", prex, 
                   ".s2.sai ", path, "fq/", prex, ".s2.fq ",
                   "|samtools view -bS -F 4 -q 37 -@ 10 |samtools sort -@ 10 >", path, 
                   "bam/", prex, ".s2.bam","\n")
  cmd_10 <- paste0("#samtools index ",path, "bam/", prex, ".s2.bam", "\n")
  cmd_11 <- paste0("samtools view -@ 10 -q 37 -bS ", path, "bam/", prex, 
                   ".s2.bam ", "|bedtools bamtobed -i - >", path, "unique/", prex, ".s2.bed", "\n")
  cmd_12 <- paste0("cat ", path, "unique/", prex, ".s2.bed |sed 's/", prex, "/n1BNA/g' |",
                   "sort -k 4b,4 >", path, "unique/", prex, ".sort.bed", "\n")
  cmd_13 <- paste0("join -t $'\\t' -j 4 ", path, "unique/", prex[c(2,4,6)], ".sort.bed ", path, 
                   "unique/", prex[c(1,3,5)], ".sort.bed >", path, "unique/", c("n1", "n2", "n3"), ".join.bed", "\n")
  for (m in cmd_01) {cat(m, file = shel, append = T)}
  for (m in cmd_02) {cat(m, file = shel, append = T)}
  for (m in cmd_03) {cat(m, file = shel, append = T)}
  for (m in cmd_04) {cat(m, file = shel, append = T)}
  for (m in cmd_05) {cat(m, file = shel, append = T)}
  for (m in cmd_06) {cat(m, file = shel, append = T)}
  for (m in cmd_07) {cat(m, file = shel, append = T)}
  for (m in cmd_08) {cat(m, file = shel, append = T)}
  for (m in cmd_09) {cat(m, file = shel, append = T)}
  for (m in cmd_10) {cat(m, file = shel, append = T)}
  for (m in cmd_11) {cat(m, file = shel, append = T)}
  for (m in cmd_12) {cat(m, file = shel, append = T)}
  for (m in cmd_13) {cat(m, file = shel, append = T)}
  print(paste0("nohup bash ", shel, " >", paste0(path,"Log/"), basename(shel), ".log", " 2>&1 &"))
}
###########################################################################
#04、分析unique模式下的MM RNA-MM DNA [结合在原位还是异位]
{
  te01 <- rtracklayer::import(con = inx3)
  te02 <- rtracklayer::import(con = inx6, format = "bed")#inx6、inx8
  for (i in c("n1","n2","n3")) {
    rd01 <- fread(file = paste0(path,"unique/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
    my01 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V2", start.field = "V3", end.field = "V4", strand.field = "V6", keep.extra.columns = T)
    my02 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V7", start.field = "V8", end.field = "V9", strand.field = "V11", keep.extra.columns = T)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = te02, subject = my01, ignore.strand = T))
    qry1 <- cbind(as.data.frame(te02[queryHits(ov01)]),as.data.frame(my01[subjectHits(ov01)]))
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = te02, subject = my02, ignore.strand = T))
    qry2 <- cbind(as.data.frame(te02[queryHits(ov02)]),as.data.frame(my02[subjectHits(ov02)]))
    tmp1 <- qry1[which(qry1$V1 %in% qry2$V1),]
    tmp2 <- qry2[which(qry2$V1 %in% qry1$V1),]
    tmp2 <- tmp2[match(x = tmp1$V1, tmp2$V1),]
    tmp1$MMPOSr <- paste0("chr",tmp1[,1],":",tmp1[,2],"-",tmp1[,3])
    tmp2$MMPOSd <- paste0("chr",tmp2[,1],":",tmp2[,2],"-",tmp2[,3])
    tmp1$tagRNA <- paste0("chr",tmp1[,8],":",tmp1[,9],"-",tmp1[,10])
    tmp2$tagDNA <- paste0("chr",tmp2[,8],":",tmp2[,9],"-",tmp2[,10])
    wr01 <- cbind(tmp1[,c(13,21,11:12,6,20)],tmp2[,c(13,21,11:12,6,20)])
    colnames(wr01) <- c("rID","tagRNA","widthRNA","strandRNA","dupID-RNA","MMPOSr",
                        "dID","tagDNA","widthDNA","strandDNA","dupID-DNA","MMPOSd")
    wr01$diff <- ifelse(wr01$`dupID-RNA`==wr01$`dupID-DNA`,wr01$diff <- "same",wr01$diff <- "diff")
    openxlsx::write.xlsx(x = wr01, file = paste0(path,"unique/",i,".MM.rna.MM.dna.r.xlsx"))
  }
}
###########################################################################
#05、MM RNA tag在哪些TE上富集 [以TE consensus长度来标准化 tags number]
{
  cons <- read.table("/Reference/aaSHY/DatabaseFiles/RepeatMasker/TE_ID_fromUCSC.txt", sep = "\t", stringsAsFactors = F, header = F)
  myta <- data.frame()
  for (i in seq_along(cons$V1)) {
    if (cons$V10[i]=="-") {
      leng <- abs(cons$V14[i]) + abs(cons$V15[i])
    }
    else if (cons$V10[i]=="+") {
      leng <- abs(cons$V15[i]) + abs(cons$V16[i])
    }
    temp <- data.frame(name = paste0(cons$V11[i],"::",cons$V13[i],"::",cons$V12[i]),
                       leng = leng)
    myta <- rbind(myta, temp)
  }
  myta$gene <- sapply(str_split(myta$name, pattern="::"),"[",1)
  #write.csv(x = myta, 
  #          file = "/Reference/aaSHY/DatabaseFiles/RepeatMasker/TE_ID_fromUCSC.csv", 
  #          quote = F, row.names = F)
  ###
  tmp1 <- read.table(file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/n1.MM.dtag.TE.txt", header = T, stringsAsFactors = F, sep = "\t")
  tmp1$colo <- "other"
  tmp1$colo[grep(tmp1$Var1, pattern="MT2_Mm|MM-int")] <- "mark"
  tmp1$leng <- NA
  for (i in seq_along(tmp1$Var1)) {
    name <- str_split(tmp1$Var1[i], pattern = "::")[[1]][1]
    leng <- myta[which(myta$gene==name),2]
    tmp1$leng[i] <- leng
  }
  tmp1$norm <- (tmp1$Freq*1000)/tmp1$leng
  tmp1 <- tmp1[order(tmp1$norm, decreasing = T),]
  tmp1$Var1 <- factor(x = tmp1$Var1, levels = tmp1$Var1)
  p_01 <- ggplot(data = tmp1) +
    geom_point(mapping = aes(x = Var1, y = log2(norm), color=colo)) +
    geom_point(data = tmp1[which(tmp1$colo=="mark"),],mapping = aes(x = Var1, y = log2(norm))) +
    theme_classic() +
    labs(y="Log2 (Normalized Tags Number)", x="TEs (subfamily)", title = "This Is Title...") +
    theme(axis.title = element_text(family = "serif", size = 12),
          plot.title = element_text(family = "serif", size = 15),
          axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 0.01),
          axis.line.y = element_line(linewidth = 1.2)) +
    geom_text_repel(data = tmp1[which(tmp1$colo=="mark"),], 
                    mapping = aes(x = Var1, y = log2(norm), label = Var1))
  p_01
  tmp1[grep(tmp1$Var1,pattern="IAPEy-int|IAPEY_LTR|MM-int|RLTR1B-int"),]
  ggsave(plot = p_01, 
         filename = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/GRIDseq-Freq-Length-TE.mt2.pdf",
         units = "cm", width = 20, height = 16)
  #write.csv(x = tmp1, 
  #          file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/GRIDseq-Freq-Length-TE.csv", 
  #          quote = F, row.names = F)
  ###
  tmp1$class <- sapply(str_split(tmp1$Var1, pattern="::"),"[",3)
  tmp2 <- tmp1[grep(tmp1$class, pattern = "LTR"),]
  tmp2 <- tmp2[order(tmp2$norm, decreasing = T),]
  tmp2$Var1 <- factor(x = tmp2$Var1, levels = tmp2$Var1)
  p_02 <- ggplot(data = tmp2) +
    geom_point(mapping = aes(x = Var1, y = log2(norm), color=colo)) +
    geom_point(data = tmp2[which(tmp2$colo=="mark"),],mapping = aes(x = Var1, y = log2(norm))) +
    theme_classic() +
    labs(y="Log2 (Normalized Tags Number)", x="ERVs (subfamily)", title = "This Is Title...") +
    theme(axis.title = element_text(family = "serif", size = 12),
          plot.title = element_text(family = "serif", size = 15),
          axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 0.01),
          axis.line.y = element_line(linewidth = 1.2)) +
    geom_text_repel(data = tmp2[which(tmp2$colo=="mark"),], 
                    mapping = aes(x = Var1, y = log2(norm), label = Var1))
  p_02
  ggsave(plot = p_02, 
         filename = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/GRIDseq-Freq-Length-ERV.mt2.pdf",
         units = "cm", width = 20, height = 16)
  #write.csv(x = tmp2, 
  #          file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/GRIDseq-Freq-Length-ERV.csv", 
  #          quote = F, row.names = F)
}





