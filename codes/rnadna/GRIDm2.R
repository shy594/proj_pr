#Write by ~~~ on 2023.09.11
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm2/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inx5 = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/GenBank/rRNA.fa"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
}
##
##
##
#03、跑命令
{
  shel <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shel)
  zemp <- list.files(path = paste0(path, "rawdata"), full.names = T)
  prex <- list.files(path = paste0(path, "rawdata"))
  cmd_01 <- paste0("parallel-fastq-dump -s ",zemp," ",#--gzip 
                   "-t 10 --split-3 -O ", path,"fq/","\n")
  cmd_02 <- paste0("rename s/fastq/fq/ ", path,"fq/*.fastq","\n")
  pair <- paste0(path,"zemp/pair.txt"); cat("", file = pair)
  for (m in paste0(path,"fq/",prex,".fq\n")) {cat(m, file =pair, append = T)}
  cmd_03 <- paste0("cat ",pair," |head -n 2 >",path,"zemp/ff01.txt","\n")
  cmd_04 <- paste0("cat ",pair," |tail -n 2 >",path,"zemp/ff02.txt","\n")
  cmd_05 <- paste0("fastuniq -i ", path, "zemp/", c("ff01", "ff02"), 
                   ".txt", " -t q -o ", path, "fq/rep",c(1,2),"gdna.s1.fq", " -p ", path,
                   "fq/rep",c(1,2), "cdna.s1.fq", "\n")
  cmd_06 <- paste0("rRNAdust ", inx5, " ", path, "fq/", prex, ".s1.fq", " -o ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_07 <- paste0("#bwa aln -t 10 -f ",path,"bam/", prex, ".s2.sai ", inx4, " ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_08 <- paste0("#bwa samse ", inx4, " ", path, "bam/", prex, ".s2.sai ", path, "fq/", prex, ".s2.fq ",
                   "|samtools view -bS -F 4 -@ 10 |samtools sort -@ 10 >", path, "bam/", prex, ".s2.bam","\n")
  cmd_09 <- paste0("#samtools index ",path, "bam/", prex, ".s2.bam", "\n")
  cmd_10 <- paste0("#samtools view -@ 10 -q 0 -bS ", path, "bam/", prex, 
                   ".s2.bam ", "|bedtools bamtobed -i - >", path, "bed/", prex, ".s2.bed", "\n")
  cmd_11 <- paste0("#cat ", path, "bed/", prex, ".s2.bed |sed 's/", prex, "/n1BNA/g' |",
                   "sort -k 4b,4 >", path, "bed/", prex, ".sort.bed", "\n")
  cmd_12 <- paste0("#join -t $'\\t' -j 4 ", path, "bed/", prex[c(2,4)], ".sort.bed ", path, 
                   "bed/", prex[c(1,3)], ".sort.bed >", path, "bed/", c("n1", "n2"), ".join.bed", "\n")
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
  print(paste0("nohup bash ", shel, " >", paste0(path,"Log/"), basename(shel), ".log", " 2>&1 &"))
}
#04、分析数据
##
##
##
te01 <- rtracklayer::import(con = inx6, format = "bed")#inx6、inx8
te02 <- rtracklayer::import(con = inx3, format = "gtf")
ge01 <- rtracklayer::import(con = inx2, format = "gtf")
ge01 <- ge01[which(ge01$type=="gene")]
#
for (i in c("n1","n2")) {
  rd01 <- fread(file = paste0(path,"bed/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
  rd01 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V2", start.field = "V3", end.field = "V4", strand.field = "V6", keep.extra.columns = T)
  ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = te01, subject = rd01, ignore.strand = T))
  rs01 <- as.data.frame(rd01[unique(subjectHits(ov01))])
  colnames(rs01) <- paste0("V", seq_along(colnames(rs01)))
  write.table(x = rs01, file = paste0(path, "specialOut/",i,".MM.dna.bed"), sep = "\t",
              quote = F, col.names = F, row.names = F)
  wr01 <- data.frame(chr = rs01$V8, start = rs01$V9, end = rs01$V10, name = ".", score = 1, strand = rs01$V12)
  write.table(x = wr01, file = paste0(path, "specialOut/",i,".MM.DNA.bed"), sep = "\t",
              quote = F, col.names = F, row.names = F)
  dtag <- unique(rs01[,c(6,8,9,10,12)])
  dtag <- makeGRangesFromDataFrame(df = dtag, seqnames.field = "V8", start.field = "V9",end.field = "V10", strand.field = "V12")
  #dtag <- dtag + 1000
  ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = te02, subject = dtag, ignore.strand = T))
  tf01 <- as.data.frame(te02[queryHits(ov02)])
  tf01$name <- paste0(tf01$gene_id, "::", tf01$family_id, "::", tf01$class_id)
  tf02 <- as.data.frame(table(tf01$name))
  tf02 <- tf02[order(tf02$Freq, decreasing = T),]
  write.table(x = tf02, file = paste0(path,"specialOut/",i,".MM.dtag.TE.txt"), sep = "\t", quote = F, row.names = F)
  tf01$name <- tf01$family_id
  tf02 <- as.data.frame(table(tf01$name))
  tf02 <- tf02[order(tf02$Freq, decreasing = T),]
  write.table(x = tf02, file = paste0(path,"specialOut/",i,".MM.dtag.family.txt"), sep = "\t", quote = F, row.names = F)
  tf01$name <- tf01$class_id
  tf02 <- as.data.frame(table(tf01$name))
  tf02 <- tf02[order(tf02$Freq, decreasing = T),]
  write.table(x = tf02, file = paste0(path,"specialOut/",i,".MM.dtag.class.txt"), sep = "\t", quote = F, row.names = F)
  rs02 <- makeGRangesFromDataFrame(df = rs01, keep.extra.columns = T, 
                                   seqnames.field = "V8", start.field = "V9", end.field = "V10", strand.field = "V12")
  ov03 <- suppressWarnings(GenomicRanges::findOverlaps(query = ge01, subject = rs02, ignore.strand = T))
  df01 <- as.data.frame(ge01[queryHits(ov03)])
  df02 <- as.data.frame(rs02[subjectHits(ov03)])
  fina <- cbind(df02, df01)
  colnames(fina) <- paste0("V", seq_along(colnames(fina)))
  out1 <- data.frame(rnaPos = paste0("chr",fina$V6,":",fina$V7,"-",fina$V8), 
                     rnaStrand = fina$V10,
                     dnaPos = paste0("chr",fina$V1,":",fina$V2,"-",fina$V3), 
                     dnaStrand = fina$V5, markID = fina$V11,
                     geneName = fina$V25, geneType = fina$V27, 
                     genePos = paste0("chr",fina$V14,":",fina$V15,"-",fina$V16))
  write.table(x = out1, file = paste0(path,"specialOut/",i,".MM.gene.txt"), sep = "\t", quote = F, row.names = F)
  inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
  wr02 <- out1[which(out1$geneName %in% inte$V4),]
  write.table(x = wr02, file = paste0(path,"specialOut/",i,".MM.2Clike.txt"), sep = "\t", quote = F, row.names = F)
  print(intersect(out1$geneName, inte$V4))
}
##2C-like 共有基因
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
twoc <- intersect(ff01$geneName, ff02$geneName)
ff01 <- ff01[which(ff01$geneName %in% twoc),]
ff02 <- ff02[which(ff02$geneName %in% twoc),]
write.table(x = ff01, file = paste0(path,"specialOut/n1.MM.2Clike.ov.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = ff02, file = paste0(path,"specialOut/n2.MM.2Clike.ov.txt"), sep = "\t", quote = F, row.names = F)
##统计TE
tf03 <- as.data.frame(te02)
se01 <- tf03 %>% dplyr::group_by(class_id) %>% summarise(length = sum(width), count = n())
se02 <- tf03 %>% dplyr::group_by(family_id) %>% summarise(length = sum(width), count = n())
se03 <- tf03 %>% dplyr::group_by(gene_id) %>% summarise(length = sum(width), count = n())
se01 <- se01[order(se01$length, decreasing = T), ]
se02 <- se02[order(se02$length, decreasing = T), ]
se03 <- se03[order(se03$length, decreasing = T), ]
write.table(x = se01, file = paste0(path,"specialOut/classID.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = se02, file = paste0(path,"specialOut/familyID.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = se03, file = paste0(path,"specialOut/repeatID.txt"), sep = "\t", quote = F, row.names = F)


