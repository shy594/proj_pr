#Write by ~~~ on 2024.02.15
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","ggforce","bedtoolsr","refGenome","yyplot","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
################################################################################
################################################################################
#############################RADICL-seq FA 2%###################################
################################################################################
################################################################################
#RADICL-seq FA 2%
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/"
  prex <- c("n1","n2","n3")
  bath <- paste0(path, "specialOut/")
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inxA = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/bwa/MuERVL"#bwa index -a is -p STR ref.fa
  inxB = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MT_Mm/bwa/MT"#bwa index -a is -p STR ref.fa
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bwa/Rn45s"#bwa index -a is -p STR ref.fa
}
################################################################################
#03、MM RNA结合的DNA tag在MM、MT2_Mm consensus序列上的map结果, 在IGV展示
for (n in prex) {
  ff01 <- read.table(file = paste0(bath, n, ".MM.dna.bed"), sep = "\t", header = F, stringsAsFactors = F)
  ff01$V6 <- gsub(x = ff01$V6, pattern = "n.*BNA", replacement = paste0(n,"DNA"))
  write.table(x = ff01$V6, 
              file = paste0(path, "tagOnConsensusTE/readsName.", n,".MMtaggeddna.txt"), quote = F, sep = "\t", col.names = F, row.names = F)
} #已完成，弃用
{
  shell <- paste0(path,"tagOnConsensusTE/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("seqtk subseq ", path, "fq/", prex, "DNA.s2.fq ",path,
                   "tagOnConsensusTE/readsName.", prex,".MMtaggeddna.txt ",
                   ">", path, "tagOnConsensusTE/seqtk.",prex,"DNA.s2.fq", "\n")
  for (m in cmd_01) {cat(m, file = shell, append = T)}
  for (x in c(inxA, inxB)) {
    cmd_02 <- paste0("bwa aln -t 10 -f ", path,"tagOnConsensusTE/", prex,basename(x),".s2.sai ", x, " ",
                     path, "tagOnConsensusTE/seqtk.",prex,"DNA.s2.fq", "\n")
    cmd_03 <- paste0("bwa samse ", x, " ", path, "tagOnConsensusTE/", prex,basename(x),".s2.sai ", 
                     path, "tagOnConsensusTE/seqtk.", prex, "DNA.s2.fq |samtools view -bS -F 4 -@ 10 |",
                     "samtools sort -@ 10 >", path, "tagOnConsensusTE/", prex,basename(x),".s2.bam","\n")
    cmd_04 <- paste0("samtools index ",path, "tagOnConsensusTE/", prex,basename(x),".s2.bam", "\n")
    cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     path, "tagOnConsensusTE/", prex,basename(x),".s2.bam -o ",
                     path, "tagOnConsensusTE/", prex,basename(x),".s2.bw", "\n")
    for (m in cmd_02) {cat(m, file = shell, append = T)}
    for (m in cmd_03) {cat(m, file = shell, append = T)}
    for (m in cmd_04) {cat(m, file = shell, append = T)}
    for (m in cmd_05) {cat(m, file = shell, append = T)}
  }
  print(paste0("nohup bash ", shell, " >", paste0(path,"tagOnConsensusTE/"), basename(shell), ".log", " 2>&1 &"))
}
################################################################################
################################################################################
#################################GRID-seq#######################################
################################################################################
################################################################################
################################################################################
#GRID-seq
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  prex <- c("n1","n2")
  bath <- paste0(path, "specialOut/")
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inxA = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/bwa/MuERVL"#bwa index -a is -p STR ref.fa
  inxB = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MT2_Mm/bwa/MT2"#bwa index -a is -p STR ref.fa
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bwa/Rn45s"#bwa index -a is -p STR ref.fa
  inxD = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
}
################################################################################
#03、MM RNA结合的DNA tag在MM、MT2_Mm consensus序列上的map结果, 在IGV展示
for (n in prex) {
  ff01 <- read.table(file = paste0(bath, n, ".MM.dna.bed"), sep = "\t", header = F, stringsAsFactors = F)
  ff01$V6 <- gsub(x = ff01$V6, pattern = "abBNA", replacement = paste0("rep",substr(n,2,2),"gdna"))
  write.table(x = ff01$V6, 
              file = paste0(path, "tagOnConsensusTE/readsName.", n,".MMtaggeddna.txt"), quote = F, sep = "\t", col.names = F, row.names = F)
} #已完成，弃用
{
  shell <- paste0(path,"tagOnConsensusTE/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("#seqtk subseq ", path, "fq/", c("rep1","rep2"), "gdna.s2.fq ",path,
                   "tagOnConsensusTE/readsName.", prex,".MMtaggeddna.txt ",
                   ">", path, "tagOnConsensusTE/seqtk.",prex,"DNA.s2.fq", "\n")
  for (m in cmd_01) {cat(m, file = shell, append = T)}
  for (x in c(inxA, inxB)) {
    cmd_02 <- paste0("bwa aln -t 10 -f ", path,"tagOnConsensusTE/", prex,basename(x),".s2.sai ", x, " ",
                     path, "tagOnConsensusTE/seqtk.",prex,"DNA.s2.fq", "\n")
    cmd_03 <- paste0("bwa samse ", x, " ", path, "tagOnConsensusTE/", prex,basename(x),".s2.sai ", 
                     path, "tagOnConsensusTE/seqtk.", prex, "DNA.s2.fq |samtools view -bS -F 4 -@ 10 |",
                     "samtools sort -@ 10 >", path, "tagOnConsensusTE/", prex,basename(x),".s2.bam","\n")
    cmd_04 <- paste0("samtools index ",path, "tagOnConsensusTE/", prex,basename(x),".s2.bam", "\n")
    cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     path, "tagOnConsensusTE/", prex,basename(x),".s2.bam -o ",
                     path, "tagOnConsensusTE/", prex,basename(x),".s2.bw", "\n")
    for (m in cmd_02) {cat(m, file = shell, append = T)}
    for (m in cmd_03) {cat(m, file = shell, append = T)}
    for (m in cmd_04) {cat(m, file = shell, append = T)}
    for (m in cmd_05) {cat(m, file = shell, append = T)}
  }
  print(paste0("nohup bash ", shell, " >", paste0(path,"tagOnConsensusTE/"), basename(shell), ".log", " 2>&1 &"))
}






