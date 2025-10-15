#Write by ~~~ on 2023.11.09
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
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa.fai"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inx9 = "/Reference/aaSHY/BED/special/spliceF.txt"
}
##
##
##
##
##
##
#04、MM Interact with MM
te01 <- rtracklayer::import(con = inx6, format = "bed")#inx6、inx8
te02 <- rtracklayer::import(con = inx3, format = "gtf")
ge01 <- rtracklayer::import(con = inx2, format = "gtf")
ge01 <- ge01[which(ge01$type=="gene")]
#
for (i in c("n1","n2","n3")) {
  rd01 <- fread(file = paste0(path,"bed/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
  #取MM与RNA tags的交集
  rd02 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V2", start.field = "V3", end.field = "V4", strand.field = "V6", keep.extra.columns = T)
  ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = te01, subject = rd02, ignore.strand = T))
  rd02 <- as.data.frame(rd02[unique(subjectHits(ov02))])
  colnames(rd02) <- c("rnaChr","rnaSta","rnaEnd","rnaWidth","rnaStrand","ID","score","dnaChr","dnaSta","dnaEnd","dnaWidth","dnaStrand")
  #取MM与DNA tags的交集
  rd03 <- makeGRangesFromDataFrame(df = rd02, seqnames.field = "dnaChr", start.field = "dnaSta", end.field = "dnaEnd", strand.field = "dnaStrand", keep.extra.columns = T)
  ov03 <- suppressWarnings(GenomicRanges::findOverlaps(query = te01, subject = rd03, ignore.strand = T))
  rd03 <- as.data.frame(rd03[unique(subjectHits(ov03))])
  #标识哪些是MM RNA-MM DNA的Interaction?
  wr01 <- rd02
  wr01$maSELF <- NA; wr01$maSELF[which(wr01$ID %in% rd03$ID)] <- "YES"
  wr01$maSELF[which(!(wr01$ID %in% rd03$ID))] <- "NOO"
  print(as.data.frame(table(wr01$maSELF))); print(sum(as.data.frame(table(wr01$maSELF))[,2]))
  wr01$maSTRAND <- NA; wr01$maSTRAND[which(wr01$rnaStrand == wr01$dnaStrand & wr01$maSELF == "YES")] <- "YES"
  wr01$maSTRAND[which(wr01$rnaStrand != wr01$dnaStrand & wr01$maSELF == "YES")] <- "NOO"
  print(as.data.frame(table(wr01$maSTRAND))); print(sum(as.data.frame(table(wr01$maSTRAND))[,2]))
  write.table(x = wr01, file = paste0(path, "MMSELF/",i,".rnaDNA.bed"), sep = "\t", quote = F, col.names = T, row.names = F)
  }
##
##
##
#05、Ssrp1 strand specific
{
  path = "/ChIP_seq_1/aaSHY/rnadna/Ssrp1Strand/"
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("CT1","CT2")
  lay1 <- list.files(path = paste0(path,"rawdata"), full.names = T)
  cmd_01 <- paste0("#parallel-fastq-dump -t 20 --split-3 --gzip -s ", lay1, 
                   " -O ",path,"fq/", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ",path,"fq/*","\n")
  lay2 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay3 <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_03 <- paste0("#trim_galore --phred33 --fastqc -j 2 --illumina ","--clip_R1 10 ",
                   "--clip_R2 10 --paired ", lay2," ",lay3, " -o ", path, "trim", "\n")
  lay2 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay3 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_04 <- paste0("#STAR --runThreadN 10 --readFilesCommand zcat --outSAMstrandField intronMotif ",
                   "--outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate ",
                   "--outFileNamePrefix ",paste0(path,"bam/",prex)," ",
                   "--outSAMprimaryFlag AllBestScore --outSAMmultNmax -1 ",
                   "--genomeDir ",inx5," --readFilesIn ",lay2," ",lay3,"\n")
  cmd_05 <- paste0("#mv ",paste0(path,"bam/",prex,"Aligned.sortedByCoord.out.bam "),paste0(path,"bam/",prex,".bam"),"\n")
  cmd_06 <- paste0("#samtools index ",paste0(path,"bam/",prex,".bam"),"\n")
  cmd_07 <- paste0("stringtie -m 100 --rf ", paste0(path,"bam/",prex,".bam")," ",
                   "-p 20 -G ", inx2, " -o ", paste0(path,"stringTie/",prex,".gtf"),"\n")
  cmd_08 <- paste0("stringtie --merge -G ",inx2, " -o ", path,"stringTie/unin.gtf ",
                   paste0(path,"stringTie/",prex,".gtf", collapse = " "))
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  for (i in cmd_03) {cat(i, file = shell, append = T)}
  for (i in cmd_04) {cat(i, file = shell, append = T)}
  for (i in cmd_05) {cat(i, file = shell, append = T)}
  for (i in cmd_06) {cat(i, file = shell, append = T)}
  for (i in cmd_07) {cat(i, file = shell, append = T)}
  for (i in cmd_08) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
  path = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/"
}
##
##
##
#06、给MM RNA和MM DNA标上链的正负信息
#基于Ssrp1
ssrp <- rtracklayer::import(con = "/ChIP_seq_1/aaSHY/rnadna/Ssrp1Strand/stringTie/unin.gtf")
ssrp <- ssrp[ssrp$type=="transcript"]
#基于TE GTF
te06 <- rtracklayer::import(con = inx3); te06 <- te06[te06$gene_id=="MM-int"]
#基于Consensus Blast
for (i in c("n1","n2","n3")) {
  rd01 <- read.table(file = paste0(path, "MMSELF/",i,".rnaDNA.bed"), sep = "\t", stringsAsFactors = F, header = T)
  wr01 <- rd01[which(rd01$maSELF=="YES"),c(1,2,3,6)]
  write.table(x = wr01, file = paste0(path, "MMSELF/", i, ".toFA.bed"), sep = "\t", quote = F, col.names = F, row.names = F)
}
{
  prex = c("n1", "n2", "n3")
  shell <- paste0(path,"src/run_b.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("bedtools getfasta -fi ", inx1, " -fo ", path, "MMSELF/", prex, ".seq.fa ",
                   "-bed ", path, "MMSELF/", prex, ".toFA.bed -nameOnly", "\n")
  cmd_10 <- paste0("blastn -db /Reference/aaSHY/BLAST/Genome/MuERVL/MM -out ", path,
                   "MMSELF/", prex, ".blast.txt -query ", path, 
                   "MMSELF/", prex, 
                   ".seq.fa -outfmt 6 -task blastn-short -word_size 15", "\n")
  for (i in cmd_09) {cat(i, file = shell, append = T)}
  for (i in cmd_10) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"Log/"),basename(shell),".log", " 2>&1 &"))
}
#补充表达量信息 [pr ChIP-seq sample 1的MM的CPM]
cpm1 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20231108/singleEND/count/TEsites.cntTable", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
cpm1 <- cpm1[,c(5,1)]
colnames(cpm1) <- c("prIP","prInput")
cpm2 <- as.data.frame(apply(cpm1, 2, function(x){x/sum(x)*1000000}))
cpm2 <- cpm2[grep(x = rownames(cpm2), pattern = "MM-int"),]
rownames(cpm2) <- sapply(str_split(rownames(cpm2), pattern = ":"), "[", 1)
#
temp <- cpm1[grep(x = rownames(cpm1), pattern = "MM-int"),]
temp["MM-int_dup1627",]
#
#开始整合信息为表
for (i in c("n1","n2","n3")) {
  print(i)
  rd01 <- read.table(file = paste0(path, "MMSELF/",i,".rnaDNA.bed"), sep = "\t", stringsAsFactors = F, header = T)
  openxlsx::write.xlsx(x = rd01, file = paste0(path, "MMSELF/",i,".rnaDNA.xlsx"), asTable = T)
  rd01 <- rd01[which(rd01$maSELF=="YES"),]
  rd01$rnaBYssrp <- NA
  rd01$rnaBYgtff <- NA ; rd01$rnaBYcons <- NA
  rd01$dnaCPMofIP <- NA; rd01$dnaCPMofInput <- NA; rd01$rnaMMName <- NA; rd01$dnaMMName <- NA
  cons <- read.table(file = paste0(path, "MMSELF/", i, ".blast.txt"), sep = "\t", header = F, stringsAsFactors = F)
  cons$strand[cons$V10 - cons$V9 >0] = "+"; cons$strand[cons$V10 - cons$V9 <0] = "-"
  for (m in seq_along(rd01$rnaChr)) {
    tmp1 <- GenomicRanges::makeGRangesFromDataFrame(df = rd01[m,], seqnames.field = "rnaChr", start.field = "rnaSta", end.field = "rnaEnd", strand.field = "rnaStrand", keep.extra.columns = T)
    ov01 <- GenomicRanges::findOverlaps(query = tmp1, subject = ssrp, ignore.strand = T)
    rd01$rnaBYssrp[m] <- paste0(unique(as.character(ssrp[unique(subjectHits(ov01))]@strand)), collapse = ",")
    ov01 <- GenomicRanges::findOverlaps(query = tmp1, subject = te06, ignore.strand = T)
    rd01$rnaBYgtff[m] <- paste0(unique(as.character(te06[unique(subjectHits(ov01))]@strand)), collapse = ",")
    if (rd01$ID[m] %in% cons$V1) {
      rd01$rnaBYcons[m] <- cons[which(cons$V1==rd01$ID[m]),13]
    }
    rd01$rnaMMName[m] <- paste0(unique(as.character(te06[unique(subjectHits(ov01))]$transcript_id)), collapse = ",")
    #
    tmp2 <- GenomicRanges::makeGRangesFromDataFrame(df = rd01[m,], seqnames.field = "dnaChr", start.field = "dnaSta", end.field = "dnaEnd", strand.field = "dnaStrand", keep.extra.columns = T)
    ov02 <- GenomicRanges::findOverlaps(query = tmp2, subject = te06, ignore.strand = T)
    rd01$dnaMMName[m] <- paste0(unique(as.character(te06[unique(subjectHits(ov02))]$transcript_id)), collapse = ",")
    rd01$dnaCPMofIP[m] <- cpm2[rd01$dnaMMName[m],1]
    rd01$dnaCPMofInput[m] <- cpm2[rd01$dnaMMName[m],2]
  }
  wr02 <- rd01
  openxlsx::write.xlsx(x = wr02, file = paste0(path, "MMSELF/",i,".rnaMM.xlsx"), asTable = T)
}


######################################################################
twoc <- read.table(file = inx8, stringsAsFactors = F, sep = "\t")
twoc <- GenomicRanges::makeGRangesFromDataFrame(df = twoc, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3", strand.field = "V6")
ov03 <- GenomicRanges::findOverlaps(query = twoc, subject = ssrp, ignore.strand = T)
temp <- as.character(ssrp[unique(subjectHits(ov03))]@strand)
table(temp)







