#write by ~~~ at 2024.01.23
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2",
              "clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/RN_try_02/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/mervl_zyw.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inxA = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/allGene.bed"
  inxB = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bowtie2/Rn45s"
  inxD = "/Reference/aaSHY/zOther/TEconsensus/MuERVL/bt2/MuERVL"
}
################################################################################
#03、跑命令
{
  for (i in c("src","trim","Log","peak","dumpROOM","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }; rm(i)
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("IP","Input")
  cmd_02 <- paste0("trim_galore --phred33 --fastqc --illumina ",path,
                   "fq/", prex, ".fq.gz -o ", path, "trim", "\n")
  cmd_04 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inx1," -U ", path, "trim/", prex, "_trimmed.fq.gz",
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_05 <- paste0("run-csem -p 10 --no-extending-reads --bam ",path,
                   "bam/",prex,".bam"," 140 ",path,
                   "bam/csem/",prex,"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",path,
                   "bam/csem/",prex,".sorted.bam -o ",path,
                   "bw/bcv/",prex,".bw","\n")
  cmd_07 <- paste0("bamCompare -p 10 --scaleFactorsMethod SES ",
                   "--sampleLength 100000 -b1 ",path,"bam/csem/",prex[1],
                   ".sorted.bam -b2 ",path,"bam/csem/",prex[2],
                   ".sorted.bam -o ", path,"bw/bcp/",prex[1],"_Input.bw","\n")
  cmd_08 <- paste0("bamCompare -p 10 ",
                   "-b1 ",path,"bam/csem/",prex[1],
                   ".sorted.bam -b2 ",path,"bam/csem/",prex[2],
                   ".sorted.bam -o ", path,"bw/bcp/",prex[1],"_Input-notSES.bw","\n")
  bw01 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-1_Input-notSES.bw"
  bw02 <- paste0(path,"bw/bcp/",prex[1],"_Input-notSES.bw")
  bw03 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/KO/bw/bcp/IP-1_Input-notSES.bw"
  tags <- c("WT","RN","KO")
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center ",
                   "-a 3000 -b 3000 --missingDataAsZero -R ",c(inx7,inx8)," -S ",
                   paste0(c(bw01,bw02,bw03),collapse = " "),
                   " -o ",path,"PLOT/dddd-",c("01","02"),".mat.gz","\n")
  cmd_10 <- paste0("plotProfile --perGroup -m ",path,
                   "PLOT/dddd-",c("01","02"),".mat.gz ",
                   "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                   "--samplesLabel ",paste0(tags,collapse = " ")," ",
                   "-o ",path,"PLOT",
                   "/WT_RN_KOH3K9me3-",gsub(".bed","",basename(c(inx7,inx8))),".pdf","\n")
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_10) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}






















