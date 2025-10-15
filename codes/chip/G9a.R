#write by ~~~ at 2024.06.18
#1、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/publicData/G9a/"
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
  inxD = "/Reference/aaSHY/zOther/TEconsensus/INDEX/MuERVL/bt2/MuERVL"
}
################################################################################
#3、跑命令
{
  for (i in c("src","fq","Log","peak","dumpROOM","DESeq2","rRNA/csem","mervl/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("IP","Input")
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 --gzip -s ",
                   list.files(path = paste0(path, "rawdata"), pattern = "sra", recursive = T, full.names = T),
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*gz", "\n")
  cmd_03 <- paste0("trim_galore --phred33 -j 4 ",
                   "-o ", path, "trim ", path, "fq/", prex, ".fq.gz", "\n")
  cmd_04 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 10 -x ",
                   inx1," -U ", path, "trim/", prex, "_trimmed.fq.gz",
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_05 <- paste0("run-csem -p 10 --no-extending-reads --bam ",
                   path,"bam/",prex,".bam"," 30 ",
                   path,"bam/csem/",prex,"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   path,"bam/csem/",prex,".sorted.bam -o ",
                   path,"bw/bcv/",prex,".bw\n")
  cmd_07 <- paste0("bamCompare -p 10 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   path,"bam/csem/",prex[1],".sorted.bam -b2 ",
                   path,"bam/csem/",prex[2],".sorted.bam -o ",
                   path,"bw/bcp/",prex[1],"_Input.bw","\n")
  cmd_08 <- paste0("bamCompare -p 10 -b1 ",
                   path,"bam/csem/",prex[1],".sorted.bam -b2 ",
                   path,"bam/csem/",prex[2],".sorted.bam -o ",
                   path,"bw/bcp/",prex[1],"_Input-notSES.bw","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (q in prex[1]) {
    cmd_09 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 5000 ",
                     "--missingDataAsZero ", "-R ", c(inx7,inx8,inx9,inxA,inxB), " -S ", path,
                     "bw/bcp/",q,"_Input.bw -o ", path,"deeptools/", 
                     basename(c(inx7,inx8,inx9,inxA,inxB)),q, ".mat.gz","\n")
    cmd_10 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),q,
                     ".mat.gz -out ",path,"PLOT/",
                     basename(c(inx7,inx8,inx9,inxA,inxB)),q,".pdf","\n")
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_10) {cat(i, append = T, file = shell)}
  }
  cmd_11 <- paste0("macs2 callpeak -f BAM -g mm -q 0.1 --keep-dup 1 -t ",
                   path,"bam/csem/",prex[1],".sorted.bam -c ",
                   path,"bam/csem/",prex[2],".sorted.bam -n ",
                   prex[1],"-N --outdir ",
                   path,"peak/", "\n")
  cmd_12 <- paste0("macs2 callpeak -f BAM -g mm -q 0.1 --broad --broad-cutoff 0.1 ",
                   "--keep-dup 1 -t ",path,"bam/csem/", prex[1],
                   ".sorted.bam -c ", path,"bam/csem/", prex[2],
                   ".sorted.bam -n ", prex[1],
                   "-B --outdir ", path,"peak/", "\n")
  cmd_13 <- paste0("TEtranscripts -t ",
                   paste0(path,"bam/csem/",prex[1],".sorted.bam", collapse = " ")," -c ",
                   paste0(path,"bam/csem/",prex[2],".sorted.bam", collapse = " ")," ",
                   "--GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project SETDB1 --outdir ",paste0(path,"count/"),"\n")
  cmd_14 <- paste0("TElocal -b ",
                   path,"bam/csem/",prex,".sorted.bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  for (i in cmd_11) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  for (kk in c("mervl","rRNA")) {
    if (kk == "mervl") {
      inxM <- inxD
    }
    else {
      inxM <- inxC
    }
    cmd_15 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 10 -x ",
                     inxM," -U ", path, "trim/", prex, "_trimmed.fq.gz |",
                     "samtools view -F 4 -b > ", path, kk, "/",prex,".bam","\n")
    cmd_16 <- paste0("run-csem -p 10 --no-extending-reads --bam ",
                     paste0(path,kk,"/",prex,".bam")," 30 ",
                     paste0(path,kk,"/csem/",prex),"\n")
    cmd_17 <- paste0("bamCoverage --normalizeUsing CPM -p 10 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     paste0(path,kk,"/csem/",prex,".sorted.bam -o "),
                     paste0(path,kk,"/",prex,".",kk,".bw\n"))
    for (i in cmd_15) {cat(i, append = T, file = shell)}
    for (i in cmd_16) {cat(i, append = T, file = shell)}
    for (i in cmd_17) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}