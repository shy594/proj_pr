#write by ~~~ at 2024.01.23
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/IAPEY_LTR::ERVK::LTR.bed"
  inx9 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/RLTR1B-int::ERV1::LTR.bed"
  inxA = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/Reference/INDEX/mm10/mm10"
  inxB = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/Reference/INDEX/MuERVL/bt/MuERVL"
  inxC = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/Reference/INDEX/MT2_Mm/bt/MT2_Mm"
  inxD = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/Reference/INDEX/IAPEY_LTR/bt/IAPEY_LTR"
  inxE = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/Reference/INDEX/RLTR1B-int/bt/RLTR1B"
  inxF = "/ChIP_seq_2/aaSHY/pr/publicData/2OMe/Reference/INDEX/rRNA_v2/bt/rRNA_v2"
}
##
##
##
#03、构建索引【弃用 已完成】
{
  fa01 <- list.files(path = "/Reference/aaSHY/zOther/TEconsensus/fa", pattern = "fai", recursive = T, full.names = T)
  fa01 <- gsub(x = fa01, pattern = ".fai", replacement = "")
  prex <- dir(path = "/Reference/aaSHY/zOther/TEconsensus/fa/")
  prex <- prex[grep(x = prex, pattern = "fa", invert = T)]
  for (i in prex) {
    if (dir.exists(paths = paste0(path,"Reference/fa/", i))==F) {
      dir.create(path = paste0(path,"Reference/fa/", i), recursive = T)
    }
  }
  for (i in prex) {
    if (dir.exists(paths = paste0(path,"Reference/INDEX/", i,"/bt"))==F) {
      dir.create(path = paste0(path,"Reference/INDEX/", i,"/bt"), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("cat ",fa01," ",path,"Reference/fa/spikein.fa |grep -v -P ^$ >",path,
                   "Reference/fa/",prex,"/",prex,".fa","\n")
  cmd_02 <- paste0("samtools faidx ",path,"Reference/fa/",prex,"/",prex,".fa","\n")
  cmd_03 <- paste0("cut -f 1,2 ",path,"Reference/fa/",prex,"/",prex,".fa.fai >",path,
                   "Reference/fa/",prex,"/chrom.sizes.txt","\n")
  cmd_04 <- paste0("bowtie-build ", path,"Reference/fa/",prex,"/",prex,".fa", " ",
                   path,"Reference/INDEX/", prex,
                   "/bt/",gsub(x = prex, pattern = "-int", replacement = ""),"\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
##
##
##
#04、跑命令
{
  for (i in c("src","trim","bam","Log","PLOT","bamCoverage","zoomROOM")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
    rm(i)
    shell <- paste0(path,"src/run_a.sh")
    cat("#!/bin/bash\n", file = shell)
    frex <- list.files(path = paste0(path, "rawdata"), recursive = T, full.names = T)
    frex <- frex[grep(frex, pattern="other|txt|mESC|adaptor", invert = T)]
    cmd_01 <- paste0("#parallel-fastq-dump -t 10 --split-3 -s ",frex,
                     " -O ", path, "fq", "\n")
    cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*fastq", "\n")
    prex <- gsub(x = basename(frex), pattern = ".sra", replacement = "")
    cmd_03 <- paste0("#fastx_trimmer -f 6 -i ", path, "fq/", prex[1:2],".fq |",
                     "fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 25 -M 10 -Q 33 ",
                     "-o ",path,"trim/",prex[1:2],"_v2.fq","\n")
    cmd_04 <- paste0("#fastx_trimmer -f 6 -i ", path, "fq/", prex[3],".fq |",
                     "fastx_clipper -a AGATCGGAAGAGCACACGTCT -l 25 -M 10 -Q 33 ",
                     "-o ",path,"trim/",prex[3],"_v2.fq","\n")
    for (i in cmd_01) {cat(i, append = T, file = shell)}
    for (i in cmd_02) {cat(i, append = T, file = shell)}
    for (i in cmd_03) {cat(i, append = T, file = shell)}
    for (i in cmd_04) {cat(i, append = T, file = shell)}
    for (inxY in c(inxB,inxC,inxD,inxE,inxF)){
      cmd_05 <- paste0("#bowtie -q -n 2 --norc -m 1 --best --strata -S -x ",
                       inxY," ",path,"trim/",prex,"_v2.fq |samtools view -F 4 -bS |",
                       "samtools sort >",path,"bam/",basename(inxY),"_",prex,".bam","\n")
      cmd_06 <- paste0("#samtools index ",path,"bam/",basename(inxY),"_",prex,".bam","\n")
      cmd_07 <- paste0("#genomeCoverageBed -bga -ibam ",path,
                       "bam/",basename(inxY),"_",prex,".bam -g ",path,
                       "Reference/fa/",
                       sapply(str_split(inxY, pattern = "/"),"[",9),
                       "/chrom.sizes.txt >",path,
                       "bamCoverage/",basename(inxY),"_",prex,".bedGraph", "\n")
      for (i in cmd_05) {cat(i, append = T, file = shell)}
      for (i in cmd_06) {cat(i, append = T, file = shell)}
      for (i in cmd_07) {cat(i, append = T, file = shell)}
    };rm(inxY)
    cmd_08 <- paste0("bowtie -q -n 2 --norc -m 1 --best --strata -S -x ",
                     inxA," ",path,"trim/",prex,"_v2.fq |samtools view -F 4 -bS |",
                     "samtools sort >",path,"bam/",basename(inxA),"_",prex,".bam","\n")
    cmd_09 <- paste0("samtools index ",path,"bam/",basename(inxA),"_",prex,".bam","\n")
    cmd_10 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     path,"bam/",basename(inxA),"_",prex,".bam -o ",
                     path,"bamCoverage/",basename(inxA),"_",prex,".bw","\n")
    for (i in cmd_08) {cat(i, append = T, file = shell)}
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}


