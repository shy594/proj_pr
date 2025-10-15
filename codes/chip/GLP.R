#write by ~~~ at 2024.03.15
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/publicData/Glp/"
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
#03、跑命令
{
  frex <- list.files(paste0(path, "rawdata"), full.names = T, pattern = "sra", recursive = T)
  prex <- c("IP","Input")
  for (i in c("src","trim","Log","peak","dumpROOM","DESeq2","rRNA/csem","mervl/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  rm(i)
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("#parallel-fastq-dump -t 10 --gzip --split-3 -s ",frex,
                   " -O ", path, "fq", "\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path, "fq/*fastq", "\n")
  #zcat merge IP-1 IP-2
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  for (i in cmd_02) {cat(i, file = shell, append = T)}
  cmd_03 <- paste0("trim_galore -j 2 --phred33 --fastqc --illumina ",path,
                   "fq/", prex, ".fq.gz -o ", path, "trim", "\n")
  cmd_04 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inx1," -U ", path, "trim/", prex, "_trimmed.fq.gz",
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_05 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"bam/",prex,".bam")," 140 ",
                   paste0(path,"bam/csem/",prex),"\n")
  cmd_06 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/csem/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (x in seq(1)) {
    zrex <- prex[c(x,x+1)]
    for (q in zrex[1]) {
      cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                       paste0(path,"bam/csem/",zrex[1],".sorted.bam -b2 "),
                       paste0(path,"bam/csem/",zrex[2],".sorted.bam -o "),
                       paste0(path,"bw/bcp/",zrex[1],"_Input.bw"),"\n")
      cmd_08 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 5000 ",
                       "--missingDataAsZero -R ", c(inx7,inx8), " -S ", path,
                       "bw/bcp/",q,"_Input.bw -o ", path,"deeptools/", 
                       c(basename(inx7),basename(inx8)),q, ".mat.gz","\n")
      cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 500 ",
                       "--missingDataAsZero -p 20 -R ", inx9, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inx9),q,".mat.gz","\n")
      cmd_10 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                       "--missingDataAsZero -p 20 -R ", inxA, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inxA),q,".mat.gz","\n")
      cmd_11 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 30000 ",
                       "--missingDataAsZero -p 20 -R ", inxB, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inxB),q,".mat.gz","\n")
      cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                       "deeptools/",c(basename(inx7),basename(inx8),basename(inx9),
                                      basename(inxA),basename(inxB)),q,
                       ".mat.gz -out ",path,"PLOT/",
                       c(basename(inx7),basename(inx8),
                         basename(inx9),basename(inxA),basename(inxB)),q,".pdf","\n")
      for (i in cmd_07) {cat(i, append = T, file = shell)}
      for (i in cmd_08) {cat(i, append = T, file = shell)}
      for (i in cmd_09) {cat(i, append = T, file = shell)}
      for (i in cmd_10) {cat(i, append = T, file = shell)}
      for (i in cmd_11) {cat(i, append = T, file = shell)}
      for (i in cmd_12) {cat(i, append = T, file = shell)}
    }
    rm(q)
    cmd_13 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                     path,"bam/csem/",zrex[1],".sorted.bam -c ",
                     path,"bam/csem/",zrex[2],".sorted.bam -n ",
                     zrex[2],"-N --outdir ",
                     path,"peak/", "\n")
    cmd_14 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                     "--keep-dup 1 -t ",path,"bam/csem/", zrex[1],
                     ".sorted.bam -c ", path,"bam/csem/", zrex[2],
                     ".sorted.bam -n ", zrex[1],
                     "-B --outdir ", path,"peak/", "\n")
    cmd_15 <- paste0("TEtranscripts -t ",
                     paste0(path,"bam/csem/",zrex[1],".sorted.bam", collapse = " ")," -c ",
                     paste0(path,"bam/csem/",zrex[2],".sorted.bam", collapse = " "),
                     " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                     "--project glp", zrex[1], " --outdir ",paste0(path,"count/"),"\n")
    cmd_16 <- paste0("TElocal -b ",
                     path,"bam/csem/",zrex,".sorted.bam ",
                     "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                     path, "count/", zrex,"\n")
    for (i in cmd_13) {cat(i, append = T, file = shell)}
    for (i in cmd_14) {cat(i, append = T, file = shell)}
    for (i in cmd_15) {cat(i, append = T, file = shell)}
    for (i in cmd_16) {cat(i, append = T, file = shell)}
  }
  rm(x)
  #print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
  #
  #
  #shell <- paste0(path,"src/run_b.sh")
  #cat("#!/bin/bash\n", file = shell)
  for (kk in c("mervl","rRNA")) {
    if (kk == "mervl") {
      inxM <- inxD
    }
    else {
      inxM <- inxC
    }
    cmd_17 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                     inxM," -U ", path, "trim/", prex, "_trimmed.fq.gz |",
                     "samtools view -F 4 -b > ", path, kk, "/",prex,".bam","\n")
    cmd_18 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                     paste0(path,kk,"/",prex,".bam")," 140 ",
                     paste0(path,kk,"/csem/",prex),"\n")
    cmd_19 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     paste0(path,kk,"/csem/",prex,".sorted.bam -o "),
                     paste0(path,kk,"/",prex,".",kk,".bw\n"))
    for (i in cmd_17) {cat(i, append = T, file = shell)}
    for (i in cmd_18) {cat(i, append = T, file = shell)}
    for (i in cmd_19) {cat(i, append = T, file = shell)}
  }
  rm(kk)
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
