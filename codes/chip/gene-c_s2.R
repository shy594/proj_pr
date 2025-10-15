#write by ~~~ at 2024.01.07
#1、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240102/single_v2/"
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
##
##
##
{
  for (i in c("src","fq","Log","peak","dumpROOM","DESeq2","rRNA/csem","mervl/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- sapply(str_split(list.files(path = paste0(path, "trim")), pattern = "_"), "[", 1)
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
  cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   paste0(path,"bam/csem/",prex[3:4],".sorted.bam -b2 "),
                   paste0(path,"bam/csem/",prex[1:2],".sorted.bam -o "),
                   paste0(path,"bw/bcp/",prex[3:4],"_Input.bw"),"\n")
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (q in prex[3:4]) {
    cmd_08 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                     "--missingDataAsZero -p 20 ", "-R ", c(inx7,inx8), " -S ", path,
                     "bw/bcp/",q,"_Input.bw -o ", path,"deeptools/", 
                     c(basename(inx7),basename(inx8)),q, ".mat.gz","\n")
    cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 500 ",
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
    for (i in cmd_08) {cat(i, append = T, file = shell)}
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    for (i in cmd_11) {cat(i, append = T, file = shell)}
    for (i in cmd_12) {cat(i, append = T, file = shell)}
  }
  cmd_13 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   path,"bam/csem/",prex[3:4],".sorted.bam -c ",
                   path,"bam/csem/",prex[1:2],".sorted.bam -n ",
                   prex[3:4],"-N --outdir ",
                   path,"peak/", "\n")
  cmd_14 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",path,"bam/csem/", prex[3:4],
                   ".sorted.bam -c ", path,"bam/csem/", prex[1:2],
                   ".sorted.bam -n ", prex[3:4],
                   "-B --outdir ", path,"peak/", "\n")
  cmd_15 <- paste0("TEtranscripts -t ",
                   paste0(path,"bam/csem/",prex[3:4],".sorted.bam", collapse = " ")," -c ",
                   paste0(path,"bam/csem/",prex[1:2],".sorted.bam", collapse = " "),
                   " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project ccPrmt --outdir ",paste0(path,"count/"),"\n")
  cmd_16 <- paste0("TElocal -b ",
                   path,"bam/csem/",prex,".sorted.bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  for (i in cmd_15) {cat(i, append = T, file = shell)}
  for (i in cmd_16) {cat(i, append = T, file = shell)}
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
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}


















