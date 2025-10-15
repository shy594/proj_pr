#write by ~~~~ at 2024.02.01
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/publicData/Suv39h1/"
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
#跑命令
{
  shell <- paste0(path,"run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- basename(list.files(path = paste0(path, "fq"), pattern = "gz", recursive = T))
  prex <- unique(gsub(pattern = "_..fq.gz", replacement = "", fixed = F, x = prex))
  one <- paste0(path,"fq/",prex,"_1.fq.gz")
  two <- paste0(path,"fq/",prex,"_2.fq.gz")
  cmd_01 <- paste0("#trim_galore --phred33 --illumina ",
                   "--paired ",one," ",two," -o ",paste0(path,"trim\n"))
  one <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  two <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_02 <- paste0("#bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inx1," -1 ", one, " -2 ",two," ",
                   "|samtools view -F 4 -b > ", path, "bam/",prex,".bam","\n")
  cmd_03 <- paste0("#run-csem -p 20 --no-extending-reads --bam ",paste0(path,"bam/",prex,".bam")," 40 ",paste0(path,"bam/csem/",prex),"\n")
  cmd_03 <- paste0("#samtools sort -@ 10 ",paste0(path,"bam/",prex,".bam"), " >",
                   paste0(path,"bam/",prex,".sorted.bam"), "\n", "samtools index ",paste0(path,"bam/",prex,".sorted.bam"), "\n")
  cmd_04 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  bed1 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n1.noChr.expand3000bp.bed")
  bed2 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n2.noChr.expand3000bp.bed")
  bed3 <- paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/intersEnrich/n3.noChr.expand3000bp.bed")
  cmd_05 <- paste0("computeMatrix reference-point -p 20 -a 3000 -b 3000 --referencePoint center ",
                   "--missingDataAsZero -R ", paste0(bed1," ",bed2, " ",bed3), " -S ",
                   paste0(path,"bw/bcv/",prex,".bw", collapse = " ")," ",
                   "-o ",path,"PLOT/fa2.n1n2n3.Suv39h1.mat.gz","\n")
  cmd_06 <- paste0("plotHeatmap --heatmapHeight 12 --perGroup --colorList ",
                   "white,red -m ", path,"PLOT/fa2.n1n2n3.Suv39h1.mat.gz ",
                   "--samplesLabel  HP1 Suv39h1 Suv39h2 ",
                   "-out ",path,"PLOT/fa2.n1n2n3.Suv39h1.pdf","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}

