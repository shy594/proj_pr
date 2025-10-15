#write by ~~~ at 2024.01.23
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2",
              "clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
################################################################################
####目标
####pr、Fbl、WT、RN、KO H3K9me3、G9a、Trim28、Setdb1
####pr、Fbl需要在MERVL上的共同富集 [fig.3 F,G]
####pr、Fbl需要在2-cell genes的TSS上富集 [fig.3 H,I]
####pr、Fbl需要在Fbl peak上有富集 [fig.3 D]
####pr、Fbl需要在1165 MERVL target genes上富集 [fig.4 E,F]
####pr、Fbl、G9a、Trim28、Setdb1的相关性热图 [fig.5 A]
####WT H3K9me3在GRID-seq MERVL Binding sites的富集 [fig.6 A]
####WT H3K9me3、Trim28、Setdb1在1165 MERVL target genes上富集 [fig.6 B-D]
####WT、RN、KO H3K9me3在MERVL上富集 [fig.6 E]
{
  dirs <- c("pr","Fbl","H3K9me3_WT","H3K9me3_RN","H3K9me3_KO","G9a","Trim28","Setdb1")
  bed1 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bed/MERVL/mervl_zyw.bed"
  bed2 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bed/MERVL/MERVL-int::ERVL::LTR.bed"
  bed3 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bed/MERVLBindingSite/mervlBindingTags.bed"
  bed4 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bed/the1165/the1165.gene.bed"
  bed5 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bed/the2C/twoC.bed"
  bed6 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bed/theFblPeak/Fbl.narrowPeak"
  #
  #
  #
  bam1 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/bam/csem/IP-1.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/bam/csem/Input-1.sorted.bam")
  bam2 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/bam/csem/IP-1.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/bam/csem/Input-1.sorted.bam")
  bam3 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bam/csem/IP-1.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bam/csem/Input-1.sorted.bam")
  bam4 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/RN/bam/csem/IP-1.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/RN/bam/csem/Input-1.sorted.bam")
  bam5 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/KO/bam/csem/IP-1.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/KO/bam/csem/Input-1.sorted.bam")
  bam6 <- c("/ChIP_seq_2/aaSHY/pr/publicData/G9a/bam/csem/IP.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/publicData/G9a/bam/csem/Input.sorted.bam")
  bam7 <- c("/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/n1/bam/csem/TRIM28.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/n1/bam/csem/Input.sorted.bam")
  bam8 <- c("/ChIP_seq_2/aaSHY/pr/publicData/SETDB1/n2/bam/csem/IP.sorted.bam",
            "/ChIP_seq_2/aaSHY/pr/publicData/SETDB1/n2/bam/csem/Input.sorted.bam")
  bams <- list(bam1,bam2,bam3,bam4,bam5,bam6,bam7,bam8)
  #
  #
  #
  para_1 <- c("RPKM","CPM","BPM")
  para_2 <- c("readCount","SES")
  para_u <- c("RPKM","CPM","BPM","readCount","SES","SES_sL2","SES_sL2_derep_Q10")
  #
  #
  #
  ###shell <- paste0(path,"srcLog/run_d.sh")
  ###cat("#!/bin/bash\n", file = shell)
  ###for (i in seq_along(dirs)) {
    ###print(i)
    ###cmd_01 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 ",
    ###                 "-b1 ",bams[[i]][1]," ",
    ###                 "-b2 ",bams[[i]][2]," -o ",path,"bw/",dirs[i],
    ###                 "/IP_Input-SES.bw","\n")
    ###cmd_01 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 200000 ",
    ###                 "-b1 ",bams[[i]][1]," ",
    ###                 "-b2 ",bams[[i]][2]," -o ",path,"bw/",dirs[i],
    ###                 "/IP_Input-SES_sL2.bw","\n")
    ###cmd_01 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 200000 ",
    ###                 "--ignoreDuplicates --minMappingQuality 10 ",
    ###                 "-b1 ",bams[[i]][1]," ",
    ###                 "-b2 ",bams[[i]][2]," -o ",path,"bw/",dirs[i],
    ###                 "/IP_Input-SES_sL2_derep_Q10.bw","\n")
    ###cmd_02 <- paste0("bamCompare -p 20 --scaleFactorsMethod readCount ",
    ###                 "-b1 ",bams[[i]][1]," ",
    ###                 "-b2 ",bams[[i]][2]," -o ",path,"bw/",dirs[i],
    ###                 "/IP_Input-readCount.bw","\n")
    ###cmd_03 <- paste0("bamCompare -p 20 --scaleFactorsMethod None ",
    ###                 "--normalizeUsing ",para_1," -b1 ",bams[[i]][1]," ",
    ###                 "-b2 ",bams[[i]][2]," -o ",path,"bw/",dirs[i],
    ###                 "/IP_Input-",para_1,".bw","\n")
    ###for (i in cmd_01) {cat(i, append = T, file = shell)}
    ###for (i in cmd_02) {cat(i, append = T, file = shell)}
    ###for (i in cmd_03) {cat(i, append = T, file = shell)}
  ###}
  ###print(paste0("nohup bash ",shell," >",path,"srcLog/run_d.log"," 2>&1 &"))
  ##
  ##
  ##
  ##
  ##
  ##
  shell <- paste0(path,"srcLog/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  for (j in para_u) {
    for (l in c("pr","Fbl")) {
      cmd_04 <- paste0("#computeMatrix scale-regions -p 20 -m 5000 -a 5000 -b ",
                       "5000 -R ", c(bed1,bed2), " -S ",path,"bw/",l,"/",
                       "IP_Input-",j,".bw ",
                       "--skipZeros -o ",path,
                       "srcLog/aaaa-",c("01","02"),".mat.gz","\n")
      cmd_05 <- paste0("#plotProfile -m ",path,
                       "srcLog/aaaa-",c("01","02"),".mat.gz ",
                       "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                       "-o ",path,"thePlot/",j,"/",l,
                       "-",gsub(".bed","",basename(c(bed1,bed2))),".pdf","\n")
      for (i in cmd_04) {cat(i, append = T, file = shell)}
      for (i in cmd_05) {cat(i, append = T, file = shell)}
    }
    for (m in c("pr","Fbl")) {
      cmd_06 <- paste0("#computeMatrix scale-regions -p 20 -m 30000 -a 5000 -b ",
                       "5000 -R ", c(bed5), " -S ",path,"bw/",m,"/",
                       "IP_Input-",j,".bw ",
                       "--skipZeros -o ",path,
                       "srcLog/bbbb-",c("01"),".mat.gz","\n")
      cmd_07 <- paste0("#plotProfile -m ",path,
                       "srcLog/bbbb-",c("01"),".mat.gz ",
                       "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                       "-o ",path,"thePlot/",j,"/",m,
                       "-",gsub(".bed","",basename(c(bed5))),".pdf","\n")
      cmd_08 <- paste0("#computeMatrix scale-regions -p 20 -m 1000 -a 2500 -b ",
                       "2500 -R ", c(bed4), " -S ",path,"bw/",m,"/",
                       "IP_Input-",j,".bw ",
                       "--skipZeros -o ",path,
                       "srcLog/bbbb-",c("01"),".mat.gz","\n")
      cmd_09 <- paste0("#plotProfile -m ",path,
                       "srcLog/bbbb-",c("01"),".mat.gz ",
                       "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                       "-o ",path,"thePlot/",j,"/",m,
                       "-",gsub(".bed","",basename(c(bed4))),".pdf","\n")
      for (i in cmd_06) {cat(i, append = T, file = shell)}
      for (i in cmd_07) {cat(i, append = T, file = shell)}
      for (i in cmd_08) {cat(i, append = T, file = shell)}
      for (i in cmd_09) {cat(i, append = T, file = shell)}
    }
    cmd_10 <- paste0("#computeMatrix reference-point -p 20 --referencePoint center ",
                     "-a 1000 -b 1000 --missingDataAsZero -R ",c(bed6)," -S ",
                     path,"bw/Fbl/","IP_Input-",j,".bw ",
                     path,"bw/pr/","IP_Input-",j,".bw ",
                     "-o ",path,"srcLog/bbbb-",c("02"),".mat.gz","\n")
    cmd_11 <- paste0("#plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                     path,"srcLog/bbbb-",c("02"),".mat.gz"," --samplesLabel ",
                     paste0(c("Fbl","pr"), collapse = " ")," -out ",
                     path,"thePlot/",j,"/",
                     gsub(".bed","",basename(c(bed6))),".pdf","\n")
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    for (i in cmd_11) {cat(i, append = T, file = shell)}
    tags <- c("pr","Fbl","G9a","Trim28","Setdb1")
    cmd_12 <- paste0("#multiBigwigSummary bins -p 20 -b ",
                     paste0(path,"bw/",tags,"/","IP_Input-",j,".bw", collapse = " ")," ",
                     "-l ",paste0(tags,collapse = " ")," ",
                     "-o ",path,"srcLog/bbbb-",c("03"),".mat.gz","\n")
    cmd_13 <- paste0("#plotCorrelation -in ",path,
                     "srcLog/bbbb-",c("03"),".mat.gz ",
                     "-c spearman -p heatmap --plotFileFormat pdf --removeOutliers ",
                     "--colorMap bwr -o ", path,"thePlot/",j,
                     "/Corre_R1_genom-spearman.pdf","\n")
    cmd_14 <- paste0("#plotCorrelation -in ",path,
                     "srcLog/bbbb-",c("03"),".mat.gz ",
                     "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                     "--colorMap bwr -o ", path,"thePlot/",j,
                     "/Corre_R1_genom-pearson.pdf","\n")
    for (i in cmd_12) {cat(i, append = T, file = shell)}
    for (i in cmd_13) {cat(i, append = T, file = shell)}
    for (i in cmd_14) {cat(i, append = T, file = shell)}
    for (n in c("H3K9me3_WT","Trim28","Setdb1")) {
      cmd_15 <- paste0("#computeMatrix scale-regions -p 20 -m 1000 -a 2500 -b ",
                       "2500 -R ", c(bed4), " -S ",path,"bw/",n,"/",
                       "IP_Input-",j,".bw ",
                       "--skipZeros -o ",path,
                       "srcLog/cccc-",c("01"),".mat.gz","\n")
      cmd_16 <- paste0("#plotProfile -m ",path,
                       "srcLog/cccc-",c("01"),".mat.gz ",
                       "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                       "-o ",path,"thePlot/",j,"/",n,
                       "-",gsub(".bed","",basename(c(bed4))),".pdf","\n")
      for (i in cmd_15) {cat(i, append = T, file = shell)}
      for (i in cmd_16) {cat(i, append = T, file = shell)}
    }
    tags <- c("H3K9me3_WT","H3K9me3_RN","H3K9me3_KO")
    cmd_17 <- paste0("#computeMatrix reference-point -p 20 --referencePoint center ",
                     "-a 3000 -b 3000 --missingDataAsZero -R ",c(bed1,bed2,bed4)," -S ",
                     paste0(path,"bw/",tags,"/","IP_Input-",j,".bw",collapse = " "),
                     " -o ",path,"srcLog/dddd-",c("01","02","03"),".mat.gz","\n")
    cmd_18 <- paste0("#plotProfile --perGroup -m ",path,
                     "srcLog/dddd-",c("01","02","03"),".mat.gz ",
                     "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                     "--samplesLabel ",paste0(tags,collapse = " ")," ",
                     "-o ",path,"thePlot/",j,
                     "/WT_RN_KOH3K9me3-",gsub(".bed","",basename(c(bed1,bed2,bed4))),"-center.pdf","\n")
    for (i in cmd_17) {cat(i, append = T, file = shell)}
    for (i in cmd_18) {cat(i, append = T, file = shell)}
    cmd_19 <- paste0("computeMatrix scale-regions -p 20 -m 1000 -a 2500 ",
                     "-b 2500 --missingDataAsZero -R ",c(bed4)," -S ",
                     paste0(path,"bw/",tags,"/","IP_Input-",j,".bw",collapse = " "),
                     " -o ",path,"srcLog/dddd-",c("04"),".mat.gz","\n")
    cmd_20 <- paste0("plotProfile --perGroup -m ",path,
                     "srcLog/dddd-",c("04"),".mat.gz ",
                     "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                     "--samplesLabel ",paste0(tags,collapse = " ")," ",
                     "-o ",path,"thePlot/",j,
                     "/WT_RN_KOH3K9me3-",gsub(".bed","",basename(c(bed4))),"-scale.pdf","\n")
    for (i in cmd_19) {cat(i, append = T, file = shell)}
    for (i in cmd_20) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"srcLog/run_a.log"," 2>&1 &"))
  ###
  ###
  ###
  ###
  ###
  ###
  shell <- paste0(path,"srcLog/run_f.sh")
  cat("#!/bin/bash\n", file = shell)
  for (e in para_u[c(5,6,7)]) {
    cmd_21 <- paste0("#computeMatrix reference-point -p 32 --referencePoint center ",
                     "-a 3000 -b 3000 --missingDataAsZero -R ",c(bed3)," -S ",
                     path,"bw/H3K9me3_WT/","IP_Input-",e,".bw ",
                     "-o ",path,"srcLog/ffff-",c("01"),".mat.gz","\n")
    cmd_22 <- paste0("#plotProfile -m ",path,
                     "srcLog/ffff-",c("01"),".mat.gz ",
                     "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                     "-o ",path,"thePlot/",e,"/H3K9me3_WT",
                     "-",gsub(".bed","",basename(c(bed3))),".pdf","\n")
    for (i in cmd_21) {cat(i, append = T, file = shell)}
    for (i in cmd_22) {cat(i, append = T, file = shell)}
  }
  for (f in para_u[c(5,6,7)]) {
    tags <- c("H3K9me3_WT","H3K9me3_RN","H3K9me3_KO")
    cmd_23 <- paste0("#computeMatrix reference-point -p 32 --referencePoint center ",
                     "-a 3000 -b 3000 --missingDataAsZero -R ",c(bed3)," -S ",
                     paste0(path,"bw/",tags,"/","IP_Input-",f,".bw",collapse = " "),
                     " -o ",path,"srcLog/ffff-",c("01"),".mat.gz","\n")
    cmd_24 <- paste0("#plotProfile --perGroup -m ",path,
                     "srcLog/ffff-",c("01"),".mat.gz ",
                     "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
                     "--samplesLabel ",paste0(tags,collapse = " ")," ",
                     "-o ",path,"thePlot/",f,
                     "/WT_RN_KOH3K9me3-",gsub(".bed","",basename(c(bed3))),".pdf","\n")
    for (i in cmd_23) {cat(i, append = T, file = shell)}
    for (i in cmd_24) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"srcLog/run_f.log"," 2>&1 &"))
}
################################################################################
#在Zscan4的9个旁系同源基因 ± 5kb [测试]
{
  ###path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/"
  ###red1 <- read.table("/Reference/aaSHY/BED/special/allGene-2.bed",header = F)
  ###red2 <- red1[grep("Zscan4",red1$V4),]
  ###write.table(red2, 
  ###            file = paste0(path,"dumpROOM/Zscan4/Zscan4.bed"),
  ###            sep = "\t", quote = F, col.names = F, row.names = F)
  ###bed1 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/dumpROOM/Zscan4/Zscan4.bed"
  ###bw01 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bw/pr/IP_Input-SES.bw"
  ###bw02 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bw/Fbl/IP_Input-SES.bw"
  ###bw03 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bw/H3K9me3_WT/IP_Input-SES.bw"
  ###bw04 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bw/H3K9me3_RN/IP_Input-SES.bw"
  ###bw05 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20250303_bcp/bw/H3K9me3_KO/IP_Input-SES.bw"
  ###tags <- c("pr","Fbl","H3K9me3_WT","H3K9me3_RN","H3K9me3_KO")
  ###cmd_01 <- paste0("computeMatrix scale-regions -p 20 -m 5000 -a 5000 -b ",
  ###                 "5000 -R ", c(bed1), " ",
  ###                 "-S ",paste0(c(bw01,bw02,bw03,bw04,bw05), collapse = " ")," ",
  ###                 "--skipZeros -o ",path,"dumpROOM/Zscan4/the.mat.gz","\n")
  ###cmd_02 <- paste0("plotProfile --perGroup -m ",path,
  ###                 "dumpROOM/Zscan4/the.mat.gz ",
  ###                 "--samplesLabel ",paste0(tags,collapse = " ")," ",
  ###                 "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
  ###                 "-o ",path,"dumpROOM/Zscan4/ChIP_signal_Zscan4_v1.pdf","\n")
  ###cmd_03 <- paste0("plotProfile -m ",path,
  ###                 "dumpROOM/Zscan4/the.mat.gz ",
  ###                 "--samplesLabel ",paste0(tags,collapse = " ")," ",
  ###                 "--plotHeight 10 --plotWidth 11 --legendLocation upper-right ",
  ###                 "-o ",path,"dumpROOM/Zscan4/ChIP_signal_Zscan4_v2.pdf","\n")
}
















