#write by ~~~ at 2024.03.15
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
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
  ff01 <- list.files(paste0(path, "fq"), full.names = T, pattern = "R1")
  ff02 <- list.files(paste0(path, "fq"), full.names = T, pattern = "R2")
  prex <- gsub(basename(ff01), pattern="_R.*", replacement="")
  path <- paste0(path,"single/")
  for (i in c("src","trim","Log","peak","dumpROOM","DESeq2","rRNA/csem","mervl/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  rm(i)
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_01 <- paste0("#zcat ",ff01," ",ff02, " |gzip >",path,"fq/",prex,".fq.gz","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  cmd_02 <- paste0("#trim_galore --phred33 --fastqc --illumina ",path,
                   "fq/", prex, ".fq.gz -o ", path, "trim", "\n")
  cmd_03 <- paste0("#bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inx1," -U ", path, "trim/", prex, "_trimmed.fq.gz",
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_04 <- paste0("#run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"bam/",prex,".bam")," 140 ",
                   paste0(path,"bam/csem/",prex),"\n")
  cmd_05 <- paste0("#bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/csem/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (x in seq(2)) {
    zrex <- prex[c(x,x+2)]
    for (q in zrex[2]) {
      cmd_07 <- paste0("#bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                       paste0(path,"bam/csem/",zrex[2],".sorted.bam -b2 "),
                       paste0(path,"bam/csem/",zrex[1],".sorted.bam -o "),
                       paste0(path,"bw/bcp/",zrex[2],"_Input.bw"),"\n")
      cmd_08 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 5000 ",
                       "--missingDataAsZero -p 20 ", "-R ", c(inx7,inx8), " -S ", path,
                       "bw/bcp/",q,"_Input.bw -o ", path,"deeptools/", 
                       c(basename(inx7),basename(inx8)),q, ".mat.gz","\n")
      cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 500 ",
                       "--missingDataAsZero -p 20 -R ", inx9, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inx9),q,".mat.gz","\n")
      cmd_10 <- paste0("#computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                       "--missingDataAsZero -p 20 -R ", inxA, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inxA),q,".mat.gz","\n")
      cmd_11 <- paste0("#computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 30000 ",
                       "--missingDataAsZero -p 20 -R ", inxB, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inxB),q,".mat.gz","\n")
      cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                       "deeptools/",c(basename(inx7),basename(inx8),basename(inx9),
                                      basename(inxA),basename(inxB)),q,
                       ".mat.gz -out ",path,"PLOT/",
                       c(basename(inx7),basename(inx8),
                         basename(inx9),basename(inxA),basename(inxB)),q,"_v2.pdf","\n")
      for (i in cmd_07) {cat(i, append = T, file = shell)}
      for (i in cmd_08) {cat(i, append = T, file = shell)}
      for (i in cmd_09) {cat(i, append = T, file = shell)}
      for (i in cmd_10) {cat(i, append = T, file = shell)}
      for (i in cmd_11) {cat(i, append = T, file = shell)}
      for (i in cmd_12) {cat(i, append = T, file = shell)}
    }
    rm(q)
    cmd_13 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                     path,"bam/csem/",zrex[2],".sorted.bam -c ",
                     path,"bam/csem/",zrex[1],".sorted.bam -n ",
                     zrex[2],"-N --outdir ",
                     path,"peak/", "\n")
    cmd_14 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                     "--keep-dup 1 -t ",path,"bam/csem/", zrex[2],
                     ".sorted.bam -c ", path,"bam/csem/", zrex[1],
                     ".sorted.bam -n ", zrex[2],
                     "-B --outdir ", path,"peak/", "\n")
    cmd_15 <- paste0("TEtranscripts -t ",
                     paste0(path,"bam/csem/",zrex[2],".sorted.bam", collapse = " ")," -c ",
                     paste0(path,"bam/csem/",zrex[1],".sorted.bam", collapse = " "),
                     " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                     "--project pr", zrex[2], " --outdir ",paste0(path,"count/"),"\n")
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
################################################################################
#04、peak Analysis (pvalue 0.05)
{
  gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
  gtf1 <- as.data.frame(gtf1)
  colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
  gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
  txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
  twoC <- read.table(inxB)
  gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
  txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
  #
  for (kk in c("allg","twoc")) {
    if (kk == "allg") {
      inxM <- txb1
    }
    else {
      inxM <- txb2
    }
    for (i in (prex[3:4])) {
      pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "single/peak/",i,"-N_peaks.narrowPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = inxM))
      pdf(file = paste0(path, "single/PLOT/",i,"-",kk,"-N_chipseekerAnno.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk01), paste0(path, "single/peak/",i,"-",kk,"-N_chipseekerAnno.csv"), quote = F, row.names = F)
    } 
  }
}
################################################################################
#05、gene费舍尔检验 [适用于ChIP-seq]
{
  bath = paste0(path, "single/")
  prex <- unique(sapply(str_split(list.files(path = paste0(bath, "fq")), pattern = "-R|.fq"), "[", 1))
  for (x in seq(2)) {
    zrex <- prex[c(x,x+2)]; rm(x)
    twoc <- read.table(file = inxB, sep = "\t")
    #Observed
    red1 <- read.table(file = paste0(bath,"count/pr_",zrex[2],".cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
    colnames(red1) <- zrex[c(2,1)]
    red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
    #Expected
    mmu1 <- read.table(file = paste0(inx3,".fai"))
    refL <- sum(mmu1$V2)
    mmu2 <- rtracklayer::import(con = inx4, format = "gtf")
    mmu3 <- as.data.frame(mmu2[which(mmu2$type=="exon")])
    for (m in seq(1)) {
      print(m)
      red3 <- red2[,c(m,m+1)]
      for (i in seq_along(rownames(red3))){
        print(i)
        red3[i,3] <- round(sum(mmu3[which(mmu3[,12] == rownames(red3)[i]),4])/refL*sum(red1[,1]),4)
        red3[i,4] <- round(sum(mmu3[which(mmu3[,12] == rownames(red3)[i]),4])/refL*sum(red1[,2]),4)
        red3[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red3[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
        red3[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red3[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
        red3[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red3[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
      }
      colnames(red3) <- c("IP","Input","expcIP","expcInput","pvalue","greate","lesser")
      red3$pMark[red3$pvalue<0.05] <- "arresting"
      red3$gMark[red3$greate<0.05] <- "arresting"
      red3$lMark[red3$lesser<0.05] <- "arresting"
      red3$geneName <- rownames(red3)
      red3$marker <- NA
      red3$marker[which(red3$geneName %in% twoc$V4)] <- "YES"
      red3$marker[which(!(red3$geneName %in% twoc$V4))] <- "NO"
      red3 <- red3[order(red3$greate, decreasing = F),]
      openxlsx::write.xlsx(x = red3, asTable = T, file = paste0(bath,"fisher/",zrex[2],"FisherGene.xlsx"))
    }
  }; rm(m)
}
################################################################################
#06、tete费舍尔检验 [适用于ChIP-seq]
{
  bath = paste0(path, "single/")
  prex <- unique(sapply(str_split(list.files(path = paste0(bath, "fq")), pattern = "-R|.fq"), "[", 1))
  for (x in seq(2)) {
    zrex <- prex[c(x,x+2)]; rm(x)
    #Observed
    red1 <- read.table(file = paste0(bath,"count/pr_",zrex[2],".cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
    colnames(red1) <- zrex[c(2,1)]
    red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
    #Expected
    mmu1 <- read.table(file = paste0(inx3,".fai"))
    refL <- sum(mmu1$V2)
    mmu2 <- read.table(file = inx5, sep = "\t", stringsAsFactors = F, header = F)
    mmu3 <- mmu2
    mmu3[,9] <- gsub("; transcript_id.*family_id |; class_id ",":",gsub("gene_id |;$","",mmu3[,9]))
    mmu3[,8] <- mmu3[,5]-mmu3[,4]+1
    for (m in seq(1)) {
      print(m)
      red3 <- red2[,c(m,m+1)]
      for (i in seq_along(rownames(red3))){
        print(i)
        red3[i,3] <- round(sum(mmu3[which(mmu3[,9] == rownames(red3)[i]),8])/refL*sum(red3[,1]),4)
        red3[i,4] <- round(sum(mmu3[which(mmu3[,9] == rownames(red3)[i]),8])/refL*sum(red3[,2]),4)
        red3[i,5] <- suppressWarnings(fisher.test(matrix(as.numeric(red3[i,1:4]),ncol = 2,byrow = F),alternative = "two.sided")$p.value)
        red3[i,6] <- suppressWarnings(fisher.test(matrix(as.numeric(red3[i,1:4]),ncol = 2,byrow = F),alternative = "greater")$p.value)
        red3[i,7] <- suppressWarnings(fisher.test(matrix(as.numeric(red3[i,1:4]),ncol = 2,byrow = F),alternative = "less")$p.value)
      }
      colnames(red3) <- c("IP","Input","expcIP","expcInput","pvalue","greate","lesser")
      #
      red3$pMark[red3$pvalue<0.05] <- "arresting"
      red3$gMark[red3$greate<0.05] <- "arresting"
      red3$lMark[red3$lesser<0.05] <- "arresting"
      red3$repeatName <- rownames(red3)
      openxlsx::write.xlsx(x = red3, asTable = T, file = paste0(bath,"fisher/",zrex[2],"FisherTE.xlsx"))
    }; rm(m)
  }
}
################################################################################
#XD、以H3K9me3 WT-1为主体, 画pr、Fbl的富集
{
  inx0 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/peak/IP-1-N_peaks.narrowPeak"
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  bath = paste0(path, "single/")
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  shell <- paste0(bath,"dumpROOM/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("computeMatrix reference-point --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -p 20 -R ", inx0,
                   " -S",
                   " /ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-1_Input.bw ",
                   paste0(bath, "bw/bcp/", c("IP-1","IP-2"),"_Input.bw", collapse = " ")," ",
                   paste0(fath, "bw/bcp/", c("IP-1","IP-2"),"_Input.bw", collapse = " "),
                   " -o ",bath,"dumpROOM/coLocation_v1.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                   bath, "dumpROOM/coLocation_v1.mat.gz --samplesLabel ",
                   paste0(c("WT-1","pr-1","pr-2","Fbl-1","Fbl-2"), collapse = " "),
                   " -out ",bath,"dumpROOM/coLocation_v1.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",bath,"dumpROOM/run_a.log"," 2>&1 &"))
}
################################################################################
#XD、画Ssrp1-1, Npm1-1, pr-1, Fbl-1, Setdb1-1, H3K9me3-WT-1, Kap1-2的相关性热图
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  bath = paste0(path, "single/")
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  shell <- paste0(bath,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  bw01 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-1_Input.bw",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-3_Input.bw")[1]
  bw02 <- list.files(paste0(bath,"bw/bcp"), full.names = T)[1]
  bw03 <- list.files(paste0(fath,"bw/bcp"), full.names = T)[1]
  bw04 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/KAP1/bw/bcp/KAP1_vs.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/bw/bcp/Trim28_vs.bw")[c(1,4)]
  bw05 <- c("/ChIP_seq_2/Ssrp1_ChIP_bp/bamCompare/D419_Ssrp1_mm10.bw",
            "/ChIP_seq_2/aaSHY/Npm1/ChIP-seq/bw/bcp/Npm1IP-Input.bw")
  tags <- c("H3K9me3-WT-1","H3K9me3-WT-2","pr-1","pr-2","Fbl-1","Fbl-2",
            "Setdb1-1","Setdb1-2","Kap1-1","Kap1-2","Ssrp1-1","Npm1-1")
  tags <- c("H3K9me3-WT-1","pr-1","Fbl-1","Setdb1-1","Kap1-2","Ssrp1-1","Npm1-1")
  cmd_01 <- paste0("multiBigwigSummary bins -p 20 -b ", 
                   paste0(c(bw01,bw02,bw03,bw04,bw05), collapse = " "),
                   " -l ",paste0(tags,collapse = " "),
                   " -o ",paste0(bath,"dumpROOM/pearsonCorre.npz","\n"))
  cmd_02 <- paste0("plotCorrelation -in ",bath,
                   "dumpROOM/pearsonCorre.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", bath,"dumpROOM/pearsonCorre_v1-xd.pdf ",
                   "--outFileCorMatrix ",bath,"dumpROOM/pearsonCorre_v1-xd.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(bath,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#XD、画Ssrp1-1, Npm1-1, pr-1, Fbl-1, Setdb1-1, H3K9me3-WT-1, Kap1-2的相关性热图 [改进]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  bath = paste0(path, "single/")
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  shell <- paste0(bath,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  bw01 <- c("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-1_Input.bw",
            "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-3_Input.bw")
  bw02 <- list.files(paste0(bath,"bw/bcp"), full.names = T)
  bw03 <- list.files(paste0(fath,"bw/bcp"), full.names = T)
  bw04 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/KAP1/bw/bcp/KAP1_vs.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/bw/bcp/Trim28_vs.bw")
  bw05 <- c("/ChIP_seq_2/Ssrp1_ChIP_bp/bamCompare/D419_Ssrp1_mm10.bw",
            "/ChIP_seq_2/aaSHY/Npm1/ChIP-seq/bw/bcp/Npm1IP-Input.bw")
  bw06 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_1_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_2_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP051217_Control_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP069149_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw")
  tags <- c("H3K9me3-WT-1","H3K9me3-WT-2","pr-1","pr-2","FBL-1","FBL-2",
            "SETDB1-1","SETDB1-2","KAP1-1","KAP1-2","SSRP1-1","NPM1-1",
            "SRP001533-1","SRP001533-2","SRP051217-1","SRP064188-1","SRP069149-1","SRP094580-1")
  cmd_01 <- paste0("multiBigwigSummary bins -p 20 -b ",
                   paste0(c(bw01,bw02,bw03,bw04,bw05,bw06), collapse = " "),
                   " -l ",paste0(tags,collapse = " "),
                   " -o ",paste0(bath,"dumpROOM/pearsonCorre.npz","\n"))
  cmd_02 <- paste0("plotCorrelation -in ",bath,
                   "dumpROOM/pearsonCorre.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", bath,"dumpROOM/pearsonCorre_v2-xd.pdf ",
                   "--outFileCorMatrix ",bath,"dumpROOM/pearsonCorre_v2-xd.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(bath,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#XD、画Ssrp1-1, Npm1-1, pr-1, Fbl-1, Setdb1-1, H3K9me3-WT-1, Kap1-2的相关性热图 [改进]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  bath = paste0(path, "single/")
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  shell <- paste0(bath,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  bw01 <- list.files(paste0(bath,"bw/bcp"), full.names = T)[1]
  bw02 <- list.files(paste0(fath,"bw/bcp"), full.names = T)[1]
  bw03 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP021942_G9a_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP069149_Kap1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP007651_Kap1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP090985_Kap1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP098642_Kap1_pseudocount_1_log2.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/KAP1/bw/bcp/KAP1_vs.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/bw/bcp/Trim28_vs.bw")
  bw04 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_TF_Hdac1_pseudocount_1_log2.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP136128_Hdac1_pseudocount_1_log2.bw",
            "/disk5/zx/ChIP_seq/Hdac4_5_HA_merge_rep1/bamcompare/Hdac4_ChIP_VS_Hdac4_input_SES.bw",
            "/disk5/zx/ChIP_seq/HDac45_AND_5HA_JAICE_rep2_merge/bamcompare/HDac4_merge_VS_HDac4_input_merge_SES.bw",
            "/disk5/zx/ChIP_seq/Hdac4_5_HA_merge_rep1/bamcompare/Hdac5_ChIP_VS_Hdac5_input_SES.bw",
            "/disk5/zx/ChIP_seq/HDac45_AND_5HA_JAICE_rep2_merge/bamcompare/HDac5_merge_VS_HDac5_input_merge_SES.bw")
  bw05 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw")
  bw06 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/Glp/bw/bcp/IP_Input.bw")
  tags <- c("pr-1","FBL-1","G9a",
            "KAP1-1","KAP1-2","KAP1-3","KAP1-4","KAP1-5","KAP1-6",
            "HDAC1-1","HDAC1-2","HDAC4-1","HDAC4-2","HDAC5-1","HDAC5-2",
            "SETDB1-1","SRP094580-1","Glp")
  cmd_01 <- paste0("multiBigwigSummary bins -p 20 -b ",
                   paste0(c(bw01,bw02,bw03,bw04,bw05,bw06), collapse = " "),
                   " -l ",paste0(tags,collapse = " "),
                   " -o ",paste0(bath,"dumpROOM/pearsonCorre.npz","\n"))
  cmd_02 <- paste0("plotCorrelation -in ",bath,
                   "dumpROOM/pearsonCorre.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", bath,"dumpROOM/pearsonCorre_v4-xd.pdf ",
                   "--outFileCorMatrix ",bath,"dumpROOM/pearsonCorre_v4-xd.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(bath,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#XD、画Ssrp1-1, Npm1-1, pr-1, Fbl-1, Setdb1-1, H3K9me3-WT-1, Kap1-2的相关性热图 [数据、参数需要在文章统一]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  bath = paste0(path, "single/")
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  shell <- paste0(bath,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  bw01 <- list.files(paste0(bath,"bw/bcp"), full.names = T)[1]
  bw02 <- list.files(paste0(fath,"bw/bcp"), full.names = T)[1]
  bw03 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP021942_G9a_pseudocount_1_log2.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/CTCF/bw/bcp/IP_Input.bw",
            "/disk5/Mouse_ESC_ChIP/bamCompare/SRP098642_Kap1_pseudocount_1_log2.bw")
  bw04 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw")
  tags <- c("pr-1","FBL-1","G9a","CTCF","KAP1-4","SETDB1-1")
  cmd_01 <- paste0("multiBigwigSummary bins -p 20 -b ",
                   paste0(c(bw01,bw02,bw03,bw04), collapse = " "),
                   " -l ",paste0(tags,collapse = " "),
                   " -o ",paste0(bath,"dumpROOM/pearsonCorre.npz","\n"))
  cmd_02 <- paste0("plotCorrelation -in ",bath,
                   "dumpROOM/pearsonCorre.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", bath,"dumpROOM/pearsonCorre_v6-xd.pdf ",
                   "--outFileCorMatrix ",bath,"dumpROOM/pearsonCorre_v6-xd.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(bath,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#XD、画pr-1, Fbl-1, Setdb1-1, H3K9me3-WT-1, Kap1-2的相关性热图 [20250210投稿前确定下来]-------------
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  bath = paste0(path, "single/")
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  shell <- paste0(bath,"src/run_c.sh")
  cat("#!/bin/bash\n", file = shell)
  ###not SES
  ###bw01 <- list.files(paste0(bath,"bw/bcp"), full.names = T)[1]
  ###bw02 <- list.files(paste0(fath,"bw/bcp"), full.names = T)[1]
  ###bw03 <- c("/ChIP_seq_2/aaSHY/pr/publicData/G9a/bw/bcp/IP_Input-notSES.bw",
  ###          "/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/n1/bw/bcp/IP_Input-notSES.bw")
  ###bw04 <- c("/ChIP_seq_2/aaSHY/pr/publicData/SETDB1/n2/bw/bcp/IP_Input-notSES.bw")
  ###SES
  bw01 <- list.files(paste0(bath,"bw/bcp"), full.names = T)[2]
  bw02 <- list.files(paste0(fath,"bw/bcp"), full.names = T)[2]
  bw03 <- c("/ChIP_seq_2/aaSHY/pr/publicData/G9a/bw/bcp/IP_Input.bw",
            "/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/n1/bw/bcp/IP_Input.bw")
  bw04 <- c("/ChIP_seq_2/aaSHY/pr/publicData/SETDB1/n2/bw/bcp/IP_Input.bw")
  tags <- c("pr-1","FBL-1","GSE46536_G9a","GSE166041_TRIM28","GSE18371_SETDB1")
  cmd_01 <- paste0("multiBigwigSummary bins -p 20 -b ",
                   paste0(c(bw01,bw02,bw03,bw04), collapse = " "),
                   " -l ",paste0(tags,collapse = " "),
                   " -o ",paste0(bath,"dumpROOM/spearmanCorre_genome.npz","\n"))
  cmd_02 <- paste0("#multiBigwigSummary BED-file -p 20 -b ",
                   paste0(c(bw01,bw02,bw03,bw04), collapse = " "),
                   " -l ",paste0(tags,collapse = " "),
                   " --BED ",inxA,
                   " -o ",paste0(bath,"dumpROOM/spearmanCorre_transcriptome.npz","\n"))
  cmd_03 <- paste0("plotCorrelation -in ",bath,
                   "dumpROOM/spearmanCorre_genome.npz ",
                   "-c spearman -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", bath,"dumpROOM/spearmanCorre_v9_genom.pdf ",
                   "--outFileCorMatrix ",bath,"dumpROOM/spearmanCorre_v9_genom.index.txt","\n")
  cmd_04 <- paste0("#plotCorrelation -in ",bath,
                   "dumpROOM/spearmanCorre_transcriptome.npz ",
                   "-c pearson -p heatmap --plotFileFormat pdf --removeOutliers ",
                   "--colorMap bwr -o ", bath,"dumpROOM/spearmanCorre_v8_trans.pdf ",
                   "--outFileCorMatrix ",bath,"dumpROOM/spearmanCorre_v8_trans.index.txt","\n")
  for (m in cmd_01){cat(m, file = shell, append = T)}
  for (m in cmd_02){cat(m, file = shell, append = T)}
  for (m in cmd_03){cat(m, file = shell, append = T)}
  for (m in cmd_04){cat(m, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(bath,"Log/"),basename(shell),".log", " 2>&1 &"))
}
################################################################################
#09、peak Anno: 看Fbl、pr promoter peak落在哪些基因上
{
  #peak Anno
  {
    gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
    gtf1 <- as.data.frame(gtf1)
    colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
    gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
    txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
    twoC <- read.table(inxB)
    gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
    txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
    #
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    for (kk in c("allg","twoc")) {
      if (kk == "allg") {
        inxM <- txb1
      }
      else {
        inxM <- txb2
      }
      for (i in (prex[3:4])) {
        pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "single/peak/",i,"-N_peaks.narrowPeak"), 
                                                          tssRegion = c(-2000, 2000), TxDb = inxM))
        pdf(file = paste0(path, "single/PLOT/",i,"-",kk,"-N_chipseekerAnno.pdf"), width = 8, height = 8)
        ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
        dev.off()
        pk01 <- as.data.frame(pk01)
        pk01$annotation <- gsub(x = pk01$annotation, pattern = ", ", replacement = "_")
        write.csv(pk01, paste0(path, "single/peak/",i,"-",kk,"-N_chipseekerAnno.csv"), quote = F, row.names = F)
      } 
    }
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  }
  #find 2-cell genes, which Peak at its Promoter [太少了, 不行]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-twoc-N_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter"),]
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-2-twoc-N_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter"),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp3 <- read.csv(file = paste0(path, "single/peak/IP-1-twoc-N_chipseekerAnno.csv"), header = T)
    tmp3 <- tmp3[grep(tmp3$annotation, pattern="Promoter"),]
    tmp4 <- read.csv(file = paste0(path, "single/peak/IP-2-twoc-N_chipseekerAnno.csv"), header = T)
    tmp4 <- tmp4[grep(tmp4$annotation, pattern="Promoter"),]
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$geneId)
    temp[["d"]] <- unique(tmp4$geneId)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "twoc genes, which Peak at its Promoter",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "pr-2","FBL-1", "FBL-2"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9","#B8E5FA"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    intersect(intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]]),temp[["d"]])
    #intersect(temp[["a"]],temp[["c"]])
    #intersect(temp[["b"]],temp[["c"]])
  }
  #find the genes, which Peak at its Promoter [交集是101基因, 继续分析]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter"),]
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-2-allg-N_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter"),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp3 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp3 <- tmp3[grep(tmp3$annotation, pattern="Promoter"),]
    tmp4 <- read.csv(file = paste0(path, "single/peak/IP-2-allg-N_chipseekerAnno.csv"), header = T)
    tmp4 <- tmp4[grep(tmp4$annotation, pattern="Promoter"),]
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$geneId)
    temp[["d"]] <- unique(tmp4$geneId)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, which Peak at its Promoter",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "pr-2","FBL-1", "FBL-2"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9","#B8E5FA"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    ne01 <- intersect(intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]]),temp[["d"]])
  }
  #交集的101基因, 在早期胚胎发育阶段的 expression pheatmap
  {
    #用2014 Deng的数据 [不行]
    cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
    cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
    prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
    prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
    colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
    cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
    cn02$geneName <- rownames(cn02)
    cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = T),]#invert = F, 那么输出TE的CPM
    df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
    df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
    df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
    df03 <- tidyr::spread(data = df02[,-4], key = "stage", value = "cpmExp")
    df03 <- df03[which(df03$geneName %in% ne01),]
    rownames(df03) <- df03$geneName; df03 <- df03[,-1]
    df03 <- df03[,c("oocyte","zygote","twoEarly","twoMid","twoLate","four","eight","blastoEarly","blastoMid","blastoLate")]
    #
    df03 <- df03[rowSums(df03)>0,]
    pheatmap(mat = df03,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             labels_col = colnames(df03), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
  }
  {
    #用2020 pacbio corresponding illumina的数据 [可以]
    stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/count/star/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    colnames(stag) <- sapply(str_split(colnames(stag), "\\."), "[", 8)
    stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))
    stag <- stag[grep(x = rownames(stag), pattern = ":", invert = T), -7]
    stag <- stag[,c("oocyte","zygote","two","four","eight","blasto")]
    df03 <- stag[ne01,]#;df03 <- df03[rowSums(df03)==0,]; rownames(df03)
    df03 <- df03[rowSums(df03)>0,]
    pdf("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/peak/prFbl.promoterGene.v1.pdf",
        width = 10, height = 15)
    pheatmap(mat = df03,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             labels_col = colnames(df03), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    dev.off()
  }
  #这些基因有多少是与MERVL MT2融合的?
  {
    ne02 <- rownames(df03)
    gf01 <- import(con = inx4)
    gf01 <- gf01[which(gf01$type=="gene")]
    tf01 <- import(con = inx5)
    tf01 <- tf01[tf01$gene_id %in% c("MERVL-int","MT2_Mm")]
    #组装的转录本中有哪些是MERVL嵌合的 [MERVL与转录本的外显子嵌合, 则算嵌合]
    {
      gf02 <- import(con = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf")
      gf02 <- gf02[which(gf02$type == "exon")]
      ov01 <- suppressWarnings(findOverlaps(query = gf02, subject = tf01, ignore.strand=T))
      mf02 <- gf02[unique(queryHits(ov01))]
      #对嵌合的转录本做注释
      gf02 <- import(con = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf")
      gf02 <- gf02[which(gf02$transcript_id %in% unique(mf02$transcript_id))]
      ov02 <- suppressWarnings(findOverlaps(query = gf02, subject = gf01, ignore.strand=T))
      theG <- gf01[unique(subjectHits(ov02))]
      intersect(ne02,unique(theG$gene_name))
    }
    #组装的转录本中有哪些是MERVL嵌合的 [MERVL与转录本嵌合, 则算嵌合]
    {
      gf02 <- import(con = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf")
      gf02 <- gf02[which(gf02$type == "transcript")]
      ov01 <- suppressWarnings(findOverlaps(query = gf02, subject = tf01, ignore.strand=T))
      mf02 <- gf02[unique(queryHits(ov01))]
      #对嵌合的转录本做注释
      gf02 <- mf02
      ov02 <- suppressWarnings(findOverlaps(query = gf02, subject = gf01, ignore.strand=T))
      theG <- gf01[unique(subjectHits(ov02))]
      intersect(ne02,unique(theG$gene_name))
    }
  }
}
#10、调试/测试Macs2 [因为Nelfa不在call出的peak里]
{
  #peak Anno
  #测试Macs2参数
  #跑macs
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
         path,"single/bam/csem/IP-1.sorted.bam -c ",
         path,"single/bam/csem/Input-1.sorted.bam -n ",
         "testDup-N --outdir ",
         path,"single/peak/", "\n")
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "single/peak/testDup-N_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pk01 <- as.data.frame(pk01)
  pk01$annotation <- gsub(x = pk01$annotation, pattern = ", ", replacement = "_")
  pk02 <- makeGRangesFromDataFrame(df = pk01, keep.extra.columns = F, ignore.strand = T)
  #
  tstG <- gtf1[which(gtf1$type == "gene" & gtf1$gene_id == "Nelfa")]
  #
  ov01 <- GenomicRanges::findOverlaps(query = pk02, subject = tstG)
  pk01[76760,]
  write.csv(pk01, paste0(path, "single/peak/testDup-N_chipseekerAnno.csv"), quote = F, row.names = F)
}
################################################################################
#11、peak Anno: 看Fbl、pr promoter peak落在哪些基因上 [改进, 具体如下]
{
  #peak Anno
  {
    gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
    gtf1 <- as.data.frame(gtf1)
    colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
    gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
    txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
    twoC <- read.table(inxB)
    gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
    txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
    #
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    for (kk in c("allg","twoc")) {
      if (kk == "allg") {
        inxM <- txb1
      }
      else {
        inxM <- txb2
      }
      for (i in (prex[3:4])) {
        pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "single/peak/",i,"-B_peaks.broadPeak"), 
                                                          tssRegion = c(-2000, 2000), TxDb = inxM))
        pdf(file = paste0(path, "single/PLOT/",i,"-",kk,"-B_chipseekerAnno.pdf"), width = 8, height = 8)
        ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
        dev.off()
        pk01 <- as.data.frame(pk01)
        pk01$annotation <- gsub(x = pk01$annotation, pattern = ", ", replacement = "_")
        write.csv(pk01, paste0(path, "single/peak/",i,"-",kk,"-B_chipseekerAnno.csv"), quote = F, row.names = F)
      } 
    }
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  }
  #find the genes, [方案01: Narrow peak, which Peak at its Promoter] [Venn]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter"),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter"),]
    tmp3 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/grid.mervlNum.gene.xlsx")
    tmp3 <- tmp3[which(tmp3$rank >3 & tmp3$biotype == "protein_coding"),]
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$name)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, which Peak at its Promoter",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1", "MERVL-tags rank >3"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    ne01 <- intersect(intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]]),temp[["d"]])
  }
  #find the genes, [方案02: Narrow peak, which Peak not the Distal、Downstream] [Venn]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    tmp3 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/grid.mervlNum.gene.xlsx")
    tmp3 <- tmp3[grep(tmp3$name, pattern="mt-", invert = T), ]
    tmp3 <- tmp3[which(tmp3$rank >7 & tmp3$biotype == "protein_coding"),]
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$name)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, which Peak not the Distal、Downstream",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1", "MERVL-tags rank >7"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    
    ne01 <- intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]])
  }
  #find the genes, [方案03: Narrow peak, which Peak at its Promoter] [GSEA]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter"),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter"),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 1000")
  }
  #find the genes, [方案04: Narrow peak, which Peak not the Distal、Downstream] [GSEA]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-N_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","rank")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 4000")
  }
  #find the genes, [方案05: Broad peak, which Peak not the Distal、Downstream] [Venn]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    tmp3 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/grid.mervlNum.gene.xlsx")
    tmp3 <- tmp3[grep(tmp3$name, pattern="mt-", invert = T), ]
    tmp3 <- tmp3[which(tmp3$rank >5 & tmp3$biotype == "protein_coding"),]
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$name)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, which Peak not the Distal、Downstream",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1", "MERVL-tags rank >5"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    
    ne01 <- intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]])
  }
  #find the genes, [方案06: Broad peak, which Peak at Promoter] [Venn]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="promoter", ignore.case = T, invert = F),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="promoter", ignore.case = T, invert = F),]
    tmp3 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/grid.mervlNum.gene.xlsx")
    tmp3 <- tmp3[grep(tmp3$name, pattern="mt-", invert = T), ]
    tmp3 <- tmp3[which(tmp3$rank >3 & tmp3$biotype == "protein_coding"),]
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$name)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, which Peak at the promoter",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1", "MERVL-tags rank >3"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    
    ne01 <- intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]])
  }
  #find the genes, [方案07: Broad peak, which Peak at its Promoter] [GSEA]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter"),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter"),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 2000")
  }
  #find the genes, [方案08: Broad peak, which Peak not the Distal、Downstream] [GSEA]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 4000")
  }
  #find the genes, [方案09: Broad peak, which Peak at its Promoter] [GSEA] [过滤非编码]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter"),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter"),]
    ne01 <- unique(as.character(intersect(tmp1$geneId, tmp2$geneId)))
    gf01 <- import(con = inx4); gf01 <- gf01[which(gf01$type=="gene")]
    ne01 <- intersect(ne01, gf01$gene_name[which(gf01$gene_biotype=="protein_coding")])
    #
    wr01 <- data.frame(theGenes = ne01)
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 2000")
  }
  #find the genes, [方案10: Broad peak, which Peak not the Distal、Downstream] [GSEA] [过滤非编码]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Distal|Downstream", ignore.case = T, invert = T),]
    ne01 <- unique(as.character(intersect(tmp1$geneId, tmp2$geneId)))
    gf01 <- import(con = inx4); gf01 <- gf01[which(gf01$type=="gene")]
    ne01 <- intersect(ne01, gf01$gene_name[which(gf01$gene_biotype=="protein_coding")])
    #
    wr01 <- data.frame(theGenes = ne01)
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01,
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 3000")
  }
  #find the genes, [方案11: Broad peak, all peaks] [Venn]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp3 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/grid.mervlNum.gene.xlsx")
    tmp3 <- tmp3[grep(tmp3$name, pattern="mt-", invert = T), ]
    tmp3 <- tmp3[which(tmp3$rank >3 & tmp3$biotype == "protein_coding"),]
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    temp[["c"]] <- unique(tmp3$name)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, use all peaks",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1", "MERVL-tags rank >3"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    
    ne01 <- intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]])
  }
  #find the genes, [方案12: Broad peak, all peaks] [GSEA]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    wr01 <- data.frame(theGenes = unique(intersect(tmp1$geneId, tmp2$geneId)))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 6000")
  }
  #find the genes, [方案13: Broad peak, all peaks] [GSEA] [过滤非编码]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    ne01 <- unique(as.character(intersect(tmp1$geneId, tmp2$geneId)))
    gf01 <- import(con = inx4); gf01 <- gf01[which(gf01$type=="gene")]
    ne01 <- intersect(ne01, gf01$gene_name[which(gf01$gene_biotype=="protein_coding")])
    #
    wr01 <- data.frame(theGenes = ne01)
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01,
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","rank")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 4000")
  }
  #find the genes, [方案14: Broad peak, all peaks] [Venn] [确定下来的]--------------------------------------
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneId)
    temp[["b"]] <- unique(tmp2$geneId)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, use all broadPeaks",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1"),
                                      fill=c("#FD763F","#23BAC5"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    pdf(file = paste0(path, "single/PLOT/broadPeaks.Venn.v1.pdf"))
    grid.draw(p_01)
    dev.off()
    ne01 <- intersect(temp[["a"]],temp[["b"]])
  }
  #find the genes, [方案15: Broad peak, which Peak at its Promoter] [GSEA] [python]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter", ignore.case = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter", ignore.case = T),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA --max-size 2000")
  }
  #find the genes, [方案16: Broad peak, which Peak at its Promoter] [GSEA] [R] [确定下来的-又丢弃]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter", ignore.case = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter", ignore.case = T),]
    myg1 <- intersect(tmp1$geneId, tmp2$geneId)
    myg2 <- suppressWarnings(clusterProfiler::bitr(wr01, fromType = "SYMBOL", 
                                                   toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
    myg3 <- data.frame(ont = rep("myGene", length(unique(myg2$ENTREZID))),
                       gen = unique(myg2$ENTREZID))
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    lis1 <- suppressWarnings(clusterProfiler::bitr(geneID = myta$name, fromType = "SYMBOL",
                                                   toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = T))
    lis2 <- dplyr::distinct(.data = lis1, SYMBOL, .keep_all = T)
    lis2$score <- apply(lis2, 1, function(x){myta[which(myta$name==x[1]),"score"]})#match(前, 后) 前在后多少位
    lis3 <- lis2[order(lis2$score, decreasing = T), c(2,3)]
    lis4 <- lis3$score; names(lis4) <- lis3$ENTREZID
    #
    gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = lis4,TERM2GENE = myg3,
                                                   maxGSSize = 2000, 
                                                   pAdjustMethod = "fdr", pvalueCutoff = 0.05))
    p_01 <- enrichplot::gseaplot2(x = gsea, 
                                  geneSetID = 1, title = "GSEA for myGene", pvalue_table = T)
    ggsave(filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA/theGenes.prerank.r8.eps",
           width = 14, height = 14, units = "cm", plot = p_01)
    pdf(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA/theGenes.prerank.r8.pdf", 
        width = 10, height = 10)
    print(p_01)
    dev.off()
  }
  #find the genes, [方案17: Broad peak, which Peak at its Promoter] [GSEA] [R] [确定下来的]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter", ignore.case = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter", ignore.case = T),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    myg1 <- intersect(tmp1$geneId, tmp2$geneId)
    myg2 <- suppressWarnings(clusterProfiler::bitr(wr01, fromType = "SYMBOL", 
                                                   toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
    myg3 <- data.frame(ont = rep("myGene", length(unique(myg2$ENTREZID))),
                       gen = unique(myg2$ENTREZID))
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene-2.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$score <- as.numeric(scale(myta$rank))
    lis1 <- suppressWarnings(clusterProfiler::bitr(geneID = myta$name, fromType = "SYMBOL",
                                                   toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop = T))
    lis2 <- dplyr::distinct(.data = lis1, SYMBOL, .keep_all = T)
    lis2$score <- apply(lis2, 1, function(x){myta[which(myta$name==x[1]),"score"]})#match(前, 后) 前在后多少位
    lis3 <- lis2[order(lis2$score, decreasing = T), c(2,3)]
    lis4 <- lis3$score; names(lis4) <- lis3$ENTREZID
    #
    gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = lis4,TERM2GENE = myg3,
                                                   maxGSSize = 2000, 
                                                   pAdjustMethod = "fdr", pvalueCutoff = 0.05))
    gsea@result
    p_01 <- enrichplot::gseaplot2(x = gsea,
                                  geneSetID = 1, title = "GSEA for myGene", pvalue_table = T)
    ggsave(filename = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA/theGenes.prerank.r9.eps",
           width = 14, height = 14, units = "cm", plot = p_01)
    pdf(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEA/theGenes.prerank.r9.pdf", 
        width = 10, height = 10)
    print(p_01)
    dev.off()
  }
  #
  #
  #
  #find the genes, [方案14: Broad peak, all peaks] [试验] [GSEA]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp1 <- tmp1[grep(tmp1$annotation, pattern="Promoter", ignore.case = T),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    tmp2 <- tmp2[grep(tmp2$annotation, pattern="Promoter", ignore.case = T),]
    wr01 <- data.frame(theGenes = intersect(tmp1$geneId, tmp2$geneId))
    wr01 <- as.data.frame(x = t(wr01))
    write.table(x = wr01, 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt",
                sep = "\t", col.names = F, quote = F)
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta <- myta[order(myta$rank, decreasing = T),]
    #myta$score <- as.numeric(scale(myta$rank))
    myta$score <- as.numeric(scale(seq(from = -10, to = 10, by=20/length(myta$name))[1:55353]))
    write.table(x = myta[,c("name","score")], 
                file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk",
                sep = "\t", row.names = F, col.names = F, quote = F)
    #
    paste0("gseapy prerank -r ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/tagsMERVL.rnk ",
           "-g /ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/theG.gmt -p 10 -o ",
           "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/enrichGSEAv2 --max-size 2000")
  }
  #find the genes, [方案14: Broad peak, all peaks] [试验] [Boxplot]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- read.csv(file = paste0(path, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
    ne01 <- unique(as.character(intersect(tmp1$geneId, tmp2$geneId)))
    #gf01 <- import(con = inx4); gf01 <- gf01[which(gf01$type=="gene")]
    #ne01 <- intersect(ne01, gf01$gene_name[which(gf01$gene_biotype=="protein_coding")])
    #
    #
    bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
    myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
    myta$mark <- NA
    myta$mark[which(myta$name %in% ne01)] <- "vennPeak"
    myta$mark[which(!(myta$name %in% ne01))] <- "other"
    myta$mark <- factor(x = myta$mark, levels = c("vennPeak","other"))
    ggplot(data = myta) +
      geom_boxplot(aes(x = mark, y = log2(rank +0.1), fill = mark)) +
      geom_rug() +
      theme_classic() +
      labs(x="Gene Set", y="Log2 (tags rank + 0.1)", title = "this is title") +
      theme(plot.margin = margin(1,1,1,1,unit = "cm"),
            plot.title = element_text(family = "serif", size = 15),
            axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
            axis.text.y = element_text(family = "serif", size = 12),
            axis.title = element_text(family = "serif", size = 12)) +
      geom_signif(mapping = aes(x = mark, y = log2(rank+0.1)),
                  comparisons = list(c("vennPeak", "other")),
                  map_signif_level=T, textsize=6, test = wilcox.test, step_increase = 0.12)
    #
    #
    #
    myta <- myta[order(myta$rank, decreasing = T),]
    myta$name <- factor(x = myta$name, levels = myta$name)
    
    p_01 <- ggplot() +
      geom_point(data = myta,mapping = aes(x = name, y = log2(rank+0.1)), color="grey40") +
      geom_point(data = myta[which(myta$mark=="vennPeak"),],
                 mapping = aes(x = name, y = log2(rank+0.1)+1), color="cyan", size =0.01) +
      geom_hline(yintercept = log2(3+0.1), linetype = 2, color = "red", linewidth = 0.8) +
      theme_classic() +
      scale_x_discrete(expand = c(0.01,0)) +
      labs(y="Log2 (n1Num/n1NumExp +0.1)", x="Genes (all)", title = "This Is Title...GRID-seq n1n2") +
      theme(axis.title = element_text(family = "serif", size = 12),
            plot.title = element_text(family = "serif", size = 15),
            axis.ticks.x = element_blank(),
            #legend.position = c(0.7,0.7),
            #legend.text = element_text(family = "serif", size = 12),
            #legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 1),
            axis.line.y = element_line(linewidth = 1.2)) +
      geom_rug(data = myta[which(myta$mark=="vennPeak"),],
               aes(x = name, y = log2(rank+0.1)),alpha = 0.2,sides = "t")
    p_01  
    #scale_color_manual(values = c("steelblue","gray"), labels = c("rank >  3","rank <=3")) +
    #guides(color = guide_legend(title = "title")) #+
    #geom_text_repel(data = tst1[which(tst1$mark =="object"),],
    #                mapping = aes(x = name, y = log2(rank+1), label = name))
  }
  #find the genes, [方案Xx: no peak, based on fisher检验] [Venn]
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
    tmp1 <- openxlsx::read.xlsx(xlsxFile = paste0(path, "single/fisher/IP-1FisherGene.xlsx"))
    tmp1 <- tmp1[which(tmp1$greate < 0.05),]
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
    tmp2 <- openxlsx::read.xlsx(xlsxFile = paste0(path, "single/fisher/IP-1FisherGene.xlsx"))
    tmp2 <- tmp2[which(tmp2$greate < 0.05),]
    #
    tmp3 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/grid.mervlNum.gene.xlsx")
    tmp3 <- tmp3[grep(tmp3$name, pattern="mt-", invert = T), ]
    tmp3 <- tmp3[which(tmp3$rank >3 & tmp3$biotype == "protein_coding"),]
    #
    temp <- list()
    temp[["a"]] <- unique(tmp1$geneName)
    temp[["b"]] <- unique(tmp2$geneName)
    temp[["c"]] <- unique(tmp3$name)
    p_01 <- VennDiagram::venn.diagram(x = temp,
                                      main = "all genes, which Peak not the Distal、Downstream",
                                      sub  = "", #a,b are samples of pr, c,d are FBL
                                      category.names = c("pr-1", "FBL-1", "MERVL-tags rank >3"),
                                      fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                      filename = NULL, 
                                      main.fontfamily = "serif", sub.fontfamily = "serif", 
                                      cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                      print.mode = "raw",#c("percent", "raw"),
                                      units = "cm", height = 12, width = 12)
    dev.off();grid.draw(p_01)
    
    ne01 <- intersect(intersect(temp[["a"]],temp[["b"]]),temp[["c"]])
  }
  #交集的101基因, 在早期胚胎发育阶段的 expression pheatmap
  {
    #用2014 Deng的数据 [不行]
    cn01 <- read.table(file = "/ChIP_seq_1/aaSHY/rna2014Deng/myRun/count/stageALL.cntTable", sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
    cn01 <- cn01[,colnames(cn01)[grep(x = colnames(cn01), pattern = "gene", invert = T)]]
    prex <- unique(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7))
    prex <- prex[c(11,10,2,3,8,5,12,1,7,4,9,6)]
    colnames(cn01) <- paste(sapply(str_split(colnames(cn01), pattern = "\\."), "[", 7),sapply(str_split(colnames(cn01), pattern = "\\."), "[", 8), sep = "-")
    cn02 <- as.data.frame(apply(cn01, 2, function(x){x/sum(x)*1000000}))
    cn02$geneName <- rownames(cn02)
    cn03 <- cn02[grep(rownames(cn02), pattern=":", invert = T),]#invert = F, 那么输出TE的CPM
    df01 <- tidyr::gather(data = cn03, key = "sample", value = "cpm", -"geneName")
    df01$stage <- sapply(str_split(df01$sample, "-"), "[", 1)
    df02 <- dplyr::reframe(.data = df01, .by = c(stage, geneName), cpmExp = mean(cpm), cpmSD = sd(cpm)/sqrt(length(cpm)))
    df03 <- tidyr::spread(data = df02[,-4], key = "stage", value = "cpmExp")
    df03 <- df03[which(df03$geneName %in% ne01),]
    rownames(df03) <- df03$geneName; df03 <- df03[,-1]
    df03 <- df03[,c("oocyte","zygote","twoEarly","twoMid","twoLate","four","eight","blastoEarly","blastoMid","blastoLate")]
    #
    df03 <- df03[rowSums(df03)>0,]
    pheatmap(mat = df03,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             labels_col = colnames(df03), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    
  }
  {
    #用2020 pacbio corresponding illumina的数据 [可以]
    stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/count/star/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    colnames(stag) <- sapply(str_split(colnames(stag), "\\."), "[", 8)
    stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))
    stag <- stag[grep(x = rownames(stag), pattern = ":", invert = T), -7]
    stag <- stag[,c("oocyte","zygote","two","four","eight","blasto")]
    df03 <- stag[ne01,]#;df03 <- df03[rowSums(df03)==0,]; rownames(df03)
    df03 <- df03[rowSums(df03)>0,]
    pdf("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/peak/prFbl.promoterGene.v1.pdf",
        width = 10, height = 15)
    pheatmap(mat = df03,
             main="the heatmap of all Genes",
             scale = "row", cellwidth = 20, 
             show_rownames = T, show_colnames = T, cluster_row = T, cluster_cols = F,
             labels_col = colnames(df03), angle_col = "45",
             breaks = seq(-2,2,by=0.01), border = F,
             color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                        colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
    dev.off()
  }
  #这些基因有多少是与MERVL MT2融合的?
  {
    ne02 <- rownames(df03)
    gf01 <- import(con = inx4)
    gf01 <- gf01[which(gf01$type=="gene")]
    tf01 <- import(con = inx5)
    tf01 <- tf01[tf01$gene_id %in% c("MERVL-int","MT2_Mm")]
    #组装的转录本中有哪些是MERVL嵌合的 [MERVL与转录本的外显子嵌合, 则算嵌合]
    {
      gf02 <- import(con = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf")
      gf02 <- gf02[which(gf02$type == "exon")]
      ov01 <- suppressWarnings(findOverlaps(query = gf02, subject = tf01, ignore.strand=T))
      mf02 <- gf02[unique(queryHits(ov01))]
      #对嵌合的转录本做注释
      gf02 <- import(con = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf")
      gf02 <- gf02[which(gf02$transcript_id %in% unique(mf02$transcript_id))]
      ov02 <- suppressWarnings(findOverlaps(query = gf02, subject = gf01, ignore.strand=T))
      theG <- gf01[unique(subjectHits(ov02))]
      intersect(ne02,unique(theG$gene_name))
    }
    #组装的转录本中有哪些是MERVL嵌合的 [MERVL与转录本嵌合, 则算嵌合]
    {
      gf02 <- import(con = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/stringTie/v_KO74.merge.gtf")
      gf02 <- gf02[which(gf02$type == "transcript")]
      ov01 <- suppressWarnings(findOverlaps(query = gf02, subject = tf01, ignore.strand=T))
      mf02 <- gf02[unique(queryHits(ov01))]
      #对嵌合的转录本做注释
      gf02 <- mf02
      ov02 <- suppressWarnings(findOverlaps(query = gf02, subject = gf01, ignore.strand=T))
      theG <- gf01[unique(subjectHits(ov02))]
      intersect(ne02,unique(theG$gene_name))
    }
  }
}
#############################################################################
#12、以FBL peak为主体, 看pr、Npm1的富集
{
  inx0 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/peak/IP-1-N_peaks.narrowPeak"
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/"
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  bw01 <- paste0(fath, "bw/bcp/", c("IP-1","IP-2"),"_Input-notSES.bw", collapse = " ")
  bw02 <- paste0(path, "bw/bcp/", c("IP-1","IP-2"),"_Input-notSES.bw", collapse = " ")
  bw03 <- "/ChIP_seq_2/aaSHY/Npm1/chipSEQ/bw/bcp/Npm1IP-Input.bw"
  shell <- paste0(path,"dumpROOM/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("computeMatrix reference-point --referencePoint center -a 1000 -b 1000 ",
                   "--missingDataAsZero -p 20 -R ", inx0,
                   " -S ",paste(bw01,bw02,bw03),
                   " -o ",path,"dumpROOM/coLocation_v4.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/coLocation_v4.mat.gz --samplesLabel ", 
                   paste0(c("FBL-1","FBL-2","pr-1","pr-2","NPM1-1"), collapse = " "),
                   " -out ",path,"dumpROOM/coLocation_v4.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"dumpROOM/run_a.log"," 2>&1 &"))
}
#############################################################################
#13、以FBL peak为主体, 看pr、H3K9me3、Setdb1、Trim28的富集
{
  inx0 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/peak/IP-1-N_peaks.narrowPeak"
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/"
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  bw01 <- paste0(fath, "bw/bcp/", c("IP-1","IP-2"),"_Input.bw", collapse = " ")
  bw02 <- paste0(path, "bw/bcp/", c("IP-1","IP-2"),"_Input.bw", collapse = " ")
  #bw03 <- paste0("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-1_Input.bw ",
  #               "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-3_Input.bw")
  bw03 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_1_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_2_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP051217_Control_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP069149_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw")
  bw04 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pr/publicData/TRIM28/bw/bcp/Trim28_vs.bw ",
                 "/ChIP_seq_2/aaSHY/pr/publicData/KAP1/bw/bcp/KAP1_vs.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP007651_Kap1_pseudocount_1_log2.bw")
  shell <- paste0(path,"dumpROOM/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("computeMatrix reference-point --referencePoint center -a 1000 -b 1000 ",
                   "--missingDataAsZero -p 20 -R ", inx0,
                   " -S ",paste(bw01,bw02,bw03,bw04,bw05),
                   " -o ",path,"dumpROOM/coLocation_v3.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/coLocation_v3.mat.gz --samplesLabel ", 
                   paste0(c("FBL-1","FBL-2","pr-1","pr-2",
                            "SRP001533-1","SRP001533-2",
                            "SRP051217","SRP064188","SRP069149","SRP094580",
                            "Setdb1-1","Setdb1-2","Trim28-1","Trim28-2","Trim28-3"),collapse = " "),
                   " -out ",path,"dumpROOM/coLocation_v3.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"dumpROOM/run_a.log"," 2>&1 &"))
}
################################################################################
#14、取所有Alu的区域，看看pr是否有富集
{
  shell <- paste0(path,"single/src/run_e.sh")
  cat("#!/bin/bash\n", file = shell)
  name <- c("MERVL-int::ERVL::LTR")
  bed1 <- list.files("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily", pattern = "Alu", full.names = T)
  bed2 <- paste0("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/",name,".bed")
  bw01 <- paste("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single//bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single//bw/bcp/IP-2_Input.bw", collapse = " ")
  cmd_01 <- paste0("computeMatrix reference-point --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -p 10 -R ", c(bed1,bed2),
                   " -S ",bw01,
                   " -o ",path,"single/dumpROOM/enrichAlu/TE.pr.",basename(c(bed1,bed2)),
                   ".mat.gz","\n")
  cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "single/dumpROOM/enrichAlu/TE.pr.",basename(c(bed1,bed2)),".mat.gz ",
                   "--samplesLabel pr-1 pr-2 ",
                   "--regionsLabel ",
                   gsub(x = basename(c(bed1,bed2)), pattern="::SINE.bed|::LTR.bed", replacement = ""),
                   " -out ",path,"single/dumpROOM/enrichAlu/TE.pr.",basename(c(bed1,bed2)),".pdf","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"single/Log/run_e.log"," 2>&1 &"))
}
################################################################################
#15、取所有MERVL-int的区域，看看pr是否有富集
{
  shell <- paste0(path,"single/src/run_e.sh")
  cat("#!/bin/bash\n", file = shell)
  name <- c("MERVL-int::ERVL::LTR")
  bed2 <- paste0("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/",name,".bed")
  bw01 <- paste("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/bw/bcp/IP-1_Input-notSES.bw",
                "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/bw/bcp/IP-2_Input-notSES.bw", collapse = " ")
  cmd_01 <- paste0("computeMatrix reference-point --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -p 10 -R ",bed2,
                   " -S ",bw01,
                   " -o ",path,"single/dumpROOM/enrichMERVL/TE.pr.",basename(c(bed2)),
                   ".mat.gz","\n")
  cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "single/dumpROOM/enrichMERVL/TE.pr.",basename(c(bed2)),".mat.gz ",
                   "--samplesLabel pr-1 pr-2 ",
                   "--regionsLabel ",
                   gsub(x = basename(c(bed1,bed2)), pattern="::SINE.bed|::LTR.bed", replacement = ""),
                   " -out ",path,"single/dumpROOM/enrichMERVL/TE.pr.",basename(c(bed1,bed2)),".pdf","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"single/Log/run_e.log"," 2>&1 &"))
}
################################################################################
#15、pr、FBL ChIP-seq在MERVL融合基因上的富集信号 [效果不理想]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/"
  bw01 <- list.files(path = path, pattern = ".bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")]
  {
    bed1 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_a.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix reference-point --referencePoint center ",
                       "-p 10 -a 5000 -b 5000 --missingDataAsZero ",
                       "-R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel FBL-1 FBL-2 pr-1 pr-2 ",
                       "-out ",path,"zempROOM/prFBL-IP-VS-",basename(inx0),".r1.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_a.log"," 2>&1 &"))
  }
}
################################################################################
#15、pr、FBL ChIP-seq在MERVL融合基因上的富集信号 [bamCompare不要SES]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/"
  shell <- paste0(path,"zempROOM/run_b.sh")
  cat("#!/bin/bash\n", file = shell)
  for (i in c("pr","FBL")) {
    cmd_01 <- paste0("bamCompare -p 10 --scaleFactorsMethod readCount -b1 ",path,i,
                     "/single/bam/csem/IP-",c(1,2),".sorted.bam ", "-b2 ",path,i,
                     "/single/bam/csem/Input-",c(1,2),".sorted.bam ",
                     "-o ",path,i,
                     "/single/bw/bcp/IP-",c(1,2),"_Input-notSES.bw","\n")
    for (i in cmd_01) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"zempROOM/run_b.log"," 2>&1 &"))
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,3)]
  {
    bed1 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 500 -b 500 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel FBL-1 pr-1 ",
                       "-out ",path,"zempROOM/prFBL-IP-VS-",basename(inx0),".r8.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Review 
#
#
#
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#01、以pr peak为主体, 看pr、Fbl的富集
{
  inx0 <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/peak/IP-1-N_peaks.narrowPeak"
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/"
  fath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/single/"
  bw01 <- paste0(fath, "bw/bcp/IP-1_Input-notSES.bw", collapse = " ")
  bw02 <- paste0(path, "bw/bcp/IP-1_Input-notSES.bw", collapse = " ")
  shell <- paste0(path,"dumpROOM/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("computeMatrix reference-point --referencePoint center -a 1000 -b 1000 ",
                   "--missingDataAsZero -p 20 -R ", inx0,
                   " -S ",paste(bw02,bw01),
                   " -o ",path,"dumpROOM/coLocation_v4.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/coLocation_v4.mat.gz --samplesLabel ", 
                   paste0(c("pr-1","FBL-1"), collapse = " "),
                   " -out ",path,"dumpROOM/review-coLocation_v1.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"dumpROOM/run_a.log"," 2>&1 &"))
}
#-------------------------------------------------------------------------------
#02、How significant is the positive correlation? P-value? [Fig.5a]
# 安装必要包
if (!require("reticulate")) install.packages("reticulate")
library(reticulate)

# 加载 Python 库
np <- import("numpy")
sp_stats <- import("scipy.stats")

# 读取您的 npz 文件
npz_file <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/single/dumpROOM/pearsonCorre_genome.npz"
npz_data <- np$load(npz_file)

# 提取矩阵
data_matrix <- npz_data[["matrix"]]

# 获取样本名称（从矩阵中提取或手动设置）
sample_names <- colnames(data_matrix)
tags <- c("pr-1","FBL-1","GSE46536_G9a","GSE166041_TRIM28","GSE18371_SETDB1")
if (is.null(sample_names)) {
  sample_names <- tags
}

# 应用异常值移除（复现 --removeOutliers）
remove_outliers <- function(x) {
  qnt <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- qnt[2] - qnt[1]
  lower <- qnt[1] - 1.5 * iqr
  upper <- qnt[2] + 1.5 * iqr
  x[x < lower | x > upper] <- NaN
  return(x)
}

# 对每列应用异常值移除
data_clean <- apply(data_matrix, 2, remove_outliers)

# 初始化结果矩阵
n_samples <- ncol(data_clean)
cor_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)
p_matrix <- matrix(NA, nrow = n_samples, ncol = n_samples)

# 计算每对样本的相关性（使用安全的提取方法）
for (i in 1:n_samples) {
  for (j in i:n_samples) {
    # 获取两列数据
    x <- data_clean[, i]
    y <- data_clean[, j]
    
    # 移除同时为NaN的行
    valid_idx <- !(is.nan(x) | is.nan(y))
    x_clean <- x[valid_idx]
    y_clean <- y[valid_idx]
    
    # 确保有足够的数据点
    if (length(x_clean) > 2 && length(y_clean) > 2) {
      # 计算Spearman相关（与plotCorrelation相同）
      result <- sp_stats$spearmanr(x_clean, y_clean)
      
      # 安全提取结果 - 处理所有可能的返回格式
      if (py_has_attr(result, "correlation")) {
        # 最新版本 (1.10+)
        cor_value <- result$correlation
        p_value <- result$pvalue
      } else if (py_has_attr(result, "statistic")) {
        # 中等版本 (0.19-1.9)
        cor_value <- result$statistic
        p_value <- result$pvalue
      } else if (py_has_attr(result, "__getitem__")) {
        # 旧版本 (元组格式)
        cor_value <- py_to_r(result[0])
        p_value <- py_to_r(result[1])
      } else {
        stop("无法识别的返回格式: ", class(result))
      }
      
      # 存储结果
      cor_matrix[i, j] <- cor_value
      cor_matrix[j, i] <- cor_value
      p_matrix[i, j] <- p_value
      p_matrix[j, i] <- p_value
    } else {
      # 数据不足时设为NA
      cor_matrix[i, j] <- NA
      cor_matrix[j, i] <- NA
      p_matrix[i, j] <- NA
      p_matrix[j, i] <- NA
    }
  }
}

# 设置行列名称
rownames(cor_matrix) <- sample_names
colnames(cor_matrix) <- sample_names
rownames(p_matrix) <- sample_names
colnames(p_matrix) <- sample_names

# 打印结果
print("Spearman相关系数矩阵 (与plotCorrelation一致):")
print(cor_matrix)

print("P值矩阵:")
print(p_matrix)

# 保存结果
#write.csv(cor_matrix, "exact_correlation_matrix.csv")
#write.csv(p_matrix, "exact_pvalue_matrix.csv")

# 可视化结果
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

# 绘制相关系数热图
pheatmap(cor_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,
         number_format = "%.2f",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Spearman Correlation (SciPy Calculation)")

# 绘制p值热图（-log10转换）
log_p_matrix <- -log10(p_matrix)
log_p_matrix[is.infinite(log_p_matrix)] <- max(log_p_matrix[is.finite(log_p_matrix)], na.rm = TRUE) + 1

# 创建显著性标记矩阵
sig_matrix <- matrix("", nrow = nrow(p_matrix), ncol = ncol(p_matrix))
sig_matrix[p_matrix < 0.001] <- "***"
sig_matrix[p_matrix >= 0.001 & p_matrix < 0.01] <- "**"
sig_matrix[p_matrix >= 0.01 & p_matrix < 0.05] <- "*"

# 绘制热图
pheatmap(log_p_matrix,
         color = colorRampPalette(c("white", "yellow", "red"))(100),
         display_numbers = sig_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Significance of Correlation (-log10(p-value))")



