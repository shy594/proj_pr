#write by ~~~ at 2024.01.23
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
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
  for (o in c("KO", "RN", "WT")) {
   path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
   path = paste0(path, o, "/"); rm(o)
   for (i in c("src","trim","Log","peak","dumpROOM","DESeq2","rRNA/csem","mervl/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
     if (dir.exists(paths = paste0(path, i))==F) {
       dir.create(path = paste0(path,i), recursive = T)
     }
   }
   rm(i)
   shell <- paste0(path,"src/run_a.sh")
   cat("#!/bin/bash\n", file = shell)
   prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R"), "[", 1))
   cmd_01 <- paste0("zcat ", path, "fq/", prex, "-R1.fq.gz ",path, "fq/", prex, "-R2.fq.gz |",
                    "gzip >", path, "fq/", prex, ".fq.gz", "\n")
   cmd_02 <- paste0("trim_galore --phred33 --fastqc --illumina ",path,
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
   for (i in cmd_01) {cat(i, append = T, file = shell)}
   for (i in cmd_02) {cat(i, append = T, file = shell)}
   for (i in cmd_04) {cat(i, append = T, file = shell)}
   for (i in cmd_05) {cat(i, append = T, file = shell)}
   for (i in cmd_06) {cat(i, append = T, file = shell)}
   for (x in seq(2)) {
     zrex <- prex[c(x,x+2)]
     for (q in zrex[2]) {
       cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
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
                      "--project h3k9me3", zrex[2], " --outdir ",paste0(path,"count/"),"\n")
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
}
################################################################################
#04、suppl
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  o="WT";path = paste0(path, o, "/"); rm(o)
  shell <- paste0(path,"src/run_suppl.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R"), "[", 1))
  prex <- prex[grep(x = prex, pattern = "gz", invert = T)]
  inxM = "/Reference/aaSHY/zOther/TEconsensus/RLTR1B/bt2/RLTR1B"
  cmd_17 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inxM," -U ", path, "trim/", prex, "_trimmed.fq.gz |",
                   "samtools view -F 4 -b > ", path, "otherTE/RLTR1B/",prex,".bam","\n")
  cmd_18 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"otherTE/RLTR1B/",prex,".bam")," 140 ",
                   paste0(path,"otherTE/RLTR1B/csem/",prex),"\n")
  cmd_19 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"otherTE/RLTR1B/csem/",prex,".sorted.bam -o "),
                   paste0(path,"otherTE/RLTR1B/",prex,".bw\n"))
  for (i in cmd_17) {cat(i, append = T, file = shell)}
  for (i in cmd_18) {cat(i, append = T, file = shell)}
  for (i in cmd_19) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_suppl.log"," 2>&1 &"))
}
################################################################################
#05、TE Binding situation [适用于ChIP-seq]
{
  for (o in c("KO", "RN", "WT")) {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    for (x in seq(2)) {
      zrex <- prex[c(x,x+2)]; rm(x)
      #Observed
      red1 <- read.table(file = paste0(path,"count/h3k9me3",zrex[2],".cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
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
        openxlsx::write.xlsx(x = red3, asTable = T, file = paste0(path,"fisher/",zrex[2],"FisherTE.xlsx"))
      }; rm(m)
    }
  }
}
################################################################################
#06、基因的Binding situation [适用于ChIP-seq]
{
  for (o in c("KO", "RN", "WT")) {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    for (x in seq(2)) {
      zrex <- prex[c(x,x+2)]; rm(x)
      twoc <- read.table(file = inxB, sep = "\t")
      #Observed
      red1 <- read.table(file = paste0(path,"count/h3k9me3",zrex[2],".cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
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
        openxlsx::write.xlsx(x = red3, asTable = T, file = paste0(path,"fisher/",zrex[2],"FisherGene.xlsx"))
      }
    }; rm(m)
  }
}
################################################################################
#07、peak anno
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
  for (o in c("KO", "RN", "WT")) {
    print(o)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    for (x in seq(2)) {
      zrex <- prex[c(x,x+2)]; rm(x)
      pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-B_peaks.broadPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb1))
      pdf(file = paste0(path, "peak/", zrex[2], "-B_peaks-Anno.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk01), paste0(path, "peak/",zrex[2],"-B_peaks-Anno.csv"), quote = F, row.names = F)
      pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-N_peaks.narrowPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb1))
      pdf(file = paste0(path, "peak/", zrex[2], "-N_peaks-Anno.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk02), paste0(path, "peak/",zrex[2],"-N_peaks-Anno.csv"), quote = F, row.names = F)
      pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-B_peaks.broadPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb2))
      pdf(file = paste0(path, "peak/", zrex[2], "-B_peaks-Anno-2C.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk03), paste0(path, "peak/",zrex[2],"-B_peaks-Anno-2C.csv"), quote = F, row.names = F)
      pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-N_peaks.narrowPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb2))
      pdf(file = paste0(path, "peak/", zrex[2], "-N_peaks-Anno-2C.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk04), paste0(path, "peak/",zrex[2],"-N_peaks-Anno-2C.csv"), quote = F, row.names = F)
    }
  }
}
################################################################################
#08、peak nearest gene
{
  gtf3 <- rtracklayer::import(con = inx4, format = "gtf")
  gtf3 <- gtf2[which(gtf2$type=="gene")]
  for (o in c("RN", "KO", "WT")) {
    print(o)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    for (x in seq(2)) {
      print(x)
      zrex <- prex[c(x,x+2)]; rm(x)
      tmp1 <- read.table(file = paste0(path, "peak/",zrex[2],"-B_peaks.broadPeak"))
      colnames(tmp1) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue")
      tmp1$nearestGen <- NA
      tmp1$theDistanc <- NA
      tmp2 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp1, keep.extra.columns = T)
      for (i in seq_along(tmp2)) {
        print(i)
        tmp3 <- tmp2[i]
        ov01 <- suppressWarnings(GenomicRanges::distanceToNearest(x = tmp3, subject = gtf3))
        tmp1$nearestGen[i] <- paste0(gtf3[unique(subjectHits(ov01))]$gene_name, collapse = ",")
        tmp1$theDistanc[i] <- ov01@elementMetadata$distance
      }; rm(i, ov01)
      write.table(x = tmp1, file = paste0(path, "peak/",zrex[2],"-B_peaks-nearGene.broadPeak.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
      tmp4 <- read.table(file = paste0(path, "peak/", zrex[2], "-N_peaks.narrowPeak"))
      colnames(tmp4) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue", "offset")
      tmp4$nearestGen <- NA
      tmp4$theDistanc <- NA
      tmp5 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp4, keep.extra.columns = T)
      for (i in seq_along(tmp5)) {
        print(i)
        tmp6 <- tmp5[i]
        ov02 <- suppressWarnings(GenomicRanges::distanceToNearest(x = tmp6, subject = gtf3))
        tmp4$nearestGen[i] <- paste0(gtf3[unique(subjectHits(ov02))]$gene_name, collapse = ",")
        tmp4$theDistanc[i] <- ov02@elementMetadata$distance
      }; rm(i, ov02)
      write.table(x = tmp4, file = paste0(path, "peak/",zrex[2],"-N_peaks-nearGene.narrowPeak.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
    }
  }
}
################################################################################
#09、取出不同样本的位置重叠的ov peak
{
  for (o in c("KO", "RN", "WT")) {
    print(o)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    tmp1 <- read.table(file = paste0(path, "peak/",prex[3],"-B_peaks.broadPeak"))
    tmp2 <- read.table(file = paste0(path, "peak/",prex[4],"-B_peaks.broadPeak"))
    tmp3 <- read.table(file = paste0(path, "peak/",prex[3],"-N_peaks.narrowPeak"))
    tmp4 <- read.table(file = paste0(path, "peak/",prex[4],"-N_peaks.narrowPeak"))
    colnames(tmp1) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue")
    colnames(tmp2) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue")
    colnames(tmp3) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue", "offset")
    colnames(tmp4) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue", "offset")
    tmp5 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp1, keep.extra.columns = T)
    tmp6 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp2, keep.extra.columns = T)
    tmp7 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp3, keep.extra.columns = T)
    tmp8 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp4, keep.extra.columns = T)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = tmp5, subject = tmp6))
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = tmp7, subject = tmp8))
    wr01 <- tmp1[unique(queryHits(ov01)),]
    wr02 <- tmp2[unique(subjectHits(ov01)),]
    wr03 <- tmp3[unique(queryHits(ov02)),]
    wr04 <- tmp4[unique(subjectHits(ov02)),]
    write.table(x = wr01, file = paste0(path, "peak/",prex[3],"-B_peaks_ov.broadPeak"), quote = F, sep = "\t", col.names = F, row.names = F)
    write.table(x = wr02, file = paste0(path, "peak/",prex[4],"-B_peaks_ov.broadPeak"), quote = F, sep = "\t", col.names = F, row.names = F)
    write.table(x = wr03, file = paste0(path, "peak/",prex[3],"-N_peaks_ov.narrowPeak"), quote = F, sep = "\t", col.names = F, row.names = F)
    write.table(x = wr04, file = paste0(path, "peak/",prex[4],"-N_peaks_ov.narrowPeak"), quote = F, sep = "\t", col.names = F, row.names = F)
  }
}
################################################################################
#10、peak anno for ov peak
{
  {
    gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
    gtf1 <- as.data.frame(gtf1)
    colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
    gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
    txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
    twoC <- read.table(inxB)
    gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
    txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
  }
  #
  for (o in c("KO", "RN", "WT")) {
    print(o)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    for (x in seq(2)) {
      zrex <- prex[c(x,x+2)]; rm(x)
      pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-B_peaks_ov.broadPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb1))
      pdf(file = paste0(path, "peak/", zrex[2], "-B_peaks_ov-Anno.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk01), paste0(path, "peak/",zrex[2],"-B_peaks_ov-Anno.csv"), quote = F, row.names = F)
      pk02 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-N_peaks_ov.narrowPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb1))
      pdf(file = paste0(path, "peak/", zrex[2], "-N_peaks_ov-Anno.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk02, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk02), paste0(path, "peak/",zrex[2],"-N_peaks_ov-Anno.csv"), quote = F, row.names = F)
      pk03 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-B_peaks_ov.broadPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb2))
      pdf(file = paste0(path, "peak/", zrex[2], "-B_peaks_ov-Anno-2C.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk03, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk03), paste0(path, "peak/",zrex[2],"-B_peaks_ov-Anno-2C.csv"), quote = F, row.names = F)
      pk04 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",zrex[2],"-N_peaks_ov.narrowPeak"), 
                                                        tssRegion = c(-2000, 2000), TxDb = txb2))
      pdf(file = paste0(path, "peak/", zrex[2], "-N_peaks_ov-Anno-2C.pdf"), width = 8, height = 8)
      ChIPseeker::plotAnnoPie(pk04, legend.position = "rightside")
      dev.off()
      write.csv(as.data.frame(pk04), paste0(path, "peak/",zrex[2],"-N_peaks_ov-Anno-2C.csv"), quote = F, row.names = F)
    }
  }
}
################################################################################
#11、peak nearest gene for ov peak
{
  for (o in c("KO", "RN", "WT")) {
    print(o)
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    path = paste0(path, o, "/"); rm(o)
    prex <- unique(sapply(str_split(list.files(path = paste0(path, "fq")), pattern = "-R|.fq"), "[", 1))
    for (x in seq(2)) {
      print(x)
      zrex <- prex[c(x,x+2)]; rm(x)
      tmp1 <- read.table(file = paste0(path, "peak/",zrex[2],"-B_peaks_ov.broadPeak"))
      colnames(tmp1) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue")
      tmp1$nearestGen <- NA
      tmp1$theDistanc <- NA
      tmp2 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp1, keep.extra.columns = T)
      for (i in seq_along(tmp2)) {
        print(i)
        tmp3 <- tmp2[i]
        ov01 <- suppressWarnings(GenomicRanges::distanceToNearest(x = tmp3, subject = gtf3))
        tmp1$nearestGen[i] <- paste0(gtf3[unique(subjectHits(ov01))]$gene_name, collapse = ",")
        tmp1$theDistanc[i] <- ov01@elementMetadata$distance
      }; rm(i, ov01)
      write.table(x = tmp1, file = paste0(path, "peak/",zrex[2],"-B_peaks_ov-nearGene.broadPeak.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
      tmp4 <- read.table(file = paste0(path, "peak/",zrex[2],"-N_peaks_ov.narrowPeak"))
      colnames(tmp4) <- c("chr", "start", "end", "name", "score", "strand", "signal", "pValue", "qValue", "offset")
      tmp4$nearestGen <- NA
      tmp4$theDistanc <- NA
      tmp5 <- GenomicRanges::makeGRangesFromDataFrame(df = tmp4, keep.extra.columns = T)
      for (i in seq_along(tmp5)) {
        print(i)
        tmp6 <- tmp5[i]
        ov02 <- suppressWarnings(GenomicRanges::distanceToNearest(x = tmp6, subject = gtf3))
        tmp4$nearestGen[i] <- paste0(gtf3[unique(subjectHits(ov02))]$gene_name, collapse = ",")
        tmp4$theDistanc[i] <- ov02@elementMetadata$distance
      }; rm(i, ov02)
      write.table(x = tmp4, file = paste0(path, "peak/",zrex[2],"-N_peaks_ov-nearGene.narrowPeak.txt"), quote = F, sep = "\t", col.names = T, row.names = F)
    }
  }
}
################################################################################
#12、H3K9me3在MERVL等TE上的富集
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  rep1 <- list.files(path = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily", full.names = T)
  shell <- paste0(path,"zempROOM/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  for (i in paste0(c("IAPEz-int","IAPEy-int","MT2_Mm","MERVL-int"),":")) {
    inx0 <- rep1[grep(x = rep1, pattern = i)]
    ff01 <- read.table(file = inx0, sep = "\t")
    cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 3000 -b 3000 -m 2000 ",
                     "--missingDataAsZero -p 20 -R ", inx0, " -S ",
                     paste0(path, c("KO","RN","WT"),"/bw/bcp/IP-1_Input-notSES.bw", collapse = " ")," ",
                     paste0(path, c("KO","RN","WT"),"/bw/bcp/IP-3_Input-notSES.bw", collapse = " ")," ",
                     "-o ",path,"zempROOM/", basename(inx0),".mat.gz","\n")
    cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "zempROOM/",basename(inx0),".mat.gz ",
                     "--samplesLabel ", paste0(c("KO","RN","WT"), "IP-1_Input", collapse = " ")," ",
                     paste0(c("KO","RN","WT"), "IP-3_Input", collapse = " ")," ",
                     "-out ",path,"zempROOM/",basename(inx0),"-r2.pdf","\n")
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_12) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"zempROOM/run_a.log"," 2>&1 &"))
}
{
  inx0 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/mervl_zyw.bed"
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 3000 -b 3000 -m 2000 ",
                   "--missingDataAsZero -p 20 -R ", inx0, " -S ",
                   paste0(path, c("KO","RN","WT"),"/bw/bcp/IP-1_Input.bw", collapse = " ")," ",
                   paste0(path, c("KO","RN","WT"),"/bw/bcp/IP-3_Input.bw", collapse = " ")," ",
                   "-o ",path,"zempROOM/", basename(inx0),".mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "zempROOM/",basename(inx0),".mat.gz ",
                   "--samplesLabel ", paste0(c("KO","RN","WT"), "IP-1_Input", collapse = " ")," ",
                   paste0(c("KO","RN","WT"), "IP-3_Input", collapse = " ")," ",
                   "-out ",path,"zempROOM/",basename(inx0),".pdf","\n")
}
################################################################################
#13、H3K9me3在MERVL融合基因上的富集
{
  inx0 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 3000 -b 3000 -m 2000 ",
                   "--missingDataAsZero -R ", inx0, " -S ",
                   paste0(path, c("KO","RN","WT"),"/bw/bcp/IP-1_Input.bw", collapse = " ")," ",
                   paste0(path, c("KO","RN","WT"),"/bw/bcp/IP-3_Input.bw", collapse = " ")," ",
                   "-o ",path,"zempROOM/", basename(inx0),".mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "zempROOM/",basename(inx0),".mat.gz ",
                   "--samplesLabel ", paste0(c("KO","RN","WT"), "IP-1_Input", collapse = " ")," ",
                   paste0(c("KO","RN","WT"), "IP-3_Input", collapse = " ")," ",
                   "-out ",path,"zempROOM/",basename(inx0),".pdf","\n")
}
################################################################################
#14、H3K9me3在1165个基因上的富集 [基于bamCoverage的]
{
  inx0 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
  cmd_09 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                   "--missingDataAsZero -R ", inx0, " -S ",
                   paste0(path, c("KO","WT"),"/bw/bcv/IP-1.bw", collapse = " ")," ",
                   paste0(path, c("KO","WT"),"/bw/bcv/IP-3.bw", collapse = " ")," ",
                   paste0(path, c("KO","WT"),"/bw/bcv/Input-1.bw", collapse = " ")," ",
                   paste0(path, c("KO","WT"),"/bw/bcv/Input-1.bw", collapse = " ")," ",
                   "-o ",path,"zempROOM/", basename(inx0),".mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "zempROOM/",basename(inx0),".mat.gz ",
                   "--samplesLabel ", 
                   paste0("KO", c("IP-1","IP-3","Input-1","Input-3"), collapse = " ")," ",
                   paste0("WT", c("IP-1","IP-3","Input-1","Input-3"), collapse = " ")," ",
                   "-out ",path,"zempROOM/K9me3-KO-WT-",basename(inx0),".pdf","\n")
}
################################################################################
#15、H3K9me3在1165个基因上的富集 [基于bamCompare的]
{
  cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/KO/bam/csem/IP-1.sorted.bam ",
                   "-b2 ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bam/csem/IP-1.sorted.bam ",
                   "-o ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/zempROOM/KO-WT-1.bw","\n")
  cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/KO/bam/csem/IP-3.sorted.bam ",
                   "-b2 ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/WT/bam/csem/IP-3.sorted.bam ",
                   "-o ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/zempROOM/KO-WT-3.bw","\n")
}
{
  inx0 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 10000 ",
                   "--missingDataAsZero -R ", inx0, " -S ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/zempROOM/KO-WT-1.bw ",
                   "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/zempROOM/KO-WT-3.bw ",
                   "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                   "--samplesLabel KO-WT-1 KO-WT-3 ",
                   "-out ",path,"zempROOM/K9me3-KO-WT-IP-VS-",basename(inx0),".pdf","\n")
}
################################################################################
#16、H3K9me3在1165个基因上的富集 [基于bamCompare的; 带SES]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "_Input.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,3,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
    #bed2 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    #bed3 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
    #bed4 <- "/Reference/aaSHY/BED/special/mervl_zyw.bed"
    #beds <- c(bed1,bed2,bed3,bed4)
    #beds <- c("/Reference/aaSHY/BED/TEsubFamily/IAPEz-int::ERVK::LTR.bed")
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_k.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".SES.r7-scale.pdf","\n")
      cmd_03 <- paste0("#computeMatrix reference-point --referencePoint TSS -p 20 ",
                       "-a 2000 -b 2000 --missingDataAsZero -R ", inx0, " ",
                       "-S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_04 <- paste0("#plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".SES.r7-tss.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
      for (i in cmd_03) {cat(i, append = T, file = shell)}
      for (i in cmd_04) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_k.log"," 2>&1 &"))
  }
}
################################################################################
#16、H3K9me3在1254个基因上的富集 [基于bamCompare的; 带SES]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "_Input.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,3,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
    #bed2 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    #bed3 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
    #bed4 <- "/Reference/aaSHY/BED/special/mervl_zyw.bed"
    #beds <- c(bed1,bed2,bed3,bed4)
    #beds <- c("/Reference/aaSHY/BED/TEsubFamily/IAPEz-int::ERVK::LTR.bed")
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_k.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".SES.r12-scale.pdf","\n")
      cmd_03 <- paste0("computeMatrix reference-point --referencePoint TSS -p 20 ",
                       "-a 2000 -b 2000 --missingDataAsZero -R ", inx0, " ",
                       "-S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_04 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".SES.r12-tss.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
      for (i in cmd_03) {cat(i, append = T, file = shell)}
      for (i in cmd_04) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_k.log"," 2>&1 &"))
  }
}
################################################################################
#17、H3K9me3 WT KO RN在MERVL-int上的富集 [bamCompare, not SES]--------------------------------------------------
#####WT KO用的bamCompare默认参数
#####WT RN用的bamCompare SES参数
#####WT RN KO画在一起，用 SES 参数----------------------------------------------center、3kb
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bwbw <- list.files(path = path, pattern = "_Input-notSES.bw", full.names = T, recursive = T)
  bw01 <- bwbw[grep(x = bwbw, pattern = "bcp")][c(1,5)]
  bw02 <- bwbw[grep(x = bwbw, pattern = "bcp")][c(3,5)]
  bw03 <- bwbw[grep(x = bwbw, pattern = "bcp")][c(1,3,5)]
  {
    bed1 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
    bed2 <- "/Reference/aaSHY/BED/special/mervl_zyw.bed"
    beds <- c(bed1,bed2)
    shell <- paste0(path,"/zempROOM/run_o.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("#computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"/zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("#plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"/zempROOM/K9me3-WT_KO-",basename(inx0),".r7.pdf","\n")
      cmd_03 <- paste0("#computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw02, collapse = " ")," ",
                       "-o ",path,"/zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_04 <- paste0("#plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel RN-1 WT-1 ",
                       "-out ",path,"/zempROOM/K9me3-WT_RN-",basename(inx0),".r7-SES.pdf","\n")
      cmd_05 <- paste0("#computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw03, collapse = " ")," ",
                       "-o ",path,"/zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_06 <- paste0("#plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"/zempROOM/K9me3-WT_RN-",basename(inx0),".r7-SES-Union.pdf","\n")
      cmd_07 <- paste0("computeMatrix reference-point --referencePoint center ",
                       "-p 20 -a 3000 -b 3000 --missingDataAsZero ",
                       "-R ", inx0, " -S ",paste(bw03, collapse = " ")," ",
                       "-o ", path,"/zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_08 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"/zempROOM/K9me3-WT_RN-",basename(inx0),".r7-notSES-Union-center.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
      for (i in cmd_03) {cat(i, append = T, file = shell)}
      for (i in cmd_04) {cat(i, append = T, file = shell)}
      for (i in cmd_05) {cat(i, append = T, file = shell)}
      for (i in cmd_06) {cat(i, append = T, file = shell)}
      for (i in cmd_07) {cat(i, append = T, file = shell)}
      for (i in cmd_08) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_s.log"," 2>&1 &"))
  }
}
################################################################################
#17、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare不要SES]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  shell <- paste0(path,"zempROOM/run_j.sh")
  cat("#!/bin/bash\n", file = shell)
  for (i in c("WT","KO","RN")) {
    cmd_01 <- paste0("bamCompare -p 10 --scaleFactorsMethod readCount -b1 ",path,i,
                     "/bam/csem/IP-",c(1,3),".sorted.bam ", "-b2 ",path,i,
                     "/bam/csem/Input-",c(1,3),".sorted.bam ",
                     "-o ",path,i,
                     "/bw/bcp/IP-",c(1,3),"_Input-notSES.bw","\n")
    for (i in cmd_01) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"zempROOM/run_j.log"," 2>&1 &"))
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
    bed2 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    bed3 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
    bed4 <- "/Reference/aaSHY/BED/special/mervl_zyw.bed"
    beds <- c(bed1,bed2,bed3,bed4)
    shell <- paste0(path,"zempROOM/run_l.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r4.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_l.log"," 2>&1 &"))
  }
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    beds <- c("/Reference/aaSHY/BED/TEsubFamily/IAPEz-int::ERVK::LTR.bed",
              "/Reference/aaSHY/BED/TEsubFamily/ETnERV-int::ERVK::LTR.bed")
    shell <- paste0(path,"zempROOM/run_m.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r4.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_m.log"," 2>&1 &"))
  }
}
################################################################################
#18、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare不要SES] [筛1165: 取top]
{
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  myta <- myta[order(myta$rank, decreasing = T),]
  myta[1:5,]
  ge01 <- import(con = inx4)
  ge01 <- ge01[which(ge01$type=="gene")]
  wr01 <- as.data.frame(ge01[which(ge01$gene_name %in% myta$name)])
  wr01 <- wr01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed", 
              sep = "\t", col.names = F, row.names = F, quote = F)
  bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.head200.gene.bed"
  bed2 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.tail200.gene.bed"
  {
    path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
    bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
    bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
    {
      beds <- c("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.head200.gene.bed",
                "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.tail200.gene.bed")
      shell <- paste0(path,"zempROOM/run_n.sh")
      cat("#!/bin/bash\n", file = shell)
      for (inx0 in beds) {
        cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                         "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                         "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
        cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                         "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                         "--samplesLabel KO-1 WT-1 ",
                         "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r4.pdf","\n")
        for (i in cmd_01) {cat(i, append = T, file = shell)}
        for (i in cmd_02) {cat(i, append = T, file = shell)}
      }
      print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
    }
  }
}
################################################################################
#19、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare要SES, 但是修改sample Length]
#效果仍然不好
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  shell <- paste0(path,"zempROOM/run_j.sh")
  cat("#!/bin/bash\n", file = shell)
  for (i in c("WT","KO","RN")) {
    cmd_01 <- paste0("bamCompare -p 10 --scaleFactorsMethod SES --sampleLength 1000 -b1 ",path,i,
                     "/bam/csem/IP-",c(1,3),".sorted.bam ", "-b2 ",path,i,
                     "/bam/csem/Input-",c(1,3),".sorted.bam ",
                     "-o ",path,i,
                     "/bw/bcp/IP-",c(1,3),"_Input-reviseSES.bw","\n")
    for (i in cmd_01) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"zempROOM/run_j.log"," 2>&1 &"))
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "reviseSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
    bed2 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    bed3 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
    bed4 <- "/Reference/aaSHY/BED/special/mervl_zyw.bed"
    beds <- c(bed1,bed2,bed3,bed4)
    shell <- paste0(path,"zempROOM/run_l.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r5.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_l.log"," 2>&1 &"))
  }
}
################################################################################
#20、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare不要SES] 
#[筛1165: CPM normalize后, 取pr或Fbl上调]
{
  red1 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/normCPM-1165.csv", row.names = 1)
  red2 <- red1
  red2 <- red2[which(rowMeans(red2[,c(1,2)])/rowMeans(red2[,c(3,4)]) >0),]
  red2 <- red2[which(rowMeans(red2[,c(5,6)])/rowMeans(red2[,c(7,8)]) >0),]
  ge01 <- import(con = inx4)
  ge01 <- ge01[which(ge01$type=="gene")]
  ge01 <- as.data.frame(ge01[which(ge01$gene_name %in% rownames(red2))])
  wr01 <- ge01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.normCPM.FblUp.prUp.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.normCPM.FblUp.prUp.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r5.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#21、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare不要SES] 
#[筛1165: CPM normalize后, 取pr上调]
{
  red1 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-fbl/normCPM-1165.csv", row.names = 1)
  red2 <- red1
  red2 <- red2[which(rowMeans(red2[,c(1,2)])/rowMeans(red2[,c(3,4)]) >0),]
  ge01 <- import(con = inx4)
  ge01 <- ge01[which(ge01$type=="gene")]
  ge01 <- as.data.frame(ge01[which(ge01$gene_name %in% rownames(red2))])
  wr01 <- ge01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.normCPM.prUp.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.normCPM.prUp.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r6.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#22、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare不要SES]
#[筛1165: pr、Fbl broad peak注释到的基因，和MERVL-RNA高频结合的基因（rank>3、rank>5）的交集]
{
  cath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  tmp1 <- read.csv(file = paste0(cath, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
  tmp1 <- unique(tmp1$geneId)
  cath = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/FBL/"
  tmp2 <- read.csv(file = paste0(cath, "single/peak/IP-1-allg-B_chipseekerAnno.csv"), header = T)
  tmp2 <- unique(tmp2$geneId)
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  tmp3 <- as.character(myta$name)
  #
  ne01 <-intersect(tmp1, tmp3)#pr Peak & 1165 = 535
  ne02 <- intersect(ne01, tmp2)#pr Peak & Fbl Peak & 1165 = 142
  #
  ge01 <- import(con = inx4)
  ge01 <- ge01[which(ge01$type=="gene")]
  ge01 <- as.data.frame(ge01[which(ge01$gene_name %in% ne02)])
  wr01 <- ge01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.PrmtFbl1Peak.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.prPeak.bed"
    bed2 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.prFbl1Peak.bed"
    beds <- c(bed1,bed2)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r5.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#23、选择1165对不对？
{
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
  myta <- myta[which(myta$rank >3),]
  myta <- myta[order(myta$uniNum, decreasing = T),]
  ne01 <- myta[1:500,"name"]
  #
  ge01 <- import(con = inx4)
  ge01 <- ge01[which(ge01$type=="gene")]
  ge01 <- as.data.frame(ge01[which(ge01$gene_name %in% ne01)])
  wr01 <- ge01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/top500.byNum.rank.Is3.tagMERVL.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/top500.byNum.rank.Is3.tagMERVL.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r5.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#24、H3K9me3在1165个基因上的富集 [基于bamCompare的; bamCompare不要SES] 
#by [center, not scale region]
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
    bed2 <- "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed"
    bed3 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
    bed4 <- "/Reference/aaSHY/BED/special/mervl_zyw.bed"
    beds <- c(bed1,bed2,bed3,bed4)
    shell <- paste0(path,"zempROOM/run_o.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix reference-point --referencePoint center ",
                       "-p 10 -a 5000 -b 5000 --missingDataAsZero ",
                       "-R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r6.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_o.log"," 2>&1 &"))
  }
}
################################################################################
#25、MERVL-融合基因、1165基因和2C基因重合
{
  fusg <- read.table(file = "/Reference/aaSHY/BED/special/mervl.2cell.fusionGene.bed", sep = "\t")
  twoc <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed", sep = "\t")
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  #
  temp <- list()
  temp[["a"]] <- unique(fusg$V4)
  temp[["b"]] <- unique(twoc$V4)
  temp[["c"]] <- unique(myta$name)
  p_01 <- VennDiagram::venn.diagram(x = temp,
                                    main = "the Genes",
                                    sub  = "", #a,b are samples of pr, c,d are FBL
                                    category.names = c("MERVL-Fus-1490", "2C-517", "MERVL-tags rank >3"),
                                    fill=c("#FD763F","#23BAC5","#B2DBB9"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)
  dev.off();grid.draw(p_01)
}
################################################################################
#26、MERVL-融合转录本、1165基因的重合 [得到159个基因]
{
  ge01 <- import(inx4)
  ge01 <- ge01[which(ge01$type == "gene")]
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  myta <- ge01[which(ge01$gene_name %in% myta$name)]
  fusg <- read.table(file = "/Reference/aaSHY/BED/special/mervl.2cell.fusionTranscript.bed", sep = "\t")
  fusg <- makeGRangesFromDataFrame(df = fusg, seqnames.field = "V1", start.field = "V2",
                                   end.field = "V3", strand.field = "V6")
  ov01 <- findOverlaps(query = ge01, subject = fusg, ignore.strand = T)
  tmp1 <- ge01[unique(queryHits(ov01))]
  length(unique(tmp1$gene_name))
  ov02 <- findOverlaps(query = myta, subject = fusg, ignore.strand = T)
  tmp2 <- as.data.frame(myta[unique(queryHits(ov02))])
  wr02 <- tmp2[,c(1,2,3,12,4,5)]
  wr02$width <- "."
  write.table(x = wr02, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/tagMERVL.mervlFus.159.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/tagMERVL.mervlFus.159.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r8.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#27、MERVL融合基因 & 1165基因
{
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.mervlNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  fusg <- read.table(file = "/Reference/aaSHY/BED/special/mervl.2cell.fusionGene.bed", sep = "\t")
  ne01 <- unique(c(myta$name, fusg$V4))
  ge01 <- import(inx4)
  ge01 <- ge01[which(ge01$type == "gene")]
  wr01 <- as.data.frame(ge01[which(ge01$gene_name %in% ne01)])
  wr01 <- wr01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, 
              file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/tagMERVL.mervlFus.2501.bed",
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/tagMERVL.mervlFus.2501.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_n.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2000 -b 2000 -m 10000 ",
                       "--missingDataAsZero -R ", inx0, " -S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),".mat.gz ",
                       "--samplesLabel KO-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r8.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#28、H3K9me3在MERVL binding tags上的富集
{
  ###red1 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/n1.mervl.DNA.bed",
  ###                   sep = "\t", header = F)
  ###red2 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/n2.mervl.DNA.bed",
  ###                   sep = "\t", header = F)
  ###reds <- distinct(rbind(red1,red2))
  ###write.table(x = reds, 
  ###            file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/mervl_tags.bed",
  ###            sep = "\t", quote = F, col.names = F, row.names = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "_Input.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,3,5)]
  bw02 <- list.files(path = path, pattern = "_Input-notSES.bw", full.names = T, recursive = T)
  bw02 <- bw02[grep(x = bw02, pattern = "bcp")][c(1,3,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/mervl_tags.bed"
    beds <- c(bed1)
    shell <- paste0(path,"zempROOM/run_y.sh")
    cat("#!/bin/bash\n", file = shell)
    for (inx0 in beds) {
      cmd_03 <- paste0("computeMatrix reference-point --referencePoint center -p 20 ",
                       "-a 3000 -b 3000 --missingDataAsZero -R ", inx0, " ",
                       "-S ",paste(bw01, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),"-1.mat.gz","\n")
      cmd_04 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),"-1.mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".SES.r1-center.pdf","\n")
      ###
      ###
      ###
      cmd_07 <- paste0("computeMatrix reference-point --referencePoint center -p 20 ",
                       "-a 3000 -b 3000 --missingDataAsZero -R ", inx0, " ",
                       "-S ",paste(bw02, collapse = " ")," ",
                       "-o ",path,"zempROOM/IP-VS-", basename(inx0),"-2.mat.gz","\n")
      cmd_08 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",path,
                       "zempROOM/IP-VS-",basename(inx0),"-2.mat.gz ",
                       "--samplesLabel KO-1 RN-1 WT-1 ",
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".notSES.r1-center.pdf","\n")
      #for (i in cmd_03) {cat(i, append = T, file = shell)}
      #for (i in cmd_04) {cat(i, append = T, file = shell)}
      for (i in cmd_07) {cat(i, append = T, file = shell)}
      for (i in cmd_08) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_y.log"," 2>&1 &"))
  }
}













