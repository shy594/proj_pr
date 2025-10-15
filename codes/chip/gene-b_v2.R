#write by ~~~ at 2024.03.15
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/ff/"
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
  cmd_01 <- paste0("zcat ",ff01," ",ff02, " |gzip >",path,"fq/",prex,".fq.gz","\n")
  for (i in cmd_01) {cat(i, file = shell, append = T)}
  cmd_02 <- paste0("trim_galore --phred33 --fastqc --illumina ",path,
                   "fq/", prex, ".fq.gz -o ", path, "trim", "\n")
  cmd_03 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inx1," -U ", path, "trim/", prex, "_trimmed.fq.gz",
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_04 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"bam/",prex,".bam")," 140 ",
                   paste0(path,"bam/csem/",prex),"\n")
  cmd_05 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
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
      cmd_07 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                       paste0(path,"bam/csem/",zrex[2],".sorted.bam -b2 "),
                       paste0(path,"bam/csem/",zrex[1],".sorted.bam -o "),
                       paste0(path,"bw/bcp/",zrex[2],"_Input.bw"),"\n")
      cmd_08 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 5000 ",
                       "--missingDataAsZero ", "-R ", c(inx7,inx8), " -S ", path,
                       "bw/bcp/",q,"_Input.bw -o ", path,"deeptools/", 
                       c(basename(inx7),basename(inx8)),q, ".mat.gz","\n")
      cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 500 ",
                       "--missingDataAsZero -R ", inx9, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inx9),q,".mat.gz","\n")
      cmd_10 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                       "--missingDataAsZero -R ", inxA, " -S ",path,
                       "bw/bcp/",q,"_Input.bw -o ",path,"deeptools/",
                       basename(inxA),q,".mat.gz","\n")
      cmd_11 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 30000 ",
                       "--missingDataAsZero -R ", inxB, " -S ",path,
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
#05、gene费舍尔检验--适用于ChIP-seq
{
  bath = paste0(path, "single/")
  prex <- unique(sapply(str_split(list.files(path = paste0(bath, "fq")), pattern = "-R|.fq"), "[", 1))
  for (x in seq(2)) {
    zrex <- prex[c(x,x+2)]; rm(x)
    twoc <- read.table(file = inxB, sep = "\t")
    #Observed
    red1 <- read.table(file = paste0(bath,"count/ff_",zrex[2],".cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
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
#06、tete费舍尔检验--适用于ChIP-seq
{
  bath = paste0(path, "single/")
  prex <- unique(sapply(str_split(list.files(path = paste0(bath, "fq")), pattern = "-R|.fq"), "[", 1))
  for (x in seq(2)) {
    zrex <- prex[c(x,x+2)]; rm(x)
    #Observed
    red1 <- read.table(file = paste0(bath,"count/ff_",zrex[2],".cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
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
#07、peak所在基因: ff & pr
{
  ff01 <- read.csv(paste0(path,"single/peak/IP-1-allg-N_chipseekerAnno.csv"))
  ff02 <- read.csv(paste0(path,"single/peak/IP-2-allg-N_chipseekerAnno.csv"))
  bath <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  pf01 <- read.csv(paste0(bath,"single/peak/IP-1-allg-N_chipseekerAnno.csv"))
  pf02 <- read.csv(paste0(bath,"single/peak/IP-2-allg-N_chipseekerAnno.csv"))
  tmp1 <- list()
  tmp1[["a"]] <- unique(ff01$geneId)
  tmp1[["b"]] <- unique(ff02$geneId)
  tmp1[["c"]] <- unique(pf01$geneId)
  tmp1[["d"]] <- unique(pf02$geneId)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "ff-IP-1,ff-IP-2,pr-IP-1,pr-IP-2", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("ff-IP-1","ff-IP-2","pr-IP-1","pr-IP-2"),
                                    fill=c("#FFFFCC","#CCFFFF","#FFCCCC","#CCFFCC"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT-ff/geneVenn.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  intersect(intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]]),tmp1[["d"]])
}
################################################################################
#08、peak所在基因: ff & pr [取各自样本1; distal intergenic: -10000; downstream: -10000]
{
  ff01 <- read.csv(paste0(path,"single/peak/IP-1-allg-N_chipseekerAnno.csv"))
  tst1 <- ff01[grep(ff01$annotation, pattern="Distal|Down"),]
  tst1 <- tst1[which(tst1$distanceToTSS > -10000 & tst1$distanceToTSS < 1000),]
  tst2 <- ff01[grep(ff01$annotation, pattern="Distal|Down", invert = T),]
  ffa <- rbind(tst1,tst2)
  #
  bath <- "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/pr/"
  pf01 <- read.csv(paste0(bath,"single/peak/IP-2-allg-N_chipseekerAnno.csv"))
  tst1 <- pf01[grep(pf01$annotation, pattern="Distal|Down"),]
  tst1 <- tst1[which(tst1$distanceToTSS > -10000 & tst1$distanceToTSS < 1000),]
  tst2 <- pf01[grep(pf01$annotation, pattern="Distal|Down", invert = T),]
  prmt <- rbind(tst1,tst2)
  #
  tmp1 <- list()
  tmp1[["a"]] <- unique(ffa$geneId)
  tmp1[["b"]] <- unique(prmt$geneId)
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "ff-IP-1,pr-IP-1", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("ff-IP-1","pr-IP-1"),
                                    fill=c("#FFFFCC","#CCFFFF"),
                                    units = "cm", height = 12, width = 12)
  pdf(file = paste0(path, "PLOT-ff/geneVenn.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  intersect(intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]]),tmp1[["d"]])
}
tst3 <- sapply(str_split(prmt$annotation, pattern = "\\("),"[",1)
table(tst3)
################################################################################
#09、取所有Alu的区域，看看ff是否有富集
{
  shell <- paste0(path,"single/src/run_e.sh")
  cat("#!/bin/bash\n", file = shell)
  name <- c("MERVL-int::ERVL::LTR")
  bed1 <- list.files("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily", pattern = "Alu", full.names = T)
  bed2 <- paste0("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/",name,".bed")
  bw01 <- paste("/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/ff/single//bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pr/ChIPseq/20240315/ff/single//bw/bcp/IP-2_Input.bw", collapse = " ")
  cmd_01 <- paste0("computeMatrix reference-point --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -p 10 -R ", c(bed1,bed2),
                   " -S ",bw01,
                   " -o ",path,"single/dumpROOM/enrichAlu/TE.ff.",basename(c(bed1,bed2)),
                   ".mat.gz","\n")
  cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "single/dumpROOM/enrichAlu/TE.ff.",basename(c(bed1,bed2)),".mat.gz ",
                   "--samplesLabel ff-1 ff-2 ",
                   "--regionsLabel ",
                   gsub(x = basename(c(bed1,bed2)), pattern="::SINE.bed|::LTR.bed", replacement = ""),
                   " -out ",path,"single/dumpROOM/enrichAlu/TE.ff.",basename(c(bed1,bed2)),".pdf","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"single/Log/run_e.log"," 2>&1 &"))
}





