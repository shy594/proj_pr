#write by ~~~ at 2024.06.07
#01、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/zssss/chipSEQ/"
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
  lay1 <- paste0(path,"fq/",prex,"_1.fq.gz")
  lay2 <- paste0(path,"fq/",prex,"_2.fq.gz")
  lay3 <- paste0(path,"trim/",prex,"_1_val_1.fq.gz")
  lay4 <- paste0(path,"trim/",prex,"_2_val_2.fq.gz")
  cmd_03 <- paste0("trim_galore --phred33 -j 2 ",
                   "--illumina --fastqc ",#"-a 'AGATCGGAAGAGC -a G{10}' -a2 'AGATCGGAAGAGC -a G{10}' --fastqc ",
                   "--clip_R1 3 --clip_R2 3 --three_prime_clip_R1 1 --three_prime_clip_R2 1 ",
                   "-o ", path, "trim --paired ",lay1," ",lay2, "\n")
  cmd_04 <- paste0("bowtie2 -q --sensitive -k 5 --no-mixed --no-discordant ",
                   "--reorder -p 20 -x ",inx1," -1 ",lay3," -2 ",lay4,
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
                   paste0(path,"bam/csem/",prex[1],".sorted.bam -b2 "),
                   paste0(path,"bam/csem/",prex[2],".sorted.bam -o "),
                   paste0(path,"bw/bcp/",prex[1],"_Input.bw"),"\n")
  cmd_08 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                   "--missingDataAsZero -R ", c(inx7,inx8,inx9,inxA,inxB),
                   " -S ", paste(path,"bw/bcp/",prex[1],"_Input.bw", sep = "", collapse = " "),
                   " -o ", path,"deeptools/", basename(c(inx7,inx8,inx9,inxA,inxB)), ".mat.gz","\n")
  cmd_09 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "deeptools/",basename(c(inx7,inx8,inx9,inxA,inxB)),".mat.gz ",
                   "-out ",path,"PLOT/",
                   basename(c(inx7,inx8,inx9,inxA,inxB)),".pdf","\n")
  cmd_10 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --keep-dup 1 -t ",
                   path,"bam/csem/",prex[1],".sorted.bam -c ",
                   path,"bam/csem/",prex[2],".sorted.bam -n ",
                   prex[1],"-N --outdir ",
                   path,"peak/", "\n")
  cmd_11 <- paste0("macs2 callpeak -f BAM -g mm -p 0.001 --broad --broad-cutoff 0.001 ",
                   "--keep-dup 1 -t ",path,"bam/csem/", prex[1],
                   ".sorted.bam -c ", path,"bam/csem/", prex[2],
                   ".sorted.bam -n ", prex[1],
                   "-B --outdir ", path,"peak/", "\n")
  cmd_12 <- paste0("TEtranscripts -t ",
                   paste0(path,"bam/csem/",prex[1],".sorted.bam", collapse = " ")," -c ",
                   paste0(path,"bam/csem/",prex[2],".sorted.bam", collapse = " "),
                   " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project zscan4 --outdir ",paste0(path,"count/"),"\n")
  cmd_13 <- paste0("TElocal -b ",
                   path,"bam/csem/",prex,".sorted.bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_10) {cat(i, append = T, file = shell)}
  for (i in cmd_11) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  for (i in cmd_13) {cat(i, append = T, file = shell)}
  for (kk in c("mervl","rRNA")) {
    if (kk == "mervl") {
      inxM <- inxD
    }
    else {
      inxM <- inxC
    }
    cmd_14 <- paste0("bowtie2 -q --sensitive -k 5 --no-mixed --no-discordant ",
                     "--reorder -p 20 -x ",inxM," -1 ", lay3," -2 ",lay4," ",
                     "|samtools view -F 4 -b > ", path, kk, "/",prex,".bam","\n")
    cmd_15 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                     paste0(path,kk,"/",prex,".bam")," 140 ",
                     paste0(path,kk,"/csem/",prex),"\n")
    cmd_16 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     paste0(path,kk,"/csem/",prex,".sorted.bam -o "),
                     paste0(path,kk,"/",prex,".",kk,".bw\n"))
    for (i in cmd_14) {cat(i, append = T, file = shell)}
    for (i in cmd_15) {cat(i, append = T, file = shell)}
    for (i in cmd_16) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
################################################################################
#04、peak Analysis (pvalue 0.05)
{
  prex <- c("IP","Input")
  gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
  gtf1 <- as.data.frame(gtf1)
  colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
  gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
  txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
  twoC <- read.table(inxB)
  gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
  txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
  #
  for (i in (prex[1])) {
    pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",i,"-N_peaks.narrowPeak"), 
                                                      tssRegion = c(-2000, 2000), TxDb = txb1))
    pdf(file = paste0(path, "PLOT/",i,"-N_chipseekerAnno.pdf"), width = 8, height = 8)
    ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
    dev.off()
    write.csv(as.data.frame(pk01), paste0(path, "peak/",i,"-N_chipseekerAnno.csv"), quote = F, row.names = F)
  }
}
################################################################################
#05、TE Binding situation [适用于ChIP-seq]
{
  #Observed
  red1 <- read.table(file = paste0(path,"count/zscan4.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- prex
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
    openxlsx::write.xlsx(x = red3, asTable = T, file = paste0(path,"fisher/",prex[1:2][m],"FisherTE.xlsx"))
  }
}
################################################################################
#06、Gene Binding situation [适用于ChIP-seq]
{
  twoc <- read.table(file = inxB, sep = "\t")
  #Observed
  red1 <- read.table(file = paste0(path,"count/zscan4.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
  colnames(red1) <- prex
  red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
  #Expected
  mmu1 <- read.table(file = paste0(inx3,".fai"))
  refL <- sum(mmu1$V2)
  mmu2 <- rtracklayer::import(con = inx4,format = "gtf")
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
    openxlsx::write.xlsx(x = red3, asTable = T, file = paste0(path,"fisher/",prex[1:2][m],"FisherGene.xlsx"))
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#07、取circRNA侧翼, 看看ZSCAN4是否有富集, 直接延伸 [早期胚胎发育的数据]
{
  devs <- c("oocyte","zygote","two","four","eight","morula","blasto")
  for (i in devs){
    ff01 <- read.table(file = paste0("/ChIP_seq_1/aaSHY/ngsSUPeR/aaMerge/theCIRC/find_circ/",i,".circCandidate.txt"))
    ff01$mark <- "."
    wr01 <- ff01[,c("V1","V2","V3","mark","V6","mark")]
    write.table(x = wr01, 
                file = paste0("/ChIP_seq_1/aaSHY/ngsSUPeR/aaMerge/theCIRC/find_circ/",i,".circCandidate.bed"),
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
  bed1 <- paste0("/ChIP_seq_1/aaSHY/ngsSUPeR/aaMerge/theCIRC/find_circ/",devs,".circCandidate.bed", collapse = " ")
  bw01 <- "/ChIP_seq_2/aaSHY/Zscan4/chipSEQ/bw/bcp/IP_Input.bw"
  cmd_08 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                   "--missingDataAsZero -p 20 ", "-R ", bed1, " -S ",bw01," -o ", path,
                   "dumpROOM/", "circ.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/circ.mat.gz --samplesLabel ZSCAN4-circRNA ",
                   "--regionsLabel ",
                   paste0(devs, collapse = " "), " -out ",
                   path,"dumpROOM/circ.zscan4.v3.pdf","\n")
}
################################################################################
#08、取circRNA侧翼, 看看ZSCAN4是否有富集, 直接延伸 [Ssrp1的数据]
{
  path = "/ChIP_seq_1/aaSHY/circrna/Ssrp1/aaReRun/"
  devs <- c("ctrlRep1","ctrlRep2","KORep1","KORep2")
  for (i in devs){
    ff01 <- read.table(file = paste0("/ChIP_seq_1/aaSHY/circrna/Ssrp1/aaReRun/circRNA/findcirc/",i,".circCandidate.txt"))
    ff01$mark <- "."
    wr01 <- ff01[,c("V1","V2","V3","mark","V6","mark")]
    write.table(x = wr01, 
                file = paste0("/ChIP_seq_1/aaSHY/circrna/Ssrp1/aaReRun/circRNA/findcirc/",i,".circCandidate.bed"),
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
  bed1 <- paste0("/ChIP_seq_1/aaSHY/circrna/Ssrp1/aaReRun/circRNA/findcirc/",devs,".circCandidate.bed", collapse = " ")
  bw01 <- "/ChIP_seq_2/aaSHY/Zscan4/chipSEQ/bw/bcp/IP_Input.bw"
  cmd_08 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                   "--missingDataAsZero -p 20 ", "-R ", bed1, " -S ",bw01," -o ", path,
                   "dumpROOM/", "circ.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/circ.mat.gz --samplesLabel ZSCAN4-circRNA ",
                   "--regionsLabel ",
                   paste0(devs, collapse = " "), " -out ",
                   path,"dumpROOM/circ.zscan4.ssrp1CIRC.v1.pdf","\n")
}
################################################################################
#09、取circRNA侧翼, 看看ZSCAN4是否有富集, 直接延伸 [GSE157788的数据 #RNase处理]
{
  path = "/ChIP_seq_1/aaSHY/circrna/circZNF827/"
  devs <- c("mESC1R","mESC1","mESC2R","mESC2")
  ff01 <- read.table(file = "/ChIP_seq_1/aaSHY/circrna/circZNF827/suppl/GSE157788_find_circ.txt", header = T)[,1:10]
  ff01$chrom <- gsub(x = ff01$chrom, pattern = "chr", fixed = T, replacement = "")
  colnames(ff01) <- c("chrom","start","end","name","score","strand",devs)
  for (i in devs) {
    wr01 <- ff01[,c("chrom","start","end","name","score","strand",i)]
    wr01 <- wr01[which(wr01[,i] != 0),1:6]
    write.table(x = wr01, 
                file = paste0("/ChIP_seq_1/aaSHY/circrna/circZNF827/dumpROOM/",i,".circCandidate.bed"),
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
  bed1 <- paste0("/ChIP_seq_1/aaSHY/circrna/circZNF827/dumpROOM/",devs,".circCandidate.bed", collapse = " ")
  bw01 <- "/ChIP_seq_2/aaSHY/Zscan4/chipSEQ/bw/bcp/IP_Input.bw"
  cmd_08 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                   "--missingDataAsZero -p 20 ", "-R ", bed1, " -S ",bw01," -o ", path,
                   "dumpROOM/", "circ.zscan4.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/circ.zscan4.mat.gz --samplesLabel ZSCAN4-circRNA ",
                   "--regionsLabel ",
                   paste0(devs, collapse = " "), " -out ",
                   path,"dumpROOM/circ.zscan4.gse157788CIRC.v1.pdf","\n")
}
################################################################################
#10、取所有Alu的区域，看看ZSCAN4是否有富集
{
  shell <- paste0(path,"src/run_i.sh")
  cat("#!/bin/bash\n", file = shell)
  name <- c("MERVL-int::ERVL::LTR")
  bed1 <- list.files("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily", pattern = "Alu", full.names = T)
  bed2 <- paste0("/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/",name,".bed")
  bw01 <- "/ChIP_seq_2/aaSHY/Zscan4/chipSEQ/bw/bcp/IP_Input.bw"
  cmd_01 <- paste0("computeMatrix reference-point --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -p 10 -R ", c(bed1,bed2),
                   " -S ",bw01,
                   " -o ",path,"dumpROOM/Alu/TE.zscan4.",basename(c(bed1,bed2)),
                   ".mat.gz","\n")
  cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/Alu/TE.zscan4.",basename(c(bed1,bed2)),".mat.gz ",
                   "--samplesLabel ZSCAN4-TE ",
                   "--regionsLabel ",
                   gsub(x = basename(c(bed1,bed2)), pattern="::SINE.bed|::LTR.bed", replacement = ""),
                   " -out ",path,"dumpROOM/Alu/TE.zscan4.",basename(c(bed1,bed2)),".pdf","\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"Log/run_i.log"," 2>&1 &"))
}
