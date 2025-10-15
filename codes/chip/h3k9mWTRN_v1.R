#write by ~~~ at 2023.11.14
#1、搭环境
{
  rm(list=ls())
  for (i in c("ReactomePA","ggplot2","stringr","ChIPseeker","DESeq2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db","TxDb.Mmusculus.UCSC.mm10.knownGene")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ChIPseq/20231120/H3K9me3/mWTRN/"
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
#3、跑命令
{
  for (i in c("src","fq","Log","peak","dumpROOM","DESeq2","rRNA/csem","mervl/csem","fisher","deeptools","bam/csem","trim","PLOT","bw/bcv","bw/bcp","count")) {
    if (dir.exists(paths = paste0(path, i))==F) {
      dir.create(path = paste0(path,i), recursive = T)
    }
  }
  shell <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  prex <- c("mWTRN-IP-1","mWTRN-IP-2","mWTRN-Input-1","mWTRN-Input-2")#; unique(sapply(str_split(list.files(path = paste0(path, "rawdata")), pattern = "74-|_"), "[", 2))
  lay1 <- list.files(path = paste0(path,"rawdata"), recursive = T, pattern = "_1.fq.gz", full.names = T)
  lay2 <- list.files(path = paste0(path,"rawdata"), recursive = T, pattern = "_2.fq.gz", full.names = T)
  lay3 <- paste0(path,"dumpROOM/",prex,"_1.fq.gz")
  lay4 <- paste0(path,"dumpROOM/",prex,"_2.fq.gz")
  cmd_01 <- paste0("zcat ", c(lay1, lay2), " |sed 's/^@/@",c(rep(1,4),rep(2,4)),"/g' |gzip >", c(lay3,lay4), "\n")
  cmd_02 <- paste0("zcat ",lay3, " ", lay4," |gzip >",paste0(path,"fq/",prex,".fq.gz"),"\n")
  cmd_03 <- paste0("fastqc -t 10 -o ", path, "fq/ ", path,"fq/",prex,".fq.gz", "\n")
  cmd_04 <- paste0("trim_galore --phred33 --fastqc --illumina ",
                   "--clip_R1 10 ",path, 
                   "fq/", prex, ".fq.gz -o ", path, "trim", "\n")
  cmd_05 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                   inx1," -U ", path, "trim/", prex, "_trimmed.fq.gz",
                   " |samtools view -F 4 -b > ", path,
                   "bam/",prex,".bam","\n")
  cmd_06 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                   paste0(path,"bam/",prex,".bam")," 140 ",
                   paste0(path,"bam/csem/",prex),"\n")
  cmd_07 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                   "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                   paste0(path,"bam/csem/",prex,".sorted.bam -o "),
                   paste0(path,"bw/bcv/",prex,".bw\n"))
  cmd_08 <- paste0("bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                   paste0(path,"bam/csem/",prex[1:2],".sorted.bam -b2 "),
                   paste0(path,"bam/csem/",prex[3:4],".sorted.bam -o "),
                   paste0(path,"bw/bcp/",prex[1:2],"VSinput.bw"),"\n")
  for (i in cmd_01) {cat(i, append = T, file = shell)}
  for (i in cmd_02) {cat(i, append = T, file = shell)}
  for (i in cmd_03) {cat(i, append = T, file = shell)}
  for (i in cmd_04) {cat(i, append = T, file = shell)}
  for (i in cmd_05) {cat(i, append = T, file = shell)}
  for (i in cmd_06) {cat(i, append = T, file = shell)}
  for (i in cmd_07) {cat(i, append = T, file = shell)}
  for (i in cmd_08) {cat(i, append = T, file = shell)}
  for (q in prex[1:2]) {
    cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                     "--missingDataAsZero -p 20 ", "-R ", c(inx7,inx8), " -S ", path,
                     "bw/bcp/",q,"VSinput.bw -o ", path,"deeptools/", 
                     c(basename(inx7),basename(inx8)),q, ".mat.gz","\n")
    cmd_10 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 500 ",
                     "--missingDataAsZero -p 20 -R ", inx9, " -S ",path,
                     "bw/bcp/",q,"VSinput.bw -o ",path,"deeptools/",
                     basename(inx9),q,".mat.gz","\n")
    cmd_11 <- paste0("computeMatrix scale-regions -p 20 -a 2000 -b 2000 -m 5000 ",
                     "--missingDataAsZero -p 20 -R ", inxA, " -S ",path,
                     "bw/bcp/",q,"VSinput.bw -o ",path,"deeptools/",
                     basename(inxA),q,".mat.gz","\n")
    cmd_12 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 -m 30000 ",
                     "--missingDataAsZero -p 20 -R ", inxB, " -S ",path,
                     "bw/bcp/",q,"VSinput.bw -o ",path,"deeptools/",
                     basename(inxB),q,".mat.gz","\n")
    cmd_13 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "deeptools/",c(basename(inx7),basename(inx8),basename(inx9),
                                    basename(inxA),basename(inxB)),q,
                     ".mat.gz -out ",path,"PLOT/",
                     c(basename(inx7),basename(inx8),
                       basename(inx9),basename(inxA),basename(inxB)),q,".pdf","\n")
    for (i in cmd_09) {cat(i, append = T, file = shell)}
    for (i in cmd_10) {cat(i, append = T, file = shell)}
    for (i in cmd_11) {cat(i, append = T, file = shell)}
    for (i in cmd_12) {cat(i, append = T, file = shell)}
    for (i in cmd_13) {cat(i, append = T, file = shell)}
  }
  cmd_14 <- paste0("macs2 callpeak -f BAM -g mm -q 0.05 --keep-dup 1 -t ",
                   path,"bam/csem/",prex[1:2],".sorted.bam -c ",
                   path,"bam/csem/",prex[3:4],".sorted.bam -n ",
                   prex[1:2],"-N --outdir ",
                   path,"peak/", "\n")
  cmd_15 <- paste0("macs2 callpeak -f BAM -g mm -q 0.05 --broad --broad-cutoff 0.05 ",
                   "--keep-dup 1 -t ",path,"bam/csem/", prex[1:2],
                   ".sorted.bam -c ", path,"bam/csem/", prex[3:4],
                   ".sorted.bam -n ", prex[1:2],
                   "-B --outdir ", path,"peak/", "\n")
  cmd_16 <- paste0("TEtranscripts -t ",
                   paste0(path,"bam/csem/",prex[1:2],".sorted.bam", collapse = " ")," -c ",
                   paste0(path,"bam/csem/",prex[3:4],".sorted.bam", collapse = " "),
                   " --GTF ",inx4," --TE ",inx5," --sortByPos --mode multi ",
                   "--project ccPrmt --outdir ",paste0(path,"count/"),"\n")
  cmd_17 <- paste0("TElocal -b ",
                   path,"bam/csem/",prex,".sorted.bam ",
                   "--GTF ",inx4," --TE ",inx6," --sortByPos --project ",
                   path, "count/", prex,"\n")
  for (i in cmd_14) {cat(i, append = T, file = shell)}
  for (i in cmd_15) {cat(i, append = T, file = shell)}
  for (i in cmd_16) {cat(i, append = T, file = shell)}
  for (i in cmd_17) {cat(i, append = T, file = shell)}
  for (kk in c("mervl","rRNA")) {
    if (kk == "mervl") {
      inxM <- inxD
    }
    else {
      inxM <- inxC
    }
    cmd_18 <- paste0("bowtie2 -q --sensitive -k 5 --reorder -p 20 -x ",
                     inxM," -U ", path, "trim/", prex, "_trimmed.fq.gz |",
                     "samtools view -F 4 -b > ", path, kk, "/",prex,".bam","\n")
    cmd_19 <- paste0("run-csem -p 20 --no-extending-reads --bam ",
                     paste0(path,kk,"/",prex,".bam")," 140 ",
                     paste0(path,kk,"/csem/",prex),"\n")
    cmd_20 <- paste0("bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 ",
                     "--samFlagExclude 256 -of bigwig -bs 20 -b ",
                     paste0(path,kk,"/csem/",prex,".sorted.bam -o "),
                     paste0(path,kk,"/",prex,".",kk,".bw\n"))
    for (i in cmd_18) {cat(i, append = T, file = shell)}
    for (i in cmd_19) {cat(i, append = T, file = shell)}
    for (i in cmd_20) {cat(i, append = T, file = shell)}
  }
  print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
##
##
##
#4、peak Analysis (pvalue 0.05)
gtf1 <- rtracklayer::import(con = inx4, format = "gtf")
gtf1 <- as.data.frame(gtf1)
colnames(gtf1)[c(10,12)] <- c("gene_name","gene_id")
gtf1 <- GenomicRanges::makeGRangesFromDataFrame(df = gtf1, keep.extra.columns = T)
txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf1, taxonomyId = 10090))
twoC <- read.table(inxB)
gtf2 <- gtf1[which(gtf1$gene_id %in% unique(twoC$V4))]
txb2 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = gtf2, taxonomyId = 10090))
#
for (i in (prex[1:2])) {
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "peak/",i,"-N_peaks.narrowPeak"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "PLOT/",i,"-N_chipseekerAnno.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, legend.position = "rightside")
  dev.off()
  write.csv(as.data.frame(pk01), paste0(path, "peak/",i,"-N_chipseekerAnno.csv"), quote = F, row.names = F)
}
##
##
##
#5、TE Binding situation [适用于ChIP-seq]
#Observed
red1 <- read.table(file = paste0(path,"count/ccPrmt.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
colnames(red1) <- prex
red2 <- red1[grep(rownames(red1), pattern = ":", invert = F),]
#Expected
mmu1 <- read.table(file = paste0(inx3,".fai"))
refL <- sum(mmu1$V2)
mmu2 <- read.table(file = inx5, sep = "\t", stringsAsFactors = F, header = F)
mmu3 <- mmu2
mmu3[,9] <- gsub("; transcript_id.*family_id |; class_id ",":",gsub("gene_id |;$","",mmu3[,9]))
mmu3[,8] <- mmu3[,5]-mmu3[,4]+1
for (m in seq(2)) {
  print(m)
  red3 <- red2[,c(m,m+2)]
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
##
##
##
#6、基因的Binding situation [适用于ChIP-seq]
twoc <- read.table(file = inxB, sep = "\t")
#Observed
red1 <- read.table(file = paste0(path,"count/ccPrmt.cntTable"), sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
colnames(red1) <- prex
red2 <- red1[grep(rownames(red1), pattern = ":", invert = T),]
#Expected
mmu1 <- read.table(file = paste0(inx3,".fai"))
refL <- sum(mmu1$V2)
mmu2 <- rtracklayer::import(con = inx4,format = "gtf")
mmu3 <- as.data.frame(mmu2[which(mmu2$type=="exon")])
for (m in seq(2)) {
  print(m)
  red3 <- red2[,c(m,m+2)]
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
##
##
##
#7、基因的差异分析
data <- read.table(paste0(path, "count/ccPrmt.cntTable"), header = T, check.names = F, stringsAsFactors = F)
gene <- data[grep(data[,1], pattern = ":", invert = T),]
tete <- data[grep(data[,1], pattern = ":", invert = F),]
rc_a <- gene
rc_b <- tete
data <- read.table(paste0(path, "count/TEsites.cntTable"), header = T, check.names = F, stringsAsFactors = F)
rc_c <- data
shy_diff <- function(mydata, id, Class){
  counts <- mydata
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  namess <- sapply(str_split(basename(colnames(counts)),pattern = ".sorted"), "[", 1)
  colnames(counts) <- namess
  yyyCOL <- data.frame(row.names = colnames(counts), Class, stringsAsFactors = T)
  ddsdds <- DESeqDataSetFromMatrix(countData = counts, colData = yyyCOL, design = ~Class)
  ddsdds <- ddsdds[rowMeans(counts(ddsdds)) > 5,]
  ddsdds <- DESeq(ddsdds)
  resres <- results(ddsdds, contrast = c("Class","IP","WT"))
  write.csv(resres, paste0(path,"DESeq2/",id,".csv"))
}
shy_diff(mydata = rc_a[,1:5], id = "FblGene", Class = factor(c("IP","IP","WT","WT")))
shy_diff(mydata = rc_b[,1:5], id = "FblTE", Class = factor(c("IP","IP","WT","WT")))
shy_diff(mydata = rc_c[,1:5], id = "FblTEsites", Class = factor(c("IP","IP","WT","WT")))
##
##
##
################################################################################
#MERVL RNA bound gene：ChIP-seq signal trend
twoc <- read.table(file = inxB, stringsAsFactors = F, sep = "\t")
twoc <- GenomicRanges::makeGRangesFromDataFrame(df = twoc, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3", strand.field = "V6")
df01 <- data.frame()
for (i in c("n1","n2","n3")) {
  merv <- read.table(file = paste0("/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/mervlSELF/",i,".rnaDNA.bed"), sep = "\t", header = T, stringsAsFactors = F)
  merv <- GenomicRanges::makeGRangesFromDataFrame(df = merv, keep.extra.columns = T, seqnames.field = "dnaChr", start.field = "dnaSta", end.field = "dnaEnd", strand.field = "dnaStrand")
  ov01 <- GenomicRanges::findOverlaps(query = twoc, subject = merv, ignore.strand = T)
  temp <- as.data.frame(table(twoc[queryHits(ov01)]$V4))
  temp$samp <- i
  df01 <- rbind(df01, temp)
}
df02 <- tidyr::spread(data = df01, key = "samp", value = "Freq")
df02[is.na(df02)] <- 0
df02 <- df02[order(as.character(df02$Var1), decreasing = F),]
#
expr <- read.table(paste0(path, "count/ccPrmt.cntTable"), header = T, check.names = F, stringsAsFactors = F, row.names = 1)
expr <- expr[grep(expr[,1], pattern = ":", invert = T),]
colnames(expr) <- sapply(str_split(basename(colnames(expr)), pattern = ".sor"), "[", 1)
cpm1 <- as.data.frame(apply(expr, 2, function(x){x/sum(x)*1000000}))
colnames(cpm1) <- paste0("cpm", colnames(cpm1))
#所有2C基因的chip signal
cpm2 <- cpm1[which(rownames(cpm1) %in% as.character(twoc$V4)),]
cpm2 <- cpm2[order(rownames(cpm2), decreasing = F),]
cpm2$geneName <- rownames(cpm2)
wr02 <- cpm2[,c(5,7,6,8)-4]
#与MERVL RNA互作的2C基因
expr <- expr[which(rownames(expr) %in% as.character(df02$Var1)),]
expr <- expr[order(rownames(expr), decreasing = F),]
cpm1 <- cpm1[which(rownames(cpm1) %in% as.character(df02$Var1)),]
cpm1 <- cpm1[order(rownames(cpm1), decreasing = F),]
df03 <- cbind(df02, expr, cpm1)
df03 <- df03[,c(1:4,5,7,6,8,c(5,7,6,8)+4)]
openxlsx::write.xlsx(x = df03, file = paste0(path, "mervl/mervlTwoLikeSignal.xlsx"), asTable = T)
























