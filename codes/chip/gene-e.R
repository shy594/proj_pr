#write by ~~~ at 2023.07.16
#1、搭环境
{
  for (i in c("ReactomePA","ggplot2","stringr","clusterProfiler","ChIPseeker","TxDb.Mmusculus.UCSC.mm10.knownGene","topGO","Rgraphviz","pathview","org.Mm.eg.db")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/publicData/KAP1/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bowtie2/mm10"
  inx2 = "/Reference/aaSHY/BED/TEsubFamily/MERVL-int::ERVL::LTR.bed"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/BED/special/rRNA.bed"
  inxx = c(inx2,inx3)
  inx4 = "/Reference/aaSHY/zOther/Rn45s/INDEX/Rn45s"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx7 = "/Reference/aaSHY/GSEAgmt/TwoGenes.gmt"
}
#2、写命令
{
prefix <- sapply(stringr::str_split(basename(list.files(paste0(path,"fq"), pattern = "gz$",full.names = T)),".fq"),"[",1)
cmd_1 <- paste0("#trim_galore --phred33 --fastqc -o ",
                paste0(path,"trim "),
                paste0(path,"fq/",prefix,".fq.gz"),"\n")
cmd_2 <- paste0("#bowtie2 -q --very-sensitive -k 5 --reorder -p 20 -x ",inx1," -U ",
                paste0(path,"trim/",prefix,"_trimmed.fq.gz"),
                " |samtools view -F 4 -b > ",
                paste0(path,"bam/",prefix,".bam"),"\n")
cmd_3 <- paste0("#samtools index ",paste0(path,"bam/",prefix,".bam"),"\n")
cmd_4 <- paste0("#run-csem -p 20 --bam ",
                paste0(path,"bam/",prefix,".bam")," --no-extending-reads 50 ",
                paste0(path,"bam/csem/",prefix),"\n")
cmd_5 <- paste0("#bamCoverage --normalizeUsing CPM -p 20 --minMappingQuality 1 --samFlagExclude 256 --outFileFormat bigwig --binSize 20 -b ",
                paste0(path,"bam/csem/",prefix,".sorted.bam -o "),
                paste0(path,"bw/bcv/",prefix,".bcv.bw\n"))
cmd_6 <- paste0("#bamCompare -p 20 --scaleFactorsMethod SES --sampleLength 100000 -b1 ",
                paste0(path,"bam/csem/",prefix[2],".sorted.bam -b2 "),
                paste0(path,"bam/csem/",prefix[1],".sorted.bam -o "),
                paste0(path,"bw/bcp/IPvsInput.bw"),"\n")
cmd_7 <- paste0("#computeMatrix reference-point --referencePoint center -p 30 -a ",
                c(5000,5000)," -b ",c(5000,5000)," --missingDataAsZero --skipZeros ",
                "-R ", inxx, " -S ",
                paste0(path,"bw/bcp/IPvsInput.bw -o "),
                paste0(path,"deeptools/",basename(inxx),".mat.gz"),"\n")
cmd_8 <- paste0("#computeMatrix scale-regions -p 30 -a 2000 -b 2000 -m 5000 ",
                "--missingDataAsZero -p 20 ",
                "-R ", inxx, " -S ",
                paste0(path,"bw/bcp/IPvsInput.bw -o "),
                paste0(path,"deeptools/",basename(inxx),".mat.gz"),"\n")
cmd_9 <- paste0("#plotProfile --plotHeight 8 --plotWidth 8 -m ",
                paste0(path,"deeptools/",basename(inxx),".mat.gz")," -out ",
                paste0(path,"deeptools/profile/",basename(inxx),".pdf"),"\n")
cmd_a <- paste0("#plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                paste0(path,"deeptools/",basename(inxx),".mat.gz")," -out ",
                paste0(path,"deeptools/heatmap/",basename(inxx),".pdf"),"\n")
cmd_b <- paste0("macs2 callpeak -f AUTO -g mm --keep-dup auto -t ",
                paste0(path,"bam/csem/",prefix[2],".sorted.bam -c "),
                paste0(path,"bam/csem/",prefix[1],".sorted.bam -n KAP1 --outdir "),
                paste0(path,"peak/"))
shell <- paste0("/ChIP_seq_2/aaSHY/pr/publicData/KAP1/run_a.sh")
cat("#!/bin/bash\n", file = shell)
for (i in cmd_1) {cat(i, append = T, file = shell)}
for (i in cmd_2) {cat(i, append = T, file = shell)}
for (i in cmd_3) {cat(i, append = T, file = shell)}
for (i in cmd_4) {cat(i, append = T, file = shell)}
for (i in cmd_5) {cat(i, append = T, file = shell)}
for (i in cmd_6) {cat(i, append = T, file = shell)}
for (i in cmd_7) {cat(i, append = T, file = shell)}
for (i in cmd_8) {cat(i, append = T, file = shell)}
for (i in cmd_9) {cat(i, append = T, file = shell)}
for (i in cmd_a) {cat(i, append = T, file = shell)}
for (i in cmd_b) {cat(i, append = T, file = shell)}
print(paste0("nohup bash ",shell," >",path,"Log/run_a.log"," 2>&1 &"))
}
##
##
##
#3、画peak
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
gene <- rtracklayer::import(con = inx5, format = "gtf")
txdb <- GenomicFeatures::makeTxDbFromGRanges(gr = gene, taxonomyId = 10090)#GenomicFeatures::seqinfo(txdb)
#
twoc <- read.table(file = inx7, header = F, sep = "\t", stringsAsFactors = F)[,-c(1,2)]
twoc <- as.data.frame(t(twoc))
gen2 <- gene[which(gene$gene_name %in% twoc$V1)]
txdb <- GenomicFeatures::makeTxDbFromGRanges(gr = gen2, taxonomyId = 10090)#GenomicFeatures::seqinfo(txdb)
#
peak <- "/ChIP_seq_2/aaSHY/pr/publicData/KAP1/peak/KAP1_peaks.narrowPeak"
x <- annotatePeak(peak, tssRegion=c(-1000, 1000), TxDb=txdb)#as.GRanges(x) %>% head(5)
write.csv(x = as.data.frame(x), file = paste0(path,"peak/chipseeker_anno.csv"))
#
peakAnno <- ChIPseeker::annotatePeak(peak, tssRegion=c(-2000, 2000), TxDb=txdb, annoDb="org.Mm.eg.db")
temp <- peakAnno@detailGenomicAnnotation
test <- as.data.frame(peakAnno@anno[as.integer(rownames(temp[which(temp$Intron==T),]))])
test$posi <- paste0("chr",test$seqnames,":",test$start,":",test$end)
write.table(x = test, file = paste0(path,"peak/IntronPeak.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
#
pdf(file = "IPme3_merge_peaks_annoPie.pdf",width = 8,height = 8)
plotAnnoPie(peakAnno,legend.position ="rightside")
dev.off()





















