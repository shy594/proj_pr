#Write by Sun Haiayng on 2023.10.24
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","ggforce","bedtoolsr","refGenome","yyplot","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/RADICLfa2/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa.fai"
  inx5 = "/ChIP_seq_1/aaSHY/rnadna/RADICL/GenBank/rRNA.fa"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inx9 = "/Reference/aaSHY/BED/special/spliceF.txt"
}
##
##
##
##
##
##
#04、小RNA
te01 <- rtracklayer::import(con = inx3, format = "gtf")
ge01 <- rtracklayer::import(con = inx2, format = "gtf")
txb1 <- suppressWarnings(GenomicFeatures::makeTxDbFromGRanges(gr = ge01, taxonomyId = 10090))
ge01 <- ge01[which(ge01$type=="gene")]
ge02 <- ge01[which(ge01$gene_biotype %in% c("snRNA", "snoRNA"))]
twoc <- read.table(file = inx7, sep = "\t")
#
for (i in c("n1","n2","n3")) {
  rd01 <- fread(file = paste0(path,"bed/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
  rd01 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V2", start.field = "V3", end.field = "V4", strand.field = "V6", keep.extra.columns = T)
  ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = ge02, subject = rd01, ignore.strand = T))
  rs01 <- as.data.frame(rd01[unique(subjectHits(ov01))])
  colnames(rs01) <- paste0("V", seq_along(colnames(rs01)))
  write.table(x = rs01, file = paste0(path, "snoRNA/",i,".snRNA-DNA.bed"), sep = "\t",
              quote = F, col.names = F, row.names = F)
  wr01 <- data.frame(chr = rs01$V8, start = rs01$V9, end = rs01$V10, name = ".", score = 1, strand = rs01$V12)
  write.table(x = wr01, file = paste0(path, "snoRNA/",i,".onlyDNA.bed"), sep = "\t",
              quote = F, col.names = F, row.names = F)
  #peak analysis
  pk01 <- suppressWarnings(ChIPseeker::annotatePeak(peak = paste0(path, "snoRNA/",i,".onlyDNA.bed"), 
                                                    tssRegion = c(-2000, 2000), TxDb = txb1))
  pdf(file = paste0(path, "snoRNA/",i,".chipseekerAnno.pdf"), width = 8, height = 8)
  ChIPseeker::plotAnnoPie(pk01, 
                          main = paste0("\n\nDistribution of Regions \nBound by snRNA/snoRNA (sample ",i,")"), 
                          legend.position = "rightside")
  dev.off()
  #sheet by gene
  wr02 <- GenomicRanges::makeGRangesFromDataFrame(df = wr01, keep.extra.columns = T)
  ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = ge01, subject = wr02, ignore.strand = T))
  ge03 <- as.data.frame(ge01[queryHits(ov02)])
  ge04 <- ge03[,c("gene_name", "width")]
  ge05 <- as.data.frame(dplyr::group_by(.data = ge04, gene_name, width) %>% dplyr::summarise(count = table(gene_name)))
  tem1 <- as.data.frame(ge01)
  tem2 <- tem1[which(tem1$gene_name %in% unique(setdiff(tem1$gene_name, ge05$gene_name))),c("gene_name","width")]
  tem2$count <- 0
  ge05 <- rbind(ge05, tem2)
  ge05$cpm <- sapply(ge05$count, function(x){x/sum(ge05$count)*1000000})
  ge05$tpm <- NA
  refL <- sum(ge05$count/ge05$width)
  for (x in seq_along(ge05$gene_name)) {
    ge05$tpm[x] <- (ge05$count[x]*1000000/ge05$width[x])/refL
  }
  ge05$marker <- NA
  ge05$marker[which(ge05$gene_name %in% twoc$V4)] <- "YES"
  ge05$marker[which(!(ge05$gene_name %in% twoc$V4))] <- "NO"
  ge05$self <- NA
  ge05$self[which(ge05$gene_name %in% ge02$gene_name)] <- "YES"
  ge05$self[which(!(ge05$gene_name %in% ge02$gene_name))] <- "NO"
  ge05 <- ge05[order(ge05$tpm, decreasing = T),]
  openxlsx::write.xlsx(x = ge05, file = paste0(path, "snoRNA/",i,".countByGene.xlsx"), asTable = T)
  #ucsc track
  wr03 <- as.data.frame(wr02[unique(subjectHits(ov02))])
  wr03 <- wr03[grep(x=wr03$seqnames, pattern="GL|JH", invert=T),]
  wr03$seqnames <- gsub(x = wr03$seqnames, pattern = "MT", replacement = "M")
  wr03$seqnames <- paste0("chr", wr03$seqnames)
  wr03 <- wr03[,c(1,2,3,6,7,5)]
  write.table(x = wr03, file = paste0(path, "snoRNA/",i,".ucscTrack.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
}
##
##
##
##
##
##
