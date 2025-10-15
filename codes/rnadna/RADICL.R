#Write by Sun Haiayng on 2023.09.04
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
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inx5 = "/ChIP_seq_1/aaSHY/rnadna/RADICL/GenBank/rRNA.fa"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT2_Mm::ERVL::LTR.bed"
  inx9 = "/Reference/aaSHY/BED/special/spliceF.txt"
}
##
##
##
#03、跑命令
{
  shel <- paste0(path,"src/run_a.sh")
  cat("#!/bin/bash\n", file = shel)
  zemp <- list.files(path = paste0(path, "rawdata"), full.names = T)
  prex <- list.files(path = paste0(path, "rawdata"))
  cmd_01 <- paste0("parallel-fastq-dump -s ",zemp," ",#--gzip 
                   "-t 10 --split-3 -O ", path,"fq/","\n")
  cmd_02 <- paste0("rename s/fastq/fq/ ", path,"fq/*.fastq","\n")
  pair <- paste0(path,"zemp/pair.txt"); cat("", file = pair)
  for (m in paste0(path,"fq/",prex,".fq\n")) {cat(m, file =pair, append = T)}
  cmd_03 <- paste0("cat ",pair," |head -n 2 >",path,"zemp/ff01.txt","\n")
  cmd_04 <- paste0("cat ",pair," |head -n 4 |tail -n 2 >",path,"zemp/ff02.txt","\n")
  cmd_05 <- paste0("cat ",pair," |tail -n 2 >",path,"zemp/ff03.txt","\n")
  cmd_06 <- paste0("fastuniq -i ", path, "zemp/", c("ff01", "ff02", "ff03"), 
                   ".txt", " -t q -o ", path, "fq/n",c(1,2,3),"DNA.s1.fq", " -p ", path,
                   "fq/n",c(1,2,3), "RNA.s1.fq", "\n")
  cmd_07 <- paste0("#rRNAdust ", inx5, " ", path, "fq/", prex, ".s1.fq", " -o ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_08 <- paste0("#bwa aln -t 10 -f ",path,"bam/", prex, ".s2.sai ", inx4, " ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_09 <- paste0("#bwa samse ", inx4, " ", path, "bam/", prex, ".s2.sai ", path, "fq/", prex, ".s2.fq ",
                   "|samtools view -bS -F 4 -@ 10 |samtools sort -@ 10 >", path, "bam/", prex, ".s2.bam","\n")
  cmd_10 <- paste0("#samtools index ",path, "bam/", prex, ".s2.bam", "\n")
  cmd_11 <- paste0("#samtools view -@ 10 -q 0 -bS ", path, "bam/", prex, 
                   ".s2.bam ", "|bedtools bamtobed -i - >", path, "bed/", prex, ".s2.bed", "\n")
  cmd_12 <- paste0("#cat ", path, "bed/", prex, ".s2.bed |sed 's/", prex, "/n1BNA/g' |",
                   "sort -k 4b,4 >", path, "bed/", prex, ".sort.bed", "\n")
  cmd_13 <- paste0("#join -t $'\\t' -j 4 ", path, "bed/", prex[c(2,4,6)], ".sort.bed ", path, 
                   "bed/", prex[c(1,3,5)], ".sort.bed >", path, "bed/", c("n1", "n2", "n3"), ".join.bed", "\n")
  for (m in cmd_01) {cat(m, file = shel, append = T)}
  for (m in cmd_02) {cat(m, file = shel, append = T)}
  for (m in cmd_03) {cat(m, file = shel, append = T)}
  for (m in cmd_04) {cat(m, file = shel, append = T)}
  for (m in cmd_05) {cat(m, file = shel, append = T)}
  for (m in cmd_06) {cat(m, file = shel, append = T)}
  for (m in cmd_07) {cat(m, file = shel, append = T)}
  for (m in cmd_08) {cat(m, file = shel, append = T)}
  for (m in cmd_09) {cat(m, file = shel, append = T)}
  for (m in cmd_10) {cat(m, file = shel, append = T)}
  for (m in cmd_11) {cat(m, file = shel, append = T)}
  for (m in cmd_12) {cat(m, file = shel, append = T)}
  for (m in cmd_13) {cat(m, file = shel, append = T)}
  print(paste0("nohup bash ", shel, " >", paste0(path,"Log/"), basename(shel), ".log", " 2>&1 &"))
}
##
##
##
#04、分析数据：2C-like共有基因、统计TE、哪些基因含有MM RNA与MM或GSAT_MM互作的位置？
te01 <- rtracklayer::import(con = inx6, format = "bed")#inx6、inx8
te02 <- rtracklayer::import(con = inx3, format = "gtf")
ge01 <- rtracklayer::import(con = inx2, format = "gtf")
ge01 <- ge01[which(ge01$type=="gene")]
#
for (i in c("n1","n2","n3")) {
  rd01 <- fread(file = paste0(path,"bed/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
  rd01 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V2", start.field = "V3", end.field = "V4", strand.field = "V6", keep.extra.columns = T)
  ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = te01, subject = rd01, ignore.strand = T))
  rs01 <- as.data.frame(rd01[unique(subjectHits(ov01))])
  colnames(rs01) <- paste0("V", seq_along(colnames(rs01)))
  write.table(x = rs01, file = paste0(path, "specialOut/",i,".MM.dna.bed"), sep = "\t",
              quote = F, col.names = F, row.names = F)
  wr01 <- data.frame(chr = rs01$V8, start = rs01$V9, end = rs01$V10, name = ".", score = 1, strand = rs01$V12)
  write.table(x = wr01, file = paste0(path, "specialOut/",i,".MM.DNA.bed"), sep = "\t",
              quote = F, col.names = F, row.names = F)
  dtag <- unique(rs01[,c(6,8,9,10,12)])
  dtag <- makeGRangesFromDataFrame(df = dtag, seqnames.field = "V8", start.field = "V9",end.field = "V10", strand.field = "V12")
  ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = te02, subject = dtag, ignore.strand = T))
  tf01 <- as.data.frame(te02[queryHits(ov02)])
  tf01$name <- paste0(tf01$gene_id, "::", tf01$family_id, "::", tf01$class_id)
  tf02 <- as.data.frame(table(tf01$name))
  tf02 <- tf02[order(tf02$Freq, decreasing = T),]
  tf02[grep(x = tf02$Var1, pattern = "RLTR1B"),]
  write.table(x = tf02, file = paste0(path,"specialOut/",i,".MM.dtag.TE.txt"), sep = "\t", quote = F, row.names = F)
  
  
  tf01$name <- tf01$family_id
  tf02 <- as.data.frame(table(tf01$name))
  tf02 <- tf02[order(tf02$Freq, decreasing = T),]
  write.table(x = tf02, file = paste0(path,"specialOut/",i,".MM.dtag.family.txt"), sep = "\t", quote = F, row.names = F)
  tf01$name <- tf01$class_id
  tf02 <- as.data.frame(table(tf01$name))
  tf02 <- tf02[order(tf02$Freq, decreasing = T),]
  write.table(x = tf02, file = paste0(path,"specialOut/",i,".MM.dtag.class.txt"), sep = "\t", quote = F, row.names = F)
  rs02 <- makeGRangesFromDataFrame(df = rs01, keep.extra.columns = T, 
                                   seqnames.field = "V8", start.field = "V9", end.field = "V10", strand.field = "V12")
  ov03 <- suppressWarnings(GenomicRanges::findOverlaps(query = ge01, subject = rs02, ignore.strand = T))
  df01 <- as.data.frame(ge01[queryHits(ov03)])
  df02 <- as.data.frame(rs02[subjectHits(ov03)])
  fina <- cbind(df02, df01)
  colnames(fina) <- paste0("V", seq_along(colnames(fina)))
  out1 <- data.frame(rnaPos = paste0("chr",fina$V6,":",fina$V7,"-",fina$V8), 
                     rnaStrand = fina$V10,
                     dnaPos = paste0("chr",fina$V1,":",fina$V2,"-",fina$V3), 
                     dnaStrand = fina$V5, markID = fina$V11,
                     geneName = fina$V25, geneType = fina$V27, 
                     genePos = paste0("chr",fina$V14,":",fina$V15,"-",fina$V16))
  write.table(x = out1, file = paste0(path,"specialOut/",i,".MM.gene.txt"), sep = "\t", quote = F, row.names = F)
  inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
  wr02 <- out1[which(out1$geneName %in% inte$V4),]
  write.table(x = wr02, file = paste0(path,"specialOut/",i,".MM.2Clike.txt"), sep = "\t", quote = F, row.names = F)
  print(intersect(out1$geneName, inte$V4))
}
#
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
twoc <- intersect(intersect(ff01$geneName, ff02$geneName),ff03$geneName)
ff01 <- ff01[which(ff01$geneName %in% twoc),]
ff02 <- ff02[which(ff02$geneName %in% twoc),]
ff03 <- ff03[which(ff03$geneName %in% twoc),]
write.table(x = ff01, file = paste0(path,"specialOut/n1.MM.2Clike.ov.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = ff02, file = paste0(path,"specialOut/n2.MM.2Clike.ov.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = ff03, file = paste0(path,"specialOut/n3.MM.2Clike.ov.txt"), sep = "\t", quote = F, row.names = F)
#
tf03 <- as.data.frame(te02)
se01 <- tf03 %>% dplyr::group_by(class_id) %>% summarise(length = sum(width), count = n())
se02 <- tf03 %>% dplyr::group_by(family_id) %>% summarise(length = sum(width), count = n())
se03 <- tf03 %>% dplyr::group_by(gene_id) %>% summarise(length = sum(width), count = n())
se01 <- se01[order(se01$length, decreasing = T), ]
se02 <- se02[order(se02$length, decreasing = T), ]
se03 <- se03[order(se03$length, decreasing = T), ]
write.table(x = se01, file = paste0(path,"specialOut/classID.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = se02, file = paste0(path,"specialOut/familyID.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = se03, file = paste0(path,"specialOut/repeatID.txt"), sep = "\t", quote = F, row.names = F)
#
myte <- te02[which(te02$gene_id=="MM-int")]
lsff <- list.files(path = paste0(path, "specialOut"), pattern = "MM.DNA.bed", full.names = T)
for (i in seq(3)) {
  print(paste0("sample n",i, " :: GSAT_MM & MM DNA tag & all Genes"))
  dtag <- read.table(file = lsff[i], header = F, sep = "\t", stringsAsFactors = F)
  dtag <- makeGRangesFromDataFrame(df = dtag, seqnames.field = "V1", start.field = "V2", end.field = "V3", strand.field = "V6")
  ov04 <- suppressWarnings(GenomicRanges::findOverlaps(query = dtag, subject = myte, ignore.strand = T))
  dtdf <- dtag[queryHits(ov04)]
  ov05 <- suppressWarnings(GenomicRanges::findOverlaps(query = dtdf, subject = ge01, ignore.strand = T))
  df01 <- as.data.frame(dtdf[queryHits(ov05)])
  df02 <- as.data.frame(ge01[subjectHits(ov05)])
  resu <- cbind(df01, df02)
  assign(paste0("int",i), unique(resu$gene_name))
  print(unique(resu$gene_name))
  print("----------------------------------------------------")
  print("----------------------------------------------------")
}
fa01 <- intersect(intersect(int1, int2), int3)
fa02 <- intersect(intersect(int1, int2), int3)
need <- intersect(fa01, fa02)
##
##
##
#05、画图：韦恩图：与2C-like基因的重叠
temp <- list()
inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
temp[["a"]] = inte$V4
temp[["b"]] = unique(ff01$geneName)
temp[["c"]] = unique(ff02$geneName)
temp[["d"]] = unique(ff03$geneName)
p_01 <- VennDiagram::venn.diagram(x = temp, filename = NULL, 
                                  main = "RADICL-seq",
                                  sub = "MM RNA targeted 2C gene", 
                                  main.fontfamily = "serif", sub.fontfamily = "serif",
                                  cex = 1.5, main.cex = 2, sub.cex = 2,
                                  print.mode = "raw",#c("percent", "raw"),
                                  category.names = c("2C-Genes", "n1Sample", "n2Sample", "n3Sample"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC","#CCCCFF"),
                                  units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
pdf(file = paste0(path, "Plot/MMGeneVenn.pdf"))
grid.draw(p_01)
dev.off()
##
##
##
#05、画图：韦恩图：第二种方法
temp <- list()
inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
temp[["a"]] = inte$V4
temp[["b"]] = unique(ff01$geneName)
temp[["c"]] = unique(ff02$geneName)
temp[["d"]] = unique(ff03$geneName)
A = data.frame(x0 = 6, y0 = 6, r0 = 5)
B = data.frame(x0 = 6, y0 = 7.5, r0 = 2.7)
C = data.frame(x0 = 7.5, y0 = 5.1, r0 = 2.7)
D = data.frame(x0 = 4.5, y0 = 5.1, r0 = 2.7)
p_02 <- ggplot() + 
        geom_circle(data = A, aes(x0=x0, y0=y0, r=r0), fill="#FFFFCC", alpha = 0.5) +
        geom_circle(data = B, aes(x0=x0, y0=y0, r=r0), fill="#CCFFFF", alpha = 0.5) +
        geom_circle(data = C, aes(x0=x0, y0=y0, r=r0), fill="#CCCCFF", alpha = 0.4) +
        geom_circle(data = D, aes(x0=x0, y0=y0, r=r0), fill="#FFCCCC", alpha = 0.4) +
        coord_fixed() +
        theme_transparent()
ggsave(plot = p_02, filename = paste0(path, "Plot/TE-MMGeneVenn-v2.pdf"), width = 16, height = 16, units = "cm")
##
##
##
#06、画图：桑葚图：与2C-like基因的重叠
temp <- list()
inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
temp[["a"]] = inte$V4
temp[["b"]] = unique(ff01$geneName)
temp[["c"]] = unique(ff02$geneName)
temp[["d"]] = unique(ff03$geneName)
nb01 <- data.frame(twoc = temp[["a"]], twocM = "Yes", n1Sample = "No", n2Sample = "No", n3Sample = "No")
nb01$n1Sample[which(temp[["a"]] %in% temp[["b"]])] <- "Yes"
nb01$n2Sample[which(temp[["a"]] %in% temp[["c"]])] <- "Yes"
nb01$n3Sample[which(temp[["a"]] %in% temp[["d"]])] <- "Yes"
nb01$Freq <- 1
nb02 <- nb01 %>% dplyr::group_by(n1Sample,n2Sample,n3Sample) %>% dplyr::summarise(.groups = "keep", Freq=sum(Freq))
nb02$n1Sample <- factor(nb02$n1Sample, levels = c("Yes","No"))
p_03 <- ggplot(nb02, aes(y = Freq, axis1 = n1Sample, axis2 = n2Sample, axis3 = n3Sample)) +
        geom_alluvium(aes(fill = n3Sample), size = 0, width = 0, knot.pos = 0, reverse = T) +
        scale_fill_manual(values = c("#FF5F5D","#3F7C85"), 
                          labels=c("Targeted by MM","Not Targeted by MM")) +
        labs(fill = "2C-like Genes", y = "2C-like Genes Count\n") +
        theme_void() +
        theme(axis.text = element_text(family = "serif", size=16),
              axis.title = element_text(family = "serif", size=16, angle = 90),
              legend.text = element_text(family = "serif", size=16),
              legend.title = element_text(family = "serif", size=16),
              plot.margin = unit (c(1,1,1,1), 'cm')) +
        geom_stratum(width = 1/6, reverse = F, color = "white", size = 1, fill = c("grey50","grey90","grey50","grey90","grey50","grey90")) +
        scale_x_continuous(breaks = 1:3, labels = c("n1Sample", "n2Sample","n3Sample"))
ggsave(plot = p_03, filename = paste0(path, "Plot/MMGeneSankey.pdf"), width = 22, height = 15, units = "cm")
##
##
##
#07、画图：点图：哪些TE富集MM RNA
pt01 <- read.table(file = paste0(path,"specialOut/repeatID.txt"), stringsAsFactors = F, sep = "\t", header = T)
pt02 <- read.table(file = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa.fai",stringsAsFactors = F, sep = "\t", header = F)
refL <- sum(pt02$V2)
pt01$frac <- pt01$length/refL
pt03 <- read.table(file = paste0(path,"specialOut/n1.MM.dtag.TE.txt"), stringsAsFactors = F, sep = "\t", header = T)
pt03 <- pt03[1:10,]
pt03$name <- sapply(str_split(pt03$Var1, pattern = "::"),"[",1)
pt03$frac <- NA
for (i in seq_along(pt03$name)) {
  pt03$frac[i] <- pt01[which(pt01$gene_id == pt03$name[i]),"frac"]
}
pt04 <- read.table(file = paste0(path,"specialOut/n2.MM.dtag.TE.txt"), stringsAsFactors = F, sep = "\t", header = T)
pt05 <- read.table(file = paste0(path,"specialOut/n3.MM.dtag.TE.txt"), stringsAsFactors = F, sep = "\t", header = T)
pt03$erro <- NA; pt03$Freq2 <- NA; pt03$Freq3 <- NA; pt03$FreqR
for (i in seq_along(pt03$name)) {
  pt03$Freq2[i] <- pt04[which(pt04$Var1 == pt03$Var1[i]), "Freq"]
  pt03$Freq3[i] <- pt04[which(pt05$Var1 == pt03$Var1[i]), "Freq"]
  pt03$FreqR[i] <- mean(c(pt03$Freq[i],pt03$Freq2[i],pt03$Freq3[i]))
  pt03$erro[i] <- sd(c(pt03$Freq[i],pt03$Freq2[i],pt03$Freq3[i]))/sqrt(3)
}
pt03$family <- sapply(str_split(pt03$Var1, pattern = "::"),"[",2)
pt03$name <- factor(x = pt03$name, levels = pt03$name)
print(pt03$name)
p_04 <- ggplot(data = pt03) + 
        geom_errorbar(mapping = aes(x = name, y = FreqR, ymin = FreqR -erro, ymax = FreqR +erro), alpha =0.5, width = 0.1) +
        geom_point(mapping = aes(x = name, y = FreqR, size = frac, color = family)) +
        theme_classic() +
        scale_size_continuous(limits = c(0,0.03), breaks = c(0.001, 0.01, 0.02)) +
        scale_y_continuous(limits = c(450,1500), breaks = seq(500,1500,100)) +
        labs(x = "", 
             y = "Interaction number of MM with TE", title = "This is TITLE...", 
             size = "Sum(TE Length)/Genome Length", color = "TE Family") +
        theme(axis.title = element_text(family = "serif", size = 12, colour = "black"),
              axis.text.y = element_text(family = "serif", size = 12, colour = "black"),
              axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
              legend.text = element_text(family = "serif", size = 12, colour = "black"),
              legend.title = element_text(family = "serif", size = 12, colour = "black"),
              legend.position = c(0.8,0.7),
              plot.margin = margin(1,1,1,1, unit = "cm"))
ggsave(plot = p_04, filename = paste0(path, "Plot/TE-Erich.pdf"), width = 16, height = 16, units = "cm")
##
##
##
##
##
##
##
##
##
#08、分析数据：MM RNA靶向哪些可变剪接基因
stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
stag <- stag[grep(x = rownames(stag), pattern = ":", invert = T), -7]
colnames(stag) <- sapply(str_split(colnames(stag), ".bam."), "[", 2)
stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))#write.table(x = stag, file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/stageExp.csv", quote = F, col.names = T, row.names = F, sep = ",")
#
sf01 <- read.table(file = inx9, stringsAsFactors = F, sep = "\t", header = F)
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
wr01 <- ff01[which(ff01$geneName %in% sf01$V1),]
wr02 <- ff02[which(ff02$geneName %in% sf01$V1),]
wr03 <- ff03[which(ff03$geneName %in% sf01$V1),]
write.table(x = wr01, file = paste0(path, "spliceFactor/n1.MM.SF.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = wr02, file = paste0(path, "spliceFactor/n2.MM.SF.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = wr03, file = paste0(path, "spliceFactor/n3.MM.SF.txt"), sep = "\t", quote = F, row.names = F)
#
sfgg <- intersect(intersect(wr01$geneName, wr02$geneName), wr03$geneName)
write.table(x = wr01[which(wr01$geneName %in% sfgg),], file = paste0(path, "spliceFactor/n1.MM.SF.ov.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = wr02[which(wr02$geneName %in% sfgg),], file = paste0(path, "spliceFactor/n2.MM.SF.ov.txt"), sep = "\t", quote = F, row.names = F)
write.table(x = wr03[which(wr03$geneName %in% sfgg),], file = paste0(path, "spliceFactor/n3.MM.SF.ov.txt"), sep = "\t", quote = F, row.names = F)
##
##
##
#09、不同样本的可变剪接基因集的韦恩关系
temp <- list()
temp[["a"]] <- sf01$V1
temp[["b"]] <- unique(wr01$geneName)
temp[["c"]] <- unique(wr02$geneName)
temp[["d"]] <- unique(wr03$geneName)
p_05 <- VennDiagram::venn.diagram(x = temp, filename = NULL, 
                                  main = "2% FA RADICL",
                                  sub = "MM RNA targeted SF Gene", 
                                  main.fontfamily = "serif", sub.fontfamily = "serif", 
                                  cex = 1.5, main.cex = 2, sub.cex = 2,
                                  print.mode = "raw",#c("percent", "raw"),
                                  category.names = c("SF Gene Set", "n1Sample", "n2Sample", "n3Sample"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC","#CCCCFF"),
                                  units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
pdf(file = paste0(path, "spliceFactor/MM-SF-Venn.pdf"))
grid.draw(p_05)
dev.off()
##
##
##
#10、柱状图：不同样本的被MM靶向的可变剪接基因的数量
pt04 <- data.frame(name = c(rep("all SF",2), rep("n1Sample",2), rep("n2Sample",2), rep("n3Sample",2)),
                   cate = rep(c("Yes","No"),4),
                   numb = c(length(sf01$V1),length(sf01$V1)-length(sf01$V1),
                            length(unique(wr01$geneName)), length(sf01$V1)-length(unique(wr01$geneName)),
                            length(unique(wr02$geneName)), length(sf01$V1)-length(unique(wr02$geneName)),
                            length(unique(wr03$geneName)), length(sf01$V1)-length(unique(wr03$geneName))))
p_06 <- ggplot(data = pt04) +
        geom_bar(mapping = aes(x = name, weight = numb, fill = cate), position = "stack", width = 0.5) +
        scale_fill_manual(values = c('#007172','#F29325')) +
        labs(title = "2% FA RADICL\nMM RNA target Splice Factor Gene",
             x="Splice Factor Gene Set", y="Splice Factor Gene Number",fill="Does Targeted \nby MM RNA?") +
        theme_classic() +
        theme(plot.margin = margin(1,1,1,1,unit="cm"),
              axis.title  = element_text(family = "serif", size = 12),
              axis.text.y = element_text(family = "serif", size = 12),
              axis.text.x = element_text(family = "serif", size = 12),
              legend.text = element_text(family = "serif", size = 12),
              legend.title = element_text(family = "serif", size = 12)) +
        geom_label(data = pt04[which(pt04$cate=="Yes"),], 
                   mapping = aes(x = name, y = numb, label=numb))
ggsave(plot = p_06, filename = paste0(path, "spliceFactor/SF-Number.pdf"), width = 16, height = 16, units = "cm")
##
##
##
#11、柱状图：不同样本的被MM靶向的可变剪接基因的数量：Top10
sf01 <- sf01
sf02 <- as.data.frame(table(wr01$geneName)); sf02 <- sf02[order(sf02$Freq, decreasing = T),]
sf03 <- as.data.frame(table(wr02$geneName)); sf03 <- sf03[order(sf03$Freq, decreasing = T),]
sf04 <- as.data.frame(table(wr03$geneName)); sf04 <- sf04[order(sf04$Freq, decreasing = T),]
pt05 <- data.frame()
for (i in seq_along(sf02$Var1[1:10])) {
  if (as.character(sf02$Var1[i]) %in% sf03$Var1) {
    fq02 = sf03[which(sf03$Var1==as.character(sf02$Var1[i])),"Freq"]
  }
  else {
    fq02 = NA
  }
  if (as.character(sf02$Var1[i]) %in% sf04$Var1) {
    fq03 = sf04[which(sf04$Var1==as.character(sf02$Var1[i])),"Freq"]
  }
  else {
    fq03 = NA
  }
  temp <- data.frame(name = as.character(sf02$Var1[i]),
                     fq01 = sf02$Freq[i], fq02 = fq02, fq03 = fq03)
  pt05 <- rbind(pt05, temp)
}
for (i in seq_along(pt05$name)) {
  if (TRUE %in% is.na(pt05[i,])) {
    pt05[i,][is.na(pt05[i,])] <- mean(as.numeric(pt05[i,2:4]), na.rm = T)
  }
}
for (i in seq_along(pt05$name)) {
  pt05$freq[i] <- rowMeans(pt05[i,2:4])
  pt05$erro[i] <- sd(pt05[i,2:4])/sqrt(3)
}
pt05 <- pt05[order(pt05$freq, decreasing = T),]
pt05$name <- factor(pt05$name, levels = pt05$name)
p_07 <- ggplot(data = pt05) + 
        geom_errorbar(mapping = aes(x = name, y = freq, ymin = freq-erro, ymax = freq+erro), linewidth=1, width = 0.2) +
        geom_bar(mapping = aes(x = name, y = freq), stat = "identity", linewidth = 1, fill = "steelblue") +
        theme_classic() +
        scale_size_continuous(limits = c(0,0.03), breaks = c(0.001, 0.01, 0.02)) +
        scale_y_continuous(expand = c(0.02,0), limits = c(0,20), breaks = seq(0,18,2)) +
        labs(x = "", y = "Interaction number of MM with SF", 
             title = "Interaction number of MM with Splice Factor (Top 10)") +
        theme(axis.title = element_text(family = "serif", size = 12, colour = "black"),
              axis.text.y = element_text(family = "serif", size = 12, colour = "black"),
              axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
              legend.text = element_text(family = "serif", size = 12, colour = "black"),
              legend.title = element_text(family = "serif", size = 12, colour = "black"),
              legend.position = c(0.8,0.7),
              plot.margin = margin(1,1,1,1, unit = "cm"))
p_07
ggsave(plot = p_07, filename = paste0(path, "spliceFactor/SF-Freq.pdf"), width = 14, height = 12, units = "cm")
##
##
##
#12、折线图：Top 10剪接基因的表达变化 (分面展示)
pt06 <- data.frame()
for (i in seq_along(pt05$name)) {
  temp <- stag[which(rownames(stag) == pt05$name[i]),]
  pt06 <- rbind(pt06, temp)
}
pt06$name <- rownames(pt06)
pt06 <- tidyr::gather(data = pt06, key = "stage", value = "expres", -"name")
pt06$name <- factor(pt06$name, levels = pt06$name[1:10])
pt06$stage <- factor(x = pt06$stage, levels = unique(pt06$stage))
p_08 <- ggplot(data = pt06) + 
        geom_line(mapping = aes(x = stage, y = expres, group = name)) +
        geom_point(mapping = aes(x = stage, y = expres)) +
        labs(x="", y="CPM Expression",
             title = "CPM: Top 10 Splice Factor (by the number of MM target)") +
        facet_wrap(.~ name, scales = "free_y") +
        theme(plot.margin = margin(1,1,1,1,unit = "cm"),
              axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
              strip.text.x = element_text(size = 12, color = "red", face = "bold"))
p_08
ggsave(plot = p_08, filename = paste0(path, "spliceFactor/SF-Top10-Expr.pdf"), width = 18, height = 15, units = "cm")
##
##
##
#13、热图：CPM heatmap of Splice Factors
stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/abRepeats/count/star/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
stag <- stag[grep(x = rownames(stag), pattern = ":", invert = T), -7]
colnames(stag) <- sapply(str_split(colnames(stag), ".bam."), "[", 2)
stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))#write.table(x = stag, file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/stageExp.csv", quote = F, col.names = T, row.names = F, sep = ",")
#
sf01 <- read.table(file = inx9, stringsAsFactors = F, sep = "\t", header = F)
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
wr01 <- ff01[which(ff01$geneName %in% sf01$V1),]
wr02 <- ff02[which(ff02$geneName %in% sf01$V1),]
wr03 <- ff03[which(ff03$geneName %in% sf01$V1),]
wr04 <- intersect(wr01$geneName, intersect(wr02$geneName, wr03$geneName))
sf02 <- stag[which(rownames(stag) %in% sf01$V1),]
#
sf01 <- sf01
wr01 <- wr01
wr02 <- wr02
wr03 <- wr03
pt07 <- data.frame()
for (i in seq_along(sf01$V1)) {
  name = sf01$V1[i]
  if (as.character(sf01$V1[i]) %in% wr01$geneName) {
    frq1 = length(grep(x = wr01$geneName, pattern = as.character(sf01$V1[i])))
  }
  else {
    frq1 = 0
  }
  if (as.character(sf01$V1[i]) %in% wr02$geneName) {
    frq2 = length(grep(x = wr02$geneName, pattern = as.character(sf01$V1[i])))
  }
  else {
    frq2 = 0
  }
  if (as.character(sf01$V1[i]) %in% wr03$geneName) {
    frq3 = length(grep(x = wr03$geneName, pattern = as.character(sf01$V1[i])))
  }
  else {
    frq3 = 0
  }
  mean <- mean(c(frq1,frq2,frq3))
  temp <- data.frame(name = name, frq1 = frq1, frq2 = frq2, frq3 = frq3, mean = mean)
  pt07 <- rbind(pt07, temp)
}
pdf(file = paste0(path, "spliceFactor/SF-Heatmap.pdf"), width = 10, height = 10)
#
tmp1 <- sf02[which(rownames(sf02) %in% unique(c(wr01$geneName,wr02$geneName,wr03$geneName))),]
tmp2 <- sf02[which(!(rownames(sf02) %in% unique(c(wr01$geneName,wr02$geneName,wr03$geneName)))),]
sf03 <- rbind(tmp1, tmp2)
num1 <- pt07[which(pt07$name %in% rownames(tmp1)),]
num2 <- pt07[which(pt07$name %in% rownames(tmp2)),]
annotation_row <- data.frame(class = c(rep("SF_MM",187), rep("SF_MM_NO",110)),
                             numberMM = c(num1$mean, num2$mean))
rownames(annotation_row) <- rownames(sf03)
pheatmap(mat = sf03,
         main="\nCPM heatmap of Splice Factors",
         scale = "row", cellwidth = 60, 
         fontsize_col = 14,
         annotation_row = annotation_row, annotation_names_row = F,
         annotation_colors = list(class = c("SF_MM" = "#007172", "SF_MM_NO"="#F29325")),
         show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
         labels_col = colnames(sf03), angle_col = "45",
         breaks = seq(-2,2, by=0.01), border = F,
         color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                    colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
#
tmp1 <- tmp1
annotation_row <- data.frame(class = c(rep("SF_MM",187)),numberMM = num1$mean)
rownames(annotation_row) <- rownames(tmp1)
pheatmap(mat = tmp1,
         main="\nCPM heatmap of Splice Factors (187)",
         scale = "row", cellwidth = 60, 
         fontsize_col = 14,
         annotation_row = annotation_row, annotation_names_row = F,
         annotation_colors = list(class = c("SF_MM" = "#007172")),
         show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
         labels_col = colnames(sf03), angle_col = "45",
         breaks = seq(-2,2, by=0.01), border = F,
         color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                    colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
#
tmp3 <- sf02[which(rownames(sf02) %in% wr04),]
tmp4 <- sf02[which(!(rownames(sf02) %in% wr04)),]
sf04 <- rbind(tmp3, tmp4)
num1 <- pt07[which(pt07$name %in% rownames(tmp3)),]
num2 <- pt07[which(pt07$name %in% rownames(tmp4)),]
annotation_row <- data.frame(class = c(rep("SF_MM_46",46), rep("SF_Other",251)),
                             numberMM = c(num1$mean, num2$mean))
rownames(annotation_row) <- rownames(sf03)
pheatmap(mat = sf04,
         main="\nCPM heatmap of Splice Factors",
         scale = "row", cellwidth = 60, 
         fontsize_col = 14,
         annotation_row = annotation_row, annotation_names_row = F,
         annotation_colors = list(class = c("SF_MM_46" = "#007172", "SF_Other"="#F29325")),
         show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
         labels_col = colnames(sf04), angle_col = "45",
         breaks = seq(-2, 2, by=0.01), border = F,
         color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                    colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
dev.off()
##
##
##
##
##
##
##
##
##
##
##
##
#14、内含子长度比较：数据获取
ge01 <- rtracklayer::import(con = inx2, format = "gtf")
ge02 <- ge01[which(ge01$type=="gene")]
ge03 <- ge01[which(ge01$type=="exon")]
my01 <- data.frame()
for (i in unique(ge01$gene_name)) {
  print(which(unique(ge01$gene_name)==i))
  tmp1 <- as.data.frame(ge02[which(ge02$gene_name==i)])[,c(1:3,12,8,5)]
  tmp1$score <- "."
  tmp1 <- bedtoolsr::bt.sort(tmp1)
  tmp1 <- bedtoolsr::bt.merge(tmp1)
  tmp2 <- as.data.frame(ge03[which(ge03$gene_name==i)])[,c(1:3,12,8,5)]
  tmp2$score <- "."
  tmp2 <- bedtoolsr::bt.sort(tmp2)
  bedtoolsr::bt.merge(tmp2)
  subt <- bedtoolsr::bt.subtract(a = tmp1, b = tmp2)
  genL <- sum(tmp1$V3-tmp1$V2)
  exoL <- sum(tmp2$V3-tmp2$V2)
  intL <- sum(subt$V3-subt$V2)
  temp <- data.frame(name = i, geneL = genL, exonL = exoL, intoL = intL)
  my01 <- rbind(my01, temp)
}
base::save(my01, file = "/Reference/aaSHY/Ensembl99/GRCm38/RData/featureL.RData")
##
##
##
#15、内含子长度比较：关于剪接因子
my02 <- my01
my02$mark <- NA
nam1 <- as.data.frame(table(c(wr01$geneName,wr02$geneName,wr03$geneName)))#所有MM结合的剪接基因
nam1 <- as.character(nam1[order(nam1$Freq, decreasing = T),1])
tmp1 <- my02[which(my02$name %in% nam1),]
tmp1$mark <- "SF_MM"
nam2 <- nam1[1:ceiling(length(nam1)*0.5)]
tmp2 <- my02[which(my02$name %in% nam2),]
tmp2$mark <- "SF_MM_TopP50"
nam3 <- sf01$V1
tmp3 <- my02[which(my02$name %in% nam3),]
tmp3$mark <- "SF_Gene"
pt08 <- rbind(tmp1, tmp2, tmp3)
nam4 <- unique(ge02[grep(ge02$gene_biotype, pattern = "protein_coding")]$gene_name)
nam4 <- setdiff(nam4, nam1)
for (i in seq(3)) {
  namX <- sample(nam4, 1000, replace = F)
  tmp4 <- my02[which(my02$name %in% namX),]
  tmp4$mark <- paste0("RandomSampleN",i)
  pt08 <- rbind(pt08, tmp4)
}
#
pt08$mark <- factor(pt08$mark, levels = unique(pt08$mark))
p_09 <- ggplot(data = pt08) + 
        geom_jitter(aes(x = mark, y = log2(intoL+1)), alpha = 0.1) +
        geom_boxplot(mapping = aes(x = mark, y = log2(intoL+1)), fill = "#F29325") + 
        theme_classic() +
        scale_y_continuous(limits = c(0, 35)) +
        labs(x="Gene Set", y="Log2 (SUM Length of Introns +1)", title = "Length Distribution of Introns of Specific Genes") +
        theme(plot.margin = margin(1,1,1,1,unit = "cm"),
              plot.title = element_text(family = "serif", size = 15),
              axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
              axis.text.y = element_text(family = "serif", size = 12),
              axis.title = element_text(family = "serif", size = 12)) +
        geom_signif(mapping = aes(x = mark, y = log2(intoL+1)),
                    comparisons = list(c("SF_MM_TopP50", "SF_Gene"),
                                       c("SF_MM_TopP50", "RandomSampleN1"),
                                       c("SF_MM_TopP50", "RandomSampleN2"),
                                       c("SF_MM_TopP50", "RandomSampleN3")),
                    map_signif_level=T, textsize=6, test = wilcox.test, step_increase = 0.12)
ggsave(plot = p_09, filename = paste0(path, "spliceFactor/SF-Box.pdf"), width = 18, height = 17, units = "cm")
##
##
##
#16、内含子长度比较：所有MM结合的基因
my03 <- my01
my03$mark <- NA
inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff03 <- read.table(file = paste0(path, "specialOut/n3.MM.2Clike.txt"), sep = "\t", stringsAsFactors = F, header = T)
nam1 <- as.data.frame(table(c(ff01$geneName,ff02$geneName,ff03$geneName)))#所有MM结合的2C基因
nam1 <- as.character(nam1[order(nam1$Freq, decreasing = T),1])
tmp1 <- my03[which(my03$name %in% nam1),]
tmp1$mark <- "TWOC_MM"
nam2 <- nam1[1:ceiling(length(nam1)*0.5)]
tmp2 <- my03[which(my03$name %in% nam2),]
tmp2$mark <- "TWOC_MM_TopP50"
nam3 <- inte$V4
tmp3 <- my03[which(my03$name %in% nam3),]
tmp3$mark <- "TWOC_Gene"
pt09 <- rbind(tmp1, tmp2, tmp3)
nam4 <- unique(ge02[grep(ge02$gene_biotype, pattern = "protein_coding")]$gene_name)
nam4 <- setdiff(nam4, nam1)
for (i in seq(3)) {
  namX <- sample(nam4, 1000, replace = F)
  tmp4 <- my03[which(my03$name %in% namX),]
  tmp4$mark <- paste0("RandomSampleN",i)
  pt09 <- rbind(pt09, tmp4)
}
#
#
ff04 <- read.table(file = paste0(path, "specialOut/n1.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff05 <- read.table(file = paste0(path, "specialOut/n2.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
ff06 <- read.table(file = paste0(path, "specialOut/n3.MM.gene.txt"), sep = "\t", stringsAsFactors = F, header = T)
nam5 <- as.data.frame(table(c(ff04$geneName,ff05$geneName,ff06$geneName)))#所有MM结合的基因
nam5 <- as.character(nam5[order(nam5$Freq, decreasing = T),1])[1:1000]
tmp5 <- my03[which(my03$name %in% nam5),]
tmp5$mark <- "Gene_MM_Top1000"
pt09 <- rbind(pt09, tmp5)
#
#
pt09$mark <- factor(pt09$mark, levels = unique(pt09$mark))
p_10 <- ggplot(data = pt09) + 
        geom_jitter(aes(x = mark, y = log2(intoL+1)), alpha = 0.1) +
        geom_boxplot(mapping = aes(x = mark, y = log2(intoL+1)), fill = "#F29325") + 
        theme_classic() +
        scale_y_continuous(limits = c(0, 35)) +
        labs(x="Gene Set", y="Log2 (SUM Length of Introns +1)", title = "Length Distribution of Introns of Specific Genes") +
        theme(plot.margin = margin(1,1,1,1,unit = "cm"),
              plot.title = element_text(family = "serif", size = 15),
              axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
              axis.text.y = element_text(family = "serif", size = 12),
              axis.title = element_text(family = "serif", size = 12)) +
        geom_signif(mapping = aes(x = mark, y = log2(intoL+1)),
                    comparisons = list(c("TWOC_MM_TopP50", "TWOC_Gene"),
                                       c("TWOC_MM_TopP50", "RandomSampleN1"),
                                       c("TWOC_MM_TopP50", "RandomSampleN2"),
                                       c("TWOC_MM_TopP50", "RandomSampleN3")),
                    map_signif_level=T, textsize=6, test = wilcox.test, step_increase = 0.12)
p_10
ggsave(plot = p_10, filename = paste0(path, "specialOut/TWOC-Box.pdf"), width = 18, height = 17, units = "cm")
#
#
transPlotR::trancriptVis(gtfFile = inx2, gene = "Inpp4b", Chr = 1, collapse = F, exonWidth = 0.5)





