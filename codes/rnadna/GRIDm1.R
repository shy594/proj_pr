#Write by ~~~ on 2023.09.11
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","VennDiagram","stringr","Biostrings","tidyr","pheatmap",
              "enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer",
              "Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2",
              "ggrepel","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  inx5 = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/GenBank/rRNA.fa"
  inx6 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  inx7 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/special/twoC.bed"
  inx8 = "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT_Mm::ERVL::LTR.bed"
}
################################################################################
#03、跑命令
{
  shel <- paste0(path,"src/run_b.sh")
  cat("#!/bin/bash\n", file = shel)
  zemp <- list.files(path = paste0(path, "rawdata"), full.names = T)
  prex <- list.files(path = paste0(path, "rawdata"))
  cmd_01 <- paste0("#parallel-fastq-dump -s ",zemp," ",#--gzip 
                   "-t 10 --split-3 -O ", path,"fq/","\n")
  cmd_02 <- paste0("#rename s/fastq/fq/ ", path,"fq/*.fastq","\n")
  pair <- paste0(path,"zemp/pair.txt"); cat("", file = pair)
  for (m in paste0(path,"fq/",prex,".fq\n")) {cat(m, file =pair, append = T)}
  cmd_03 <- paste0("#cat ",pair," |head -n 2 >",path,"zemp/ff01.txt","\n")
  cmd_04 <- paste0("#cat ",pair," |tail -n 2 >",path,"zemp/ff02.txt","\n")
  cmd_05 <- paste0("#fastuniq -i ", path, "zemp/", c("ff01", "ff02"), 
                   ".txt", " -t q -o ", path, "fq/rep",c(1,2),"cdna.s1.fq", " -p ", path,
                   "fq/rep",c(1,2), "gdna.s1.fq", "\n")
  cmd_06 <- paste0("#rRNAdust ", inx5, " ", path, "fq/", prex, ".s1.fq", " -o ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_07 <- paste0("#bwa aln -t 10 -f ",path,"bam/", prex, ".s2.sai ", inx4, " ", path, "fq/", prex, ".s2.fq", "\n")
  cmd_08 <- paste0("#bwa samse ", inx4, " ", path, "bam/", prex, ".s2.sai ", path, "fq/", prex, ".s2.fq ",
                   "|samtools view -bS -F 4 -@ 10 |samtools sort -@ 10 >", path, "bam/", prex, ".s2.bam","\n")
  cmd_09 <- paste0("#samtools index ",path, "bam/", prex, ".s2.bam", "\n")
  cmd_10 <- paste0("#samtools view -@ 10 -q 0 -bS ", path, "bam/", prex, 
                   ".s2.bam ", "|bedtools bamtobed -i - >", path, "bed/", prex, ".s2.bed", "\n")
  cmd_11 <- paste0("cat ", path, "bed/", prex, ".s2.bed |sed 's/", prex, "/abBNA/g' |",
                   "sort -k 4b,4 >", path, "bed/", prex, ".sort.bed", "\n")
  cmd_12 <- paste0("join -t $'\\t' -j 4 ", path, "bed/", prex[c(1,3)], ".sort.bed ", path, 
                   "bed/", prex[c(2,4)], ".sort.bed >", path, "bed/", c("n1", "n2"), ".join.bed", "\n")
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
  print(paste0("nohup bash ", shel, " >", paste0(path,"Log/"), basename(shel), ".log", " 2>&1 &"))
}
################################################################################
#04、分析数据
{
  te01 <- rtracklayer::import(con = inx6, format = "bed")#inx6、inx8
  te02 <- rtracklayer::import(con = inx3, format = "gtf")
  ge01 <- rtracklayer::import(con = inx2, format = "gtf")
  ge01 <- ge01[which(ge01$type=="gene")]
  #
  for (i in c("n1","n2")) {
    rd01 <- fread(file = paste0(path,"bed/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
    rd01 <- makeGRangesFromDataFrame(df = rd01, seqnames.field = "V2", 
                                     start.field = "V3", end.field = "V4", 
                                     strand.field = "V6", keep.extra.columns = T)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = te01, 
                                                         subject = rd01, ignore.strand = T))
    rs01 <- as.data.frame(rd01[unique(subjectHits(ov01))])
    colnames(rs01) <- paste0("V", seq_along(colnames(rs01)))
    write.table(x = rs01, 
                file = paste0(path, "specialOut/",i,".MM.dna.bed"), sep = "\t",
                quote = F, col.names = F, row.names = F)
    wr01 <- data.frame(chr = rs01$V8, start = rs01$V9, end = rs01$V10, name = ".", score = 1, strand = rs01$V12)
    write.table(x = wr01, 
                file = paste0(path, "specialOut/",i,".MM.DNA.bed"), sep = "\t",
                quote = F, col.names = F, row.names = F)
    dtag <- unique(rs01[,c(6,8,9,10,12)])
    dtag <- makeGRangesFromDataFrame(df = dtag, seqnames.field = "V8", start.field = "V9",end.field = "V10", strand.field = "V12")
    #dtag <- dtag + 1000
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = te02, subject = dtag, ignore.strand = T))
    tf01 <- as.data.frame(te02[queryHits(ov02)])
    tf01$name <- paste0(tf01$gene_id, "::", tf01$family_id, "::", tf01$class_id)
    tf02 <- as.data.frame(table(tf01$name))
    tf02 <- tf02[order(tf02$Freq, decreasing = T),]
    write.table(x = tf02, 
                file = paste0(path,"specialOut/",i,".MM.dtag.TE.txt"), sep = "\t", quote = F, row.names = F)
    tf01$name <- tf01$family_id
    tf02 <- as.data.frame(table(tf01$name))
    tf02 <- tf02[order(tf02$Freq, decreasing = T),]
    write.table(x = tf02, 
                file = paste0(path,"specialOut/",i,".MM.dtag.family.txt"), sep = "\t", quote = F, row.names = F)
    tf01$name <- tf01$class_id
    tf02 <- as.data.frame(table(tf01$name))
    tf02 <- tf02[order(tf02$Freq, decreasing = T),]
    write.table(x = tf02, 
                file = paste0(path,"specialOut/",i,".MM.dtag.class.txt"), sep = "\t", quote = F, row.names = F)
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
    write.table(x = out1, file = paste0(path,"specialOut/",i,".MM.gene.txt"), 
                sep = "\t", quote = F, row.names = F)
    inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
    wr02 <- out1[which(out1$geneName %in% inte$V4),]
    write.table(x = wr02, file = paste0(path,"specialOut/",i,".MM.2Clike.txt"), 
                sep = "\t", quote = F, row.names = F)
    print(intersect(out1$geneName, inte$V4))
  }
}
##2C-like 共有基因
{
  ff01 <- read.table(file = paste0(path, "specialOut/n1.MM.2Clike.txt"), sep = "\t", 
                     stringsAsFactors = F, header = T)
  ff02 <- read.table(file = paste0(path, "specialOut/n2.MM.2Clike.txt"), sep = "\t", 
                     stringsAsFactors = F, header = T)
  twoc <- intersect(ff01$geneName, ff02$geneName)
  ff01 <- ff01[which(ff01$geneName %in% twoc),]
  ff02 <- ff02[which(ff02$geneName %in% twoc),]
  write.table(x = ff01, 
              file = paste0(path,"specialOut/n1.MM.2Clike.ov.txt"), 
              sep = "\t", quote = F, row.names = F)
  write.table(x = ff02, 
              file = paste0(path,"specialOut/n2.MM.2Clike.ov.txt"), 
              sep = "\t", quote = F, row.names = F)
}
##统计TE
{
  tf03 <- as.data.frame(te02)
  se01 <- tf03 %>% dplyr::group_by(class_id) %>% summarise(length = sum(width), count = n())
  se02 <- tf03 %>% dplyr::group_by(family_id) %>% summarise(length = sum(width), count = n())
  se03 <- tf03 %>% dplyr::group_by(gene_id) %>% summarise(length = sum(width), count = n())
  se01 <- se01[order(se01$length, decreasing = T), ]
  se02 <- se02[order(se02$length, decreasing = T), ]
  se03 <- se03[order(se03$length, decreasing = T), ]
  write.table(x = se01, 
              file = paste0(path,"specialOut/classID.txt"), sep = "\t", quote = F, row.names = F)
  write.table(x = se02, 
              file = paste0(path,"specialOut/familyID.txt"), sep = "\t", quote = F, row.names = F)
  write.table(x = se03, 
              file = paste0(path,"specialOut/repeatID.txt"), sep = "\t", quote = F, row.names = F)
}
################################################################################
#画图：韦恩图：与2C-like基因的重叠
{
  temp <- list()
  inte <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
  temp[["x"]] = inte$V4
  temp[["y"]] = ff01$geneName
  temp[["z"]] = ff02$geneName
  p_01 <- VennDiagram::venn.diagram(x = temp, filename = NULL, 
                                    main = "GRID-seq",
                                    sub = "MM RNA targeted 2C gene", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("2ClikeGenes", "n1 Sample", "n2 Sample"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "Plot/MMGeneVenn.pdf"))
  grid.draw(p_01)
  dev.off()
}
################################################################################
#画图：点图
{
  pt01 <- read.table(file = paste0(path,"specialOut/repeatID.txt"), stringsAsFactors = F, sep = "\t", header = T)
  pt02 <- read.table(file = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa.fai",stringsAsFactors = F, sep = "\t", header = F)
  refL <- sum(pt02$V2)
  pt01$frac <- pt01$length/refL
  pt03 <- read.table(file = paste0(path,"specialOut/n1.MM.dtag.TE.txt"), 
                     stringsAsFactors = F, sep = "\t", header = T)
  pt03 <- pt03[1:10,]
  pt03$name <- sapply(str_split(pt03$Var1, pattern = "::"),"[",1)
  pt03$frac <- NA
  for (i in seq_along(pt03$name)) {
    pt03$frac[i] <- pt01[which(pt01$gene_id == pt03$name[i]),"frac"]
  }
  pt04 <- read.table(file = paste0(path,"specialOut/n2.MM.dtag.TE.txt"), 
                     stringsAsFactors = F, sep = "\t", header = T)
  pt03$erro <- NA; pt03$Freq2 <- NA
  for (i in seq_along(pt03$name)) {
    pt03$Freq2[i] <- pt04[which(pt04$Var1 == pt03$Var1[i]), "Freq"]
    pt03$erro[i] <- sd(c(pt03$Freq[i],pt03$Freq2[i]))/sqrt(2)
    pt03$FreqR[i] <- mean(c(pt03$Freq[i],pt03$Freq2[i]))
  }
  pt03$family <- sapply(str_split(pt03$Var1, pattern = "::"),"[",2)
  pt03$name <- factor(x = pt03$name, levels = pt03$name)
  print(pt03$name)
  p_02 <- ggplot(data = pt03) + 
    geom_errorbar(aes(x = name, y = Freq, ymin = Freq -erro, ymax = Freq +erro), alpha =0.5, width = 0.1) +
    geom_point(mapping = aes(x = name, y = Freq, size = frac, color = family)) +
    theme_classic() +
    scale_size_continuous(limits = c(0.001,0.03), breaks = c(0.001, 0.01, 0.02)) +
    scale_y_continuous(limits = c(1000,7000), breaks = seq(1000,7000,1000)) +
    labs(x = "", 
         y = "Interaction number of MM with TE", title = "This is TITLE...", 
         size = "Sum(TE Length)/Genome Length", color = "TE Family") +
    theme(axis.title   = element_text(family = "serif", size = 12),
          axis.text.y  = element_text(family = "serif", size = 12),
          axis.text.x  = element_text(family = "serif", size = 12, 
                                      angle = 60, vjust = 1, hjust = 1),
          legend.text  = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          legend.position = c(0.8,0.7),
          plot.margin = margin(1,1,1,1, unit = "cm"))
  ggsave(plot = p_02, 
         filename = paste0(path, "Plot/TE-Erich.pdf"), width = 16, height = 16, units = "cm")
}
################################################################################
#
#
#
#
#
#
#
#
#
###########################7SK, Rpph1
#
#
#
#
#
#
#
#
#
#
#
#
################################################################################
#####整理出目标RNA的DNA Interaction的信息
{
  gf_1 <- rtracklayer::import(con = inx2, format = "gtf")
  gf_2 <- gf_1[which(gf_1$type=="gene")]
  tf_1 <- rtracklayer::import(con = inx3, format = "gtf")
  #####want <- "Rpph1"
  #####aim1 <- gf_1[which(gf_1$gene_name=="Rpph1" & gf_1$type=="exon")]
  want <- "7SK"
  aim1 <- tf_1[which(tf_1$gene_id == want)]
  #
  #
  #
  #
  for (i in c("n1","n2")) {
    print(i)
    rd_1 <- fread(file = paste0(path,"bed/",i,".join.bed"), sep = "\t", stringsAsFactors = F)
    rd_1 <- makeGRangesFromDataFrame(df = rd_1, 
                                     seqnames.field = "V2", start.field = "V3", 
                                     end.field = "V4", strand.field = "V6", keep.extra.columns = T)
    ov_1 <- suppressWarnings(findOverlaps(query = aim1, subject = rd_1, ignore.strand = T))
    rd_2 <- as.data.frame(rd_1[unique(subjectHits(ov_1))])
    colnames(rd_2) <- paste0("V", seq_along(colnames(rd_2)))
    
    wr01 <- rd_2
    write.table(x = wr01, 
                quote = F, col.names = F, row.names = F, sep = "\t",
                file = paste0(path, "specialOut_",want,"/",i,"_",want,"-is-RNA_DNA.bed"))
    wr02 <- data.frame(chr = rd_2$V8, start = rd_2$V9, 
                       end = rd_2$V10, name = ".", score = 1, strand = rd_2$V12)
    write.table(x = wr02, 
                quote = F, col.names = F, row.names = F, sep = "\t",
                file = paste0(path, "specialOut_",want,"/",i,"_",want,"-is-DNA.bed"))
    #
    #
    #
    #
    #
    #
    dtag <- makeGRangesFromDataFrame(df = unique(rd_2[,c(6,8,9,10,12)]), 
                                     seqnames.field = "V8", start.field = "V9",
                                     end.field = "V10",   strand.field = "V12")
    ##dtag <- dtag + 1000
    vv_1 <- suppressWarnings(findOverlaps(query = tf_1, subject = dtag, ignore.strand = T))
    tf_2 <- as.data.frame(tf_1[queryHits(vv_1)])
    wr03 <- tf_2
    write.table(x = wr03, sep = "\t", quote = F, row.names = F,
                file = paste0(path, "specialOut_",want,"/",i,"_",want,"-is-TE.DNA.Posi.bed"))
    
    
    tf_2$name <- paste0(tf_2$gene_id, ":", tf_2$family_id, ":", tf_2$class_id)
    tf_3 <- as.data.frame(table(tf_2$name))
    tf_3 <- tf_3[order(tf_3$Freq, decreasing = T),]
    wr04 <- tf_3
    write.table(x = wr04, sep = "\t", quote = F, row.names = F,
                file = paste0(path, "specialOut_",want,"/",i,"_",want,"-is-TE.DNA.Info.bed"))
    #
    #
    #
    #
    #
    #
    rd_d <- makeGRangesFromDataFrame(df = rd_2, keep.extra.columns = T, 
                                     seqnames.field = "V8",start.field = "V9", 
                                     end.field = "V10",  strand.field = "V12")
    vv_2 <- suppressWarnings(findOverlaps(query = gf_2, subject = rd_d, ignore.strand = T))
    df01 <- as.data.frame(gf_2[queryHits(vv_2)])
    df02 <- as.data.frame(rd_d[subjectHits(vv_2)])
    fina <- cbind(df02, df01)
    colnames(fina) <- paste0("V", seq_along(colnames(fina)))
    twos <- read.table(file = inx7, sep = "\t", stringsAsFactors = F)
    fina$Mark <- "not_2C_gene"
    fina$Mark[which(fina$V25 %in% twos$V4)] <- "2C_gene"
    wr05 <- data.frame(RNA_Posi = paste0("chr",fina$V6,":",fina$V7,"-",fina$V8), 
                       RNA_Strand = fina$V10, 
                       tag_ID = fina$V11,
                       DNA_Posi = paste0("chr",fina$V1,":",fina$V2,"-",fina$V3), 
                       DNA_Strand = fina$V5,
                       Gene_Name = fina$V25, Gene_Mark = fina$Mark,
                       Gene_Type = fina$V27, 
                       Gene_Posi = paste0("chr",fina$V14,":",fina$V15,"-",fina$V16))
    write.table(x = wr05, sep = "\t", quote = F, row.names = F,
                file = paste0(path, "specialOut_",want,"/",i,"_",want,"-is-Gene.Info.bed"))
  }
}
################################################################################
#####计算相应DNA TAGs在所有基因 ± 2 kb区域的数量
{
  #want <- "7SK"
  want <- "Rpph1"
  bed1 <- read.table(file = paste0(path, "specialOut_",want,"/n1_",want,"-is-DNA.bed"), sep = "\t")
  colnames(bed1) <- c("chr","start","end","name","score","strand")
  bed1 <- makeGRangesFromDataFrame(df = bed1, ignore.strand = F, keep.extra.columns = T)
  bed1Num <- length(bed1$name)
  bed2 <- read.table(file = paste0(path, "specialOut_",want,"/n2_",want,"-is-DNA.bed"), sep = "\t")
  colnames(bed2) <- c("chr","start","end","name","score","strand")
  bed2 <- makeGRangesFromDataFrame(df = bed2, ignore.strand = F, keep.extra.columns = T)
  bed2Num <- length(bed2$name)
  #
  #
  #
  #
  gf_1 <- rtracklayer::import(con = inx2, format = "gtf")
  gf_2 <- gf_1[which(gf_1$type=="gene")]
  gf_3 <- as.data.frame(gf_2 +2000)
  
  tst1 <- read.table(file = "/Reference/aaSHY/Ensembl99/GRCm38/fa/chrom.size", row.names = 1)
  seqs <- tst1$V2
  names(seqs) <- rownames(tst1)
  gf_3 <- suppressWarnings(makeGRangesFromDataFrame(df = gf_3, 
                                                    keep.extra.columns = T, seqinfo = seqs))
  gf_3 <- GenomicRanges::trim(gf_3, use.names = T)
  sums <- read.table(file = paste0(inx1,".fai"))
  sums <- sum(sums$V2)
  myta <- data.frame();wwww <- 0
  for (i in unique(gf_3$gene_name)) {#c("Duxf3","Setdb1","Zscan4c","Ssrp1","pr","Dot1l","Trim28")
    wwww <- wwww +1
    print(wwww)
    #i="Ptk6"
    gf_4 <- gf_3[which(gf_3$gene_name==i)]
    type <- paste0(unique(gf_4$gene_biotype), collapse = ",")
    bed4 <- bedtoolsr::bt.merge(gf_4)
    colnames(bed4) <- c("chr","start","end")
    gr_4 <- makeGRangesFromDataFrame(df = bed4, ignore.strand = T)
    wide <- sum(gr_4@ranges@width)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = gr_4, subject = bed1, ignore.strand = T))
    num1 <- length(unique(subjectHits(ov01)))
    exp1 <- round((wide/sums)*bed1Num, 6)
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = gr_4, subject = bed2, ignore.strand = T))
    num2 <- length(unique(subjectHits(ov02)))
    exp2 <- round((wide/sums)*bed2Num, 6)
    temp <- data.frame(name = i, biotype = type,
                       plus2kb = wide, 
                       n1Num = num1, n1NumExp = exp1, n2Num = num2, n2NumExp = exp2)
    myta <- rbind(myta, temp)
  }
  myta$rank1 <- myta$n1Num/myta$n1NumExp
  myta$rank2 <- myta$n2Num/myta$n2NumExp
  openxlsx::write.xlsx(x = myta, 
                       file = paste0(path, "specialOut_",want,"/",want,"-Num.Gene.xlsx"))
}
#####计算相应DNA TAGs在所有基因 ± 2 kb区域的数量, 作图
{
  #want <- "7SK"
  want <- "Rpph1"
  myta <- openxlsx::read.xlsx(paste0(path, "specialOut_",want,"/",want,"-Num.Gene.xlsx"))
  myta$uniNum <- rowSums(myta[,c(4,6)])
  myta$uniNumExp <- rowSums(myta[,c(5,7)]); myta$rank <- myta$uniNum/myta$uniNumExp
  {
    test <- myta[which(myta$biotype=="protein_coding"),]
    test$mark <- "other"
    test$mark[which(test$name %in% c("Duxf3","Nelfa","Zscan4d"))] <- "object"
    tst1 <- test[order(test$rank, decreasing = T),]
    tst1$colo <- NA
    tst1$colo[which(tst1$rank > 3)] <- "color-1"
    tst1$colo[which(tst1$rank <=3)] <- "color-2"
    #
    tst1$name <- factor(x = tst1$name, levels = tst1$name)
    p_01 <- ggplot(data = tst1) +
      geom_point(mapping = aes(x = name, y = log2(rank+0.1), color=colo)) +
      geom_point(data = tst1[which(tst1$mark=="object"),],mapping = aes(x = name, y = log2(rank+0.1))) +
      geom_hline(yintercept = log2(3+0.1), linetype = 2, color = "red", linewidth = 0.8) +
      theme_classic() +
      scale_x_discrete(expand = c(0.01,0)) +
      labs(y="Log2 (n1Num/n1NumExp +0.1)", x="Genes (ProteinCoding)", 
           title = paste0("GRID-Seq, ",want," RNA Interaction (Gene ± 2kb)")) +
      theme(axis.title = element_text(family = "serif", size = 18),
            plot.title = element_text(family = "serif", size = 20),
            axis.ticks.x = element_blank(),
            plot.margin = margin(unit = "cm",1,1,1,1),
            legend.position = c(0.7,0.7),
            legend.text = element_text(family = "serif", size = 16),
            legend.title = element_text(family = "serif", size = 16),
            axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 1),
            axis.line.y = element_line(linewidth = 1.2)) +
      scale_color_manual(values = c("steelblue","gray"), labels = c("rank >  3","rank <=3")) +
      guides(color = guide_legend(title = "Legend")) +
      geom_text_repel(data = tst1[which(tst1$mark =="object"),],
                      mapping = aes(x = name, y = log2(rank+0.1), label = name))
    p_01
    ggsave(plot = p_01, units = "cm", width = 18, height = 16,
           filename = paste0(path, "specialOut_",want,"/",want,"-Num.Gene-Score-v1.pdf"))
  }
}
################################################################################
#####画个点图, 被bound的TE, 按照RNA-DNA tags数量排序取top 10, 作图
{
  #want <- "7SK"
  want <- "Rpph1"
  
  
  pt01 <- read.table(paste0(path,"specialOut/repeatID.txt"), 
                     stringsAsFactors = F, sep = "\t", header = T)
  pt02 <- read.table(paste0(inx1,".fai"), 
                     stringsAsFactors = F, sep = "\t", header = F)
  refL <- sum(pt02$V2)
  pt01$frac <- pt01$length/refL
  pt03 <- read.table(paste0(path, "specialOut_",want,"/n1_",want,"-is-TE.DNA.Info.bed"), 
                     stringsAsFactors = F, sep = "\t", header = T)
  pt03 <- pt03[1:10,]
  pt03$name <- sapply(str_split(pt03$Var1, pattern = ":"),"[",1)
  pt03$frac <- NA
  for (i in seq_along(pt03$name)) {
    pt03$frac[i] <- pt01[which(pt01$gene_id == pt03$name[i]),"frac"]
  }
  pt04 <- read.table(paste0(path, "specialOut_",want,"/n2_",want,"-is-TE.DNA.Info.bed"),
                     stringsAsFactors = F, sep = "\t", header = T)
  pt03$erro <- NA; pt03$Freq2 <- NA
  for (i in seq_along(pt03$name)) {
    pt03$Freq2[i] <- pt04[which(pt04$Var1 == pt03$Var1[i]), "Freq"]
    pt03$erro[i] <- sd(c(pt03$Freq[i],pt03$Freq2[i]))/sqrt(2)
    pt03$FreqR[i] <- mean(c(pt03$Freq[i],pt03$Freq2[i]))
  }
  pt03$family <- sapply(str_split(pt03$Var1, pattern = ":"),"[",2)
  pt03$name <- factor(x = pt03$name, levels = pt03$name)
  print(pt03$name)
  p_02 <- ggplot(data = pt03) + 
    geom_errorbar(aes(x = name, y = Freq, ymin = Freq -erro, ymax = Freq +erro), alpha =0.5, width = 0.1) +
    geom_point(mapping = aes(x = name, y = Freq, size = frac, color = family)) +
    theme_classic() +
    scale_size_continuous(limits = c(0.001,0.03), breaks = c(0.001, 0.01, 0.02)) +
    scale_y_continuous(limits = c(0,1000), breaks = seq(0,1000,200)) +
    labs(x = "", y = paste0("Interaction number of ",want," with TE"), 
         title = paste0("GRID-Seq, ",want," RNA Interaction (On TE) (Top 10 Most enriched)"), 
         size = "Sum(TE Length)/Genome Length", color = "TE Family") +
    theme(axis.title   = element_text(family = "serif", size = 12),
          axis.text.y  = element_text(family = "serif", size = 12),
          axis.text.x  = element_text(family = "serif", size = 12, 
                                      angle = 60, vjust = 1, hjust = 1),
          legend.text  = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12),
          legend.position = c(0.8,0.7),
          plot.margin = margin(1,1,1,1, unit = "cm"))
  p_02
  ggsave(plot = p_02, width = 16, height = 16, units = "cm", 
         filename = paste0(path, "specialOut_",want,"/",want,"-TE.Enrich.pdf"))
}
################################################################################
#####画个点图, 被bound的TE, 按照RNA-DNA tags数量做Normalize
{
  #want <- "7SK"
  want <- "Rpph1"
  
  
  
  tf_1 <- rtracklayer::import(con = inx3, format = "gtf")
  name <- unique(paste0(tf_1$gene_id,":",tf_1$family_id,":",tf_1$class_id))
  gsiz <- read.table(file = paste0(inx1,".fai"), sep = "\t", stringsAsFactors = F)
  gsiz <- sum(gsiz$V2)
  red1 <- data.frame(tete = name, gsiz = gsiz)
  red1[,c("tsiz","n1","n2")] <- NA
  for (i in red1$tete) {
    print(i)
    tf_2 <- tf_1[which(tf_1$gene_id==str_split(i, pattern = ":")[[1]][1])]
    red1$tsiz[which(red1$tete==i)] <- sum(tf_2@ranges@width)
    for (N in c("n1","n2")) {
      ff_1 <- read.table(paste0(path, "specialOut_",want,"/",N,"_",want,"-is-TE.DNA.Info.bed"),
                         stringsAsFactors = F, sep = "\t", header = T)
      if (i %in% ff_1$Var1) {
        red1[which(red1$tete==i),N] <- ff_1[which(ff_1$Var1==i),"Freq"]
      }
      else {
        red1[which(red1$tete==i),N] <- 0
      }
    }
  }
  red1$nomz <- rowMeans(red1[,c("n1","n2")])/(red1$tsiz/red1$gsiz)
  wr_1 <- red1[order(red1$nomz, decreasing = T),]
  write.csv(x = wr_1, quote = F, row.names = F, 
            paste0(path, "specialOut_",want,"/",want,".TE.Norm.Nums.csv"))
}
################################################################################
#####画个点图, 被bound的TE, 按照RNA-DNA tags数量做Normalize, 作图
{
  want <- "7SK"
  #want <- "Rpph1"
  
  
  red1 <- read.csv(paste0(path, "specialOut_",want,"/",want,".TE.Norm.Nums.csv"))
  red1 <- red1[which(red1$nomz !=0),]
  red1$tete <- factor(red1$tete, levels = unique(red1$tete))
  red2 <- red1[grep("MM-int|MT_MmMT_Mm", red1$tete),]
  red2$labe <- sapply(str_split(red2$tete,":"),"[",1)
  
  ggplot() +
    geom_point(data = red1, aes(x = tete, y = log2(nomz)), color = "#d1c2d3") +
    geom_point(data = red2, aes(x = tete, y = log2(nomz)), color = "green") +
    geom_text_repel(data = red2, 
              aes(x = tete, y = log2(nomz), label = paste0(labe," (Position: ",rownames(red2),")"))) +
    scale_y_continuous(limits = c(floor(log2(min(red1$nomz))),
                                  ceiling(log2(max(red1$nomz))))) +
    scale_x_discrete(expand = c(0.02,1)) +
    theme_classic() +
    labs(x = "Sorted TE", y = "Log2 (Normalized Tags Number)", 
         title = paste0("GRID-Seq, ",want," RNA Interaction (On TE)")) +
    theme(axis.title   = element_text(family = "serif", size = 12),
          axis.text.y  = element_text(family = "serif", size = 12),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(1,1,1,1, unit = "cm"))
}
################################################################################
#####画个韦恩图, 被bound的基因是多少事2-cell基因
{
  want <- "Rpph1"
  red1 <- data.frame()
  for (i in c("n1","n2")) {
    tmp1 <- read.table(paste0(path, "specialOut_",want,"/",i,"_",want,"-is-Gene.Info.bed"),
                       sep = "\t", header = T)[,6:8]
    tmp1$samp <- i
    red1 <- rbind(red1,tmp1)
  }
  twos <- read.table(inx7, header = F)
  #
  #
  vens <- list()
  vens[["x"]] = twos$V4
  vens[["y"]] = unique(red1$Gene_Name[which(red1$samp=="n1")])
  vens[["z"]] = unique(red1$Gene_Name[which(red1$samp=="n2")])
  p_01 <- VennDiagram::venn.diagram(x = vens, filename = NULL, 
                                    main = "------------",
                                    sub = paste0("GRID-Seq -- ",want," Interact Info"), 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("2C_gene", "n1 Sample", "n2 Sample"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  #pdf(file = paste0(path, "Plot/venn.pdf"))
  grid.draw(p_01)
  #dev.off()
}
################################################################################
#####画个柱状图, 按tags number算, 被bound的top 10基因是哪些
{
  want <- "Rpph1"
  red1 <- data.frame()
  for (i in c("n1","n2")) {
    tmp1 <- read.table(paste0(path, "specialOut_",want,"/",i,"_",want,"-is-Gene.Info.bed"),
                       sep = "\t", header = T)[,6:8]
    tmp1$samp <- i
    red1 <- rbind(red1,tmp1)
  }
  red2 <- red1 %>%
    dplyr::group_by(Gene_Name, Gene_Mark, Gene_Type) %>% 
    dplyr::summarise(cou1 = n()/2, .groups = "drop")
  red2 <- red2[order(red2$cou1, decreasing = T),]
  #
  #
  #
  #
  twos <- grep("not_",red2$Gene_Mark,invert = T)
  red3 <- red2[c(1:10,if(length(twos) >= 10) {twos[1:10]} else {twos}),]
  red3$Gene_Name <- factor(red3$Gene_Name, levels = unique(red3$Gene_Name))
  red3$Gene_Mark <- factor(red3$Gene_Mark, levels = unique(red3$Gene_Mark))
  ggplot(data = red3) +
    geom_bar(aes(x = Gene_Name, y = cou1), stat = "identity") +
    facet_wrap(~Gene_Mark, scales = "free") +
    theme_classic() +
    labs(x = "", y = "tags number", 
         title = paste0("top 10 genes (by tags number)__bound by ",want)) +
    theme(plot.title   = element_text(family = "serif", size = 18, vjust = 1),
          axis.title   = element_text(family = "serif", size = 16),
          axis.text.y  = element_text(family = "serif", size = 16),
          axis.text.x  = element_text(family = "serif", size = 16, 
                                      angle = 60, vjust = 1, hjust = 1),
          plot.margin = margin(1,1,1,1, unit = "cm"),
          strip.background = element_blank(),
          strip.text = element_text(family = "serif", size = 16))
}













