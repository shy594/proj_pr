#Write by Sun Haiayng on 2024.04.19
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","ggforce","bedtoolsr","refGenome","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)#"yyplot",
}
#02、给索引
{
  prex <- c("n1","n2","n3")
  bath <- paste0(path, "specialOut/")
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/fa/GRCm38.dna.primary_assembly.fa"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/bwa/mm10"
  regs_01 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MM-int::ERVL::LTR.bed"
  regs_02 <- "/Reference/aaSHY/Ensembl99/GRCm38/BED/TEsubFamily/MT_Mm::ERVL::LTR.bed"
  ### inxA = "/Reference/aaSHY/zOther/TEconsensus/MuERVL/bwa/MuERVL"#bwa index -a is -p STR ref.fa
  ### inxB = "/Reference/aaSHY/zOther/TEconsensus/MT_Mm/bwa/MT2"#bwa index -a is -p STR ref.fa
  ### inxC = "/Reference/aaSHY/zOther/Rn45s/INDEX/bwa/Rn45s"#bwa index -a is -p STR ref.fa
}
################################################################################
#03、MM RNA tags bed file convert to BigWig file
#基于R [慢] [×]
{
  bed1 <- read.table(file = paste0(bath,"n1.MM.DNA.bed"), sep = "\t")
  colnames(bed1) <- c("chr","start","end","name","score","strand")
  bed2 <- dplyr::group_by(.data = bed1, chr, start, end, name, strand) %>% summarise(.groups = "keep", score = sum(score))
  bed2 <- as.data.frame(bed2)
  bed2 <- GenomicRanges::makeGRangesFromDataFrame(df = bed2, keep.extra.columns = T)
  bed3 <- bedtoolsr::bt.merge(bed2)
  bed3$V4 <- "."; bed3$V5 <- "."; bed3$V6 <- "."
  colnames(bed3) <- c("chr","start","end","name","score","strand")
  for (i in seq_along(bed3$chr)) {
    print(i)
    temp <- GenomicRanges::makeGRangesFromDataFrame(df = bed3[i,], keep.extra.columns = T)
    ov01 <- GenomicRanges::findOverlaps(query = temp, subject = bed2, ignore.strand=T)
    bed3$score[i] <- sum(bed2[unique(subjectHits(ov01))]$score)
  }
  bed3$score <- as.numeric(bed3$score)
  bed3 <- GenomicRanges::makeGRangesFromDataFrame(df = bed3, keep.extra.columns = T)
  tst1 <- read.table(file = "/Reference/aaSHY/Ensembl99/GRCm38/fa/chrom.size", row.names = 1)
  seqs <- tst1$V2
  names(seqs) <- rownames(tst1)
  bed3 <- GenomicRanges::makeGRangesFromDataFrame(df = bed3, keep.extra.columns = T, seqinfo = seqs)
  rtracklayer::export.bw(object = bed3, con = paste0(bath,"n1.MM.DNA.bw"),format="bigWig")
}
#基于Linux [快] [single] [×]
{
  size <- "/Reference/aaSHY/Ensembl99/GRCm38/fa/chrom.size"
  bath <- paste0(path, "specialOut/")
  shel <- paste0(path, "src/run_w.sh")
  cat("#!/bin/bash\n", file = shel)
  for (kk in c("n1","n2","n3")) {
    cmd_01 <- paste0("sort -k 1,1 -k 2,2n -k 3,3n ",bath,kk,".MM.DNA.bed ",
                     "> ",bath,kk,".MM.DNA.sort.bed","\n")
    cmd_02 <- paste0("bedtools genomecov -i ",bath,kk,".MM.DNA.sort.bed ",
                     "-g ",size," -bg >",bath,kk,"MM.DNA.sort.bedgraph","\n")
    cmd_03 <- paste0("bedGraphToBigWig ",bath,kk,"MM.DNA.sort.bedgraph ",size,
                     " ",bath,kk,".MM.DNA.sort.bw","\n")
    for (i in cmd_01) {cat(i, file = shel, append = T)}
    for (i in cmd_02) {cat(i, file = shel, append = T)}
    for (i in cmd_03) {cat(i, file = shel, append = T)}
  }
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
#基于Linux [快] [merge] [OK]
{
  size <- "/Reference/aaSHY/Ensembl99/GRCm38/fa/chrom.size"
  bath <- paste0(path, "specialOut/")
  shel <- paste0(path, "src/run_w.sh")
  cat("#!/bin/bash\n", file = shel)
  cmd_01 <- paste0("sort -k 1,1 -k 2,2n -k 3,3n ",
                   "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.bed ",
                   "> /ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bed","\n")
  cmd_02 <- paste0("bedtools genomecov -i ",
                   "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bed ",
                   "-g ",size," -bg >",
                   "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bedgraph","\n")
  cmd_03 <- paste0("bedGraphToBigWig ",
                   "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bedgraph ",
                   size," /ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bw","\n")
  for (i in cmd_01) {cat(i, file = shel, append = T)}
  for (i in cmd_02) {cat(i, file = shel, append = T)}
  for (i in cmd_03) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#04、MM RNA tags在2-cell基因上的富集
#[single] [×]
{
  shel <- paste0(path, "src/run_f.sh")
  cat("#!/bin/bash\n", file = shel)
  twoc <- "/Reference/aaSHY/BED/special/TWOC_gene.bed"
  bw01 <- paste0(bath,c("n1","n2","n3"),".MM.DNA.sort.bw", collapse = " ")
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", twoc, " -S ",bw01," -o ",path,
                   "PLOT/MM.twoc.mat.gz","\n")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", twoc, " -S ",bw01," -o ",path,
                   "PLOT/MM.twoc.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/MM.twoc.mat.gz --samplesLabel ", "n1 n2 n3",
                   " -out ",path,"PLOT/MM.twoc.mat.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
#[merge] [OK]
{
  shel <- paste0(path, "src/run_f.sh")
  cat("#!/bin/bash\n", file = shel)
  twoc <- "/Reference/aaSHY/BED/special/TWOC_gene.bed"
  bw01 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bw"
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", twoc, " -S ",bw01," -o ",path,
                   "PLOT/MM.twoc.mat.gz","\n")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", twoc, " -S ",bw01," -o ",path,
                   "PLOT/MM.twoc.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/MM.twoc.mat.gz --samplesLabel merge-n1n2 ",
                   "-out ",path,"PLOT/MM.twoc.mat.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#05、MM RNA tags在2-cell stage的MM Fusion基因上的富集
#[筛选出组装的融合转录本] [外显子水平的重叠] [OK] [被融合的基因 X]
{
  gf01 <- import(inx2); gf01 <- gf01[which(gf01$type=="gene")]
  tf01 <- import(inx3); tf01 <- tf01[which(tf01$gene_id %in% c("MM-int","MT_Mm"))]
  tf02 <- import("/ChIP_seq_1/aaSHY/pacbioIllumina/aaMerge/stringtie/star/two.gtf")
  #
  tf03 <- tf02[which(tf02$type == "transcript")]
  tf04 <- tf02[which(tf02$type == "exon")]
  ov01 <- suppressWarnings(findOverlaps(query = tf04, subject = tf01, ignore.strand = T))
  tf04 <- tf04[unique(queryHits(ov01))]
  tf03 <- tf03[which(tf03$transcript_id %in% unique(tf04$transcript_id))]
  #
  ov02 <- suppressWarnings(findOverlaps(query = gf01, subject = tf03, ignore.strand = T))
  ne02 <- as.data.frame(gf01[unique(queryHits(ov02))])
  wr02 <- data.frame(chr = ne02$seqnames, start = ne02$start, end = ne02$end, 
                     name = ne02$gene_name, score = ".", strand = ne02$strand)
  write.table(x = wr02, file = "/Reference/aaSHY/BED/special/MM.2cell.fusionGene.bed", 
              quote = F, sep = "\t", col.names = F, row.names = F)
  #
  ne02 <- as.data.frame(tf03)
  wr02 <- data.frame(chr = ne02$seqnames, start = ne02$start, end = ne02$end, 
                     name = ".", score = ".", strand = ne02$strand)
  write.table(x = wr02, file = "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed", 
              quote = F, sep = "\t", col.names = F, row.names = F)
  #
  ne03 <- tf02[which(tf02$transcript_id %in% tf03$transcript_id)]
  wr03 <- ne03
  rtracklayer::export(object = wr03, format = "gtf",
                      con = "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.gtf")
  #sed修改gene_id为gene_name
}
#[做富集] [merge] [OK]
{
  shel <- paste0(path, "src/run_f.sh")
  cat("#!/bin/bash\n", file = shel)
  twoc <- "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed"
  bw01 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/merge-n1n2.MM.DNA.sort.bw"
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", twoc, " -S ",bw01," -o ",path,
                   "PLOT/MM.twoc.mat.gz","\n")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 5000 -b 5000 ",
                   "--missingDataAsZero -R ", twoc, " -S ",bw01," -o ",path,
                   "PLOT/MM.twoc.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/MM.twoc.mat.gz --samplesLabel merge-n1n2 ",
                   "-out ",path,"PLOT/MM.twoc.mat.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#06、计算MM DNA tags在所有基因 ± 2 kb区域的数量
{
  bed1 <- read.table(file = paste0(bath,"n1.MM.DNA.bed"), sep = "\t")
  colnames(bed1) <- c("chr","start","end","name","score","strand")
  bed1 <- makeGRangesFromDataFrame(df = bed1, ignore.strand = F, keep.extra.columns = T)
  bed1Num <- length(bed1$name)
  bed2 <- read.table(file = paste0(bath,"n2.MM.DNA.bed"), sep = "\t")
  colnames(bed2) <- c("chr","start","end","name","score","strand")
  bed2 <- makeGRangesFromDataFrame(df = bed2, ignore.strand = F, keep.extra.columns = T)
  bed2Num <- length(bed2$name)
  #
  gf01 <- import(inx2)
  gf01 <- gf01[which(gf01$type=="gene")]
  gf02 <- gf01 +2000
  gf02 <- as.data.frame(gf02)
  tst1 <- read.table(file = "/Reference/aaSHY/Ensembl99/GRCm38/fa/chrom.size", row.names = 1)
  seqs <- tst1$V2
  names(seqs) <- rownames(tst1)
  gf02 <- suppressWarnings(GenomicRanges::makeGRangesFromDataFrame(df = gf02, keep.extra.columns = T, seqinfo = seqs))
  gf02 <- GenomicRanges::trim(gf02, use.names = T)
  sums <- read.table(file = paste0(inx1,".fai"))
  sums <- sum(sums$V2)
  myta <- data.frame();wwww <- 0
  for (i in unique(gf02$gene_name)) {#c("Duxf3","Setdb1","Zscan4c","Ssrp1","pprr","Dot1l","Trim28")
    wwww <- wwww +1
    print(wwww)
    #i="Ptk6"
    gf03 <- gf02[which(gf02$gene_name==i)]
    type <- paste0(unique(gf03$gene_biotype), collapse = ",")
    bed3 <- bedtoolsr::bt.merge(gf03)
    colnames(bed3) <- c("chr","start","end")
    gr03 <- makeGRangesFromDataFrame(df = bed3, ignore.strand = T)
    wide <- sum(gr03@ranges@width)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = gr03, subject = bed1, ignore.strand = T))
    num1 <- length(unique(subjectHits(ov01)))
    exp1 <- round((wide/sums)*bed1Num, 6)
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = gr03, subject = bed2, ignore.strand = T))
    num2 <- length(unique(subjectHits(ov02)))
    exp2 <- round((wide/sums)*bed2Num, 6)
    temp <- data.frame(name = i, biotype = type,
                       plus2kb = wide, 
                       n1Num = num1, n1NumExp = exp1, n2Num = num2, n2NumExp = exp2)
    myta <- rbind(myta, temp)
  }
  myta$rank1 <- myta$n1Num/myta$n1NumExp
  myta$rank2 <- myta$n2Num/myta$n2NumExp
  openxlsx::write.xlsx(x = myta, file = paste0(bath, "grid.MMNum.gene.xlsx"))
}
#
#
#
#
#下面的方案01、方案02, 分别是指"基于单个的n1、n2样本"、"合并n1、n2样本"
################################################################################
##########################*******方案01**********###############################
##########################*******方案01**********###############################
##########################*******方案01**********###############################
################################################################################
#07、MM DNA tags在所有基因±2 kb区域的数量的分布图
{
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta$rank1 <- myta$n1Num/myta$n1NumExp
  myta$rank2 <- myta$n2Num/myta$n2NumExp
  {
    test <- myta[which(myta$biotype=="protein_coding"),]
    test$mark <- "other"
    test$mark[which(test$name %in% c("Duxf3","Nelfa"))] <- "object"
    tst1 <- test[order(test$rank1, decreasing = T),]
    tst1$colo <- NA
    tst1$colo[which(tst1$rank1 > 3)] <- "color-1"
    tst1$colo[which(tst1$rank1 <=3)] <- "color-2"
    tst2 <- test[order(test$rank2, decreasing = T),]
    tst2$colo <- NA
    tst2$colo[which(tst2$rank2 > 3)] <- "color-1"
    tst2$colo[which(tst2$rank2 <=3)] <- "color-2"
    #
    tst1$name <- factor(x = tst1$name, levels = tst1$name)
    p_01 <- ggplot(data = tst1) +
      geom_point(mapping = aes(x = name, y = log2(rank1+1), color=colo)) +
      geom_point(data = tst1[which(tst1$mark=="object"),],mapping = aes(x = name, y = log2(rank1+1))) +
      geom_hline(yintercept = log2(3+1), linetype = 2, color = "red", linewidth = 0.8) +
      theme_classic() +
      scale_x_discrete(expand = c(0.01,0)) +
      labs(y="Log2 (n1Num/n1NumExp +1)", x="Genes (ProteinCoding)", title = "This Is Title...GRID-seq n1") +
      theme(axis.title = element_text(family = "serif", size = 12),
            plot.title = element_text(family = "serif", size = 15),
            axis.ticks.x = element_blank(),
            legend.position = c(0.7,0.7),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 1),
            axis.line.y = element_line(linewidth = 1.2)) +
      scale_color_manual(values = c("steelblue","gray"), labels = c("rank >  3","rank <=3")) +
      guides(color = guide_legend(title = "title")) +
      geom_text_repel(data = tst1[which(tst1$mark =="object"),],
                      mapping = aes(x = name, y = log2(rank1+1), label = name))
    p_01
    tst2$name <- factor(x = tst2$name, levels = tst2$name)
    p_02 <- ggplot(data = tst2) +
      geom_point(mapping = aes(x = name, y = log2(rank2+1), color=colo)) +
      geom_point(data = tst2[which(tst2$mark=="object"),],mapping = aes(x = name, y = log2(rank2+1))) +
      geom_hline(yintercept = log2(3+1), linetype = 2, color = "red", linewidth = 0.8) +
      theme_classic() +
      scale_x_discrete(expand = c(0.01,0)) +
      labs(y="Log2 (n2Num/n1NumExp +1)", x="Genes (ProteinCoding)", title = "This Is Title...GRID-seq n2") +
      theme(axis.title = element_text(family = "serif", size = 12),
            plot.title = element_text(family = "serif", size = 15),
            axis.ticks.x = element_blank(),
            legend.position = c(0.7,0.7),
            legend.text  = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 1),
            axis.line.y = element_line(linewidth = 1.2)) +
      scale_color_manual(values = c("steelblue","gray"), labels = c("rank >  3","rank <=3")) +
      guides(color = guide_legend(title = "title")) +
      geom_text_repel(data = tst2[which(tst2$mark =="object"),],
                      mapping = aes(x = name, y = log2(rank2+1), label = name))
    p_02
    plot <- p_01 +p_02 +plot_layout(ncol = 1)
    ggsave(plot = plot, filename = paste0(path, "specialOut/grid.MMNum.gene.protein.pdf"), 
           units = "cm", width = 18, height = 26)
  }
}
################################################################################
#08、n1、n2样本rank >3的基因的Venn图 [得到535个基因] [rank值是标准化的tags数量]
{
  tmp1 <- list()
  tmp1[["a"]] = tst1[which(tst1$rank1 >3),"name"]
  tmp1[["b"]] = tst2[which(tst2$rank2 >3),"name"]
  p_03 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "n1 rank >3 & n2 rank >3", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("n1", "n2"), fill=c("#FD763F","#23BAC5"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "specialOut/grid.MMNum.gene.protein.venn.pdf"))
  grid.draw(p_03)
  dev.off()
}
################################################################################
#09、pprr KO、ff KD中535个基因的表达热图
{
  ne01 <- base::intersect(tmp1[["a"]],tmp1[["b"]])
  pprr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/count/pprrKO.cpm.onlyGene.xlsx")
  pprr <- pprr[which(pprr$theName %in% ne01),]
  rownames(pprr) <- pprr$theName; pprr <- pprr[,1:4]
  ffa <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202306/count/ff.cpm.onlyGene.xlsx")
  ffa <- ffa[which(ffa$theName %in% ne01),]
  rownames(ffa) <- ffa$theName; ffa <- ffa[,1:4]
  cnt1 <- cbind(pprr, ffa)
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nthe heatmap of the venn 535 Genes",
           scale = "row", cellwidth = 30, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#10、MM DNA tags在535个基因上的富集 [弱富集]
{
  {
    ne01 <- base::intersect(tmp1[["a"]],tmp1[["b"]])
    wr01 <- as.data.frame(gf01[which(gf01$gene_name %in% ne01)])
    wr01 <- wr01[,c("seqnames","start","end","gene_name","score","strand")]
    wr01$score <- "."
    write.table(x = wr01, file = paste0(path, "specialOut/grid.MMNum.gene.protein.venn.bed"), 
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
  {
    shel <- paste0(path, "src/run_g.sh")
    cat("#!/bin/bash\n", file = shel)
    theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.bed")
    bw01 <- paste0(bath,c("n1","n2"),".MM.DNA.sort.bw", collapse = " ")
    cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 2500 -b 2500 ",
                     "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                     "PLOT/MM.theG.mat.gz","\n")
    cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                     "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                     "PLOT/MM.theG.mat.gz","\n")
    cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "PLOT/MM.theG.mat.gz --samplesLabel ", "n1 n2",
                     " -out ",path,"PLOT/MM.theG.mat.pdf","\n")
    for (i in cmd_09) {cat(i, file = shel, append = T)}
    for (i in cmd_12) {cat(i, file = shel, append = T)}
    print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
  }
}
################################################################################
#11、pprr、ff ChIP-seq在535个基因上的富集 [有富集]
{
  shel <- paste0(path, "src/run_h.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.bed")
  bw01 <- paste("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/pprr/single/bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/pprr/single/bw/bcp/IP-2_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/ff/single/bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/ff/single/bw/bcp/IP-2_Input.bw")
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                   "PLOT/pprrff.theG.mat.gz","\n")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                   "PLOT/pprrff.theG.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/pprrff.theG.mat.gz --samplesLabel ", 
                   "pprr-1 pprr-2 ff-1 ff-2 ",
                   "-out ",path,"PLOT/pprrff.theG.mat.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#12、H3K9me3，Setdb1, Kap1 ChIP-seq在535个基因上的富集 [仅H3K9me3无富集]
{
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.bed")
  bw03 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw")
  bw04 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("SRP094580 ","Setdb1-2 ","Trim28-1 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r3",c("-1","-2","-3"),".pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#13、H3K9me3在IAP上的富集 [有富集]
{
  inx0 <- c("/Reference/aaSHY/BED/TEsubFamily/IAPEy-int::ERVK::LTR.bed /Reference/aaSHY/BED/TEsubFamily/IAPEY_LTR::ERVK::LTR.bed")
  bw03 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_1_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_H3K9me3_2_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP051217_Control_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP069149_H3K9me3_pseudocount_1_log2.bw ",
                 "/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw")
  shell <- paste0(path,"dumpROOM/run_a.sh")
  cat("#!/bin/bash\n", file = shell)
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -p 20 -R ", inx0,
                   " -S ",paste(bw03),
                   " -o ",path,"dumpROOM/testH3K9me3.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",
                   path, "dumpROOM/testH3K9me3.mat.gz --samplesLabel ", 
                   paste0(c("SRP001533-1","SRP001533-2",
                            "SRP051217","SRP064188","SRP069149","SRP094580"),collapse = " "),
                   " -out ",path,"dumpROOM/testH3K9me3.mat.pdf","\n")
  for (i in cmd_09) {cat(i, append = T, file = shell)}
  for (i in cmd_12) {cat(i, append = T, file = shell)}
  print(paste0("nohup bash ",shell," >",path,"dumpROOM/run_a.log"," 2>&1 &"))
}
################################################################################
#14、H3K9me3 bamCoverage在535个基因上的富集 [不需要]
{
  shel <- paste0(path, "src/run_j.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.bed")
  bw01 <- paste0("/disk5/Mouse_ESC_ChIP/bamCoverage/SRP094580_histone_H3K9me3.bw")
  bw02 <- paste0("/disk5/Mouse_ESC_ChIP/bamCoverage/SRP094580_histone_Control.bw")
  bw03 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/bw/bcv/TRIM28.bcv.bw")
  bw04 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/bw/bcv/Input.bcv.bw")
  bw05 <- paste0("/disk5/Mouse_ESC_ChIP/bamCoverage/SRP064188_SETDB1.bw")
  bw06 <- paste0("/disk5/Mouse_ESC_ChIP/bamCoverage/SRP064188_Control.bw")
  bwbw <- c(bw01,bw02,bw03,bw04,bw05,bw06)
  for (e in c(1,3,5)) {
    labs <- c("SRP094580IP","SRP094580Input","Kap1IP","Kap1Input","Setdb1IP","Setdb1Input")
    labs <- labs[e:c(e+1)]
    cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                     "--missingDataAsZero -R ", theG, " -S ",
                     paste(bwbw[e:c(e+1)], collapse = " ")," -o ",path,
                     "PLOT/H3K9me3Setdb1Kap1.theG.mat.gz","\n")
    cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "PLOT/H3K9me3Setdb1Kap1.theG.mat.gz ","--samplesLabel ", 
                     paste(labs, collapse = " ")," -out ",
                     path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r5.",c("-1","-2","-3")[(e+1)/2],".pdf","\n")
    for (i in cmd_09) {cat(i, file = shel, append = T)}
    for (i in cmd_12) {cat(i, file = shel, append = T)}
  }
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
##########################*******方案02**********###############################
##########################*******方案02**********###############################
##########################*******方案02**********###############################
################################################################################
#15、MM DNA tags在所有基因±2 kb区域的数量的分布图, 并取rank >3, 得到1165个基因
{
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
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
      labs(y="Log2 (n1Num/n1NumExp +0.1)", x="Genes (ProteinCoding)", title = "This Is Title...GRID-seq uni") +
      theme(axis.title = element_text(family = "serif", size = 12),
            plot.title = element_text(family = "serif", size = 15),
            axis.ticks.x = element_blank(),
            legend.position = c(0.7,0.7),
            legend.text = element_text(family = "serif", size = 12),
            legend.title = element_text(family = "serif", size = 12),
            axis.text.x = element_blank(), axis.line.x = element_line(linewidth = 1),
            axis.line.y = element_line(linewidth = 1.2)) +
      scale_color_manual(values = c("steelblue","gray"), labels = c("rank >  3","rank <=3")) +
      guides(color = guide_legend(title = "title")) +
      geom_text_repel(data = tst1[which(tst1$mark =="object"),],
                      mapping = aes(x = name, y = log2(rank+0.1), label = name))
    p_01
    ggsave(plot = p_01, filename = paste0(path, "specialOut/grid.MMNum.gene.protein.uni.r2.pdf"), 
           units = "cm", width = 18, height = 16)
  }
}
################################################################################
#16、1165个基因的表达热图
#####在MM KD RNA-seq
{
  ne01 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed")
  ne01 <- ne01$V4
  mmmm <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/ask/MM/MMKD/count/MM-KD.gene.CPM.xlsx")
  mmmm <- mmmm[which(mmmm$geneName %in% ne01),]
  rownames(mmmm) <- mmmm$geneName; mmmm <- mmmm[,c(6,7,2,3)]
  cnt1 <- mmmm
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nthe heatmap of the venn 1165 Genes",
           scale = "row", cellwidth = 30, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#16、1254个基因的表达热图
#####在MM KD RNA-seq
{
  ne01 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed")
  ne01 <- ne01$V4
  mmmm <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/ask/MM/MMKD/count/MM-KD.gene.CPM.xlsx")
  mmmm <- mmmm[which(mmmm$geneName %in% ne01),]
  rownames(mmmm) <- mmmm$geneName; mmmm <- mmmm[,c(6,7,2,3)]
  cnt1 <- mmmm
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nMM KD____1254 genes",
           scale = "row", cellwidth = 30, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
#####在MT2__ KD RNA-seq
{
  ne01 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed")
  ne01 <- ne01$V4
  mmmm <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/ask/MM/MMKD/count/MM-KD.gene.CPM.xlsx")
  mmmm <- mmmm[which(mmmm$geneName %in% ne01),]
  rownames(mmmm) <- mmmm$geneName; mmmm <- mmmm[,c(6,7,4,5)]
  cnt1 <- mmmm
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nMT2 KD____1254 genes",
           scale = "row", cellwidth = 30, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
#####在MT2__ Ac RNA-seq
{
  ne01 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed")
  ne01 <- ne01$V4
  mmmm <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/ask/MM/MM-Activation/count/MM-Act.onlyGene.xlsx")
  mmmm <- mmmm[which(mmmm$theName %in% ne01),]
  rownames(mmmm) <- mmmm$theName; mmmm <- mmmm[,c(5,6,1,3)]
  cnt1 <- mmmm
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nMT2 Activation____1254 genes",
           scale = "row", cellwidth = 30,
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#16、pprr KO、ff KD中1165个基因的表达热图
{
  ne01 <- as.character(tst1[which(tst1$rank >3 & tst1$biotype=="protein_coding"), "name"])
  pprr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/count/pprrKO.cpm.onlyGene.xlsx")
  pprr <- pprr[which(pprr$theName %in% ne01),]
  rownames(pprr) <- pprr$theName; pprr <- pprr[,1:4]
  ffa <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202306/count/ff.cpm.onlyGene.xlsx")
  ffa <- ffa[which(ffa$theName %in% ne01),]
  rownames(ffa) <- ffa$theName; ffa <- ffa[,1:4]
  cnt1 <- cbind(pprr, ffa)
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  pheatmap(mat = cnt1,
           main="\nthe heatmap of the venn 1165 Genes",
           scale = "row", cellwidth = 30, 
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(cnt1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
}
################################################################################
#17、pprr KO、ff KD中1165个基因的表达热图 [设WT CPM为1]
{
  ne01 <- as.character(tst1[which(tst1$rank >3 & tst1$biotype=="protein_coding"), "name"])
  pprr <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/count/pprrKO.cpm.onlyGene.xlsx")
  pprr <- pprr[which(pprr$theName %in% ne01),]
  rownames(pprr) <- pprr$theName; pprr <- pprr[,1:4]
  ffa <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202306/count/ff.cpm.onlyGene.xlsx")
  ffa <- ffa[which(ffa$theName %in% ne01),]
  rownames(ffa) <- ffa$theName; ffa <- ffa[,1:4]
  cnt1 <- cbind(pprr, ffa)
  cnt1 <- cnt1[rowSums(cnt1) >0,]
  cnt2 <- cnt1
  cnt2 <- cnt2 +0.001
  cnt3 <- as.data.frame(apply(cnt2[,c(1:4)], 2, function(x){x/rowMeans(cnt2[,c(3,4)])}))
  cnt4 <- as.data.frame(apply(cnt2[,c(5:8)], 2, function(x){x/rowMeans(cnt2[,c(7,8)])}))
  dfa1 <- cbind(cnt3,cnt4)
  dfa1 <- log2(dfa1 +0.001)
  dfa1 <- dfa1[which(rowSums(is.na(dfa1)) <1),]
  #wr01 <- dfa1
  #wr01$geneName <- rownames(wr01)
  #write.csv(x = wr01, file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/normCPM-1165.csv", quote = F)
  dfa1[1,]
  pheatmap(mat = dfa1,
           main="\nthe heatmap of the venn 1165 Genes",
           scale = "none", cellwidth = 30,
           show_rownames = F, show_colnames = T, cluster_row = T, cluster_cols = F,
           labels_col = colnames(dfa1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  class(dfa2)
}
################################################################################
#18、pprr KO、ff KD中1165个基因里共同上调的基因, 和rank >5、rank >7的基因做交集 [CPM FC >0]
{
  red1 <- read.csv(file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/normCPM-1165.csv", row.names = 1)
  red2 <- red1
  red2 <- red2[which(rowMeans(red2[,c(1,2)])/rowMeans(red2[,c(3,4)]) >0),]
  red2 <- red2[which(rowMeans(red2[,c(5,6)])/rowMeans(red2[,c(7,8)]) >0),]
  #
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  ne01 <- as.character(myta$name)
  temp <- list()
  temp[["a"]] <- unique(red2$geneName)
  temp[["b"]] <- unique(ne01)
  #setdiff(temp[["a"]],temp[["b"]]); intersect(temp[["a"]],temp[["b"]])
  p_01 <- VennDiagram::venn.diagram(x = temp,
                                    main = "356 genes co-UP......genes rank >3",
                                    sub  = "", #a,b are samples of pprr, c,d are ff
                                    category.names = c("356 genes", "MM-tags rank >3"),
                                    fill=c("#FD763F","#23BAC5"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)
  dev.off();grid.draw(p_01)
}
################################################################################
#18-2、pprr KO、ff KD中1165个基因里共同上调的基因, 和rank >5、rank >7的基因做交集 [FC >1.5 & padj <.05]
{
  red1 <- read.csv(file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/DESeq2/pprrKO.geneOnly.csv", row.names = 1)
  red1 <- rownames(red1)[which(red1$log2FoldChange >log2(1.5) & red1$padj <0.05)]
  red2 <- read.csv(file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202306/DESeq2/ff.geneOnly.csv", row.names = 1)
  red2 <- rownames(red2)[which(red2$log2FoldChange >log2(1.5) & red2$padj <0.05)]
  #
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  ne01 <- as.character(myta$name)
  #
  temp <- list()
  temp[["a"]] <- unique(red1)
  temp[["b"]] <- unique(red2)
  temp[["c"]] <- unique(ne01)
  #setdiff(temp[["a"]],temp[["b"]]); intersect(temp[["a"]],temp[["b"]])
  p_01 <- VennDiagram::venn.diagram(x = temp,
                                    main = "pprr-Up...ff-Up...1165 rank >3",
                                    sub  = "", #a,b are samples of pprr, c,d are ff
                                    category.names = c("pprr-Up","ff-Up","MM-tags rank >3"),
                                    fill=c("#FD763F","#23BAC5","red"),
                                    filename = NULL, 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 1.5,
                                    print.mode = "raw",#c("percent", "raw"),
                                    units = "cm", height = 12, width = 12)
  dev.off();grid.draw(p_01)
}
################################################################################




################################################################################
#19、筛选出共同上调的基因, 用GSEA的方法展示是否富集在rank较大值的基因 [×]
{
  red1 <- read.csv(file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/normCPM-1165.csv", row.names = 1)
  red2 <- red1
  red2 <- red2[which(rowMeans(red2[,c(1,2)])/rowMeans(red2[,c(3,4)]) >0),]
  red2 <- red2[which(rowMeans(red2[,c(5,6)])/rowMeans(red2[,c(7,8)]) >0),]
  wr01 <- data.frame(theGenes = red2$geneName)
  wr01 <- as.data.frame(x = t(wr01))
  write.table(x = wr01, 
              file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/theG.gmt",
              sep = "\t", col.names = F, quote = F)
  #
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[order(myta$rank, decreasing = T),]
  myta$score <- as.numeric(scale(myta$rank))
  write.table(x = myta[,c("name","score")], 
              file = "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/tagsMM.rnk",
              sep = "\t", row.names = F, col.names = F, quote = F)
  #
  paste0("gseapy prerank -r ",
         "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/tagsMM.rnk ",
         "-g /ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/theG.gmt -p 10 -o ",
         "/ChIP_seq_2/aaSHY/pprr/RNAseq/202104/aapprrKO74/PLOT-ff/enrichGSEA")
}
################################################################################
#20、MM DNA tags在1165个基因上的富集
{
  {
    ne01 <- as.character(tst1[which(tst1$rank >3 & tst1$biotype=="protein_coding"), "name"])
    wr01 <- as.data.frame(gf01[which(gf01$gene_name %in% ne01)])
    wr01 <- wr01[,c("seqnames","start","end","gene_name","score","strand")]
    wr01$score <- "."
    write.table(x = wr01, file = paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed"), 
                quote = F, col.names = F, row.names = F, sep = "\t")
  }
  {
    shel <- paste0(path, "src/run_g.sh")
    cat("#!/bin/bash\n", file = shel)
    theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed")
    bw01 <- paste0(bath,c("n1","n2"),".MM.DNA.sort.bw", collapse = " ")
    cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 2500 -b 2500 ",
                     "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                     "PLOT/MM.theG.mat.gz","\n")
    cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                     "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                     "PLOT/MM.theG.mat.gz","\n")
    cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                     "PLOT/MM.theG.mat.gz --samplesLabel ", "n1 n2",
                     " -out ",path,"PLOT/MM.theG.mat.pdf","\n")
    for (i in cmd_09) {cat(i, file = shel, append = T)}
    for (i in cmd_12) {cat(i, file = shel, append = T)}
    print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
  }
}
################################################################################
#21、pprr、ff ChIP-seq在1165个基因上的富集 [×]
{
  shel <- paste0(path, "src/run_h.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed")
  bw01 <- paste("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/pprr/single/bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/pprr/single/bw/bcp/IP-2_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/ff/single/bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/ff/single/bw/bcp/IP-2_Input.bw")
  cmd_09 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                   "PLOT/pprrff.theG.mat.gz","\n")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                   "PLOT/pprrff.theG.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/pprrff.theG.mat.gz --samplesLabel ", 
                   "pprr-1 pprr-2 ff-1 ff-2 ",
                   "-out ",path,"PLOT/pprrff.theG.mat.r4.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
#22、pprr、ff ChIP-seq在1165个基因上的富集 [确定下来的] [挑选了样本]--------------------------------
{
  shel <- paste0(path, "src/run_h.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed")
  bw01 <- paste("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/pprr/single/bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/pprr/single/bw/bcp/IP-2_Input.bw")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                   "PLOT/pprrff.theG.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/pprrff.theG.mat.gz --samplesLabel ", 
                   "pprr-1 pprr-2 ",
                   "-out ",path,"PLOT/pprr.theG.mat.r1.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
  #
  #
  shel <- paste0(path, "src/run_h.sh")
  cat("#!/bin/bash\n", file = shel)
  bw01 <- paste("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/ff/single/bw/bcp/IP-1_Input.bw",
                "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240315/ff/single/bw/bcp/IP-2_Input.bw")
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bw01," -o ",path,
                   "PLOT/pprrff.theG.mat.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/pprrff.theG.mat.gz --samplesLabel ", 
                   "ff-1 ff-2 ",
                   "-out ",path,"PLOT/ff.theG.mat.r1.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#23、H3K9me3，Setdb1, Kap1, CTCF ChIP-seq在1165个基因上的富集 [前人跑出的bw]
{
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed")
  bw03 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw")
  bw04 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_12 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("SRP094580 ","Setdb1-2 ","Trim28-1 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r7",c("-1","-2","-3"),".pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
#23、H3K9me3，Setdb1, Kap1, CTCF ChIP-seq在1165个基因上的富集 [我的跑出的bw] [最后用到]---------------
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed")
  bw03 <- c("/ChIP_seq_2/aaSHY/pprr/publicData/H3K9me3/n2/bw/bcp/IP_Input-notSES.bw")
  bw04 <- c("/ChIP_seq_2/aaSHY/pprr/publicData/SETDB1/n2/bw/bcp/IP_Input-notSES.bw")
  #bw04 <- c("/disk5/Mouse_ESC_ChIP/bamCompare/SRP001533_SetDB1_pseudocount_1_log2.bw")
  bw05 <- c("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input-notSES.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2500 -b 2500 -m 1000 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("GSE90893_H3K9me3 ","GSE18371_SETDB1 ","GSE166041_TRIM28 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r8",c("-1","-2","-3"),".pdf","\n")
  for (i in cmd_01) {cat(i, file = shel, append = T)}
  for (i in cmd_02) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}
################################################################################
#23、CTCF ChIP-seq在1165个基因、MM tags上的富集
{
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed")
  theG <- paste0(path, "specialOut/mergeTwo/MM.DNA.sort.bed")
  bw06 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/CTCF/bw/bcp/IP_Input.bw")
  bwbw <- c(bw06)
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat-4.gz","\n")
  cmd_12 <- paste0("plotHeatmap --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat-4.gz",
                   " --samplesLabel CTCF ", 
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r3-4.pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >", paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
#24、H3K9me3，Setdb1, Kap1 ChIP-seq在MM fusion基因上的富集
{
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed"
  bw03 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP094580_histone_H3K9me3_pseudocount_1_log2.bw")
  bw04 <- paste0("/disk5/Mouse_ESC_ChIP/bamCompare/SRP064188_SETDB1_pseudocount_1_log2.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_09 <- paste0("computeMatrix scale-regions -p 20 -a 2500 -b 2500 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_12 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("SRP094580 ","Setdb1-2 ","Trim28-1 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.MMfus.theG.mat.r1",c("-1","-2","-3"),".pdf","\n")
  for (i in cmd_09) {cat(i, file = shel, append = T)}
  for (i in cmd_12) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log"," 2>&1 &"))
}
################################################################################
##########################*******方案03**********###############################
##########################*******方案03**********###############################
##########################*******方案03**********###############################
################################################################################
#方案03见代码"pprr_v2.R"的章节11, 具体步骤如下
#step 1: 注释ff、pprr的broad peak基因, 画二者的Venn
#step 2: 取二者的交集基因, 只保留Peak at Promoter的基因, 整理成gmt格式
#step 3: 根据前面计算的MM DNA tags normalized value (rank), 降序排序基因, rank取z-score, 整理成rnk格式
#step 4: 做GSEA [R clusterProfiler]
################################################################################
#25、MM DNA tags在所有基因±2 kb区域的数量的分布图, 并取rank >3, 得到1165个基因,
#####MM-KD、MM-Ac RNA-seq在1165个基因上的富集
{
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  ge01 <- import(con = inx2)
  ge01 <- ge01[which(ge01$type=="gene")]
  ge01 <- as.data.frame(ge01[which(ge01$gene_name %in% myta$name)])
  wr01 <- ge01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed", 
              sep = "\t", col.names = F, row.names = F, quote = F)
  bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
  #
  #关于如何合并重复样本的bw
  {
    ac01 <- list.files("/ChIP_seq_2/aaSHY/pprr/ask/MM/MM-Activation/bamCoverage", pattern = "bw", full.names = T)[c(1,3)]
    wt01 <- list.files("/ChIP_seq_2/aaSHY/pprr/ask/MM/MM-Activation/bamCoverage", pattern = "bw", full.names = T)[c(5,6)]
    cmd_01 <- paste0("bigWigMerge  bw file ", "bedGraph file","\n")
    cmd_02 <- paste0("sort -k1,1 -k2,2n ", "bedGraph file > sorted bedGraph file","\n")
    cmd_03 <- paste0("bedGraphToBigWig ","sorted bedGraph file ",
                     "/Reference/aaSHY/Ensembl99/GRCm38/fa/chrom.size ",
                     "bw file","\n")
  }
  #
  {
    bw01 <- list.files("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/dumpROOM/MMKDAC/KD", pattern = "bw", full.names = T)
    bw01 <- paste(bw01, collapse = " ")
    cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed1, " -S ",bw01," -o ", path,
                     "dumpROOM/", "the1165.gene.KD.mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/the1165.gene.KD.mat.gz --samplesLabel ",
                     "WT MM-KD ",
                     "--regionsLabel the1165.gene ",
                     "-out ",path,"dumpROOM/the1165.gene.KD.v1.pdf","\n")
    #
    #
    ####在1254基因上的图
    bed2 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
    cmd_03 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed2, " -S ",bw01," -o ", path,
                     "dumpROOM/", "the1254.gene.KD.mat.gz","\n")
    cmd_04 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/the1254.gene.KD.mat.gz --samplesLabel ",
                     "WT MM-KD ",
                     "--regionsLabel the1254.gene ",
                     "-out ",path,"dumpROOM/the1254.gene.KD.v1.pdf","\n")
  }
  #
  {
    bw02 <- list.files("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/dumpROOM/MMKDAC/AC", pattern = "bw", full.names = T)
    bw02 <- paste(bw02, collapse = " ")
    cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed1, " -S ",bw02," -o ", path,
                     "dumpROOM/", "the1165.gene.AC.mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/the1165.gene.AC.mat.gz --samplesLabel ",
                     "MT2-AC WT ",
                     "--regionsLabel the1165.gene ",
                     "-out ",path,"dumpROOM/the1165.gene.AC.v1.pdf","\n")
  }
  #
  {
    bw02 <- list.files("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/dumpROOM/pprrKO74", pattern = "bw", full.names = T)
    bw02 <- paste(bw02, collapse = " ")
    cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed1, " -S ",bw02," -o ", path,
                     "dumpROOM/", "the1165.gene.pprr.mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/the1165.gene.pprr.mat.gz --samplesLabel ",
                     "KO-74 WT ",
                     "--regionsLabel the1165.gene ",
                     "-out ",path,"dumpROOM/the1165.gene.pprr.v1.pdf","\n")
  }
  #
  {
    bw02 <- list.files("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/dumpROOM/ffKD", pattern = "bw", full.names = T)
    bw02 <- paste(bw02, collapse = " ")
    cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed1, " -S ",bw02," -o ", path,
                     "dumpROOM/", "the1165.gene.ff.mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/the1165.gene.ff.mat.gz --samplesLabel ",
                     "ff-KD WT ",
                     "--regionsLabel the1165.gene ",
                     "-out ",path,"dumpROOM/the1165.gene.ff.v1.pdf","\n")
  }
  #
  bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1165.gene.bed"
  {
    mm = "WT"
    bw02 <- list.files(paste0("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240124/toRUN/",mm), recursive = T, pattern = ".bw", full.names = T)
    bw02 <- bw02[grep(bw02, pattern="bcp")]
    bw02 <- paste(bw02, collapse = " ")
    cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 0 -b 0 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed1, " -S ",bw02," -o ", path,
                     "dumpROOM/", "the1165.gene.H3K9me3-",mm,".mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/","the1165.gene.H3K9me3-",mm,".mat.gz ",
                     "--samplesLabel ","IP-1 IP-3 ",
                     "--regionsLabel the1165.gene ",
                     "-out ",path,"dumpROOM/the1165.gene.H3K9me3-",mm,".v1.pdf","\n")
  }
  {
    bw02 <- list.files(paste0("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240124/toRUN/"), recursive = T, pattern = ".bw", full.names = T)
    bw02 <- bw02[grep(bw02, pattern="bcp")]
    bw02 <- paste(bw02, collapse = " ")
    cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 5000 -b 5000 -m 10000 ",
                     "--missingDataAsZero ", "-R ", bed1, " -S ",bw02," -o ", path,
                     "dumpROOM/", "the1165.gene.H3K9me3.mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --perGroup --heatmapHeight 12 --colorList white,red -m ",
                     path, "dumpROOM/","the1165.gene.H3K9me3.mat.gz ",
                     "--samplesLabel ","KO-IP-1 KO-IP-3 RN-IP-1 RN-IP-3 WT-IP-1 WT-IP-3 ",
                     "--regionsLabel the1165.gene ",
                     "-out ",path,"dumpROOM/the1165.gene.H3K9me3.v1.pdf","\n")
  }
}
################################################################################
#26、MM DNA tags在所有基因±2 kb区域的数量的分布图, 并取rank >3, 得到1165个基因, 
#26、MM KD AC Boxplot
{
  bath = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  #
  mt2c <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/ask/MM/MM-Activation/count/MM-Act.onlyGene.xlsx")
  mt2c <- mt2c[which(mt2c$theName %in% myta$name), c(7,1,3,5,6)]
  mt2c$MM-Act <- rowMeans(mt2c[,c(2,3)])
  mt2c$ctrl <- rowMeans(mt2c[,c(4,5)])
  mt2c <- mt2c[,c(1,6,7)]
  mt2c <- tidyr::gather(data = mt2c, key = "class", value = "cpm", -theName)
  ggplot(data = mt2c) +
    geom_boxplot(aes(x = class, y = log2(cpm +0.1))) +
    theme_classic() +
    #scale_y_continuous(limits = c(0, 35)) +
    labs(x="", y="Log2 (CPM +0.1)", title = "Expression of the 1165 genes") +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(family = "serif", size = 15),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          axis.text.y = element_text(family = "serif", size = 12),
          axis.title = element_text(family = "serif", size = 12)) +
    geom_signif(aes(x = class, y = log2(cpm +0.1)),
                comparisons = list(c("ctrl", "MM-Act")),
                test.args = c("less"),
                map_signif_level=T, textsize=6, test = wilcox.test, step_increase = 0.12)
  wilcox.test(x = mt2c$cpm[which(mt2c$class=="ctrl")],
              y = mt2c$cpm[which(mt2c$class=="MM-Act")], alternative = "less")
  #
  merv <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_2/aaSHY/pprr/ask/MM/MMKD/count/MMKD.onlyGene.xlsx")
  merv <- merv[which(merv$theName %in% myta$name), c(7,1,2,5,6)]
  merv$MMKD <- rowMeans(merv[,c(2,3)])
  merv$ctrl <- rowMeans(merv[,c(4,5)])
  merv <- merv[,c(1,6,7)]
  merv <- tidyr::gather(data = merv, key = "class", value = "cpm", -theName)
  ggplot(data = merv) +
    geom_boxplot(aes(x = class, y = log2(cpm +0.1))) +
    theme_classic() +
    #scale_y_continuous(limits = c(0, 35)) +
    labs(x="", y="Log2 (CPM +0.1)", title = "Expression of the 1165 genes") +
    theme(plot.margin = margin(1,1,1,1,unit = "cm"),
          plot.title = element_text(family = "serif", size = 15),
          axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
          axis.text.y = element_text(family = "serif", size = 12),
          axis.title = element_text(family = "serif", size = 12)) +
    geom_signif(aes(x = class, y = log2(cpm +0.1)),
                comparisons = list(c("ctrl", "MMKD")),
                test.args = c("great"),
                map_signif_level=T, textsize=6, test = wilcox.test, step_increase = 0.12)
  wilcox.test(x = merv$cpm[which(merv$class=="ctrl")],
              y = merv$cpm[which(merv$class=="MMKD")], alternative = "great")
}
################################################################################
################################################################################
#27、计算MM DNA tags在所有基因 ± 2 kb区域的数量 [###1170]
{
  bath <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/"
  bed1 <- read.table(file = paste0(bath,"n1.MM.DNA.bed"), sep = "\t")
  colnames(bed1) <- c("chr","start","end","name","score","strand")
  bed1 <- makeGRangesFromDataFrame(df = bed1, ignore.strand = F, keep.extra.columns = T)
  bed1Num <- length(bed1$name)
  bed2 <- read.table(file = paste0(bath,"n2.MM.DNA.bed"), sep = "\t")
  colnames(bed2) <- c("chr","start","end","name","score","strand")
  bed2 <- makeGRangesFromDataFrame(df = bed2, ignore.strand = F, keep.extra.columns = T)
  bed2Num <- length(bed2$name)
  sums <- read.table(file = paste0(inx1,".fai"))
  sums <- sum(sums$V2)
  #
  ge01 <- import(inx2)
  ge01 <- ge01[which(ge01$type=="gene")]
  gf01 <- as.data.frame(ge01)
  te01 <- import(inx3)
  te01 <- te01[grep(x = te01$gene_id, pattern = "MM-int|MT_Mm")] +1
  te02 <- bedtoolsr::bt.merge(i = te01)
  te02$name <- paste0("MM-", seq_along(te02$V1)); te02$strand <- "."
  te02 <- makeGRangesFromDataFrame(df = te02, keep.extra.columns = T,
                                   seqnames.field = "V1", start.field = "V2", end.field = "V3")
  tf02 <- as.data.frame(te02)
  fus1 <- read.table(file = "/Reference/aaSHY/BED/special/MM.2cell.fusionTranscript.bed", sep = "\t")
  fus1 <- makeGRangesFromDataFrame(df = fus1, seqnames.field = "V1", start.field = "V2",
                                   end.field = "V3", strand.field = "V6")
  fusf <- as.data.frame(fus1); fusf$fusID <- paste0("fus-", seq_along(fusf$seqnames))
  #为fusion gene注释gene
  ov01 <- findOverlaps(query = ge01, subject = fus1, ignore.strand = T)
  tmp1 <- cbind(fusf[subjectHits(ov01),c(1,2,3,6)], gf01[queryHits(ov01),c(1,2,3,12)])
  colnames(tmp1) <- c("fus_chr","fus_sta","fus_end","fus_ID","ge_chr","ge_sta","ge_end","ge_nam")
  dta1 <- dplyr::group_by(.data = tmp1, fus_chr,fus_sta,fus_end,fus_ID) %>% 
    dplyr::summarise(.groups = "keep",
                     ge_chr = paste0(ge_chr, collapse = ","),
                     ge_sta = paste0(ge_sta, collapse = ","),
                     ge_end = paste0(ge_end, collapse = ","),
                     ge_nam = paste0(ge_nam, collapse = ","))
  ov02 <- findOverlaps(query = te02, subject = fus1, ignore.strand = T)
  tmp2 <- cbind(fusf[subjectHits(ov02),c(1,2,3,6)], tf02[queryHits(ov02),c(1,2,3,6)])
  colnames(tmp2) <- c("fus_chr","fus_sta","fus_end","fus_ID","te_chr","te_sta","te_end","te_nam")
  dta2 <- dplyr::group_by(.data = tmp2, fus_chr,fus_sta,fus_end,fus_ID) %>% 
    dplyr::summarise(.groups = "keep",
                     te_chr = paste0(te_chr, collapse = ","),
                     te_sta = paste0(te_sta, collapse = ","),
                     te_end = paste0(te_end, collapse = ","),
                     te_nam = paste0(te_nam, collapse = ","))
  dta3 <- merge(x = dta1, y = dta2, by = "fus_ID")[,-c(9:11)]
  dta3[1,]
  dta3[4,]
  dta3[5,]
  dta4 <- as.data.frame(tidyr::separate_rows(data = dta3, ge_chr, ge_sta, ge_end, ge_nam, sep = ","))
  dta5 <- as.data.frame(tidyr::separate_rows(data = dta4, te_chr, te_sta, te_end, te_nam, sep = ","))
  dta5$maxPos <- NA
  for (i in seq_along(dta5$fus_ID)) {
    dta5$minPos[i] <- min(as.numeric(dta5[i,c("ge_sta","ge_end","te_sta","te_end")]))
    dta5$maxPos[i] <- max(as.numeric(dta5[i,c("ge_sta","ge_end","te_sta","te_end")]))
  }
  dta5$wide <- dta5$maxPos - dta5$minPos +1
  dta6 <- data.frame()
  for (i in unique(dta5$ge_nam)) {
    www1 <- dta5[which(dta5$ge_nam == i),-c(1:4)]
    www2 <- www1[which(www1$wide==min(www1$wide)),][1,]
    dta6 <- rbind(dta6, www2)
  }
  dta6[1,]
  #
  dta7 <- data.frame()
  nums <- 0
  for (i in seq_along(dta6$ge_nam)) {
    nums <- nums +1
    print(nums)
    hit1 <- dta6[i,]
    hit1 <- makeGRangesFromDataFrame(df = hit1, ignore.strand = T,
                                     seqnames.field = "ge_chr", start.field = "minPos",
                                     end.field = "maxPos", keep.extra.columns = T)
    wide <- sum(hit1@ranges@width)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = hit1, subject = bed1, ignore.strand = T))
    num1 <- length(unique(subjectHits(ov01)))
    exp1 <- round((wide/sums)*bed1Num, 6)
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = hit1, subject = bed2, ignore.strand = T))
    num2 <- length(unique(subjectHits(ov02)))
    exp2 <- round((wide/sums)*bed2Num, 6)
    type <- ge01$gene_biotype[which(ge01$gene_name==dta6$ge_nam[i])][1]
    temp <- data.frame(name = dta6$ge_nam[i], 
                       biotype = type,
                       wideL = wide, 
                       n1Num = num1, n1NumExp = exp1, n2Num = num2, n2NumExp = exp2)
    dta7 <- rbind(dta7, temp)
  }
  dta7$rank1 <- dta7$n1Num/dta7$n1NumExp
  dta7$rank2 <- dta7$n2Num/dta7$n2NumExp
  dta7$rank <- rowMeans(dta7[,c("rank1","rank2")])
  #
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta$mark <- NA
  for (i in seq_along(dta7$name)) {
    myta[which(myta$name == dta7$name[i]),"rank"] <- dta7[which(dta7$name == dta7$name[i]),"rank"]
    myta$mark[which(myta$name == dta7$name[i])] <- "fusg"
  }
  myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  myta <- myta[order(myta$rank, decreasing = T),]
  wr01 <- as.data.frame(ge01[which(ge01$gene_name %in% myta$name)])
  wr01 <- wr01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1170.gene.bed", 
              sep = "\t", col.names = F, row.names = F, quote = F)
}
#28、计算MM DNA tags在所有基因 ± 2 kb区域的数量 [###1254]
{
  wraa <- dta5
  lian <- as.data.frame(ge01[which(ge01$gene_name %in% unique(dta5$ge_nam))])
  wraa$lian <- NA
  for (i in seq_along(wraa$fus_ID)) {
    wraa$lian[i] <- as.character(lian$strand[which(lian$gene_name==wraa$ge_nam[i])][1])
  }
  wraa$fusPos <- paste0(wraa$fus_chr.x, ":", wraa$fus_sta.x,"-",wraa$fus_end.x)
  wraa$genePos <- paste0(wraa$ge_chr,":",wraa$ge_sta,"-",wraa$ge_end)
  wraa$MMPos <- paste0(wraa$te_chr,":",wraa$te_sta,"-",wraa$te_end)
  wraa <- wraa[,c("fus_ID","fusPos","ge_nam","lian","genePos","te_nam","MMPos")]
  wraa[1,]
  openxlsx::write.xlsx(x = wraa, 
                       file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/fus-gene-MM.xlsx")
  #吕老师返回了/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/MM-fus.xlsx
  eta1 <- openxlsx::read.xlsx(xlsxFile = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/MM-fus.xlsx")
  eta2 <- dplyr::group_by(.data = eta1, ge_nam) %>% 
    dplyr::summarise(.groups = "keep",
                     chr = paste0(chr, collapse = ","), start = paste0(start, collapse = ","),
                     end = paste0(end, collapse = ","))
  eta2$seqs <- NA; eta2$minPos <- NA; eta2$maxPos <- NA
  for (i in seq_along(eta2$ge_nam)) {
    eta2$seqs[i] <- str_split(eta2$chr[i], pattern = ",")[[1]][1]
    eta2$minPos[i] <- min(as.numeric(str_split(eta2$start[i], pattern = ",")[[1]]))
    eta2$maxPos[i] <- max(as.numeric(str_split(eta2$end[i], pattern = ",")[[1]]))
  }
  eta2 <- as.data.frame(eta2[,c(1,5,6,7)])
  openxlsx::write.xlsx(x = eta2, 
                       file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/MM-fus-2.xlsx")
  #
  eta3 <- data.frame()
  nums <- 0
  for (i in seq_along(eta2$ge_nam)) {
    nums <- nums +1
    print(nums)
    hit1 <- eta2[i,]
    hit1 <- makeGRangesFromDataFrame(df = hit1, ignore.strand = T,
                                     seqnames.field = "seqs", start.field = "minPos",
                                     end.field = "maxPos", keep.extra.columns = T)
    wide <- sum(hit1@ranges@width)
    ov01 <- suppressWarnings(GenomicRanges::findOverlaps(query = hit1, subject = bed1, ignore.strand = T))
    num1 <- length(unique(subjectHits(ov01)))
    exp1 <- round((wide/sums)*bed1Num, 6)
    ov02 <- suppressWarnings(GenomicRanges::findOverlaps(query = hit1, subject = bed2, ignore.strand = T))
    num2 <- length(unique(subjectHits(ov02)))
    exp2 <- round((wide/sums)*bed2Num, 6)
    type <- ge01$gene_biotype[which(ge01$gene_name==eta2$ge_nam[i])][1]
    temp <- data.frame(name = eta2$ge_nam[i], 
                       biotype = type,
                       wideL = wide, 
                       n1Num = num1, n1NumExp = exp1, n2Num = num2, n2NumExp = exp2)
    eta3 <- rbind(eta3, temp)
  }
  eta3$rank1 <- eta3$n1Num/eta3$n1NumExp
  eta3$rank2 <- eta3$n2Num/eta3$n2NumExp
  eta3$rank <- rowMeans(eta3[,c("rank1","rank2")])
  myta <- openxlsx::read.xlsx(xlsxFile = paste0(bath, "grid.MMNum.gene.xlsx"))
  myta$mark <- NA
  for (i in seq_along(eta3$name)) {
    myta[which(myta$name == eta3$name[i]),c(4,5,6,7,8,9,12)] <- eta3[which(eta3$name == eta3$name[i]),c(4,5,6,7,8,9,10)]
    myta[which(myta$name == eta3$name[i]),c(3,10,11)] <- "fusg"
    myta$mark[which(myta$name == eta3$name[i])] <- "fusg"
  }
  #openxlsx::write.xlsx(x = myta, file = paste0(bath, "grid.MMNum.gene-2.xlsx"))
  #
  #myta <- myta[which(myta$biotype=="protein_coding" & myta$rank >3),]
  eee1 <- myta[which(myta$name %in% c("Duxf3","AW822073")),]
  myta <- myta[which(myta$rank >2 & myta$n1Num >4 & myta$n2Num >4),]
  myta <- rbind(myta, eee1)
  myta <- myta[order(myta$rank, decreasing = T),]
  wr01 <- as.data.frame(ge01[which(ge01$gene_name %in% myta$name)])
  wr01 <- wr01[,c(1,2,3,12,4,5)]
  wr01$width <- "."
  write.table(x = wr01, file = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed", 
              sep = "\t", col.names = F, row.names = F, quote = F)
}
{
  path = "/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240124/toRUN/"
  bw01 <- list.files(path = path, pattern = "notSES.bw", full.names = T, recursive = T)
  bw01 <- bw01[grep(x = bw01, pattern = "bcp")][c(1,5)]
  {
    bed1 <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
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
                       "-out ",path,"zempROOM/K9me3-IP-VS-",basename(inx0),".r11-2.pdf","\n")
      for (i in cmd_01) {cat(i, append = T, file = shell)}
      for (i in cmd_02) {cat(i, append = T, file = shell)}
    }
    print(paste0("nohup bash ",shell," >",path,"zempROOM/run_n.log"," 2>&1 &"))
  }
}
################################################################################
#29、H3K9me3, Setdb1, Kap1的ChIP-seq在1254个基因上的富集 [对上面的23的修改]
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
  bw03 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/H3K9me3/n1/bw/bcp/IP_Input-notSES.bw")
  bw04 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/SETDB1/n1/bw/bcp/IP_Input-notSES.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input-notSES.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2500 -b 2500 -m 1000 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("GSE73432-H3K9me3 ","GSE73432-SETDB1 ","GSE166041-TRIM28 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r6",c("-1","-2","-3"),".pdf","\n")
  for (i in cmd_01) {cat(i, file = shel, append = T)}
  for (i in cmd_02) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
  bw03 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/H3K9me3/n2/bw/bcp/IP_Input-notSES.bw")
  bw04 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/SETDB1/n1/bw/bcp/IP_Input-notSES.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input-notSES.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2500 -b 2500 -m 1000 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("SRP094580-H3K9me3 ","GSE73432-SETDB1 ","GSE166041-TRIM28 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r6",c("-1-2","-2","-3"),".pdf","\n")
  for (i in cmd_01) {cat(i, file = shel, append = T)}
  for (i in cmd_02) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  theG <- "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/PLOT/the1254.gene.bed"
  bw03 <- paste0("/ChIP_seq_2/aaSHY/pprr/ChIPseq/20240124/toRUN/WT/bw/bcp/IP-1_Input-notSES.bw")
  bw04 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/SETDB1/n1/bw/bcp/IP_Input-notSES.bw")
  bw05 <- paste0("/ChIP_seq_2/aaSHY/pprr/publicData/TRIM28/n1/bw/bcp/IP_Input-notSES.bw")
  bwbw <- c(bw03,bw04,bw05)
  cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2500 -b 2500 -m 1000 ",
                   "--missingDataAsZero -R ", theG, " -S ",bwbw," -o ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz","\n")
  cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 --colorList white,red -m ",path,
                   "PLOT/H3K9me3Setdb1Kap1.theG.mat",c("-1","-2","-3"),".gz",
                   " --samplesLabel ", 
                   c("WT-H3K9me3 ","GSE73432-SETDB1 ","GSE166041-TRIM28 "),
                   "-out ",path,"PLOT/H3K9me3Setdb1Kap1.theG.mat.r6",c("-1-3","-2","-3"),".pdf","\n")
  for (i in cmd_01) {cat(i, file = shel, append = T)}
  for (i in cmd_02) {cat(i, file = shel, append = T)}
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
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
#------------------------------拒稿/审稿阶段
#
#
#
#
#
#
#
#
################################################################################
#01、画下表达激活的标记, 比如H3K9ac, H3K4me3, H3K27ac
#### 这是实验室之前的人跑出的bw
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  beds = c(paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed"),
           regs_01,regs_02,
           paste0(path, "specialOut/mergeTwo/MM.DNA.sort.bed"))
  beds_nads = c("MM_target_1165","MM","MT_Mm","MM_target")
  
  bwbw = c("/disk5/zyw/fromChIPseq2/zyw/E14_histone/bamcompare/H3K4me3_VS_Control__H3K4me1_2_3_H3K36me3_SES.bw",
           "/disk5/zyw/fromChIPseq2/zyw/E14_histone/bamcompare/H3K9ac_VS_Control__H3K14ac_H3K9ac_SES.bw",
           "/disk5/zyw/fromChIPseq2/zyw/E14_histone/bamcompare/H3K27ac_VS_Control__H3K27ac_H2BK20ac_SES.bw")
  bwbw_nads = c("H3K4me3","H3K9ac","H3K27ac")
  
  
  for(G in seq_along(beds_nads)) {
    if (beds_nads[G] == "MM_target_1165") {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2500 -b 2500 -m 1000 ",
                       "--missingDataAsZero -R ", beds[G], " -S ",bwbw," -o ",path,
                       "PLOT/tmps-",seq(length(bwbw)),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 ",
                       "--colorList white,red -m ",path,
                       "PLOT/tmps-",seq(length(bwbw)),".mat.gz ",
                       "--samplesLabel ",bwbw_nads," ",
                       "-out ",path,
                       "PLOT/",beds_nads[G],"_by_",bwbw_nads,".pdf","\n")
    } else {
      cmd_01 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000  ",
                       "--missingDataAsZero -R ", beds[G], " -S ",bwbw," -o ",path,
                       "PLOT/tmps-",seq(length(bwbw)),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 ",
                       "--colorList white,red -m ",path,
                       "PLOT/tmps-",seq(length(bwbw)),".mat.gz ",
                       "--samplesLabel ",bwbw_nads," ",
                       "-out ",path,
                       "PLOT/",beds_nads[G],"_by_",bwbw_nads,".pdf","\n")
    }
    for (i in cmd_01) {cat(i, file = shel, append = T)}
    for (i in cmd_02) {cat(i, file = shel, append = T)}
  }
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}
#02、画下表达激活的标记, 比如H3K9ac, H3K4me3, H3K27ac
#### 这是我跑出的bw, GSE61874
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_j.sh")
  cat("#!/bin/bash\n", file = shel)
  beds = c(paste0(path, "specialOut/grid.MMNum.gene.protein.venn.uni.bed"),
           regs_01,regs_02,
           paste0(path, "specialOut/mergeTwo/MM.DNA.sort.bed"))
  beds_nads = c("MM_target_1165","MM","MT_Mm","MM_target")
  
  bwbw = list.files(paste0("/sc/aaSHY/epiGenetic/mesc_markers/sets-01/bw/bcp"),full.names = T)
  bwbw_nads = gsub("_Input-SES.bw","",basename(bwbw))
  
  for(G in seq_along(beds_nads)) {
    if (beds_nads[G] == "MM_target_1165") {
      cmd_01 <- paste0("computeMatrix scale-regions -p 10 -a 2500 -b 2500 -m 1000 ",
                       "--missingDataAsZero -R ", beds[G], " -S ",bwbw," -o ",path,
                       "PLOT/j-tmps-",seq(length(bwbw)),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 ",
                       "--colorList white,red -m ",path,
                       "PLOT/j-tmps-",seq(length(bwbw)),".mat.gz ",
                       "--samplesLabel ",bwbw_nads," ",
                       "-out ",path,
                       "PLOT/",beds_nads[G],"_by_",bwbw_nads,".pdf","\n")
    } else {
      cmd_01 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000  ",
                       "--missingDataAsZero -R ", beds[G], " -S ",bwbw," -o ",path,
                       "PLOT/j-tmps-",seq(length(bwbw)),".mat.gz","\n")
      cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 ",
                       "--colorList white,red -m ",path,
                       "PLOT/j-tmps-",seq(length(bwbw)),".mat.gz ",
                       "--samplesLabel ",bwbw_nads," ",
                       "-out ",path,
                       "PLOT/",beds_nads[G],"_by_",bwbw_nads,".pdf","\n")
    }
    for (i in cmd_01) {cat(i, file = shel, append = T)}
    for (i in cmd_02) {cat(i, file = shel, append = T)}
  }
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}
################################################################################
#03、将MM target sites分为是否在重复序列的情况
{
  site_01 <- read.table(paste0(path, "specialOut/mergeTwo/MM.DNA.sort.bed"), header = F)
  colnames(site_01) <- c("chr","start","end","name","score","strand")
  site_02 <- makeGRangesFromDataFrame(df = site_01, keep.extra.columns = T)
  reps_01 <- import(inx3)
  reps_02 <- reps_01#[which(reps_01$class_id == "Satellite")]
  ovov <- suppressWarnings(findOverlaps(query = site_02, subject = reps_02, ignore.strand =T))
  site_02_reps <- as.data.frame(site_02[unique(queryHits(ovov))])[,c(1:3,6,7,5)]
  nots <- setdiff(as.numeric(rownames(site_01)), unique(queryHits(ovov)))
  site_02_nots <- as.data.frame(site_02[nots])[,c(1:3,6,7,5)]
  write.table(x = site_02_reps,
              file = paste0(path, "specialOut/mergeTwo/MM-targets-reps.bed"),
              col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(x = site_02_nots,
              file = paste0(path, "specialOut/mergeTwo/MM-targets-nots.bed"),
              col.names = F, row.names = F, quote = F, sep = "\t")
}
#---------~~~---#@-----
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_i.sh")
  cat("#!/bin/bash\n", file = shel)
  beds = c(paste0(path, "specialOut/mergeTwo/MM-targets-reps.bed"),
           paste0(path, "specialOut/mergeTwo/MM-targets-nots.bed"))
  beds_nads = c("MM_YesReps","MM_NotReps")
  
  bwbw = c("/disk5/zyw/fromChIPseq2/zyw/E14_histone/bamcompare/H3K4me3_VS_Control__H3K4me1_2_3_H3K36me3_SES.bw",
           "/disk5/zyw/fromChIPseq2/zyw/E14_histone/bamcompare/H3K9ac_VS_Control__H3K14ac_H3K9ac_SES.bw",
           "/disk5/zyw/fromChIPseq2/zyw/E14_histone/bamcompare/H3K27ac_VS_Control__H3K27ac_H2BK20ac_SES.bw")
  bwbw_nads = c("H3K4me3","H3K9ac","H3K27ac")
  for(G in seq_along(beds_nads)) {
    cmd_01 <- paste0("computeMatrix reference-point -p 20 ",
                     "--referencePoint center -a 5000 -b 5000 ",
                     "--missingDataAsZero -R ", beds[G], " -S ",bwbw," -o ",path,
                     "PLOT/activeMarker/i-isReps/tmps-",seq(length(bwbw)),".mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 ",
                     "--colorList white,red -m ",path,
                     "PLOT/activeMarker/i-isReps/tmps-",seq(length(bwbw)),".mat.gz ",
                     "--samplesLabel ",bwbw_nads," ",
                     "-out ",path,
                     "PLOT/activeMarker/i-isReps/",beds_nads[G],"_by_",bwbw_nads,".pdf","\n")
    for (i in cmd_01) {cat(i, file = shel, append = T)}
    for (i in cmd_02) {cat(i, file = shel, append = T)}
  }
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}
{
  path = "/ChIP_seq_1/aaSHY/rnadna/GRIDm1/"
  shel <- paste0(path, "src/run_j-2.sh")
  cat("#!/bin/bash\n", file = shel)
  beds = c(paste0(path, "specialOut/mergeTwo/MM-targets-reps.bed"),
           paste0(path, "specialOut/mergeTwo/MM-targets-nots.bed"))
  beds_nads = c("MM_YesReps","MM_NotReps")
  
  bwbw = list.files(paste0("/sc/aaSHY/epiGenetic/mesc_markers/sets-01/bw/bcp"),full.names = T)
  bwbw_nads = gsub("_Input-SES.bw","",basename(bwbw))
  for(G in seq_along(beds_nads)) {
    cmd_01 <- paste0("computeMatrix reference-point -p 20 --referencePoint center -a 5000 -b 5000  ",
                     "--missingDataAsZero -R ", beds[G], " -S ",bwbw," -o ",path,
                     "PLOT/activeMarker/j-isReps/tmps-",seq(length(bwbw)),".mat.gz","\n")
    cmd_02 <- paste0("plotHeatmap --zMin -0.5 --zMax 0.5 --heatmapHeight 12 ",
                     "--colorList white,red -m ",path,
                     "PLOT/activeMarker/j-isReps/tmps-",seq(length(bwbw)),".mat.gz ",
                     "--samplesLabel ",bwbw_nads," ",
                     "-out ",path,
                     "PLOT/activeMarker/j-isReps/",beds_nads[G],"_by_",bwbw_nads,".pdf","\n")
    for (i in cmd_01) {cat(i, file = shel, append = T)}
    for (i in cmd_02) {cat(i, file = shel, append = T)}
  }
  print(paste0("nohup bash ",shel, " >",paste0(path,"Log/"),basename(shel),".log", " 2>&1 &"))
}























