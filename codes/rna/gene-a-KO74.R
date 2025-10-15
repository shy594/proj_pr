#write by ~~~ at 2024.04.16
#01、搭环境"myenrichplot",
{
  for (i in c("data.table","scales","stringr","Biostrings","tidyr","pheatmap","enrichplot","enrichplot","DESeq2","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/"
  inx1 = "/Reference/aaSHY/Ensembl99/GRCm38/INDEX/STARv279a"
  inx2 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/BED/special/rRNA.bed"
  inx3 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/Mus_musculus.GRCm38.99.gtf"
  inx4 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf"
  inx5 = "/Reference/aaSHY/Ensembl99/GRCm38/gtf/GRCm38_Ensembl_rmsk_TE.gtf.locInd"
}
#######################################################################
#04、DESeq2：Gene、TE、TE sites
{
  shy_diff <- function(mark, geneFile, siteFile, Class, VS){
    tmp1 <- read.table(geneFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp1) <- sapply(str_split(basename(colnames(tmp1)),pattern = ".bam"), "[", 1)
    dfclas <- data.frame(row.names = colnames(tmp1), Class, stringsAsFactors = T)
    colLen <- ceiling(length(colnames(tmp1))/2)
    geneAA <- tmp1[grep(rownames(tmp1), pattern = ":", invert = T),]
    teOnly <- tmp1[grep(rownames(tmp1), pattern = ":", invert = F),]
    teGene <- tmp1
    #
    geneAA <- DESeq2::DESeqDataSetFromMatrix(countData = geneAA, colData = dfclas, design = ~Class)
    geneAA <- geneAA[rowSums(BiocGenerics::counts(geneAA) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(geneAA)) > 1,]
    geneAA <- DESeq2::DESeq(geneAA,)
    geneAA <- DESeq2::results(geneAA, contrast = VS)
    write.csv(geneAA, paste0(path,"DESeq2/",mark, ".geneOnly.csv"))
    upup <- subset(geneAA, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(geneAA, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".geneOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".geneOnly.down.csv"))
    #
    teOnly <- DESeq2::DESeqDataSetFromMatrix(countData = teOnly, colData = dfclas, design = ~Class)
    teOnly <- teOnly[rowSums(BiocGenerics::counts(teOnly) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teOnly)) > 1,]
    teOnly <- DESeq2::DESeq(teOnly,)
    teOnly <- DESeq2::results(teOnly, contrast = VS)
    write.csv(teOnly, paste0(path,"DESeq2/",mark, ".teOnly.csv"))
    upup <- subset(teOnly, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teOnly, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".teOnly.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".teOnly.down.csv"))
    #
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.upup.csv"))
    write.csv(down,paste0(path,"DESeq2/",mark, ".te.GeneBackGround.down.csv"))
    #
    #
    tmp2 <- read.table(siteFile, header = T, check.names = F, stringsAsFactors = F, row.names = 1)
    colnames(tmp2) <- sapply(str_split(basename(colnames(tmp2)),pattern = ".bam"), "[", 1)
    colLen <- ceiling(length(colnames(tmp2))/2)
    #
    teGene <- tmp2
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    teGene <- teGene[grep(rownames(teGene), pattern = ":", invert = F),]
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.GeneBackGround.down.csv"))
    #
    teGene <- tmp2[grep(rownames(tmp2), pattern = ":", invert = F),]
    teGene <- DESeq2::DESeqDataSetFromMatrix(countData = teGene, colData = dfclas, design = ~Class)
    teGene <- teGene[rowSums(BiocGenerics::counts(teGene) > 2) >= colLen, ]#[rowMeans(BiocGenerics::counts(teGene)) > 1,]
    teGene <- DESeq2::DESeq(teGene,)
    teGene <- DESeq2::results(teGene, contrast = VS)
    write.csv(teGene, paste0(path,"DESeq2/",mark, ".te.site.teOnly.csv"))
    upup <- subset(teGene, log2FoldChange > log2(1.5) & padj < 0.05)
    down <- subset(teGene, log2FoldChange < -log2(1.5)& padj < 0.05)
    write.csv(upup, paste0(path,"DESeq2/",mark, ".te.site.teOnly.upup.csv"))
    write.csv(down, paste0(path,"DESeq2/",mark, ".te.site.teOnly.down.csv"))
  }
  shy_diff(mark = "prKO", 
           geneFile = paste0(path, "count/pr.cntTable"),
           siteFile = paste0(path, "count/pr.site.cntTable"),
           Class = factor(c("KO","KO","WT","WT")), VS = c("Class","KO","WT"))
}
#######################################################################
#05、CPM Excel
{
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/count/pr.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("prKD-1","prKD-2","WT-1","WT-2")
  tmp1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = F),]#TE
  tmp2 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),]#Gene
  tmp3 <- cnt1#All
  tmp1 <- as.data.frame(apply(tmp1, 2, function(x){x/sum(x)*1000000}))
  tmp2 <- as.data.frame(apply(tmp2, 2, function(x){x/sum(x)*1000000}))
  tmp3 <- as.data.frame(apply(tmp3, 2, function(x){x/sum(x)*1000000}))
  tmp1$theName <- rownames(tmp1)
  tmp2$theName <- rownames(tmp2)
  tmp3$theName <- rownames(tmp3)
  openxlsx::write.xlsx(x = tmp1, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/count/prKO.cpm.onlyTE.xlsx")
  openxlsx::write.xlsx(x = tmp2, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/count/prKO.cpm.onlyGene.xlsx")
  openxlsx::write.xlsx(x = tmp3, file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/count/prKO.cpm.allGeneTE.xlsx")
}
#######################################################################
#06、pr KO和ff KD做上、下调的基因的交集
{
  {
    ff01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/DESeq2/prKO.geneOnly.csv", row.names = 1)
    ff01$mark <- NA
    ff01$mark[which(ff01$log2FoldChange > log2(1.5) & ff01$padj < 0.05)] <- "upup"
    ff01$mark[which(ff01$log2FoldChange < -log2(1.5) & ff01$padj < 0.05)] <- "down"
    #
    ff02 <- read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/202306/DESeq2/ff.geneOnly.csv", row.names = 1)
    ff02$mark <- NA
    ff02$mark[which(ff02$log2FoldChange > log2(1.5) & ff02$padj < 0.05)] <- "upup"
    ff02$mark[which(ff02$log2FoldChange < -log2(1.5) & ff02$padj < 0.05)] <- "down"
    #
    twoc <- read.table(file = "/Reference/aaSHY/BED/special/TWOC_gene.bed")
  }
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(ff01[which(ff01$mark=="upup"),])
  tmp1[["b"]] = rownames(ff02[which(ff02$mark=="upup"),])
  p_00 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "pr KO74 upup & ff KD upup", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr-UP", "ff-UP"), fill=c("#FD763F","#23BAC5"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT-ff/geneVenn.upup.r2.pdf"))
  grid.draw(p_00)
  dev.off()
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(ff01[which(ff01$mark=="upup"),])
  tmp1[["b"]] = rownames(ff02[which(ff02$mark=="upup"),])
  tmp1[["c"]] = twoc$V4
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "pr KO74 upup & ff KD upup", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr", "ff","twoC"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT-ff/geneVenn.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  #
  tmp2 <- list()
  tmp2[["a"]] = rownames(ff01[which(ff01$mark=="down"),])
  tmp2[["b"]] = rownames(ff02[which(ff02$mark=="down"),])
  tmp2[["c"]] = twoc$V4
  p_02 <- VennDiagram::venn.diagram(x = tmp2, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "pr KO74 down & ff KD down", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr", "ff","twoC"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT-ff/geneVenn.down.pdf"))
  grid.draw(p_02)
  dev.off()
  intersect(base::intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]])
  intersect(base::intersect(tmp2[["a"]],tmp2[["b"]]),tmp2[["c"]])
}
################################################################################
#07、pr---ff共同上调的566个基因--------25个2C-marker [pheatmap]
{
  pdf("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/PLOT-ff/geneVenn.upup.twoc.r2.pdf", width = 8, height = 10)
  gg01 <- intersect(intersect(tmp1[["a"]], tmp1[["b"]]),tmp1[["c"]])
  cnt1 <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/count/pr.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt1) <- c("prKO-1","prKO-2","WT-1","WT-2")
  hta1 <- cnt1[grep(x = rownames(cnt1), pattern=":", invert = T),c(3,4,1,2)]#onlyGene
  hta1 <- as.data.frame(apply(hta1, 2, function(x){x/sum(x)*1000000}))
  hta1 <- hta1[which(rownames(hta1) %in% gg01),]
  pheatmap(mat = hta1,
           main="\nexpression of selected genes -- pr Data",
           scale = "row", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = F, cluster_cols = F,
           row_names_justify = "left",
           labels_col = colnames(hta1), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  #
  cnt2 <- read.table("/ChIP_seq_2/aaSHY/pr/RNAseq/202306/count/ccff.cntTable", header = T, sep = "\t", row.names = 1)
  colnames(cnt2) <- c("ffKD-1","ffKD-2","WT-1","WT-2")
  hta2 <- cnt2[grep(x = rownames(cnt2), pattern=":", invert = T),c(3,4,1,2)]#onlyGene
  hta2 <- as.data.frame(apply(hta2, 2, function(x){x/sum(x)*1000000}))
  hta2 <- hta2[which(rownames(hta2) %in% gg01),]
  pheatmap(mat = hta2,
           main="\nexpression of selected genes -- ff Data",
           scale = "row", cellwidth = 20, 
           show_rownames = T, show_colnames = T, cluster_row = F, cluster_cols = F,
           row_names_justify = "left",
           labels_col = colnames(hta2), angle_col = "45",
           breaks = seq(-2,2,by=0.01), border = F,
           color  = c(colorRampPalette(colors = c("blue","white"))(length(seq(-2,2,by=0.01))/2),
                      colorRampPalette(colors = c("white","red"))(length(seq(-2,2,by=0.01))/2)))
  dev.off()
}
################################################################################
#07、GSEA: 在共同上调的基因上做
{
  ne01 <- base::intersect(tmp1[["a"]],tmp1[["b"]])
  res1 <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/DESeq2/prKD.geneOnly.csv"
  res1 <- "/ChIP_seq_2/aaSHY/pr/RNAseq/202306/DESeq2/ff.geneOnly.csv"
  res1 <- read.csv(file = res1, header = T, stringsAsFactors = F)[,c(1,3)]
  res2 <- suppressWarnings(clusterProfiler::bitr(geneID = res1$X, OrgDb = "org.Mm.eg.db",
                                                 fromType = "SYMBOL", toType = "ENTREZID", drop = T))
  res3 <- dplyr::distinct(.data = res2, SYMBOL, .keep_all = T)
  res3$rank <- apply(res3, 1, function(x){res1[which(res1$X==x[1]),2]})
  res4 <- res3[order(res3$rank, decreasing = T),c(2,3)]
  res5 <- res4$rank; names(res5) <- res4$ENTREZID
  #
  theG <- data.frame(term = "the_Genes", gene = ne01)
  theG <- suppressWarnings(clusterProfiler::bitr(theG$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db"))
  theG <- data.frame(ont = rep("the_Genes", length(unique(theG$ENTREZID))), 
                   gen = unique(theG$ENTREZID))
  gsea <- suppressWarnings(clusterProfiler::GSEA(geneList = res5, 
                                                 maxGSSize = 1000, TERM2GENE = theG, pvalueCutoff = 1))
  rm(res1,res2,res3,res4,res5,theG)
  #
  p_01 <- enrichplot::gseaplot2(x = gsea, geneSetID = 1, 
                                title = "GSEA for the 656 Genes", pvalue_table = T)
  pdf(file = "path", width = 12, height = 10)
  print(p_01)
  dev.off()
}
#######################################################################
#08、Venn图加上MM RNA互做的基因：ff和pr结合的基因也有MM RNA结合吗？
{
  ff01 <- read.table("/ChIP_seq_1/aaSHY/rnadna/GRIDm1/specialOut/n2.MM.DNA.bed")
  ff01 <- GenomicRanges::makeGRangesFromDataFrame(df = ff01, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
  gf01 <- import(inx4)
  gf01 <- gf01[which(gf01$type=="gene")]
  ov01 <- suppressWarnings(findOverlaps(query = gf01, subject = ff01, ignore.strand=T))
  theG <- gf01[queryHits(ov01)]$gene_name
  theG <- as.data.frame(table(theG))
  theG <- theG[order(theG$Freq, decreasing = T),]
  #
  ff01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/DESeq2/prKD.geneOnly.csv", row.names = 1)
  ff01$mark <- NA
  ff01$mark[which(ff01$log2FoldChange > log2(1.5) & ff01$pvalue < 0.05)] <- "upup"
  ff01$mark[which(ff01$log2FoldChange < -log2(1.5) & ff01$pvalue < 0.05)] <- "down"
  #
  ff02 <- read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/202306/DESeq2/ff.geneOnly.csv", row.names = 1)
  ff02$mark <- NA
  ff02$mark[which(ff02$log2FoldChange > log2(1.5) & ff02$pvalue < 0.05)] <- "upup"
  ff02$mark[which(ff02$log2FoldChange < -log2(1.5) & ff02$pvalue < 0.05)] <- "down"
  #
  gen1 <- as.character(theG$theG[which(theG$Freq>1)])
  gen1 <- theG$theG
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(ff01[which(ff01$mark=="upup"),])
  tmp1[["b"]] = rownames(ff02[which(ff02$mark=="upup"),])
  tmp1[["c"]] = gen1
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "pr KO74 upup & ff KD upup & MMTagGene2", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr", "ff","MMTagGene2"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT-ff/geneVenn.upup.pdf"))
  grid.draw(p_01)
  dev.off()
}
#######################################################################
#09、ff和pr共同结合的基因, 如果和MM KD/activation-的225个基因, overlap还剩多少?
{
  ge01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/MMKD/DESeq2/onlyGene.MM.csv", row.names = 1)#不带TE背景
  ge01$mark <- NA
  ge01$mark[which(ge01$log2FoldChange >  log2(1.5) & ge01$padj < 0.05)] <- "upup"
  ge01$mark[which(ge01$log2FoldChange < -log2(1.5) & ge01$padj < 0.05)] <- "down"
  #
  ge02 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/ask/MM/mt2Activation/DESeq2/rmSampAc2Ctrl1/MT2.geneOnly.csv", row.names = 1)
  ge02$mark <- NA
  ge02$mark[which(ge02$log2FoldChange >  log2(1.5) & ge02$padj < 0.05)] <- "upup"
  ge02$mark[which(ge02$log2FoldChange < -log2(1.5) & ge02$padj < 0.05)] <- "down"
  tmp4 <- list()
  tmp4[["a"]] = rownames(ge01[which(ge01$mark=="down"),])
  tmp4[["b"]] = rownames(ge02[which(ge02$mark=="upup"),])
  gen1 <- base::intersect(tmp4[["a"]],tmp4[["b"]])
  #
  ff01 <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/aaprKO74/DESeq2/prKD.geneOnly.csv", row.names = 1)
  ff01$mark <- NA
  ff01$mark[which(ff01$log2FoldChange > log2(1.5) & ff01$pvalue < 0.05)] <- "upup"
  ff01$mark[which(ff01$log2FoldChange < -log2(1.5) & ff01$pvalue < 0.05)] <- "down"
  #
  ff02 <- read.csv("/ChIP_seq_2/aaSHY/pr/RNAseq/202306/DESeq2/ff.geneOnly.csv", row.names = 1)
  ff02$mark <- NA
  ff02$mark[which(ff02$log2FoldChange > log2(1.5) & ff02$pvalue < 0.05)] <- "upup"
  ff02$mark[which(ff02$log2FoldChange < -log2(1.5) & ff02$pvalue < 0.05)] <- "down"
  #
  tmp1 <- list()
  tmp1[["a"]] = rownames(ff01[which(ff01$mark=="upup"),])
  tmp1[["b"]] = rownames(ff02[which(ff02$mark=="upup"),])
  tmp1[["c"]] = gen1
  p_01 <- VennDiagram::venn.diagram(x = tmp1, filename = NULL, 
                                    main = "gene Overlap",
                                    sub  = "pr KO74 upup & ff KD upup & kdUP-acDOWN", 
                                    main.fontfamily = "serif", sub.fontfamily = "serif", 
                                    cex = 1.5, main.cex = 2, sub.cex = 2,
                                    print.mode = "raw",#c("percent", "raw"),
                                    category.names = c("pr", "ff","kdUP-acDOWN"), fill=c("#FFFFCC","#CCFFFF","#FFCCCC"),
                                    units = "cm", height = 12, width = 12)#c('#FFFFCC','#CCFFFF',"#FFCCCC")#c('#FFFFCC','#CCFFFF',"#FFCCCC","#CCCCFF", "#CCFFCC")
  pdf(file = paste0(path, "PLOT-ff/geneVenn.upup.pdf"))
  grid.draw(p_01)
  dev.off()
  intersect(intersect(tmp1[["a"]],tmp1[["b"]]),tmp1[["c"]])
}
#######################################################################
#10、基因火山图
{
  shy_volcano_gene <- function(resFile, saveFile){
    resdata <- as.data.frame(read.csv(resFile, header = T, stringsAsFactors = F))
    resdata <- na.omit(resdata)
    resdata$threshold[resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05] = "Upregulated genes"
    resdata$threshold[resdata$log2FoldChange < -log2(1.5) & resdata$padj <0.05] = "Downregulated genes"
    resdata$threshold[!(resdata$log2FoldChange > log2(1.5) & resdata$padj <0.05) & !(resdata$log2FoldChange < -log2(1.5) & resdata$padj < 0.05)] = "Other genes"
    resdata$threshold <- factor(resdata$threshold, levels = c("Upregulated genes","Downregulated genes","Other genes"))
    p_1 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(padj), colour=threshold)) +
      labs(col="", x="Log2FoldChange",y="-Log10(adjusted p-value)") + 
      geom_point(size=0.8) +
      scale_x_continuous(limits = c(-6, 9), breaks = seq(-6, 9, 2)) +
      geom_hline(aes(yintercept=-log10(0.05)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=-log2(1.5)), linetype="dashed", colour="royalblue") +
      geom_vline(aes(xintercept=log2(1.5)), linetype="dashed", colour="royalblue") +
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      scale_color_manual(values = c("red","blue","gray"),
                         labels = c(paste0(names(table(resdata$threshold)[1])," (",unname(table(resdata$threshold)[1]),")"),
                                    paste0(names(table(resdata$threshold)[2])," (",unname(table(resdata$threshold)[2]),")"),
                                    paste0(names(table(resdata$threshold)[3])," (",unname(table(resdata$threshold)[3]),")"))) +
      theme(legend.position = c(0.2,0.9), 
            legend.background = element_rect(fill = NA),
            legend.spacing.x = unit(0,"mm"),
            legend.spacing.y = unit(0,"mm"),
            legend.key = element_blank(),
            legend.text = element_text(size = 13),
            axis.title = element_text(size = 13)) +
      guides(color = guide_legend(override.aes = list(size = 3)))
    ggsave(plot = p_1, saveFile, units = "cm", width = 18, height = 16)
  }
  shy_volcano_gene(resFile  = paste0(path, "DESeq2/prKD.geneOnly.csv"),
                   saveFile = paste0(path, "PLOT/prKD.geneOnly.volcano.pdf"))
}






