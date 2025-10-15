#write by ~~~ at 2023.07.02
#1、搭环境
{
  rm(list=ls())
  for (i in c("homologene","Hmisc","aPEAR","GOplot","DESeq2","rtracklayer",
              "ggrepel","ReactomePA","stringr","tidyr","Biostrings",
              "IRanges","GenomicRanges","ggplot2","ggpubr","ggsci","ggthemes",
              "clusterProfiler","topGO","Rgraphviz","pathview","org.Mm.eg.db")) {library(i, character.only = T)};rm(i)
}
#2、给索引
{
  path = "/ChIP_seq_2/aaSHY/pr/ask/MS"
}
#3、筛选数据
{
  #mafi <- read.table(file = paste0(path, "/MS-data-2.txt"), header = T, stringsAsFactors = F, sep = "\t")
  #mafi[mafi=="-"] <- 0.5
  #mafi[,2:7] <- as.data.frame(apply(mafi[,2:7],2,as.numeric))
  #tmp1 <- dplyr::group_by(.data = mafi, by=Genes) %>% 
  #  dplyr::summarise(Peptides_Control=sum(Peptides_Control),
  #                   Peptides_MM_S=sum(Peptides_MM_S),
  #                   Peptides_MM_AS=sum(Peptides_MM_AS))
  #mafi <- as.data.frame(tmp1)
  #colnames(mafi)[1] <- "Gene"
  #mafi$Gene <- gsub(x = mafi$Gene, pattern = " ", replacement = "")
  #mafi$Gene <- toupper(mafi$Gene)
  #tmp2 <- homologene::homologene(genes = mafi$Gene, inTax = 9606, outTax = 10090)
  #colnames(tmp2) <- c("hsa","mmu","hsaID","mmuID")
  #mafi$geneName <- NA
  #for (i in seq_along(mafi$Gene)) {
  #  print(i)
  #  if (mafi$Gene[i] %in% tmp2$hsa) {
  #    mafi$geneName[i] <- tmp2[which(tmp2$hsa==mafi$Gene[i]),2][1]
  #  }
  #  else {
  #    mafi$geneName[i] <- Hmisc::capitalize(tolower(mafi$Gene[i]))
  #  }
  #}
  #mafi$theY <- log2((mafi$Peptides_MM_S+mafi$Peptides_MM_AS)/mafi$Peptides_Control +1)
  #mafi <- mafi[order(mafi$theY, decreasing = T),]
  #mafi$theX <- seq_along(mafi$theY)
}
######cpms <- read.table(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/count/TEtranscripts_geneName/ccpr.cntTable", header = T, row.names = 1)
######cpms <- cpms[,c(5,6)]
######cpms <- as.data.frame(apply(cpms, 2, function(x){x/sum(x)*1000000}))
######cpms$meanCPM <- rowMeans(cpms); cpms$name <- rownames(cpms)
######mafi$theY <- NA
######for (i in seq_along(mafi$Gene)) {
######  print(i)
######  if (mafi$Gene[i] %in% cpms$name) {
######    mafi$theY[i] <- cpms[which(cpms$name==mafi$Gene[i]),3]
######  }
######}
######mafi$theY <- log10(mafi$theY+1)
######纵轴
######mafi$mark[which(mafi$Peptides_Control==0 | mafi$Peptides_MM_S/mafi$Peptides_Control>=4)] <- "signif"
######mafi$mark[which(!(mafi$Peptides_Control==0 | mafi$Peptides_MM_S/mafi$Peptides_Control>=4))] <- "zNot_fuji"
######write.csv(x = mafi, file = paste0(path,"/MS-data.mark.csv"), row.names = F, quote = F)
######4、基于control和MM拉下的肽段数量的分布
mafi <- read.table(file = paste0(path, "/MS-data-3.txt"), header = T, stringsAsFactors = F, sep = "\t")
a1 <- c("pr","ff","Nop56","Nop58","Nop2","Gemin5"); labe <- mafi[which(mafi$Gene %in% a1),]
tmp1 <- mafi[which(mafi$defScoreLog>=1),]
tmp2 <- mafi[which(mafi$defScoreLog<1),]
p_1 <- ggplot() + 
       geom_point(data = tmp1, aes(x=theRank, y = defScoreLog), size=2, color="#5086C4") +
       geom_point(data = tmp2, aes(x=theRank, y = defScoreLog), size=2, color="gray") +
       geom_hline(aes(yintercept = log2(2)), 
                  linetype="dashed", linewidth=1,colour="red") +
       #geom_text_repel(data = labe, 
       #                aes(x= theRank, y = defScoreLog, label = Gene),
       #                color="indianred",
       #                force=20,size=4,point.padding = 0.5,hjust = -5,
       #                #arrow = arrow(length = unit(0.01, "npc"), type = "open", ends = "last"),
       #                segment.color="black",segment.size=0.5,segment.alpha=0.8,nudge_y=1) +
       labs(title = "MM RNA Pull Down", x="Rank of Positive Hits", y="log2(FC of S/(C +1))") +
       theme_classic() +
       theme(axis.title = element_text(family = "serif", size = 12, color="black"),
             legend.position.inside = c(0.8,0.2),
             axis.line = element_line(linewidth = 1.2),
             axis.text  = element_text(family = "serif", size = 12, color="black"),
             text = element_text(family = "serif", size = 12, color="black")) #+
       #scale_color_manual(values=c("red", "gray")) #+
       #scale_y_continuous(breaks = c(seq(0,30,by=5))) +
       #scale_x_continuous(limits = c(0,15), breaks = c(seq(0,15,by=1)))
p_1
ggsave(plot = p_1, filename = "/ChIP_seq_2/aaSHY/pr/ask/MS/MS-pulldown-noLAbel.r.pdf", units = "cm", width = 16, height = 14)
##
##
##
##
##
################################################################################
#富集分析
{
  mafi <- read.table(file = paste0(path, "/MS-data-3.txt"), header = T, stringsAsFactors = F, sep = "\t")
  names <- unique(mafi$Gene[which(mafi$defScoreLog>=1)])
  geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(names), keytype="SYMBOL", column="ENTREZID"))
  gobp <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, 
                                    keyType = "ENTREZID", ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
  gomf <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, 
                                    keyType = "ENTREZID", ont = "MF",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
  kegg <- clusterProfiler::enrichKEGG(gene = geyy, organism = "mmu", 
                                      keyType = "kegg", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  kegg <- DOSE::setReadable(kegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  mybp <- gobp@result[which(gobp@result$pvalue<0.05),]
  mymf <- gomf@result[which(gomf@result$pvalue<0.05),]
  mykg <- kegg@result[which(kegg@result$pvalue<0.05),]
  mykg$Description <- sapply(str_split(mykg$Description, pattern = " - Mus"), "[", 1)
  mybp$geneID <- gsub(x = mybp$geneID, pattern = "/", replacement = ",")
  mymf$geneID <- gsub(x = mymf$geneID, pattern = "/", replacement = ",")
  mykg$geneID <- gsub(x = mykg$geneID, pattern = "/", replacement = ",")
  openxlsx::write.xlsx(x = mybp, file = "/ChIP_seq_2/aaSHY/pr/ask/MS/signif.gobp.xlsx")
  openxlsx::write.xlsx(x = mymf, file = "/ChIP_seq_2/aaSHY/pr/ask/MS/signif.gomf.xlsx")
  openxlsx::write.xlsx(x = mykg, file = "/ChIP_seq_2/aaSHY/pr/ask/MS/signif.kegg.xlsx")
}
##
##
##
##
##
################################################################################
#上述KEGG结果的可视化
{
  mykg[4,]
  dta1 <- mykg[c(1,2,4,5,6,8,9,10,11,13),]
  #
  tmp1 <- aPEAR::findPathClusters(dta1,cluster = "hier",minClusterSize = 1)
  p_01 <- plotPathClusters(enrichment = dta1,sim = tmp1$similarity,
                           clusters = tmp1$clusters,fontSize = 4,colorBy = "pvalue",
                           nodeSize = "Count",drawEllipses = TRUE,outerCutoff = 0.01)+
    scale_color_gradientn(colours = c("#1A5592",'white','#B83D3D'),name = "pvalue")
  #
  dta1$Description <- factor(dta1$Description, levels = rev(dta1$Description))
  #p_02 <- ggplot(dta1)+
  #        geom_point(aes(x=-log10(pvalue), y=Description,color=-log10(pvalue),size=Count)) +
  #        theme_classic() +
  #        labs(y="",title = "KEGG: genes that score >= 1") +
  #        scale_color_gradient2(low = "white", mid = "red3", high = "red3", midpoint = 14, 
  #                              limits = c(2, 40), breaks = c(2,20,40)) +
  #        geom_segment(aes(x=0,xend = -log10(pvalue),y=Description,yend=Description), 
  #                     color="red3", linetype=2) +
  #        theme(text = element_text(size = 12,family = "serif",color = "black"),
  #              axis.title = element_text(size = 12,family = "serif",color = "black"),
  #              axis.text  = element_text(size = 12,family = "serif",color = "black"))
  #ggsave(plot = p_02,
  #       units = "cm", height = 16, width = 16, filename = "/ChIP_seq_2/aaSHY/pr/ask/MS/signif.kegg.pdf")
  p_02 <- ggplot(dta1)+
    geom_bar(aes(x=-log10(pvalue), y=Description), stat = "identity", fill="#1ba784", alpha=1.0) +
    theme_few() +
    geom_vline(aes(xintercept = 5), color = "white", linewidth  = 0.8) +
    geom_vline(aes(xintercept = 10), color = "white", linewidth = 0.8) +
    geom_vline(aes(xintercept = 15), color = "white", linewidth = 0.8) +
    geom_vline(aes(xintercept = 20), color = "white", linewidth = 0.8) +
    geom_vline(aes(xintercept = 25), color = "white", linewidth = 0.8) +
    geom_vline(aes(xintercept = 30), color = "white", linewidth = 0.8) +
    labs(y="",title = "KEGG: genes that score >= 1") +
    theme(text = element_text(size = 12,family = "serif",color = "black"),
          panel.border = element_rect(linewidth = 1),
          plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
          axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
          axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  p_02
  ggsave(plot = p_02,
         units = "cm", height = 10, width = 20, filename = "/ChIP_seq_2/aaSHY/pr/ask/MS/signif.kegg.r3.pdf")
}
##
##
##
##
##
################################################################################
#老师提供的 GO BP, 用DAVID在线网站分析的
{
  ff01 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ask/MS/MS.GO.BP.DAVID-selectedByLu.txt",
                     header = T, sep = "\t", stringsAsFactors = F)
  ff01 <- tidyr::separate(data = ff01, col = "Term", into = c("GOID","term"),sep = "~")
  ff01[1,]
  ff01$term <- factor(x = ff01$term, levels = rev(ff01$term))
  ff01$PValue
  p_03 <- ggplot(ff01)+
          geom_bar(aes(x=-log10(PValue), y=term), stat = "identity", fill="#1ba784", alpha=1.0) +
          theme_few() +
          geom_vline(aes(xintercept = 05), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 10), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 15), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 20), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 25), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 30), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 35), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 40), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 45), color = "white", linewidth = 0.8) +
          geom_vline(aes(xintercept = 50), color = "white", linewidth = 0.8) +
          labs(y="",title = "GO BP: genes that score >= 1") +
          theme(text = element_text(size = 12,family = "serif",color = "black"),
                panel.border = element_rect(linewidth = 1),
                plot.title = element_text(size = 12,family = "serif",color = "#1ba784",face = "bold"),
                axis.title = element_text(size = 12,family = "serif",color = "#1ba784"),
                axis.text  = element_text(size = 12,family = "serif",color = "#1ba784"))
  
  ggsave(plot = p_03,units = "cm", 
         height = 10, width = 20, filename = "/ChIP_seq_2/aaSHY/pr/ask/MS/MS.GO.BP.DAVID-selected.r1.pdf")
}
##
##
##
##
################################################################################
#和弦图
#5、突出GO term与ff等key Gene的联系
{
  a1a1 <- c("pr","ff","Nop56","Nop58","Nop2","Gemin5")
  mafi$mark <- NA
  mafi$mark[which(mafi$theY > 2.3)] <- "signif"
  mafi$mark[which(mafi$theY <=2.3)] <- "nosignif"
  names <- mafi$geneName[which(mafi$mark=="signif")]
  geyy <- na.omit(mapIds(x=org.Mm.eg.db, keys=unique(names), keytype="SYMBOL", column="ENTREZID"))
  gogo <- clusterProfiler::enrichGO(gene = geyy, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",pvalueCutoff = 0.05,qvalueCutoff = 0.05, readable = T)
  term <- gogo@result[1:6,]
  term$geneID <- gsub(x = term$geneID, pattern = "/", replacement = ",")
  term$category <- "Biological Process"
  term <- term[,c(1,2,10,6,8)]
  colnames(term)<-c('ID', 'term','category','adj_pval','genes')
  set1 <- term
  View(set1)
  #
  temp <- read.csv(file = "/ChIP_seq_2/aaSHY/pr/RNAseq/202104/DESeq2/pr.KO74.gene.csv", header = T, stringsAsFactors = F)
  colnames(temp) <- c("ID","baseMean","logFC","lfcSE","stat","pvalue","padj")
  set2 <- temp
  View(set2)
  #
  circ <- GOplot::circle_dat(set1, set2)
  circ$genes <- Hmisc::capitalize(tolower(circ$genes))
  circ <- circ[which(circ$genes %in% mafi[which(mafi$mark=="signif"),"geneName"]),]
  chrd <- GOplot::chord_dat(circ, genes = distinct(distinct(circ[,c("genes","logFC")])), process = unique(circ$term))
  GOplot::GOChord(chrd, gene.order = 'logFC')
  #
  chrd <- as.data.frame(chrd)
  dat1 <- mafi[which(mafi$mark=="signif"),]
  dat1$fcfc <- dat1$Peptides_MM_S/dat1$Peptides_Control
  for (i in rownames(chrd)){
    chrd[i,"logFC"] <- dat1[which(dat1$geneName==i),"fcfc"]
  }
  #chrd$logFC[is.infinite(chrd$logFC)] <- 10
  chrd$logFC <- scale(chrd$logFC)
  #chrd[which(rownames(chrd)=="ff"),]
  pt_1 <- GOplot::GOChord(chrd, 
                          gene.order = 'logFC',
                          gene.space = 0.25, gene.size = 4,
                          ribbon.col = brewer.pal(6, "Set3"))
  pdf(file = paste0(path,"/goplot.pdf"), width = 10, height = 11, family = "serif")
  pt_1
  dev.off()
}







