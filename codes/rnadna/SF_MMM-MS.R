#Write by Sun Haiayng on 2023.09.18
###########################Environment和Index###################################
#01、搭环境
{
  rm(list=ls())
  for (i in c("data.table","stringr","ggforce","bedtoolsr","refGenome","yyplot","Biostrings","tidyr","pheatmap","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","Mfuzz","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","pathview","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、给索引
{
  path = "/ChIP_seq_1/aaSHY/pacbioIllumina/withMS/"
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
stag <- read.table(file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/ccStages.cntTable", sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
stag <- stag[grep(x = rownames(stag), pattern = ":", invert = T), -7]
colnames(stag) <- sapply(str_split(colnames(stag), ".bam."), "[", 2)
stag <- as.data.frame(apply(stag, 2, function(x){x/sum(x)*1000000}))#write.table(x = stag, file = "/ChIP_seq_1/aaSHY/pacbioIllumina/count/tetoolkit/stageExp.csv", quote = F, col.names = T, row.names = F, sep = ",")
#
sf01 <- read.table(file = "/Reference/aaSHY/BED/special/spliceF.txt", stringsAsFactors = F, sep = "\t", header = F)
ms01 <- read.table(file = "/ChIP_seq_2/aaSHY/pr/ask/MS/MS-data.txt", sep = "\t", header = T, stringsAsFactors = F)
inte <- intersect(sf01$V1, unique(ms01$Gene))
my01 <- stag[which(rownames(stag) %in% inte),]
my01$name <- rownames(my01)
my01 <- tidyr::gather(my01, key = "stage", value = "expres", -"name")
my02 <- ms01[which(ms01$Gene %in% inte),]
#
my01$stage <- factor(x = my01$stage, levels = unique(my01$stage))
p_01 <- ggplot(data = my01) + 
        geom_line(mapping = aes(x = stage, y = expres, group = name)) +
        geom_point(mapping = aes(x = stage, y = expres)) +
        labs(x="", y="CPM Expression",
             title = "CPM: 20 SF which occured in MM MS") +
        facet_wrap(.~ name, scales = "free_y") +
        theme(plot.margin = margin(1,1,1,1,unit = "cm"),
              axis.text.x = element_text(family = "serif", size = 12, angle = 60, vjust = 1, hjust = 1, colour = "black"),
              strip.text.x = element_text(size = 12, color = "red", face = "bold"))
p_01
View(my02)
##
##
##
#MM RNA pull down的蛋白, 有没有前面预测过的MM剪接位点的蛋白?
site <- c("Sox8","Sox9","Sox10","Tbx15","Tbx7","Tbx1","Tgif1","Pknox1","Tgif2lx","Arid3a",
          "Ybx1","Srsf3","Srsf5","Mbnl1","Khsrp","Hnrnpf","Hnrnp3","Hnrnp2","Hnrnp1",
          "Khdrbs1","Khdrbs2")

inte <- intersect()
stag[which(rownames(stag)=="Rplp0"),]

