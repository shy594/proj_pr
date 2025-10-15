#1、write by ~~~ on 2023.02.03
target = c("ENSMUST00000207370.1", "ENSMUST00000207659.1", "ENSMUST00000107843.10", "ENSMUST00000208829.1", "ENSMUST00000045325.13", "ENSMUST00000209124.1", "ENSMUST00000208897.1", "ENSMUST00000208778.1", "ENSMUST00000207735.1", "ENSMUST00000208312.1", "ENSMUST00000208938.1", "ENSMUST00000209056.2", "ENSMUST00000207522.1", "ENSMUST00000207702.1")

#2、使用featureCounts统计reads
{
bam <- list.files(path = "/ChIP_seq_1/sht/asJianJie/pacbio/illumina/bam_hisat2",pattern = "*bam$", full.names = T)
input <- paste0(bam,collapse = " ")
paste0("featureCounts -F 'GTF' -p -T 30 -g transcript_id -a ",
       "/Reference/zyw_reference/zyw_reference/gtf/mouse/gencode.vM21.primary_assembly.annotation.gtf -o ",
       "/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/featureCounts/stages_T.cnt ",
       input)
}

#3、处理featureCounts的结果
{
vta = read.table("/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/featureCounts/stages_T.cnt", 
                 header = TRUE, skip = 1, row.names = 1, check.names = F, stringsAsFactors = F)
vta = vta[,-c(1,2,3,4,5)]
colnames(vta) <- sapply(str_split(basename(colnames(vta)),pattern = "\\."), "[",1)
vta <- apply(vta,2,function(x){x/sum(x)*1000000})
data = data.frame("oocyte"=rowMeans(vta[,c(12,13)]),
                  "zygote"=rowMeans(vta[,c(14,15,16)]),
                  "x2cell"=rowMeans(vta[,c(1,2,3)]),
                  "x4cell"=rowMeans(vta[,c(4,5)]),
                  "x8cell"=rowMeans(vta[,c(6,7,8)]),
                  "blasto"=rowMeans(vta[,c(9,10)]))
}

#4、作图
{
for (i in target) {
#  i="ENSMUST00000207370.1"
  bpp=as.data.frame(t(data[i,]))
  bpp$qw = rownames(bpp)
  bpp$qw <- factor(bpp$qw, levels = bpp$qw)
  p_1 <- ggplot(data = bpp) +
    geom_bar(mapping = aes(x=qw,y=bpp[,1]),stat = "identity", fill="indianred4") +
    labs(title = i) + ylab(label ="expression level (CPM)") +
    geom_label(mapping = aes(x=qw,y=bpp[,1]),label=bpp[,1])
  p_1
  ggsave(plot = p_1, filename = paste0("/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/featureCounts/pr/fea_T_",
                                       i,".pdf"))
}
}





