#1、write by ~~~ on 2023.02.03
target = c("pr", "pp8", "Ybx1", "Ybx3", "Nop2", "Angel2", "Fbl",
           "Nop58","Nop56","Syncrip","Hnrnpu","EWSR1","Mecp2","Dhx9",
           "Cnot1","Cnot2","Cnot3","Cnot4","Cnot7","Cnot8")
target = c("Oga","Abce1")
gtf <- import("/Reference/aaSHT/appTools/TEtoolkit/yesCHR/gencode.vM25.annotation.gtf")
gtf <- gtf[which(gtf$type=="gene"),]
gtf <- as.data.frame(gtf)
gtf[which(gtf$gene_name=="4930594M22Rik"),]

#2、使用featureCounts统计reads
bam <- list.files(path = "/ChIP_seq_1/sht/asJianJie/pacbio/illumina/bam_hisat2",pattern = "*bam$", full.names = T)
input <- paste0(bam,collapse = " ")
paste0("featureCounts -F 'GTF' -p -T 30 -g gene_name -a ",
       "/Reference/zyw_reference/zyw_reference/gtf/mouse/gencode.vM21.primary_assembly.annotation.gtf -o ",
       "/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/featureCounts/stages.cnt ",
       input)

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

for (i in target) {
  i="Abce1"
  bpp=as.data.frame(t(data[i,]))
  bpp$qw = rownames(bpp)
  bpp$qw <- factor(bpp$qw, levels = bpp$qw)
  p_1 <- ggplot(data = bpp) +
    geom_bar(mapping = aes(x=qw,y=bpp[,1]),stat = "identity", fill="indianred4") +
    labs(title = i) + ylab(label ="expression level (CPM)") +
    geom_label(mapping = aes(x=qw,y=bpp[,1]),label=bpp[,1])
  p_1
  ggsave(plot = p_1, filename = paste0("/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/featureCounts/pr/fea_",
                                       i,".pdf"))
}





#3、使用TEtranscripts统计reads

vta = read.table("/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/stages.cntTable", 
                 header = TRUE, row.names = 1, check.names = F, stringsAsFactors = F)
colnames(vta) <- sapply(str_split(basename(colnames(vta)),pattern = "\\."), "[",1)
vta <- apply(vta,2,function(x){x/sum(x)*1000000})
data = data.frame("oocyte"=rowMeans(vta[,c(12,13)]),
                  "zygote"=rowMeans(vta[,c(14,15,16)]),
                  "x2cell"=rowMeans(vta[,c(1,2,3)]),
                  "x4cell"=rowMeans(vta[,c(4,5)]),
                  "x8cell"=rowMeans(vta[,c(6,7,8)]),
                  "blasto"=rowMeans(vta[,c(9,10)]))
data[1,]


for (i in target) {
  i="Abce1"
  id = gtf[which(gtf$gene_name==i),"gene_id"]
  bpp=as.data.frame(t(data[id,]))
  bpp$qw = rownames(bpp)
  bpp$qw <- factor(bpp$qw, levels = bpp$qw)
  p_1 <- ggplot(data = bpp) +
    geom_bar(mapping = aes(x=qw,y=bpp[,1]),stat = "identity", fill="indianred4") +
    labs(title = i) + ylab(label ="expression level (CPM)") +
    geom_label(mapping = aes(x=qw,y=bpp[,1]),label=bpp[,1])
  p_1
  ggsave(plot = p_1, filename = paste0("/ChIP_seq_1/sht/asJianJie/pacbio/illumina/count/basedonHisat2/featureCounts/pr/pic_",
                                       i,".pdf"))
}










