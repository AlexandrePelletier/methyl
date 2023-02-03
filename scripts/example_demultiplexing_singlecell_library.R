#D1#
out<-"outputs/droso"
dir.create(out)

library(Seurat)
library(data.table)


cr_dir<-"~/RUN/Run_737_single-cell/ouput/multi/D1/outs/multi/count/raw_feature_bc_matrix/"

umis<-Read10X(cr_dir)

#filter cells with enough CMO and mRNA
#mRNA
n_mRNAs<-colSums(umis$`Gene Expression`)
summary(n_mRNAs)
n_umis_lib<-data.table(bc=colnames(umis$`Gene Expression`),
           nCount_RNA=colSums(umis$`Gene Expression`),
           nCount_CMO=colSums(umis$`Multiplexing Capture`))
n_umis_libf<-n_umis_lib[nCount_RNA>0&nCount_CMO>0]
n_umis_libf #160k
n_umis_libf[,bc_rna_rank:=rank(-nCount_RNA)]
ggplot(n_umis_libf[order(bc_rna_rank)])+
  geom_line(aes(x=bc_rna_rank,y=nCount_RNA))+geom_hline(yintercept =500 )+
  scale_x_log10()+scale_y_log10()

ggplot(n_umis_libf[order(bc_rna_rank)])+
  geom_boxplot(aes(x=nCount_RNA>500,y=nCount_CMO))+scale_y_log10()

n_umis_libfr<-n_umis_libf[nCount_RNA>500]
n_umis_libfr#11k
n_umis_libfr[,bc_cmo_rank:=rank(-nCount_CMO)]

ggplot(n_umis_libfr[order(bc_cmo_rank)])+
  geom_line(aes(x=bc_cmo_rank,y=nCount_CMO))+
  geom_hline(yintercept =1600)+
  scale_x_log10()+scale_y_log10()
n_umis_libfrc<-n_umis_libfr[nCount_CMO>1600]

n_umis_libfrc #9.7k

lapply(umis,dim)
#demultiplexing
d1<-CreateSeuratObject(counts = umis$`Gene Expression`[,n_umis_libfrc$bc])


rowSums(umis$`Multiplexing Capture`[,n_umis_libfrc$bc])#have all CMO
#need keep only CM01-4
cmos<-c("CMO301","CMO302","CMO303")
d1[["CMO"]] <- CreateAssayObject(counts = umis$`Multiplexing Capture`[cmos,n_umis_libfrc$bc])

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
d1 <- NormalizeData(d1, assay = "CMO", normalization.method = "CLR")
d1 <- HTODemux(d1, assay = "CMO",positive.quantile = 0.95)

head(d1@meta.data)

table(d1$CMO_classification.global)
 # Doublet Negative  Singlet 
 #    1042       57     8619 


table(d1$hash.ID)


saveRDS(d1,fp(out,"d1.rds"))

d1<-subset(d1,CMO_classification.global=="Singlet")

saveRDS(d1,fp(out,"d1_singlet.rds"))


d1$sample<-sapply(d1$hash.ID,function(x){
  if(x=='CMO301') return('ztll2')
  if(x=='CMO302')return('ztll6')
    if(x=='CMO303')return('ztll10')

})

#rq: comment demultiplexer
d1_demux<-SplitObject(d1,split.by = "sample")


#and save only the matrix
d1_demux$ztll10@assays$RNA@counts

saveRDS(d1,fp(out,"d1_singlet.rds"))


"/disks/PROJECT/PhD_AlexandrePelletier/methyl/scripts/example_demultiplexing_singlecell_library.R"
