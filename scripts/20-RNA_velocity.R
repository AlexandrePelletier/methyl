out<-"outputs/20-RNA_velocity"
dir.create(out)

source("scripts/utils/new_utils.R")

#infer RNA velocity using scVelo to 1) find if LGA HSC have a differentiation bias ?
#2) study which genes (IEGs ?) are freshly transcribed (++ premRNA/ mRNA) 
#and 3) if there is transcription differences between LGA and control

####need first install python####
renv::snapshot()
#rm python files

renv::use_python()

library(reticulate)

reticulate::py_install(packages = c("numpy" ,"scipy", "cython" ,"numba" ,"matplotlib", "scikit-learn", "h5py", "click"))
reticulate::py_install(packages ="umap-learn")
reticulate::py_install(packages ="pysam",pip=T)
reticulate::py_install(packages = c("velocyto"),pip = T)
reticulate::py_install(packages ="scvelo",pip=T)

#need downgrade numba to 0.52
reticulate::py_install(packages ="numba==0.52")
reticulate::py_run_string("
import scvelo as scv
scv.logging.print_version()
                          ") ##Running scvelo 0.2.4 (python 3.7.3) on 2022-04-07 13:23

renv::install("bioc::LoomExperiment")
renv::snapshot()

#need also samtools
system("../singlecell/tools/samtools-1.12/samtools")


#run velocyto counting pipeline
  #first, need dl genome repeat sequence to mask here : https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Human&db=0&hgta_group=allTracks&hgta_track=rmsk&hgta_table=rmsk&hgta_regionType=genome&position=&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf
dir.create("outputs/20-RNA_velocity/velocyto_counts")
#run 20A- 


#then use scVelo (https://scvelo.readthedocs.io/getting_started.html) in python to generate the RNA velocity matrix
# calculate and formulate the anndata

#test for 1 : CBP2
reticulate::py_run_string("
import scvelo as scv
filename='outputs/20-RNA_velocity/velocyto_counts/cbp2/possorted_genome_bam_97NV3.loom'
out='outputs/20-RNA_velocity'
#load loom file
adata = scv.read(filename, cache=True)

#adata.var_names_make_unique()
#scv.pl.proportions(adata)

# #basic preprocessing
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

#recover dynmics
scv.tl.recover_dynamics(adata)

scv.tl.velocity(adata, mode='dynamical')
adata.write(out+'/cbp2.h5ad', compression='gzip')

scv.tl.velocity_graph(adata)

adata.write(out+'/cbp2.h5ad', compression='gzip')

                          ")

#get velocity matrix 
#run 20B
velo_mat<-fread("outputs/20-RNA_velocity/cbp2_velocity_vec_matrix.csv")
sum(!is.na(as.matrix(velo_mat)))

#for all HTO data
#need first merge all velocito object #see https://github.com/basilkhuder/Seurat-to-RNA-Velocity#integrating-loom-file-and-meta-data
#with cells kept in seurat analysis
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")

mtd[,cell_id:=paste(str_extract(bc,"[ATCG]+"),orig.ident,sep="_")]
fwrite(mtd,"outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
#update reticulat because of bug
renv::install("rstudio/reticulate")

#get umap coord
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
umap_coord<-data.table(cbps@reductions$ref.umap@cell.embeddings,keep.rownames = "bc")
umap_coord<-merge(umap_coord,mtd[,.(bc,cell_id)])
fwrite(umap_coord,"outputs/06-integr_singlecell_cbps/umap_cbps.csv")

table(mtd[hto==T]$batch)

#run 20B-scvelo
#get velocity matrix
library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5ad", dest = "h5seurat", overwrite = TRUE)
cbps_velo <- LoadH5Seurat("outputs/20-RNA_velocity/cbps_hto_velo_anndata.h5seurat") #Error: Missing required datasets 'levels' and 'values'

cbps<-readRDS("outputs/06-integr_singlecell_cbps/cbps_filtered.rds")
head(cbps[[]])

velo<-fread("outputs/20-RNA_velocity/cbps_hto_velocity_matrix.csv")
velo[1:10,.(cell_id)]
velo<-as.matrix(t(data.frame(velo,row.names = "cell_id")))

mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdf<-mtd[cell_id%in%colnames(velo)]
cbps_h<-cbps[,mtdf$bc]
rm(cbps)
cbps_h<-RenameCells(cbps_h,new.names = mtdf$cell_id )

cbps_h[["velocity"]]<-CreateAssayObject(data =velo[2:nrow(velo),] )

mtdvelo<-fread("outputs/20-RNA_velocity/cbps_hto_velocity_metadata.csv")
mtdvelo[,cell_id:=V1]
cbps_h<-AddMetaData(cbps_h,metadata = data.frame(mtdvelo[,-c("V1","batch")],row.names = "cell_id"))
saveRDS(cbps_h,fp(out,"cbps_hto_with_velocity_assay.rds"))

#next : analyze diff bias, diff gene regulation by lineage (IEGs upreg in HSC ?), reduce upregulation IEGs in LGA HSC ?
cbps_h<-readRDS(fp(out,"cbps_hto_with_velocity_assay.rds"))

#diff bias####
#diff "intensity" by lineage
mtd<-data.table(cbps_h@meta.data,keep.rownames = "cell_id")
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")
mtdl<-mtd[lineage_hmap%in%lins]
mtdl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdl[,velo_len_avg:=mean(velocity_length),by=.(lineage_hmap,sample)]
ggplot(mtdl)+
  geom_boxplot(aes(x=group,y=velocity_length,fill=group,group=sample),outlier.shape = NA)+
  facet_wrap("lineage_hmap")+coord_cartesian(ylim = c(0,100))

ggplot(mtdl)+
  geom_boxplot(aes(x=group,y=velocity_length,fill=group),outlier.shape = NA)+
  facet_wrap("lineage_hmap")+coord_cartesian(ylim = c(0,65))

ggplot(unique(mtdl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=velo_len_avg,fill=group))

#Differentiation Analysis : ponderate average of the pseudotime of the predicted cells transition
#need=> transition graph : correl change of expression between 2 cells and predict change based on the RNA velocity
# need get velocity_graph
trans<-fread("outputs/20-RNA_velocity/cbps_hto_transition_matrix.csv")
trans<-as.matrix(t(data.frame(trans,row.names = "cell_id")))
dim(trans) #12684 12683
trans[1:10,1:10]
trans<-trans[2:nrow(trans),]

mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdf<-mtd[cell_id%in%colnames(trans)]
#add pseudotime
pseudo<-fread("outputs/12-Pseudotime/metadata_pseudotime.csv")
pseudo[,cell_id:=ps(str_extract(bc,"[ATCG]+"),orig.ident,sep='_')]
mtdf<-merge(mtdf,pseudo[,.(cell_id,pseudotime)])
mtdff<-mtdf[pseudotime!=Inf]
summary(mtdff$pseudotime)
mtdff[,pseudo_bias:=sapply(1:.N,function(i){
  cell<-cell_id[i]
  return(sum(pseudotime*trans[cell,cell_id]))
  })]

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group,group=sample),outlier.shape = NA)+
  facet_wrap("lineage_hmap")+
  coord_cartesian(ylim = c(0,75))

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group),outlier.shape = NA)+
  facet_wrap("lineage_hmap")+
  coord_cartesian(ylim = c(0,50))

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

mtdffl<-mtdff[lineage_hmap%in%lins]
mtdffl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdffl[,avg_pseudobias:=mean(pseudo_bias),by=.(lineage_hmap,sample)]

ggplot(unique(mtdffl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=avg_pseudobias,fill=group))
unique(mtdffl[lineage_hmap=="MPP/LMPP"],by=c("sample","lineage_hmap"))

#unsupervised velocity analysis####
#1) which genes are velo markers in HSC ? 2) which genes are diff velo in LGA HSC
#1) which genes are velo markers in HSC
#=> lineage upreg genes 
library(Seurat)
cbps_h<-readRDS(fp(out,"cbps_hto_with_velocity_assay.rds"))
DefaultAssay(cbps_h)<-"velocity"
Idents(cbps_h)<-"lineage_hmap"
m_lins<-FindAllMarkers(cbps_h,min.pct = 0)
m_lins<-data.table(m_lins)
fwrite(m_lins,fp(out,"velo_genes_lineage_markers.csv.gz"))
#add avg_diff and avg_mean
m_lins<-fread(fp(out,"velo_genes_lineage_markers.csv.gz"))
m_lins2<-Reduce(rbind,lapply(unique(m_lins$cluster),function(lin)merge(m_lins[cluster==lin],data.table(FoldChange(cbps_h,features =m_lins[cluster==lin]$gene ,ident.1 = lin,mean.fxn = rowMeans,fc.name ="avg_diff"),keep.rownames = "gene")[,.(gene,avg_diff)][,cluster:=lin])))
m_lins3<-Reduce(rbind,lapply(unique(m_lins$cluster),function(lin)merge(m_lins2[cluster==lin],melt(data.table(AverageExpression(cbps_h,assays = "velocity",features = m_lins2[cluster==lin]$gene)$velocity,
                                                                                                             keep.rownames = "gene")[,.SD,.SDcols=c("gene",lin)],
                                                                                                  id.vars="gene",
                                                                                                  measure.vars=lin,
                                                                                                  variable.name="cluster",
                                                                                                  value.name = "avg_cluster"))))

fwrite(m_lins3,fp(out,"velo_genes_lineage_markers.csv.gz"))
m_lins3<-fread(fp(out,"velo_genes_lineage_markers.csv.gz"))
m_lins3[,top10:=p_val_adj<0.001&avg_diff>0.1&rank(avg_diff)>(.N-10),.(cluster)]
ggplot(m_lins3[cluster%in%lins&avg_cluster>0.1],aes(x=log2(avg_diff),y=-log10(p_val_adj),col=p_val_adj<0.001&avg_diff>0.1))+
  geom_point()+
  facet_wrap("cluster")+
  geom_label_repel(aes(label=ifelse(top10,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#velo upreg lineage markers 
table(m_lins3[p_val_adj<0.001&avg_diff>0.1&avg_cluster>0.1]$cluster)

   # B cell          DC Erythro-Mas         HSC      LT-HSC    Lymphoid       Mk/Er    MPP/LMPP     Myeloid      T cell 
   #     1235         423        1178         961          34         859         457         609         809        1468 

#gost analysis by lineage
library(gprofiler2)
velo_markers_list<-split(m_lins3[p_val_adj<0.001&avg_diff>0.1&avg_cluster>0.1]$gene,m_lins3[p_val_adj<0.001&avg_diff>0.1&avg_cluster>0.1]$cluster)
res_gost<-gost(velo_markers_list,organism = "hsapiens",exclude_iea = T,sources=c('GO','KEGG','TF'),
               significant = T,domain_scope = "annotated")

saveRDS(res_gost,fp(out,"res_gost_lineage_velo_markers.rds"))
dt_gost<-data.table(res_gost$result)
fwrite(dt_gost,fp(out,"res_gost_lineage_velo_markers.csv"))

lapply(split(dt_gost,by = c("query","source")),function(x)head(x[order(p_value)]$term_name,20))
#top20 terms by lineage
# $`HSC.GO:BP`
#  [1] "regulation of macromolecule metabolic process"                          
#  [2] "regulation of nitrogen compound metabolic process"                      
#  [3] "regulation of metabolic process"                                        
#  [4] "regulation of primary metabolic process"                                
#  [5] "regulation of cellular metabolic process"                               
#  [6] "positive regulation of nitrogen compound metabolic process"             
#  [7] "negative regulation of macromolecule metabolic process"                 
#  [8] "positive regulation of macromolecule metabolic process"                 
#  [9] "negative regulation of metabolic process"                               
# [10] "positive regulation of cellular metabolic process"                      
# [11] "positive regulation of nucleobase-containing compound metabolic process"
# [12] "positive regulation of cellular process"                                
# [13] "negative regulation of cellular metabolic process"                      
# [14] "positive regulation of metabolic process"                               
# [15] "negative regulation of nitrogen compound metabolic process"             
# [16] "regulation of gene expression"                                          
# [17] "regulation of nucleobase-containing compound metabolic process"         
# [18] "negative regulation of biological process"                              
# [19] "negative regulation of cellular process"                                
# [20] "macromolecule metabolic process"                                        
# 
# $HSC.KEGG
#  [1] "Hepatitis B"                                     "Pathways in cancer"                             
#  [3] "Th17 cell differentiation"                       "Chronic myeloid leukemia"                       
#  [5] "Osteoclast differentiation"                      "Prostate cancer"                                
#  [7] "MAPK signaling pathway"                          "Human T-cell leukemia virus 1 infection"        
#  [9] "Colorectal cancer"                               "Growth hormone synthesis, secretion and action" 
# [11] "Thyroid hormone signaling pathway"               "Human cytomegalovirus infection"                
# [13] "TGF-beta signaling pathway"                      "Hepatitis C"                                    
# [15] "Acute myeloid leukemia"                          "Kaposi sarcoma-associated herpesvirus infection"
# [17] "Neurotrophin signaling pathway"                  "FoxO signaling pathway"                         
# [19] "Adherens junction"                               "ErbB signaling pathway"                         
# 
# $HSC.TF
#  [1] "Factor: ETF; motif: GVGGMGG; match class: 1"                  
#  [2] "Factor: ETF; motif: GVGGMGG"                                  
#  [3] "Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS; match class: 1"   
#  [4] "Factor: TAFII250; motif: RARRWGGCGGMGGNGR; match class: 1"    
#  [5] "Factor: MOVO-B; motif: GNGGGGG"                               
#  [6] "Factor: Churchill; motif: CGGGNN; match class: 1"             
#  [7] "Factor: BTEB3; motif: CCNNSCCNSCCCCKCCCCC; match class: 1"    
#  [8] "Factor: Sp1; motif: NGGGGCGGGGN; match class: 1"              
#  [9] "Factor: E2F-3; motif: GGCGGGN; match class: 1"                
# [10] "Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT; match class: 1"
# [11] "Factor: sp4; motif: NNGNARGRGGCGGRGCNNRR"                     
# [12] "Factor: TIEG1; motif: NCCCNSNCCCCGCCCCC; match class: 1"      
# [13] "Factor: sp4; motif: NNGNARGRGGCGGRGCNNRR; match class: 1"     
# [14] "Factor: MAZ; motif: GGGMGGGGSSGGGGGGGGGGGG; match class: 1"   
# [15] "Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1"            
# [16] "Factor: E2F-1; motif: NGGGCGGGARV; match class: 1"            
# [17] "Factor: BTEB1; motif: GGGGGCGGGGCNGSGGGNGS"                   
# [18] "Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS"                   
# [19] "Factor: SP1:SP3; motif: CCSCCCCCYCC"                          
# [20] "Factor: SP2; motif: GNNGGGGGCGGGGSN"



head(m_lins3[cluster=="HSC"][order(p_val_adj)][p_val_adj==0],100)
genes_of_int<-c("SOCS3","KLF2","KLF4","EGR1","JUNB","IER2","IER3","GADD45B","ID3")
m_lins3[cluster=="HSC"&gene%in%genes_of_int] #only EGR1 and KLF2 are ++ upreg in HSC
head(m_lins3[cluster=="HSC"&avg_diff>0&],100) #LDHB, RP prot, SLC38A9,IGFBP7, MAP2K6.. are downreg ++ in HSC

#2) which genes are diff velo in LGA HSC ?
#=>LGA vs Ctrl differential kinetics in HSC 
dk_lga_hsc<-data.table(FindMarkers(cbps_h,min.pct = 0,logfc.threshold = 0,
                                   ident.1 = "lga",ident.2 = "ctrl",mean.fxn = rowMeans,fc.name ="avg_diff" ,
                                   group.by = "group",subset.ident = "HSC"),
                       keep.rownames = "gene")
dk_lga_hsc[p_val_adj<0.001] #UBC ++ downreg, UBC gene is one of the two stress-regulated polyubiquitin genes (UBB and UBC) in mammals
head(dk_lga_hsc[p_val_adj<0.001],100)
dk_lga_hsc[gene%in%genes_of_int]

fwrite(dk_lga_hsc,fp(out,"diff_kinetics_lga_vs_ctrl_hsc.csv.gz"))

dk_lga_hsc<-fread(fp(out,"diff_kinetics_lga_vs_ctrl_hsc.csv.gz"))


dk_lga_hsc[,top20:=p_val_adj<0.001&abs(avg_diff)>0.1&rank(abs(avg_diff))>(.N-20)]
ggplot(dk_lga_hsc,aes(x=avg_diff,y=-log10(p_val_adj),col=p_val_adj<0.001&abs(avg_diff)>0.1))+
  geom_point()+
  geom_label_repel(aes(label=ifelse(top20,gene,"")),
                    size=3,
                   max.overlaps = 3000)+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")


#add mean regulation in HSC
avg_regul<-AverageExpression(cbps_h,assays = "velocity",features = dk_lga_hsc[p_val_adj<0.001]$gene)
head(avg_regul$velocity)
avg_velo_hsc<-data.table(gene=rownames(avg_regul$velocity),avg_velo=avg_regul$velocity[,"HSC"])
dk_lga_hsc_sig<-merge(dk_lga_hsc,avg_velo_hsc)
fwrite(dk_lga_hsc_sig,fp(out,"diff_kinetics_lga_vs_ctrl_hsc_with_avg_velo_padj0.001.csv.gz"))
dk_lga_hsc_sig<-fread(fp(out,"diff_kinetics_lga_vs_ctrl_hsc_with_avg_velo_padj0.001.csv.gz"))

dk_lga_hsc_sig[gene%in%genes_of_int]
head(dk_lga_hsc_sig[order(p_val_adj)],100)
dk_lga_hsc_sig[p_val_adj<10^-25&avg_velo<0&abs(avg_diff)>0.1][order(p_val_adj)]
dk_lga_hsc_sig[p_val_adj<10^-25&avg_velo>0&abs(avg_diff)>0.1][order(p_val_adj)]

dk_lga_hsc_sig[,regul:=ifelse(avg_diff<0,"down","up")]
velo_genes_hsc_lga_list<-split(dk_lga_hsc_sig[p_val_adj<0.001&abs(avg_velo)>0.1]$gene,dk_lga_hsc_sig[p_val_adj<0.001&abs(avg_velo)>0.1]$regul)
res_gost_hsc_lga<-gost(velo_genes_hsc_lga_list,organism = "hsapiens",exclude_iea = T,sources=c('GO','KEGG','TF'),
               significant = T,domain_scope = "annotated")

saveRDS(res_gost_hsc_lga,fp(out,"res_gost_lga_vs_ctrl_hsc_velo_padj0.001_absvelo0.1.rds"))
dt_gost<-data.table(res_gost_hsc_lga$result)
fwrite(dt_gost,fp(out,"res_gost_lga_vs_ctrl_hsc_velo_padj0.001_absvelo0.1.csv"))

lapply(split(dt_gost,by = c("query","source")),function(x)head(x[order(p_value)]$term_name,20))
# $`down.GO:BP`
#  [1] "cytoplasmic translation"                                       
#  [2] "cellular macromolecule metabolic process"                      
#  [3] "cellular metabolic process"                                    
#  [4] "nitrogen compound metabolic process"                           
#  [5] "macromolecule metabolic process"                               
#  [6] "primary metabolic process"                                     
#  [7] "nucleic acid metabolic process"                                
#  [8] "cellular nitrogen compound metabolic process"                  
#  [9] "regulation of primary metabolic process"                       
# [10] "RNA processing"                                                
# [11] "cellular process"                                              
# [12] "regulation of nitrogen compound metabolic process"             
# [13] "regulation of cellular metabolic process"                      
# [14] "RNA metabolic process"                                         
# [15] "nucleobase-containing compound metabolic process"              
# [16] "regulation of nucleobase-containing compound metabolic process"
# [17] "organic substance metabolic process"                           
# [18] "metabolic process"                                             
# [19] "gene expression"                                               
# [20] "regulation of macromolecule metabolic process"                 
# 

# $`down.GO:MF`
#  [1] "RNA binding"                              "heterocyclic compound binding"           
#  [3] "organic cyclic compound binding"          "nucleic acid binding"                    
#  [5] "protein binding"                          "binding"                                 
#  [7] "transcription coregulator activity"       "mRNA binding"                            
#  [9] "enzyme binding"                           "structural constituent of ribosome"      
# [11] "protein serine/threonine kinase activity" "histone binding"                         
# [13] "demethylase activity"                     "methylation-dependent protein binding"   
# [15] "methylated histone binding"               "transcription factor binding"            
# [17] "protein demethylase activity"             "histone demethylase activity"            
# 
# $down.KEGG
# [1] "Coronavirus disease - COVID-19"    "Ribosome"                          "Apoptosis"                        
# [4] "Thyroid hormone signaling pathway"
# 
# $down.TF
#  [1] "Factor: ZF5; motif: GGSGCGCGS; match class: 1"                     
#  [2] "Factor: ZF5; motif: GSGCGCGR; match class: 1"                      
#  [3] "Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN; match class: 1"           
#  [4] "Factor: E2F-2; motif: GCGCGCGCGYW"                                 
#  [5] "Factor: E2F-4; motif: SNGGGCGGGAANN; match class: 1"               
#  [6] "Factor: E2F-2; motif: GCGCGCGCGYW; match class: 1"                 
#  [7] "Factor: ZF5; motif: GGSGCGCGS"                                     
#  [8] "Factor: E2F; motif: GGCGSG"                                        
#  [9] "Factor: ETF; motif: GVGGMGG"                                       
# [10] "Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1"                 
# [11] "Factor: ZF5; motif: GSGCGCGR"                                      
# [12] "Factor: E2F-1:HES-7; motif: GGCRCGTGSYNNWNGGCGCSM; match class: 1" 
# [13] "Factor: E2F; motif: GGCGSG; match class: 1"                        
# [14] "Factor: E2F-3:HES-7; motif: NNNSGCGCSNNNNNCRCGYGNN; match class: 1"
# [15] "Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT; match class: 1"     
# [16] "Factor: ETF; motif: GVGGMGG; match class: 1"                       
# [17] "Factor: E2F-4; motif: SNGGGCGGGAANN"                               
# [18] "Factor: IRX-1; motif: NACRYNNNNNNNNRYGNN; match class: 1"          
# [19] "Factor: E2F-1; motif: NGGGCGGGARV"                                 
# [20] "Factor: Egr-1; motif: GCGCATGCG; match class: 1"                   
# 
# $`up.GO:BP`
#  [1] "cellular protein metabolic process"        "cellular macromolecule metabolic process" 
#  [3] "macromolecule metabolic process"           "cellular metabolic process"               
#  [5] "protein metabolic process"                 "cellular catabolic process"               
#  [7] "metabolic process"                         "primary metabolic process"                
#  [9] "regulation of cellular catabolic process"  "regulation of catabolic process"          
# [11] "organic substance metabolic process"       "catabolic process"                        
# [13] "protein modification process"              "cellular protein modification process"    
# [15] "nitrogen compound metabolic process"       "regulation of cellular metabolic process" 
# [17] "macromolecule catabolic process"           "macromolecule modification"               
# [19] "intracellular transport"                   "organonitrogen compound metabolic process"
# 
# $`up.GO:CC`
#  [1] "intracellular anatomical structure"       "intracellular organelle"                 
#  [3] "intracellular membrane-bounded organelle" "cytosol"                                 
#  [5] "nucleus"                                  "organelle"                               
#  [7] "membrane-bounded organelle"               "nucleoplasm"                             
#  [9] "nuclear lumen"                            "cytoplasm"                               
# [11] "protein-containing complex"               "intracellular organelle lumen"           
# [13] "membrane-enclosed lumen"                  "organelle lumen"                         
# [15] "catalytic complex"                        "nuclear body"                            
# [17] "nuclear protein-containing complex"       "intracellular protein-containing complex"
# [19] "nuclear speck"                            "transferase complex"                     
# 
# $`up.GO:MF`
#  [1] "RNA binding"                                 "protein binding"                            
#  [3] "enzyme binding"                              "catalytic activity, acting on a protein"    
#  [5] "heterocyclic compound binding"               "binding"                                    
#  [7] "nucleic acid binding"                        "organic cyclic compound binding"            
#  [9] "ubiquitin-like protein transferase activity" "ubiquitin-protein transferase activity"     
# [11] "catalytic activity"                          "transferase activity"                       
# [13] "mRNA binding"                                "ubiquitin-like protein binding"             
# [15] "kinase binding"                              "ubiquitin protein ligase binding"           
# [17] "ubiquitin-like protein ligase binding"       "cadherin binding"                           
# [19] "modification-dependent protein binding"      "protein kinase binding"                     
# 
# $up.KEGG
#  [1] "Ubiquitin mediated proteolysis"              "Autophagy - animal"                         
#  [3] "Endocytosis"                                 "Nucleocytoplasmic transport"                
#  [5] "Phosphatidylinositol signaling system"       "Shigellosis"                                
#  [7] "Amyotrophic lateral sclerosis"               "Neurotrophin signaling pathway"             
#  [9] "Chronic myeloid leukemia"                    "Salmonella infection"                       
# [11] "Inositol phosphate metabolism"               "Colorectal cancer"                          
# [13] "Endometrial cancer"                          "Longevity regulating pathway"               
# [15] "Bacterial invasion of epithelial cells"      "Spinocerebellar ataxia"                     
# [17] "Mitophagy - animal"                          "Sphingolipid signaling pathway"             
# [19] "Protein processing in endoplasmic reticulum" "Cellular senescence"                        
# 
# $up.TF
#  [1] "Factor: ETF; motif: GVGGMGG"                                       
#  [2] "Factor: ZF5; motif: GGSGCGCGS; match class: 1"                     
#  [3] "Factor: E2F-3:HES-7; motif: NNNSGCGCSNNNNNCRCGYGNN; match class: 1"
#  [4] "Factor: ETF; motif: GVGGMGG; match class: 1"                       
#  [5] "Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN; match class: 1"           
#  [6] "Factor: E2F; motif: GGCGSG"                                        
#  [7] "Factor: ZF5; motif: GSGCGCGR; match class: 1"                      
#  [8] "Factor: E2F-4; motif: SNGGGCGGGAANN; match class: 1"               
#  [9] "Factor: E2F-2; motif: GCGCGCGCGYW; match class: 1"                 
# [10] "Factor: E2F-2; motif: GCGCGCGCGYW"                                 
# [11] "Factor: E2F-1:HES-7; motif: GGCRCGTGSYNNWNGGCGCSM; match class: 1" 
# [12] "Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1"                 
# [13] "Factor: IRX-1; motif: NACRYNNNNNNNNRYGNN; match class: 1"          
# [14] "Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT"                     
# [15] "Factor: ZF5; motif: GGSGCGCGS"                                     
# [16] "Factor: E2F1; motif: GSGCGGGAAN"                                   
# [17] "Factor: E2F-1; motif: NGGGCGGGARV"                                 
# [18] "Factor: E2F-4; motif: NNTTCCCGCCNN"                                
# [19] "Factor: E2F-1; motif: NNNNGGCGGGAARN"                              
# [20] "Factor: E2F-3:HES-7; motif: NNNSGCGCSNNNNNCRCGYGNN"  
# 

#==> IN LGA HSC, downregulation of RNA velocity / decrease active transcription of methylation-dependent protein binding 
##and gene target of Factor: Egr-1; motif: GCGCATGCG; match class: 1

#more stringeant analysis

dk_lga_hsc_sig[,regul:=ifelse(avg_diff<0,"down","up")]
velo_genes_hsc_lga_list<-split(dk_lga_hsc_sig[p_val_adj<0.001&abs(avg_diff)>0.1]$gene,dk_lga_hsc_sig[p_val_adj<0.001&abs(avg_diff)>0.1]$regul)
lapply(velo_genes_hsc_lga_list, length)
# $down
# [1] 823
# 
# $up
# [1] 991
res_gost_hsc_lga<-gost(velo_genes_hsc_lga_list,organism = "hsapiens",exclude_iea = T,sources=c('GO','KEGG','TF'),
               significant = T,domain_scope = "annotated")

saveRDS(res_gost_hsc_lga,fp(out,"res_gost_lga_vs_ctrl_hsc_velo_padj0.001_absdiff0.1.rds"))
dt_gost<-data.table(res_gost_hsc_lga$result)
fwrite(dt_gost,fp(out,"res_gost_lga_vs_ctrl_hsc_velo_padj0.001_absdiff0.1.csv"))

lapply(split(dt_gost,by = c("query","source")),function(x)head(x[order(p_value)]$term_name,20))
# $`down.GO:BP`
#  [1] "cytoplasmic translation"                                                 "gene expression"                                                        
#  [3] "cellular macromolecule metabolic process"                                "macromolecule metabolic process"                                        
#  [5] "cellular process"                                                        "cytoskeleton organization"                                              
#  [7] "macromolecule biosynthetic process"                                      "RNA metabolic process"                                                  
#  [9] "cellular nitrogen compound metabolic process"                            "cellular nitrogen compound biosynthetic process"                        
# [11] "cellular protein metabolic process"                                      "nitrogen compound metabolic process"                                    
# [13] "translation"                                                             "peptide biosynthetic process"                                           
# [15] "regulation of macromolecule metabolic process"                           "regulation of metabolic process"                                        
# [17] "negative regulation of nucleobase-containing compound metabolic process" "organelle organization"                                                 
# [19] "peptide metabolic process"                                               "cellular metabolic process"                                             
# 
# $`down.GO:CC`
#  [1] "intracellular anatomical structure"           "intracellular organelle"                      "cytosol"                                     
#  [4] "non-membrane-bounded organelle"               "intracellular non-membrane-bounded organelle" "organelle"                                   
#  [7] "cytosolic ribosome"                           "nucleoplasm"                                  "protein-containing complex"                  
# [10] "nucleus"                                      "nuclear lumen"                                "membrane-bounded organelle"                  
# [13] "intracellular membrane-bounded organelle"     "polysome"                                     "ribonucleoprotein complex"                   
# [16] "cytoplasm"                                    "membrane-enclosed lumen"                      "organelle lumen"                             
# [19] "intracellular organelle lumen"                "cell-substrate junction"                     
# 
# $`down.GO:MF`
#  [1] "RNA binding"                              "structural constituent of ribosome"       "DNA-binding transcription factor binding"
#  [4] "nucleic acid binding"                     "protein binding"                          "transcription coregulator activity"      
#  [7] "heterocyclic compound binding"            "organic cyclic compound binding"          "transcription factor binding"            
# [10] "cadherin binding"                         "mRNA binding"                             "binding"                                 
# 
# $down.KEGG
# [1] "Coronavirus disease - COVID-19"                  "Ribosome"                                        "Epstein-Barr virus infection"                   
# [4] "Kaposi sarcoma-associated herpesvirus infection" "Notch signaling pathway"                         "Hepatitis B"                                    
# [7] "Thyroid hormone signaling pathway"              
# 
# $down.TF
#  [1] "Factor: ETF; motif: GVGGMGG; match class: 1"                       "Factor: ETF; motif: GVGGMGG"                                      
#  [3] "Factor: ZF5; motif: GGSGCGCGS; match class: 1"                     "Factor: MOVO-B; motif: GNGGGGG"                                   
#  [5] "Factor: E2F-2; motif: GCGCGCGCGYW"                                 "Factor: Churchill; motif: CGGGNN; match class: 1"                 
#  [7] "Factor: E2F; motif: GGCGSG; match class: 1"                        "Factor: E2F-1:HES-7; motif: GGCRCGTGSYNNWNGGCGCSM; match class: 1"
#  [9] "Factor: TAFII250; motif: RARRWGGCGGMGGNGR; match class: 1"         "Factor: E2F-3; motif: GGCGGGN"                                    
# [11] "Factor: E2F-4; motif: SNGGGCGGGAANN; match class: 1"               "Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN; match class: 1"          
# [13] "Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT; match class: 1"     "Factor: ZF5; motif: GSGCGCGR; match class: 1"                     
# [15] "Factor: E2F-2; motif: GCGCGCGCGYW; match class: 1"                 "Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS; match class: 1"       
# [17] "Factor: Sp1; motif: GGGGCGGGGC"                                    "Factor: E2F-3; motif: GGCGGGN; match class: 1"                    
# [19] "Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1"                 "Factor: Sp1; motif: NGGGGCGGGGN"                                  
# 
# $`up.GO:BP`
#  [1] "cellular macromolecule metabolic process"          "macromolecule metabolic process"                   "cellular metabolic process"                       
#  [4] "cellular protein metabolic process"                "protein metabolic process"                         "regulation of cellular metabolic process"         
#  [7] "organic substance metabolic process"               "metabolic process"                                 "regulation of metabolic process"                  
# [10] "primary metabolic process"                         "nitrogen compound metabolic process"               "regulation of nitrogen compound metabolic process"
# [13] "regulation of macromolecule metabolic process"     "regulation of primary metabolic process"           "regulation of catabolic process"                  
# [16] "organonitrogen compound metabolic process"         "cellular protein modification process"             "protein modification process"                     
# [19] "cellular catabolic process"                        "regulation of cellular catabolic process"         
# 
# $`up.GO:CC`
#  [1] "intracellular anatomical structure"       "cytosol"                                  "intracellular membrane-bounded organelle"
#  [4] "intracellular organelle"                  "nucleus"                                  "nucleoplasm"                             
#  [7] "cytoplasm"                                "nuclear lumen"                            "membrane-bounded organelle"              
# [10] "organelle"                                "protein-containing complex"               "membrane-enclosed lumen"                 
# [13] "intracellular organelle lumen"            "organelle lumen"                          "catalytic complex"                       
# [16] "nuclear body"                             "nuclear protein-containing complex"       "intracellular protein-containing complex"
# [19] "nuclear speck"                            "cytoplasmic ribonucleoprotein granule"   
# 
# $`up.GO:MF`
#  [1] "RNA binding"                                 "enzyme binding"                              "protein binding"                            
#  [4] "heterocyclic compound binding"               "organic cyclic compound binding"             "catalytic activity, acting on a protein"    
#  [7] "nucleic acid binding"                        "binding"                                     "kinase binding"                             
# [10] "ubiquitin-like protein transferase activity" "enzyme regulator activity"                   "transcription factor binding"               
# [13] "catalytic activity"                          "ubiquitin-protein transferase activity"      "mRNA binding"                               
# [16] "transferase activity"                        "protein kinase binding"                      "protein domain specific binding"            
# [19] "protein serine/threonine kinase activity"    "cadherin binding"                           
# 
# $up.KEGG
#  [1] "Ubiquitin mediated proteolysis"            "Autophagy - animal"                        "Shigellosis"                              
#  [4] "Sphingolipid signaling pathway"            "Phosphatidylinositol signaling system"     "Longevity regulating pathway"             
#  [7] "Neurotrophin signaling pathway"            "Fc epsilon RI signaling pathway"           "MAPK signaling pathway"                   
# [10] "Salmonella infection"                      "Phospholipase D signaling pathway"         "Regulation of actin cytoskeleton"         
# [13] "Acute myeloid leukemia"                    "Platelet activation"                       "Bacterial invasion of epithelial cells"   
# [16] "Prostate cancer"                           "Insulin signaling pathway"                 "EGFR tyrosine kinase inhibitor resistance"
# [19] "B cell receptor signaling pathway"         "Nucleocytoplasmic transport"              
# 
# $up.TF
#  [1] "Factor: ETF; motif: GVGGMGG"                                        "Factor: ETF; motif: GVGGMGG; match class: 1"                       
#  [3] "Factor: E2F-2; motif: GCGCGCGCGYW; match class: 1"                  "Factor: E2F-3:HES-7; motif: NNNSGCGCSNNNNNCRCGYGNN; match class: 1"
#  [5] "Factor: ZF5; motif: GGSGCGCGS; match class: 1"                      "Factor: ZF5; motif: GSGCGCGR; match class: 1"                      
#  [7] "Factor: E2F; motif: GGCGSG"                                         "Factor: E2F-2; motif: GCGCGCGCGYW"                                 
#  [9] "Factor: Egr-1; motif: GCGCATGCG; match class: 1"                    "Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS; match class: 1"        
# [11] "Factor: RNF96; motif: BCCCGCRGCC"                                   "Factor: E2F-4; motif: SNGGGCGGGAANN; match class: 1"               
# [13] "Factor: Egr-1; motif: GCGCATGCG"                                    "Factor: E2F-1; motif: NNNSSCGCSAANN"                               
# [15] "Factor: ZF5; motif: NRNGNGCGCGCWN; match class: 1"                  "Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1"                 
# [17] "Factor: E2F-1; motif: NNNNGGCGGGAARN"                               "Factor: sp4; motif: NNGNARGRGGCGGRGCNNRR"                          
# [19] "Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN; match class: 1"            "Factor: E2F-3; motif: GGCGGGN; match class: 1"       

dt_gost[source=="TF"&query=="down"][1:100][,.(p_value,term_name)]
#          p_value                                                            term_name
#   1: 1.442577e-27                          Factor: ETF; motif: GVGGMGG; match class: 1
#   2: 5.612617e-26                                          Factor: ETF; motif: GVGGMGG
#   3: 1.914500e-24                        Factor: ZF5; motif: GGSGCGCGS; match class: 1
#   4: 1.631698e-22                                       Factor: MOVO-B; motif: GNGGGGG
#   5: 2.532861e-21                                    Factor: E2F-2; motif: GCGCGCGCGYW
#   6: 1.163400e-20                     Factor: Churchill; motif: CGGGNN; match class: 1
#   7: 1.674011e-20                           Factor: E2F; motif: GGCGSG; match class: 1
#   8: 2.316826e-20    Factor: E2F-1:HES-7; motif: GGCRCGTGSYNNWNGGCGCSM; match class: 1
#   9: 4.910691e-20            Factor: TAFII250; motif: RARRWGGCGGMGGNGR; match class: 1
#  10: 6.138178e-20                                        Factor: E2F-3; motif: GGCGGGN
#  11: 9.871933e-20                  Factor: E2F-4; motif: SNGGGCGGGAANN; match class: 1
#  12: 1.469573e-19              Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN; match class: 1
#  13: 1.525609e-19        Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT; match class: 1
#  14: 3.453983e-19                         Factor: ZF5; motif: GSGCGCGR; match class: 1
#  15: 4.852621e-19                    Factor: E2F-2; motif: GCGCGCGCGYW; match class: 1
#  16: 5.448752e-19           Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS; match class: 1
#  17: 6.673311e-19                                       Factor: Sp1; motif: GGGGCGGGGC
#  18: 1.083631e-18                        Factor: E2F-3; motif: GGCGGGN; match class: 1
#  19: 1.372239e-18                    Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1
#  20: 3.018652e-18                                      Factor: Sp1; motif: NGGGGCGGGGN
#  21: 3.864190e-18             Factor: sp4; motif: NNGNARGRGGCGGRGCNNRR; match class: 1
#  22: 1.528971e-17                   Factor: E2F-1; motif: WWTGGCGCCAAA; match class: 1
#  23: 1.957637e-17                Factor: TCF-1; motif: ACATCGRGRCGCTGW; match class: 1
#  24: 4.349381e-17             Factor: IRX-1; motif: NACRYNNNNNNNNRYGNN; match class: 1
#  25: 5.305464e-17                                     Factor: RNF96; motif: BCCCGCRGCC
#  26: 5.720913e-17                                           Factor: E2F; motif: GGCGSG
#  27: 1.253312e-16                                  Factor: KROX; motif: CCCGCCCCCRCCCC
#  28: 1.290958e-16                    Factor: E2F-1; motif: NGGGCGGGARV; match class: 1
#  29: 1.416476e-16                       Factor: MOVO-B; motif: GNGGGGG; match class: 1
#  30: 1.440242e-16                                         Factor: ZF5; motif: GSGCGCGR
#  31: 1.816918e-16           Factor: MAZ; motif: GGGMGGGGSSGGGGGGGGGGGG; match class: 1
#  32: 1.918698e-16 Factor: AP-2gamma:Elk-1; motif: NGCCKNRGGSGRCGGAAGTG; match class: 1
#  33: 2.075431e-16                      Factor: Egr-1; motif: GCGCATGCG; match class: 1
#  34: 2.628345e-16                                    Factor: Sp1; motif: NNGGGGCGGGGNN
#  35: 2.785192e-16                             Factor: sp4; motif: NNGNARGRGGCGGRGCNNRR
#  36: 2.865972e-16                            Factor: BTEB3; motif: CCNNSCCNSCCCCKCCCCC
#  37: 4.528106e-16                                    Factor: E2F-1; motif: NGGGCGGGARV
#  38: 5.935704e-16                                        Factor: ZF5; motif: GGSGCGCGS
#  39: 7.707035e-16                               Factor: E2F-2; motif: NWTTTGGCGCCAWWNN
#  40: 8.380152e-16                  Factor: SP1:SP3; motif: CCSCCCCCYCC; match class: 1
#  41: 1.180900e-15                                     Factor: Egr-1; motif: GCGGGGGCGG
#  42: 1.195042e-15                                  Factor: E2F-1; motif: NNNSSCGCSAANN
#  43: 1.533039e-15   Factor: E2F-3:HES-7; motif: NNNSGCGCSNNNNNCRCGYGNN; match class: 1
#  44: 1.576342e-15                             Factor: IRX-1; motif: NACRYNNNNNNNNRYGNN
#  45: 2.023542e-15                                      Factor: EGR; motif: CGCCCCCGCNN
#  46: 2.873654e-15                                    Factor: Sp1; motif: NGGGGGCGGGGYN
#  47: 3.014528e-15              Factor: TIEG1; motif: NCCCNSNCCCCGCCCCC; match class: 1
#  48: 3.724149e-15                     Factor: E2F-4; motif: NTTTCSCGCC; match class: 1
#  49: 3.893226e-15             Factor: WT1; motif: RGGNGGGGGAGGRGGNGGRG; match class: 1
#  50: 7.931774e-15                                 Factor: CTCF; motif: NCCRSTAGGGGGCGC
#  51: 8.646630e-15                                  Factor: SP2; motif: GNNGGGGGCGGGGSN
#  52: 9.165508e-15                                     Factor: Churchill; motif: CGGGNN
#  53: 1.087766e-14                             Factor: KLF3; motif: NNNNNNGGGCGGGGCNNGN
#  54: 1.566566e-14                     Factor: RNF96; motif: BCCCGCRGCC; match class: 1
#  55: 1.713455e-14                                  Factor: SP1:SP3; motif: CCSCCCCCYCC
#  56: 2.080461e-14                                  Factor: E2F-4; motif: SNGGGCGGGAANN
#  57: 2.559053e-14               Factor: E2F-2; motif: NWTTTGGCGCCAWWNN; match class: 1
#  58: 2.647208e-14                                       Factor: WT1; motif: CGCCCCCNCN
#  59: 2.877408e-14                                   Factor: Sp1; motif: GGNDGGRGGCGGGG
#  60: 3.669186e-14                                  Factor: Sp2; motif: NYSGCCCCGCCCCCY
#  61: 3.681001e-14            Factor: BTEB3; motif: CCNNSCCNSCCCCKCCCCC; match class: 1
#  62: 3.809121e-14                           Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS
#  63: 4.912620e-14                        Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT
#  64: 5.157622e-14                                      Factor: Egr-1; motif: GCGCATGCG
#  65: 5.329684e-14                                     Factor: WT1; motif: GNGGGGGCGGGG
#  66: 9.197673e-14                             Factor: CPBP; motif: GNNRGGGHGGGGNNGGGRN
#  67: 1.073754e-13                                   Factor: SP1; motif: NRGKGGGCGGGGCN
#  68: 1.110999e-13                            Factor: TAFII250; motif: RARRWGGCGGMGGNGR
#  69: 1.202244e-13           Factor: MAZ; motif: GGGGGAGGGGGNGRGRRRGNRG; match class: 1
#  70: 1.323255e-13                     Factor: Egr-1; motif: GCGGGGGCGG; match class: 1
#  71: 1.361425e-13                                Factor: TCF-1; motif: ACATCGRGRCGCTGW
#  72: 1.478339e-13                    Factor: E2F-1:HES-7; motif: GGCRCGTGSYNNWNGGCGCSM
#  73: 1.490297e-13                                     Factor: E2F-4; motif: NTTTCSCGCC
#  74: 1.536193e-13                                    Factor: SP1; motif: GGCCCCGCCCCCN
#  75: 1.576931e-13                       Factor: WT1; motif: CGCCCCCNCN; match class: 1
#  76: 1.808850e-13                         Factor: p300; motif: ACNTCCG; match class: 1
#  77: 1.840476e-13                 Factor: E2F-1; motif: NNNNGGCGGGAARN; match class: 1
#  78: 2.159829e-13                 Factor: ZIC4; motif: NNCCNCCCRYNGYGN; match class: 1
#  79: 2.222785e-13             Factor: GCMa:Erg; motif: ATGCGGGCGGAARKG; match class: 1
#  80: 2.247058e-13                                   Factor: E2F-1; motif: WWTGGCGCCAAA
#  81: 2.856986e-13                       Factor: Sp1; motif: CCCCGCCCCN; match class: 1
#  82: 2.879733e-13          Factor: PATZ; motif: GGGGNGGGGGMKGGRRNGGNRN; match class: 1
#  83: 3.000669e-13                                       Factor: Sp1; motif: CCCCGCCCCN
#  84: 3.250043e-13                                  Factor: E2F-1:DP-1; motif: TTTCSCGC
#  85: 3.753701e-13                                    Factor: E2F-2; motif: GCGCGCGCNCS
#  86: 3.838223e-13                     Factor: WT1; motif: GNGGGGGCGGGG; match class: 1
#  87: 4.013774e-13              Factor: ZNF383; motif: SSNGGGMGGNGSNGGS; match class: 1
#  88: 4.119295e-13             Factor: CPBP; motif: GNNRGGGHGGGGNNGGGRN; match class: 1
#  89: 5.357700e-13                         Factor: E2F; motif: TTTSGCGS; match class: 1
#  90: 6.588175e-13                      Factor: Sp1; motif: NGGGGCGGGGN; match class: 1
#  91: 1.179991e-12                   Factor: Kaiso; motif: GCMGGGRGCRGS; match class: 1
#  92: 1.495160e-12             Factor: RUNX2; motif: WRACCGCANWAACCGCAN; match class: 1
#  93: 2.231909e-12                             Factor: WT1; motif: RGGNGGGGGAGGRGGNGGRG
#  94: 2.547394e-12                    Factor: ZF5; motif: NRNGNGCGCGCWN; match class: 1
#  95: 3.075408e-12                      Factor: E2F1; motif: GSGCGGGAAN; match class: 1
#  96: 3.559844e-12                      Factor: ETF; motif: CCCCGCCCCYN; match class: 1
#  97: 3.860973e-12                Factor: TCF-1; motif: ACATCGRGRCGCTGW; match class: 1
#  98: 4.878634e-12                                 Factor: ZIC4; motif: NNCCNCCCRYNGYGN
#  99: 5.912921e-12                                    Factor: ZF5; motif: NRNGNGCGCGCWN
# 100: 6.088935e-12                           Factor: Sp1; motif: NGGGGGCGGGGCCNGGGGGGGG
#diff 0.5
velo_genes_hsc_lga_list<-split(dk_lga_hsc_sig[p_val_adj<0.001&abs(avg_diff)>0.1]$gene,dk_lga_hsc_sig[p_val_adj<0.001&abs(avg_diff)>0.5]$regul)
lapply(velo_genes_hsc_lga_list, length)
# $down
# [1] 109
# 
# $up
# [1] 194
res_gost_hsc_lga<-gost(velo_genes_hsc_lga_list,organism = "hsapiens",exclude_iea = T,sources=c('GO','KEGG','TF'),
               significant = T,domain_scope = "annotated")

saveRDS(res_gost_hsc_lga,fp(out,"res_gost_lga_vs_ctrl_hsc_velo_padj0.001_absdiff0.5.rds"))
dt_gost<-data.table(res_gost_hsc_lga$result)
fwrite(dt_gost,fp(out,"res_gost_lga_vs_ctrl_hsc_velo_padj0.001_absdiff0.5.csv"))

lapply(split(dt_gost,by = c("query","source")),function(x)head(x[order(p_value)]$term_name,20))
# $`down.GO:BP`
#  [1] "cytoplasmic translation"                         "translation"                                     "peptide biosynthetic process"                   
#  [4] "peptide metabolic process"                       "gene expression"                                 "amide biosynthetic process"                     
#  [7] "macromolecule biosynthetic process"              "cellular nitrogen compound biosynthetic process" "macromolecule metabolic process"                
# [10] "cellular amide metabolic process"                "cellular nitrogen compound metabolic process"    "cellular macromolecule biosynthetic process"    
# [13] "metabolic process"                               "cellular biosynthetic process"                   "nitrogen compound metabolic process"            
# [16] "organic substance biosynthetic process"          "organic substance metabolic process"             "biosynthetic process"                           
# [19] "cellular macromolecule metabolic process"        "primary metabolic process"                      
# 
# $`down.GO:CC`
#  [1] "cytosolic ribosome"                           "ribosomal subunit"                            "ribosome"                                    
#  [4] "ribonucleoprotein complex"                    "intracellular organelle"                      "polysome"                                    
#  [7] "non-membrane-bounded organelle"               "intracellular non-membrane-bounded organelle" "cytosolic large ribosomal subunit"           
# [10] "nucleus"                                      "protein-containing complex"                   "cytosolic small ribosomal subunit"           
# [13] "organelle"                                    "cell-substrate junction"                      "intracellular anatomical structure"          
# [16] "nucleoplasm"                                  "small ribosomal subunit"                      "large ribosomal subunit"                     
# [19] "focal adhesion"                               "nuclear lumen"                               
# 
# $`down.GO:MF`
# [1] "structural constituent of ribosome" "structural molecule activity"       "nucleic acid binding"               "RNA binding"                       
# [5] "heterocyclic compound binding"      "organic cyclic compound binding"   
# 
# $down.KEGG
# [1] "Ribosome"                                        "Coronavirus disease - COVID-19"                  "Shigellosis"                                    
# [4] "Kaposi sarcoma-associated herpesvirus infection"
# 
# $down.TF
#  [1] "Factor: WT1; motif: SMCNCCNSC; match class: 1"                 "Factor: RNF96; motif: BCCCGCRGCC; match class: 1"             
#  [3] "Factor: E2F-4; motif: GCGGGAAANA; match class: 1"              "Factor: E2F; motif: TTTSGCGS; match class: 1"                 
#  [5] "Factor: E2F; motif: GGCGSG; match class: 1"                    "Factor: MOVO-B; motif: GNGGGGG"                               
#  [7] "Factor: TAFII250; motif: RARRWGGCGGMGGNGR; match class: 1"     "Factor: ZF5; motif: GGSGCGCGS; match class: 1"                
#  [9] "Factor: DB1; motif: GGRRRRGRRGGAGGGGGNGRRR"                    "Factor: E2F-1:Elk-1; motif: SGCGCSNNAMCGGAAGT; match class: 1"
# [11] "Factor: ZF5; motif: GSGCGCGR; match class: 1"                  "Factor: MAZ; motif: GGGGGAGGGGGNGRGRRRGNRG; match class: 1"   
# [13] "Factor: MAZ; motif: GGGGAGGG; match class: 1"                  "Factor: WT1; motif: SMCNCCNSC"                                
# [15] "Factor: E2F-2; motif: GCGCGCGCNCS; match class: 1"             "Factor: CKROX; motif: SCCCTCCCC"                              
# [17] "Factor: CKROX; motif: SCCCTCCCC; match class: 1"               "Factor: E2F-2; motif: NWTTTGGCGCCAWWNN"                       
# [19] "Factor: E2F; motif: TTTCGCGC"                                  "Factor: Zbtb37; motif: NYACCGCRNTCACCGCR; match class: 1"     
# 
# $`up.GO:BP`
#  [1] "protein metabolic process"                                               "positive regulation of nitrogen compound metabolic process"             
#  [3] "cellular protein metabolic process"                                      "cellular protein modification process"                                  
#  [5] "protein modification process"                                            "positive regulation of cellular metabolic process"                      
#  [7] "positive regulation of macromolecule metabolic process"                  "regulation of nitrogen compound metabolic process"                      
#  [9] "macromolecule modification"                                              "positive regulation of metabolic process"                               
# [11] "cellular macromolecule metabolic process"                                "organonitrogen compound metabolic process"                              
# [13] "positive regulation of nucleobase-containing compound metabolic process" "regulation of macromolecule metabolic process"                          
# [15] "regulation of cellular metabolic process"                                "positive regulation of RNA metabolic process"                           
# [17] "regulation of primary metabolic process"                                 "macromolecule metabolic process"                                        
# [19] "positive regulation of cellular process"                                 "nitrogen compound metabolic process"                                    
# 
# $`up.GO:CC`
#  [1] "cytosol"                                  "intracellular membrane-bounded organelle" "intracellular anatomical structure"      
#  [4] "nucleus"                                  "nucleoplasm"                              "nuclear lumen"                           
#  [7] "cytoplasm"                                "intracellular organelle"                  "membrane-bounded organelle"              
# [10] "organelle"                                "organelle lumen"                          "membrane-enclosed lumen"                 
# [13] "intracellular organelle lumen"            "intracellular protein-containing complex" "membrane coat"                           
# [16] "coated membrane"                          "catalytic complex"                        "clathrin coat"                           
# [19] "protein-containing complex"               "transferase complex"                     
# 
# $`up.GO:MF`
# [1] "enzyme binding"                          "catalytic activity, acting on a protein" "protein binding"                        
# 
# $up.KEGG
# [1] "Ubiquitin mediated proteolysis"              "Sphingolipid signaling pathway"              "Platelet activation"                        
# [4] "Long-term potentiation"                      "Protein processing in endoplasmic reticulum"
# 
# $up.TF
#  [1] "Factor: ETF; motif: GVGGMGG; match class: 1"                "Factor: ETF; motif: GVGGMGG"                               
#  [3] "Factor: TIEG1; motif: NCCCNSNCCCCGCCCCC; match class: 1"    "Factor: ZXDL; motif: GSGSCNNGGGMRGCNCCGGGS; match class: 1"
#  [5] "Factor: TCF-1; motif: ACATCGRGRCGCTGW; match class: 1"      "Factor: E2F-1; motif: NNNNGGCGGGAARN; match class: 1"      
#  [7] "Factor: Sp1; motif: CCCCGCCCCN"                             "Factor: E2F-4; motif: NNTTCCCGCCNN"                        
#  [9] "Factor: RUNX3; motif: NRACCGCANWAACCRCAN; match class: 1"   "Factor: E2F-3; motif: GGCGGGN; match class: 1"             
# [11] "Factor: AP-2; motif: MKCCCSCNGGCG; match class: 1"          "Factor: BTEB3; motif: CCNNSCCNSCCCCKCCCCC"                 
# [13] "Factor: E2F-1; motif: NNNSSCGCSAANN"                        "Factor: KROX; motif: CCCGCCCCCRCCCC"                       
# [15] "Factor: BTEB3; motif: CCNNSCCNSCCCCKCCCCC; match class: 1"  "Factor: E2F-1; motif: NNNNGGCGGGAARN"                      
# [17] "Factor: E2F-1:DP-2; motif: TTTSSCGC; match class: 1"        "Factor: AP-2; motif: GSCCSCRGGCNRNRNN; match class: 1"     
# [19] "Factor: E2F3; motif: NNRGMKGGAR"                            "Factor: RNF96; motif: BCCCGCRGCC"     
# 

#IEGs/EGRn gene velocity analysis ?####
#1) def IEGs/EGRn,  2) more active than others in HSC ?, 3) compared to all lineage. 4) Which are diff velo in LGA ?
#1) def IEGs/EGRn
#IEGs/EGRn = all genes of the EGRn (all altered regulons in LGA)
tftargets<-fread("outputs/16-GRN_final/tf_target_interactions.csv")
altered_regulons<-c("KLF2","KLF4","EGR1","JUN","ARID5A","FOSB","KLF10")
egrn_genes<-unique(c(tftargets[tf%in%altered_regulons]$target,altered_regulons))
all_genes<-unique(c(tftargets$target,tftargets$tf))

length(egrn_genes) #518
length(all_genes) #2459

cbps_h<-readRDS(fp(out,"cbps_hto_with_velocity_assay.rds"))
#2) more active than others in HSC ?
#increased velocity in HSC
DefaultAssay(cbps_h)<-"velocity"
velo<-data.table(t(as.matrix(cbps_h@assays$velocity@data[intersect(all_genes,rownames(cbps_h)),])),keep.rownames = "cell")
velo<-melt(velo,id.vars = "cell",variable.name ="gene",value.name = "velocity" )

mtd<-data.table(cbps_h@meta.data,keep.rownames = "cell")
velo_dt<-merge(velo,mtd)
velo_dt
velo_dt[,egrn_gene:=gene%in%egrn_genes]
velo_dt[,velo_mean:=mean(velocity),by=.(cell,egrn_gene)]
ggplot(unique(velo_dt,by=c("cell","egrn_gene"))[lineage_hmap=="HSC"&group=="ctrl"])+geom_boxplot(aes(x=egrn_gene,y=velo_mean))

#! because egrn genes are genes express in HSC ?
#=> kept only genes express in more than 10% of HSC
DefaultAssay(cbps_h)<-"RNA"
hscs<-colnames(cbps_h)[cbps_h$lineage_hmap=="HSC"]
genes_hsc<-rownames(cbps_h)[rowSums(cbps_h@assays$RNA@counts[,hscs]>0)>0.1*length(hscs)]
length(genes_hsc) #8k
head(genes_hsc)

velo_dt[,express_in_hsc:=gene%in%genes_hsc]
velo_dt[(express_in_hsc),velo_mean:=mean(velocity),by=.(cell,egrn_gene)]

ggplot(unique(velo_dt,by=c("cell","egrn_gene"))[lineage_hmap=="HSC"&group=="ctrl"&(express_in_hsc)])+
  geom_boxplot(aes(x=egrn_gene,y=velo_mean))


#3) compared to all lineage.
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

velo_dt[,express_in_lin:={
  cells<-colnames(cbps_h)[cbps_h$lineage_hmap==lineage_hmap[1]]
  genes<-rownames(cbps_h)[rowSums(cbps_h@assays$RNA@counts[,cells]>0)>0.1*length(cells)]
  gene%in%genes
  },by="lineage_hmap"]
velo_dt[,ngenes.express:=length(unique(gene[express_in_lin==T])),by="lineage_hmap"]
unique(velo_dt[,.(lineage_hmap,ngenes.express)])
#     lineage_hmap ngenes.express
#  1:     MPP/LMPP           2284
#  2:          HSC           2354
#  3:  Erythro-Mas           2342
#  4:     Lymphoid           2094
#  5:       B cell           1604
#  6:      Myeloid           2386
#  7:       T cell           1546
#  8:           DC           1621
#  9:       LT-HSC           2254
# 10:        Mk/Er           1586

velo_dt[,velo_mean:=mean(velocity[express_in_lin==T]),by=.(cell,egrn_gene,lineage_hmap)]



ggplot(unique(velo_dt[lineage_hmap%in%lins],by=c("cell","egrn_gene"))[group=="ctrl"])+
  geom_boxplot(aes(x=egrn_gene,y=velo_mean))+facet_wrap("lineage_hmap")

#4) Which are diff velo in LGA ?
ggplot(unique(velo_dt[lineage_hmap%in%lins],by=c("cell","egrn_gene")))+
  geom_boxplot(aes(x=egrn_gene,y=velo_mean,fill=group))+facet_wrap("lineage_hmap")
#not really for HSC egrn genes, but seems increase velo for non egrn genes in HSC,
#and more importantly in MPP while a decrease of egrn genes velo.
#signif ?
velo_dtc<-unique(velo_dt,by=c("cell","egrn_gene"))
velo_dtc[,pval:=wilcox.test(velo_mean[group=="lga"],velo_mean[group=="ctrl"])$p.value,by=c("lineage_hmap","egrn_gene")]
unique(velo_dtc[,.(egrn_gene,lineage_hmap,pval)])

#     egrn_gene lineage_hmap         pval
#  1:      TRUE     MPP/LMPP 1.461996e-11
#  2:     FALSE     MPP/LMPP 5.888737e-41
#  3:      TRUE          HSC 3.640636e-01
#  4:     FALSE          HSC 5.959674e-38
#  5:      TRUE  Erythro-Mas 1.511895e-01
#  6:     FALSE  Erythro-Mas 8.828127e-16
#  7:      TRUE     Lymphoid 9.740312e-01
#  8:     FALSE     Lymphoid 2.638874e-11
#  9:      TRUE       B cell 2.623415e-04
# 10:     FALSE       B cell 1.691000e-01
# 11:      TRUE      Myeloid 2.017183e-06
# 12:     FALSE      Myeloid 1.916581e-24
# 13:      TRUE       T cell 2.245428e-04
# 14:     FALSE       T cell 5.752957e-03
# 15:      TRUE           DC 8.506311e-01
# 16:     FALSE           DC 8.506311e-01
# 17:      TRUE       LT-HSC 4.271997e-03
# 18:     FALSE       LT-HSC 5.781954e-03
# 19:      TRUE        Mk/Er 7.491968e-02
# 20:     FALSE        Mk/Er 1.661697e-01

#most velo bias egrn genes ? [TODO]

#doesnt see significant differences (at least at sample level)
#- due to bad velocity modelling ? ####
#=> dynamical modelling
# see https://scvelo.readthedocs.io/DynamicalModeling/
#run 20C
velo<-fread("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_matrix.csv")
velo[1:10,.(cell_id)]
velo<-as.matrix(t(data.frame(velo,row.names = "cell_id")))

cbps_h[["velocity"]]<-CreateAssayObject(data =velo[2:nrow(velo),] )

mtdvelo<-fread("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_metadata.csv")
mtdvelo[,cell_id:=V1]
cbps_h<-AddMetaData(cbps_h,metadata = data.frame(mtdvelo[,-c("V1","batch")],row.names = "cell_id"))
saveRDS(cbps_h@assays$velocity,fp(out,"cbps_hto_dynamical_velocity_assay.rds"))

#velo length bias ?
mtd<-data.table(cbps_h@meta.data,keep.rownames = "cell_id")
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")
mtdl<-mtd[lineage_hmap%in%lins]
mtdl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdl[,velo_len_avg:=mean(velocity_length),by=.(lineage_hmap,sample)]
ggplot(mtdl)+
  geom_boxplot(aes(x=group,y=velocity_length,fill=group,group=sample),outlier.shape = NA)+
  facet_wrap("lineage_hmap")+coord_cartesian(ylim = c(0,100))

ggplot(mtdl)+
  geom_boxplot(aes(x=group,y=velocity_length,fill=group),outlier.shape = NA)+
  facet_wrap("lineage_hmap")+coord_cartesian(ylim = c(0,65))

ggplot(unique(mtdl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=velo_len_avg,fill=group))

#Differentiation Analysis with pseudotime####
trans<-fread("outputs/20-RNA_velocity/cbps_hto_dynamical_transition_matrix.csv")
trans<-as.matrix(t(data.frame(trans,row.names = "cell_id")))
dim(trans) #12684 12683
trans[1:10,1:10]
trans<-trans[2:nrow(trans),]

mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtdf<-mtd[cell_id%in%colnames(trans)]

#add pseudotime
pseudo<-fread("outputs/12-Pseudotime/metadata_pseudotime.csv")
pseudo[,cell_id:=ps(str_extract(bc,"[ATCG]+"),orig.ident,sep='_')]
mtdf<-merge(mtdf,pseudo[,.(cell_id,pseudotime)])
mtdff<-mtdf[pseudotime!=Inf]
summary(mtdff$pseudotime)
mtdff[,pseudo_pred:=sapply(1:.N,function(i){
  cell<-cell_id[i]
  return(sum(pseudotime*trans[cell,cell_id]))
  })]

mtdff[,pseudo_bias:=pseudo_pred-pseudotime]


fwrite(mtdff,fp(out,"pseudo_bias_rna_velo_based_lineages.csv.gz"))
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")


mtdff[,lineage_hmap:=factor(lineage_hmap,levels=c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas"))]
ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group,group=sample),outlier.shape = NA)+
  facet_wrap("lineage_hmap")

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))

ggplot(mtdff[lineage_hmap%in%lins])+
  geom_boxplot(aes(x=group,y=pseudo_bias,fill=group),outlier.shape = NA)+
  facet_wrap("lineage_hmap")

lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")

mtdffl<-mtdff[lineage_hmap%in%lins]
mtdffl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdffl[,avg_pseudobias:=mean(pseudo_bias),by=.(lineage_hmap,sample)]

ggplot(unique(mtdffl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=avg_pseudobias,fill=group))
unique(mtdffl[lineage_hmap=="MPP/LMPP"],by=c("sample","lineage_hmap"))

#latent time (~pseudotime) analysis ####
#what lineage/ct are the bigger latent time
mtdvelo<-fread("outputs/20-RNA_velocity/cbps_hto_dynamical_velocity_metadata.csv")
mtdvelo[,cell_id:=V1]

mtdcbps<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")
mtd<-merge(mtdvelo[,-"batch"],mtdcbps)
ggplot(mtd)+geom_boxplot(aes(x=cell_type_hmap,y=latent_time))
fwrite(mtd,fp(out,"metadata_cbps_hto_dynamical_merged.csv"))

#keys Genes following main dynamics [todo, see tuto]####
#Top-likelihood genes
#=> assume Driver genes because display pronounced dynamic behavior

genes_mtd<-fread("outputs/20-RNA_velocity/cbps_hto_dynamical_genes_metadata.csv")
genes_mtd<-genes_mtd[!is.na(fit_likelihood)]
genes_mtd#1742
head(genes_mtd[order(-fit_likelihood)],20)
head(genes_mtd[order(-fit_likelihood),.(Gene,fit_likelihood)],50)
#             Gene fit_likelihood
#  1:        RPL41      0.7024196
#  2:         FTH1      0.4409182
#  3:        RPS21      0.4193488
#  4:        RPL26      0.4011444
#  5:   AC002454.1      0.3567726
#  6:        SERF2      0.3304627
#  7:         TMA7      0.3147169
#  8:        MYADM      0.3128698
#  9:        SNHG5      0.3084447
# 10:         CPA3      0.3069082
# 11:        SSBP1      0.3048434
# 12:       SH2D4B      0.3040029
# 13:         MYL6      0.2988126
# 14:         AREG      0.2983997
# 15:       RPL27A      0.2977637
# 16:          FOS      0.2969011
# 17:         CD99      0.2942007
# 18:        PPM1G      0.2932404
# 19:      SEPTIN6      0.2924571
# 20:       TIPARP      0.2923209
# 21:         FOSB      0.2871259
# 22: ADAMTSL4-AS1      0.2866284
# 23:       ATP5MG      0.2844941
# 24:      POU2AF1      0.2822485
# 25:       IGFBP7      0.2806364
# 26:       ADGRG6      0.2798706
# 27:        COPS9      0.2797676
# 28:    LINC01220      0.2766412
# 29:        DUSP1      0.2734219
# 30:        MLLT3      0.2731625
# 31:        RPS4X      0.2726536
# 32:   AC008440.1      0.2719394
# 33:         KLF2      0.2718106
# 34:       NDUFS5      0.2700264
# 35:         ACTB      0.2697666
# 36:       AKAP12      0.2694474
# 37:       COX6B1      0.2688626
# 38:        RBPMS      0.2686618
# 39:       SPINK2      0.2678411
# 40:         ATF3      0.2665551
# 41:        SNHG7      0.2664279
# 42:        MYBL2      0.2653777
# 43:        TSTD1      0.2645509
# 44:       RNF130      0.2644003
# 45:       PET100      0.2629900
# 46:        SNRPF      0.2615482
# 47:   SPACA6P-AS      0.2603058
# 48:        IFFO2      0.2587123
# 49:         NASP      0.2567542
# 50:         GYPC      0.2563279
#             Gene fit_likelihood

#likelihood by lineage (lineage specific driver genes)
dyna_genes<-fread("outputs/20-RNA_velocity/dynamical_genes_by_lineage.csv",header = T)
dyna_genes<-melt(dyna_genes[,-"V1"],value.name = "gene",variable.name = "lineage",measure.vars = colnames(dyna_genes)[-1])
dyna_genes[lineage=="HSC"]
#     lineage         gene
#   1:     HSC        RPL41
#   2:     HSC        RPS21
#   3:     HSC   AC002454.1
#   4:     HSC         CD99
#   5:     HSC      SEPTIN6
#   6:     HSC       ATP5MG
#   7:     HSC         TMA7
#   8:     HSC       SPINK2
#   9:     HSC        MLLT3
#  10:     HSC         NASP
#  11:     HSC        RBPMS
#  12:     HSC         RHOH
#  13:     HSC        SNHG5
#  14:     HSC        MYADM
#  15:     HSC          TKT
#  16:     HSC          FOS
#  17:     HSC        SNHG7
#  18:     HSC         FTH1
#  19:     HSC       COX6B1
#  20:     HSC        SNRPF
#  21:     HSC       TMSB10
#  22:     HSC       EPSTI1
#  23:     HSC       RPL27A
#  24:     HSC       RNF130
#  25:     HSC      PIP4K2A
#  26:     HSC         GBP2
#  27:     HSC        YPEL5
#  28:     HSC         PTMA
#  29:     HSC       SMIM24
#  30:     HSC         GYPC
#  31:     HSC         ATF3
#  32:     HSC         KLF4 *
#  33:     HSC         ZNF3
#  34:     HSC        ZFP36
#  35:     HSC         KLF2 *
#  36:     HSC ADAMTSL4-AS1
#  37:     HSC         FOSB
#  38:     HSC     SERPINB1
#  39:     HSC         DTD1
#  40:     HSC        STX11
#  41:     HSC         BEX1
#  42:     HSC         JUND
#  43:     HSC        H3F3B
#  44:     HSC         AREG
#  45:     HSC       TUBA1A
#  46:     HSC         YBX3
#  47:     HSC        SOCS2
#  48:     HSC         MCL1
#  49:     HSC        DUSP1
#  50:     HSC       DYNLL1
#  51:     HSC         FDX1
#  52:     HSC        RPL26
#  53:     HSC         EGR1 *
#  54:     HSC        PTGS2
#  55:     HSC         MAFF
#  56:     HSC          IDS
#  57:     HSC       TIPARP
#  58:     HSC        KLF10
#  59:     HSC         SRGN
#  60:     HSC         BTG2
#  61:     HSC        BUD31
#  62:     HSC      TSC22D2
#  63:     HSC         YRDC
#  64:     HSC       NFKBIA
#  65:     HSC       CHST11
#  66:     HSC      CHORDC1
#  67:     HSC         EMP3
#  68:     HSC       PTGER4
#  69:     HSC        NR4A1
#  70:     HSC     PPP1R15A
#  71:     HSC      SERTAD3
#  72:     HSC         LMNA
#  73:     HSC        HSPA8
#  74:     HSC         IER2
#  75:     HSC       DNAJB1
#  76:     HSC         NEU1
#  77:     HSC        RPS4X
#  78:     HSC        RUNX3
#  79:     HSC        IFFO2
#  80:     HSC        HSPH1
#  81:     HSC         EREG
#  82:     HSC          AVP
#  83:     HSC          ID2
#  84:     HSC        ARL4A
#  85:     HSC          PNP
#  86:     HSC          KIN
#  87:     HSC        CXCL8
#  88:     HSC        UBE2S
#  89:     HSC       MAGED2
#  90:     HSC         TOB2
#  91:     HSC      CCDC173
#  92:     HSC       ARRDC2
#  93:     HSC        PPM1N
#  94:     HSC        NR4A2
#  95:     HSC   AL021155.5
#  96:     HSC         MLF1
#  97:     HSC         SKIL
#  98:     HSC       NDUFS5
#  99:     HSC         ICA1
# 100:     HSC         CD74
#      lineage         gene


#gene level velocity analysis [todo, redo]####
#IEGs/ EGRns
#IEGs / EGRns LGAs bias


#Velocity on control cells (sans HTO)  ####
#est ce que la dynamique de velocity est inversée ? (va vers les cellules differenciées?)
#run 20D
# => non

#Top-likelihood genes compared to HTO cells





