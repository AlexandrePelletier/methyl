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
#run 20A- [to get from singlecell project scripts]


#then use scVelo (https://scvelo.readthedocs.io/getting_started.html) in python to generate the RNA velocity matrix
# calculate and formulate the anndata

#test for 1 : CBP2
reticulate::py_run_string("
import scvelo as scv
filename='../singlecell/outputs/09-Velocity/velocyto_counts/cbp2/possorted_genome_bam_97NV3.loom'
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

#diff bias
#diff "intensity" by lineage
mtd<-data.table(cbps_h@meta.data,keep.rownames = "cell_id")
lins<-c("LT-HSC","HSC","MPP/LMPP","Myeloid","Lymphoid","Erythro-Mas")
mtdl<-mtd[lineage_hmap%in%lins]
mtdl[,lineage_hmap:=factor(lineage_hmap,levels = lins)]
mtdl[,velo_len_avg:=mean(velocity_length),by=.(lineage_hmap,sample)]

ggplot(unique(mtdl,by=c("sample","lineage_hmap")))+geom_boxplot(aes(x=lineage_hmap,y=velocity_length,fill=group))

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

ggplot(mtdff)+geom_boxplot(aes(x=lineage_hmap,y=pseudo_bias,fill=group))

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


#add mean regulation in HSC
avg_regul<-AverageExpression(cbps_h,assays = "velocity",features = dk_lga_hsc[p_val_adj<0.001]$gene)
head(avg_regul$velocity)
avg_velo_hsc<-data.table(gene=rownames(avg_regul$velocity),avg_velo=avg_regul$velocity[,"HSC"])
dk_lga_hsc_sig<-merge(dk_lga_hsc,avg_velo_hsc)
fwrite(dk_lga_hsc_sig,fp(out,"diff_kinetics_lga_vs_ctrl_hsc_with_avg_velo_padj0.001.csv.gz"))
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


#IEGs/EGRn gene velocity analysis ?####
#1) def IEGs/EGRn,  2) check in HSC, 3) compared to all lineage. 4) Which are diff velo in LGA ?


#doesnt see significant differences (at least at sample level)
#- due to bad velocity modelling ? ####
#=> dynamical modelling
# see https://scvelo.readthedocs.io/DynamicalModeling/
#run 20C
#due to not enough control cells (sans HTO) ?####



