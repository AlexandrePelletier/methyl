### Project Setup ==================================================================================
out<-"outputs/16-GRN_final"
dir.create(out)
source("scripts/utils/new_utils.R")
library(Seurat)
library(Signac)

####Functions####
GetMotifIDs<-function(object,motif.names,assay=NULL,return_dt=FALSE){
  if(is.null(assay))assay<-DefaultAssay(object)
  idx<-match(motif.names,object@assays[[assay]]@motifs@motif.names)
  if(return_dt){
    return(
      data.table(motif.name=motif.names,
                 motif.id=names(object@assays[[assay]]@motifs@motif.names[idx]))
      )
    }else{
  return(names(object@assays[[assay]]@motifs@motif.names[idx]))
    }
  
  }
CheckMotif<-function(object,peaks,motif.name,assay = NULL,return.peaks=FALSE){
  require("Signac")
  if(is.null(assay))assay<-DefaultAssay(object)
  motif<-GetMotifID(object,motif.name,assay=assay)
  motif.all <- GetMotifData(
    object = object, assay = assay, slot = "data"
  )
  
  motifs_peaks_tf <- motif.all[peaks,motif , drop = FALSE]
  if(return.peaks){
    motifs_peaks_tf<-rownames(motifs_peaks_tf)[as.vector(motifs_peaks_tf==1)]
    return(motifs_peaks_tf)
  }else{
    motifs_peaks_tf_vec<-as.vector(motifs_peaks_tf==1)
    names(motifs_peaks_tf_vec)<-rownames(motifs_peaks_tf)
    return(motifs_peaks_tf_vec)
  }
  
 
}



  
### Analysis =======================================================================================


#clean regulons list based on atac
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")

atacs<-readRDS("outputs/14-DMCs_atac_integr/cbps_atacs.rds")
atacs[["lin_peaks"]]<-readRDS("outputs/14-DMCs_atac_integr/cbps_lin_spe_peaks_assay.rds")
atacs@assays$lin_peaks@motifs<-readRDS("outputs/14-DMCs_atac_integr/atacs_cbps_lin_peaks_motif_object.rds")
DefaultAssay(atacs)<-"lin_peaks"

#for EGR1
peaks_hsc_genes<-fread("outputs/14-DMCs_atac_integr/peaks_hsc_genes_anno.csv.gz")

peaks_close_EGR1_target<-peaks_hsc_genes[gene_name%in%regulons_list$EGR1]$query_region
egr1_peaks<-CheckMotif(atacs,
                       peaks =peaks_close_EGR1_target ,
                       motif.name = "EGR1",
                       return.peaks = TRUE)
length(egr1_peaks)/length(peaks_close_EGR1_target) #30%
egr1_genes<-intersect(peaks_hsc_genes[query_region%in%egr1_peaks]$gene_name,
                      regulons_list$EGR1)
length(egr1_genes)/length(regulons_list$EGR1) #96% ! (25/26)

#egr1_extended
peaks_close_EGR1_target<-peaks_hsc_genes[gene_name%in%regulons_list$EGR1e]$query_region
egr1_peaks<-CheckMotif(atacs,
                       peaks =peaks_close_EGR1_target ,
                       motif.name = "EGR1",
                       return.peaks = TRUE)
length(egr1_peaks)/length(peaks_close_EGR1_target) #17%
egr1_genes<-intersect(peaks_hsc_genes[query_region%in%egr1_peaks]$gene_name,
                      regulons_list$EGR1e)
length(egr1_genes)/length(regulons_list$EGR1e) #93% (393/424)

#for all
tfs_scenic<-unique(str_remove(names(regulons_list),"e$"))
regulons_tf_atac<-unlist(atacs@assays$lin_peaks@motifs@motif.names[atacs@assays$lin_peaks@motifs@motif.names%in%tfs_scenic])
length(regulons_tf_atac)/length(tfs_scenic)#107/157

regulons_atac_list<-regulons_list[str_remove(names(regulons_list),"e$")%in%regulons_tf_atac]
length(regulons_atac_list) #174
length(regulons_list) #250


regulons_atac_listf<-lapply(names(regulons_atac_list), function(regulon_name){
  targets<-regulons_atac_list[[regulon_name]]
  motif_name<-str_remove(regulon_name,"e$")
  peaks_close_targets<-peaks_hsc_genes[gene_name%in%targets]$query_region
  tf_peaks<-CheckMotif(atacs,
                         peaks =peaks_close_targets ,
                         motif.name = motif_name,
                         return.peaks = TRUE)
  filtered_targets<-intersect(peaks_hsc_genes[query_region%in%tf_peaks]$gene_name,
                        targets)
  return(filtered_targets)
  })

names(regulons_atac_listf)<-names(regulons_atac_list)

cat(unlist(lapply(1:length(regulons_atac_listf),
       function(i)paste(names(regulons_atac_listf)[i],"=",round(length(regulons_atac_listf[[i]])/length(regulons_atac_list[[i]])*100),"%"))),sep = "\n")

mean(unlist(lapply(1:length(regulons_atac_listf),
       function(i)length(regulons_atac_listf[[i]])/length(regulons_atac_list[[i]])))) #59%


#make a df of interactions tf > targets
regulons<-Reduce(rbind,lapply(names(regulons_atac_listf), function(t)data.table(tf=rep(t,length(regulons_atac_listf[[t]])),target=regulons_atac_listf[[t]])))
regulons[,extended:=str_detect(tf,"e$")]
regulons[,tf:=str_remove(tf,"e$")]
regulons[(extended)] #25397 tf> target interaction
regulons[(!extended)] #4808 tf> target interaction with high confidence
fwrite(regulons,fp(out,"tf_target_interactions.csv"))
regulons<-fread(fp(out,"tf_target_interactions.csv"))

#start build network only with tf> interact with high conf
regf<-regulons[(!extended)]

length(unique(regf$tf)) #72 tfs
regf[,n.target:=.N,by="tf"]
summary(unique(regf,by="tf")$n.target)
   # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 4.00   16.50   32.50   66.78   75.50  462.00

#show network using ggnet
#renv::install("briatte/ggnet")
library(ggnet)
library(network)
library(sna)

?network

regf<-regf[!is.na(target)]
net<-as.network(regf[,.(tf,target)],loops = T,directed = T)
net
 # Network attributes:
 #  vertices = 1802 
 #  directed = TRUE 
 #  hyper = FALSE 
 #  loops = TRUE 
 #  multiple = FALSE 
 #  bipartite = FALSE 
 #  total edges= 4808 
 #    missing edges= 0 
 #    non-missing edges= 4808 
 # 
 # Vertex attribute names: 
 #    vertex.names 
 # 
 # Edge attribute names not shown 

saveRDS(net,fp(out,"network_tf_target_hi_conf.rds"))

#only with tf of interest
egr1_modul<-c("KLF2","EGR1","KLF4")
reg_egr1<-regf[tf%in%c(egr1_modul)] #add only targets of the tfs altered

fwrite(reg_egr1,fp(out,"egr1_KLF2_KLF4_network_tf_target_interactions.csv"))

net_genes<-union(reg_egr1$tf,reg_egr1$target)
reg_egr1r1<-unique(rbind(reg_egr1,regf[target%in%egr1_modul&tf%in%net_genes])) #add also regulators of this tfs in this newtwork
fwrite(reg_egr1r1,fp(out,"egr1_network_plus_tf_regulators_tf_target_interactions.csv"))

# reg_egr1r2<-regf[tf%in%c(egr1_modul)|target%in%egr1_modul] #add upstream regulators of egr1_modul
# fwrite(reg_egr1r2,fp(out,"egr1_network_plus_tf_regulators.csv"))

#reg_egr1<-unique(rbind(reg_egr1,regf[tf%in%net_genes&target%in%net_genes])) #add all interactions of this genes presents
tfs<-unique(reg_egr1r1$tf)
net_egr1<-as.network(reg_egr1r1[,.(tf,target)],loops = T,directed = T)
net_egr1
#  Network attributes:
#   vertices = 123 
#   directed = TRUE 
#   hyper = FALSE 
#   loops = TRUE 
#   multiple = FALSE 
#   bipartite = FALSE 
#   total edges= 161 
#     missing edges= 0 
#     non-missing edges= 161 
# 
#  Vertex attribute names: 
#     vertex.names 
# 
# No edge attributes

#add a vertex attributes wich indicates if the gene is a tf or not
net_egr1 %v% "type" = ifelse(network.vertex.names(net_egr1) %in% regf$tf, "tf", "gene")

#add methyl info
res_m<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
res_m[gene_score_add>500,meth:="cadetblue"]
res_m[gene_score_add<=500,meth:="cornsilk3"]

net_egr1 %v% "meth" = sapply(res_m[network.vertex.names(net_egr1),on="gene"]$meth,function(x)ifelse(is.na(x),"cornsilk3",x))


#add expr info 
res_e<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]

res_e[padj>0.05,deg:="black"]
res_e[padj<=0.05&log2FoldChange>0.5,deg:="red"]
res_e[padj<=0.05&log2FoldChange<(-0.5),deg:="blue"]

net_egr1 %v% "deg" = res_e[network.vertex.names(net_egr1),on="gene"]$deg


#add atac info 
#on vertice
#need add target info
res_a<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_hsc_logFC0.csv.gz")
res_a<-res_a[!str_detect(peak,"chr[XY]")]
peaks_hsc_genes[,peak:=query_region]
res_at<-merge(res_a,peaks_hsc_genes,by="peak")
res_at[,target:=gene_name]

res_at[p_val_adj<0.001&avg_log2FC>0.25,da:="red"]
res_at[p_val_adj<0.001&avg_log2FC<(-0.25),da:="blue"]
res_at[is.na(da),da:="black"]

net_egr1 %v% "da" = res_at[network.vertex.names(net_egr1),on="target"]$da

#on edge : df tf-target link with if peak with motif found, FC / pval of the change
#need merge network df with res_atac df
#add TF info on res_atac => for each peak, merge with tf(of the network)-peak dt 
tfs<-unique(reg_egr1r1[,.(tf,target)]$tf) #"KLF4" "EGR1" "KLF2" "ATF4" "ATF3" "JUN"  "FOS"  "JUNB"
peaks<-unique(res_at[target%in%reg_egr1r1$target]$peak)
length(peaks)#690

motif.all <- GetMotifData(
    object = atacs, assay = "lin_peaks", slot = "data"
  )

motifs_peaks_tfs <- motif.all[peaks,GetMotifIDs(atacs,tfs) , drop = FALSE]
tf_peak_dt<-melt(data.table(data.frame(as.matrix(motifs_peaks_tfs==1)),keep.rownames = "peak"),id.vars = "peak",variable.name ="motif.id" ,value.name = "is_present")
tf_peak_dt<-merge(tf_peak_dt,GetMotifIDs(atacs,tfs,return_dt=TRUE))
tf_peak_dt<-tf_peak_dt[is_present==TRUE]

res_at_tf<-merge(res_at,tf_peak_dt,by="peak")
res_at_tf[,tf:=motif.name]

#merge with network df
reg_egr1r1_peaks<-merge(reg_egr1r1,
                        res_at_tf[,.(tf,target,peak,p_val,p_val_adj,avg_log2FC,pct.1,pct.2,type,da)],
                        by = c("tf","target"),
                        all.x = T)

reg_egr1r1_peaks[,n.tf.target.peaks:=.N,by=.(tf,target)]
reg_egr1r1_peaks[,biggest_change:=p_val==min(p_val),.(tf,target)]
reg_egr1r1_peaks[(biggest_change)]
unique(reg_egr1r1_peaks[(biggest_change)],by=c("tf","target"))#ok


reg_egr1r1_peaks[,da.peak:=p_val_adj<0.001&abs(avg_log2FC)>0.25]

reg_egr1r1_peaks[,n.da.tf.target.peaks:=sum(da.peak),.(tf,target)]

reg_egr1r1_peaks[da.peak==T] #9

reg_egr1r1_peak1<-reg_egr1r1_peaks[(biggest_change)]

#ADD edge atttibute (tf> target) : color depend of if atac based tf> target interact is altered by chromatin change
#if tf-gene peak dn : blue if tf-gene peak up :  red

net_egr1_a<-network(reg_egr1r1_peak1[,-c("extended","biggest_change")],loops = T,directed = T)
list.edge.attributes(net_egr1_a)
as.matrix(net_egr1_a,attrname='da')

set.edge.attribute(net_egr1, "color", ifelse(net_egr1 %e% "dap" > 1, "black", "grey75"))
ggnet2(bip, shape = "mode", edge.label = "weights", edge.label.color = "color")


#genes_of_interest<-union(res_e[padj<=0.05&abs(log2FoldChange)>0.5]$gene,union(res_m[gene_score_add>500]$gene,unique(reg_egr1$tf)))


#GRN sans selection, label only fgenes of interest
ggnet2(net_egr1,
       color = "meth",
       label = genes_of_interest,label.color = "deg",label.size = 3,
       size = "type" ,size.palette = c("tf"=3,"gene"=1),
       shape = "type",
       edge.alpha = 0.6,
       arrow.size = 5, arrow.gap =0.01) +
  theme(panel.background = element_rect(fill = "white"))
ggsave(fp(out,"all_KLF2_FOSB_JUN_EGR1_KLF4_ARID5A_network_label_only_altered_regulons_genes_and_TF_upstream.pdf"),width = 15,height = 10)

ggnet2(net_egr1,
       color = "meth",
       label = T,label.color = "deg",label.size = 3,
       size = "degree" ,size.min = 2,size.cut = 4,
       shape = "type",
       edge.alpha = 0.6,
       arrow.size = 5, arrow.gap =0.02) +
  theme(panel.background = element_rect(fill = "white"))
ggsave(fp(out,"all_KLF2_FOSB_JUN_EGR1_KLF4_ARID5A_network_all_label_min2Con.pdf"),width = 10,height = 6)

summary(res_m[gene_score_add>250])
res_anno<-fread("outputs/02-gene_score_calculation_and_validation/res_anno.csv.gz")

res_anno[,ncpg.region:=.N,by=.(region_type,gene)]
res_anno[,ncpg.prom:=sum(region_type=="promoter"),by="gene"]
res_anno[,ncpg.enh:=sum(region_type=="other"),by="gene"]

resg<-unique(res_anno[order(gene,pval)],by="gene")
resg[gene_score_add>700&pval>0.01]
plot(density(resg$gene_score_add))

#with only egr1 modul
egr1_modul<-c("KLF2","FOSB","JUN","EGR1","KLF4","ARID5A","KLF10","JUNB")
reg_egr1_o<-regf[tf%in%egr1_modul&target%in%egr1_modul] 

net_egr1o<-as.network(reg_egr1_o[,.(tf,target)],loops = T,directed = T)


ggnet2(net_egr1o,label=T,
       edge.alpha = 0.6,
       arrow.size = 5, arrow.gap =0.02) +
  theme(panel.background = element_rect(fill = "white"))

ggsave(fp(out,"network_altered_regulons_and_genes.pdf"))

ggnet2(net_egr1f,
       color = "meth",
       label = T,label.color = "deg",label.size = 3,
       size = "degree",size.min = 2,size.cut = 4,
       shape = "type" ,
       edge.alpha = 0.7,
       arrow.size = 5, arrow.gap =0.02) +
  theme(panel.background = element_rect(fill = "white"))


ggsave(fp(out,"network_altered_regulons_and_genes_mincon2.pdf"))




### Complete =======================================================================================
message("Success!", appendLF = TRUE)

