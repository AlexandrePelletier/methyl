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

#% EGR1 in all peaks ?
egr1_in_peaks<-CheckMotif(atacs,
                       peaks =rownames(atacs) ,
                       motif.name = "EGR1",
                       return.peaks = FALSE)
sum(egr1_in_peaks)/nrow(atacs) #12%

#% EGR1 in HSC peaks ?
length(intersect(peaks_hsc_genes$query_region,rownames(atacs))) / nrow(atacs) #56% of all peaks are in HSC peaks
egr1_in_hscpeaks<-CheckMotif(atacs,
                       peaks =peaks_hsc_genes$query_region ,
                       motif.name = "EGR1",
                       return.peaks = FALSE)
sum(egr1_in_hscpeaks)/length(unique(peaks_hsc_genes$query_region)) #19% 


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

#for KLF2

peaks_close_KLF2_target<-peaks_hsc_genes[gene_name%in%regulons_list$KLF2]$query_region
klf2_peaks<-CheckMotif(atacs,
                       peaks =peaks_close_KLF2_target ,
                       motif.name = "KLF2",
                       return.peaks = TRUE)
length(klf2_peaks)/length(peaks_close_KLF2_target) #32%
klf2_genes<-intersect(peaks_hsc_genes[query_region%in%klf2_peaks]$gene_name,
                      regulons_list$KLF2)
length(klf2_genes)/length(regulons_list$KLF2) #90% !

#% KLF2 in all peaks ?
klf2_in_peaks<-CheckMotif(atacs,
                       peaks =rownames(atacs) ,
                       motif.name = "KLF2",
                       return.peaks = FALSE)
sum(klf2_in_peaks)/nrow(atacs) #17%

#% KLF2 in HSC peaks ?
klf2_in_hscpeaks<-CheckMotif(atacs,
                       peaks =peaks_hsc_genes$query_region ,
                       motif.name = "KLF2",
                       return.peaks = FALSE)
sum(klf2_in_hscpeaks)/length(unique(peaks_hsc_genes$query_region)) #24% 

#for KLF4

peaks_close_KLF4_target<-peaks_hsc_genes[gene_name%in%regulons_list$KLF4]$query_region
klf4_peaks<-CheckMotif(atacs,
                       peaks =peaks_close_KLF4_target ,
                       motif.name = "KLF4",
                       return.peaks = TRUE)
length(klf4_peaks)/length(peaks_close_KLF4_target) #32%
klf4_genes<-intersect(peaks_hsc_genes[query_region%in%klf4_peaks]$gene_name,
                      regulons_list$KLF4)
length(klf4_genes)/length(regulons_list$KLF4) #95% !

#% KLF4 in all peaks ?
klf4_in_peaks<-CheckMotif(atacs,
                       peaks =rownames(atacs) ,
                       motif.name = "KLF4",
                       return.peaks = FALSE)
sum(klf4_in_peaks)/nrow(atacs) #24%

#% KLF2 in HSC peaks ?
klf4_in_hscpeaks<-CheckMotif(atacs,
                       peaks =peaks_hsc_genes$query_region ,
                       motif.name = "KLF4",
                       return.peaks = FALSE)
sum(klf4_in_hscpeaks)/length(unique(peaks_hsc_genes$query_region)) #32% 


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
regulons[(extended)] #25397 tf > target interaction
regulons[(!extended)] #4808 tf > target interaction with high confidence
fwrite(regulons,fp(out,"tf_target_interactions.csv"))
regulons<-fread(fp(out,"tf_target_interactions.csv"))

#%TF-target conserved
regulonsf<-fread(fp(out,"tf_target_interactions.csv"))[!(extended)]
regulons_old<-fread("outputs/10-SCENIC/regulons.csv")
res_conserved<-sapply(unique(regulons_old$tf), function(t)length(intersect(regulons_old[tf==t]$gene,regulonsf[tf==t]$target))/nrow(regulons_old[tf==t]))
res_conserved[c("EGR1","KLF2","KLF4")]
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
res_m[gene_score_add>500,meth:="darkred"]
res_m[gene_score_add<=500,meth:="black"]

net_egr1 %v% "meth" = sapply(res_m[network.vertex.names(net_egr1),on="gene"]$meth,function(x)ifelse(is.na(x),"cornsilk3",x))


#add expr info 
res_e<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]

res_e[padj>0.05,deg:="cornsilk3"]
res_e[padj<=0.05&log2FoldChange>0,deg:="coral2"]
res_e[padj<=0.05&log2FoldChange>0.5,deg:="coral3"]
res_e[padj<=0.05&log2FoldChange<(0),deg:="cadetblue3"]
res_e[padj<=0.05&log2FoldChange<(-0.5),deg:="cadetblue4"]
res_e[padj<=0.05&log2FoldChange<(-0.25)]
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
res_at[is.na(da),da:="grey75"]

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
reg_egr1r1_peaks[(biggest_change)|is.na(biggest_change)]
unique(reg_egr1r1_peaks[(biggest_change)|is.na(biggest_change)],by=c("tf","target"))#ok


reg_egr1r1_peaks[,da.peak:=p_val_adj<0.001&abs(avg_log2FC)>0.25]

reg_egr1r1_peaks[,n.da.tf.target.peaks:=sum(da.peak),.(tf,target)]

reg_egr1r1_peaks[da.peak==T] #9

#reg_egr1r1_peak1<-reg_egr1r1_peaks[(biggest_change)|is.na(biggest_change)]

#add DMCs infos on edge
#need merge peaks DMCs df with reg_egr1 df
peaks_cpgs<-fread("outputs/14-DMCs_atac_integr/cpgs_in_lin_OCRs.csv.gz")
peaks_meth<-merge(peaks_cpgs,fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz"))
peaks_meth[,peak:=peaks]
peaks_meth_hsc<-peaks_meth[peak%in%peaks_hsc_genes$peak]

reg_egr1r1_peak_meth<-merge(reg_egr1r1_peaks,
                             peaks_meth_hsc[,.(peak,cpg_id,logFC,P.Value,adj.P.Val)],
                             by="peak",
                             all.x=T)
reg_egr1r1_peak_meth[,n.cpg.peak:=.N,by=.(peak,tf)]
reg_egr1r1_peak_meth[,biggest_meth_change:=P.Value==min(P.Value),.(peak,tf)]
reg_egr1r1_peak_meth[(biggest_meth_change)|is.na(biggest_meth_change)] #ok


reg_egr1r1_peak_meth[,dmcs:=P.Value<0.001&abs(logFC)>25]

reg_egr1r1_peak_meth[,n.dmcs.peak:=sum(dmcs),.(peak,tf)]

reg_egr1r1_peak_meth[dmcs==T] #36

reg_egr1r1_peak1_meth1<-reg_egr1r1_peak_meth[biggest_change|is.na(biggest_change)|(biggest_meth_change)|is.na(biggest_meth_change)]
reg_egr1r1_peak1_meth1[,da:=unique(da[biggest_change==T]),.(tf,target)]
reg_egr1r1_peak1_meth1[,n.dmcs:=sum(unique(n.dmcs.peak[biggest_meth_change==T]),na.rm = T),.(tf,target)]

reg1<-unique(reg_egr1r1_peak1_meth1[,-c("extended","biggest_change","biggest_meth_change","peak")])
reg1<-unique(reg1,by=c('tf',"target"))
reg1<-reg1[!is.na(da)]

#ADD edge atttibute (tf> target) : color depend of if atac based tf> target interact is altered by chromatin change
#if tf-gene peak dn : blue if tf-gene peak up :  red

net_egr1_a<-network(reg1,loops = T,directed = T)
list.edge.attributes(net_egr1_a)
as.matrix(net_egr1_a,attrname='da')
net_egr1_a %e% "da"
net_egr1_a %e% "da"=sapply(net_egr1_a %e% "da",function(x)ifelse(x=="grey75","darkgrey",x))

net_egr1_a %e% "dmc_line"=net_egr1_a %e% "n.dmcs"+1
net_egr1_a %e% "dmc_line"=sapply(net_egr1_a %e% "dmc_line",function(x)ifelse(is.na(x),1,x))

#set.edge.attribute(net_egr1, "color", ifelse(net_egr1 %e% "dap" > 1, "black", "grey75"))
#add vertices attributes
net_egr1_a %v% "type" = ifelse(network.vertex.names(net_egr1_a) %in% regf$tf, "tf", "gene")
net_egr1_a %v% "deg" = res_e[network.vertex.names(net_egr1_a),on="gene"]$deg
net_egr1_a %v% "deg" = sapply(net_egr1_a %v% "deg",function(x)ifelse(is.na(x),"cornsilk3",x))

net_egr1_a %v% "meth" = sapply(res_m[network.vertex.names(net_egr1_a),on="gene"]$meth,function(x)ifelse(is.na(x),"black",x))

#genes_of_interest<-union(res_e[padj<=0.05&abs(log2FoldChange)>0.5]$gene,union(res_m[gene_score_add>500]$gene,unique(reg_egr1$tf)))


#GRN sans selection,label all genes
ggnet2(net_egr1_a,
       color = "deg",
       label = T,label.color = "meth",label.size = 2,
       size = "type" ,size.palette = c("tf"=3,"gene"=1),
       shape = "type",
       edge.alpha = 0.8,
       edge.size=0.5,
       edge.color = "da",
       arrow.size = 5,
       edge.lty = "dmc_line",
       arrow.gap =0.02) +
  theme(panel.background = element_rect(fill = "white"))

ggsave(fp(out,"final_network_EGR1_KLF2_KLF4_tf_targets_2.pdf"),width = 10,height = 10)
reg_egr1r1_peak1_meth1[n.dmcs.peak>1]

#compare network epigen alter with other network
#make complete tf - target epigen anno datafram
#1)chrine change regulons anno
res_a<-fread("outputs/15-chromatin_change_LGA_vs_Ctrl/differential_peaks_accessibility_lga_vs_ctrl_hsc_logFC0.csv.gz")
res_a<-res_a[!str_detect(peak,"chr[XY]")]
peaks_hsc_genes[,peak:=query_region]
res_at<-merge(res_a,peaks_hsc_genes,by="peak")
res_at[,target:=gene_name]


peaks<-unique(res_at[target%in%regf$target]$peak)
length(peaks)#12880

motif.all <- GetMotifData(
    object = atacs, assay = "lin_peaks", slot = "data"
  )

motifs_peaks_tfs <- motif.all[peaks,GetMotifIDs(atacs,unique(regf$tf)) , drop = FALSE]
tf_peak_dt<-melt(data.table(data.frame(as.matrix(motifs_peaks_tfs==1)),keep.rownames = "peak"),id.vars = "peak",variable.name ="motif.id" ,value.name = "is_present")
tf_peak_dt<-merge(tf_peak_dt,GetMotifIDs(atacs,unique(regf$tf),return_dt=TRUE))
tf_peak_dt<-tf_peak_dt[is_present==TRUE]

res_at_tf<-merge(res_at,tf_peak_dt,by="peak")
res_at_tf[,tf:=motif.name]

#merge with network df
reg_tf_peaks<-merge(regf,
                        res_at_tf[,.(tf,target,peak,p_val,p_val_adj,avg_log2FC,pct.1,pct.2,type)],
                        by = c("tf","target"),
                        all.x = T)

reg_tf_peaks[,n.tf.target.peaks:=.N,by=.(tf,target)]
reg_tf_peaks[,biggest_change:=p_val==min(p_val),.(tf,target)]
reg_tf_peaks[(biggest_change)|is.na(biggest_change)]
unique(reg_tf_peaks[(biggest_change)|is.na(biggest_change)],by=c("tf","target"))#ok


reg_tf_peaks[,da.peak:=p_val_adj<0.001&abs(avg_log2FC)>0.25]
reg_tf_peaks[,n.da.tf.target.peaks:=sum(da.peak),.(tf,target)]
reg_tf_peaks[da.peak==T] #59
fwrite(reg_tf_peaks,fp(out,"regulons_chrine_change_anno.csv"))

#2) meth change regulons infos
#need merge peaks DMCs df with regf df
peaks_cpgs<-fread("outputs/14-DMCs_atac_integr/cpgs_in_lin_OCRs.csv.gz")
peaks_meth<-merge(peaks_cpgs,fread("outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz"))
peaks_meth[,peak:=peaks]
peaks_meth_hsc<-peaks_meth[peak%in%peaks_hsc_genes$peak]

reg_tf_peaks_meth<-merge(reg_tf_peaks,
                             peaks_meth_hsc[,.(peak,cpg_id,logFC,P.Value,adj.P.Val)],
                             by="peak",
                             all.x=T)
reg_tf_peaks_meth[,n.cpg.peak:=.N,by=.(peak,tf)]
reg_tf_peaks_meth[,biggest_meth_change:=P.Value==min(P.Value),.(peak,tf)]
reg_tf_peaks_meth[(biggest_meth_change)|is.na(biggest_meth_change)] #ok


reg_tf_peaks_meth[,dmcs:=P.Value<0.001&abs(logFC)>25]

reg_tf_peaks_meth[,n.dmcs.peak:=sum(dmcs),.(peak,tf)]

reg_tf_peaks_meth[dmcs==T] #485
fwrite(reg_tf_peaks_meth,fp(out,"regulons_chrine__and_meth_change_anno.csv.gz"))
reg_tf_peaks_meth<-fread(fp(out,"regulons_chrine__and_meth_change_anno.csv.gz"))

#% genes of EGR1 / KLF2 KLF4 network epigen altered
#for peak containing interactions
reg_tf_peaks_methf<-reg_tf_peaks_meth[!is.na(peak)]
pct.alter.netint<-length(unique(reg_tf_peaks_methf[tf%in%c("EGR1","KLF2","KLF4")&((da.peak)|(dmcs))]$target))/length(unique(reg_tf_peaks_methf[tf%in%c("EGR1","KLF2","KLF4")]$target))#22%
pct.alter.netint#23%

unique(reg_tf_peaks_methf[tf%in%c("EGR1","KLF2","KLF4")&((da.peak)|(dmcs))])
nonplotted<-setdiff(unique(reg_tf_peaks_methf[tf%in%c("EGR1","KLF2","KLF4")&((da.peak)|(dmcs))]$target),
         c("ARRDC2","KLF13","AHNAK","RHOC","HEXIM1","TSC22D1","INTS6","TINAGL1","CD151","PTGER4","DDIT","TOB1","ETS1",
           "SOCS3","JUNB","CDKN1C","PRNP","GADD45B","LMNA","MYADM","DDIT4"))
nonplotted#"ID3"     "UBC"     "EIF5"    "TGIF1"   "NFATC1"  "PHACTR2" "IFRD1"  
reg_tf_peaks_methf[tf%in%tfsint&target%in%nonplotted&((da.peak)|(dmcs))]
#compared to randomly peaks 3 TFs
rand_alters<-sapply(1:1000,function(i){
  set.seed(i)
  tfs<-sample(unique(reg_tf_peaks_methf$tf),3)
  pct.alter<-length(unique(reg_tf_peaks_methf[tf%in%tfs&((da.peak)|(dmcs))]$target))/length(unique(reg_tf_peaks_methf[tf%in%tfs]$target))
  return(pct.alter)
  })
mean(rand_alters) #10%
median(rand_alters) #9.6%

sum(pct.alter.netint<rand_alters) #0.5% (p=0.005)

#compared to randomly peaks n tf-gene pair
tfsint<-c("EGR1","KLF2","KLF4")
tf_genes<-unique(reg_tf_peaks_methf[,.(tf,target)])
ntfgenes_int<-nrow(unique(reg_tf_peaks_methf[tf%in%tfsint][,.(tf,target)])) #154
rand_alters<-sapply(1:1000,function(i){
  set.seed(i)
  tf_genes_s<-tf_genes[sample(1:.N,ntfgenes_int)]
  reg_tf_peaks_meth_s<-merge(reg_tf_peaks_methf,tf_genes_s,by=c("tf","target"))
  pct.alter<-nrow(unique(reg_tf_peaks_meth_s[((da.peak)|(dmcs))][,.(tf,target)]))/ntfgenes_int
  return(pct.alter)
  })
mean(rand_alters) #10%
median(rand_alters) #9%

sum(pct.alter.netint<rand_alters) #0% (p<0.001)

#same for DEGs
#need add DEGs infos
res_e<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]
res_e[,degs:=padj<0.05&abs(log2FoldChange)>0.5]
res_e[,down:=padj<0.05&log2FoldChange<0]

res_e[,target:=gene]

reg_tf_peaks_meth_expr<-merge(reg_tf_peaks_meth,res_e,by="target")
fwrite(reg_tf_peaks_meth_expr,fp(out,"regulons_chrine_meth_and_expr_change_anno.csv.gz"))

reg_tf_peaks_meth_exprf<-reg_tf_peaks_meth_expr[!is.na(peak)]

pct.alter.netint<-length(unique(reg_tf_peaks_meth_exprf[tf%in%c("EGR1","KLF2","KLF4")&(down)]$target))/length(unique(reg_tf_peaks_meth_exprf[tf%in%c("EGR1","KLF2","KLF4")]$target))


pct.alter.netint#22%

#compared to randomly peaks 3 TFs
rand_alters<-sapply(1:1000,function(i){
  set.seed(i)
  tfs<-sample(unique(reg_tf_peaks_meth_exprf$tf),3)
  pct.alter<-length(unique(reg_tf_peaks_meth_exprf[tf%in%tfs&(down)]$target))/length(unique(reg_tf_peaks_meth_exprf[tf%in%tfs]$target))
  return(pct.alter)
  })
mean(rand_alters) #6%
median(rand_alters) #5%

sum(pct.alter.netint<rand_alters) #0.2% (p=0.002)

#compared to randomly peaks n tf-gene pair
tfsint<-c("EGR1","KLF2","KLF4")
tf_genes<-unique(reg_tf_peaks_meth_exprf[,.(tf,target)])
ntfgenes_int<-nrow(unique(reg_tf_peaks_meth_exprf[tf%in%tfsint][,.(tf,target)])) #154
rand_alters<-sapply(1:1000,function(i){
  set.seed(i)
  tf_genes_s<-tf_genes[sample(1:.N,ntfgenes_int)]
  reg_tf_peaks_meth_exprf_s<-merge(reg_tf_peaks_meth_exprf,tf_genes_s,by=c("tf","target"))
  pct.alter<-nrow(unique(reg_tf_peaks_meth_exprf_s[(down)][,.(tf,target)]))/ntfgenes_int
  return(pct.alter)
  })
mean(rand_alters) #5%
median(rand_alters) #5%

sum(pct.alter.netint<rand_alters) #0% (p<0.001)


