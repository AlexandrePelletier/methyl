### Project Setup ==================================================================================
library(here)
out <- here("outputs", "13-GRN_integr")
dir.create(out, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  source("scripts/utils/new_utils.R")
})


### Tables and Figures Theme =======================================================================
theme_set(theme_light())


### Functions ======================================================================================


### Analysis =======================================================================================
regulons_list<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/regulons_list.rds")

#make a df of interactions tf > targets
regulons<-Reduce(rbind,lapply(names(regulons_list), function(t)data.table(tf=rep(t,length(regulons_list[[t]])),target=regulons_list[[t]])))
regulons[,extended:=str_detect(tf,"e$")]
regulons[,tf:=str_remove(tf,"e$")]
regulons[(extended)] #59955 tf> target interaction
regulons[(!extended)] #12017 tf> target interaction with high confidence
fwrite(regulons,fp(out,"tf_target_interactions.csv"))
regulons<-fread(fp(out,"tf_target_interactions.csv"))

#start build netork only with tf> interact with high conf
regf<-regulons[(!extended)]

length(unique(regf$tf)) #106 tfs
regf[,n.target:=.N,by="tf"]
summary(unique(regf,by="tf")$n.target)
  #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 10.00   20.25   52.50  113.37  130.00  896.00 

#show network using ggnet
renv::install("briatte/ggnet")
library(ggnet)
library(network)
library(sna)

?network

regf<-regf[!is.na(target)]
net<-as.network(regf[,.(tf,target)],loops = T,directed = T)
net
 # Network attributes:
 #  vertices = 2477 
 #  directed = TRUE 
 #  hyper = FALSE 
 #  loops = TRUE 
 #  multiple = FALSE 
 #  bipartite = FALSE 
 #  total edges= 12017 
 #    missing edges= 0 
 #    non-missing edges= 12017 
 # 
 # Vertex attribute names: 
 #    vertex.names 
 # 
 # Edge attribute names not shown 

saveRDS(net,fp(out,"network_tf_target_hi_conf.rds"))

#only with tf of interest
egr1_modul<-c("KLF2","FOSB","JUN","EGR1","KLF4","ARID5A")
reg_egr1<-regf[tf%in%c(egr1_modul)] #add only targets of the tfs altered

fwrite(reg_egr1,fp(out,"egr1_network_tf_target_interactions.csv"))

net_genes<-union(reg_egr1$tf,reg_egr1$target)
reg_egr1r<-unique(rbind(reg_egr1,regf[target%in%egr1_modul&tf%in%net_genes])) #add also regulators of this tfs in this newtwork
fwrite(reg_egr1r,fp(out,"egr1_network_plus_tf_regulators_tf_target_interactions.csv"))

#reg_egr1<-regf[tf%in%c(egr1_modul)|target%in%egr1_modul] #add upstream regulators of egr1_modul

#reg_egr1<-unique(rbind(reg_egr1,regf[tf%in%net_genes&target%in%net_genes])) #add all interactions of this genes presents
tfs<-unique(reg_egr1$tf)
net_egr1<-as.network(reg_egr1[,.(tf,target)],loops = T,directed = T)
net_egr1
 # Network attributes:
 #  vertices = 208 
 #  directed = TRUE 
 #  hyper = FALSE 
 #  loops = TRUE 
 #  multiple = FALSE 
 #  bipartite = FALSE 
 #  total edges= 382 
 #    missing edges= 0 
 #    non-missing edges= 382 

#add a vertex attributes wich indicates if the gene is a tf or not
net_egr1 %v% "type" = ifelse(network.vertex.names(net_egr1) %in% regf$tf, "tf", "gene")

#add methyl / expression info
res_m<-fread("outputs/02-gene_score_calculation_and_validation/res_genes.csv.gz")
res_e<-fread("outputs/09-LGA_vs_Ctrl_Activated/res_pseudobulkDESeq2_by_lineage.csv.gz")[lineage=="HSC"]

res_e[padj>0.05,deg:="black"]
res_e[padj<=0.05&log2FoldChange>0.5,deg:="red"]
res_e[padj<=0.05&log2FoldChange<(-0.5),deg:="blue"]

net_egr1 %v% "deg" = res_e[network.vertex.names(net_egr1),on="gene"]$deg

res_m[gene_score_add>500,meth:="cadetblue"]
res_m[gene_score_add<=500,meth:="cornsilk3"]

net_egr1 %v% "meth" = sapply(res_m[network.vertex.names(net_egr1),on="gene"]$meth,function(x)ifelse(is.na(x),"cornsilk3",x))
genes_of_interest<-union(res_e[padj<=0.05&abs(log2FoldChange)>0.5]$gene,union(res_m[gene_score_add>500]$gene,unique(reg_egr1$tf)))


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

