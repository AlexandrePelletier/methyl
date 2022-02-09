#stimulation response charact at regulon level
#do a heatmap with col = regulon line = lineage or gene set

source("scripts/utils/new_utils.R")
out<-"outputs/08-HTO_signature/regulon_level"
dir.create(out)


regulons<-fread("outputs/10-SCENIC/regulons.csv")
regn<-unique(regulons$tf)
length(regn) #106
#regulon order ~ gene overlap
overlap_reg<-sapply(names(regulons_list),function(reg){
  reg1<-regulons_list[[reg]]
  return(sapply(regulons_list,function(reg2)length(intersect(reg1,reg2))/length(union(reg1,reg2))))
  })
clus<-hclust(as.dist(1-overlap_reg))
reg_order<-rownames(overlap_reg)[clus$order]

regul_act<-data.table(data.frame(as.matrix(readRDS("outputs/10-SCENIC/TF_AUC_assay.rds")@data)),
                      keep.rownames = "regulon")

regul_act<-melt(regul_act,id.vars = "regulon",variable.name = "bc",value.name = "activity")
mtd<-fread("outputs/06-integr_singlecell_cbps/metadata_cbps_filtered.csv.gz")                   
regul_act_mtd<-merge(regul_act,mtd[,bc:=str_replace(bc,"-",'.')])
fwrite(regul_act_mtd,fp(out,"regul_activity_by_cells_mtd.csv.gz"))
regul_act_mtd<-fread(fp(out,"regul_activity_by_cells_mtd.csv.gz"))

#moy acti by lineage for regulon non extended
regulf_act<-regul_act_mtd[regulon%in%regn]
remove(regul_act,regul_act_mtd)
regulf_act[,activity.lin:=mean(activity),.(regulon,lineage_hmap)]
rlin<-unique(regulf_act,by=c("regulon","lineage_hmap"))
rlinf<-rlin[differentiated==F]

ml<-as.matrix(data.frame(dcast(rlinf[,.(regulon,lineage_hmap,activity.lin)],formula = regulon~lineage_hmap),
                         row.names = "regulon"))
clus<-hclust(dist(ml))

# reg_order<-rownames(ml)[clus$order]
lin_order<-c("LT.HSC","HSC","MPP.LMPP","Erythro.Mas","Myeloid","Lymphoid")

break_cols_up<-c(1:10*3/100)
pheatmap::pheatmap(t(ml)[lin_order,reg_order],
                   cluster_rows = F,cluster_cols = F,
                   color = colorRampPalette(c("white", "red"))(length(break_cols_up)+1),
                   breaks =break_cols_up)



#regulon activation by lineage
regul_lin<-fread("outputs/10-SCENIC/regulon_activity_HTO_vs_not_by_lineage.csv.gz")
regul_lin<-regul_lin[!str_detect(regulon,'e$')]
unique(regul_lin$regulon) #106

#redo regulon activation calcul
regul_act_mtd<-fread(fp(out,"regul_activity_by_cells_mtd.csv.gz"))
regul_act_mtdf<-regul_act_mtd[sample%in%c("ctrlM555","ctrlM518","ctrlM537")]
#regul_act_mtdf[,avg_log2FC:=]
plot(density(scale(exp(regul_act_mtdf$activity))))
#regul_lin[,activation_score:=sign(avg_log2FC)*-log10(p_val_adj)]
#regul_lin[activation_score==Inf,activation_score:=sign(activation_score)*303]

SCENIC::runSCENIC_3_scoreCells
AUCell::AUCell_calcAUC


regul_lin[,activation_score:=ifelse(p_val_adj<0.05,avg_log2FC,0)]
#regul_lin[,activation_score_scaled:=scale(activation_score,center = F),by="lineage"]

regul_mat<-dcast(regul_lin[,.(regulon,lineage,activation_score)],formula = regulon~lineage)
regul_mat<-regul_mat[,.(regulon,`LT-HSC`,HSC,`MPP/LMPP`,`Erythro-Mas`,Myeloid,Lymphoid)]

regul_mat<-t(as.matrix(data.frame(regul_mat,row.names = "regulon")))
# clus<-hclust(dist(t(regul_mat)))
# reg_order<-colnames(regul_mat)[clus$order]

rownames(regul_mat)<-paste0(rownames(regul_mat),"_activation")
#regul_mat_o<-regul_mat[,order(regul_mat["HSC",],decreasing = T)]
break_cols_up<-c(1:100*2/1000)
break_cols<-sort(c(-(break_cols_up),break_cols_up))
pheatmap::pheatmap(regul_mat,
                   cluster_rows = F,cluster_cols = T,
                   show_colnames = T,
                  color = colorRampPalette(c("darkblue","white", "red"))(length(break_cols)+1),
         breaks =break_cols)  


#combine the 2 heatmap
library(ComplexHeatmap)
# first heatmap 
library(circlize)
col_fun1 = colorRamp2(c(0, 0.4), c("white" ,"darkblue"))
#col_fun1(seq(0,0.3))
h1 <- Heatmap(t(ml)[lin_order,reg_order], 
              col = col_fun1,
              cluster_rows =FALSE,
              cluster_columns=FALSE,
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              name="regulon activity",
              column_title = "Regulon",
              row_title = "Lineage")

# second heatmap
col_fun2 = colorRamp2(c(-0.15,0, 0.15), c("chartreuse3","white", "darkorange2"))
h2 <-   Heatmap(regul_mat[paste0(lin_order,"_activation"),reg_order], 
                col=col_fun2,
                cluster_rows = FALSE,
                cluster_columns=FALSE,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                name="regulon activation",
                column_title = "Regulon",
                row_title = "Lineage Activation")

h1 %v% h2


# gene set to test enrichment for regulon

#goprofiler
renv::install("gprofiler2")
library(gprofiler2)
regulons_list<-split(regulons$gene,regulons$tf)
reg_enrich<-gost(query = regulons_list,significant = F,organism = "hsapiens",correction_method = "gSCS")
saveRDS(reg_enrich,fp(out,"reg_enrich_gprofil.rds"))
head(reg_enrich)
head(reg_enrich$result)

re<-data.table(reg_enrich$result)
split(re[significant==T&source=="GO:BP"]$term_name,re[significant==T&source=="GO:BP"]$query)

re[query=="EGR1"&significant==T&source=="GO:BP"]
re[query=="JUNB"&significant==T&source=="GO:BP"]$term_name

terms_of_interest<-c("cellular response to cytokine stimulus"   ,
                     "cellular response to chemical stimulus",
                     "cell activation" ,
                     "cellular response to external stimulus" ,
                     "response to stress"    ,
                     "regulation of cellular response to stress" ,
                     "regulation of response to stimulus" ,
                     "cellular response to stress",
                     "regulation of cellular response to stress",
                     "cellular response to growth factor stimulus",
                     "regulation of cellular response to growth factor stimulus" ,
                     "cellular response to extracellular stimulus" ,
                     "regulation of response to external stimulus" ,
                     "negative regulation of response to stimulus" ,
                     "positive regulation of response to cytokine stimulus" ,
                     "positive regulation of response to stimulus",
                     "regulation of stress-activated MAPK cascade"  ,
                     "inactivation of MAPK activity" ,
                     "regulation of cell activation",
                     "negative regulation of cell activation" ,
                     "positive regulation of cell activation",
                     "leukocyte activation" ,
                     "lymphocyte activation" ,
                      "regulation of cell differentiation",
                     "positive regulation of cell differentiation" ,
                     "negative regulation of cell differentiation",
                     "cell differentiation"  ,
                     "regulation of cell differentiation" ,
                     "stem cell differentiation",
                     "hemopoiesis"  ,
                     "embryonic hemopoiesis" ,
                     "primitive hemopoiesis"  ,
                     "positive regulation of hemopoiesis",
                     "regulation of hemopoiesis" ,
                     "negative regulation of hemopoiesis" ,
                     "hematopoietic progenitor cell differentiation"   ,
                     "lymphoid progenitor cell differentiation" ,
                     "regulation of lymphocyte differentiation" ,
                     "positive regulation of lymphocyte differentiation"  ,
                     "negative regulation of erythrocyte differentiation",
                     "positive regulation of erythrocyte differentiation" ,
                     "positive regulation of myeloid cell differentiation",
                     "negative regulation of myeloid cell differentiation" ,
                     "cell population proliferation",
                     "regulation of cell population proliferation",
                     "negative regulation of cell population proliferation" ,
                     "cell cycle"  ,
                     "regulation of cell cycle",
                     "mitotic cell cycle" ,
                     "regulation of mitotic cell cycle" ,
                     "negative regulation of cell cycle",
                     "negative regulation of cell cycle process" ,
                     "positive regulation of cell cycle",
                     "positive regulation of cell cycle process"  )
ref<-re[term_name%in%terms_of_interest]

fwrite(ref,fp(out,"res_regulons_enrich_for_BP_of_interest.csv"))
ref<-fread(fp(out,"res_regulons_enrich_for_BP_of_interest.csv"))


#hsc quiescence prolif geneset (from https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0020301)
hsc_dir<-"ref/hsc_signatures/hsc_act_sign/"
hsc_gs<-lapply(list.files(hsc_dir),function(f)fread(fp(hsc_dir,f),
                                                   select = c(2,5:7),
                                                   col.names = c("gene_mouse","log2FC","day","pval")))
names(hsc_gs)<-str_remove(list.files(hsc_dir),".csv")

hsc_gs<-Reduce(rbind,lapply(names(hsc_gs), function(x)hsc_gs[[x]][,gene_set:=x]))

gene_trad<-convertMouseGeneList(unique(hsc_gs$gene_mouse))
gene_trad[,gene_mouse:=MGI.symbol][,gene:=HGNC.symbol]
hsc_gs<-merge(hsc_gs,gene_trad[,.(gene_mouse,gene)])

fwrite(hsc_gs,fp(out,"hsc_state_signatures_geneset.csv"))

g<-readRDS("../singlecell/outputs/05-SCENIC/cbps0-8_clean/int/1.1_genesKept.Rds")

hsc_list<-split(hsc_gs$gene,hsc_gs$gene_set)
res_hsc<-OR2(query = regulons_list,
              terms_list = hsc_list,
              size_universe = length(g))

fwrite(res_hsc,fp(out,"res_regulons_enrich_hsc_state_signatures_geneset.csv"))
res_hsc<-fread(fp(out,"res_regulons_enrich_hsc_state_signatures_geneset.csv"))

#degs_hto
degs<-fread("outputs/08-HTO_signature/by_lineage/res_pseudobulk_DESeq2_3replicates.csv.gz")

degs<-degs[padj<0.05&abs(log2FoldChange)>0.5&lineage%in%c("LT-HSC","HSC","MPP/LMPP","Erythro-Mas","Myeloid","Lymphoid")]
degs[log2FoldChange>0,gene_set:=paste0(lineage,"_up")]
degs[log2FoldChange<0,gene_set:=paste0(lineage,"_dn")]

degs_list<-split(degs$gene,degs$gene_set)

res_degs<-OR2(query = regulons_list,
              terms_list = degs_list,
              size_universe = length(g))

fwrite(res_degs,fp(out,"res_regulons_enrich_degs_hto_by_lineage.csv"))
res_degs<-fread(fp(out,"res_regulons_enrich_degs_hto_by_lineage.csv"))

#rbind mat
res_bp<-fread(fp(out,"res_regulons_enrich_for_BP_of_interest.csv"),select = c(1,11,3),col.names = c("regulon","term","padj"))
res_hsc<-fread(fp(out,"res_regulons_enrich_hsc_state_signatures_geneset.csv"),select = c(9,1,8),col.names = c("regulon","term","padj"))
res_degs<-fread(fp(out,"res_regulons_enrich_degs_hto_by_lineage.csv"),select=c(9,1,8),col.names = c("regulon","term","padj"))

enrich_mat<-as.matrix(data.frame(dcast(Reduce(rbind,list(res_degs,res_bp,res_hsc)),formula = term~regulon),row.names = "term"))
enrich_mat<-(-log10(enrich_mat))



# third heatmap
#filter go bp
terms_of_interest<-c("cell activation" ,
                    "regulation of cell activation",
                     "cellular response to external stimulus" ,
                      "regulation of response to external stimulus" ,
                     "cellular response to stress",
                    "regulation of cellular response to stress" ,
                     "cellular response to growth factor stimulus",
                     "regulation of cellular response to growth factor stimulus" ,
                      "regulation of cell differentiation",
                     "stem cell differentiation",
                     "regulation of hemopoiesis" ,
                     "hematopoietic progenitor cell differentiation"   ,
                    "cell population proliferation",
                     "regulation of cell population proliferation",
                     "regulation of cell cycle" )

col_fun3 = colorRamp2(c(0,10, 50), c("white", "chartreuse4","black"))
h3 <-   Heatmap(enrich_mat[unique(res_degs$term),reg_order], 
                col=col_fun3,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                cluster_rows = TRUE,
                cluster_columns=FALSE,
                name="DEGs enrichment",
                column_title = "Regulon",
                row_title = "DEGs")

col_fun4 = colorRamp2(c(0,5, 25), c("white", "chartreuse4","black"))
Heatmap(enrich_mat[terms_of_interest,reg_order], 
                col=col_fun3,
        column_names_gp = grid::gpar(fontsize = 8),
  row_names_gp = grid::gpar(fontsize = 8),)

h4 <-   Heatmap(enrich_mat[terms_of_interest,reg_order], 
                col=col_fun4,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                cluster_rows = TRUE,
                cluster_columns=FALSE,
                name="GO BP enrichment",
                column_title = "Regulon",
                row_title = "Biological Process")



h5 <-   Heatmap(enrich_mat[unique(res_hsc$term),reg_order], 
                col=col_fun4,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                cluster_rows = TRUE,
                cluster_columns=FALSE,
                name="HSC signatures enrichment",
                column_title = "Regulon",
                row_title = "HSC state signature")


h4 %v% h5
#h1 %v% h2 %v% h3 %v% h4 %v% h5


min(res_degs$padj)
pdf(fp(out,"regulon_annotation.pdf"),height = 12)
 Heatmap(t(enrich_mat[c(unique(res_hsc$term),terms_of_interest),reg_order]), 
                col=col_fun4,
                row_names_gp = grid::gpar(fontsize = 8),
                column_names_gp = grid::gpar(fontsize = 8),
                cluster_rows = TRUE,
                cluster_columns=TRUE,
                name="Gene sets enrichment",
                row_title = "Regulon",
                column_title = "Gene sets")

dev.off()

pdf(fp(out,"regulon_activation.pdf"),width = 15)
h2 %v% h3
dev.off()
