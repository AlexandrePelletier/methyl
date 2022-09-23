#sirnas experiments

out<-"outputs/25-siRNAs_KLF2_cbl"
dir.create(out)
library(Seurat)
source("scripts/utils/new_utils.R")

#siCTRL
mat<-Read10X("~/RUN/Run_727_single-cell/output/count/humain/sict/outs/filtered_feature_bc_matrix/")

sict<-CreateSeuratObject(mat,project = "siCTRL_cbls")
sict #1919 cells
sict[["percent.mt"]]<-PercentageFeatureSet(sict,pattern = "MT-")
VlnPlot(object = sict, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

sict<-subset(sict,nFeature_RNA>2500&nCount_RNA>10000&percent.mt<8)
sict#1533 cells

#siKLF2
mat<-Read10X("~/RUN/Run_727_single-cell/output/count/humain/siklf2/outs/filtered_feature_bc_matrix/")

sikl<-CreateSeuratObject(mat,project = "siKLF2_cbls")
sikl #2836 cells
sikl[["percent.mt"]]<-PercentageFeatureSet(sikl,pattern = "MT-")
VlnPlot(object = sikl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

sikl<-subset(sikl,nFeature_RNA>2500&nCount_RNA>10000&percent.mt<8)
sikl#2116  cells

#merge
sict$sirna<-"CTRL"
sikl$sirna<-"KLF2"

sicbl<-merge(sict,sikl)

#norm
sicbl_list<-SplitObject(sicbl,split.by = "orig.ident")

sicbl_list<-lapply(sicbl_list,"SCTransform")

#KLF2 downregulation?
sicbl<-merge(sicbl_list[[1]],sicbl_list[[2]])
sicbl#3649 feature

VlnPlot(sicbl,features = "KLF2")

#clustering
VariableFeatures(sicbl)<-intersect(VariableFeatures(sicbl_list[[1]]),VariableFeatures(sicbl_list[[2]]))
sicbl<-RunPCA(sicbl)
sicbl<-FindNeighbors(sicbl,dims = 1:30)
sicbl<-RunUMAP(sicbl,dims = 1:30)
sicbl<-FindClusters(sicbl,resolution = 0.6)

sicbl <- CellCycleScoring(sicbl,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)
DimPlot(sicbl,label=T)#cluster ~ cell cycle

head(sicbl[[]])
sicbl_list<-SplitObject(sicbl,split.by = "orig.ident")

sicbl_list<-lapply(sicbl_list, SCTransform,vars.to.regress=c("S.Score","G2M.Score","percent.mt"))
sicbl<-merge(sicbl_list[[1]],sicbl_list[[2]])

#clustering
VariableFeatures(sicbl)<-intersect(VariableFeatures(sicbl_list[[1]]),VariableFeatures(sicbl_list[[2]]))
sicbl<-RunPCA(sicbl)
sicbl<-FindNeighbors(sicbl,dims = 1:30)
sicbl<-RunUMAP(sicbl,dims = 1:30)
DimPlot(sicbl,label=T)#ok
FeaturePlot(sicbl,features = "percent.mt")#ok

sicbl<-FindClusters(sicbl,resolution = 0.6)
DimPlot(sicbl,label=T)#
FeaturePlot(sicbl,features =c("GATA1","VPREB1","LTB","MPO","SOCS3","AVP"))#ok

#map on hmap
library(parallel)
hmap<-readRDS("outputs/05-make_hematomap/hematomap_ctrls_sans_stress.rds")

DefaultAssay(hmap)<-"integrated"
hmap[["pca.annoy.neighbors"]] <- LoadAnnoyIndex(object = hmap[["pca.annoy.neighbors"]], file = "outputs/05-make_hematomap/reftmp.idx")

sicbl_list<-lapply(sicbl_list,function(x){
    x<-SCTransform(x)
    x <- CellCycleScoring(x,s.features = cc.genes$s.genes,
                          g2m.features = cc.genes$g2m.genes,
                          set.ident = TRUE,
                          search=TRUE)


  
 # x$CC.Difference <- x$S.Score - x$G2M.Score
  return(x)
  },mc.cores = 2)

sicbl_list<-lapply(sicbl_list, SCTransform,vars.to.regress=c("percent.mt","S.Score","G2M.Score"),
                  return.only.var.genes=F)


anchors <- list()
for (i in 1:length(sicbl_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = hmap,
    query = sicbl_list[[i]],
    k.filter = NA,
    reference.reduction = "pca", 
    reference.neighbors = "pca.annoy.neighbors", 
    dims = 1:50
  )
}


for (i in 1:length(sicbl_list)) {
  sicbl_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = sicbl_list[[i]],
    reference = hmap, 
    refdata = list(
      cell_type = "cell_type", 
      lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "ref.umap"
  )
}


# Merge the queries
sicbl <- merge(sicbl_list[[1]], sicbl_list[2:length(sicbl_list)],merge.dr = c("ref.pca","ref.umap"))
DimPlot(sicbl, reduction = "ref.umap", group.by =  "predicted.cell_type", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()


DimPlot(sicbl, reduction = "ref.umap", group.by =  "predicted.lineage", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

DimPlot(sicbl, reduction = "ref.umap", group.by =  c("Phase","orig.ident")) 

#KLF2 reudction in MPP ?
VlnPlot(sicbl,"KLF2",group.by = "predicted.lineage")

VlnPlot(sicbl,"KLF2",group.by = "predicted.lineage")


#diff of lin ?
ggplot(sicbl@meta.data)+geom_bar(aes(x=sirna,fill=predicted.lineage),position = "fill")
ggplot(sicbl@meta.data)+geom_bar(aes(x=sirna,fill=predicted.cell_type),position = "fill")

#DEGs between sirnas ?
Idents(sicbl)<-"predicted.lineage"
degs<-FindMarkers(sicbl,group.by = "sirna",ident.1 = "KLF2",ident.2 = "CTRL",subset.ident = "MPP/LMPP")
degs<-data.table(degs,keep.rownames = "gene")
degs[p_val_adj<0.001] #762
degs_list<-split(degs[p_val_adj<0.001]$gene,degs[p_val_adj<0.001]$avg_log2FC>0)
names(degs_list)<-c("down","up")
lapply(degs_list, head)
# $down
# [1] "MTRNR2L12" "MIF"       "SELENOH"   "HNRNPH1"   "LRRC75A"   "RPS7"     
# 
# $up
# [1] "CD81"  "JDP2"  "MBD2"  "VEGFB"

#MBD2 : Methyl-CpG Binding Domain Protein, 
# Binds CpG islands in promoters where the DNA is methylated at position 5
#of cytosine within CpG dinucleotides. Binds hemimethylated DNA as well. 
#Recruits histone deacetylases and DNA methyltransferases.
# Acts as transcriptional repressor and plays a role in gene silencing. 
# Functions as a scaffold protein, targeting GATAD2A and GATAD2B to chromatin 
#to promote repression.
# May enhance the activation of some unmethylated cAMP-responsive promoters

#JPD2:Jun Dimerization Protein 2
#Component of the AP-1 transcription factor that represses transactivation mediated
#by the Jun family of proteins. 
# Involved in a variety of transcriptional responses associated with AP-1 
#such as UV-induced apoptosis, cell differentiation, tumorigenesis and antitumogeneris. 
# Can also function as a repressor by recruiting histone deacetylase 3/HDAC3 to the promoter region of JUN
#=> downreg of KLF2 => upreg of epigenetic machinery CpG Methylation binder / gene silencer + JUN family activity regulators
# +=> downreg of ++ genes  (due to MBD2 ?)
#which genes are downreg ?

library(gprofiler2)
res_gos<-gost(query = degs_list$down,organism = "hsapiens")
gprofiler2::gostplot(res_gos)
resgos<-data.table(res_gos$result)
lapply(split(resgos,by = "source"),head)

lapply(split(resgos,by = "source"),function(x)head(x[p_value<0.001]$term_name,20))

# $`GO:BP`
#  [1] "cellular macromolecule metabolic process"       
#  [2] "ATP metabolic process"                          
#  [3] "organonitrogen compound metabolic process"      
#  [4] "biosynthetic process"                           
#  [5] "regulation of mRNA metabolic process"           
#  [6] "organic substance biosynthetic process"         
#  [7] "cellular protein metabolic process"             
#  [8] "oxidative phosphorylation"                      
#  [9] "cellular biosynthetic process"                  
# [10] "regulation of cellular metabolic process"       
# [11] "RNA localization"                               
# [12] "cellular nitrogen compound biosynthetic process"
# [13] "peptide metabolic process"                      
# [14] "organonitrogen compound biosynthetic process"   
# [15] "cellular metabolic process"                     
# [16] "translation"                                    
# [17] "cellular localization"                          
# [18] "protein metabolic process"                      
# [19] "organelle organization"                         
# [20] "peptide biosynthetic process"                   
# 
                     
# 
# $`GO:MF`
#  [1] "protein binding"                                                
#  [2] "RNA binding"                                                    
#  [3] "mRNA binding"                                                   
#  [4] "cadherin binding"                                               
#  [5] "enzyme binding"                                                 
#  [6] "ubiquitin protein ligase binding"                               
#  [7] "ubiquitin-like protein ligase binding"                          
#  [8] "electron transfer activity"                                     
#  [9] "oxidoreduction-driven active transmembrane transporter activity"
# [10] "unfolded protein binding"                                       
# [11] "ribonucleoprotein complex binding"                              
# [12] "heterocyclic compound binding"                                  
# [13] "cell adhesion molecule binding"                                 
# [14] "organic cyclic compound binding"                                
# [15] "identical protein binding"                                      
# [16] "protein tag"                                                    
# [17] "nucleic acid binding"                                           
# [18] "oxidoreductase activity"                                        
# [19] "primary active transmembrane transporter activity"              
# [20] "protein-containing complex binding"                             
# 

# $KEGG
#  [1] "Amyotrophic lateral sclerosis"                    
#  [2] "Parkinson disease"                                
#  [3] "Prion disease"                                    
#  [4] "Huntington disease"                               
#  [5] "Spliceosome"                                      
#  [6] "Pathways of neurodegeneration - multiple diseases"
#  [7] "Oxidative phosphorylation"                        
#  [8] "Alzheimer disease"                                
#  [9] "Diabetic cardiomyopathy"                          
# [10] "Chemical carcinogenesis - reactive oxygen species"
# [11] "Thermogenesis"                                    
# [12] "Proteasome"                                       
# [13] "mRNA surveillance pathway"                        
# [14] "Nucleocytoplasmic transport"                      
# [15] "Non-alcoholic fatty liver disease"                
# 
# $REAC
#  [1] "Metabolism of RNA"                                                                                                  
#  [2] "Processing of Capped Intron-Containing Pre-mRNA"                                                                    
#  [3] "mRNA Splicing - Major Pathway"                                                                                      
#  [4] "mRNA Splicing"                                                                                                      
#  [5] "Cellular response to chemical stress"                                                                               
#  [6] "Cellular responses to stress"                                                                                       
#  [7] "Cellular responses to stimuli"                                                                                      
#  [8] "Host Interactions of HIV factors"                                                                                   
#  [9] "Cytoprotection by HMOX1"                                                                                            
# [10] "Negative regulation of NOTCH4 signaling"                                                                            
# [11] "Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins."
# [12] "Stabilization of p53"                                                                                               
# [13] "HIV Infection"                                                                                                      
# [14] "Degradation of beta-catenin by the destruction complex"                                                             
# [15] "ER-Phagosome pathway"                                                                                               
# [16] "Infectious disease"                                                                                                 
# [17] "Regulation of HMOX1 expression and activity"                                                                        
# [18] "AUF1 (hnRNP D0) binds and destabilizes mRNA"                                                                        
# [19] "Regulation of RUNX3 expression and activity"                                                                        
# [20] "GLI3 is processed to GLI3R by the proteasome"                                                                       
# 
# $TF
#  [1] "Factor: E2F-3:HES-7; motif: NNNSGCGCSNNNNNCRCGYGNN; match class: 1"
#  [2] "Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN; match class: 1"           
#  [3] "Factor: ZF5; motif: GGSGCGCGS; match class: 1"                     
#  [4] "Factor: E2F-4; motif: SNGGGCGGGAANN"                               
#  [5] "Factor: ER81; motif: RCCGGAARYN; match class: 1"                   
#  [6] "Factor: E2F-1; motif: NGGGCGGGARV"                                 
#  [7] "Factor: E2F; motif: GGCGSG"                                        
#  [8] "Factor: ZF5; motif: GSGCGCGR; match class: 1"                      
#  [9] "Factor: ZF5; motif: GGSGCGCGS"                                     
# [10] "Factor: ER81; motif: RCCGGAARYN"                                   
# [11] "Factor: E2F-4; motif: SNGGGCGGGAANN; match class: 1"               
# [12] "Factor: E2F4; motif: YCCCGCCNCNNSSNNSNN"                           
# [13] "Factor: ZF5; motif: GSGCGCGR"                                      
# [14] "Factor: E2F-2; motif: NWTTTGGCGCCAWWNN"                            
# [15] "Factor: Egr-1; motif: GCGCATGCG; match class: 1"                   
# [16] "Factor: E2F-4; motif: NTTTCSCGCC"                                  
# [17] "Factor: E2F-1; motif: NNNSSCGCSAANN"                               
# [18] "Factor: p300; motif: ACNTCCG; match class: 1"                      
# [19] "Factor: E2F-4:DP-1; motif: TTTSGCGC"                               
# [20] "Factor: EHF; motif: CSCGGAARTN"                                    
# 
# $WP
#  [1] "mRNA processing"                                                 
#  [2] "Electron transport chain: OXPHOS system in mitochondria"         
#  [3] "Cholesterol synthesis disorders"                                 
#  [4] "Oxidative phosphorylation"                                       
#  [5] "VEGFA-VEGFR2 signaling pathway"                                  
#  [6] "Cholesterol biosynthesis pathway"                                
#  [7] "Cholesterol metabolism with Bloch and Kandutsch-Russell pathways"
#  [8] "Proteasome degradation"                                          
#  [9] "Sterol regulatory element-binding proteins (SREBP) signaling"    
# [10] "Nonalcoholic fatty liver disease" 

#save
saveRDS(sicbl,fp(out,"sicbls.rds"))
fwrite(degs,fp(out,"degs_siklf2_vs_ctrl_mpplmpp.csv.gz"))
fwrite(resgos,fp(out,"res_gost_degs_siklf2_vs_ctrl_mpplmpp.csv.gz"))

#degs MPP
#DEGs between sirnas ?
Idents(sicbl)<-"predicted.cell_type"
degs<-FindMarkers(sicbl,group.by = "sirna",ident.1 = "KLF2",ident.2 = "CTRL",subset.ident = "MPP")
degs<-data.table(degs,keep.rownames = "gene")
degs[p_val_adj<0.001] #101
degs_list<-split(degs[p_val_adj<0.001]$gene,degs[p_val_adj<0.001]$avg_log2FC>0)
names(degs_list)<-c("down","up")
lapply(degs_list, head)
# $down
# [1] "MTRNR2L12" "LRRC75A"   "RPL30"     "PLCG2"    
# [5] "HNRNPH1"   "RPS27A"   
# 
# $up
# [1] "CD81"

res_gos<-gost(query = degs_list$down,organism = "hsapiens")
gprofiler2::gostplot(res_gos)
resgos<-data.table(res_gos$result)
lapply(split(resgos,by = "source"),head)

lapply(split(resgos,by = "source"),function(x)head(x[p_value<0.001]$term_name,20))
# $CORUM
# [1] "Ribosome, cytoplasmic"
# 
# $`GO:BP`
#  [1] "cytoplasmic translation"                            
#  [2] "ATP metabolic process"                              
#  [3] "peptide metabolic process"                          
#  [4] "translation"                                        
#  [5] "peptide biosynthetic process"                       
#  [6] "cellular amide metabolic process"                   
#  [7] "cellular macromolecule biosynthetic process"        
#  [8] "aerobic respiration"                                
#  [9] "amide biosynthetic process"                         
# [10] "energy derivation by oxidation of organic compounds"
# [11] "oxidative phosphorylation"                          
# [12] "generation of precursor metabolites and energy"     
# [13] "organonitrogen compound metabolic process"          
# 
# 
# $`GO:MF`
# [1] "RNA binding"     "mRNA binding"   
# [3] "protein binding" "protein tag"    
# 
# 
# $KEGG
# [1] "Spliceosome"                  
# [2] "Shigellosis"                  
# [3] "Ribosome"                     
# [4] "Amyotrophic lateral sclerosis"
# 
# 
# $REAC
#  [1] "Metabolism of RNA"                                                           
#  [2] "SRP-dependent cotranslational protein targeting to membrane"                 
#  [3] "Peptide chain elongation"                                                    
#  [4] "Viral mRNA Translation"                                                      
#  [5] "Eukaryotic Translation Termination"                                          
#  [6] "Selenocysteine synthesis"                                                    
#  [7] "Eukaryotic Translation Elongation"                                           
#  [8] "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)"
#  [9] "Infectious disease"                                                          
# [10] "Regulation of expression of SLITs and ROBOs"                                 
# [11] "Response of EIF2AK4 (GCN2) to amino acid deficiency"                         
# [12] "Formation of a pool of free 40S subunits"                                    
# [13] "mRNA Splicing - Major Pathway"                                               
# [14] "mRNA Splicing"                                                               
# [15] "L13a-mediated translational silencing of Ceruloplasmin expression"           
# [16] "GTP hydrolysis and joining of the 60S ribosomal subunit"                     
# [17] "Influenza Infection"                                                         
# [18] "Nonsense-Mediated Decay (NMD)"                                               
# [19] "Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)"   
# [20] "Selenoamino acid metabolism"                                                 
# 
# $TF
# [1] "Factor: E2F-1; motif: NGGGCGGGARV"
# 
#save
fwrite(degs,fp(out,"degs_siklf2_vs_ctrl_mpp.csv.gz"))
fwrite(resgos,fp(out,"res_gost_degs_siklf2_vs_ctrl_mpp.csv.gz"))

