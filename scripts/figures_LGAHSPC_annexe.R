#2B-HTO_signature
#scDEGs
res_dup<-fread("../singlecell/outputs/03-HTOs_stimulation/duplicates/scdegs_hto_vs_not_wilcoxon.csv")
#volcano
genes_of_interest<-c("SOCS3","HES1","JUN","FOS","JUNB","ZFP36","EGR1",
                      "DUSP2","DUSP1","FOSB","SOCS1","KLF2","KLF4",
                       "PLK2","PLK3","ID1","MYC","","ID2","IDS","RGCC")

ggplot(res_dup,aes(x=avg_log2FC,y=-log10(p_val_adj),col=p_val_adj<0.01&abs(avg_log2FC)>0.4))+
  geom_point()+ 
  geom_label_repel(aes(label = ifelse(p_val_adj<0.01&
                                        abs(avg_log2FC)>0.4&
                                        gene%in%genes_of_interest,gene,"")),
                   max.overlaps = 3000,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_manual(values = c("grey","red")) +
  theme_minimal() +
  theme(legend.position = "bottom")+facet_wrap("sample")
ggsave(fp(out,"2B-volcano_scDEGs_antibody_activation_signature_in_duplicates.pdf"))

table(res_dup[p_val_adj<0.01&abs(avg_log2FC)>0.4]$sample)
# ctrlM518 ctrlM537 ctrlM555 
#      572      675      331 

res_dup[,degs_in_3:=sum(p_val_adj<0.01&abs(avg_log2FC)>0.4)==3,by="gene"]
length(unique(res_dup[degs_in_3==T]$gene)) #210
cat(paste(unique(res_dup[degs_in_3==T]$gene),collapse = "\n"))

#supp : signature enrichment by lineage 
source("scripts/utils/new_utils.R")
res_sign_kegg<-fread("outputs/08-HTO_signature/by_lineage/res_hto_signature_kegg_by_lineage.csv")

res_sign_kegg[,n.enriched:=Count]
res_sign_kegg[,n_gene_set:=as.numeric(str_extract(BgRatio,"^[0-9]+"))]
res_sign_kegg[,pct.enriched:=n.enriched/n_gene_set]
ggplot(res_sign_kegg[lineage%in%c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")&p.adjust<0.05])+
  geom_point(aes(x=Description,y=pct.enriched,size=n.enriched,col=-log10(p.adjust)))+
  facet_grid(~lineage,scales = "free_x",space = "free")+
  theme(axis.text.x = element_text(angle = 60,vjust =1,hjust = 1))
ggsave("outputs/figures_epi_response/figure2/2Bsupp-kegg_dotplot_antibody_activation_signature_by_lineage.pdf")



res_sign_go_list<-readRDS("outputs/08-HTO_signature/by_lineage/res_hto_signature_go_bp_by_lineage.rds")
res_sign_go<-Reduce(function(x,y)rbind(x,y,fill=TRUE),lapply(names(res_sign_go_list), function(lin)data.table(as.data.frame(res_sign_go_list[[lin]]))[,lineage:=lin]))

res_sign_go[,n.enriched:=Count]
res_sign_go[,n_gene_set:=as.numeric(str_extract(BgRatio,"^[0-9]+"))]
res_sign_go[,pct.enriched:=n.enriched/n_gene_set]
res_sign_go[,top10:=p.adjust<=max(sort(p.adjust)[1:10],na.rm = T),by="lineage"]
ggplot(res_sign_go[top10==T&lineage%in%c("HSC","MPP/LMPP","Erythro-Mas","Lymphoid","Myeloid")&p.adjust<0.05])+
  geom_point(aes(x=Description,y=pct.enriched,size=n.enriched,col=-log10(p.adjust)))+
  facet_grid(~lineage,scales = "free_x",space = "free")+
  theme(axis.text.x = element_text(angle = 66,vjust =1,hjust = 1))
ggsave("outputs/figures_epi_response/figure2/2Bsupp-go_dotplot_antibody_activation_signature_by_lineage.pdf")


#percentage of cells below each pseudotime
mtd
n_bins<-30
pseudotime_thr<-1:n_bins*(max(mtdf$pseudotime,na.rm = T)/n_bins)
pseudo_bins<-data.table(bin=0:(n_bins-1),
                        pseudotime_thr=pseudotime_thr)
mtdf[,bin:=sum(pseudotime>=pseudotime_thr),by="bc"]
summary(mtdf$bin)

pseudo_bins_mtd<-merge(pseudo_bins,mtdf)

pseudo_bins_mtd[,n.cells:=.N,by="sample_hto"]
pseudo_bins_mtd[,n.cells.bin:=.N,by=.(sample_hto,bin)]
pseudo_bins_mtd[,pct.cells.bin:=.N/n.cells,by=.(sample_hto,bin)]
pseudo_bins_<-unique(pseudo_bins_mtd,by=c("sample_hto","bin"))
pseudo_bins_[,n.cells.below.thr:=NA]
pseudo_bins_[,n.cells.below.thr:=as.numeric(n.cells.below.thr)]

for(binn in 0:n_bins){
  pseudo_bins_[bin<=binn,n.cells.below.thr:=ifelse(is.na(n.cells.below.thr),sum(n.cells.bin),n.cells.below.thr),by="sample_hto"]

  }

pseudo_bins_[,pct.cells.below.thr:=n.cells.below.thr/n.cells]
ggplot(pseudo_bins_,aes(x=factor(bin),y=pct.cells.below.thr,fill=group))+
  geom_boxplot()+
facet_wrap("hto")
ggsave(fp(out,"boxplot_cumulated_pseudotime_by_group.pdf"))

ggplot(pseudo_bins_[hto==T],aes(x=factor(bin),y=pct.cells.below.thr,fill=group))+
  geom_boxplot()+
facet_wrap("hto")

plot(mtd$pseudotime)


# pseudotime influencé par group:hto ?
glm.bins<-stats::glm(n.cells.bin~n.cells+pseudotime_thr*group*hto+batch,family=poisson(),data = pseudo_bins_)
summary(glm.bins)

lme4::lmer(formula = time ~ 1 + group + XXX + (1|ID), ...)

ggplot(mtd[(hto)])+
  # geom_density(aes(x=pseudotime,group=sample)) +
  geom_density(aes(x=pseudotime,col = group))


# Call:
# stats::glm(formula = n.cells.bin ~ n.cells + pseudotime_thr * 
#     group * hto + batch, family = poisson(), data = pseudo_bins_)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -19.0806   -5.8031   -0.6116    3.5736   22.0908  
# 
# Coefficients: (1 not defined because of singularities)
#                                   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                      1.743e+00  8.509e-02  20.480  < 2e-16 ***
# n.cells                          5.848e-04  1.315e-05  44.468  < 2e-16 ***
# pseudotime_thr                  -3.750e-03  9.145e-04  -4.101 4.12e-05 ***
# grouplga                         2.534e-02  2.621e-02   0.967    0.334    
# htoTRUE                          1.282e+00  8.635e-02  14.842  < 2e-16 ***
# batchcbp0_lga                    1.375e+00  6.412e-02  21.445  < 2e-16 ***
# batchcbp2                        1.627e-01  2.682e-02   6.068 1.29e-09 ***
# batchcbp4                       -6.519e-03  2.744e-02  -0.238    0.812    
# batchcbp6a                       1.330e+00  5.935e-02  22.408  < 2e-16 ***
# batchcbp6b                       1.413e+00  7.076e-02  19.975  < 2e-16 ***
# batchcbp6c                       1.475e+00  6.540e-02  22.562  < 2e-16 ***
# batchcbp7a                       1.399e+00  6.287e-02  22.253  < 2e-16 ***
# batchcbp7b                       1.276e+00  5.611e-02  22.743  < 2e-16 ***
# batchcbp7c                       1.405e+00  5.504e-02  25.526  < 2e-16 ***
# batchcbp8                               NA         NA      NA       NA    
# pseudotime_thr:grouplga         -1.091e-02  1.332e-03  -8.187 2.67e-16 ***
# pseudotime_thr:htoTRUE          -1.280e-02  2.150e-03  -5.953 2.63e-09 ***
# grouplga:htoTRUE                -2.021e-01  4.510e-02  -4.482 7.38e-06 ***
# pseudotime_thr:grouplga:htoTRUE  3.034e-02  2.799e-03  10.841  < 2e-16 ***


#=>  counts par pseudotime est influencé par group_hto.

#test counts par pseudotime influencé par group in stimulated cond ?
glm.bins_hto<-stats::glm(n.cells.bin~n.cells+pseudotime_thr*group+batch,family=poisson(),data = pseudo_bins_[hto==T])
summary(glm.bins_hto)

# Call:
# stats::glm(formula = n.cells.bin ~ n.cells + pseudotime_thr * 
#     group + batch, family = poisson(), data = pseudo_bins_[hto == 
#     T])
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -12.8282   -4.4758   -0.5995    2.1113   14.6747  
# 
# Coefficients:
#                           Estimate Std. Error z value Pr(>|z|)    
# (Intercept)              2.9826825  0.0402223  74.155  < 2e-16 ***
# n.cells                  0.0007477  0.0000200  37.377  < 2e-16 ***
# pseudotime_thr          -0.0167324  0.0019449  -8.603  < 2e-16 ***
# grouplga                -0.2820240  0.0406433  -6.939 3.95e-12 ***
# batchcbp4               -0.0278552  0.0287041  -0.970    0.332    
# batchcbp8               -0.0393944  0.0292333  -1.348    0.178    
# pseudotime_thr:grouplga  0.0191529  0.0024593   7.788 6.82e-15 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
#     Null deviance: 10728  on 317  degrees of freedom
# Residual deviance:  7650  on 311  degrees of freedom
# AIC: 9177.9
# 
# Number of Fisher Scoring iterations: 5

#=>  counts par pseudotime est influencé par lga in stim condition .


#calc AUC by samples [TO UPDATE WITH COMPLETE DATA]
library(zoo)

x <- 1:10
y <- 3*x+25
id <- order(x)

AUC <- sum(diff(x[id])*rollmean(y[id],2))

#for 1
x <- pseudo_bins_[sample_hto=="ctrlF547FALSE"]$pseudotime_thr
y <- pseudo_bins_[sample_hto=="ctrlF547FALSE"]$pct.cells.below.thr
id <- order(x)

AUC <- sum(diff(x[id])*rollmean(y[id],2))

#for all
#need first calculate pct.cumul for each pseudotime [TO DO]
df<-data.table(pseudotime=pseudotime_thr)
pseudo_bins_[,AUC:=sum(diff(pseudotime_thr[order(pseudotime_thr)])*rollmean(pct.cells.below.thr[order(pseudotime_thr)],2)),by="sample_hto"]

pseudo_bins_s<-unique(pseudo_bins_,by="sample_hto")
ggplot(pseudo_bins_s)+geom_boxplot(aes(x=group,y=AUC,fill=sex))+facet_wrap("hto")


pseudo_bins_s[,AUC.mean:=mean(AUC),by="group_hto"]
pseudo_bins_s[,AUC.sem:= sd(AUC)/sqrt(.N),by="group_hto"]


ggplot(unique(pseudo_bins_s,by="group_hto"))+geom_col(aes(x=group,y=AUC.mean))+facet_wrap("hto")

plot(density(pseudo_bins_s$AUC))
pseudo_bins_s[AUC<=max(boxplot.stats(pseudo_bins_s$AUC)$out)]
t.test(pseudo_bins_s[group=='lga'&hto==T]$AUC,pseudo_bins_s[group=='ctrl'&hto==T]$AUC) #not sig



#visuellement le pique de l'influence semble être à bin = 7 wich corresponf to HSC cells

#so, is there less cells with pseudotim <7 in stimulated LGA compared to stimulated control ?

glm.bins_7<-stats::glm(n.cells.below.thr~n.cells+group+batch,family=poisson(),data = pseudo_bins_[hto==T&bin==7])
summary(glm.bins_7)
# Call:
# stats::glm(formula = n.cells.below.thr ~ n.cells + group + batch, 
#     family = poisson(), data = pseudo_bins_[hto == T & bin == 
#         7])
# 
# Deviance Residuals: 
#    Min      1Q  Median      3Q     Max  
# -9.700  -4.184  -0.584   3.858   9.342  
# 
# Coefficients:
#               Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  4.591e+00  7.603e-02  60.384  < 2e-16 ***
# n.cells      5.995e-04  5.147e-05  11.646  < 2e-16 ***
# grouplga    -5.096e-01  5.567e-02  -9.155  < 2e-16 ***
# batchcbp4    1.698e-01  6.526e-02   2.602  0.00927 ** 
# batchcbp8    3.564e-01  6.674e-02   5.340  9.3e-08 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for poisson family taken to be 1)
# 
#     Null deviance: 540.77  on 11  degrees of freedom
# Residual deviance: 362.26  on  7  degrees of freedom
# AIC: 453.61
# 
# Number of Fisher Scoring iterations: 5


sapply(unique(pseudo_bins_$bin), function(b){
  tryCatch(
    expr = {
  glm.bins_<-stats::glm(n.cells.below.thr~n.cells+group+batch,family=poisson(),data = pseudo_bins_[hto==T&bin==b])
  res_glm.bins_<-summary(glm.bins_)
  res_glm.bins_$coefficients["grouplga",c("Estimate","Pr(>|z|)")]
  },
  error=function(c)NA
  
  )
  }
  ) #signif from bin 2 - 18, ++bin7

#wilcoxon test stat : 
wilcox.test(pseudo_bins_[hto==T&bin==7&group=="lga"]$pct.cells.below.thr,
            pseudo_bins_[hto==T&bin==7&group=="ctrl"]$pct.cells.below.thr)
#p-value = 0.04113

wilcox.test(pseudo_bins_[hto==F&bin==7&group=="lga"]$pct.cells.below.thr,
            pseudo_bins_[hto==F&bin==7&group=="ctrl"]$pct.cells.below.thr)
#p-value = 0.3

#with lines :
pseudo_bins_[,med.pct.cells:=median(pct.cells.below.thr),by=c("bin","group_hto")]
pseudo_bins_[,mad.pct.cells:=mad(pct.cells.below.thr),by=c("bin","group_hto")]

summary(pseudo_bins_$pct.cells.below.thr)

pseudo_bins_g<-unique(pseudo_bins_,by=c("group_hto","bin"))
ggplot(pseudo_bins_g,aes(x=pseudotime_thr,y=med.pct.cells,ymin=-mad.pct.cells, ymax=+mad.pct.cells,fill=group,linetype=group))+
  geom_line()+
  geom_ribbon(alpha=0.5)+
facet_wrap("hto")
