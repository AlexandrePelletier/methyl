
library(data.table)
library(stringr)
dir.create("outputs/rnaseq_cbps")
test<-fread("~/RUN/Run_595_RNA/Output/RSEM/cbl-ctrl1_Ensembl-104.genes_genename.results")

rsem_out<-"~/RUN/Run_595_RNA/Output/RSEM/"

files_rsem_out<-list.files(path = rsem_out,pattern = "cbp-?[0-9]+_Ensembl-104.genes_genename.results")
length(files_rsem_out)
sapply(paths_rsem_out, function(path)str_remove(str_extract(path,"cbp-[0-9]+"),"-"))
  
cbps_list<-lapply(files_rsem_out, function(file){
  cbp_name<-str_remove(str_extract(file,"cbp-?[0-9]+"),"-")
  path<-file.path(rsem_out,file)
  return(fread(path)[,.(gene_id.1,expected_count)][,(cbp_name):=expected_count][,-"expected_count"])
  })

cbps_list
cbps_mat<-Reduce(merge,cbps_list)

test[,hgnc_symbol:=external_gene_name]
cbps_mat<-merge(cbps_mat,test[,.(gene_id.1,gene_id,hgnc_symbol)])
fwrite(cbps_mat,"outputs/rnaseq_cbps/bulk_rna_seq_gene_count_matrix_cbps170621.csv")
