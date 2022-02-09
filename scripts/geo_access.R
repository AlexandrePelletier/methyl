#geo data
out<-"outputs/geo_access"
dir.create(out)

#methyldata
mtd<-fread("datasets/cd34/metadata_cl_190421.csv")
methyl_files<-list.files("datasets/cd34/HPA",full.names = T,pattern = ,paste0(c("mspi",mtd$sample),collapse = "|"))


tar("outputs/geo_access/methyl_counts.tar.gz",methyl_files,compression = "gzip")

#scrnaseq data
h5_dirs<-c(cbp2="~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp2b_tri/single_cell_barcode_539_HTO_cbp2b/outs/",
           cbp4="~/RUN/Run_539_10x_standard/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/",
           cbp3="~/RUN/Run_505_10x_standard/Output/cellranger_count/run_505_10xm_standard_CBP3_10x-CBP3/outs/",
           cbp6a="~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-a/outs/",
           cbp6b="~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-b/outs",
           cbp6c="~/RUN/Run_538_10x_standard/Output/cellranger_count/run_538_10x-cbp6-c/outs/",
           cbp7a="~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-a/outs/",
           cbp7b="~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-b/outs/",
           cbp7c="~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp7-c/outs/",
           cbp8="~/RUN/Run_554_10x_standard/Output/cellranger_count/single_cell_barcode_run_554_10xcbp8/outs/")
h5_files<-paste0(h5_dirs,"filtered_feature_bc_matrix.h5")
raw_files<-c(h5_files,cbp3_hto="~/RUN/Run_505_10x_hto/Output/CiteSeq/HTO_CBP3/umi_count/")

tar("outputs/geo_access/scRNAseq_counts_files.tar.gz",raw_files,compression = "gzip")
