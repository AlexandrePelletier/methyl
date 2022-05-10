#put the binary directory in the path
old_path <- Sys.getenv("PATH")

Sys.setenv(PATH = paste(old_path, "/disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12","/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin", sep = ":"))
Sys.getenv("PATH")
system("velocyto")
dir.create("outputs/20-RNA_velocity/velocyto_counts/cbp2/",recursive = T)

#with in bash :

system("velocyto run -b  ~/RUN/run_539_single_cell/Output/cellranger_count_cbp2b_tri/single_cell_barcode_539_HTO_cbp2b/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o outputs/20-RNA_velocity/velocyto_counts/cbp2/ -m ref/hg38_rmsk.gtf ~/RUN/run_539_single_cell/Output/cellranger_count_cbp2b_tri/single_cell_barcode_539_HTO_cbp2b/outs/possorted_genome_bam.bam ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

