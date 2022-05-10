#put the binary directory in the path
old_path <- Sys.getenv("PATH")

Sys.setenv(PATH = paste(old_path, "/disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12","/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin", sep = ":"))
Sys.getenv("PATH")

dir.create("outputs/20-RNA_velocity/velocyto_counts/cbp4/",recursive = T)

system("velocyto run -b  ~/RUN/run_539_single_cell/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o outputs/20-RNA_velocity/velocyto_counts/cbp4/ -m ref/hg38_rmsk.gtf ~/RUN/run_539_single_cell/Output/cellranger_count_cbp4_tri/single_cell_barcode_539_HTO_cbp4b/outs/possorted_genome_bam.bam ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

