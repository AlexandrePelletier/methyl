#put the binary directory in the path
old_path <- Sys.getenv("PATH")

Sys.setenv(PATH = paste(old_path, "/disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12","/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin", sep = ":"))
Sys.getenv("PATH")
system("velocyto")

dir.create("outputs/20-RNA_velocity/velocyto_counts/cbp3/",recursive = T)

system("velocyto run -b  ~/RUN/run_505_10xm_standard/Output/cellranger_count/run_505_10xm_standard_CBP3_10x-CBP3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o outputs/20-RNA_velocity/velocyto_counts/cbp3/ -m ref/hg38_rmsk.gtf ~/RUN/run_505_10xm_standard/Output/cellranger_count/run_505_10xm_standard_CBP3_10x-CBP3/outs/possorted_genome_bam.bam ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf")


