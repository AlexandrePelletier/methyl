#put the binary directory in the path
source("../methyl/scripts/utils/new_utils.R")
out<-"outputs/20-RNA_velocity/velocyto_counts/cbp8"
dir.create(out,recursive = T)
cellranger_outs<-"~/RUN/run_554_RNA/Output/cellranger_count/single_cell_barcode_run_554_10xcbp8/outs"


old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/disks/DATATMP/PhD_AlexandrePelletier/tools/samtools-1.12","/disks/DATATMP/PhD_AlexandrePelletier/singlecell/python/r-reticulate/bin", sep = ":"))
Sys.getenv("PATH")
system("velocyto")

barcodes_path<-fp(cellranger_outs,"filtered_feature_bc_matrix/barcodes.tsv.gz")
bam_path<-fp(cellranger_outs,"possorted_genome_bam.bam")

cmd_string<-paste("velocyto run -b",barcodes_path,"-o",out,"-m ref/hg38_rmsk.gtf",bam_path,"ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
system(cmd_string)
