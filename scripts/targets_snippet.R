f <- function(input, out) {
  dir.create(out)
  
  
  fread(input,sep="\t",
             select = c("cpg_id","P.Value","adj.P.Val","AveExpr","logFC"),
             col.names = c("cpg_id","pval","padj","avg.meth","meth.change"))
}


list(
  tar_target(file1_path, "outputs/01-lga_vs_ctrl_limma_DMCs_analysis/res_limma.tsv.gz", format = "file"),
  tar_target(file1_dt, command = f(input = file1_path, out = "")),
  tar_target(file1_plot, command = plot(file1_dt)),
  tar_target(file1_g, command = g(file1_dt))
)
)

mtd<-tar_read(file1_dt)
g <- function(x) ...
g(mtd)
