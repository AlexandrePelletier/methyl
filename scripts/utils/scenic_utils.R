GetRegulonActivity<-function(analysis.dir ){
  require("SCENIC")
  require("stringr")
  require("withr")
  withr::with_dir(analysis.dir,
      tf_mat <- loadInt(readRDS("int/scenicOptions.Rds"),"aucell_regulonAUC")@assays@data$AUC
  )
  
  regulon_names<-rownames(tf_mat)
  rownames(tf_mat)<-str_replace(str_remove(regulon_names," \\(.*\\)"),"_extended","e")
  return(tf_mat)
}


AddGeneTargetInfo<-function(seurat_object,analysis.dir,assay="SCENIC"){
  require("SCENIC")
  require("Seurat")
  require(withr)
  require("stringr")

  withr::with_dir(analysis.dir,{
  scenic_options<-readRDS("int/scenicOptions.Rds")
  regulon_targets_infos <- loadInt(scenic_options, "regulonTargetsInfo")
  }

  )
    regulon_targets_infos[highConfAnnot==T,regulon:=TF]
  regulon_targets_infos[highConfAnnot==F,regulon:=paste0(TF,"e")]
  targets<-unique(regulon_targets_infos[,targets:=paste(gene,collapse = "/"),by="regulon"][,.(regulon,targets)])
  setwd(here())
  seurat_object[[assay]]<-AddMetaData(object=seurat_object[[assay]],
                                      metadata=data.frame(targets,row.names = "regulon")[rownames(seurat_object),"targets"],
                                      col.name="targets")
  
  return(seurat_object)

}





