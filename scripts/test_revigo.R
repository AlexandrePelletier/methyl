renv::install("bioc::rrvgo")

library(rrvgo)
go_analysis<-fread("outputs/08-HTO_signature/res_hto_signature_go_bp_up.csv")
go_analysis_sig<-go_analysis[p.adjust<0.001]
simMatrix <- calculateSimMatrix(go_analysis_sig$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")


scores <- setNames(-log10(go_analysis_sig$qvalue), go_analysis_sig$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)

scatterPlot(simMatrix, reducedTerms,size = "score")


treemapPlot(reducedTerms)
