if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MAGeCKFlute")
install.packages("ggplot2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

library(MAGeCKFlute)
library(ggplot2)
library(clusterProfiler)
library("gprofiler2")
library("xlsx")
setwd("/Volumes/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data/1. batch/mageck-tables")

if(length(list.files(".",pattern = "contrasts",r=F))>=2) {
  stop("There are two files containing the word contrasts, please rename one" )
} else if(length(list.files(".",pattern = "contrasts",r=F))==0) {
  stop("The contrast file is missing in your folder directory")
}

contrasts <- read.delim(list.files(".",pattern = "contrasts",r=F)[1], header=F,sep = "\t")
contrasts 
gdata = ReadBeta("plasmid.gene_summary.txt")

dir.create("Mageck-graphs")
dir.create("Mageck-graphs/Scatter")
dir.create("Mageck-graphs/Ranks")
dir.create("Mageck-graphs/9square")
dir.create("Mageck-graphs/Consistency")
dir.create("Mageck-graphs/Pathway_analysis")
dir.create("Mageck-graphs/Density")

for(i in 1:nrow(contrasts)){
  
  jpeg(paste0("Mageck-graphs/Density/DensityView_",contrasts$V1[i],"_",contrasts$V2[i],".jpeg"))
  p1 = DensityView(gdata, samples=c(contrasts$V1[i],contrasts$V2[i] ))
  print(p1)
  dev.off()
  
  jpeg(paste0("Mageck-graphs/Scatter/ScatterView_",contrasts$V1[i],"_",contrasts$V2[i],".jpeg"))
  p2 = ScatterView(gdata, auto_cut_diag = TRUE, display_cut = TRUE, x=contrasts$V1[i],label = "Gene",
                   y=contrasts$V2[i],
                   groups = c("top", "bottom"))
  print(p2)
  dev.off()
  
  
  ########Pathway analysis##########
  datasources <- c("KEGG", "REAC")
  min_DEG_pathway = 2
  for(g in c("top","bottom")) {
  gostres <- gost(query=as.character(p2$data[p2$data$group == g,]$Gene),
                  organism="hsapiens",
                  significant=TRUE,
                  correction_method="fdr",
                  sources=datasources,
                  evcodes=TRUE,
                  user_threshold=0.05)
  
  # Make data frame of gost result
  pathway_gostres <- gostres$result
  # Select only significantly enriched pathways (according to adjusted p-value)
  pathway_gostres <- as.data.frame(pathway_gostres[which(pathway_gostres$significant==TRUE),])
  # Select only pathways with a min. number of DEG
  pathway_gostres <- pathway_gostres[which(pathway_gostres$intersection_size>=min_DEG_pathway),]
  pathway_gostres_table = pathway_gostres
  pathway_gostres_table$parents <- NULL
  
 # write.table(pathway_gostres_table, file=paste0("Mageck-graphs/Pathway_analysis/Pathway_analysis_",contrasts$V1[i],"_",
  #                                         contrasts$V2[i],"_",g,".txt")
   #           ,sep="\t", quote = F, col.names = T, row.names = F)
  write.xlsx(pathway_gostres_table, file=paste0("Mageck-graphs/Pathway_analysis/Pathway_analysis_",contrasts$V1[i],"_",
                                                contrasts$V2[i],"_",g,".xlsx"),
             sheetName = "pathway", append = FALSE)
  }
  
  jpeg(paste0("Mageck-graphs/Consistency/ConsistencyView_",contrasts$V1[i],"_",contrasts$V2[i],".jpeg"))
  p3 = ConsistencyView(gdata, contrasts$V1[i], contrasts$V2[i])
  print(p3)
  dev.off()
  
  ctrlname <- contrasts$V1[i]
  treatname <- contrasts$V2[i]
  
  #gdata_cc[apply(gdata_cc[,-1], 1, function(x) !all(x>10)),]
  
  #gdata_cc = NormalizeBeta(gdata, id = "Gene", samples=c(ctrlname, treatname), 
  #                        method="cell_cycle")
  gdata_cc <- gdata
  head(gdata_cc)
  #
  # gdata_cc[apply(gdata_cc[,-1], 1, function(x) !all(x>10)),]
  gdata_cc$Control = rowMeans(gdata_cc[,ctrlname,drop=FALSE])
  gdata_cc$Treatment = rowMeans(gdata_cc[,treatname,drop=FALSE])
  
  gdata_cc$Diff = gdata_cc$Treatment - gdata_cc$Control
  gdata_cc$Rank = rank(gdata_cc$Diff)
  
  jpeg(paste0("Mageck-graphs/Ranks/Rank_names_",contrasts$V1[i],"_",contrasts$V2[i],".jpeg"))
  p1 = ScatterView(gdata_cc, x="Rank", y="Diff",label="Gene",
                   top=5,auto_cut_y=TRUE, groups=c("top","bottom"))
  print(p1)
  dev.off()
  
  jpeg(paste0("Mageck-graphs/9square/ScatterView_grid_names_",contrasts$V1[i],"_",contrasts$V2[i],".jpeg"))
  p1 = ScatterView(gdata_cc, x=ctrlname, y=treatname, label = "Gene",
                   model = "ninesquare", top=5, display_cut=TRUE,
                   x_cut = c(-0.5,0.5), y_cut = c(-0.5,0.5)
  )
  print(p1)
  dev.off()
  
  #Square9 = p1$data
  #idx=Square9$group=="topcenter"
  #geneList = Square9$Diff[idx]
  #names(geneList) = Square9$Gene[idx]
  #universe = Square9$Gene
  
  # Enrichment analysis
  #kegg1 = EnrichAnalyzer(geneList = geneList, universe = universe)
  
  #jpeg(paste0("EnrichedView_",contrasts$V1[i],"_",contrasts$V2[i]))
  #EnrichedView(kegg1, top = 6, bottom = 0)
  #dev.off()
}


