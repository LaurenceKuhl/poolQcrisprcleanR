Sys.setenv(LANG = "en")
install.packages("devtools",repos = "http://cran.us.r-project.org")
library(devtools)
install_github("francescojm/CRISPRcleanR")
library(CRISPRcleanR)
library("dplyr")
library("magrittr")
rm(list = ls())
setwd("/Volumes/groups/Neuro-Oncology/Clinical Neuro-Oncology/data backup/GBM_Screens/Broad screening/Knockout_synthetic_lethal_raw_data/1batch")
data("Brunello_Library")
#create the construct barcode present in the poolq output in 
#the brunello table
Brunello <- Brunello_Library
Brunello$Construct.Barcode <- 
  substr(Brunello$Target.Context.Sequence,5,24)
Brunello <- Brunello %>%
  select(Construct.Barcode,CODE,GENES)

#nb reads of the sequenced plasmids
if(length(list.files(".",pattern = "CP0041_Brunello",r=F))>=2) {
  stop("There are two files starting with CP0041_Brunello, please rename one" )
}

control_counts <- read.delim(list.files(".",pattern = "CP0041_Brunello")[1], 
                             header=T,sep = ",")
control_counts <- control_counts %>%
  select(Construct.Barcode,sum.pXPR_003)

#read conditions per well#######################################################
conditionslist <- list.files(path=".",pattern = "_conditions", recursive="true")
conditions <- read.delim(conditionslist[1], header=F,sep = ",")
for(i in 2:length(conditions)) {
  table <- read.delim(conditionslist[i], header=F,sep = ",")
  conditions <- rbind(conditions, table)
  conditions <- conditions[!(is.na(conditions$V2) | conditions$V2==""), ]
}
conditions <- unique(conditions$V2)
#################################################################################

#number of plates
plates <- sub(" ","_", dir(".",pattern="Plate"))

#read count files per plate output by poolq########################
files <- list.files(path=".",pattern = "^counts*", recursive="true")

for(i in 1:length(plates)) {
  table <- read.delim(files[i], header=T,sep="\t")
  assign(plates[i],table)
}
###################################################################

#####################Big table summing up all counts for all replicates and conditions#######

counts <- purrr::reduce(lapply(ls(pattern = "Plate*"), get), dplyr::left_join,
                      by = c('Construct.Barcode',"Construct.IDs"))

counts <- counts  %>% left_join(control_counts) %>% 
                  left_join(Brunello)

conditions <- c("GS.9_DMSO_REP_A"  , "GS.9_DMSO_REP_B"  , "GS.9_Abema_REP_A", 
                "GS.9_Abema_REP_B",  "GS.9_Rego_REP_A",   "GS.9_Rego_REP_B",  
                "GS.9_Evero_REP_A",  "GS.9_Evero_REP_B" , "LN229_DMSO_REP_A", 
                "LN229_DMSO_REP_B" , "LN229_Abema_REP_A" ,"LN229_Abema_REP_B",
                "LN229_Evero_REP_A" ,"LN229_Evero_REP_B" ,"LN229_Rego_REP_A")
conditions <- unique(substr(conditions,1,nchar(conditions)-2))

for(i in conditions) {
  counts[paste(i)] <- counts %>% select(starts_with(i)) %>% rowSums()
}

counts <- counts  %>%
  select(CODE,GENES,sum.pXPR_003,conditions)

#############################################################################
merge_table <- counts  %>%
  select(CODE,GENES,sum.pXPR_003)
colnames(merge_table) <- c("sgRNA","gene","CP0041")

############################CRISPRcleanR#####################################
dir.create("crisprcleanR-graphs")
for(i in conditions){
  input <- counts  %>%
    select(CODE,GENES,sum.pXPR_003,i)
  colnames(input) <- c("sgRNA","gene","CP0041",i)
  normANDfcs <- ccr.NormfoldChanges(Dframe=input,
                                    saveToFig = TRUE,
                                    libraryAnnotation=Brunello_Library,
                                    EXPname=i)
  
  write.csv(normANDfcs$norm_counts, file=paste0(i,"_norm_counts.csv"),row.names=FALSE)
  write.csv(normANDfcs$logFCs, file=paste0(i,"_logFCs.csv"),row.names=FALSE)
  
  gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brunello_Library)
  correctedFCs <- ccr.GWclean(gwSortedFCs,display=FALSE,label=i)
  correctedCounts <- ccr.correctCounts(i,
                                       normANDfcs$norm_counts,
                                       correctedFCs,
                                       Brunello_Library,
                                       minTargetedGenes=3,
                                       OutDir='./')
  
  merge_table <- merge_table %>% right_join(correctedCounts)
}

write.table(merge_table, file="all_norm_counts.tsv",row.names=FALSE, sep = "\t")


