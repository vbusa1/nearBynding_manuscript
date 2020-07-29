library(biomaRt)
library(dplyr)

file_info<-read.csv("file_information.csv", stringsAsFactors = F)

## ANNOTATE THE DATA
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
proteins <- file_info$protein %>% unique()
protein_domains <- getBM(attributes=c("hgnc_symbol","interpro","interpro_description"),
                         filters="hgnc_symbol",
                         values=proteins,
                         mart=ensembl)
RBP_annotation<-as.data.frame(matrix(NA, nrow = length(unique(file_info$protein)), ncol = 12))
colnames(RBP_annotation)<-c("helicase", "zinc_finger", "RRM", "KH", "dsRBD", "AB", "armadillo", "P_loop", "G_patch", "SAP", "WD40", "winged_helix")
row.names(RBP_annotation)<-RBP_annotation$Protein<-unique(file_info$protein)
for (protein in unique(file_info$protein)){
  RBP_annotation[protein, "helicase"]<- sum(grepl("helicase", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                              "interpro_description"],
                                                  ignore.case = T)) >0
  RBP_annotation[protein, "KH"]<-sum(grepl("K Homology", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                         "interpro_description"],
                                           ignore.case = T)) >0
  RBP_annotation[protein, "zinc_finger"]<-sum(grepl("zinc finger", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                                   "interpro_description"],
                                                    ignore.case = T)) >0
  RBP_annotation[protein, "RRM"]<-sum(grepl("RNA recognition motif", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                                     "interpro_description"],
                                            ignore.case = T)) >0
  RBP_annotation[protein, "dsRBD"]<-sum(grepl("Double-stranded RNA-binding", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                                             "interpro_description"])) >0
  RBP_annotation[protein, "AB"]<-sum(grepl("alpha-beta plait", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                               "interpro_description"],
                                           ignore.case = T)) >0
  RBP_annotation[protein, "armadillo"]<-sum(grepl("Armadillo", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                               "interpro_description"],
                                                  ignore.case = T)) >0
  RBP_annotation[protein, "P_loop"]<-sum(grepl("P-loop", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                         "interpro_description"],
                                               ignore.case = T)) >0
  RBP_annotation[protein, "G_patch"]<-sum(grepl("G-patch", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                           "interpro_description"],
                                                ignore.case = T)) >0
  RBP_annotation[protein, "SAP"]<-sum(grepl("SAP", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                   "interpro_description"],
                                            ignore.case = T)) >0
  RBP_annotation[protein, "WD40"]<-sum(grepl("WD", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                   "interpro_description"],
                                             ignore.case = T)) >0
  RBP_annotation[protein, "winged_helix"]<-sum(grepl("Winged helix", protein_domains[which(protein_domains$hgnc_symbol == protein),
                                                                                     "interpro_description"],
                                                     ignore.case = T)) >0
}

write.csv(RBP_annotation, "RBP_annotation.csv", row.names = F)
