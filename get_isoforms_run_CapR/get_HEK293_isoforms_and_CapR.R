## data from GEO: GSE122425 Sun et al 2018 Epigenomics

library(dplyr)
library(nearBynding)

data1<-read.delim("SRR8181090.gtf", sep = "\t", stringsAsFactors = F, header = F) %>% filter(V3 == "transcript")
data2<-read.delim("SRR8181091.gtf", sep = "\t", stringsAsFactors = F, header = F) %>% filter(V3 == "transcript")
data3<-read.delim("SRR8181092.gtf", sep = "\t", stringsAsFactors = F, header = F) %>% filter(V3 == "transcript")

data1_cut<-lapply(data1$V9, function(x){strsplit(x, ";")})
data1_cut_df <- data.frame(matrix(unlist(data1_cut), nrow=length(data1_cut), byrow=T), stringsAsFactors = F)
colnames(data1_cut_df)<-c("gene_id", "transcript_id", "ref_gene_name", "cov", "FPKM", "TPM")
data1_cut_df<-sapply(data1_cut_df, function(x){gsub(".* ", "", x)}) %>% as.data.frame()
data1_cut_df$TPM<-as.numeric(as.character(data1_cut_df$TPM))
# only transcripts with expression
data1_df<-data1_cut_df %>% filter(TPM > 0)
# take the highest expressed isoform of every gene if multiple isoforms are expressed
data1_df<-aggregate(TPM ~ gene_id, data = data1_df, max)
data1_df$transcript_id<-""
for(t in 1:nrow(data1_df)){
  data1_df[t,"transcript_id"]<-data1_cut_df[which(as.character(data1_cut_df$gene_id) == 
                                                      as.character(data1_df[t, "gene_id"]) &
                                                  data1_cut_df$TPM == data1_df[t, "TPM"]),2] %>% 
      as.character()
}

data2_cut<-lapply(data2$V9, function(x){strsplit(x, ";")})
data2_cut_df <- data.frame(matrix(unlist(data2_cut), nrow=length(data2_cut), byrow=T), stringsAsFactors = F)
colnames(data2_cut_df)<-c("gene_id", "transcript_id", "ref_gene_name", "cov", "FPKM", "TPM")
data2_cut_df<-sapply(data2_cut_df, function(x){gsub(".* ", "", x)}) %>% as.data.frame()
data2_cut_df$TPM<-as.numeric(as.character(data2_cut_df$TPM))
data2_df<-data2_cut_df %>% filter(TPM > 0)
data2_df<-aggregate(TPM ~ gene_id, data = data2_df, max)
data2_df$transcript_id<-""
for(t in 2:nrow(data2_df)){
  data2_df[t,"transcript_id"]<-data2_cut_df[which(as.character(data2_cut_df$gene_id) == 
                                                      as.character(data2_df[t, "gene_id"]) &
                                                    data2_cut_df$TPM == data2_df[t, "TPM"]),2] %>% 
      as.character()
}

data3_cut<-lapply(data3$V9, function(x){strsplit(x, ";")})
data3_cut_df <- data.frame(matrix(unlist(data3_cut), nrow=length(data3_cut), byrow=T), stringsAsFactors = F)
colnames(data3_cut_df)<-c("gene_id", "transcript_id", "ref_gene_name", "cov", "FPKM", "TPM")
data3_cut_df<-sapply(data3_cut_df, function(x){gsub(".* ", "", x)}) %>% as.data.frame()
data3_cut_df$TPM<-as.numeric(as.character(data3_cut_df$TPM))
data3_df<-data3_cut_df %>% filter(TPM > 0)
data3_df<-aggregate(TPM ~ gene_id, data = data3_df, max)
data3_df$transcript_id<-""
for(t in 3:nrow(data3_df)){
  data3_df[t,"transcript_id"]<-data3_cut_df[which(as.character(data3_cut_df$gene_id) == 
                                                      as.character(data3_df[t, "gene_id"]) &
                                                    data3_cut_df$TPM == data3_df[t, "TPM"]),2] %>% 
      as.character()
}

## get list of transcripts present in all three samples
HEK293_transcripts<-data3_df$transcript_id[data3_df$transcript_id %in% 
                                               data1_df$transcript_id[data1_df$transcript_id %in% 
                                                                          data2_df$transcript_id]]
writeLines(HEK293_transcripts, "isoforms_HEK293.txt")

# map to 3'UTR
GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                          out_chain_name = "HEK293_3UTR.chain",
                          RNA_fragment = "three_prime_utr",
                          transcript_list = HEK293_transcripts,
                          alignment = "hg38")
getChainChrSize(chain = "HEK293_3UTR.chain",
                out_chr = "HEK293_3UTR.size")
ExtractTranscriptomeSequence(transcript_list = HEK293_transcripts,
                             ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                             genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                             RNA_fragment = "three_prime_utr",
                             exome_prefix = "HEK293_3UTR")
runCapR("HEK293_3UTR.fa")
processCapRout(CapR_outfile = "HEK293_3UTR.out",
                 output_prefix = "HEK293_3UTR",
                 chrom_size = "HEK293_3UTR.size",
                 genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                 RNA_fragment = "three_prime_utr",
                 chain = "HEK293_3UTR.chain")
