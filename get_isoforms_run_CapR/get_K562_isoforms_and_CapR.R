# determine the most highly expressed transcripts for every gene of K562 cells
# data from ENCODE

library(nearBynding)
library(dplyr)

data1<-read.table("ENCFF509AOR.tsv", header = T, stringsAsFactors = F)
data2<-read.table("ENCFF745GPL.tsv", header = T, stringsAsFactors = F)

# only transcripts with expression
data1_with_reads<-data1[which(data1[,6] > 0),]
genes1<-data1_with_reads$gene_id %>% unique()
# take the highest expressed isoform of every gene if multiple isoforms are expressed
isoform_max1<-data.frame(gene_id=genes1,
                    max_transcript=unlist(lapply(genes1, function(x){
                      transcript<-data1_with_reads[which(data1_with_reads$TPM==
                                                             max(filter(data1_with_reads, gene_id==x)$TPM) &
                                                           data1_with_reads$gene_id==x),
                                                    "transcript_id"] %>% unlist()
                      if(length(transcript)>0){ # if transcripts have same TPM
                        transcript<-transcript[1]
                      }
                      #have to remove numbers after decimal
                      transcript<-gsub("\\..*", "", as.character(transcript))
                      return(transcript)
                    })))

data2_with_reads<-data2[which(data2[,6] > 0),]
genes2<-data2_with_reads$gene_id %>% unique()
isoform_max2<-data.frame(gene_id=genes2,
                    max_transcript=unlist(lapply(genes2, function(x){
                      transcript<-data2_with_reads[which(data2_with_reads$TPM==
                                                             max(filter(data2_with_reads, gene_id==x)$TPM) &
                                                           data2_with_reads$gene_id==x),
                                                   "transcript_id"] %>% unlist()
                      if(length(transcript)>0){
                        transcript<-transcript[1]
                      }
                      transcript<-gsub("\\..*", "", as.character(transcript))
                      return(transcript)
                    })))
K562_transcripts<-intersect(isoform_max1$max_transcript, 
                        isoform_max2$max_transcript) %>% 
    as.vector()

writeLines(K562_transcripts, "isoforms_K562.txt")

# map to 3'UTR
GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                         out_chain_name = "K562_3UTR.chain",
                         RNA_fragment = "three_prime_utr",
                         transcript_list = K562_transcripts,
                         alignment = "hg38")
getChainChrSize(chain = "K562_3UTR.chain",
                out_chr = "K562_3UTR.size")
ExtractTranscriptomeSequence(transcript_list = K562_transcripts,
                             ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                             genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                             RNA_fragment = "three_prime_utr",
                             exome_prefix = "K562_3UTR")
runCapR("K562_3UTR.fa")
processCapRout(CapR_outfile = "K562_3UTR.out",
               output_prefix = "K562_3UTR",
               chrom_size = "K562_3UTR.size",
               genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
               RNA_fragment = "three_prime_utr",
               chain = "K562_3UTR.chain")

# map to exome
GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                         out_chain_name = "K562_exome.chain",
                         transcript_list = K562_transcripts,
                         alignment = "hg38")
getChainChrSize(chain = "K562_exome.chain",
                out_chr = "K562_exome.size")
ExtractTranscriptomeSequence(transcript_list = K562_transcripts,
                             ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                             genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                             exome_prefix = "K562_exome")
runCapR("K562_exome.fa")
processCapRout(CapR_outfile = "K562_exome.out",
               output_prefix = "K562_exome",
               chrom_size = "K562_exome.size",
               genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
               RNA_fragment = "exon",
               chain = "K562_exome.chain")
