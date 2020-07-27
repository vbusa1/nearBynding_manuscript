# determine the most highly expressed transcripts for every gene of HepG2 cells
# data from ENCODE

library(nearBynding)
library(dplyr)

data1<-read.table("ENCFF550BJR.tsv", header = T, stringsAsFactors = F)
data2<-read.table("ENCFF957EGM.tsv", header = T, stringsAsFactors = F)

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
HepG2_transcripts<-intersect(isoform_max1$max_transcript, 
                        isoform_max2$max_transcript) %>% 
    as.vector()

writeLines(HepG2_transcripts, "isoforms_HepG2.txt")

# map to 3'UTR
GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                         out_chain_name = "HepG2_3UTR.chain",
                         RNA_fragment = "three_prime_utr",
                         transcript_list = HepG2_transcripts,
                         alignment = "hg38")
getChainChrSize(chain = "HepG2_3UTR.chain",
                out_chr = "HepG2_3UTR.size")
ExtractTranscriptomeSequence(transcript_list = HepG2_transcripts,
                             ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                             genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                             RNA_fragment = "three_prime_utr",
                             exome_prefix = "HepG2_3UTR")
runCapR("HepG2_3UTR.fa")
processCapRout(CapR_outfile = "HepG2_3UTR.out",
               output_prefix = "HepG2_3UTR",
               chrom_size = "HepG2_3UTR.size",
               genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
               RNA_fragment = "three_prime_utr",
               chain = "HepG2_3UTR.chain")

# map to exome
GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                         out_chain_name = "HepG2_exome.chain",
                         transcript_list = HepG2_transcripts,
                         alignment = "hg38")
getChainChrSize(chain = "HepG2_exome.chain",
                out_chr = "HepG2_exome.size")
ExtractTranscriptomeSequence(transcript_list = HepG2_transcripts,
                             ref_genome = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                             genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                             exome_prefix = "HepG2_exome")
runCapR("HepG2_exome.fa")
processCapRout(CapR_outfile = "HepG2_exome.out",
               output_prefix = "HepG2_exome",
               chrom_size = "HepG2_exome.size",
               genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
               RNA_fragment = "exon",
               chain = "HepG2_exome.chain")
