# NONO binding relative to G4 in K562 cells

library(nearBynding)
library(Rsamtools)
library(rtracklayer)
library(dplyr)

transcript_list<-transcript_list<-readLines("get_isoforms_run_CapR/isoforms_K562.txt")

GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                          out_chain_name = "K562_exome.chain",
                          RNA_fragment = "exon",
                          transcript_list = transcript_list,
                          alignment = "hg38")
getChainChrSize(chain = "K562_exome.chain",
                out_chr = "K562_exome.size")

############ NONO
###### rep 1
in_bam<-sortBam("ENCFF351ATM.bam", "ENCFF351ATM_sorted")
CleanBAMtoBG(in_bam = in_bam,
             unwanted_chromosomes = "EBV",
             out_bedGraph = "ENCFF351ATM.bedGraph")
liftOverToExomicBG(input = "ENCFF351ATM.bedGraph",
                   chain = "K562_exome.chain",
                   chrom_size = "K562_exome.size",
                   output_bg = "ENCFF351ATM_liftOver.bedGraph")
###### rep 2
in_bam<-sortBam("ENCFF888NYK.bam", "ENCFF888NYK_sorted")
CleanBAMtoBG(in_bam = in_bam,
             unwanted_chromosomes = "EBV",
             out_bedGraph = "ENCFF888NYK.bedGraph")
liftOverToExomicBG(input = "ENCFF888NYK.bedGraph",
                   chain = "K562_exome.chain",
                   chrom_size = "K562_exome.size",
                   output_bg = "ENCFF888NYK_liftOver.bedGraph")
###### input
in_bam<-sortBam("ENCFF701APQ.bam", "ENCFF701APQ_sorted")
CleanBAMtoBG(in_bam = in_bam,
             unwanted_chromosomes = "EBV",
             out_bedGraph = "ENCFF701APQ.bedGraph")
liftOverToExomicBG(input = "ENCFF701APQ.bedGraph",
                   chain = "K562_exome.chain",
                   chrom_size = "K562_exome.size",
                   output_bg = "ENCFF701APQ_liftOver.bedGraph")

###### G4 positions
# must first convert to hg38
G4_hg19<-read.table("GSE77282_K_hits.bed", header = F, stringsAsFactors = F)
G4_hg19 <- G4_hg19[,-c(7:length(G4_hg19))]
colnames(G4_hg19)<-c('chr','start','end','id','score','strand')
G4_hg19_GRanges <- with(G4_hg19, GRanges(chr, IRanges(start, end), 
                                         id=id, score=score, strand=strand))
hg19_hg38<-import.chain("hg19ToHg38.over.chain")
G4_hg38_GRanges<-liftOver(G4_hg19_GRanges, hg19_hg38) %>% unlist()
#export(G4_hg38_GRanges, "GSE77282_K_hg38.bed", "BED")
# separate by F and R strand
G4_hg38_GRanges_F<-G4_hg38_GRanges %>% filter(strand == "+")
G4_hg38_GRanges_R<-G4_hg38_GRanges %>% filter(strand == "-")
export(G4_hg38_GRanges_F, "GSE77282_K_hg38_F.bedgraph", "BedGraph")
export(G4_hg38_GRanges_R, "GSE77282_K_hg38_R.bedgraph", "BedGraph")

liftOverToExomicBG(input = c("GSE77282_K_hg38_F.bedgraph",
                             "GSE77282_K_hg38_R.bedgraph"),
                   chain = "K562_exome.chain",
                   chrom_size = "K562_exome.size",
                   output_bg = "GSE77282_K_hg38_liftOver.bedGraph")

## correlate
write_config(chrom_size = "K562_exome.size")
runStereogene(name_config = "config.cfg",
              track_files = c("GSE77282_K_hg38_liftOver.bedGraph",
                              "ENCFF351ATM_liftOver.bedGraph")) # protein track second
runStereogene(name_config = "config.cfg",
              track_files = c("GSE77282_K_hg38_liftOver.bedGraph",
                              "ENCFF888NYK_liftOver.bedGraph"))
runStereogene(name_config = "config.cfg",
              track_files = c("GSE77282_K_hg38_liftOver.bedGraph",
                              "ENCFF701APQ_liftOver.bedGraph"))

## visualize binding
visualizeStereogene(context_file = "GSE77282_K_hg38_liftOver",
                    protein_file = c("ENCFF351ATM_liftOver",
                                     "ENCFF888NYK_liftOver"),
                    protein_file_input = "ENCFF701APQ_liftOver",
                    out_file = "NONO_G4_K_heatmap",
                    heatmap = T)


