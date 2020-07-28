# G3BP1 binding relative to CapR predicted structure across 3'UTR in HepG2 cells

library(nearBynding)
library(Rsamtools)
library(rtracklayer)
library(dplyr)

transcript_list<-transcript_list<-readLines("get_isoforms_run_CapR/isoforms_HepG2.txt")

GenomeMappingToChainFile(genome_gtf = "Homo_sapiens.GRCh38.98.gtf",
                         out_chain_name = "HepG2_3UTR.chain",
                         RNA_fragment = "three_prime_utr",
                         transcript_list = transcript_list,
                         alignment = "hg38")
getChainChrSize(chain = "HepG2_3UTR.chain",
                out_chr = "HepG2_3UTR.size")

## input
in_bam<-sortBam("ENCFF700YZD.bam", "ENCFF700YZD_sorted")
CleanBAMtoBG(in_bam = in_bam,
             unwanted_chromosomes = "EBV",
             out_bedGraph = "ENCFF700YZD.bedGraph")
liftOverToExomicBG(input = "ENCFF700YZD.bedGraph",
                   chain = "HepG2_3UTR.chain",
                   chrom_size = "HepG2_3UTR.size",
                   output_bg = "ENCFF700YZD_liftOver.bedGraph")
runStereogeneOnCapR(chrom_size = "HepG2_3UTR.size",
                    name_config = "HepG2_3UTR.cfg",
                    input_prefix = "HepG2_3UTR",
                    protein_file = "ENCFF700YZD_liftOver.bedGraph")

## rep1
in_bam<-sortBam("ENCFF395MJO.bam", "ENCFF395MJO_sorted")
CleanBAMtoBG(in_bam = in_bam,
             unwanted_chromosomes = "EBV",
             out_bedGraph = "ENCFF395MJO.bedGraph")
liftOverToExomicBG(input = "ENCFF395MJO.bedGraph",
                   chain = "HepG2_3UTR.chain",
                   chrom_size = "HepG2_3UTR.size",
                   output_bg = "ENCFF395MJO_liftOver.bedGraph")
runStereogeneOnCapR(chrom_size = "HepG2_3UTR.size",
                    name_config = "HepG2_3UTR.cfg",
                    input_prefix = "HepG2_3UTR",
                    protein_file = "ENCFF395MJO_liftOver.bedGraph")

## rep2
in_bam<-sortBam("ENCFF288LEG.bam", "ENCFF288LEG_sorted")
CleanBAMtoBG(in_bam = in_bam,
             unwanted_chromosomes = "EBV",
             out_bedGraph = "ENCFF288LEG.bedGraph")
liftOverToExomicBG(input = "ENCFF288LEG.bedGraph",
                   chain = "HepG2_3UTR.chain",
                   chrom_size = "HepG2_3UTR.size",
                   output_bg = "ENCFF288LEG_liftOver.bedGraph")
runStereogeneOnCapR(chrom_size = "HepG2_3UTR.size",
                    name_config = "HepG2_3UTR.cfg",
                    input_prefix = "HepG2_3UTR",
                    protein_file = "ENCFF288LEG_liftOver.bedGraph")

## visualize binding
visualizeCapRStereogene(CapR_prefix = "HepG2_3UTR",
                        protein_file = c("ENCFF395MJO_liftOver", 
                                         "ENCFF288LEG_liftOver"),
                        protein_file_input = "ENCFF700YZD_liftOver",
                        out_file = "HepG2_3UTR_CapR_line",
                        legend = F)
visualizeCapRStereogene(CapR_prefix = "HepG2_3UTR",
                        protein_file = c("ENCFF395MJO_liftOver", 
                                         "ENCFF288LEG_liftOver"),
                        protein_file_input = "ENCFF700YZD_liftOver",
                        out_file = "HepG2_3UTR_CapR_heatmap",
                        heatmap = T)
