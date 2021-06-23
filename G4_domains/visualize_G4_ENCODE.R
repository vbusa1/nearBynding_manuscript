library(tidyverse)
library(plyr)
library(nearBynding)
library(ggplot2)
library(ggpubr)
library(lsr)
library(gplots)
library(ggforce)
setwd("/Users/veronica/Documents/GitHub/RBP_binding/hold_files/all_dist_files/")

files<-list.files()
files<-files[grepl("_GSE77282", files)]
Accession<-lapply(files, function(x){substr(x, nchar(x)-24, nchar(x)-14)}) %>% 
    unlist()
file_info<-read.csv("../file_information.csv", stringsAsFactors = F)

# for(test_protein in unique(file_info$protein)){
#     cell_lines<- filter(file_info, protein == test_protein)[, "cell"] %>% unique()
#     for(test_cell in cell_lines){
#         replicates<-filter(file_info, protein == test_protein,
#                            cell == test_cell,
#                            input == F)[, "file"]
#         if(length(replicates) != 2){next}
#         input<-filter(file_info, protein == test_protein,
#                       cell == test_cell,
#                       input == T)[, "file"]
#         if(length(input) != 1){next}
#         visualizeStereogene(out_file = paste0("../all_G4_visuals/", test_protein,
#                                               "_", test_cell, "_unspliced_line"),
#                             context_file = paste0(test_cell, "_GSE77282_K_hg38_liftOver"),
#                             protein_file = c(paste0(replicates[1], "_liftOver"),
#                                              paste0(replicates[2], "_liftOver")),
#                             protein_file_input = paste0(input, "_liftOver"),
#                             legend = F)
#         visualizeStereogene(out_file = paste0("../all_G4_visuals/", test_protein,
#                                               "_", test_cell, "_unspliced_heatmap"),
#                             context_file = paste0(test_cell, "_GSE77282_K_hg38_liftOver"),
#                             protein_file = c(paste0(replicates[1], "_liftOver"),
#                                              paste0(replicates[2], "_liftOver")),
#                             protein_file_input = paste0(input, "_liftOver"),
#                             heatmap = T)
#     }
# }

range = c(-100, 100)
amplitude<-as.data.frame(matrix(NA, nrow = nrow(file_info), ncol = 5))
n<-1
for(test_protein in unique(file_info$protein)){
    cell_lines<- filter(file_info, protein == test_protein)[, "cell"] %>% unique()
    for(test_cell in cell_lines){
        replicates<-filter(file_info, protein == test_protein,
                           cell == test_cell,
                           input == F)[, "file"]
        if(length(replicates) != 2){next}
        input<-filter(file_info, protein == test_protein,
                      cell == test_cell,
                      input == T)[, "file"]
        if(length(input) != 1){next}
        dist_1<-read.table(paste0("GSE77282_K_hg38_", test_cell, "_liftOver~", replicates[1], "_liftOver.dist"), header = T)%>% 
            filter(range[1] <= x, x <= range[2])
        dist_2<-read.table(paste0("GSE77282_K_hg38_", test_cell, "_liftOver~", replicates[2], "_liftOver.dist"), header = T)%>% 
            filter(range[1] <= x, x <= range[2])
        dist_input<- read.table(paste0("GSE77282_K_hg38_", test_cell, "_liftOver~", input, "_liftOver.dist"), header = T)%>% 
            filter(range[1] <= x, x <= range[2])
        
        dist<-rowMeans(cbind(dist_1$Fg, dist_2$Fg)) - dist_input$Fg
        
        get_signal<-max(abs(max(dist)), abs(min(dist)))
        get_pos<-(range[1]:range[2])[which(abs(dist) == get_signal)]
        
        amplitude[n,]<-c(test_protein, test_cell, max(dist), min(dist), get_pos)
        n<-n+1
    }
}
amplitude<-amplitude[complete.cases(amplitude),]

amplitude$signal <- lapply(1:nrow(amplitude), function(x){
    if(abs(as.numeric(amplitude[x,"V3"])) > abs(as.numeric(amplitude[x,"V4"]))){
        return(amplitude[x,"V3"])
    }else{
        return(amplitude[x,"V4"])
    }
}) %>% unlist()
amplitude$signal<-as.numeric(amplitude$signal)
hist(amplitude$signal, breaks = 100)


amplitude_new<-as.data.frame(matrix(NA, nrow = nrow(file_info), ncol = 5))
n<-1
for(test_protein in unique(file_info$protein)){
    cell_lines<- filter(file_info, protein == test_protein)[, "cell"] %>% unique()
    for(test_cell in cell_lines){
        replicates<-filter(file_info, protein == test_protein,
                           cell == test_cell,
                           input == F)[, "file"]
        if(length(replicates) != 2){next}
        input<-filter(file_info, protein == test_protein,
                      cell == test_cell,
                      input == T)[, "file"]
        if(length(input) != 1){next}
        dist_1<-read.table(paste0(test_cell, "_GSE77282_K_hg38_liftOver~", replicates[1], "_liftOver.dist"), header = T)%>% 
            filter(range[1] <= x, x <= range[2])
        dist_2<-read.table(paste0(test_cell, "_GSE77282_K_hg38_liftOver~", replicates[2], "_liftOver.dist"), header = T)%>% 
            filter(range[1] <= x, x <= range[2])
        dist_input<- read.table(paste0(test_cell, "_GSE77282_K_hg38_liftOver~", input, "_liftOver.dist"), header = T)%>% 
            filter(range[1] <= x, x <= range[2])
        
        dist<-rowMeans(cbind(dist_1$Fg, dist_2$Fg)) - dist_input$Fg
        
        get_signal<-max(abs(max(dist)), abs(min(dist)))
        get_pos<-(range[1]:range[2])[which(abs(dist) == get_signal)]
        
        amplitude_new[n,]<-c(test_protein, test_cell, max(dist), min(dist), get_pos)
        n<-n+1
    }
}
amplitude_new<-amplitude_new[complete.cases(amplitude_new),]

amplitude_new$signal <- lapply(1:nrow(amplitude_new), function(x){
    if(abs(as.numeric(amplitude_new[x,"V3"])) > abs(as.numeric(amplitude_new[x,"V4"]))){
        return(amplitude_new[x,"V3"])
    }else{
        return(amplitude_new[x,"V4"])
    }
}) %>% unlist()
amplitude_new$signal<-as.numeric(amplitude_new$signal)
hist(amplitude_new$signal, breaks = 100)

# RBP_annotation<-read.csv("../RBP_annotation.csv")
# domains_amplitudes<-merge(amplitude, RBP_annotation,
#                           by.x = "V1", by.y = "Protein")
# RG_proteins<-read.delim("../RG_proteins.txt", header = F)
# domains_amplitudes$RG<-domains_amplitudes$V1 %in% RG_proteins$V3
# 
# cohensD(domains_amplitudes[which(domains_amplitudes$RG == T), "signal"],
#         domains_amplitudes[which(domains_amplitudes$RG == F), "signal"])
# t.test(domains_amplitudes[which(domains_amplitudes$RG == T), "signal"],
#        domains_amplitudes[which(domains_amplitudes$RG == F), "signal"])$p.value
# 
# stats_tests<-matrix(NA, nrow = 13, ncol = 4)
# colnames(stats_tests)<-c("domain", "G4 signal difference", "effect size", "p-value")
# for(n in 6:18){
#     p<-t.test(domains_amplitudes[which(domains_amplitudes[,n] == T), "signal"],
#               domains_amplitudes[which(domains_amplitudes[,n] == F), "signal"])$p.value %>%
#         signif(3)
#     m<-t.test(domains_amplitudes[which(domains_amplitudes[,n] == T), "signal"],
#                  domains_amplitudes[which(domains_amplitudes[,n] == F), "signal"])$estimate
#     d<-cohensD(domains_amplitudes[which(domains_amplitudes[,n] == T), "signal"],
#                domains_amplitudes[which(domains_amplitudes[,n] == F), "signal"]) %>%
#         signif(3)
#     c<-colnames(domains_amplitudes[n])
#     stats_tests[n-5,]<-c(c, signif(m[1]-m[2], 3), d, p)
# }
# write.csv(stats_tests, "../stats_domains_G4.csv", row.names = F, quote = F)
#
# gather_domains<-domains_amplitudes %>% gather("domain", "present", 6:18)
# symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), 
#                     symbols = c("****", "***", "**", "*", ".", "ns"))
# ggplot(gather_domains, aes(x = domain, y = as.numeric(signal),
#                            color = paste(present, V2))) +
#     geom_sina(size = 1, maxwidth = .7) +
#     theme_classic() +
#     labs(y = "G4 Density Signal", x = "Protein Domain") +
#     stat_compare_means(aes(domain, as.numeric(signal), group = present),
#                        label = "p.signif",
#                        method = "t.test",
#                        symnum.args = symnum.args) +
#     ylim(-5, 13) +
#     theme(axis.text.x = element_text(angle = 65, vjust = .5, hjust=.4, size = 12),
#           text = element_text(size = 14),
#           axis.text.y = element_text(size = 12)) +
#     scale_color_manual(values = c("red3", "brown1", "blue3", "dodgerblue"))
# ggsave("../visualize_G4_domains_unspliced.pdf")
# 
# gather_domains_HepG2<-filter(gather_domains, V2 == "HepG2")
# ggplot(gather_domains_HepG2, aes(x = domain, y = as.numeric(signal),
#                            color = present)) +
#     geom_sina() +
#     theme_classic() +
#     labs(y = "G4 Density Signal", x = "Protein Domain") +
#     stat_compare_means(aes(domain, as.numeric(signal), group = present),
#                        label = "p.signif",
#                        method = "t.test",
#                        symnum.args = symnum.args) +
#     ylim(-5, 13) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.4, size = 12),
#           text = element_text(size = 14),
#           axis.text.y = element_text(size = 12)) +
#     scale_color_manual(values = c("red", "blue")) +
#     ggtitle("HepG2 domain analysis")
# ggsave("../visualize_G4_domains_unspliced_HepG2.pdf")
# 
# gather_domains_K562<-filter(gather_domains, V2 == "K562")
# ggplot(gather_domains_K562, aes(x = domain, y = as.numeric(signal),
#                                  color = present)) +
#     geom_sina() +
#     theme_classic() +
#     labs(y = "G4 Density Signal", x = "Protein Domain") +
#     stat_compare_means(aes(domain, as.numeric(signal), group = present),
#                        label = "p.signif",
#                        method = "t.test",
#                        symnum.args = symnum.args) +
#     ylim(-5, 13) +
#     theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.4, size = 12),
#           text = element_text(size = 14),
#           axis.text.y = element_text(size = 12)) +
#     scale_color_manual(values = c("red", "blue")) +
#     ggtitle("K562 domain analysis")
# ggsave("../visualize_G4_domains_unspliced_K562.pdf")
# 
# 
# heatmap_plot<-as.data.frame(matrix(NA, nrow = nrow(file_info), ncol = 203))
# n<-1
# range = c(-100, 100)
# 
# for(test_protein in unique(file_info$protein)){
#     cell_lines<- filter(file_info, protein == test_protein)[, "cell"] %>% unique()
#     for(test_cell in cell_lines){
#         replicates<-filter(file_info, protein == test_protein,
#                            cell == test_cell,
#                            input == F)[, "file"]
#         if(length(replicates) != 2){next}
#         input<-filter(file_info, protein == test_protein,
#                       cell == test_cell,
#                       input == T)[, "file"]
#         if(length(input) != 1){next}
#         dist_1<-read.table(paste0(test_cell, "_GSE77282_K_hg38_liftOver~", replicates[1], "_liftOver.dist"), header = T)%>% 
#             filter(range[1] <= x, x <= range[2])
#         dist_2<-read.table(paste0(test_cell, "_GSE77282_K_hg38_liftOver~", replicates[2], "_liftOver.dist"), header = T)%>% 
#             filter(range[1] <= x, x <= range[2])
#         dist_input<- read.table(paste0(test_cell, "_GSE77282_K_hg38_liftOver~", input, "_liftOver.dist"), header = T)%>% 
#             filter(range[1] <= x, x <= range[2])
#         
#         dist<-rowMeans(cbind(dist_1$Fg, dist_2$Fg)) - dist_input$Fg
#         
#         heatmap_plot[n,]<-c(test_protein, test_cell, dist)
#         n<-n+1
#     }
# }
# 
# heatmap_plot<-merge(amplitude, heatmap_plot, 
#                     by.x = c("V1", "V2"), by.y = c("V1", "V2"))
# heatmap_plot<-heatmap_plot[complete.cases(heatmap_plot),]
# heatmap_plot$signal<-as.numeric(heatmap_plot$signal)
# heatmap_plot_sort<-heatmap_plot[order(heatmap_plot$signal),]
# 
# G4_lit<-read.csv("../G4_protein_PMID.csv")
# G4_lit_T<-G4_lit$protein[G4_lit$PMID != ""]
# 
# heatmap_plot_sort_HepG2<-heatmap_plot_sort %>% filter(V2 == "HepG2")
# heatmap_final_HepG2<-heatmap_plot_sort_HepG2[,6:ncol(heatmap_plot_sort_HepG2)] %>%
#     sapply(as.numeric)
# row.names(heatmap_final_HepG2)<-heatmap_plot_sort_HepG2$V1
# col<-unlist(lapply(row.names(heatmap_final_HepG2), function(x){
#     if(x %in% G4_lit_T){return("red")}else{return("white")}
# }))
# jpeg("../visualize_G4_HepG2_unspliced.jpeg", height = 14, width = 6, units = "in", res = 200)
# heatmap.2(heatmap_final_HepG2,
#           Rowv = NA, Colv = NA, col = redblue(256),
#           keysize = .1, RowSideColors = col,
#           dendrogram = "none", labCol = FALSE, trace = "none",
#           density.info = "none", symkey = FALSE, key = F,
#           colsep = c(100, 101), sepcolor = "lightgrey",
#           margins = c(1, 5),
#           breaks = seq(-max(heatmap_plot$signal), 
#                        max(heatmap_plot$signal), 
#                        length.out = 257))
# dev.off()
# 
# 
# heatmap_plot_sort_K562<-heatmap_plot_sort %>% filter(V2 == "K562")
# heatmap_final_K562<-heatmap_plot_sort_K562[,6:ncol(heatmap_plot_sort_K562)] %>%
#     sapply(as.numeric)
# row.names(heatmap_final_K562)<-heatmap_plot_sort_K562$V1
# col<-unlist(lapply(row.names(heatmap_final_K562), function(x){
#     if(x %in% G4_lit_T){return("red")}else{return("white")}
# }))
# jpeg("../visualize_G4_K562_unspliced.jpeg", height = 14, width = 6, units = "in", res = 200)
# heatmap.2(heatmap_final_K562,
#           Rowv = NA, Colv = NA, col = redblue(256),
#           keysize = .1, RowSideColors = col,
#           dendrogram = "none", labCol = FALSE, trace = "none",
#           density.info = "none", symkey = FALSE, key = F,
#           colsep = c(100, 101), sepcolor = "lightgrey",
#           margins = c(1, 5),
#           breaks = seq(-max(heatmap_plot$signal), 
#                        max(heatmap_plot$signal), 
#                        length.out = 257))
# dev.off()


all_amplitude<-merge(amplitude, amplitude_new, by= c("V1", "V2"))

ggplot(all_amplitude, aes(x = signal.x,
                          y =  signal.y)) +
    geom_abline(slope=1, intercept = 0, color = "grey") +
    geom_point(aes(color = V2), alpha = .7, size = 3) +
    theme_classic() +
    xlim(-5, 12) + ylim(-5, 12) +
    labs(y = "maximal G4 unspliced correlation",
         x = "maximal G4 spliced correlation") +
    theme(text = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size = 16)) +
    scale_color_manual(values = c("red", "blue"))
# ggsave("../G4/G4_spliced_unspliced_cell.jpeg")

# R2 = 0.81

# protein_functions<-read.csv("../13059_2020_1982_MOESM2_ESM.csv")
# all_info<-merge(all_amplitude, protein_functions, by.x = "V1", by.y = "protein")
# ggplot(all_info, aes(x = signal.x,
#                           y =  signal.y)) +
#     geom_abline(slope=1, intercept = 0, color = "grey") +
#     geom_point(aes(color = factor(Splicing.regulation)), alpha = .7, size = 3) +
#     theme_classic() +
#     xlim(-5, 12) + ylim(-5, 12) +
#     labs(y = "maximal G4 unspliced correlation",
#          x = "maximal G4 spliced correlation") +
#     theme(text = element_text(size = 12))
intron<-c("BCCIP", "CSTF2", "CSTF2T", "EWSR1", "FAM120A", "FUBP3", "FUS", "HNRNPA1", "HNRNPC",
          "HNRNPK", "HNRNPL", "HNRNPM", "HNRNPU", "HNRNPUL1", "KHDRBS1", "KHSRP", "MATR3", "NCBP2",
          "NONO", "PCBP2", "QKI", "RBFOX2", "SAFB", "SAFB2", "SFPQ", "SUGP2", "TAF15", "TIA1", "TIAL1",
          "AQR", "BUD13", "CSTF2T", "EFTUD2", "EWSR1", "FAM120A", "KHSRP", "PRPF4", "PRPF8",
          "RBFOX2", "RBM22", "SF3B4", "TIA1", "TIAL1", "U2AF1", "U2AF2", "SLTM", "GTF2F1", "PTBP1", "POLR2G",
          "ILF3", "EXOSC5", "FKBP4", "NSUN2", "XPO5", "XRCC6", "EIF4G2", "TARDBP", "DDX52", "AGGF1", "YWHAG",
          "STAU2", "AKAP8L", "NKRF", "DDX59", "LSM11", "CDC40", "DROSHA", "DDX42", "XRN2", "FASSTKD2", 
          "DDX51", "SRSF7", "DHX30", "ZRANB2", "WRN", "TBRG4", "SMNDC1", "PIL4", "FTO")
all_amplitude$intronic<-all_amplitude$V1 %in% intron
all_amplitude$intronic<-mapvalues(all_amplitude$intronic, from = c(TRUE, FALSE), to = c("intronic", "non-intronic"))
ggplot(all_amplitude, aes(x = signal.x,
                          y =  signal.y)) +
    facet_wrap(~intronic) +
    geom_abline(slope=1, intercept = 0, color = "grey") +
    geom_point(aes(color = intronic), size = 3,
               show.legend = F) +
    theme_classic() +
    xlim(-5, 12) + ylim(-5, 12) +
    labs(y = "maximal G4 unspliced correlation",
         x = "maximal G4 spliced correlation") +
    theme(text = element_text(size = 20))+
    stat_regline_equation(label.y = -4, label.x = 6, 
                          aes(label = ..rr.label..),
                          size = 5) +
    scale_color_manual(values = c("red", "blue"))
# ggsave("../G4/G4_spliced_unspliced_bindingsite.jpeg")
ggplot(all_amplitude, aes(x = signal.x,
                          y =  signal.y)) +
    geom_abline(slope = 1, color = "grey") +
    geom_point(aes(color = intronic),
               alpha = .7, size = 3) +
    theme_classic() +
    xlim(-5, 12) + ylim(-5, 12) +
    labs(y = "maximal G4 unspliced correlation",
         x = "maximal G4 spliced correlation") +
    theme(text = element_text(size = 20),
          legend.title = element_blank(),
          legend.text = element_text(size = 16))+
    stat_regline_equation(label.y = -4, label.x = 6, 
                          aes(label = ..rr.label..),
                          size = 5) +
    scale_color_manual(values = c("red", "blue"))
# ggsave("../G4/G4_spliced_unspliced_bindingsite_combined.jpeg")

all_amplitude$residual<-lm(all_amplitude$signal.y~all_amplitude$signal.x)$residuals
ggplot(all_amplitude, aes(x = intronic, y = residual^2)) +
    geom_sina(aes(color = intronic), pch = 1, show.legend = F, size = 2) +
    geom_boxplot(alpha = 0) +
    theme_classic() +
    xlab(NULL) +
    ylab(expression(residual^2)) +
    theme(text = element_text(size = 16)) +
    stat_compare_means(aes(label = ..p.signif..), method = "t.test",
                       comparisons = list(c("intronic", "non-intronic"))) +
    scale_color_manual(values = c("red", "blue"))
# ggsave("../G4/G4_spliced_unspliced_residual.jpeg")

ggplot(filter(all_amplitude, signal.x > 2.5),
       aes(x = as.numeric(V5.x),
           y = as.numeric(V5.y))) +
    facet_wrap(~intronic) +
    geom_abline(slope=1, intercept = 0, color = "grey") +
    geom_point(aes(color = intronic),
               show.legend = F) +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          text = element_text(size = 14)) +
    stat_regline_equation(label.y = 5, label.x = 30, 
                          aes(label = ..rr.label..),
                          size = 5) +
    ylim(0, 70) + xlim (0, 70) +
    labs(y = "position max G4 unspliced corr",
         x = "position max G4 spliced corr") +
    scale_color_manual(values = c("red", "blue"))
# ggsave("../G4/G4_spliced_unspliced_signal_pos.jpeg")

ggplot(all_amplitude,
       aes(x = as.numeric(V5.x),
           y = as.numeric(V5.y))) +
    geom_abline(slope=1, intercept = 0, color = "grey") +
    geom_point(aes(color = intronic)) +
    theme_classic() +
    theme(axis.text = element_text(size = 16),
          legend.title = element_blank(),
          text = element_text(size = 14)) +
    stat_regline_equation(label.y = -95, label.x = -50, 
                          aes(label = ..rr.label..),
                          size = 5) +
    #ylim(0, 70) + xlim (0, 70) +
    labs(y = "position max G4 unspliced corr",
         x = "position max G4 spliced corr") +
    scale_color_manual(values = c("red", "blue"))
# ggsave("../G4/G4_spliced_unspliced_signal_pos_all.jpeg")
