library(tidyverse)
library(nearBynding)
library(ggpubr)
library(gplots)
library(lsr)

files<-list.files() # list all files in directory containing all *dist files
files<-files[grepl("GSE77282", files)] # only grab G4-correlated files
Accession<-lapply(files, function(x){substr(x, nchar(x)-24, nchar(x)-14)}) %>% 
    unlist()
file_info<-read.csv("file_information.csv", stringsAsFactors = F)

amplitude<-as.data.frame(matrix(NA, nrow = nrow(file_info), ncol = 4))
n<-1
range = c(-100, 100)

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
        
        amplitude[n,]<-c(test_protein, test_cell, max(dist), min(dist))
        n<-n+1
    }
}


amplitude<-amplitude[complete.cases(amplitude),]
# get maximal amplitude
amplitude$signal <- lapply(1:nrow(amplitude), function(x){
    if(abs(as.numeric(amplitude[x,"V3"])) > abs(as.numeric(amplitude[x,"V4"]))){
        return(amplitude[x,"V3"])
    }else{
        return(amplitude[x,"V4"])
    }
}) %>% unlist()

amplitude$signal<-as.numeric(amplitude$signal)

RBP_annotation<-read.csv("RBP_annotation.csv")
domains_amplitudes<-merge(amplitude, RBP_annotation,
                          by.x = "V1", by.y = "Protein")
RG_proteins<-read.delim("RG_proteins.txt", header = F)
domains_amplitudes$RG<-domains_amplitudes$V1 %in% RG_proteins$V3

stats_tests<-matrix(NA, nrow = 13, ncol = 4)
colnames(stats_tests)<-c("domain", "G4 signal difference", "effect size", "p-value")
for(n in 6:18){
    p<-t.test(domains_amplitudes[which(domains_amplitudes[,n] == T), "signal"],
              domains_amplitudes[which(domains_amplitudes[,n] == F), "signal"])$p.value %>%
        signif(3)
    m<-t.test(domains_amplitudes[which(domains_amplitudes[,n] == T), "signal"],
                 domains_amplitudes[which(domains_amplitudes[,n] == F), "signal"])$estimate
    d<-cohensD(domains_amplitudes[which(domains_amplitudes[,n] == T), "signal"],
               domains_amplitudes[which(domains_amplitudes[,n] == F), "signal"]) %>%
        signif(3)
    c<-colnames(domains_amplitudes[n])
    stats_tests[n-5,]<-c(c, signif(m[1]-m[2], 3), d, p)
}
write.csv(stats_tests, "stats_domains_G4.csv", row.names = F, quote = F)

gather_domains<-domains_amplitudes %>% gather("domain", "present", 6:18)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", "ns"))
ggplot(gather_domains, aes(x = domain, y = as.numeric(signal),
                           color = present)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "G4 Density Signal", x = "Protein Domain") +
    stat_compare_means(aes(domain, as.numeric(signal), group = present),
                       label = "p.signif",
                       method = "t.test",
                       symnum.args = symnum.args) +
    ylim(-5, 13) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.4, size = 12),
          text = element_text(size = 14),
          axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = c("red", "blue")) +
    ggtitle("Pooled domain analysis")
ggsave("visualize_G4_domains.pdf")

gather_domains_HepG2<-filter(gather_domains, V2 == "HepG2")
ggplot(gather_domains_HepG2, aes(x = domain, y = as.numeric(signal),
                           color = present)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "G4 Density Signal", x = "Protein Domain") +
    stat_compare_means(aes(domain, as.numeric(signal), group = present),
                       label = "p.signif",
                       method = "t.test",
                       symnum.args = symnum.args) +
    ylim(-5, 13) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.4, size = 12),
          text = element_text(size = 14),
          axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = c("red", "blue")) +
    ggtitle("HepG2 domain analysis")
ggsave("visualize_G4_domains_HepG2.pdf")

gather_domains_K562<-filter(gather_domains, V2 == "K562")
ggplot(gather_domains_K562, aes(x = domain, y = as.numeric(signal),
                                 color = present)) +
    geom_boxplot() +
    theme_classic() +
    labs(y = "G4 Density Signal", x = "Protein Domain") +
    stat_compare_means(aes(domain, as.numeric(signal), group = present),
                       label = "p.signif",
                       method = "t.test",
                       symnum.args = symnum.args) +
    ylim(-5, 13) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=.4, size = 12),
          text = element_text(size = 14),
          axis.text.y = element_text(size = 12)) +
    scale_color_manual(values = c("red", "blue")) +
    ggtitle("K562 domain analysis")
ggsave("visualize_G4_domains_K562.pdf")


heatmap_plot<-as.data.frame(matrix(NA, nrow = nrow(file_info), ncol = 203))
n<-1
range = c(-100, 100)

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
        
        heatmap_plot[n,]<-c(test_protein, test_cell, dist)
        n<-n+1
    }
}

heatmap_plot<-merge(amplitude, heatmap_plot, 
                    by.x = c("V1", "V2"), by.y = c("V1", "V2"))
heatmap_plot<-heatmap_plot[complete.cases(heatmap_plot),]
heatmap_plot$signal<-as.numeric(heatmap_plot$signal)
heatmap_plot_sort<-heatmap_plot[order(heatmap_plot$signal),]

G4_lit<-read.csv("G4_protein_PMID.csv")
G4_lit_T<-G4_lit$protein[!is.na(G4_lit$PMID)]

heatmap_plot_sort_HepG2<-heatmap_plot_sort %>% filter(V2 == "HepG2")
heatmap_final_HepG2<-heatmap_plot_sort_HepG2[,6:ncol(heatmap_plot_sort_HepG2)] %>%
    sapply(as.numeric)
row.names(heatmap_final_HepG2)<-heatmap_plot_sort_HepG2$V1
col<-unlist(lapply(row.names(heatmap_final_HepG2), function(x){
    if(x %in% G4_lit_T){return("red")}else{return("white")}
}))
jpeg("visualize_G4_HepG2.jpeg", height = 14, width = 6, units = "in", res = 200)
heatmap.2(heatmap_final_HepG2,
          Rowv = NA, Colv = NA, col = redblue(256),
          keysize = .1, RowSideColors = col,
          dendrogram = "none", labCol = FALSE, trace = "none",
          density.info = "none", symkey = FALSE, key = F,
          colsep = c(100, 101), sepcolor = "lightgrey",
          margins = c(1, 5),
          breaks = seq(-max(heatmap_plot$signal), 
                       max(heatmap_plot$signal), 
                       length.out = 257))
dev.off()


heatmap_plot_sort_K562<-heatmap_plot_sort %>% filter(V2 == "K562")
heatmap_final_K562<-heatmap_plot_sort_K562[,6:ncol(heatmap_plot_sort_K562)] %>%
    sapply(as.numeric)
row.names(heatmap_final_K562)<-heatmap_plot_sort_K562$V1
col<-unlist(lapply(row.names(heatmap_final_K562), function(x){
    if(x %in% G4_lit_T){return("red")}else{return("white")}
}))
jpeg("visualize_G4_K562.jpeg", height = 14, width = 6, units = "in", res = 200)
heatmap.2(heatmap_final_K562,
          Rowv = NA, Colv = NA, col = redblue(256),
          keysize = .1, RowSideColors = col,
          dendrogram = "none", labCol = FALSE, trace = "none",
          density.info = "none", symkey = FALSE, key = F,
          colsep = c(100, 101), sepcolor = "lightgrey",
          margins = c(1, 5),
          breaks = seq(-max(heatmap_plot$signal), 
                       max(heatmap_plot$signal), 
                       length.out = 257))
dev.off()
