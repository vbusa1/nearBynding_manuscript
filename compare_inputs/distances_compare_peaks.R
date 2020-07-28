library(dplyr)
library(nearBynding)
library(ggplot2)

matched_names<-read.csv("match_CLIPper_Piranha_file_names.csv", stringsAsFactors = F) %>%
  mutate(Piranha_name = paste0(protein, "_", cell, "_rep", Piranha_rep))

pair_proteins<-read.csv("distances_compare_peaks.csv") # 3-col file with sample 1, sample 2, and distance

# make df where first 2 col are sample and category
proteins_categories<-data.frame(sample = c(matched_names$Piranha_name, matched_names$Accession, matched_names$file),
                                category = rep(paste0(matched_names$protein, "_", matched_names$cell), 3),
                                cell = rep(matched_names$cell, 3),
                                protein = rep(matched_names$protein, 3))

## plot HepG2
HepG2_proteins<-proteins_categories %>% filter(cell == "HepG2")
hold<-assessGrouping(pair_proteins, HepG2_proteins, measurement = "mean", output = "plot")
p <- assessGrouping(pair_proteins, HepG2_proteins, measurement = "mean", output = "KS.pvalue")
hold + annotate(geom = "text", y= 0.1, x = 6, label = paste("p-value = ", signif(p, 3)), hjust = "right") +
  theme(legend.position = c(0.8, 0.5), legend.text = element_text(size = 10), axis.text = element_text(size = 10))
ggsave("KS_mean_HepG2.pdf")


## plot K562
K562_proteins<-proteins_categories %>% filter(cell == "K562")
hold<-assessGrouping(pair_proteins, K562_proteins, measurement = "mean", output = "plot")
p <- assessGrouping(pair_proteins, K562_proteins, measurement = "mean", output = "KS.pvalue")
hold + annotate(geom = "text", y= 0.1, x = 6, label = paste("p-value = ", signif(p, 3)), hjust = "right") +
  theme(legend.position = c(0.8, 0.5), legend.text = element_text(size = 10), axis.text = element_text(size = 10))
ggsave("KS_mean_K562.pdf")

## plot MDS
pair_proteins <- filter(pair_proteins, distance > 0)
pair_proteins<-pair_proteins[order(pair_proteins[,2], pair_proteins[,1]), ]
pair_proteins<-reshape(pair_proteins, idvar = "Var1", timevar = "Var2", v.names = "distance", direction = "wide")
pair_proteins<-rbind(pair_proteins[nrow(pair_proteins),], pair_proteins[1:(nrow(pair_proteins)-1),])

rownames(pair_proteins)<-lapply(pair_proteins[,1], function(x){
    if(x %in% matched_names$file){
        peak_caller<-"track"
        i<-match(x, matched_names$file)
    }else if(x %in% matched_names$Accession){
        peak_caller<-"CLIPper"
        i<-match(x, matched_names$Accession)
    }else{
        peak_caller<-"Piranha"
        i<-match(x, matched_names$Piranha_name)
    }
    cell<-matched_names[i, "cell"]
    protein<-matched_names[i, "protein"]
    return(paste0(protein, "_", cell, "_", peak_caller))
}) %>% unlist()
cell<-lapply(pair_proteins[,1], function(x){
    if(x %in% matched_names$file){
        cell<-matched_names[match(x, matched_names$file), "cell"]
    }else if(x %in% matched_names$Accession){
        cell<-matched_names[match(x, matched_names$Accession), "cell"]
    }else{
        cell<-matched_names[match(x, matched_names$Piranha_name), "cell"]
    }
    return(cell)
}) %>% unlist()
protein<-lapply(pair_proteins[,1], function(x){
    if(x %in% matched_names$file){
        protein<-matched_names[match(x, matched_names$file), "protein"]
    }else if(x %in% matched_names$Accession){
        protein<-matched_names[match(x, matched_names$Accession), "protein"]
    }else{
        protein<-matched_names[match(x, matched_names$Piranha_name), "protein"]
    }
    return(protein)
}) %>% unlist()
peak_caller<-lapply(pair_proteins[,1], function(x){
    if(x %in% matched_names$file){
        peak_caller<-"track"
    }else if(x %in% matched_names$Accession){
        peak_caller<-"CLIPper"
    }else{
        peak_caller<-"Piranha"
    }
    return(peak_caller)
}) %>% unlist()
pair_proteins<-pair_proteins[,2:ncol(pair_proteins)]
colnames(pair_proteins)<-rownames(pair_proteins)
pair_proteins<-as.dist(pair_proteins)
fit <- cmdscale(pair_proteins, k=2) # k is the number of dim
plot<-data.frame(x = fit[,1],
                 y = fit[,2],
                 cell = cell,
                 protein = protein,
                 peak_caller = peak_caller)
# make sure colors are diverse
toMatch<-c("blue", "red", "gold", "orange", "magenta", "orchid", "green",
           "pink", "violet", "brown", "aquamarine", "cyan")
Color = colors()[sample(grep(paste(toMatch,collapse="|"), grDevices::colors()), length(unique(protein)))]

# plot
ggplot(plot, aes(x = x, y = y, shape = peak_caller, fill = protein)) +
    facet_wrap(~cell, nrow = 2, scales = "free_y") +
    scale_shape_manual(values = c(22, 23, 21)) +
    geom_point(alpha = .7, size = 3, color = "black") +
    theme_classic() +
    scale_fill_manual(values = Color, guide = F) +
    labs(x = "MDS_1", y = "MDS_2") +
    theme(strip.text = element_text(size = 12))
ggplot(plot, aes(x = x, y = y, shape = peak_caller, color = protein)) +
    facet_wrap(~cell, nrow = 2) +
    scale_shape_manual(values = c(15, 18, 16)) +
    geom_point(alpha = .7, size = 3) +
    theme_classic() +
    scale_color_manual(values = Color)
