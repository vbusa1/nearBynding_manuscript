setwd("/Users/veronica/Documents/GitHub/RBP_binding/hold_files/compare_peaks/")
library(nearBynding)
library(tidyverse)
library(GGally)

# # check for completeness of files, and run incomplete files again
# temp<-Sys.glob("CLIPper_peaks/*[A-Z].bedGraph")
# temp<-lapply(temp, function(x){substr(x, nchar(x)-19, nchar(x)-9)}) %>% unlist()
#
# num_dist<-table(unlist(lapply(temp, function(x){
#   return(rep(x, sum(grepl(x, Sys.glob("CLIPper_peaks/*.dist")))))})))
# rerun_dist<-filter(as.data.frame(num_dist), Freq < 6)
#
# cat("#!bin/bash",
#     file = "submit_to_JHPCE.sh",
#     sep = "\n")
# for(file in rerun_dist$Var1){
#   cat(paste0("qsub -cwd ", file, ".sh"),
#       file = "submit_to_JHPCE.sh",
#       sep = "\n",
#       append = T)
# }
#
# temp<-Sys.glob("Piranha_peaks/*[1-2].bedGraph")
# temp<-lapply(temp, function(x){substr(x, 15, nchar(x)-9)}) %>% unlist()
#
# num_dist<-table(unlist(lapply(temp, function(x){
#   return(rep(x, sum(grepl(x, Sys.glob("Piranha_peaks/*.dist")))))})))
# rerun_dist<-filter(as.data.frame(num_dist), Freq < 6)
#
# for(file in rerun_dist$Var1){
#   cat(paste0("qsub -cwd ", file, ".sh"),
#       file = "submit_to_JHPCE.sh",
#       sep = "\n",
#       append = T)
# }

ENCODE_info<-read.csv("../file_information.csv")
matched_names<-read.csv("match_CLIPper_Piranha_file_names.csv") %>%
  mutate(Piranha_name = paste0(protein, "_", cell, "_rep", Piranha_rep))
all_files<-c(as.character(matched_names$file)[c(TRUE, FALSE)], # actually only taking every other file
             as.character(matched_names$Accession)[c(TRUE, FALSE)],
             as.character(matched_names$Piranha_name)[c(TRUE, FALSE)])
pair_files<-expand.grid(all_files, all_files, stringsAsFactors = F) %>% mutate(distance = 0)

for(n in 1:nrow(pair_files)){
  if(pair_files[n, "Var1"] %in% matched_names$file){ #no peak calling
    dir_stereogene_output<-"../all_dist_files"
    m<-match(pair_files[n, "Var1"], matched_names$file)
    CapR_prefix<-paste0(matched_names[m, "cell"], "_3UTR")
    protein_file<-c(paste0(matched_names[m, "file"], "_liftOver"), paste0(matched_names[m+1, "file"], "_liftOver"))
    protein_file_input<-paste0(ENCODE_info[which(ENCODE_info$protein == as.character(matched_names[m, "protein"]) &
                                                 ENCODE_info$cell == as.character(matched_names[m, "cell"]) &
                                                 ENCODE_info$input == TRUE), "file"], "_liftOver")
  }else if(pair_files[n, "Var1"] %in% matched_names$Accession){ #CLIPper
    dir_stereogene_output<-"CLIPper_peaks"
    CapR_prefix<-paste0(matched_names[match(pair_files[n, "Var1"], matched_names$Accession), "cell"], "_3UTR")
    m<-match(pair_files[n, "Var1"], matched_names$Accession)
    protein_file<-c(paste0(matched_names[m, "Accession"], "_liftOver"), paste0(matched_names[m+1, "Accession"], "_liftOver"))
    protein_file_input<-NULL
  }else{ #Piranha
    dir_stereogene_output<-"Piranha_peaks"
    CapR_prefix<-paste0(matched_names[match(pair_files[n, "Var1"], matched_names$Piranha_name), "cell"], "_3UTR")
    m<-match(pair_files[n, "Var1"], matched_names$Piranha_name)
    protein_file<-c(paste0(matched_names[m, "Piranha_name"], "_liftOver"), paste0(matched_names[m+1, "Piranha_name"], "_liftOver"))
    protein_file_input<-NULL
  }
  if(pair_files[n, "Var2"] %in% matched_names$file){ #no peak calling
    dir_stereogene_output_2<-"../all_dist_files"
    CapR_prefix_2<-paste0(matched_names[match(pair_files[n, "Var2"], matched_names$file), "cell"], "_3UTR")
    l<-match(pair_files[n, "Var2"], matched_names$file)
    protein_file_2<-c(paste0(matched_names[l, "file"], "_liftOver"), paste0(matched_names[l+1, "file"], "_liftOver"))
    protein_file_input_2<-paste0(ENCODE_info[which(ENCODE_info$protein == as.character(matched_names[l, "protein"]) &
                                                   ENCODE_info$cell == as.character(matched_names[l, "cell"]) &
                                                   ENCODE_info$input == TRUE), "file"], "_liftOver")
  }else if(pair_files[n, "Var2"] %in% matched_names$Accession){ #CLIPper
    dir_stereogene_output_2<-"CLIPper_peaks"
    CapR_prefix_2<-paste0(matched_names[match(pair_files[n, "Var2"], matched_names$Accession), "cell"], "_3UTR")
    l<-match(pair_files[n, "Var2"], matched_names$Accession)
    protein_file_2<-c(paste0(matched_names[l, "Accession"], "_liftOver"), paste0(matched_names[l+1, "Accession"], "_liftOver"))
    protein_file_input_2<-NULL
  }else{ #Piranha
    dir_stereogene_output_2<-"Piranha_peaks"
    CapR_prefix_2<-paste0(matched_names[match(pair_files[n, "Var2"], matched_names$Piranha_name), "cell"], "_3UTR")
    l<-match(pair_files[n, "Var2"], matched_names$Piranha_name)
    protein_file_2<-c(paste0(matched_names[l, "Piranha_name"], "_liftOver"), paste0(matched_names[l+1, "Piranha_name"], "_liftOver"))
    protein_file_input_2<-NULL
    }
  distance<-bindingContextDistanceCapR(dir_stereogene_output = dir_stereogene_output,
                                            CapR_prefix = CapR_prefix,
                                            protein_file = protein_file,
                                            protein_file_input = protein_file_input,
                                            dir_stereogene_output_2 = dir_stereogene_output_2,
                                            CapR_prefix_2 = CapR_prefix_2,
                                            protein_file_2 = protein_file_2,
                                            protein_file_input_2 = protein_file_input_2)
  if(is.numeric(distance)){pair_files[n, "distance"]<-distance}
  if(n %% 100 == 0){
    write.csv(pair_files, "distances_compare_peaks.csv", row.names = F)
  }
}
write.csv(pair_files, "distances_compare_peaks.csv", row.names = F)

pair_files<-read.csv("distances_compare_peaks.csv")

pair_proteins <- filter(pair_files, distance > 0)
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

toMatch<-c("blue", "red", "gold", "orange", "magenta", "orchid", "green", "pink", "violet", "brown", "aquamarine", "cyan")
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

## pair technical replicates
expts<-matched_names[,c("protein", "cell")] %>% unique() %>%
  mutate(aligned = 0, Piranha = 0, CLIPper = 0)
for(i in 1:nrow(expts)){
  aligned<-matched_names[which(matched_names$protein == expts[i, "protein"] &
                                 matched_names$cell == expts[i, "cell"]), "file"]
  Piranha<-matched_names[which(matched_names$protein == expts[i, "protein"] &
                                 matched_names$cell == expts[i, "cell"]), "Piranha_name"]
  CLIPper<-matched_names[which(matched_names$protein == expts[i, "protein"] &
                                 matched_names$cell == expts[i, "cell"]), "Accession"]
  expts[i, "aligned"]<-bindingContextDistanceCapR(dir_stereogene_output = "../all_dist_files",
                                                  CapR_prefix = paste0(expts[i, "cell"], "_3UTR"),
                                                  protein_file = paste0(aligned[1], "_liftOver"),
                                                  CapR_prefix_2 = paste0(expts[i, "cell"], "_3UTR"),
                                                  protein_file_2 = paste0(aligned[2], "_liftOver"))
  expts[i, "Piranha"]<-bindingContextDistanceCapR(dir_stereogene_output = "Piranha_peaks",
                                                  CapR_prefix = paste0(expts[i, "cell"], "_3UTR"),
                                                  protein_file = paste0(Piranha[1], "_liftOver"),
                                                  CapR_prefix_2 = paste0(expts[i, "cell"], "_3UTR"),
                                                  protein_file_2 = paste0(Piranha[2], "_liftOver"))
  expts[i, "CLIPper"]<-bindingContextDistanceCapR(dir_stereogene_output = "CLIPper_peaks",
                                                  CapR_prefix = paste0(expts[i, "cell"], "_3UTR"),
                                                  protein_file = paste0(CLIPper[1], "_liftOver"),
                                                  CapR_prefix_2 = paste0(expts[i, "cell"], "_3UTR"),
                                                  protein_file_2 = paste0(CLIPper[2], "_liftOver"))
}

expts_plot<-gather(expts, "data", "distance", -cell, -protein)
ggplot(expts_plot, aes(x = data, y = distance)) +
  geom_boxplot() +
  theme_classic()
ggparcoord(expts, columns = 3:5, groupColumn = c(6),
           showPoints = T,
           alphaLines = .6,
           scale = "globalminmax") +
  theme_classic() +
  theme(legend.position="none",
        text = element_text(size=14)) +
  scale_color_manual(values= rep("black", 40)) +
  xlab(NULL) + ylab("Distance")

expts_mutate<-expts
#expts_mutate$Piranha<-expts_mutate$CLIPper - expts_mutate$Piranha
#expts_mutate$CLIPper<-expts_mutate$aligned - expts_mutate$CLIPper
#expts_mutate$aligned<-0
expts_mutate_plot<-gather(expts_mutate, "data", "distance", -cell, -protein)
ggplot(filter(expts_mutate_plot, data != "aligned"), aes(x = data, y = distance)) +
  geom_boxplot() +
  xlab(NULL) + 
  ylab("Aligned reads minus peak distance") +
  theme_classic() +
  theme(text = element_text(size=14))
ggsave("relative_distance.jpg")

## identify min distance
expts$min<-""
for(i in 1:nrow(expts)){
  if(expts[i, "aligned"] < expts[i, "Piranha"] &
     expts[i, "aligned"] < expts[i, "CLIPper"]){
    expts[i, "min"] <-"aligned"
  }
  else if(expts[i, "Piranha"] < expts[i, "aligned"] &
     expts[i, "Piranha"] < expts[i, "CLIPper"]){
    expts[i, "min"] <-"Piranha"
  }
  else if(expts[i, "CLIPper"] < expts[i, "aligned"] &
     expts[i, "CLIPper"] < expts[i, "Piranha"]){
    expts[i, "min"] <-"CLIPper"
  }
}
ggparcoord(expts, columns = 3:5, groupColumn = 7,
           showPoints = T,
           alphaLines = .6,
           scale = "globalminmax") +
  theme_classic() +
  scale_color_discrete(name = "Method with\nminimum distance") +
  theme(text = element_text(size=14)) +
  xlab(NULL) + ylab("Distance between replicates")
ggsave("parallel_coord_compare.jpg")
pie(c(sum(expts$min == "aligned")/40, 
      sum(expts$min == "CLIPper")/40, 
      sum(expts$min == "Piranha")/40),
    labels = c("aligned reads", "CLIPper peaks", "Piranha peaks"))






