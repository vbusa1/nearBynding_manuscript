library(nearBynding)

file_info<-read.csv("file_information.csv", stringsAsFactors = F)

HepG2_samples<-filter(file_info, cell == "HepG2")
K562_samples<-filter(file_info, cell == "K562")

complete_HepG2<-HepG2_samples %>% dplyr::count(protein) %>%
    filter(n == 3) %>%
    dplyr::select("protein") %>%
    unique() %>% unlist()

complete_K562<-K562_samples %>% dplyr::count(protein) %>%
    filter(n == 3) %>%
    dplyr::select("protein") %>%
    unique() %>% unlist()

## HepG2 cells

pair_proteins<-expand.grid(complete_HepG2, complete_HepG2, stringsAsFactors = F)
pair_proteins<-pair_proteins[apply(pair_proteins, 1, 
                                   function(x){x[1] != x[2]}),] %>% 
    mutate(distance = 0)
CapR_prefix<-CapR_prefix_2<-"HepG2_3UTR"
dir_stereogene_output<-dir_stereogene_output_2<-"all_dist_files"
for(n in 1:nrow(pair_proteins)){
    protein_file<-paste0(HepG2_samples[which(HepG2_samples$protein == pair_proteins_stem[n,1] &
                                                 HepG2_samples$input == F), "file"], "_liftOver")
    protein_file_input<-paste0(HepG2_samples[which(HepG2_samples$protein == pair_proteins_stem[n,1] &
                                                       HepG2_samples$input == T), "file"], "_liftOver")
    protein_file_2<-paste0(HepG2_samples[which(HepG2_samples$protein == pair_proteins_stem[n,2] &
                                                   HepG2_samples$input == F), "file"], "_liftOver")
    protein_file_input_2<-paste0(HepG2_samples[which(HepG2_samples$protein == pair_proteins_stem[n,2] &
                                                         HepG2_samples$input == T), "file"], "_liftOver")
    distance<-bindingContextDistanceCapR(dir_stereogene_output = dir_stereogene_output,
                                          CapR_prefix = CapR_prefix,
                                          protein_file = protein_file,
                                          protein_file_input = protein_file_input,
                                          dir_stereogene_output_2 = dir_stereogene_output_2,
                                          CapR_prefix_2 = CapR_prefix_2,
                                          protein_file_2 = protein_file_2,
                                          protein_file_input_2 = protein_file_input_2)
    pair_proteins[n, "distance"]<-distance
}
pair_proteins <- filter(pair_proteins, distance > 0)
saveRDS(pair_proteins, "wasserstein_distances_all_HepG2.rds")

## K562 cells

pair_proteins<-expand.grid(complete_K562, complete_K562, stringsAsFactors = F)
pair_proteins<-pair_proteins[apply(pair_proteins, 1, 
                                   function(x){x[1] != x[2]}),] %>% 
    mutate(distance = 0)
CapR_prefix<-CapR_prefix_2<-"K562_3UTR"
dir_stereogene_output<-dir_stereogene_output_2<-"all_dist_files"
for(n in 1:nrow(pair_proteins)){
    protein_file<-paste0(K562_samples[which(K562_samples$protein == pair_proteins_stem[n,1] &
                                                 K562_samples$input == F), "file"], "_liftOver")
    protein_file_input<-paste0(K562_samples[which(K562_samples$protein == pair_proteins_stem[n,1] &
                                                       K562_samples$input == T), "file"], "_liftOver")
    protein_file_2<-paste0(K562_samples[which(K562_samples$protein == pair_proteins_stem[n,2] &
                                                   K562_samples$input == F), "file"], "_liftOver")
    protein_file_input_2<-paste0(K562_samples[which(K562_samples$protein == pair_proteins_stem[n,2] &
                                                         K562_samples$input == T), "file"], "_liftOver")
    distance<-bindingContextDistanceCapR(dir_stereogene_output = dir_stereogene_output,
                                         CapR_prefix = CapR_prefix,
                                         protein_file = protein_file,
                                         protein_file_input = protein_file_input,
                                         dir_stereogene_output_2 = dir_stereogene_output_2,
                                         CapR_prefix_2 = CapR_prefix_2,
                                         protein_file_2 = protein_file_2,
                                         protein_file_input_2 = protein_file_input_2)
    pair_proteins[n, "distance"]<-distance
}
pair_proteins <- filter(pair_proteins, distance > 0)
saveRDS(pair_proteins, "wasserstein_distances_all_K562.rds")

## compare cell lines
complete_list<-c(paste0(complete_HepG2[complete_HepG2 %in% complete_K562], "_H"),
                 paste0(complete_K562[complete_K562 %in% complete_HepG2], "_K"))
pair_proteins<-expand.grid(complete_list, complete_list, stringsAsFactors = F)
pair_proteins<-pair_proteins[apply(pair_proteins, 1, function(x){x[1] != x[2]}),] %>%
    mutate(distance = 0)
dir_stereogene_output<-dir_stereogene_output_2<-"all_dist_files"
for(n in 1:nrow(pair_proteins)){
    if(substr(pair_proteins[n,1], nchar(pair_proteins[n,1]), nchar(pair_proteins[n,1])) == "H"){
        protein_file<-paste0(HepG2_samples[which(HepG2_samples$protein == substr(pair_proteins[n,1], 1, nchar(pair_proteins[n,1])-2) &
                                                     HepG2_samples$input == F), "file"], "_liftOver")
        protein_file_input<-paste0(HepG2_samples[which(HepG2_samples$protein == substr(pair_proteins[n,1], 1, nchar(pair_proteins[n,1])-2) &
                                                           HepG2_samples$input == T), "file"], "_liftOver")
        CapR_prefix<-"HepG2_3UTR"
    } else{
        protein_file<-paste0(K562_samples[which(K562_samples$protein == substr(pair_proteins[n,1], 1, nchar(pair_proteins[n,1])-2) &
                                                    K562_samples$input == F), "file"], "_liftOver")
        protein_file_input<-paste0(K562_samples[which(K562_samples$protein == substr(pair_proteins[n,1], 1, nchar(pair_proteins[n,1])-2) &
                                                          K562_samples$input == T), "file"], "_liftOver")
        CapR_prefix<-"K562_3UTR"
    }
    if(substr(pair_proteins[n,2], nchar(pair_proteins[n,2]), nchar(pair_proteins[n,2])) == "H"){
        protein_file_2<-paste0(HepG2_samples[which(HepG2_samples$protein == substr(pair_proteins[n,2], 1, nchar(pair_proteins[n,2])-2) &
                                                       HepG2_samples$input == F), "file"], "_liftOver")
        protein_file_input_2<-paste0(HepG2_samples[which(HepG2_samples$protein == substr(pair_proteins[n,2], 1, nchar(pair_proteins[n,2])-2) &
                                                             HepG2_samples$input == T), "file"], "_liftOver")
        CapR_prefix_2<-"HepG2_3UTR"
    } else{
        protein_file_2<-paste0(K562_samples[which(K562_samples$protein == substr(pair_proteins[n,2], 1, nchar(pair_proteins[n,2])-2) &
                                                      K562_samples$input == F), "file"], "_liftOver")
        protein_file_input_2<-paste0(K562_samples[which(K562_samples$protein == substr(pair_proteins[n,2], 1, nchar(pair_proteins[n,2])-2) &
                                                            K562_samples$input == T), "file"], "_liftOver")
        CapR_prefix_2<-"K562_3UTR"
    }
    distance<-try(DistributionDistanceWasserstein(dir_stereogene_output = dir_stereogene_output,
                                                  CapR_prefix = CapR_prefix,
                                                  protein_file = protein_file,
                                                  protein_file_input = protein_file_input,
                                                  dir_stereogene_output_2 = dir_stereogene_output_2,
                                                  CapR_prefix_2 = CapR_prefix_2,
                                                  protein_file_2 = protein_file_2,
                                                  protein_file_input_2 = protein_file_input_2),
                  silent = T)
    if(is.numeric(distance)){pair_proteins[n, "distance"]<-distance}
}
pair_proteins <- filter(pair_proteins, distance > 0)
for(protein in as.character(complete_HepG2[complete_HepG2 %in% complete_K562])){
    if(length(pair_proteins[which(pair_proteins$Var1 == paste0(protein, "_H") &
                                  pair_proteins$Var2 == paste0(protein, "_K")),3])>0){next}
    else{
        pair_proteins<-filter(pair_proteins,
                              Var1 != paste0(protein, "_K"),
                              Var1 != paste0(protein, "_H"),
                              Var2 != paste0(protein, "_K"),
                              Var2 != paste0(protein, "_H"))
    }
}
saveRDS(pair_proteins, "wasserstein_dist_both_cells.rds")


