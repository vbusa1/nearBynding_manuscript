library(ggplot2)

pair_proteins_HepG2<-readRDS("wasserstein_distances_all_HepG2.rds")
pair_proteins_K562<-readRDS("wasserstein_distances_all_K562.rds")
file_info<-read.csv("file_information.csv", stringsAsFactors = F)

pair_proteins_HepG2$protein1<-lapply(pair_proteins_HepG2$Var1, function(x){
    file_info[which(file_info$file == x), "protein"]
}) %>% unlist()
pair_proteins_HepG2$protein2<-lapply(pair_proteins_HepG2$Var2, function(x){
    file_info[which(file_info$file == x), "protein"]
}) %>% unlist()
pair_rank_HepG2<-pair_proteins_HepG2 %>% group_by(Var1) %>%
    mutate(rank = order(order(distance)))
pair_same_protein_HepG2<-filter(pair_rank_HepG2, protein1 == protein2)
ggplot(pair_same_protein_HepG2, aes(x = rank)) +
    geom_histogram(bins = 150) +
    theme_classic() +
    labs(y = NULL, x = NULL) +
    theme(axis.text = element_text(size = 12)) +
    ggtitle("HepG2 rank")
ggsave("visual_hist_HepG2_samples_rank.pdf")

sum(pair_same_protein_HepG2$rank == 1)

pair_proteins_K562$protein1<-lapply(pair_proteins_K562$Var1, function(x){
    file_info[which(file_info$file == x), "protein"]
}) %>% unlist()
pair_proteins_K562$protein2<-lapply(pair_proteins_K562$Var2, function(x){
    file_info[which(file_info$file == x), "protein"]
}) %>% unlist()
pair_rank_K562<-pair_proteins_K562 %>% group_by(Var1) %>%
    mutate(rank = order(order(distance)))
pair_same_protein_K562<-filter(pair_rank_K562, protein1 == protein2)
ggplot(pair_same_protein_K562, aes(x = rank)) +
    geom_histogram(bins = 120) +
    theme_classic() +
    labs(y = NULL, x = NULL) +
    theme(axis.text = element_text(size = 12)) +
    ggtitle("K562 rank")
ggsave("visual_hist_K562_samples_rank.pdf")

sum(pair_same_protein_K562$rank == 1)

## compare across cell lines

pair_proteins<-readRDS("wasserstein_dist_both_cells.rds")

pair_proteins_dist<- pair_proteins %>% separate(Var1, c("protein_HepG2", "cell1"), "_") %>%
    separate(Var2, c("protein_K562", "cell2"), "_") %>%
    filter(cell1 == "H" & cell2 == "K") %>% dplyr::select(-cell1, -cell2)

pair_proteins_rank<-pair_proteins_dist %>% group_by(protein_K562) %>%
    mutate(rank = order(order(distance)))

top_rank<-filter(pair_proteins_rank, rank == 1)
sum(top_rank$protein_HepG2 == top_rank$protein_K562) # 15 of 73
sum(top_rank$protein_HepG2 == top_rank$protein_K562)/nrow(top_rank) #.2055

topish_rank<-filter(pair_proteins_rank, rank == 1 | rank == 2 | rank == 3)
sum(topish_rank$protein_HepG2 == topish_rank$protein_K562) # 21 of 73
sum(topish_rank$protein_HepG2 == topish_rank$protein_K562)/nrow(top_rank) # .2877

pair_same_protein<-filter(pair_proteins_rank, protein_HepG2 == protein_K562)
ggplot(pair_same_protein, aes(x = rank)) +
    geom_histogram(bins = 65) +
    theme_classic() +
    labs(y = NULL, x = NULL) +
    theme(axis.text = element_text(size = 12)) +
    ggtitle("across cells rank")
ggsave("visual_hist_cell_lines_rank.pdf")

