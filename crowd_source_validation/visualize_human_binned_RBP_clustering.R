library(tidyverse)

data<-read.csv("Bin RBP binding (Responses) - Form Responses.csv", header = F)[2:26] %>%
  t() %>% as.data.frame()
data$protein<-c("same", "different", "same", "same", "different", "different",
                "same", "different", "same", "different", "different", "same", 
                "same", "different", "same", "different", "same", "same", "different",
                "same", "same", "different", "same", "different", "same")
data$rank<-c("different", "same", "same", "different", "different", "same",
             "same", "different", "different", "same", "different", "different",
             "same", "different", "same", "same", "different", "same", "different",
             "same", "different", "same", "same", "different", "same")
data_gather<-gather(data, "bin", "count",-V1, -protein, -rank) %>% select(-bin)

protein_change<-c(`same` = "same protein",
           `different` = "different proteins")
rank_change<-c(`same` = "closely clustered",
        `different` = "distantly clustered")

ggplot(data_gather, aes(x = as.factor(V1), fill = count)) +
  geom_bar() +
  facet_grid(~ protein + rank, scale = "free",
             labeller = labeller(protein = as_labeller(protein_change),
                                 rank = as_labeller(rank_change))) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(legend.title=element_blank(),
        legend.text=element_text(size = 14),
        text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  scale_fill_manual(values = c("red", "blue"))
ggsave("human_binned_RBP_clustering.pdf", width =9, height =4)
