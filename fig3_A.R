#of species
source("/share/data1/zhangy2/scripts/R_my_functions/zy_alpha_diversity.R")
sample_map = read.table("../APs.group", sep="\t", header=T, check.names=F)
sample_map1 = read.table("../APd.group", sep="\t", header=T, check.names=F)
dt = read.table("AP_sig.tpm", sep="\t", header=T, row.names=1, check.names=F)

#APs
species_counts <- colSums(dt > 0) 
df <- data.frame(Sample = names(species_counts), Species = species_counts) %>%
  left_join(sample_map, by = "Sample")
group_means <- df %>%
  group_by(Group) %>%
  summarise(Mean_Species = mean(Species)) %>%
  mutate(FC = log10(Mean_Species / Mean_Species[Group == "HC"]))

p_values <- list()
groups <- unique(df$Group)
for (group in groups[groups != "HC"]) {
  test <- wilcox.test(df$Species[df$Group == group], df$Species[df$Group == "HC"])
  p_values[[group]] <- test$p.value
}

plot_data <- group_means %>%
  filter(Group != "HC") %>%  
  mutate(P_value = sapply(Group, function(g) p_values[[g]])) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

#APd
df <- data.frame(Sample = names(species_counts), Species = species_counts) %>%
  left_join(sample_map1, by = "Sample")

group_means <- df %>%
  group_by(Group) %>%
  summarise(Mean_Species = mean(Species)) %>%
  mutate(FC = log10(Mean_Species / Mean_Species[Group == "HC"])) 

p_values <- list()
groups <- unique(df$Group)
for (group in groups[groups != "HC"]) {
  test <- wilcox.test(df$Species[df$Group == group], df$Species[df$Group == "HC"])
  p_values[[group]] <- test$p.value
}

plot_data1 <- group_means %>%
  filter(Group != "HC") %>%  
  mutate(P_value = sapply(Group, function(g) p_values[[g]])) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

result_plot = rbind(plot_data,plot_data1)
result_plot = rbind(plot_data,plot_data1)
result_plot$Group <- factor(
  result_plot$Group, 
  levels = rev(unique(result_plot$Group)) 
)
p1 <- ggplot(result_plot, aes(x = FC, y = Group)) +
  geom_bar(stat = "identity", fill = "gray", width = 0.7) +  
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) + 
  labs(
    x = "log10(Fold change) of species",
    y = "Comparison with HC",
    title = "# of species comparison"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.major.x = element_blank() 
  )

p1


#shanon
#APs
library(tidyverse)
library(vegan)
sample_map = read.table("../APs.group", sep="\t", header=T, check.names=F)
sample_map1 = read.table("../APd.group", sep="\t", header=T, check.names=F)
#sample_map$Group = factor(sample_map$Group, levels = group_order)
dt = read.table("AP_sig.tpm", sep="\t", header=T, row.names=1, check.names=F)

shannon_index <- vegan::diversity(t(dt),index="shannon") 
shannon_df <- data.frame(Sample = names(shannon_index), Shannon = shannon_index) %>%
  left_join(sample_map, by = "Sample")

group_means <- shannon_df %>%
  group_by(Group) %>%
  summarise(Mean_Shannon = mean(Shannon)) %>%
  mutate(FC = log10(Mean_Shannon / Mean_Shannon[Group == "HC"])) 

p_values <- list()
groups <- unique(shannon_df$Group)
for (group in groups[groups != "HC"]) {
  test <- wilcox.test(shannon_df$Shannon[shannon_df$Group == group], 
                      shannon_df$Shannon[shannon_df$Group == "HC"])
  p_values[[group]] <- test$p.value
}

plot_data <- group_means %>%
  filter(Group != "HC") %>%  
  mutate(P_value = sapply(Group, function(g) p_values[[g]])) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

#APd
shannon_index <- vegan::diversity(t(dt),index="shannon") 
shannon_df <- data.frame(Sample = names(shannon_index), Shannon = shannon_index) %>%
  left_join(sample_map1, by = "Sample")

group_means <- shannon_df %>%
  group_by(Group) %>%
  summarise(Mean_Shannon = mean(Shannon)) %>%
  mutate(FC = log10(Mean_Shannon / Mean_Shannon[Group == "HC"])) 

p_values <- list()
groups <- unique(shannon_df$Group)
for (group in groups[groups != "HC"]) {
  test <- wilcox.test(shannon_df$Shannon[shannon_df$Group == group], 
                      shannon_df$Shannon[shannon_df$Group == "HC"])
  p_values[[group]] <- test$p.value
}

plot_data1 <- group_means %>%
  filter(Group != "HC") %>% 
  mutate(P_value = sapply(Group, function(g) p_values[[g]])) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

result_plot = rbind(plot_data,plot_data1)
result_plot$Group <- factor(
  result_plot$Group, 
  levels = rev(unique(result_plot$Group))  
)

p2 <- ggplot(result_plot, aes(x = FC, y = Group)) +
  geom_col(fill = "gray", width = 0.7) +  
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  labs(
    x = "log10(Fold change) of Shannon",
    y = "Comparison with HC",
    title = "Shannon index comparison"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

p2

#shimpon
#APs
library(tidyverse)
library(vegan)
sample_map = read.table("../APs.group", sep="\t", header=T, check.names=F)
sample_map1 = read.table("../APd.group", sep="\t", header=T, check.names=F)
dt = read.table("AP_sig.tpm", sep="\t", header=T, row.names=1, check.names=F)

simpson_index <- vegan::diversity(t(dt), index = "simpson") 
simpson_df <- data.frame(Sample = names(simpson_index), Simpson = simpson_index) %>%
  left_join(sample_map, by = "Sample")

group_means <- simpson_df %>%
  group_by(Group) %>%
  summarise(Mean_Simpson = mean(Simpson)) %>%
  mutate(FC = log10(Mean_Simpson / Mean_Simpson[Group == "HC"])) 

p_values <- list()
groups <- unique(simpson_df$Group)
for (group in groups[groups != "HC"]) {
  test <- wilcox.test(simpson_df$Simpson[simpson_df$Group == group], 
                      simpson_df$Simpson[simpson_df$Group == "HC"])
  p_values[[group]] <- test$p.value
}

plot_data <- group_means %>%
  filter(Group != "HC") %>% 
  mutate(P_value = sapply(Group, function(g) p_values[[g]])) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ))

#APd
simpson_index <- vegan::diversity(t(dt), index = "simpson") 
simpson_df <- data.frame(Sample = names(simpson_index), Simpson = simpson_index) %>%
  left_join(sample_map1, by = "Sample")

group_means <- simpson_df %>%
  group_by(Group) %>%
  summarise(Mean_Simpson = mean(Simpson)) %>%
  mutate(FC = log10(Mean_Simpson / Mean_Simpson[Group == "HC"]))

p_values <- list()
groups <- unique(simpson_df$Group)
for (group in groups[groups != "HC"]) {
  test <- wilcox.test(simpson_df$Simpson[simpson_df$Group == group], 
                      simpson_df$Simpson[simpson_df$Group == "HC"])
  p_values[[group]] <- test$p.value
}

plot_data1 <- group_means %>%
  filter(Group != "HC") %>% 
  mutate(P_value = sapply(Group, function(g) p_values[[g]])) %>%
  mutate(Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ))
result_plot = rbind(plot_data,plot_data1)
result_plot$Group <- factor(
  result_plot$Group, 
  levels = rev(unique(result_plot$Group)) 
)

p3 <- ggplot(result_plot, aes(x = FC, y = Group)) +
  geom_col(fill = "gray", width = 0.7) + 
  geom_text(aes(label = Significance), vjust = -0.5, size = 5) +
  labs(
    x = "log10(Fold change) of simpson",
    y = "Comparison with HC",
    title = "simpson index comparison"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

p3

pa <- p1+p2+p3+
  plot_layout(ncol = 3)
pa
ggsave("fig3_A_diversity.pdf",pa, width=12,height=5)







