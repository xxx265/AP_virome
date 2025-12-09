library(vegan)    
library(ggplot2)  
library(ggpubr)   
library(tidyverse) 

#APd
abundance <- read.table("../00_data/AP.tpm", check.names=F, row.names = 1, sep="\t", header=T)
metadata <- read.table("../APd.group", sep="\t", header=T, check.names=F)
abundance <- abundance[, metadata$Sample]
rownames(metadata) <- metadata$Sample
print(unique(metadata$Group)) 
abundance_t <- t(abundance)
#Bray-Curtis
dist_matrix <- vegdist(abundance_t, method="bray")
adonis_result <- adonis2(dist_matrix ~ Group, data=metadata, permutations=999)
print(adonis_result)

#pairwise
pairwise_adonis <- function(dist_matrix, group_var, permutations=999) {
  groups <- unique(group_var)
  result <- data.frame()
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      subset_idx <- group_var %in% c(groups[i], groups[j])
      sub_dist <- as.dist(as.matrix(dist_matrix)[subset_idx, subset_idx])
      res <- adonis2(sub_dist ~ group_var[subset_idx], permutations=permutations)
      result <- rbind(result, data.frame(
        Group1 = groups[i],
        Group2 = groups[j],
        R2 = res$R2[1],
        p.value = res$`Pr(>F)`[1]
      ))
    }
  }
  result$p.adj <- p.adjust(result$p.value, method="BH") 
  return(result)
}

pairwise_results <- pairwise_adonis(dist_matrix, metadata$Group)
print(pairwise_results)
#write.table(pairwise_results, "pairwise_results_etiology.tsv", sep="\t", row.names = FALSE,quote = FALSE)

groups <- unique(metadata$Group)
full_matrix <- expand.grid(Group1 = groups, Group2 = groups) %>%
  left_join(pairwise_results, by = c("Group1", "Group2")) %>%
  left_join(pairwise_results, by = c("Group1" = "Group2", "Group2" = "Group1")) %>%
  mutate(
    R2 = ifelse(is.na(R2.x), R2.y, R2.x),
    p.adj = ifelse(is.na(p.adj.x), p.adj.y, p.adj.x)
  ) %>%
  select(Group1, Group2, R2, p.adj) %>%
  mutate(R2 = ifelse(Group1 == Group2, NA, R2))
full_matrix$Group1 <- factor(full_matrix$Group1, levels = groups)
full_matrix$Group2 <- factor(full_matrix$Group2, levels = groups)

# plot
p <- ggplot(full_matrix, aes(x = Group1, y = Group2, fill = R2)) +
  geom_tile(color = "grey70", linewidth = 0.2) +
  scale_fill_gradient2(low = "white", high = "#A2C9A1", 
                       midpoint = 0, limit = c(0, max(full_matrix$R2, na.rm = TRUE)),
                       name = expression(R^2), na.value = "white") +
  geom_text(aes(label = ifelse(is.na(R2), "", 
                               ifelse(p.adj < 0.05, sprintf("%.4f*", R2), 
                                      sprintf("%.4f", R2)))), 
color = "black", size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text( vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank()) +
  coord_fixed() +
  ggtitle("Pairwise Effect Size (R²)") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(full_matrix$Group2)))
p
ggsave("adonis_etiology.pdf",p, width=4, height=4)


#APs
abundance <- read.table("../00_data/AP.tpm", check.names=F, row.names = 1, sep="\t", header=T)
metadata <- read.table("../APs.group", sep="\t", header=T, check.names=F)
abundance <- abundance[, metadata$Sample]
rownames(metadata) <- metadata$Sample
print(unique(metadata$Group)) 
abundance_t <- t(abundance)
#Bray-Curtis
dist_matrix <- vegdist(abundance_t, method="bray")
adonis_result <- adonis2(dist_matrix ~ Group, data=metadata, permutations=999)
print(adonis_result)

# pairwise
pairwise_adonis <- function(dist_matrix, group_var, permutations=999) {
  groups <- unique(group_var)
  result <- data.frame()
  for (i in 1:(length(groups)-1)) {
    for (j in (i+1):length(groups)) {
      subset_idx <- group_var %in% c(groups[i], groups[j])
      sub_dist <- as.dist(as.matrix(dist_matrix)[subset_idx, subset_idx])
      res <- adonis2(sub_dist ~ group_var[subset_idx], permutations=permutations)
      result <- rbind(result, data.frame(
        Group1 = groups[i],
        Group2 = groups[j],
        R2 = res$R2[1],
        p.value = res$`Pr(>F)`[1]
      ))
    }
  }
  result$p.adj <- p.adjust(result$p.value, method="BH")  
  return(result)
}


pairwise_results <- pairwise_adonis(dist_matrix, metadata$Group)
print(pairwise_results)
#write.table(pairwise_results, "pairwise_results_severity.tsv", sep="\t", row.names = FALSE,quote = FALSE)

groups <- unique(metadata$Group)
full_matrix <- expand.grid(Group1 = groups, Group2 = groups) %>%
  left_join(pairwise_results, by = c("Group1", "Group2")) %>%
  left_join(pairwise_results, by = c("Group1" = "Group2", "Group2" = "Group1")) %>%
  mutate(
    R2 = ifelse(is.na(R2.x), R2.y, R2.x),
    p.adj = ifelse(is.na(p.adj.x), p.adj.y, p.adj.x)
  ) %>%
  select(Group1, Group2, R2, p.adj) %>%
  mutate(R2 = ifelse(Group1 == Group2, NA, R2))

full_matrix$Group1 <- factor(full_matrix$Group1, levels = groups)
full_matrix$Group2 <- factor(full_matrix$Group2, levels = groups)

# plot
p <- ggplot(full_matrix, aes(x = Group1, y = Group2, fill = R2)) +
  geom_tile(color = "grey70", linewidth = 0.2) +
  scale_fill_gradient2(low = "white", high = "#A2C9A1", 
                       midpoint = 0, limit = c(0, max(full_matrix$R2, na.rm = TRUE)),
                       name = expression(R^2), na.value = "white") +
  geom_text(aes(label = ifelse(is.na(R2), "", 
                               ifelse(p.adj < 0.05, sprintf("%.4f*", R2), 
                                      sprintf("%.4f", R2)))), 
            color = "black", size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text( vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank()) +
  coord_fixed() +
  ggtitle("Pairwise Effect Size (R²)") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(levels(full_matrix$Group2)))

p
ggsave("adonis_severity.pdf",p, width=4, height=4)
