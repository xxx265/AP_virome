library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggforce)
library(rstatix)
library(ggpubr)

#data
#APs
sample_map2 = read.table("../AP.group", sep="\t", header=T, check.names=F)
sample_map = read.table("../APs.group", sep="\t", header=T, check.names=F)
sample_map1 = read.table("../APd.group", sep="\t", header=T, check.names=F)
dt = read.table("AP_sig.tpm", sep="\t", header=T, row.names=1, check.names=F)
pv = read.table("AP_sig.py_pvalue", sep="\t", header=T, check.names=F)

#APs
# data prepare ----------------------------------------------------------------
selected_votus_depleted <- pv %>%
  filter(enriched == "Control") %>%
  pull(name) 
selected_votus_enriched <- pv %>%
  filter(enriched == "Disease") %>%
  pull(name) 
enriched_abundance <- dt %>%
  filter(rownames(.) %in% selected_votus_enriched) %>%
  colSums() %>%
  enframe(name = "Sample", value = "Abundance") %>%
  mutate(Type = "AP-enriched")
depleted_abundance <- dt %>%
  filter(rownames(.) %in% selected_votus_depleted) %>%
  colSums() %>%
  enframe(name = "Sample", value = "Abundance") %>%
  mutate(Type = "AP-depleted")
combined_data <- bind_rows(enriched_abundance, depleted_abundance) %>%
  left_join(sample_map, by = "Sample") %>%
  mutate(Group = factor(Group, levels = c("HC","MAP", "MSAP", "SAP")))
combined_data <- combined_data %>%
  mutate(Group = factor(Group, levels = c("HC","MAP", "MSAP", "SAP")))
head(combined_data)
# stat.test ----------------------------------------------------------------
stat.test <- combined_data %>%
  group_by(Type) %>%
  wilcox_test(Abundance ~ Group, p.adjust.method = "BH") %>%
  add_xy_position(x = "Group", dodge = 0.8)
stat.test_filtered <- stat.test %>%
  filter(group1 == "HC")
# --------------------------------------------------------------
color_palette <- c("HC" = "#A2C9A1", "MAP" = "#8E7FB8", "MSAP" = "#DCE125", "SAP" = "#3FA0C0")

# plot --------------------------------------------------------------
p <- ggplot(combined_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1.5, show.legend = FALSE) +
  facet_wrap(~Type, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = color_palette) +
  labs(y = "Gross relative abundance (%)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(2, "lines")
  ) +
  stat_pvalue_manual(
    stat.test,
    label = "{p.adj.signif}",  
    tip.length = 0.005,                 
    bracket.size = 0.6,                  
    label.size = 4,                    
    step_increase = 0.15,               
    position = position_dodge(0.8),      
    inherit.aes = FALSE,                 
    coord.flip = FALSE                 
  )
p
stat.test_filtered$groups <- sapply(stat.test_filtered$groups, function(x) paste(x, collapse=", "))
ggsave("fig3_E_severity.pdf",p, width=5,height=5)

#APd
# data prepare ----------------------------------------------------------------
selected_votus_depleted <- pv %>%
  filter(enriched == "Control") %>%
  pull(name) 
selected_votus_enriched <- pv %>%
  filter(enriched == "Disease") %>%
  pull(name) 
enriched_abundance <- dt %>%
  filter(rownames(.) %in% selected_votus_enriched) %>%
  colSums() %>%
  enframe(name = "Sample", value = "Abundance") %>%
  mutate(Type = "AP-enriched")
depleted_abundance <- dt %>%
  filter(rownames(.) %in% selected_votus_depleted) %>%
  colSums() %>%
  enframe(name = "Sample", value = "Abundance") %>%
  mutate(Type = "AP-depleted")
combined_data <- bind_rows(enriched_abundance, depleted_abundance) %>%
  left_join(sample_map1, by = "Sample") %>%
  mutate(Group = factor(Group, levels = c("HC","ABP", "AHP", "APN", "Other")))
combined_data <- combined_data %>%
  mutate(Group = factor(Group, levels = c("HC","ABP", "AHP", "APN", "Other")))
head(combined_data)
# stat.test ----------------------------------------------------------------
stat.test <- combined_data %>%
  group_by(Type) %>%
  wilcox_test(Abundance ~ Group, p.adjust.method = "BH") %>%
  add_xy_position(x = "Group", dodge = 0.8)

#--------------------------------------------------------------
color_palette <- c("ABP" = "#A2C9A1", "AHP" = "#8E7FB8", "APN" = "#DCE125", "HC" = "#3FA0C0", "Other" = "#F28080")

# plot --------------------------------------------------------------
p <- ggplot(combined_data, aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1.5, show.legend = FALSE) +
  facet_wrap(~Type, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = color_palette) +
  labs(y = "Gross relative abundance (%)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing = unit(2, "lines")
  ) +
  stat_pvalue_manual(
    stat.test,
    label = "{p.adj.signif}",  
    tip.length = 0.005,                  
    bracket.size = 0.6,                 
    label.size = 4,                     
    step_increase = 0.15,               
    position = position_dodge(0.8),       
    inherit.aes = FALSE,                 
    coord.flip = FALSE                    
  )
p
stat.test_filtered$groups <- sapply(stat.test_filtered$groups, function(x) paste(x, collapse=", "))
ggsave("fig3_E_etiology.pdf",p, width=5,height=5)




























