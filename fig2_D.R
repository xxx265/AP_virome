setwd("/share/data3/zhangyt/AP/01_statistical")
rm(list = ls())
library(ggplot2)
library(dplyr)
#数据处理
#得出dis_sig.py_pvalue


##-----host select-------
dt <- read.table("AP_virus.py_pvalue", sep="\t", header=TRUE, row.names=1)
dt1 = read.table("AP_rf_imp.sorted", header=T, check.names=F, row.names=1)
matching_rows <- dt[dt$name %in% row.names(dt1), ]

dis <- matching_rows[match(row.names(dt1), matching_rows$name), ]
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dis$qvalue = p.adjust(dis$pvalue, method="BH")
dis$log2FC = ifelse(dis$enriched == "Control",log2(dt$fold_change), -log2(dt$fold_change))
dis$log2FC[is.infinite(dis$log2FC) & dt$log2FC == Inf] <- 10
dis$log2FC[is.infinite(dis$log2FC) & dt$log2FC == -Inf] <- -10
dis <- dis %>%  mutate(color = ifelse((qvalue<0.05 & fold_change > 1.2),enriched, "na"))
dis_sig.py_pvalue <- dis[, c("name", "enriched", "fold_change", "qvalue", "log2FC", "color")]

write.table(dis_sig.py_pvalue, "04_value/dis_sig.py_pvalue", sep="\t", row.names = FALSE)

##-----host plot---------
dis_sig <- read.table("08_host/virus_sig", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
family <- read.table("/share/data1/zhangy2/00.pub_data/10.China_vir_catalog_20241125/vOTU.family", sep = "\t", stringsAsFactors = FALSE)

colnames(family) <- c("name", "V2", "V3", "family")
family <- tibble::as_tibble(family)

dis_sig <- dis_sig %>%
  left_join(dplyr::select(family, name, family_column = family), by = "name") %>%
  mutate(family = ifelse(is.na(family_column), "unclassified", family_column)) %>%
  dplyr::select(-family_column) 

virus_host <- read.table("08_host/virus_host", header = FALSE) 
virus_host <- virus_host %>%
  rename(
    name = V1,        
    host = V2,    
    genus = V3    
  )

colnames(virus_host) <- c("name", "host", "genus")
dis_sig_host <- dis_sig %>%
  left_join(virus_host, by = "name") %>%
  mutate(
    host = ifelse(is.na(host), "unknown", host),
    genus = ifelse(is.na(genus), "unknown", genus)
  ) %>%
  dplyr::select(name, enriched, family, host, genus)  %>%
  arrange(name)  

df_summary <- dis_sig_host %>%
  group_by(family, genus, enriched) %>%
  summarise(virus_count = n(), .groups = "drop")

top_60_genus <- df_summary %>%
  group_by(genus) %>%
  summarise(total_virus_count = sum(virus_count)) %>%
  arrange(desc(total_virus_count)) %>%
  slice_head(n = 60) %>%
  pull(genus)

df_summary_top60 <- df_summary %>%
  filter(genus %in% top_60_genus)

all_families <- unique(df_summary_top60$family)
all_genus <- unique(df_summary_top60$genus)

p <- ggplot(df_summary_top60, aes(x = family, y = genus, size = virus_count, color = family)) +
  geom_point(alpha = 0.7) +  
  facet_wrap(~enriched) + 
  #scale_size(trans = "log10") +  
  scale_size_area(max_size = 5, trans = "sqrt") + 
  scale_x_discrete(limits = all_families) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "Virus-Host Enrichment Dot Plot",
    x = "Virus Family",
    y = "Host Genus",
    size = "Virus Count",
    color = "Virus Family"
  )
p

ggsave("fig2_D_host_botplot.pdf", p, width = 6, height = 8)

#-----Supplementary Table S4-------
write.table(dis_sig_host, file = "Supplementary_Table_S3.tsv", sep = "\t", row.names = FALSE)







