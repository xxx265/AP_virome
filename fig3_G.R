library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggforce)

#data prepare
#significant vOTUs
SAP_dt = read.table("16_new_new/SAP/rf_imp.sorted", sep="\t", header=T, row.names=1, check.names=F)
SAP <- SAP_dt[1:30, ]
MSAP_dt = read.table("16_new_new/MSAP/rf_imp.sorted", sep="\t", header=T, row.names=1, check.names=F)
MSAP <- MSAP_dt[1:30, ]
MAP_dt = read.table("16_new_new/MAP/rf_imp.sorted", sep="\t", header=T, row.names=1, check.names=F)
MAP <- MAP_dt[1:30, ]
file_list <- list(SAP, MSAP, MAP)
all_viruses <- Reduce(union, lapply(file_list, rownames))
all_viruses <- union(rownames(SAP), union(rownames(MSAP), rownames(MAP)))
pv = read.table("04_value/AP_sig.py_pvalue", sep="\t", header=T, check.names=F)
pv$name <- tolower(trimws(pv$name))
all_viruses <- tolower(trimws(all_viruses))
pv_filtered <- pv %>%
  filter(!is.na(name), name %in% all_viruses)
write.table(pv_filtered, "16_new_new/heatmap/filtered_pvalue.txt", sep="\t", row.names=FALSE)
dt = read.table("04_value/AP_sig.tpm", sep="\t", header=T, row.names=1, check.names=F)
MAP_columns <- which(sample_map$Group %in% c("HC", "MAP"))
MAP_map <- sample_map[MAP_columns,]
dt1 <- dt[rownames(dt) %in% all_viruses, ]
write.table(MAP_map, "15_new_heatmap/MAP.group", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(dt1, "16_new_new/heatmap/sel.tpm", sep = "\t", row.names = TRUE, quote = FALSE)

#APs
source("/share/data1/zhangy2/scripts/R_my_functions/zy_alpha_diversity.R")
sample_map = read.table("../APs.group", sep="\t", header=T, check.names=F)
sample_map1 = read.table("../APd.group", sep="\t", header=T, check.names=F)
sample_map$Group = factor(sample_map$Group, levels = group_order)
dt = read.table("16_new_new/heatmap/sel.tpm", sep="\t", header=T, row.names=1, check.names=F)

pvalue_data <- read.table("16_new_new/heatmap/filtered_pvalue.txt", sep="\t", header=TRUE)
filtered_virus_sorted <- pvalue_data[order(pvalue_data$enriched), ]
selected_votus <- filtered_virus_sorted$name 
#write.table(selected_votus, "selected_votus.tsv", row.names = FALSE)

control_group <- "HC"
disease_groups <- unique(sample_map$Group[sample_map$Group != control_group]) 

log2fc_matrix <- matrix(nrow = length(selected_votus), ncol = length(disease_groups),
                        dimnames = list(selected_votus, disease_groups))
pval_matrix <- matrix(nrow = length(selected_votus), ncol = length(disease_groups),
                      dimnames = list(selected_votus, disease_groups))

for (group in disease_groups) {
  samples_group <- sample_map %>% filter(Group == group) %>% pull(Sample)
  samples_hc <- sample_map %>% filter(Group == control_group) %>% pull(Sample)
  fc <- log2(
    (rowMeans(dt[selected_votus, samples_group, drop = FALSE] + 1e-6)) / 
      (rowMeans(dt[selected_votus, samples_hc, drop = FALSE] + 1e-6))
  )
  log2fc_matrix[, group] <- fc
  pvals <- sapply(selected_votus, function(votu) {
    wilcox.test(
      as.numeric(dt[votu, samples_group]),
      as.numeric(dt[votu, samples_hc])
    )$p.value
  })
  pval_matrix[, group] <- p.adjust(pvals, method = "BH")  # FDR校正
}

#APd
control_group <- "HC"
disease_groups <- unique(sample_map1$Group[sample_map1$Group != control_group])  

log2fc_matrix1 <- matrix(nrow = length(selected_votus), ncol = length(disease_groups),
                        dimnames = list(selected_votus, disease_groups))
pval_matrix1 <- matrix(nrow = length(selected_votus), ncol = length(disease_groups),
                      dimnames = list(selected_votus, disease_groups))

for (group in disease_groups) {
  samples_group <- sample_map1 %>% filter(Group == group) %>% pull(Sample)
  samples_hc <- sample_map1 %>% filter(Group == control_group) %>% pull(Sample)
  
  cat("Group:", group, 
      "\n  Samples in group:", length(samples_group),
      "\n  Samples in HC:", length(samples_hc), "\n")
  
  if (length(samples_group) == 0 | length(samples_hc) == 0) {
    cat("Skipping group", group, "due to zero samples\n")
    next
  }
  fc <- log2(
    (rowMeans(dt[selected_votus, samples_group, drop = FALSE] + 1e-6)) / 
      (rowMeans(dt[selected_votus, samples_hc, drop = FALSE] + 1e-6))
  )
  log2fc_matrix1[, group] <- fc
  
  pvals <- sapply(selected_votus, function(votu) {
    wilcox.test(
      as.numeric(dt[votu, samples_group]),
      as.numeric(dt[votu, samples_hc])
    )$p.value
  })
  pval_matrix1[, group] <- p.adjust(pvals, method = "BH")  
}

log2fc_matrix_f <- cbind(log2fc_matrix, log2fc_matrix1)
pval_matrix_f <- cbind(pval_matrix, pval_matrix1)

###plot data
transposed_log2fc <- t(log2fc_matrix_f)  
transposed_pval <- t(pval_matrix_f)
group_annotation <- data.frame(
  Group = c(rep("Control", 37), rep("Disease", ncol(transposed_log2fc) - 37))
)
rownames(group_annotation) <- colnames(transposed_log2fc)  

combined_data <- transposed_log2fc
trend_symbols <- matrix("", nrow = nrow(combined_data), ncol = ncol(combined_data))
rownames(trend_symbols) <- rownames(combined_data)
colnames(trend_symbols) <- colnames(combined_data)

for (i in 1:nrow(transposed_pval)) {
  for (j in 1:ncol(transposed_pval)) {
    p_value <- transposed_pval[i, j]
    
    # 根据 p 值填充符号
    if (p_value <= 0.001) {
      trend_symbols[i, j] <- "***"
    } else if (p_value <= 0.01) {
      trend_symbols[i, j] <- "**"
    } else if (p_value <= 0.05) {
      trend_symbols[i, j] <- "*"
    } else {
      trend_symbols[i, j] <- ""  
    }
  }
}
display_func <- function(mat) {
  return(trend_symbols)  
}

min_value <- min(combined_data, na.rm = TRUE)
max_value <- max(combined_data, na.rm = TRUE)
min_value_rounded <- floor(min_value)
max_value_rounded <- round(max_value) 


group_annotation <- data.frame(
  Group = c(rep("Control", 37), rep("Disease", ncol(transposed_log2fc) - 37))
)
rownames(group_annotation) <- colnames(transposed_log2fc)  
#family
host<- read.table("08_host/virus_sig_host", sep="\t", header=T, check.names=F)
dt_votus <- filtered_virus_sorted[, "name", drop = FALSE]
result <- merge(dt_votus, host, by = "name", all = FALSE)
result <- result[order(result$enriched), ]
family <- unique(result[, c("name", "family")])
group_annotation <- merge(family[, c("name", "family")], group_annotation,  by.x = "name", by.y = "row.names", all.x = TRUE)
rownames(group_annotation) <- group_annotation$name
group_annotation <- group_annotation[, -which(names(group_annotation) == "name")]

#type
type <- read.table("/share/data1/yangw/IMID_416/bacphlip_filter", sep="\t", header=F, check.names=F)
result1 <- merge(dt_votus, type, by.x = 1, by.y = 1, all.x = TRUE)
colnames(result1)[4] <- "type"
result1$type <- ifelse(is.na(result1$type), "unknown", result1$type)
group_annotation <- merge(result1[, c("name","type")], group_annotation, by.x = "name", by.y = "row.names", all.x = TRUE)
rownames(group_annotation) <- group_annotation$name
group_annotation <- group_annotation[, -which(names(group_annotation) == "name")]
group_annotation <- group_annotation[, c("Group", "family", "type")]

#host
host<- read.table("08_host/virus_sig_host_new", sep="\t", header=T, check.names=F)
dt_votus <- filtered_virus_sorted[, "name", drop = FALSE]
result <- merge(dt_votus, host, by = "name", all = FALSE)
result <- result[order(result$enriched), ]
filtered_result_c <- result[result$enriched == "Control", ]
unique_values_c <- unique(filtered_result_c[, c("name", "host_family")])
geuns_count_c <- as.data.frame(table(unique_values_c$host_family))
filtered_result_d <- result[result$enriched == "Disease", ]
unique_values_d <- unique(filtered_result_d[, c("name", "host_family")])
geuns_count_d <- as.data.frame(table(unique_values_d$host_family))
geuns_count_c <- geuns_count_c[order(-geuns_count_c$Freq), ]
geuns_count_d <- geuns_count_d[order(-geuns_count_d$Freq), ]
top_c <- head(geuns_count_c$Var1, 5)
c_host <- data.frame(matrix(ncol = 5, nrow = 46))
colnames(c_host) <- top_c
rownames(c_host) <- rownames(group_annotation)
for (i in 1:nrow(result)) {
  virus_name <- result$name[i]
  host_family <- result$host_family[i]
    if (virus_name %in% rownames(c_host) && host_family %in% colnames(c_host)) {
     c_host[virus_name, host_family] <- result$enriched[i]
  } else {
    c_host[virus_name, host_family] <- NA
  }
}
c_host[is.na(c_host)] <- "others"
top_c_levels <- as.character(top_c)
c_host <- c_host[, colnames(c_host) %in% top_c_levels]
top_d <- head(geuns_count_d$Var1, 3)
d_host <- data.frame(matrix(ncol = 3, nrow = 46))
colnames(d_host) <- top_d
rownames(d_host) <- rownames(group_annotation)
for (i in 1:nrow(result)) {
  virus_name <- result$name[i]
  host_family <- result$host_family[i]
  
  if (virus_name %in% rownames(d_host) && host_family %in% colnames(d_host)) {
    d_host[virus_name, host_family] <- result$enriched[i]
  } else {
    d_host[virus_name, host_family] <- NA
  }
}
d_host[is.na(d_host)] <- "others"
top_d_levels <- as.character(top_d)
d_host <- d_host[, colnames(d_host) %in% top_d_levels]
group_annotation <- merge(group_annotation, c_host, by = "row.names", all.x = TRUE)
rownames(group_annotation) <- group_annotation$Row.names
group_annotation <- group_annotation[, -which(names(group_annotation) == "Row.names")]
group_annotation <- merge(group_annotation, d_host, by = "row.names", all.x = TRUE)
rownames(group_annotation) <- group_annotation$Row.names
group_annotation <- group_annotation[, -which(names(group_annotation) == "Row.names")]

group_annotation <- group_annotation[, c("Group", "family", "type", "Lachnospiraceae", "Oscillospiraceae", "Erysipelotrichaceae", 
                                         "Ruminococcaceae", "Acutalibacteraceae", "Enterobacteriaceae", "Veillonellaceae", "no-host")]
#color
annotation_colors <- list(
  type = c("Virulent" = "#fff591",
           "Temperate" = "#f6bed6",
           "unknown" = "#ADD8E6"),
  family = c("Inoviridae" = "#C9EFBE",
             "Retroviridae" = "#F0CFEA",     
             "unclassified" = "#CAC8EF",    
             "Peduoviridae" = "#9BDCFC"),  
  Group = c("Control" = "#66C466",  # 绿色
            "Disease" = "#7F5FBF"), # 紫色
  Lachnospiraceae = c("Control" = "#1976D2",  
                      "Disease" = "#E65100",
                      "others" = "white"),
  Oscillospiraceae = c("Control" = "#2196F3",  
                       "Disease" = "#FF5722",
                       "others" = "white"),
  Erysipelotrichaceae = c("Control" = "#42A5F5",  
                      "Disease" = "#FF7043",
                      "others" = "white"),
  Ruminococcaceae = c("Control" = "#64B5F6",  
                          "Disease" = "#FF8A65",
                          "others" = "white"),
  Acutalibacteraceae = c("Control" = "#BBDEFB",  
                         "Disease" = "#FFCC80",
                         "others" = "white"),
  Enterobacteriaceae = c("Control" = "#2E7D32", 
                         "Disease" = "#C62828",
                         "others" = "white"),
  Veillonellaceae = c("Control" = "#388E3C",  
                "Disease" = "#D32F2F",
                "others" = "white"),
  `no-host` = c("Control" = "#4CAF50",  
                     "Disease" = "#F44336",
                     "others" = "white")
)

library(scales)
my_palette <- gradient_n_pal(
  colours = c("#4681B8", "white", "#EF4F2F"),
  values = c(-5, 0, 3)  # 直接映射数据值到颜色
)(seq(-5, 3, length.out = 200))

#colnames()
p <- pheatmap(
  mat = combined_data,
  color = my_palette,  
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_col = 37, 
  gaps_row = 3,
  main = "Differential Abundance of vOTUs (log2FC vs HC)",
  fontsize_row = 8,
  fontsize_col = 6,
  display_numbers = trend_symbols,  
  fontsize_number = 6,
  cellwidth = 8,
  cellheight = 15,
  angle_col = 45,
  legend = TRUE,  
  annotation_col = group_annotation,
  annotation_names_col = TRUE,
  annotation_colors = annotation_colors,  
  legend_breaks = c(min_value, max_value),  
  legend_labels = c(sprintf("%d", min_value_rounded), sprintf("%d", max_value_rounded))  
)
p

library(showtext)
showtext_auto() 

ggsave("fig3_G_heatmap.pdf", p, 
       width = 20, height = 25,
       device = cairo_pdf)






