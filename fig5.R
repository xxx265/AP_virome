load("sample_color.RData")
load("group_order.RData")
library(dplyr)

#-----fig 5 A_ROC-----
source("/share/data1/zhangy2/scripts/R_my_functions/zy_plot_ROC.R")

dt_virus <- read.table("virus/virus_rf_re", skip = 1, sep = "\t", header = F)
colnames(dt_virus) <- c("Sample", "split", "Control", "Disease", "Actual", "nspecies", "seed")
dt_virus$nspecies <- 1049
dt_bac <- read.table("bac/bac_rf_re", skip = 1, sep = "\t", header = F)
colnames(dt_bac) <- c("Sample", "split", "Control", "Disease", "Actual", "nspecies", "seed")
dt_bac$nspecies <- 68
dt_merge <- read.table("merge/merge_rf_re", skip = 1, sep = "\t", header = F)
colnames(dt_merge) <- c("Sample", "split", "Control", "Disease", "Actual", "nspecies", "seed")
dt_merge$nspecies <- 68

x_virus <- plot_roc(dt_virus, pred="Disease", true="Actual", cols = "#ff7314")
x_bac <- plot_roc(dt_bac, pred="Disease", true="Actual",cols = "#56B4E9")
x_merge <- plot_roc(dt_merge, pred="Disease", true="Actual",cols = "#64cd64")

roc_virus <- x_virus$ROC[[1]]
roc_bac <- x_bac$ROC[[1]]
roc_merge <- x_merge$ROC[[1]]

labels_virus <- x_virus$labels
labels_bac <- x_bac$labels
labels_merge <- x_merge$labels

# plot
p_combined <- ggroc(list(virus = roc_virus, bac = roc_bac, merge = roc_merge)) +
  scale_color_manual(values = c("virus" = "#ff7314", "bac" = "#56B4E9", "merge" = "#64cd64")) +
  theme_bw() +
  geom_segment(data = data.frame(x = 0, y = 1),
               aes(x = x, y = y, xend = 1, yend = 0),
               color = "#d9d9d9", lwd = .4, inherit.aes = FALSE) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid = element_line(linetype = "dashed", color = "black", linewidth = 0.2),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black")
  )

p_combined <- p_combined + 
  annotate("segment", x = 0.5, xend = 0.55, y = 0.2, yend = 0.2, color = "#ff7314") +
  annotate("text", x = 0.49, y = 0.19, label = "virus: ", hjust = 0, vjust = 0) +
  annotate("text", x = 0.39, y = 0.19, label = labels_virus, hjust = 0, vjust = 0) +
  annotate("segment", x = 0.5, xend = 0.55, y = 0.15, yend = 0.15, color = "#56B4E9") +
  annotate("text", x = 0.49, y = 0.14, label = "bac: ", hjust = 0, vjust = 0) +
  annotate("text", x = 0.39, y = 0.14,label = labels_bac, hjust = 0, vjust = 0) +
  annotate("segment", x = 0.55, xend = 0.60, y = 0.1, yend = 0.1, color = "#64cd64") +
  annotate("text", x = 0.54, y = 0.09, label = "virus+bac: ", hjust = 0, vjust = 0) +
  annotate("text", x = 0.39, y = 0.09, label = labels_merge, hjust = 0, vjust = 0)
print(p_combined)
ggsave("ROC.pdf",p_combined, width=8,height = 5)

#-----fig 5 BCD_importance-------
#virus
dt_virus = read.table("virus/virus_rf_imp", sep="\t", header=T, row.names=1)
colnames(dt_virus)[1] <- "name"
dt_virus$y = rownames(dt_virus) 
map_virus = read.table("virus/AP_virus_sig.py_pvalue", sep="\t", header=T)
dtm_virus = dt_virus %>% merge(., map_virus,  by.y="name", all.x = T)
#plot
p <- dtm_virus %>%
  mutate(y = factor(y, levels = rev(y))) %>% 
  arrange(desc(Mean_Importance)) %>%  
  head(20) %>%  
  ggplot(aes(y = reorder(y, Mean_Importance), x = Mean_Importance, fill = enriched)) +  
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = sample_color) +  
  ggtitle("Virus Importance")
p
ggsave("virus/virus_sig_imp_plot.pdf", plot = p, width = 4.5, height = 6, dpi = 300)

#bac
dt_bac = read.table("bac/bac_rf_imp", sep="\t", header=T, row.names=1)
colnames(dt_bac)[1] <- "name"
dt_bac$y = rownames(dt_bac)  
map_bac = read.table("bac/AP_bac_sig.py_pvalue", sep="\t", header=T)
dtm_bac = dt_bac %>% merge(., map_bac,  by.y="name", all.x = T)
#plot
p <- dtm_bac %>%
  mutate(y = factor(y, levels = rev(y))) %>% 
  arrange(desc(Mean_Importance)) %>%  
  head(20) %>%  
  ggplot(aes(y = reorder(y, Mean_Importance), x = Mean_Importance, fill = enriched)) +  
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = sample_color) + 
  ggtitle("Bacteria Importance")  
p
ggsave("bac/bac_sig_imp_plot.pdf", plot = p, width = 6, height = 6, dpi = 300)

#merge
dt_merge = read.table("merge/merge_rf_imp", sep="\t", header=T, row.names=1)
colnames(dt_merge)[1] <- "name"
dt_merge$y = rownames(dt_merge)  
map_merge = read.table("merge/AP_merge_sig.py_pvalue", sep="\t", header=T)
dtm_merge = dt_merge %>% merge(., map_merge,  by.y="name", all.x = T)

p <- dtm_merge %>%
  mutate(y = factor(y, levels = rev(y))) %>% 
  arrange(desc(Mean_Importance)) %>%  
  head(20) %>%  
  ggplot(aes(y = reorder(y, Mean_Importance), x = Mean_Importance, fill = enriched)) +  
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = sample_color) +  
  ggtitle("Virus+Bacteria Importance")  
p
ggsave("merge/merge_sig_imp_plot_new.pdf", plot = p, width = 4.5, height = 6, dpi = 300)

#-----------fig 5 E_ROC-----------------
source('/share/data1/zhangy2/scripts/R_my_functions/zy_plot_ROC.R')
source('/share/data1/zhangy2/scripts/R_my_functions/zy_randomForest.R')
get_plot_data <- function(in_f, pred){
  dt = read.table(in_f, sep="\t", header=T)
  dt$new_group = paste(dt$seed, dt$nspecies, sep="|")
  data = calc_auc(dt, pred=pred,true="Actual", group="new_group")$table
  data$name  = rownames(data)
  data <- data %>%
    extract(new_group,c("seed","nspecies"),"(.*)\\|(.*)")
  data$nspecies = as.numeric(data$nspecies)
  plot_data <- data %>%
    group_by(nspecies) %>%
    summarise(avg.auc = mean(auc),
              sd.auc = sd(auc))
  return(plot_data)
}

get_plot = function(in_bac, in_virus, in_merge, pred, title){
  
  plot_data_bac <- get_plot_data(in_bac, pred)
  plot_data_virus <- get_plot_data(in_virus, pred)
  plot_data_merge <- get_plot_data(in_merge, pred)
  
  all_data <- rbind(
    data.frame(nspecies = plot_data_bac$nspecies, avg.auc = plot_data_bac$avg.auc, sd.auc = plot_data_bac$sd.auc, group = 'Bacteria'),
    data.frame(nspecies = plot_data_virus$nspecies, avg.auc = plot_data_virus$avg.auc, sd.auc = plot_data_virus$sd.auc, group = 'Virus'),
    data.frame(nspecies = plot_data_merge$nspecies, avg.auc = plot_data_merge$avg.auc, sd.auc = plot_data_merge$sd.auc, group = 'Merged')
  )
  max_auc_data <- all_data %>%
    group_by(group) %>%
    dplyr::filter(avg.auc == max(avg.auc)) %>%
    dplyr::select(nspecies, avg.auc, group)
  
  ggplot(all_data, aes(x = nspecies, y = avg.auc, color = group, group = group)) +
    geom_line(aes(color = group), size = 1) +
    geom_errorbar(aes(ymax = avg.auc + sd.auc, ymin = avg.auc - sd.auc, color = group), alpha = 0.5) +
    geom_point(data = max_auc_data, aes(x = nspecies, y = avg.auc), color = "black", size = 2) +
    geom_segment(data = max_auc_data, aes(x = nspecies, xend = nspecies, y = 75, yend = avg.auc, color = group), linetype = "dashed") +
    geom_segment(data = max_auc_data, aes(x = 0, xend = nspecies, y = avg.auc, yend = avg.auc, color = group), linetype = "dashed") +
    scale_color_manual(values = c("Virus" = "#ff7314", "Bacteria" = "#56B4E9", "Merged" = "#64cd64")) +
    labs(title = title, x = "nspecies", y = "avg_auc") +
    theme_minimal() +
    coord_cartesian(ylim = c(75, NA)) +
    theme(legend.position = "right") +
    annotate("segment", x = 0, xend = Inf, y = 75, yend = 75, color = "black", linetype = "solid") +
    annotate("segment", x = 0, xend = 0, y = 75, yend = Inf, color = "black", linetype = "solid") +
    coord_cartesian(xlim = c(0, max(all_data$nspecies) + 10), ylim = c(75, 100))+
    geom_text(data = max_auc_data, aes(x = nspecies, y = 90, label = paste(round(nspecies, 2)), color = group), 
              vjust = 27, hjust = 0.5, size = 3, fontface = "bold")+
    geom_text(data = max_auc_data[max_auc_data$group == "Virus", ], aes(x = 10, y = avg.auc - 0.5, label = paste(round(avg.auc, 2)), color = group), 
              vjust = -0.5, hjust = 1.5, size = 3, fontface = "bold") +
    geom_text(data = max_auc_data[max_auc_data$group == "Bacteria", ], aes(x = 10, y = avg.auc - 0.7, label = paste(round(avg.auc, 2)), color = group), 
              vjust = -1, hjust = 1.5, size = 3, fontface = "bold") +
    geom_text(data = max_auc_data[max_auc_data$group == "Merged", ], aes(x = 10, y = avg.auc - 0.1, label = paste(round(avg.auc, 2)), color = group), 
              vjust = -0.5, hjust = 1.5, size = 3, fontface = "bold") +
    theme(legend.position = "none") +   
    annotate("segment", x = 500, xend = 540, y = 80, yend = 80, color = "#ff7314", linetype = "solid") +
    annotate("text", x = 550, y = 80, label = "virus", hjust = 0, vjust = 0.3, size = 4) +
    annotate("segment", x = 500, xend = 540, y = 79, yend = 79, color = "#56B4E9", linetype = "solid") +
    annotate("text", x = 550, y = 79, label = "Bacteria", hjust = 0, vjust = 0.3, size = 4) +
    annotate("segment", x = 500, xend = 540, y = 78, yend = 78, color = "#64cd64", linetype = "solid") +
    annotate("text", x = 550, y = 78, label = "Virus+Bacteria", hjust = 0, vjust = 0.3, size = 4) 
}

p1 = get_plot("bac/bac_merge.tsv", "virus/virus_merge.tsv", "merge/merge_merge.tsv", "Disease", "Optimal AUC")

p1
ggsave("/00_data/lineplot_AUC.pdf", plot = p1, width = 10, height = 5, dpi = 300)


