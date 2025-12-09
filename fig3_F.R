load("sample_color.RData")
load("group_order.RData")

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

get_plot = function(in_MAP, in_MSAP, in_SAP, pred1, pred2, pred3, title){
  
  plot_data_MAP <- get_plot_data(in_MAP, pred1)
  plot_data_MSAP <- get_plot_data(in_MSAP, pred2)
  plot_data_SAP <- get_plot_data(in_SAP, pred3)
  
  all_data <- rbind(
    data.frame(nspecies = plot_data_MAP$nspecies, avg.auc = plot_data_MAP$avg.auc, sd.auc = plot_data_MAP$sd.auc, group = 'MAP'),
    data.frame(nspecies = plot_data_MSAP$nspecies, avg.auc = plot_data_MSAP$avg.auc, sd.auc = plot_data_MSAP$sd.auc, group = 'MSAP'),
    data.frame(nspecies = plot_data_SAP$nspecies, avg.auc = plot_data_SAP$avg.auc, sd.auc = plot_data_SAP$sd.auc, group = 'SAP')
  )
  max_auc_data <- all_data %>%
    group_by(group) %>%
    dplyr::filter(avg.auc == max(avg.auc)) %>%
    dplyr::select(nspecies, avg.auc, group)
  
  ggplot(all_data, aes(x = nspecies, y = avg.auc, color = group, group = group)) +
    geom_line(aes(color = group), size = 1) +
    geom_errorbar(aes(ymax = avg.auc + sd.auc, ymin = avg.auc - sd.auc, color = group), alpha = 0.5) +
    geom_point(data = max_auc_data, aes(x = nspecies, y = avg.auc), color = "black", size = 2) +
    geom_segment(data = max_auc_data, aes(x = nspecies, xend = nspecies, y = 60, yend = avg.auc, color = group), linetype = "dashed") +
    geom_segment(data = max_auc_data, aes(x = 0, xend = nspecies, y = avg.auc, yend = avg.auc, color = group), linetype = "dashed") +
    scale_color_manual(values = c("MAP" = "#8E7FB8", "MSAP" = "#DCE125", "SAP" = "#3FA0C0")) +
    labs(title = title, x = "nspecies", y = "avg_auc") +
    theme_minimal() +
    coord_cartesian(ylim = c(50, NA)) +
    theme(legend.position = "right") +
    annotate("segment", x = 0, xend = Inf, y = 75, yend = 75, color = "black", linetype = "solid") +
    annotate("segment", x = 0, xend = 0, y = 60, yend = Inf, color = "black", linetype = "solid") +
    coord_cartesian(xlim = c(0, max(all_data$nspecies) + 10), ylim = c(75, 100))+
    geom_text(data = max_auc_data, aes(x = nspecies, y = 90, label = paste(round(nspecies, 2)), color = group), 
              vjust = 25, hjust = 0.5, size = 3, fontface = "bold")+
    geom_text(data = max_auc_data[max_auc_data$group == "MAP", ], aes(x = 4, y = avg.auc, label = paste(round(avg.auc, 2)), color = group), 
              vjust = 0, hjust = 1.5, size = 3, fontface = "bold") +
    geom_text(data = max_auc_data[max_auc_data$group == "MSAP", ], aes(x = 4, y = avg.auc - 0.7, label = paste(round(avg.auc, 2)), color = group), 
              vjust = 0, hjust = 1.5, size = 3, fontface = "bold") +
    geom_text(data = max_auc_data[max_auc_data$group == "SAP", ], aes(x = 4, y = avg.auc, label = paste(round(avg.auc, 2)), color = group), 
              vjust = 0, hjust = 1.5, size = 3, fontface = "bold") +
    theme(legend.position = "none") +   
    annotate("segment", x = 235, xend = 245, y = 82, yend = 82, color = "#8E7FB8", linetype = "solid") +
    annotate("text", x = 250, y = 82, label = "MAP", hjust = 0, vjust = 0.3, size = 4) +
    annotate("segment", x = 235, xend = 245, y = 80, yend = 80, color = "#DCE125", linetype = "solid") +
    annotate("text", x = 250, y = 80, label = "MSAP", hjust = 0, vjust = 0.3, size = 4) +
    annotate("segment", x = 235, xend = 245, y = 78, yend = 78, color = "#3FA0C0", linetype = "solid") +
    annotate("text", x = 250, y = 78, label = "SAP", hjust = 0, vjust = 0.3, size = 4) 
}


p1 = get_plot("MAP/KF_merge.tsv", "MSAP/KF_merge.tsv", "SAP/KF_merge.tsv", "MAP", "MSAP", "SAP", "Optimal AUC")

p1
ggsave("fig3_F_AUC.pdf", plot = p1, width = 10, height = 5, dpi = 300)


