library(dplyr)

# data
kegg_dt <- read.table("/share/data1/zhangy2/00.pub_data/10.China_vir_catalog_20241125/vOTU.function.kegg.m8.meta", sep="\t", header=F, check.names=F)
kegg_fu <- read.table("/share/data1/Database/KEGG/20241226/ko_list", sep="\t", header=T, check.names=F)

enrich = read.table("04_value/AP_sig.py_pvalue", sep="\t", header=T, check.names=F)

kegg_dt <- kegg_dt %>%
  mutate(V1_prefix = sub("_.*", "", V1)) 
kegg_dt <- kegg_dt %>%
  as_tibble() %>%
  left_join(enrich, by = c("V1_prefix" = "name")) %>%
  filter(!is.na(enriched)) %>%
  dplyr::select(V1, V2, enriched)
unique_values <- unique(kegg_dt[[2]])

#fisher
unique_ke <- unique(kegg_dt$V2)

total_disease <- sum(enrich$enriched == "Disease")
total_control <- sum(enrich$enriched == "Control")

dt <- data.frame(
  kegg = unique_ke,
  stringsAsFactors = FALSE
)

dt$disease_enrich <- sapply(unique_ke, function(ko) {
  sum(kegg_dt$enriched[kegg_dt$V2 == ko] == "Disease")
})

dt$control_enrich <- sapply(unique_ke, function(ko) {
  sum(kegg_dt$enriched[kegg_dt$V2 == ko] == "Control")
})

dt$no_disease_enrich <- total_disease - dt$disease_enrich
dt$no_control_enrich <- total_control - dt$control_enrich

#Fisher 
res <- data.frame(
  kegg = character(),
  p_value = numeric()
)

for (i in 1:nrow(dt)) {
  table <- matrix(c(
    dt$control_enrich[i], dt$disease_enrich[i],
    dt$no_control_enrich[i], dt$no_disease_enrich[i]
  ), nrow = 2, byrow = TRUE)
  
  fisher_result <- fisher.test(table)
  
  res <- rbind(res, data.frame(
    kegg = dt$kegg[i],
    p_value = fisher_result$p.value
  ))
}

# fisher function
run_fisher <- function(row) {
  mat <- matrix(c(row['disease_enrich'], row['control_enrich'],
                  row['no_disease_enrich'], row['no_control_enrich']),
                nrow = 2, byrow = TRUE)
  test <- fisher.test(mat)
  return(c(odds_ratio = test$estimate, p_value = test$p.value))
}


results <- t(apply(dt[, -1], 1, run_fisher))
final_df <- cbind(dt, results)

final_df$disease <- (dt$disease_enrich/total_disease) *100
final_df$control <- (dt$control_enrich/total_control) * 100
final_df$total<- final_df$disease + final_df$control

# plot data
fit <- final_df %>%
  filter((disease > 5 | control > 5) & p_value < 0.05)
fit <- fit %>%
  mutate(group = ifelse(control > disease, "Control", "Disease"))

dt_sorted <- fit %>%
  arrange(
    group,  
    ifelse(group == "Control", -control, disease) 
  ) %>%
  select(-group) 

plot_dt <- dt_sorted[, c("kegg", "disease", "control", "p_value")]

final_df$stars <- ifelse(final_df$p_value < 0.001, "***",
                         ifelse(final_df$p_value < 0.01, "**",
                                ifelse(final_df$p_value < 0.05, "*", "")))

plot_dt_long <- reshape2::melt(plot_dt, id.vars = "kegg", variable.name = "type", value.name = "value")
plot_dt_long <- plot_dt_long %>%
  left_join(final_df %>%
              dplyr::select(kegg, stars), by = "kegg") 
plot_dt_long <- plot_dt_long %>%
  mutate(stars = if_else(type != "control", NA_character_, stars))
plot_dt_long <- plot_dt_long %>%
  filter(type != "p_value")
plot_dt_long <- plot_dt_long %>%
  mutate(type = factor(type, levels = c("control_enrich", "other_type1", "other_type2"))) %>%
  filter(type == "control_enrich") %>%
  bind_rows(
    plot_dt_long %>% filter(type != "control_enrich") 
  ) %>%
  mutate(kegg = factor(kegg, levels = rev(unique(kegg)))) 

# plot
p1 <- ggplot(plot_dt_long, aes(kegg, y=value, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", size = 0.5) +  
  coord_flip() + 
  scale_fill_manual(
    values = c("disease" = "#fb8072", "control" = "#80b1d3"),
    labels = c("disease" = "disease_enrich vOTUs", "control" = "control_enrich vOTUs")
  )+
  theme(
    panel.background = element_rect(fill = "white", color = "white"),  
    plot.background = element_rect(fill = "white", color = "white"),  
    panel.grid.major = element_line(color = "lightgray", size = 0.25), 
    panel.grid.minor = element_line(color = "lightgray", size = 0.25),  
    axis.title = element_blank(),  
    axis.text = element_text(color = "black"), 
    axis.ticks = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title.x = element_text(size = 13, face = "bold"),  
    axis.title.y = element_text(size = 13, face = "bold"),  
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),  
    legend.position = c(1, 0.2)  
  ) +
  labs(
    y = "% Occurrence rate",  
    x = "KEGG",  
    title = "KEGG Enrichment Analysis"  
  ) +
  geom_text(aes(x = kegg, y = as.numeric(value) + 2, label = stars), vjust = 1, size = 5, check_overlap = TRUE)  
p1
ggsave("fig2_E_kegg.pdf",p1, width=6,height=6)



