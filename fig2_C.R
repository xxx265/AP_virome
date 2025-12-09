setwd("/share/data3/zhangyt/AP/01_statistical")
rm(list = ls())
load("sample_color.RData")
load("group_order.RData")

# --- significant species select-------
dt = read.table("/share/data3/zhangyt/AP/00_data/AP.tpm", check.names=F, row.names = 1, sep="\t", header=T)
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dt = dt[rowMeans(dt)>0.01,]
x <- zy_pvalue(dt, sample_map, group="Group", ID="Sample",p.method="wilcox.test")
x$qvalue = p.adjust(x$pvalue, method = "BH")
write.table(x, "AP_pvalue.tsv", sep="\t", quote = FALSE)


#------ plot --------
dt = read.table("AP_pvalue.tsv", sep="\t", header=T, row.names=1, check.names=F)
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dt$qvalue = p.adjust(dt$pvalue, method="BH")
dt$log2FC = ifelse(dt$enriched == "Control",log2(dt$fold_change), -log2(dt$fold_change))
dt$log2FC[is.infinite(dt$log2FC) & dt$log2FC == Inf] <- 10
dt$log2FC[is.infinite(dt$log2FC) & dt$log2FC == -Inf] <- -10

dt <- dt %>%  mutate(color = ifelse((qvalue<0.05 & fold_change > 1.2 & (count1 > 3 | count2 > 3)),enriched, "na"))
temp_label <- dt %>%
  select(name, log2FC, qvalue, enriched) %>%
  filter(qvalue < 0.01) %>%
  arrange(qvalue) 

control_top5 <- temp_label %>%
  filter(enriched == "Control") %>%
  head(5)
disease_top5 <- temp_label %>%
  filter(enriched == "Disease") %>%
  head(5)
result <- bind_rows(control_top5, disease_top5)
dt = dt %>% mutate(label=ifelse(name %in% result$name, name,NA))

p <- ggplot(dt, aes(x=log2FC, y=-log10(qvalue), color=color)) +
  geom_point(alpha=0.8) +
  theme_bw() +
  xlim(c(-11, 11)) +
  geom_hline(yintercept=-c(log10(0.05), log10(0.01)), colour="red", linetype="dashed") +
  geom_vline(xintercept=c(log(1.2, 2), -log(1.2, 2)), colour="red", linetype="dashed") +
  scale_color_manual(values=sample_color, na.value = "#c0c0c0") +
  geom_text_repel(aes(label=label), point.padding = 0.8, min.segment.length=0,
                  box.padding = 0.8, segment.color="black", color="black")  # 标签颜色设置为黑色
p
ggsave("fig2_C_AP_volcano.pdf",p, width=6, height=6)

#-----Supplementary Table S3-------
dt = read.table("04_value/AP_sig_pvalue.tsv", sep="\t", header=T, row.names=1, check.names=F)
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dt$qvalue = p.adjust(dt$pvalue, method="BH")
dt$log2FC = ifelse(dt$enriched == "Control",log2(dt$fold_change), -log2(dt$fold_change))
dt$log2FC[is.infinite(dt$log2FC) & dt$log2FC == Inf] <- 10
dt$log2FC[is.infinite(dt$log2FC) & dt$log2FC == -Inf] <- -10
write.table(dt, "Supplementary_Table_S3.tsv", sep="\t", quote = FALSE)
