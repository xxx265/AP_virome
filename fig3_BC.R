library(patchwork)
library(readxl)
load("sample_color.RData")
load("group_order.RData")

source("/share/data1/zhangy2/scripts/R_my_functions/zy_PCoA.R")
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)

#------fig3 B PCoA etiology-------
dt = read.table("AP_sig.tpm", check.names=F, row.names = 1, sep="\t", header=T)
sample_map = read.table("../APd.group", sep="\t", header=T, check.names=F)
p0 <- zy_pcoa(dt, sample_map, group="Group", ID="Sample", sample.color = sample_color, levels=0.8, star_plot=F, title = "PCoA of virus")
p0$plot
ggsave("PCoA_etiology.pdf", p0$plot, width = 6, height = 5)

#------fig3 C severity-------
dt = read.table("AP_sig.tpm", check.names=F, row.names = 1, sep="\t", header=T)
sample_map = read.table("../APs.group", sep="\t", header=T, check.names=F)
p0 <- zy_pcoa(dt, sample_map, group="Group", ID="Sample", sample.color = sample_color, levels=0.8, star_plot=F, title = "PCoA of virus")
p0$plot
ggsave("PCoA_severity.pdf", p0$plot, width = 6, height = 5)














