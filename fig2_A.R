load("sample_color.RData")
load("group_order.RData")

source("/share/data1/zhangy2/scripts/R_my_functions/zy_compositions.R")
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)

taxo_color = c("#E8F7EF", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#ff6464",
               "#64cd64", "#ffb03f", "#6464ff", "#34597b",
               "#43b095", "#d35372", "#eeb872")

# virus family###生成family的tpm
dt = read.table("03_composition/AP_family.tpm", sep="\t", header=T, row.names=1, check.names=F)
p1 <- zy_group_compositions(dt, sample_map, ID="Sample",group="Group", taxo.color = taxo_color, top_N = 11,
                            title="Composition of virus at family level")+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_blank())+
  coord_cartesian(ylim=c(0,100))+
  scale_y_continuous(breaks=c(25,50,75,100))
p1

ggsave(p1, filename="fig2_A_composition.pdf", width=9, height = 6)

