library(dplyr)

###data collect

virus_data <- read.table("vb_control_corr", sep="\t", header=T, check.names=F)
bac_data <- read.table("vs_control_corr", sep="\t", header=T, check.names=F)
serum_data <- read.table("bs_control_corr", sep="\t", header=T, check.names=F)

control_corr <- bind_rows(virus_data, bac_data, serum_data)
control_corr$group <- "bac_virus_serum"
control_corr$qvalue <- p.adjust(control_corr$pval, method="fdr")
control_corr$source <- "control"


virus_data1 <- read.table("vb_disease_corr", sep="\t", header=T, check.names=F)
bac_data1 <- read.table("vs_disease_corr", sep="\t", header=T, check.names=F)
serum_data1 <- read.table("bs_disease_corr", sep="\t", header=T, check.names=F)

disease_corr <- bind_rows(virus_data1, bac_data1, serum_data1)
disease_corr$group <- "bac_virus_serum"
disease_corr$qvalue <- p.adjust(disease_corr$pval, method="fdr")
disease_corr$source <- "disease"

all_conn <- bind_rows(control_corr, disease_corr)

write.table(all_conn, "all_conn", sep = "\t")


###corrlation data(fig4 AB)
zy_node_attr = read.table("all_node",sep="\t",comment.char = "", header=T,check.names = F) 
zy_node_attr$id = seq(1:nrow(zy_node_attr))

zy_conn = read.table("all_conn",sep="\t",comment.char = "", header=T, check.names=F)
colnames(zy_conn)[1:2] <- c("from", "to")

zy_conn <- zy_conn %>%
  filter(corr > 0.4 | corr < (-0.4), qvalue < 0.05)
zy_conn = zy_conn %>% filter(source == "disease")#1132
#zy_conn = zy_conn %>% filter(source == "control")#1105
#write.table(zy_conn, "conrtol_conn.tsv", sep = "\t")
write.table(zy_conn, "disease_conn.tsv", sep = "\t")

# 构建关系数据
zy_d1 = data.frame(from='origin', to=unique(zy_node_attr$order))
zy_d1 <- zy_d1 %>%
  arrange(to) 
zy_d2 = data.frame(from=zy_node_attr$order, to=zy_node_attr$name)
zy_d2 <- zy_d2 %>%
  arrange(to) 
zy_tree = rbind(zy_d1, zy_d2)

# 构建每个节点的属性
temp_dt = as.data.frame(matrix(NA,
                               nrow=nrow(zy_d1)+1, ncol=ncol(zy_node_attr),
                               dimnames=list(NULL,colnames(zy_node_attr))
))
temp_dt$name = c("origin",zy_d1$to)
zy_node_attr <- rbind(temp_dt, zy_node_attr)

##############################
temp_color <- zy_node_attr %>%
  dplyr::select(color,temp_taxo) %>%
  unique() %>%
  dplyr::filter(!is.na(temp_taxo))
mycolors = temp_color$color
names(mycolors) = temp_color$temp_taxo
##############################

rownames(zy_node_attr) = zy_node_attr$name
zy_graph_data <- graph_from_data_frame(zy_tree, vertices = zy_node_attr)

zy_from = match(zy_conn$from, zy_node_attr$name)
zy_to   = match(zy_conn$to, zy_node_attr$name)
zy_conn$width = abs(zy_conn$corr)
# zy_conn$color = ifelse(zy_conn$corr>0, "gray", "red")
#zy_conn$color = zy_node_attr[match(zy_conn$from, zy_node_attr$name),"enriched"]
zy_conn$color = ifelse(zy_conn$corr > 0,"#9FDFCF", "#4B2C20" ) 


##plot
p<- ggraph(graph = zy_graph_data, layout='dendrogram', circular=T) +
  geom_conn_bundle(data = get_con(from = zy_from, to = zy_to, corr=zy_conn$corr, edge_color=zy_conn$color, width=zy_conn$width), 
                   aes(edge_colour=edge_color, edge_width=width),
                   edge_width=0.25,# 这个会覆盖aes里面的参数
                   tension=0.6,
                   alpha=0.3) +
  geom_node_point(aes(filter = leaf, fill = temp_taxo, shape = type, size=type)) + #filter 过滤
  geom_node_text(aes(x = x*1.04, y=y*1.04, filter = leaf, label=label, 
                     colour = enriched,
                     angle = atan(y/x)*360/(2*pi), 
                     hjust='outward'),
                 size=3, alpha=1)+
  scale_edge_colour_manual(values=structure(c("#9FDFCF",'#4B2C20'), names=c("#9FDFCF","#2F1B16")))+
  scale_shape_manual(values=c(21,24,22))+
  scale_fill_manual(values=mycolors) +
  expand_limits(x = c(-1.1, 1.1), y = c(-1.1, 1.1)) +
  scale_size_manual(values=c(2.8,2,2)) +
  coord_fixed()+ # 设置为标准圆形+
  theme_void()+
  expand_limits(x = c(-1.1, 1.1), y = c(-1.1, 1.1))+
  #ggtitle("Collection Network in disease \n (n = 1953)")
  ggtitle("Collection Network in control \n (n = 1105)") # 添加标题

p

ggsave("test.pdf", p, width = 40, height = 40, limitsize = FALSE)
ggsave("circle_disease.pdf", p, width = 40, height = 40, limitsize = FALSE)
ggsave("circle_control.pdf", p, width = 40, height = 40, limitsize = FALSE)


# plot legend
pa<- ggplot(temp_color,aes(x=1, y=temp_taxo, fill=temp_taxo))+
  geom_point( shape=21, size=8)+
  geom_text(aes(label=temp_taxo),x=1.1, hjust=0)+
  scale_fill_manual(values=mycolors)+
  xlim(c(1,4))+
  theme_void()

ggsave("legend.pdf",pa, width=5,height=5)

##plot most species(fig 4 CD)
merged_counts <- zy_conn %>%
  pivot_longer(cols = c(from, to), names_to = "direction", values_to = "value") %>%
  count(value) %>%
  mutate(type = ifelse(grepl("^s_", value), "bac", 
                       ifelse(grepl("^v", value), "virus", "other"))) %>%
  arrange(desc(n))

top_20 <- merged_counts %>%
  top_n(20, n)

p <- ggplot(top_20, aes(y = reorder(value, n), x = n, fill = type)) +
  geom_bar(stat = "identity") +
  labs(title = "Top 20 Most Frequent Values", x = "Count", y = "Value") +
  scale_fill_manual(values = c("bac" = "#4E79A7", "virus" = "#E15759", "other" = "grey")) 
  theme(axis.text.y = element_text(angle = 0, hjust = 1))  
p
ggsave("disease_most.pdf",p, width=5,height=5)
#ggsave("control_most.pdf",p, width=5,height=5)