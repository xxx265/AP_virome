library(patchwork)
load("sample_color.RData")
load("group_order.RData")

#plot
#diversity plot function
zy_alpha = function(dt=NA, sample_map=NA, group="Group", ID="Sample", 
                    index="shannon", 
                    sample.color=NA, 
                    box_width=0.5, 
                    title="alpha diversity", 
                    violin = F
){
  ## colors 
  if (any(is.na(sample.color))){
    sample.color = c(1:length(unique(sample_map[,group])))
  }
  message(paste(length(sample.color), "of groups to plot"))
  
  ## align dt and group
  dt = dt[,sample_map[,ID]]
  dt = dt[rowSums(dt)!=0,]
  
  #alpha
  if(tolower(index) == "obs"){
    alpha = data.frame(alpha=colSums((dt>0)+0))
  }else{
    alpha = data.frame(alpha = vegan::diversity(t(dt),index=index))
  }
  
  dm = merge(alpha,sample_map, by.x='row.names', by.y=ID)
  comp = combn(as.character(unique(dm[,group])),2,list)
  
  p = ggplot(dm, aes(x=.data[[group]], y=alpha,fill=.data[[group]]))
  if(isTRUE(violin)){
    p <- p+
      geom_violin()+
      geom_boxplot(width=box_width, fill="white",
                   position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.colour = NA)
  }else{
    p <- p+ 
      geom_boxplot(position = position_dodge2(preserve = 'single')
                   ,outlier.shape = 21,outlier.fill=NA, outlier.color="#c1c1c1")
  }
  
  ylabs = structure(c("Number of OTUs","Shannon index", "1 - Simpson index", "Invsimpson index"),
                    names=c("obs", "shannon", "simpson","invsimpson"))
  ylab = ylabs[tolower(index)]
  
  
  p <- p+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_fill_manual(values=sample.color)+
    #geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=numFunc)+
    geom_signif(comparisons =comp,test='wilcox.test',test.args=list(exact=F),step_increase = 0.1)+
    labs(title=title, y = ylab, x=NULL)
  
  p
}

sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
sample_map$Group = factor(sample_map$Group, levels = group_order)
dt = read.table("../00_data/AP.tpm", sep="\t", header=T, row.names=1, check.names=F)
p0 <- zy_nspecies(dt, sample_map, sample.color = sample_color)
p1 <- zy_alpha(dt, sample_map = sample_map, index="obs", title = "vOTU obs", sample.color = sample_color)
p2 <- zy_alpha(dt, sample_map = sample_map, index="shannon", title = "vOTU shannon", sample.color = sample_color)
p3 <- zy_alpha(dt, sample_map = sample_map, index="simpson", title = "vOTU simpson", sample.color = sample_color)

pa <- p0+p1+p2+p3+
  plot_layout(ncol = 4)
pa

ggsave("fig1_abcd_diversity.pdf",pa, width=12,height=5)


#PCoA plot function
source("/share/data1/zhangy2/scripts/R_my_functions/zy_PCoA.R")
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)

dt = read.table("../00_data/AP.tpm", check.names=F, row.names = 1, sep="\t", header=T)
p0 <- zy_pcoa(dt, sample_map, group="Group", ID="Sample", sample.color = sample_color, levels=0.8, star_plot=F, title = "PCoA of virus")
p0$plot
ggsave("fig1_e_PCoA.pdf", p0$plot, width = 6, height = 5)
































