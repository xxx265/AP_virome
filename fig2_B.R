library(ggplot2)
library(ggforce)

sample_map = read.table("../AP.group", sep="\t", header=T)
sample_map $Group = factor(sample_map$Group, levels=c("Control", "Disease"))
load("group_order.RData")
load("sample_color.RData")

#-----significant family-------
###select function
zy_pvalue = function(dt=NA, sample_map=NA, group="Group", ID="Sample", p.method="wilcox.test"){

  # ID -> ID columns name
  # gorup -> how to group data
  # dt -> profile
  # sample_map -> mapping file
  
  
  test.arg = c(wilcox.test, t.test)
  names(test.arg) = c("wilcox.test", "t.test")
  my_test = test.arg[[p.method]]
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }else{
    message("\033[31mInfo\t数据和分组完全匹配\033[0m")
  }
  dt = dt[, sample_map[,ID]]
  
  raw_ncol = ncol(dt)
  raw_nrow = nrow(dt)
  dt = dt[rowSums(dt)!=0,]
  f_ncol = ncol(dt)
  f_nrow = nrow(dt)
  
  if(f_ncol != raw_ncol || raw_nrow != f_nrow){
    message(paste("delete all items is 0 -> columns:", raw_ncol - f_ncol, " rows:" ,raw_nrow - f_nrow, sep=""))
  }
  
  grps = unique(sample_map[,group])
  com = t(combn(grps,2))
  nspecies = nrow(dt)
  names_ = rownames(dt)
  # Avg -> 平均数
  # Avg.weighted.g1 -> 这个分组的加权平均数
  result = data.frame(matrix(NA,nrow = nrow(com)*nspecies, ncol = 19,
                             dimnames = list(NULL,c("name","g1","g2","Avg.g1","Avg.g2","fold_change","enriched",
                                                    "all.avg","all.var","pvalue",
                                                    "count1","count2","total_count1","total_count2",
                                                    "rank1.avg", "rank2.avg","method","var1","var2"))))
  nr = 1
  for (n in 1:nspecies){
    cat("\r",n, " / ", nspecies)
    temp_dt = dt[n,]
    for(c in 1:nrow(com)){
      g1 = com[c,1]
      g2 = com[c,2]
      
      g1s = sample_map[which(sample_map[,group] == g1), ID] # group
      g2s = sample_map[which(sample_map[,group] == g2), ID]
      
      dt1 = as.matrix(temp_dt[,g1s]) # data
      dt2 = as.matrix(temp_dt[,g2s])
      
      c1 = sum(dt1 != 0, na.rm=T)  # count !0
      c2 = sum(dt2 != 0, na.rm=T)
      
      tc1 = length(dt1) # total count
      tc2 = length(dt2)
      
      m1 = mean(dt1, na.rm=T) # mean data
      m2 = mean(dt2, na.rm=T)
      
      var_1 = var(as.numeric(dt1), na.rm=T) # var data
      var_2 = var(as.numeric(dt2), na.rm=T)
      
      am = mean(c(dt1,dt2), na.rm=T) # total mean
      a_var=var(c(dt1,dt2), na.rm=T) # total var
      p = my_test(dt1,dt2)$p.value # pvalue
      fold_change = ifelse(m1>m2, m1/m2, m2/m1) # fold change
      enriched = ifelse(m1>m2, g1,g2) # enriched
      
      m = sample_map[which(sample_map[,group] %in% c(g1,g2)), ID]
      all_rank = rank(temp_dt[,m])
      
      rank1 = all_rank[colnames(dt1)] # rank
      rank2 = all_rank[colnames(dt2)]
      
      rank1.avg = mean(rank1) # mean rank
      rank2.avg = mean(rank2)
      
      result[nr,] = c(names_[n], g1, g2, m1, m2, 
                      fold_change, enriched, am, a_var,
                      p, c1, c2, tc1, tc2,rank1.avg, rank2.avg,
                      method=p.method, var1=var_1, var2=var_2)
      nr = nr+1
    }
  }
  result[,c(4:6,8:16,18,19)] = lapply(result[,c(4:6,8:16,18,19)], as.numeric)
  result
}

##select
dt = read.table("AP_family.tpm", check.names=F, row.names = 1, sep="\t", header=T)
dt = dt[rowMeans(dt)>0.01,]
y <- zy_pvalue(dt, sample_map, group="Group", ID="Sample",p.method="wilcox.test")
y$qvalue = p.adjust(y$pvalue, method = "BH")
write.table(y, "AP_family_pvalue.tsv", sep="\t")

data_virus_fanliy <- read.table("AP_family_pvalue.tsv", sep="\t", header=TRUE, row.names=1)
filtered_fanliy <- data_virus_fanliy[data_virus_fanliy$fold_change >= 1.2 & data_virus_fanliy$qvalue < 0.05 & (data_virus_fanliy$count1 > 3 | data_virus_fanliy$count2 > 3), ] #12
AP_family_sig_tpm <- dt[filtered_fanliy$name, ]
write.table(AP_family_sig_tpm, "AP_famliy_sig.tpm", sep="\t")

first_famliy <- filtered_fanliy[, 1, drop = FALSE]
write.table(first_famliy, "AP_family_sig.tsv", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#-----box-----
sigFunc = function(x){
  if(x < 0.001){formatC(x, digits = 1, width = 1, format = "e", flag = "0")}
  else if(x<0.05){formatC(x, digits = 3, width = 1, format = "f", flag = "0")}
  else{NA}}

getplot = function(x,y){
  dt = read.table(x, sep="\t", header=T, row.names=1, check.names=F)
  mark = read.table(y, sep="\t", header=F)
  dt = dt[mark$V1, ]
  taxo_order = sort(rowMeans(dt), decreasing=T)
  dl = melt(as.matrix(dt))
  dm = merge(dl,sample_map, by.x="Var2", by.y="Sample")
  dm$Var1 = factor(dm$Var1, levels=names(taxo_order))
  p <- ggplot(dm, aes(y=value, x=Group, fill=Group))+
    geom_boxplot(outlier.shape = 21, outlier.color = "grey", outlier.fill = NA)+
    scale_y_sqrt()+
    theme_bw()+
    scale_fill_manual(values=sample_color)+
    facet_wrap(~ Var1, scales = "free_y", ncol = length(unique(dm$Var1)))+
    geom_signif(comparisons =list(c("Control", "Disease")),test='wilcox.test',test.args=list(exact=F),step_increase = 0.1,map_signif_level=sigFunc)+
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      
    )+
    ylab("% Abundance")
  p
}

p = getplot("AP_family.tpm","AP_family_sig.tsv")
p

ggsave("fig2_B_family_box.pdf",p, width=14, height=5)

