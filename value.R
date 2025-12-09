setwd("/share/data3/zhangyt/AP/01_statistical")
rm(list = ls())


#---virus---
#/share/data1/zhangy2/scripts/calc_pvalue.py -d /share/data3/zhangyt/AP/00_data/APAC.tpm -g /share/data3/zhangyt/AP/01_statistical/sample.group -n "Sample" -G "Group" -o APACvirus.py_pvalue

#all_pvalue
dt = read.table("/share/data3/zhangyt/AP/00_data/AP.tpm", check.names=F, row.names = 1, sep="\t", header=T)
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dt = dt[rowMeans(dt)>0.01,]
x <- zy_pvalue(dt, sample_map, group="Group", ID="Sample",p.method="wilcox.test")
x$qvalue = p.adjust(x$pvalue, method = "BH")
write.table(x, "04_value/AP_pvalue.tsv", sep="\t", quote = FALSE)



#-----bac------
#all_pvalue
dt = read.table("/share/data3/zhangyt/AP/00_data/AP_bac.txt", check.names=F, row.names = 1, sep="\t", header=T)
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dt = dt[rowMeans(dt)>0.01,]
colnames(dt) <- paste("pancreatitis_", colnames(dt), sep = "")
x <- zy_pvalue(dt, sample_map, group="Group", ID="Sample",p.method="wilcox.test")
x$qvalue = p.adjust(x$pvalue, method = "BH")
write.table(x, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_pvalue.tsv", sep="\t", quote = FALSE)

#-----serum------
#all_pvalue
dt = read.table("/share/data3/zhangyt/AP/01_statistical/10_corr/AP_serum.tpm", check.names=F, row.names = 1, sep="\t", header=T)
sample_map = read.table("../AP.group", sep="\t", header=T, check.names=F)
dt = dt[rowMeans(dt)>0.01,]
colnames(dt) <- paste("pancreatitis_", colnames(dt), sep = "")
x <- zy_pvalue(dt, sample_map, group="Group", ID="Sample",p.method="wilcox.test")
x$qvalue = p.adjust(x$pvalue, method = "BH")
write.table(x, "04_value/AP_serum_pvalue.tsv", sep="\t", quote = FALSE)



#sig_tpm
######virus
AP <- read.table("04_value/AP_pvalue.tsv", sep="\t", header=TRUE, row.names=1)
filtered_virus <- AP[AP$fold_change >= 1.2 & AP$qvalue < 0.05 & (AP$count1 > 3 | AP$count2 > 3), ]#1049
AP_sig_tpm <- dt[filtered_virus$name, ]
write.table(AP_sig_tpm, "04_value/AP_sig.tpm", sep="\t", quote = FALSE)
first_virus <- filtered_virus[, 1, drop = FALSE]
write.table(first_virus, "04_value/AP_sig.tsv", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(filtered_virus, "04_value/AP_sig_pvalue.tsv", sep="\t", quote = FALSE)


#bac
AP <- read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_pvalue.tsv", sep="\t", header=TRUE, row.names=1)
filtered_virus <- AP[AP$fold_change >= 1.2 & AP$qvalue < 0.05 & (AP$count1 > 3 | AP$count2 > 3), ]#68
AP_sig_tpm <- dt[filtered_virus$name, ]
write.table(AP_sig_tpm, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_sig.tpm", sep="\t",quote = FALSE)
first_virus <- filtered_virus[, 1, drop = FALSE]
write.table(first_virus, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_sig.tsv", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(first_virus, 
            "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_sig.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

#serum
AP <- read.table("/share/data3/zhangyt/AP/01_statistical/04_value/AP_serum_pvalue.tsv", sep="\t", header=TRUE, row.names=1)
filtered_virus <- AP#15
AP_sig_tpm <- dt[filtered_virus$name, ]
write.table(AP_sig_tpm, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_sig.tpm", sep="\t",quote = FALSE)
first_virus <- filtered_virus[, 1, drop = FALSE]
write.table(first_virus, "/share/data3/zhangyt/AP/01_statistical/04_value/AP_serum_sig.tsv", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(first_virus, 
            "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_sig.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)

#sig_pvalue---rf才需要
##virus
dt <- read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/virus/AP_virus.py_pvalue", sep="\t", header=TRUE, row.names=1)
dt1 = read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/virus/virus_rf_imp.sorted", header=T, check.names=F, row.names=1)
matching_rows <- dt[dt$name %in% row.names(dt1), ]

virus <- matching_rows[match(row.names(dt1), matching_rows$name), ]
sample_map = read.table("/share/data3/zhangyt/AP/AP.group", sep="\t", header=T, check.names=F)
virus$qvalue = p.adjust(virus$pvalue, method="BH")
virus$log2FC = ifelse(virus$enriched == "Control",log2(dt$fold_change), -log2(dt$fold_change))
virus$log2FC[is.infinite(virus$log2FC) & dt$log2FC == Inf] <- 10
virus$log2FC[is.infinite(virus$log2FC) & dt$log2FC == -Inf] <- -10
virus <- virus %>%  mutate(color = ifelse((qvalue<0.05 & fold_change > 1.2),enriched, "na"))
virus_sig.py_pvalue <- virus[, c("name", "enriched", "fold_change", "qvalue", "log2FC", "color")]

write.table(virus_sig.py_pvalue, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/virus/AP_virus_sig.py_pvalue", sep="\t", row.names = FALSE)

##bac
dt <- read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac.py_pvalue", sep="\t", header=TRUE, row.names=1)
dt1 = read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/bac_rf_imp.sorted", header=T, check.names=F, row.names=1)
matching_rows <- dt[dt$name %in% row.names(dt1), ]

bac <- matching_rows[match(row.names(dt1), matching_rows$name), ]
sample_map = read.table("/share/data3/zhangyt/AP/AP.group", sep="\t", header=T, check.names=F)
bac$qvalue = p.adjust(bac$pvalue, method="BH")
bac$log2FC = ifelse(bac$enriched == "Control",log2(dt$fold_change), -log2(dt$fold_change))
bac$log2FC[is.infinite(bac$log2FC) & dt$log2FC == Inf] <- 10
bac$log2FC[is.infinite(bac$log2FC) & dt$log2FC == -Inf] <- -10
bac <- bac %>%  mutate(color = ifelse((qvalue<0.05 & fold_change > 1.2),enriched, "na"))
bac_sig.py_pvalue <- bac[, c("name", "enriched", "fold_change", "qvalue", "log2FC", "color")]

write.table(bac_sig.py_pvalue, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/bac/AP_bac_sig.py_pvalue", sep="\t", row.names = FALSE)


##merge
dt <- read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/merge/AP_merge.py_pvalue", sep="\t", header=TRUE, row.names=1)
dt1 = read.table("/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/merge/merge_rf_imp.sorted", header=T, check.names=F, row.names=1)
matching_rows <- dt[dt$name %in% row.names(dt1), ]

merge <- matching_rows[match(row.names(dt1), matching_rows$name), ]
sample_map = read.table("/share/data3/zhangyt/AP/AP.group", sep="\t", header=T, check.names=F)
merge$qvalue = p.adjust(merge$pvalue, method="BH")
merge$log2FC = ifelse(merge$enriched == "Control",log2(dt$fold_change), -log2(dt$fold_change))
merge$log2FC[is.infinite(merge$log2FC) & dt$log2FC == Inf] <- 10
merge$log2FC[is.infinite(merge$log2FC) & dt$log2FC == -Inf] <- -10
merge <- merge %>%  mutate(color = ifelse((qvalue<0.05 & fold_change > 1.2),enriched, "na"))
merge_sig.py_pvalue <- merge[, c("name", "enriched", "fold_change", "qvalue", "log2FC", "color")]

write.table(merge_sig.py_pvalue, "/share/data3/zhangyt/AP/01_statistical/07_randomforest/00_data/merge/AP_merge_sig.py_pvalue", sep="\t", row.names = FALSE)








#--- virus family
dt = read.table("03_composition/AP_family.tpm", check.names=F, row.names = 1, sep="\t", header=T)
dt = dt[rowMeans(dt)>0.01,]
y <- zy_pvalue(dt, sample_map, group="Group", ID="Sample",p.method="wilcox.test")
y$qvalue = p.adjust(y$pvalue, method = "BH")
write.table(y, "04_value/AP_family_pvalue.tsv", sep="\t")

data_virus_fanliy <- read.table("04_value/AP_family_pvalue.tsv", sep="\t", header=TRUE, row.names=1)
filtered_fanliy <- data_virus_fanliy[data_virus_fanliy$fold_change >= 1.2 & data_virus_fanliy$qvalue < 0.05 & (data_virus_fanliy$count1 > 3 | data_virus_fanliy$count2 > 3), ] #12
AP_family_sig_tpm <- dt[filtered_fanliy$name, ]
write.table(AP_family_sig_tpm, "04_value/AP_famliy_sig.tpm", sep="\t")

first_famliy <- filtered_fanliy[, 1, drop = FALSE]
write.table(first_famliy, "04_value/AP_family_sig.tsv", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


























zy_pvalue = function(dt=NA, sample_map=NA, group="Group", ID="Sample", p.method="wilcox.test"){
  # 如果有多个分组，all.avg只代表当前两个分组的均值，min_avg, min_fd...也一样
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
  # dt = dt[, colSums(dt)!=0]这一步不应该有
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
