library(dplyr)
library(data.table)

myfun = function(file,out) {
  res = fread(file, fill = T, sep = '\t', header = F)
  res$temp = rep(c(1,2,3,4),61)
  sample = res %>% filter(temp == 1) %>% select(V1)
  mapped = res %>% filter(temp == 2) %>% select(V2)
  multimapped = res %>% filter(temp == 3) %>% select(V2)
  unmapped = res %>% filter(temp == 4) %>% select(V2)
  resNew = data.frame(sample, mapped, multimapped, unmapped)
  colnames(resNew) = c('samples', 'uniMapped', 'multiMapped', 'unMapped')
  write.table(resNew, paste0('/Users/hziliang/Documents/广州实验室/project/新冠鼻拭子单细胞/results/',out,'.txt'), quote = F, sep = '\t',row.names = F)
}

# for map
myfun('/Users/hziliang/Documents/广州实验室/project/新冠鼻拭子单细胞/results/map/results.txt', 'map')
myfun('/Users/hziliang/Documents/广州实验室/project/新冠鼻拭子单细胞/results/mapRm/results.txt', 'mapRm')
myfun('/Users/hziliang/Documents/广州实验室/project/新冠鼻拭子单细胞/results/mapRmPseudomonasncbi/results.txt', 'mapRmPseudomonasncbi')
