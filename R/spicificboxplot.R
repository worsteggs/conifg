#!/usr/bin/env Rscript

# 1. 程序功能描述和主要步骤----

# 程序功能：Anova组间统计和箱线图展示
# Anova Script functions: Calculate pvalue of groups by aov and TukeyHSD
# Main steps: 
# - Reads data table input.txt
# - Calculate pvalue and save in output.txt
# - Draw boxplot and save in output.pdf

# 程序使用示例
# USAGE
# Default
# anova.r   -i data_table.txt
#                       -o otuput filename prefix for output directory name 

# 参数说明
# Options
# -i/--input    输入数据表文件 input.txt
# -o/--output   输出结果文件名前缀 output_prefix, 通常会有统计表txt和矢量图pdf
options(warn = -1)


# 2. 依赖关系检查、安装和加载----
# See whether these packages exist on comp. If not, install.
package_list <- c("optparse","reshape2","ggplot2","ggpubr","dplyr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}
# 另两种常见R包安装方法
if (FALSE){
  # Bioconductor安装
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("reshape2"))
  # Github安装
  install.packages("devtools", repo="http://cran.us.r-project.org")
  library(devtools)
  install_github("kassambara/ggpubr")
}

# 清理工作环境 clean enviroment object
rm(list=ls()) 

# 加载依赖关系 Load essential packages
library(optparse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)

# 解析命令行
if (TRUE){
  option_list <- list(
    make_option(c("-i", "--input"), type="character", default="../tRF.txt",
                help="Input table file to read [default %default]"),
    make_option(c("-g", "--group"), type="character", default="../标本分组.csv",
                help="Input table file to read [default %default]"),
    make_option(c("-o", "--output"), type="character", default="./out",
                help="output directory or prefix [default %default]"),
    make_option(c("-R", "--rna"), type="character", default=c("hsa-miR-130a-3p", "hsa-miR-29b-3p", "hsa-miR-3613-5p", "hsa-miR-24-3p", "hsa-miR-335-5p", "hsa-miR-181a-5p",'tRF-3-tRNA-Val-CAC-3-1','tRF-3-tRNA-Asp-GTC-2-7'),
                help="output directory or prefix [default %default]")
  )
  opts <- parse_args(OptionParser(option_list=option_list))

  # 显示输入输出确认是否正确
  print(paste("The input file is ", opts$input,  sep = ""))
  print(paste("The gruop file is ", opts$group,  sep = ""))
  print(paste("The output file prefix is ", opts$output, sep = ""))
  print(paste("The specific rna is ", opts$rna,  sep = ""))
}

# 3. 读取输入文件----
# 需要使用哪种方式，将其设置为TRUE，其它为FALSE即可

# 产生两组数据比较
if (FALSE){
  # 产生两种各10个正态分布数值，分别以均值为1和2，标准差为均值的0.5倍，colnames修改组名称
  set.seed(123)
  A = rnorm(10, mean = 1, sd = 0.5)
  B = rnorm(10, mean = 2, sd = 1)
  dat=as.data.frame(cbind(A,B))
  colnames(dat)=c("GroupA","GroupB")
  # 保存数据文件方便下面演示
  write.table(dat, file="opts$input", append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = F, col.names = T)
}
# 从文件中读取
if (TRUE){
  counts = read.table(opts$input, header = T, sep = '')#导入tRF数据
  patient_info <- read.csv(opts$group,header = T) # 读入样本分组信息
}
# 弹出窗口选择文件
if (FALSE){
  dat = read.table(file.choose(), header=T, row.names = NULL, sep="\t")
}
# 4. 数据处理----
#构建矩阵
eset = counts[,2:16]
rownames(eset) = counts[,1]
eset<-eset[which(rowSums(eset)>1),] 
#log
logTPM <- log2(eset+1)

#去掉单样本病人
logTPM <- logTPM[,!colnames(logTPM) %in% c('GYL.2','GYL_2')]

design <- model.matrix(~0+(patient_info$clinic_type))
colnames(design) <- c('healthy','moderate','severe')
rownames(design) <- patient_info$sample #构建design，用于后续相关性分析
design <- design[match(colnames(logTPM),rownames(design)),]#对齐design和eset(不一定需要，只是保险，但是样本名的顺序必须要一致！！！否则会出错！！！！)

#拟合统计模型
fit <- lmFit(logTPM, design) 
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(mh = moderate-healthy,
                                 sh = severe-healthy, 
                                 sm = severe-moderate,
                                 # ch = (severe+moderate)-healthy,
                                 levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
mh.data <- topTable(fit2,coef=1,adjust.method="BH",n=Inf,sort.by="P", resort.by="logFC")
sh.data <- topTable(fit2,coef=2,adjust.method="BH",n=Inf,sort.by="P", resort.by="logFC")
sm.data <- topTable(fit2,coef=3,adjust.method="BH",n=Inf,sort.by="P", resort.by="logFC")

#构建显著性标注表格
p.data = c()
ymix = c()
for (j in opts$rna) {
  temp = data.frame(P.Value = c(mh.data[j,]$P.Value, sh.data[j,]$P.Value, sm.data[j,]$P.Value))
  temp$signif <- ifelse(temp$P.Value > 0.05, "NA", ifelse(temp$P.Value > 0.01,"*", "**"))
  p.data = rbind(p.data, temp)
  maxx = max(logTPM[j,])
  minx = min(logTPM[logTPM > 0,][j,])
  gap = maxx - minx
  ytemp = c(maxx+0.1, maxx+0.1+0.15*gap, maxx+0.1+0.15*gap/2)
  ymix = c(ymix, ytemp)
}
p.data$start = rep(c('healthy', 'healthy', 'moderate'), length(opts$rna))
p.data$end = rep(c('moderate', 'severe', 'severe'), length(opts$rna))
p.data$y = ymix

p.data$sRNA = rep(opts$rna, each=3)

to_plot_exp <- logTPM
#构建长数据
t_to_plot_exp<-t(to_plot_exp)
melt_t_to_plot_exp<-melt(t_to_plot_exp)
colnames(melt_t_to_plot_exp)<-c("sample","sRNA","Log2TPM")
#添加分组信息
melt_t_to_plot_exp$clinic_type <- patient_info$clinic_type[match(melt_t_to_plot_exp$sample, patient_info$sample)]


# 5. 统计与绘图----
if (TRUE){
  p = ggboxplot(melt_t_to_plot_exp[melt_t_to_plot_exp$sRNA %in% opts$rna,],"clinic_type","Log2TPM",
            fill="clinic_type",color = "clinic_type",palette = "Paired",
            font.label = c("bold", 11),
            font.legend = "bold",
            font.x=c(11,"bold","black"),
            font.y=c(11,"bold","black"),
            font.xtickslab=c(11.5,"bold","black"),
            width = 0.5,
            add="mean",
            add.params = list(color = "red"),
            ggtheme=ggplot2::theme_minimal(),
            bxp.errorbar=T,
            bxp.errorbar.width=0.2) +
    geom_jitter(width = .1, alpha = 0.5) +
    theme_bw() +
    theme(axis.text.x=element_text(angle = 0,colour = "black"),
          panel.background=element_rect(colour = "grey"),
          strip.background = element_rect(fill = "lightyellow",size = 1),
          # strip.placement = "inside",
          strip.text = element_text(face = "bold",family = "sans"),
          panel.grid = element_blank()) +
    geom_signif(data=p.data,
                aes(xmin=start, xmax=end, annotations=signif, y_position=y),
                # textsize = 3,
                # vjust = -0.2,
                manual=TRUE,
    textsize = 4,
    vjust = 1.3,
    # y_position = c(15.5,16),
    tip_length = c(0), #设置显著性那条横线两头向下的长度
    # map_signif_level = T, #设置是否标记显著性的*号，还是直接标记数值
    test = wilcox.test) +
    facet_wrap(.~sRNA,scales="free", strip.position = "top")
  print(p)
  ggsave(paste0(opts$output,"/boxplot.pdf"),width = 1.5*length(opts$rna),height = 1.5*length(opts$rna)*2/3)
}



# 6. 保存图表----
if (FALSE){
  # 保存一个制表符，解决存在行名时，列名无法对齐的问题
  write.table("\t", file=paste(opts$output,".txt",sep=""),append = F, quote = F, eol = "", row.names = F, col.names = F)
  # 保存统计结果，有waring正常
  write.table(Tukey_HSD_table, file=paste(opts$output,".txt",sep=""), append = T, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)
  print(paste("The output table is ", opts$output, ".txt",  sep = ""))
  # 保存图片至文件，pdf方便AI修改成出版级图片
  ggsave(file=paste(opts$output,".pdf",sep=""), p, width = 5, height = 3)
  print(paste("The output figure is ", opts$output, ".pdf",  sep = ""))
}