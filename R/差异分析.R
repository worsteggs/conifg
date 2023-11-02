ibrary(edgeR)
library(dplyr)
library(limma)
library(ggplot2)
library(stringr)

# load data ----
eset = count

# 标准化文库
eset = normalizeBetweenArrays(eset)

# 构建样本信息
f <- factor(c("FGR", "FGR", "FGR", "FGR", "FGR", "FGR", "FGR", "FGR","control", "control", "control", "control", "control", "control", "control", "control"))
design <- model.matrix(~ 0 + f)
design
colnames(design) <- levels(f)

# 构建比较矩阵
contrast <- c("FGR-control")
cont.matrix <- makeContrasts(contrasts = contrast, levels=design)
cont.matrix

# 拟合
fit <- lmFit(eset, design) 
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 导出结果
res <- topTable(fit2,coef=1,adjust.method="BH",n=Inf,sort.by="P", resort.by="logFC")

# 标注上下调信息
res<-arrange(res,desc(B))
res$status<-as.factor(
  ifelse((res$logFC > 1& res$adj.P.Val<0.05),"Log2FC=Up&P",
         ifelse((res$logFC < (-1)&res$adj.P.Val<0.05),"Log2FC=Down&P",
                ifelse(res$logFC<1&res$logFC>(-1)&res$adj.P.Val>0.05,"NS","adj.P.Val>0.05"))))
table(res$status)