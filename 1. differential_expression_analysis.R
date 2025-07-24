# 加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2"))

library(limma)
library(edgeR)
library(DESeq2)

# 读取数据
counts <- read.table("TCGA-BRCA.star_counts.tsv", header=TRUE, row.names=1, sep="\t")

# 数据预处理
# 1. 过滤低表达基因（更严格的过滤标准）
keep <- rowMeans(counts) >= 10
counts <- counts[keep,]

# 2. 创建分组信息
sample_type <- ifelse(grepl(".11", colnames(counts)), "Normal", "Tumor")
design <- model.matrix(~0 + factor(sample_type))

colnames(design) <- c("Normal", "Tumor")
contrast.matrix <- makeContrasts(Tumor-Normal, levels=design)

# 3. 使用edgeR进行数据标准化
dge <- DGEList(counts = counts, group = sample_type)
dge <- calcNormFactors(dge, method="TMM")

# 4. 使用voom进行数据转换（添加权重）
v <- voom(dge, design, plot=TRUE, normalize.method="quantile")

# 5. 差异表达分析
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit, trend=TRUE)

# 6. 提取差异表达结果
results <- topTable(fit, coef=1, number=Inf, adjust.method="BH")
results$gene_id <- rownames(results)

# 6. 设置显著性阈值
results$significant <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1, 
                            ifelse(results$logFC > 1, "Up", "Down"), "NS")

# 8. 保存结果
write.csv(results, "differential_expression_results.csv", row.names=FALSE)

# 9. 生成火山图
pdf("volcano_plot2.pdf", width=10, height=8)
# 设置马卡龙配色
colors <- c("Up" ="#C1565E",    # 马卡龙红色  #ff4757
           "Down" ="#7EA4D1", # 马卡龙蓝色  #546de5
           "NS" = "#807C7D")    # 马卡龙灰色 #d2dae2

# 创建基础图形
plot(results$logFC, -log10(results$adj.P.Val),
     pch=20, 
     main="Volcano Plot",
     xlab="log2 Fold Change",
     ylab="-log10 adjusted p-value",
     col=colors[results$significant],
     cex=0.8)

# 添加阈值线
abline(h=-log10(0.05), col="#666666", lty=2, lwd=1.5)
abline(v=c(-1,1), col="#666666", lty=2, lwd=1.5)

# 添加图例
legend("topright",
       legend=c("Up-regulated", "Down-regulated", "Not significant"),
       col=c("#C1565E","#7EA4D1", "#807C7D"), #c("#FF9999", "#99CCFF", "#CCCCCC"),
       pch=20,
       cex=0.8,
       bty="n")

# 添加标题
title(main="Differential Expression Analysis",
      #sub="|log2FC| > 1 & adj.P.Val < 0.05",
      #cex.main=1.2,
      #cex.sub=0.8)
)

dev.off()

# 10. 输出统计信息
cat("\n差异表达基因统计：\n")
cat("上调基因数量：", sum(results$significant == "Up"), "\n")
cat("下调基因数量：", sum(results$significant == "Down"), "\n")
cat("非显著差异基因数量：", sum(results$significant == "NS"), "\n")

# 11. 保存显著差异基因列表
sig_genes <- results[results$significant != "NS",]
write.csv(sig_genes, "significant_genes.csv", row.names=FALSE)

# 12. 输出数据质量信息
cat("\n数据质量信息：\n")
cat("原始基因数量：", nrow(counts), "\n")
cat("过滤后基因数量：", sum(keep), "\n")
cat("样本数量：", ncol(counts), "\n")
cat("肿瘤样本数量：", sum(sample_type == "Tumor"), "\n")
cat("正常样本数量：", sum(sample_type == "Normal"), "\n") 