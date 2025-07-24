# 加载必要的包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", "gridExtra"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(gridExtra)

# 读取基因列表
genes <- readLines("genes.txt")

# 设置马卡龙配色
macaron_colors <- c(
    "#D95319",  # 马卡龙红色
    "#6BAED6",  # 马卡龙蓝色
    "#77AC30",  # 马卡龙绿色
    "#FFCC99",  # 马卡龙橙色
    "#CC99FF"   # 马卡龙紫色
)

# GO富集分析
go_bp <- enrichGO(gene = genes,
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

go_cc <- enrichGO(gene = genes,
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

go_mf <- enrichGO(gene = genes,
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)

# KEGG富集分析
kegg <- enrichKEGG(gene = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
                  organism = 'hsa',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

# 准备数据
bp_data <- as.data.frame(go_bp)[1:10,]
cc_data <- as.data.frame(go_cc)[1:10,]
mf_data <- as.data.frame(go_mf)[1:10,]

# 添加类别标签
bp_data$Category <- "Biological Process"
cc_data$Category <- "Cellular Component"
mf_data$Category <- "Molecular Function"

# 合并数据
combined_data <- rbind(bp_data, cc_data, mf_data)

# 创建整合的GO富集分析图
pdf("GO_enrichment_combined2.pdf",width = 10, height = 10)# res=100)
ggplot(combined_data, aes(x=reorder(Description, Count), y=Count, fill=Category)) +
    geom_bar(stat="identity") +
    coord_flip() +
    facet_wrap(~Category, scales="free_y", ncol=1) +
    scale_fill_manual(values=c("Biological Process"=macaron_colors[1],
                              "Cellular Component"=macaron_colors[2],
                              "Molecular Function"=macaron_colors[3])) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size=rel(1)),
          axis.text.x = element_text(size=rel(1)),
          plot.title = element_text(size=rel(1.2), face="bold"),
          legend.position = "none",
          strip.text = element_text(size=rel(1.1), face="bold"),
          strip.background = element_rect(fill="white", color="gray"),
          axis.title.x = element_text(size=rel(1)),
          axis.title.y = element_text(size=rel(1))) +
    labs(title="GO Enrichment Analysis",
         x="GO Terms",
         y="Gene Count")
dev.off()

# 绘制KEGG富集分析结果
pdf("KEGG_enrichment.pdf", width=10, height=10)
dotplot(kegg, showCategory=15, color="p.adjust") +
    scale_color_gradient(low='#E0F7FA', high='#4A148C') +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size=rel(1)),
          axis.text.x = element_text(size=rel(1)),
          plot.title = element_text(size=rel(1.2), face="bold"),
          legend.text = element_text(size=rel(1)),
          legend.title = element_text(size=rel(1)),
          axis.title.x = element_text(size=rel(1)),
          axis.title.y = element_text(size=rel(1))) +
    labs(title="KEGG Pathway Enrichment")
dev.off()

# 保存富集分析结果
write.csv(as.data.frame(go_bp), "GO_BP_results.csv")
write.csv(as.data.frame(go_cc), "GO_CC_results.csv")
write.csv(as.data.frame(go_mf), "GO_MF_results.csv")
write.csv(as.data.frame(kegg), "KEGG_results.csv")

# 输出统计信息
cat("\n富集分析统计信息：\n")
cat("GO生物过程富集项数量：", nrow(as.data.frame(go_bp)), "\n")
cat("GO细胞组成富集项数量：", nrow(as.data.frame(go_cc)), "\n")
cat("GO分子功能富集项数量：", nrow(as.data.frame(go_mf)), "\n")
cat("KEGG通路富集项数量：", nrow(as.data.frame(kegg)), "\n") 