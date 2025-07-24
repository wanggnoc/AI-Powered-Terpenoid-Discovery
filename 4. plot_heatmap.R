# 加载必要的包
if (!require("pheatmap")) install.packages("pheatmap")
if (!require("RColorBrewer")) install.packages("RColorBrewer")
library(pheatmap)
library(RColorBrewer)

# 读取数据
data <- read.csv("outAll.csv")

# 将数据重塑为宽格式
data_wide <- reshape(data, 
                    idvar = "Chem", 
                    timevar = "PDB", 
                    direction = "wide")

# 提取Bestaffinity列
affinity_matrix <- data_wide[, grep("Bestaffinity", names(data_wide))]
rownames(affinity_matrix) <- data_wide$Chem

# 重命名列名（去掉"Bestaffinity."前缀）
colnames(affinity_matrix) <- gsub("Bestaffinity.", "", colnames(affinity_matrix))

# 排除2g3r
affinity_matrix <- affinity_matrix[, !grepl("2g3r", colnames(affinity_matrix))]

# 替换蛋白质名称
protein_names <- c("2pjl" = "ESRRA",
                  "7sj3" = "CDK4",
                  "2or9" = "MYC",
                  "1g5m" = "BCL2",
                  "2enq" = "PIK3CA",
                  "6d3o" = "VEGFA")

colnames(affinity_matrix) <- protein_names[colnames(affinity_matrix)]

# 替换化合物名称
compound_names <- c("Paclitaxel" = "TKC200866",
                   "adenosine" = "TKC122482")

rownames(affinity_matrix) <- ifelse(rownames(affinity_matrix) %in% names(compound_names),
                                  compound_names[rownames(affinity_matrix)],
                                  rownames(affinity_matrix))

# 去掉Holacanthone行
affinity_matrix <- affinity_matrix[rownames(affinity_matrix) != "Holacanthone", ]
affinity_matrix2=affinity_matrix[c(-3,-42),]
# 创建马卡龙配色方案
macaron_colors <- colorRampPalette(c("#FF9999", "#FFCC99", "#FFFF99", "#99FF99", "#99CCFF"))(100)
macaron_colors <- colorRampPalette(c("#49006a", "#7a0177", "#ae017e",
                                  "#dd3497", "#f768a1", "#fa9fb5",
                                  "#fcc5c0", "#d1e5f0", "#92c5de",
                                  "#4393c3", "#2166ac", "#053061"))(100)

# 设置超出范围的颜色为灰色
affinity_matrix[affinity_matrix < -10] <- -10
affinity_matrix[affinity_matrix > 0] <- 0

# 绘制热图
pheatmap(affinity_matrix,
         color = macaron_colors,
         breaks = seq(-10, 0, length.out = 101),
         cluster_rows = FALSE,  # 关闭行聚类
         cluster_cols = FALSE,  # 关闭列聚类
         display_numbers = TRUE,
         fontsize_row = 14,     # 增大行标签字体
         fontsize_col = 14,     # 增大列标签字体
         fontsize = 14,         # 增大图例字体
         main = "Binding Affinity Heatmap",
         filename = "affinity_heatmap2.pdf",
         width = 12,
         height = 10) 


