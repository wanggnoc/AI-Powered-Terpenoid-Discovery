# 加载必要的包
library(VennDiagram)
library(showtext)

# 添加中文字体支持
font_add("SimSun", "SimSun.ttf")  # 添加宋体
showtext_auto()  # 自动使用showtext渲染

# 设置参数
venn.plot <- draw.pairwise.venn(
  area1 = 937,                # 萜类化合物靶标数量
  area2 = 813,                # 卵巢癌候选靶点数量
  cross.area = 109,           # 交集数量
  #category = c("萜类化合物靶标", "卵巢癌候选靶点"),
  fill = c("#A0D0D0", "#A7C0DF"),  # 马卡龙风格的填充颜色
  alpha = 0.6,                # 透明度
  lwd = 2,                    # 线条宽度
  cex = 1.5,                  # 文字大小
  cat.cex = 1.2,             # 类别标签大小
  cat.pos = c(0, 0),         # 类别标签位置
  cat.dist = 0.05,           # 类别标签距离
  #cat.fontfamily = "SimSun", # 使用宋体
  #fontfamily = "SimSun",     # 使用宋体
  #main = "萜类化合物靶标与卵巢癌候选靶点的交集分析",
  main.cex = 1.5,            # 标题大小
  #main.fontfamily = "SimSun" # 使用宋体
)

# 保存图片
pdf("venn_diagram2.pdf", width = 8, height = 8)
grid.draw(venn.plot)
dev.off()

# 同时保存PNG格式（通常PNG格式对中文支持更好）
png("venn_diagram.png", width = 800, height = 800, res = 100)
grid.draw(venn.plot)
dev.off()

# 打印统计信息
cat("韦恩图统计信息：\n")
cat("萜类化合物靶标数量：", 937, "\n")
cat("卵巢癌候选靶点数量：", 813, "\n")
cat("交集数量：", 109, "\n")
cat("仅在萜类化合物靶标中的数量：", 937 - 109, "\n")
cat("仅在卵巢癌候选靶点中的数量：", 813 - 109, "\n") 