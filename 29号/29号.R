# --- 依赖包 ---
# 运行此脚本前，请确保您已安装 'data.table' 和 'readxl' 包。
# 您可以使用以下命令在 R/RStudio 控制台中安装它们:
#
# install.packages("data.table")
# install.packages("readxl")
#

# 加载所需的库
# 我们将 'suppressPackageStartupMessages' 放在这里，以保持控制台输出整洁
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(readxl))

# --- 配置 ---
# 假设此脚本与 'TPM.gene.annotated.xls' 文件在同一个目录中
# 如果不在同一目录，请修改为文件的绝对路径
# 例如: "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/3.Quantity/TPM.gene.annotated.xls"
# 注意: 在 R 中，路径最好使用正斜杠 "/"
file_path <- "TPM.gene.annotated.xls"

# --- 脚本开始 ---

cat(paste("--- 正在探查文件:", file_path, "---\n"))

# 检查文件是否存在
if (!file.exists(file_path)) {
  # 使用 stop() 会在错误时停止脚本执行
  stop(paste("错误: 文件未找到:", file_path, "\n请确保文件路径正确，或者脚本与文件在同一目录下。"))
}

# 初始化一个变量来存储数据
df <- NULL
read_success <- FALSE

# 1. 尝试按 TSV (制表符分隔) 读取
# data.table::fread 非常快，并且能很好地自动检测分隔符，
# 但我们明确指定 sep="\t" 来优先尝试 TSV
cat("步骤 1: 尝试按 TSV (制表符分隔) 格式读取 (使用 data.table::fread)...\n")
tryCatch({
  # data.table = FALSE 使其返回一个标准的 data.frame，而不是 data.table
  df <- fread(file_path, sep = "\t", data.table = FALSE) 
  read_success <- TRUE
  cat("文件读取成功 (TSV 格式)！\n")
}, error = function(e_tsv) {
  # 如果出错，打印错误信息
  cat(paste("按 TSV 读取失败:", e_tsv$message, "\n"))
})

# 2. 如果 TSV 失败，尝试按 Excel (.xls) 格式读取
if (!read_success) {
  cat("\n步骤 2: 尝试按 Excel (.xls) 格式读取 (使用 readxl::read_excel)...\n")
  tryCatch({
    # read_excel 可以自动处理 .xls 和 .xlsx
    df <- read_excel(file_path)
    # read_excel 返回的是 tibble，我们将其转换为标准的 data.frame 以保持一致性
    df <- as.data.frame(df) 
    read_success <- TRUE
    cat("文件读取成功 (Excel 格式)！\n")
  }, error = function(e_xls) {
    cat(paste("按 Excel 读取也失败了:", e_xls$message, "\n"))
    stop("无法读取文件。请检查文件是否损坏，或是否已正确安装 'readxl' 包。")
  })
}

# 3. 打印文件基本信息
if (read_success && !is.null(df)) {
  cat("\n--- 1. 数据概览 (前 6 行) ---\n")
  # R 中的 head() 默认显示 6 行
  print(head(df))
  
  cat("\n--- 2. 数据维度 (行数, 列数) ---\n")
  # 使用 R 的 dim() 或 nrow()/ncol()
  cat(paste("Shape: [", nrow(df), ",", ncol(df), "]\n"))
  cat(paste("总共有", nrow(df), "个基因/转录本 (行)\n"))
  cat(paste("总共有", ncol(df), "个样本/注释 (列)\n"))
  
  cat("\n--- 3. 所有列的名称 ---\n")
  cat("列名列表:\n")
  print(colnames(df))
  
  cat("\n--- 4. 数据结构概览 (str) ---\n")
  # str() 提供了类似于 Python .info() 的列类型概览
  # 我们使用 capture.output 和 cat 来控制 str 的输出格式，防止它打印过多行
  cat(capture.output(str(df, list.len = 10)), sep = "\n")
  
  cat("\n--- 探查完毕 ---\n")
  cat("请根据上面的输出，判断哪些列是基因ID，哪些是注释，哪些是样本的TPM值。\n")
} else {
  cat("未能成功读取数据。\n")
}



# -------------------------------------------------------------------------
# 最终火山图绘制脚本 (V6 改编版 - 适应 KD vs NC 数据)
# -------------------------------------------------------------------------

# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 
# 加载必要的包
# install.packages(c("readr", "ggplot2", "ggrepel", "ggnewscale"))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggnewscale)) # 用于实现双色谱

# --- 2. 定义文件路径 ---
# (输入) 差异分析结果文件
input_file <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"

# (输出) 结果保存目录
output_dir <- "G:/XY/29/analysis"
# (输出) 图片文件名
output_png <- file.path(output_dir, "KD_vs_NC_volcano_plot.png")
output_pdf <- file.path(output_dir, "KD_vs_NC_volcano_plot.pdf")

# 自动创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste("已创建目录:", output_dir, "\n"))
}

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取文件:", input_file, "\n"))
# 使用 read_tsv (因为它是 .xls 后缀的制表符分隔文件)
data <- read_tsv(input_file, col_types = cols(.default = "c")) 
cat("文件读取完毕。\n")

# --- 3.1 数据清洗和类型转换 ---
# 转换为标准 data.frame
data <- as.data.frame(data)
# 使用 'gene' 列作为行名 (如果需要的话，但后续我们直接用列名更安全)
# rownames(data) <- data$gene 

# 关键：将我们需要的列从 'character' 转换为 'numeric'
# 我们在读取时指定了所有列为 'c' (character)，现在手动转换
cols_to_numeric <- c("logFC.KD-vs-NC", "PValue.KD-vs-NC", "FDR.KD-vs-NC")
for (col in cols_to_numeric) {
  data[, col] <- as.numeric(data[, col])
}

# --- 3.2 列名映射 (改为 R 风格的列，更易于后续处理) ---
data$log2FoldChange <- data$`logFC.KD-vs-NC`
data$pvalue_raw     <- data$`PValue.KD-vs-NC`
data$pvalue_adj     <- data$`FDR.KD-vs-NC`

# 过滤掉NA值
data <- data[!is.na(data$log2FoldChange) & !is.na(data$pvalue_raw) & !is.na(data$pvalue_adj), ]

# 计算 Y 轴 (-log10 原始P值)
data$neg_log10_pvalue <- -log10(data$pvalue_raw)
# 处理 pvalue=0 导致的 Infinite 值
infinite_rows <- is.infinite(data$neg_log10_pvalue)
if (any(infinite_rows)) {
  cat(paste("警告：有", sum(infinite_rows), "个基因的 PValue 为 0。\n"))
  # 将 Infinite 值替换为一个比最大非 Infinite 值稍大的数
  max_finite_p <- max(data$neg_log10_pvalue[!infinite_rows], na.rm = TRUE)
  data$neg_log10_pvalue[infinite_rows] <- max_finite_p * 1.1 
  cat("已将 Infinite 值替换为:", max_finite_p * 1.1, "\n")
}
plot_data <- data

# --- 4. 定义阈值并拆分数据 ---
logfc_threshold <- 2  # log2FC 阈值 (代表 4 倍差异) - 已根据您的要求从 1 调整
fdr_threshold   <- 0.01 # 显著性阈值 (使用 FDR) - 已根据您的要求从 0.05 调整

# 筛选:
# 1. 不显著
data_ns <- subset(plot_data, abs(log2FoldChange) <= logfc_threshold | pvalue_adj >= fdr_threshold)
# 2. 显著
data_sig <- subset(plot_data, abs(log2FoldChange) > logfc_threshold & pvalue_adj < fdr_threshold)
# 3. 显著上调
data_up <- subset(data_sig, log2FoldChange > 0)
# 4. 显著下调
data_down <- subset(data_sig, log2FoldChange < 0)

# 存储统计数字用于标注 (在过滤前计算总数)
n_up <- nrow(data_up)
n_down <- nrow(data_down)

cat(paste("筛选完毕：总计", (n_up + n_down), "个差异基因\n"))
cat(paste("上调:", n_up, "个\n"))
cat(paste("下调:", n_down, "个\n"))

# --- 4.1 过滤极高P值的点 (根据用户要求) ---
# 移除 Y 轴 > 320 的点，以便在图上不再显示它们
y_limit_display <- 320
plot_data <- subset(plot_data, neg_log10_pvalue <= y_limit_display)
data_ns     <- subset(data_ns,   neg_log10_pvalue <= y_limit_display)
data_up     <- subset(data_up,   neg_log10_pvalue <= y_limit_display)
data_down   <- subset(data_down, neg_log10_pvalue <= y_limit_display)
cat(paste("已过滤掉 neg_log10(PValue) >", y_limit_display, "的点，不再用于绘图。\n"))


# --- 5. 基因标签处理 (已移除) ---
# 根据您的要求，本节已删除，不再生成基因标签。

# --- 6. 开始绘图 (使用您的 V6 模板) ---
cat("开始绘制火山图...\n")
p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  
  # 第1层: 不显著的点 (灰色)
  geom_point(data = data_ns, color = "grey", size = 1.2, alpha = 0.6, stroke = 0) +
  
  # 第2层: 下调的显著点 (蓝色系渐变)
  geom_point(data = data_down, aes(color = neg_log10_pvalue, size = neg_log10_pvalue), alpha = 0.8, stroke = 0) +
  scale_color_gradient(low = "#b5dbe6", high = "#333aab", name = "-log10(p-value) Down") +
  
  # 开始一个新的颜色标度
  new_scale_color() +
  
  # 第3层: 上调的显著点 (红色系渐变)
  geom_point(data = data_up, aes(color = neg_log10_pvalue, size = neg_log10_pvalue), alpha = 0.8, stroke = 0) +
  scale_color_gradient(low = "#fe8264", high = "#a73336", name = "-log10(p-value) Up") +

  # --- 其他设置 ---
  # 阈值线
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "#999999") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "#999999") + # 注意: 水平线基于 FDR 阈值
  
  scale_size_continuous(range = c(1, 4)) +
  scale_y_continuous(expand = c(0, 0)) +  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  guides(size = "none") + # 隐藏大小图例
  
  # --- 6.1 添加统计标注 (根据用户要求) ---
  annotate("text",
           x = 18, # X 坐标 (右上角)
           y = 300, # Y 坐标
           label = paste("Up:", n_up),
           size = 5, color = "#a73336", fontface = "bold") +
  annotate("text",
           x = -18, # X 坐标 (左上角)
           y = 300, # Y 坐标
           label = paste("Down:", n_down),
           size = 5, color = "#333aab", fontface = "bold") +
  
  # 基因标签 (已移除)
  # geom_text_repel(...) 已根据您的要求删除
  
  # 标题和坐标轴
  xlab("Log2 Fold Change (KD vs NC)") +
  ylab("-log10(P-VALUE)") +
  ggtitle("Volcano Plot: KD vs NC") +
  
  # 坐标轴范围 (已根据您的要求添加)
  # 设置 X 轴范围为 -20 到 20
  # 设置 Y 轴范围为 0 到 320
  coord_cartesian(xlim = c(-20, 20), ylim = c(0, 320))

# --- 7. 显示并保存图片 ---
cat("绘图完毕，正在保存...\n")

# 保存为 PNG
ggsave(output_png, plot = p, width = 10, height = 8, units = "in", dpi = 300)
# 保存为 PDF (矢量图)
ggsave(output_pdf, plot = p, width = 10, height = 8, units = "in")

cat(paste("成功保存图片到:\n", output_png, "\n", output_pdf, "\n"))

# 在 RStudio Viewer 中显示图片
print(p)


# -------------------------------------------------------------------------
# 差异基因热图 (Heatmap) 脚本
# -------------------------------------------------------------------------

# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# 检查和安装核心包
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(RColorBrewer))

# --- 2. 定义文件路径 ---
# (输入) 差异分析结果文件
input_file <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"

# (输出) 结果保存目录
output_dir <- "G:/XY/29/analysis"
# (输出) 图片文件名
output_png <- file.path(output_dir, "KD_vs_NC_DEG_Heatmap.png")
output_pdf <- file.path(output_dir, "KD_vs_NC_DEG_Heatmap.pdf")

# 自动创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste("已创建目录:", output_dir, "\n"))
}

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取文件:", input_file, "\n"))
data <- read_tsv(input_file, col_types = cols(.default = "c")) 
cat("文件读取完毕。\n")

# 转换为标准 data.frame
data <- as.data.frame(data)

# 关键：将我们需要的列从 'character' 转换为 'numeric'
# 我们只需要差异统计列 和 6个样本的TPM列
cols_to_numeric <- c("logFC.KD-vs-NC", "PValue.KD-vs-NC", "FDR.KD-vs-NC", 
                     "con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
for (col in cols_to_numeric) {
  data[, col] <- as.numeric(data[, col])
}

# 列名映射 (用于筛选)
data$log2FoldChange <- data$`logFC.KD-vs-NC`
data$pvalue_adj     <- data$`FDR.KD-vs-NC`

# 过滤掉NA值
data <- data[!is.na(data$log2FoldChange) & !is.na(data$pvalue_adj), ]

# --- 4. 筛选差异基因 (使用与火山图相同的标准) ---
logfc_threshold <- 2    # log2FC 阈值 (4倍)
fdr_threshold   <- 0.01 # 显著性阈值

# 筛选显著差异的基因
data_sig <- subset(data, abs(log2FoldChange) > logfc_threshold & pvalue_adj < fdr_threshold)

cat(paste("根据 FDR <", fdr_threshold, "且 |log2FC| >", logfc_threshold, "标准\n"))
cat(paste("共筛选到", nrow(data_sig), "个差异基因用于绘制热图。\n"))

# --- 5. 准备表达矩阵 ---
# 提取6个样本的TPM值
sample_cols <- c("con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
heatmap_matrix <- data_sig[, sample_cols]

# 检查是否有 'SYMBOL' 列，并用它作为行名
if ("SYMBOL" %in% colnames(data_sig)) {
  # 检查是否有重复的 SYMBOL，这在热图中很重要
  if (any(duplicated(data_sig$SYMBOL))) {
    cat("警告：存在重复的基因 'SYMBOL'。将使用 make.unique() 处理。\n")
    rownames(heatmap_matrix) <- make.unique(data_sig$SYMBOL)
  } else {
    rownames(heatmap_matrix) <- data_sig$SYMBOL
  }
} else {
  # 如果没有 SYMBOL，使用 'gene' 列
  rownames(heatmap_matrix) <- data_sig$gene
}

# 关键步骤 1: Log2 转换 (处理 0 值并缩小范围)
# log2(TPM + 1)
heatmap_matrix_log2 <- log2(heatmap_matrix + 1)

# --- 6. 创建样本 (列) 注释 ---
annotation_col <- data.frame(
  Group = factor(c(rep("con", 3), rep("siha", 3)))
)
rownames(annotation_col) <- sample_cols

# 定义分组的颜色
ann_colors <- list(
  Group = c(con = "#333aab", siha = "#a73336") # 沿用火山图的蓝/红
)

# --- 7. "Nature" 风格颜色定义 ---
# 创建一个 蓝-白-红 的颜色梯度
# Z-score标准化后，负值=蓝色(下调)，0=白色(均值)，正值=红色(上调)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
# 或者使用 "navy", "white", "firebrick3"
# heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)


# --- 8. 绘制热图 ---
cat("开始绘制热图并保存...\n")

p_heatmap <- pheatmap(
  heatmap_matrix_log2,         # 输入: log2(TPM+1) 矩阵
  color = heatmap_colors,      # 颜色: 蓝-白-红 梯度
  
  # 关键步骤 2: Z-score 标准化 (非常重要)
  scale = "row",               # 对每一行(基因)进行 Z-score 标准化
  
  # 聚类
  cluster_rows = TRUE,         # 聚类行 (基因)
  cluster_cols = TRUE,         # 聚类列 (样本)
  clustering_distance_rows = "euclidean", # 聚类方法
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  
  # 样本注释
  annotation_col = annotation_col, # 添加列注释 (分组)
  annotation_colors = ann_colors,  # 列注释的颜色
  
  # 外观
  show_rownames = FALSE,       # 不显示行名 (基因太多)
  show_colnames = TRUE,        # 显示列名 (样本名)
  border_color = NA,           # 去除单元格边框
  main = paste0("Heatmap of ", nrow(data_sig), " DEGs (FDR<", fdr_threshold, ", |log2FC|>", logfc_threshold, ")"),
  fontsize = 10,
  
  # 保存到文件
  filename = output_png,       # 保存为 PNG
  width = 8,
  height = 10
)

# 再次绘制并保存为 PDF (矢量图)
pheatmap(
  heatmap_matrix_log2,
  color = heatmap_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  border_color = NA,
  main = paste0("Heatmap of ", nrow(data_sig), " DEGs (FDR<", fdr_threshold, ", |log2FC|>", logfc_threshold, ")"),
  fontsize = 10,
  filename = output_pdf, # 保存为 PDF
  width = 8,
  height = 10
)

cat(paste("成功保存热图到:\n", output_png, "\n", output_pdf, "\n"))


# -------------------------------------------------------------------------
# 差异基因相关性网络 脚本 (V5 - 仅可视化 Top 9 亚群 + 保存基因列表)
# -------------------------------------------------------------------------

# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# 检查和安装核心包
packages_to_check <- c("Hmisc", "igraph", "ggraph", "readr", "dplyr", "RColorBrewer", "tidyr")
for (pkg in packages_to_check) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# --- 2. 定义文件路径 ---
# (输入) 差异分析结果文件
input_file <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"

# (输出) 结果保存目录 (新：网络分析专用文件夹)
output_dir <- "G:/XY/29/analysis/3_Network"
# (新) 亚群基因列表的保存目录
output_modules_dir <- file.path(output_dir, "Top9_Modules_GeneLists")

# (输出) 文件名
output_png <- file.path(output_dir, "KD_vs_NC_DEG_CoNetwork_Top9_Clusters.png")
output_pdf <- file.path(output_dir, "KD_vs_NC_DEG_CoNetwork_Top9_Clusters.pdf")
output_edges_csv <- file.path(output_dir, "KD_vs_NC_DEG_CoNetwork_ALL_edges.csv")
output_nodes_csv <- file.path(output_dir, "KD_vs_NC_DEG_CoNetwork_ALL_nodes.csv")

# 自动创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste("已创建目录:", output_dir, "\n"))
}
# (新) 自动创建亚群列表目录
if (!dir.exists(output_modules_dir)) {
  dir.create(output_modules_dir, recursive = TRUE)
  cat(paste("已创建目录:", output_modules_dir, "\n"))
}

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取文件:", input_file, "\n"))
data <- read_tsv(input_file, col_types = cols(.default = "c")) 
data <- as.data.frame(data)

cols_to_numeric <- c("logFC.KD-vs-NC", "PValue.KD-vs-NC", "FDR.KD-vs-NC", 
                     "con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
for (col in cols_to_numeric) {
  data[, col] <- as.numeric(data[, col])
}

data$log2FoldChange <- data$`logFC.KD-vs-NC`
data$pvalue_adj     <- data$`FDR.KD-vs-NC`
data <- data[!is.na(data$log2FoldChange) & !is.na(data$pvalue_adj), ]

# --- 4. 筛选差异基因 ---
logfc_threshold <- 2    # log2FC 阈值 (4倍)
fdr_threshold   <- 0.01 # 显著性阈值

data_sig <- subset(data, abs(log2FoldChange) > logfc_threshold & pvalue_adj < fdr_threshold)
cat(paste("共筛选到", nrow(data_sig), "个差异基因用于构建网络。\n"))

# --- 5. 准备表达矩阵 ---
sample_cols <- c("con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
heatmap_matrix <- data_sig[, sample_cols]
# (改进) 检查 'SYMBOL' 列是否有效 (非NA, 非"-")
if ("SYMBOL" %in% colnames(data_sig) && !all(is.na(data_sig[,"SYMBOL"])) && !all(data_sig$SYMBOL == "-")) {
  gene_names <- data_sig$SYMBOL
  cat("使用 'SYMBOL' 作为基因名。\n")
} else {
  gene_names <- data_sig$gene
  cat("警告: 'SYMBOL' 列无效, 将使用 'gene' (GeneID) 作为基因名。\n")
}
gene_names <- make.unique(gene_names)
rownames(heatmap_matrix) <- gene_names
data_sig$gene_name_unique <- gene_names
heatmap_matrix_log2 <- log2(heatmap_matrix + 1)

# --- 6. 计算相关性矩阵 ---
cat("正在计算基因间的 Pearson 相关性...\n")
cor_matrix_input <- t(heatmap_matrix_log2)
cor_results <- rcorr(cor_matrix_input, type = "pearson")
cor_r_matrix <- cor_results$r
cor_p_matrix <- cor_results$P

# --- 7. 扁平化矩阵并过滤 ---
cat("正在构建边列表 (Edge List)...\n")
edges_r <- as.data.frame(as.table(cor_r_matrix), stringsAsFactors = FALSE)
names(edges_r) <- c("from", "to", "correlation")
edges_p <- as.data.frame(as.table(cor_p_matrix), stringsAsFactors = FALSE)
names(edges_p) <- c("from_p", "to_p", "pvalue")
edges_full <- cbind(edges_r, pvalue = edges_p$pvalue)

edges_filtered <- edges_full %>%
  filter(as.character(from) < as.character(to)) %>%
  filter(!is.na(correlation) & !is.na(pvalue))

# --- 8. 定义阈值并筛选边 ---
cat("正在对 P-values 进行多重假设检验校正 (FDR)...\n")
edges_filtered$fdr <- p.adjust(edges_filtered$pvalue, method = "BH")

corr_threshold <- 0.95 # |相关系数| 必须大于 0.95
fdr_threshold  <- 0.05 # FDR 必须小于 0.05

edges_final <- edges_filtered %>%
  filter(abs(correlation) > corr_threshold & fdr < fdr_threshold) %>%
  mutate(
    direction = ifelse(correlation > 0, "Positive", "Negative"),
    direction = factor(direction, levels = c("Positive", "Negative")) 
  )

cat(paste("在 |r| >", corr_threshold, "且 FDR <", fdr_threshold, "的阈值下，共找到", nrow(edges_final), "条显著连接。\n"))

# --- 8.2 安全检查 (如果边仍然过多) ---
max_edges_for_plotting <- 3000
if(nrow(edges_final) > max_edges_for_plotting) {
  cat(paste("警告：筛选后的边 (", nrow(edges_final), ") 仍然过多，不利于可视化。\n"))
  cat(paste("将只提取相关性最强 (Top", max_edges_for_plotting, ") 的边进行绘图。\n"))
  
  edges_final <- edges_final %>%
    arrange(desc(abs(correlation))) %>%
    head(max_edges_for_plotting)
  
  cat(paste("已缩减至", nrow(edges_final), "条边。\n"))
}

if(nrow(edges_final) == 0) {
  cat("警告：未找到符合阈值的边。网络为空。\n")
  stop("无法继续，网络为空。")
}

# --- 9. 准备节点列表 (Node List) ---
# (已修复) 移除 select(...) 中的 'regulation = stat.KD-vs-NC' 来避免列名错误
nodes_final <- data_sig %>%
  select(gene_name_unique, log2FoldChange) %>% 
  mutate(regulation = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
  rename(name = gene_name_unique) 

# --- 10. 创建 igraph 网络对象 ---
graph_full <- graph_from_data_frame(d = edges_final, vertices = nodes_final, directed = FALSE)
isolated_nodes <- V(graph_full)[degree(graph_full) == 0]
graph_final <- delete_vertices(graph_full, isolated_nodes)

cat(paste("最终网络包含", vcount(graph_final), "个基因 (节点) 和", ecount(graph_final), "条连接 (边)。\n"))

# --- 10.1 社区检测 (亚群划分) ---
cat("正在进行 Louvain 社区检测...\n")
set.seed(123)
community_louvain <- cluster_louvain(graph_final)
# 将社区(亚群)信息添加回节点
V(graph_final)$community <- as.factor(community_louvain$membership)
nodes_df <- as.data.frame(vertex_attr(graph_final)) # 获取包含 community 的完整节点列表

cat(paste("已将网络划分为", length(unique(community_louvain$membership)), "个亚群。\n"))

# --- 10.2 筛选 Top 9 亚群 ---
cat("正在筛选 Top 9 最大的亚群...\n")

community_sizes <- sizes(community_louvain)
community_summary <- data.frame(
  community = names(community_sizes),
  size = as.numeric(community_sizes)
) %>%
  arrange(desc(size))

top_9_communities <- head(community_summary, 9)$community

# (新) 打印总结表格到控制台
cat("--- Top 9 亚群总结 ---\n")
top_9_summary <- nodes_df %>%
  filter(community %in% top_9_communities) %>%
  group_by(community) %>%
  summarise(
    total_genes = n(),
    up_genes = sum(regulation == "Up"),
    down_genes = sum(regulation == "Down")
  ) %>%
  # 重新排序，以匹配 community_summary 的顺序
  mutate(community = factor(community, levels = top_9_communities)) %>%
  arrange(community)

print(top_9_summary)
cat("-------------------------\n")


# (新) 从图中删除所有 "非 Top 9" 的节点
nodes_to_remove <- V(graph_final)[!(V(graph_final)$community %in% top_9_communities)]
graph_top9 <- delete_vertices(graph_final, nodes_to_remove)

cat(paste("已筛选 Top 9 亚群，剩余", vcount(graph_top9), "个基因和", ecount(graph_top9), "条连接用于绘图。\n"))

# --- 10.3 (新) 保存 Top 9 亚群的基因列表 ---
cat(paste("正在保存 Top 9 亚群的基因列表到:", output_modules_dir, "\n"))
for (comm_id in top_9_communities) {
  # 筛选出该亚群的基因
  module_genes <- nodes_df %>%
    filter(community == comm_id) %>%
    select(GeneSymbol = name, log2FoldChange, Regulation = regulation) %>%
    arrange(desc(log2FoldChange)) # 按 logFC 排序
  
  # 定义文件名
  file_name <- file.path(output_modules_dir, paste0("Module_", comm_id, "_genes (n=", nrow(module_genes), ").csv"))
  
  # 保存
  write.csv(module_genes, file_name, row.names = FALSE)
}
cat("所有亚群基因列表保存完毕。\n")


# --- 11. 可视化网络 (新：仅 Top 9 亚群) ---
cat("开始绘制 Top 9 亚群分面网络图...\n")

color_up   <- "#a73336"
color_down <- "#333aab"

set.seed(123)
gg_network <- ggraph(graph_top9, layout = 'fr') + 
  
  # 1. 绘制边 (连接)
  geom_edge_link(aes(color = direction, width = abs(correlation)), alpha = 0.5) +
  scale_edge_color_manual(values = c("Positive" = "red", "Negative" = "blue"), name = "Correlation") +
  scale_edge_width(range = c(0.2, 1.5), name = "|Correlation|") +
  
  # 2. 绘制节点 (基因)
  geom_node_point(aes(color = regulation, size = abs(log2FoldChange)), alpha = 0.8) +
  scale_color_manual(values = c("Up" = color_up, "Down" = color_down), name = "Regulation") + 
  scale_size(range = c(2, 8), name = "|log2FC|") +
  
  # 3. 按社区(亚群)分面 - 现在只会显示 9 个窗格
  facet_nodes(~community, scales = "free", ncol = 3) + # 指定为 3x3 网格
  
  # 4. 主题 (移除基因标签)
  theme_graph(base_family = 'sans', background = "white", border = TRUE) +
  ggtitle(paste0("Top 9 DEG Co-expression Modules (", vcount(graph_top9), " Nodes, ", ecount(graph_top9), " Edges)"),
          subtitle = paste0("Threshold: |r| > ", corr_threshold, " & FDR < ", fdr_threshold, " | Clustered by Louvain")) +
  theme(
    legend.position = "right",
    # 增加分面标签的字体大小
    strip.text = element_text(size = 14, face = "bold")
  )

# --- 12. 保存结果 ---
cat("绘图完毕，正在保存...\n")

# 保存图片
ggsave(output_png, plot = gg_network, width = 14, height = 12, units = "in", dpi = 300)
ggsave(output_pdf, plot = gg_network, width = 14, height = 12, units = "in")

# 保存数据 (我们仍然保存完整的(未过滤的)数据，以便在 Cytoscape 中分析)
nodes_to_save <- as.data.frame(vertex_attr(graph_final)) # 保存所有节点信息
write.csv(nodes_to_save, output_nodes_csv, row.names = FALSE)
edges_to_save <- as.data.frame(as_edgelist(graph_final)) # 保存所有边信息
names(edges_to_save) <- c("from", "to")
edges_to_save <- cbind(edges_to_save, edge_attr(graph_final)) 
write.csv(edges_to_save, output_edges_csv, row.names = FALSE)

cat(paste("成功保存网络图和数据文件到:\n", output_dir, "\n"))
print(gg_network)
# -------------------------------------------------------------------------
# 导出差异基因列表 (用于 GO/KEGG 平台) 脚本
# -------------------------------------------------------------------------

# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# 检查和安装核心包
packages_to_check <- c("readr", "dplyr")
for (pkg in packages_to_check) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# --- 2. 定义文件路径 ---
# (输入) 差异分析结果文件
input_file <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"

# (输出) 结果保存目录 (新：富集分析列表专用文件夹)
output_dir <- "G:/XY/29/analysis/4_Enrichment_Lists"

# (输出) 文件名
output_all_degs <- file.path(output_dir, "KD_vs_NC_DEGs_All.csv")
output_up_degs  <- file.path(output_dir, "KD_vs_NC_DEGs_Up.csv")
output_down_degs <- file.path(output_dir, "KD_vs_NC_DEGs_Down.csv")

# 自动创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste("已创建目录:", output_dir, "\n"))
}

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取文件:", input_file, "\n"))
data <- read_tsv(input_file, col_types = cols(.default = "c")) 
data <- as.data.frame(data)

cols_to_numeric <- c("logFC.KD-vs-NC", "PValue.KD-vs-NC", "FDR.KD-vs-NC")
for (col in cols_to_numeric) {
  data[, col] <- as.numeric(data[, col])
}

data$log2FoldChange <- data$`logFC.KD-vs-NC`
data$pvalue_adj     <- data$`FDR.KD-vs-NC`
data <- data[!is.na(data$log2FoldChange) & !is.na(data$pvalue_adj), ]

# --- 4. 筛选差异基因 (使用我们统一的标准) ---
logfc_threshold <- 2    # log2FC 阈值 (4倍)
fdr_threshold   <- 0.01 # 显著性阈值

data_sig <- subset(data, abs(log2FoldChange) > logfc_threshold & pvalue_adj < fdr_threshold)
data_up  <- subset(data_sig, log2FoldChange > 0)
data_down <- subset(data_sig, log2FoldChange < 0)

cat(paste("共筛选到", nrow(data_sig), "个差异基因。\n"))
cat(paste("上调:", nrow(data_up), "个\n"))
cat(paste("下调:", nrow(data_down), "个\n"))

# --- 5. 准备基因列表 (使用 SYMBOL 或 gene ID) ---
# (改进) 检查 'SYMBOL' 列是否有效 (非NA, 非"-")
# 平台分析需要唯一的、有效的基因标识符
get_gene_list <- function(df) {
  if ("SYMBOL" %in% colnames(df) && !all(is.na(df[,"SYMBOL"])) && !all(df$SYMBOL == "-")) {
    cat("正在提取 'SYMBOL' 作为基因名。\n")
    gene_list <- df$SYMBOL
  } else {
    cat("警告: 'SYMBOL' 列无效, 将使用 'gene' (GeneID) 作为基因名。\n")
    gene_list <- df$gene
  }
  # 确保基因名是唯一的
  gene_list <- unique(gene_list)
  # 转换为数据框以便保存
  return(data.frame(GeneSymbol = gene_list))
}

# 提取列表
list_all_degs <- get_gene_list(data_sig)
list_up_degs  <- get_gene_list(data_up)
list_down_degs <- get_gene_list(data_down)

# --- 6. 保存到文件 ---
# 我们将文件保存为 .csv, 包含一个列头 "GeneSymbol"
# 很多平台也接受 .txt 文件，您可以根据需要修改后缀
cat(paste("正在保存基因列表到:", output_dir, "\n"))

write.csv(list_all_degs, file = output_all_degs, row.names = FALSE, quote = FALSE)
write.csv(list_up_degs,  file = output_up_degs,  row.names = FALSE, quote = FALSE)
write.csv(list_down_degs, file = output_down_degs, row.names = FALSE, quote = FALSE)

cat("--- 所有基因列表保存完毕 ---\n")
cat(paste("总差异基因 (", nrow(list_all_degs), "):", output_all_degs, "\n"))
cat(paste("上调基因 (", nrow(list_up_degs), "):", output_up_degs, "\n"))
cat(paste("下调基因 (", nrow(list_down_degs), "):", output_down_degs, "\n"))

# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# 检查和安装核心包
packages_to_check <- c("readr", "dplyr", "stringr")
for (pkg in packages_to_check) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr)) # 用于更方便的文本搜索

# --- 2. 定义文件路径 ---
# (输入) 差异分析结果文件
input_file <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"

# (输出) 新的结果保存目录
output_dir <- "G:/XY/29/analysis/5_Specific_Gene_Analysis"
# (输出) 筛选出的特定基因列表文件
output_specific_genes_csv <- file.path(output_dir, "Specific_DEGs_Inflammation_Macrophage.csv")

# 自动创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste("已创建新目录:", output_dir, "\n"))
}

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取文件:", input_file, "\n"))
data <- read_tsv(input_file, col_types = cols(.default = "c")) 
data <- as.data.frame(data)

# 转换数值列 (与您之前的脚本一致)
cols_to_numeric <- c("logFC.KD-vs-NC", "PValue.KD-vs-NC", "FDR.KD-vs-NC", 
                     "con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
# 检查列是否存在
cols_exist <- cols_to_numeric %in% colnames(data)
if (!all(cols_exist)) {
  cat("警告：以下列在文件中不存在：\n", paste(cols_to_numeric[!cols_exist], collapse = "\n"), "\n")
}
# 仅转换存在的列
for (col in cols_to_numeric[cols_exist]) {
  data[, col] <- as.numeric(data[, col])
}

# 列名映射
data$log2FoldChange <- data$`logFC.KD-vs-NC`
data$pvalue_adj     <- data$`FDR.KD-vs-NC`

# 过滤NA
data <- data[!is.na(data$log2FoldChange) & !is.na(data$pvalue_adj), ]
cat("数据读取和预处理完毕。\n")

# --- 4. 筛选差异基因 (DEG) ---
logfc_threshold <- 2    # log2FC 阈值 (4倍)
fdr_threshold   <- 0.01 # 显著性阈值

data_sig <- subset(data, abs(log2FoldChange) > logfc_threshold & pvalue_adj < fdr_threshold)
cat(paste("共筛选到", nrow(data_sig), "个差异基因 (DEGs)。\n"))

# --- 5. 筛选特定功能的基因 ---
cat("正在从 DEGs 中筛选特定功能基因...\n")

# 定义搜索关键词 (不区分大小写)
search_terms <- c(
  "inflammatory",   # 炎症
  "inflammation",   # 炎症
  "cytokine",       # 细胞因子
  "chemokine",      # 趋化因子
  "interleukin",    # 白介素
  "macrophage",     # 巨噬细胞
  "polarization"    # 极化
)

# 构建搜索模式 (用 | 连接所有关键词)
search_pattern <- paste(search_terms, collapse = "|")
cat(paste("使用的搜索关键词:", search_pattern, "\n"))

# *** 已修正 ***
# 确定要在哪些列中搜索 (根据您提供的列名)
potential_cols <- c("DESCRIPTION", "GO", "PATH")
search_cols <- intersect(potential_cols, colnames(data_sig))

if (length(search_cols) == 0) {
  stop(paste("错误：在文件中未找到任何预期的注释列 (如", paste(potential_cols, collapse = ", "), ")。\n这不应该发生，请检查原始文件。"))
}
cat(paste("将在以下列中搜索:", paste(search_cols, collapse = ", "), "\n"))

# 创建一个逻辑向量 (matrix)，标记每一行是否在任一列中匹配
match_matrix <- sapply(data_sig[, search_cols, drop = FALSE], function(col) {
  str_detect(col, regex(search_pattern, ignore_case = TRUE))
})
# 处理NA值 (str_detect 对 NA 返回 NA，我们将其视为 FALSE)
match_matrix[is.na(match_matrix)] <- FALSE

# 只要任何一列匹配 (rowSums > 0)，就保留该行
keep_rows <- rowSums(match_matrix) > 0

data_specific_genes <- data_sig[keep_rows, ]

cat(paste("--- 筛选完毕 --- \n共找到", nrow(data_specific_genes), "个与关键词匹配的差异基因。\n"))

# --- 6. 保存筛选结果 ---
if (nrow(data_specific_genes) > 0) {
  # 排序 (按 logFC 绝对值降序)
  data_specific_genes <- data_specific_genes %>%
    arrange(desc(abs(log2FoldChange)))
  
  # 保存到 CSV
  write.csv(data_specific_genes, file = output_specific_genes_csv, row.names = FALSE)
  
  cat(paste("成功保存特定基因列表到:\n", output_specific_genes_csv, "\n"))
  
  cat("\n--- 文件内容预览 (前 5 个基因的 SYMBOL 和 DESCRIPTION) ---\n")
  # 优先显示 SYMBOL，如果不存在则显示 gene
  gene_col_name <- ifelse("SYMBOL" %in% colnames(data_specific_genes), "SYMBOL", "gene")
  annot_col_name <- "DESCRIPTION"
  
  preview_df <- data_specific_genes[1:min(5, nrow(data_specific_genes)), c(gene_col_name, "log2FoldChange", annot_col_name)]
  print(preview_df)
  
} else {
  cat("警告：未筛选到任何符合条件的基因。无法保存文件。\n")
  cat("请检查您的关键词 (search_terms)。\n")
}

# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# 检查和安装核心包
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(readr))       
suppressPackageStartupMessages(library(RColorBrewer)) 

# --- 2. 定义文件路径 ---
# (输入) 上一步筛选出的特定基因列表
input_csv <- "G:/XY/29/analysis/5_Specific_Gene_Analysis/Specific_DEGs_Inflammation_Macrophage.csv"

# (输出) 结果保存目录
output_dir <- "G:/XY/29/analysis/5_Specific_Gene_Analysis"
# (输出) 图片文件名
output_png <- file.path(output_dir, "Top30_Specific_DEGs_Heatmap.png")
output_pdf <- file.path(output_dir, "Top30_Specific_DEGs_Heatmap.pdf")

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取特定基因文件:", input_csv, "\n"))
# 检查文件是否存在
if (!file.exists(input_csv)) {
  stop(paste("错误: 输入文件未找到:", input_csv, "\n请确保上一步已成功运行。"))
}

# 使用 read_csv 读取
data_specific <- read_csv(input_csv, col_types = cols(.default = "c")) 

# *** 关键修正 ***
# 将 read_csv 读入的 'tibble' 转换回 'data.frame'，以解决 'list' 错误
data_specific <- as.data.frame(data_specific)
# *** 修正完毕 ***

cat(paste("成功读取", nrow(data_specific), "个基因。\n"))

# 转换数值列
sample_cols <- c("con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
for (col in sample_cols) {
  data_specific[, col] <- as.numeric(data_specific[, col])
}
cat("样本TPM值已成功转换为数值。\n")

# --- 4. 筛选 Top 30 并准备矩阵 ---
n_top <- 30
data_top30 <- head(data_specific, n_top)
cat(paste("已提取 Top", n_top, "个基因用于绘图。\n"))

# 提取表达矩阵 (TPM 值)
heatmap_matrix_raw <- data_top30[, sample_cols]

# --- 4.1 设置行名 (使用 SYMBOL) ---
if ("SYMBOL" %in% colnames(data_top30) && !all(is.na(data_top30$SYMBOL))) {
  gene_names <- data_top30$SYMBOL
  cat("使用 'SYMBOL' 作为热图的行名。\n")
} else {
  gene_names <- data_top30$gene 
  cat("警告: 'SYMBOL' 列无效, 将使用 'gene' (GeneID) 作为行名。\n")
}

# 检查是否有重复的基因名，并处理
if (any(duplicated(gene_names))) {
  cat("警告：存在重复的基因名，将使用 make.unique() 处理。\n")
  rownames(heatmap_matrix_raw) <- make.unique(gene_names)
} else {
  rownames(heatmap_matrix_raw) <- gene_names
}

# --- 4.2 关键转换: Log2(TPM + 1) ---
heatmap_matrix_log2 <- log2(heatmap_matrix_raw + 1)

# --- 5. "Nature" 风格设置 (注释和颜色) ---
# 5.1 样本 (列) 注释
annotation_col <- data.frame(
  Group = factor(c(rep("con", 3), rep("siha", 3)))
)
rownames(annotation_col) <- sample_cols

# 5.2 分组颜色 (沿用您火山图的蓝/红)
ann_colors <- list(
  Group = c(con = "#333aab", siha = "#a73336") 
)

# 5.3 热图颜色 (蓝-白-红 梯度)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# --- 6. 绘制热图 ---
cat("开始绘制热图并保存...\n")

row_fontsize <- 8 

p_heatmap <- pheatmap(
  heatmap_matrix_log2,         
  color = heatmap_colors,      
  scale = "row",               
  cluster_rows = TRUE,         
  cluster_cols = TRUE,         
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,  
  show_rownames = TRUE,        
  show_colnames = TRUE,        
  border_color = NA,           
  main = paste0("Top ", n_top, " Inflammation & Macrophage DEGs Heatmap"),
  fontsize = 10,
  fontsize_row = row_fontsize, 
  filename = output_png,       
  width = 8,                   
  height = 10                  
)

# 再次绘制并保存为 PDF (矢量图)
pheatmap(
  heatmap_matrix_log2,
  color = heatmap_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  main = paste0("Top ", n_top, " Inflammation & Macrophage DEGs Heatmap"),
  fontsize = 10,
  fontsize_row = row_fontsize,
  filename = output_pdf, 
  width = 8,
  height = 10
)

cat(paste("成功保存 Top 30 基因热图到:\n", output_png, "\n", output_pdf, "\n"))


# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# 检查和安装核心包
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(readr))       
suppressPackageStartupMessages(library(RColorBrewer)) 

# --- 2. 定义文件路径 ---
# (输入) 上一步筛选出的特定基因列表
input_csv <- "G:/XY/29/analysis/5_Specific_Gene_Analysis/Specific_DEGs_Inflammation_Macrophage.csv"

# (输出) 结果保存目录
output_dir <- "G:/XY/29/analysis/5_Specific_Gene_Analysis"
# (输出) 图片文件名
output_png <- file.path(output_dir, "Top30_Specific_DEGs_Heatmap.png")
output_pdf <- file.path(output_dir, "Top30_Specific_DEGs_Heatmap.pdf")

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取特定基因文件:", input_csv, "\n"))
# 检查文件是否存在
if (!file.exists(input_csv)) {
  stop(paste("错误: 输入文件未找到:", input_csv, "\n请确保上一步已成功运行。"))
}

# 使用 read_csv 读取
data_specific <- read_csv(input_csv, col_types = cols(.default = "c")) 

# *** 关键修正 ***
# 将 read_csv 读入的 'tibble' 转换回 'data.frame'，以解决 'list' 错误
data_specific <- as.data.frame(data_specific)
# *** 修正完毕 ***

cat(paste("成功读取", nrow(data_specific), "个基因。\n"))

# 转换数值列
sample_cols <- c("con_1", "con_2", "con_3", "siha_1", "siha_2", "siha_3")
for (col in sample_cols) {
  data_specific[, col] <- as.numeric(data_specific[, col])
}
cat("样本TPM值已成功转换为数值。\n")

# --- 4. 筛选 Top 30 并准备矩阵 ---
n_top <- 30
data_top30 <- head(data_specific, n_top)
cat(paste("已提取 Top", n_top, "个基因用于绘图。\n"))

# 提取表达矩阵 (TPM 值)
heatmap_matrix_raw <- data_top30[, sample_cols]

# --- 4.1 设置行名 (使用 SYMBOL) ---
if ("SYMBOL" %in% colnames(data_top30) && !all(is.na(data_top30$SYMBOL))) {
  gene_names <- data_top30$SYMBOL
  cat("使用 'SYMBOL' 作为热图的行名。\n")
} else {
  gene_names <- data_top30$gene 
  cat("警告: 'SYMBOL' 列无效, 将使用 'gene' (GeneID) 作为行名。\n")
}

# 检查是否有重复的基因名，并处理
if (any(duplicated(gene_names))) {
  cat("警告：存在重复的基因名，将使用 make.unique() 处理。\n")
  rownames(heatmap_matrix_raw) <- make.unique(gene_names)
} else {
  rownames(heatmap_matrix_raw) <- gene_names
}

# --- 4.2 关键转换: Log2(TPM + 1) ---
heatmap_matrix_log2 <- log2(heatmap_matrix_raw + 1)

# --- 5. "Nature" 风格设置 (注释和颜色) ---
# 5.1 样本 (列) 注释
annotation_col <- data.frame(
  Group = factor(c(rep("con", 3), rep("siha", 3)))
)
rownames(annotation_col) <- sample_cols

# 5.2 分组颜色 (沿用您火山图的蓝/红)
ann_colors <- list(
  Group = c(con = "#333aab", siha = "#a73336") 
)

# 5.3 热图颜色 (蓝-白-红 梯度)
heatmap_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# --- 6. 绘制热图 ---
cat("开始绘制热图并保存...\n")

row_fontsize <- 8 

p_heatmap <- pheatmap(
  heatmap_matrix_log2,         
  color = heatmap_colors,      
  scale = "row",               
  cluster_rows = TRUE,         
  cluster_cols = TRUE,         
  annotation_col = annotation_col, 
  annotation_colors = ann_colors,  
  show_rownames = TRUE,        
  show_colnames = TRUE,        
  border_color = NA,           
  main = paste0("Top ", n_top, " Inflammation & Macrophage DEGs Heatmap"),
  fontsize = 10,
  fontsize_row = row_fontsize, 
  filename = output_png,       
  width = 8,                   
  height = 10                  
)

# 再次绘制并保存为 PDF (矢量图)
pheatmap(
  heatmap_matrix_log2,
  color = heatmap_colors,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  main = paste0("Top ", n_top, " Inflammation & Macrophage DEGs Heatmap"),
  fontsize = 10,
  fontsize_row = row_fontsize,
  filename = output_pdf, 
  width = 8,
  height = 10
)

cat(paste("成功保存 Top 30 基因热图到:\n", output_png, "\n", output_pdf, "\n"))


# --- 1. 环境准备 ---
# 清理环境
rm(list = ls()) 

# --- 1.1 安装 Bioconductor 核心包 ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
# 'org.Hs.eg.db' 是【人类】的注释数据库 (Hs = Homo sapiens)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}


# --- 1.2 加载包 ---
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db)) # 人类数据库
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble)) # 用于处理 data.frame

# --- 2. 定义文件路径 ---
# (输入1) 差异分析总表 (用于提取“背景基因”)
input_all <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"
# (输入2) 上一步筛选的特定基因 (用于提取“前景基因”)
input_specific <- "G:/XY/29/analysis/5_Specific_Gene_Analysis/Specific_DEGs_Inflammation_Macrophage.csv"

# (输出) 新的富集分析目录
output_dir <- "G:/XY/29/analysis/6_Enrichment_Analysis"
# (输出) 保存转换后的 RData 文件
output_rdata <- file.path(output_dir, "Enrichment_Input_Data.RData")

# 自动创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(paste("已创建新目录:", output_dir, "\n"))
}

# --- 3. 读取基因列表 ---
cat("正在读取前景基因 (Specific DEGs)...\n")
data_specific <- read_csv(input_specific, col_types = cols(.default = "c"))
data_specific <- as.data.frame(data_specific)

cat("正在读取背景基因 (Universe)...\n")
data_all <- read_tsv(input_all, col_types = cols(.default = "c"))
data_all <- as.data.frame(data_all)

# --- 3.1 提取 SYMBOL 列表 ---
genes_specific <- unique(data_specific$SYMBOL)
genes_universe <- unique(data_all$SYMBOL)
genes_specific <- genes_specific[!is.na(genes_specific) & genes_specific != "-"]
genes_universe <- genes_universe[!is.na(genes_universe) & genes_universe != "-"]

cat(paste("前景基因 (Specific): 提取到", length(genes_specific), "个唯一的 SYMBOL。\n"))
cat(paste("背景基因 (Universe): 提取到", length(genes_universe), "个唯一的 SYMBOL。\n"))

# --- 4. 将 SYMBOL 转换为 ENTREZID ---
cat("正在将 SYMBOL 转换为 ENTREZID (这可能需要一点时间)...\n")

# 转换前景基因
ids_specific_df <- bitr(genes_specific, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Hs.eg.db")

# 转换背景基因
ids_universe_df <- bitr(genes_universe, 
                        fromType = "SYMBOL", 
                        toType = "ENTREZID", 
                        OrgDb = "org.Hs.eg.db")

cat("--- ID 转换总结 ---\n")
cat(paste("前景基因:", length(genes_specific), "个 SYMBOL 成功转换为", nrow(ids_specific_df), "个 ENTREZID。\n"))
cat(paste("背景基因:", length(genes_universe), "个 SYMBOL 成功转换为", nrow(ids_universe_df), "个 ENTREZID。\n"))

# --- 5. 准备最终的向量并保存 ---
genes_specific_entrez <- unique(ids_specific_df$ENTREZID)
genes_universe_entrez <- unique(ids_universe_df$ENTREZID)

save(genes_specific_entrez, genes_universe_entrez, file = output_rdata)

cat("-------------------\n")
cat(paste("成功将转换后的 ID 列表保存到:\n", output_rdata, "\n"))
cat("准备工作完成！\n")


# --- 1. 环境准备 ---
# (与之前相同)
rm(list = ls()) 
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))

# --- 2. 定义文件路径 ---
# (与之前相同)
input_rdata <- "G:/XY/29/analysis/6_Enrichment_Analysis/Enrichment_Input_Data.RData"
output_dir <- "G:/XY/29/analysis/6_Enrichment_Analysis"

# --- 3. 加载准备好的数据 ---
# (与之前相同)
cat(paste("正在加载已转换的基因ID文件:", input_rdata, "\n"))
if (!file.exists(input_rdata)) {
  stop("错误: RData 文件未找到。请先运行脚本3。")
}
load(input_rdata)
cat("数据加载完毕。\n")


# --- 4. GO (Gene Ontology) 富集分析 ---
# (分别运行 BP, CC, MF)
cat("\n--- 正在运行 GO-BP (Biological Process) ... ---\n")
go_bp_results <- enrichGO(gene          = genes_specific_entrez,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID",
                          ont           = "BP", 
                          pAdjustMethod = "BH",
                          universe      = genes_universe_entrez,
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

cat("\n--- 正在运行 GO-CC (Cellular Component) ... ---\n")
go_cc_results <- enrichGO(gene          = genes_specific_entrez,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID",
                          ont           = "CC", 
                          pAdjustMethod = "BH",
                          universe      = genes_universe_entrez,
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

cat("\n--- 正在运行 GO-MF (Molecular Function) ... ---\n")
go_mf_results <- enrichGO(gene          = genes_specific_entrez,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = "ENTREZID",
                          ont           = "MF", 
                          pAdjustMethod = "BH",
                          universe      = genes_universe_entrez,
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

cat("GO 分析完毕。\n")

# --- 4.1 简化 GO 结果 ---
go_bp_simple <- NULL
if (!is.null(go_bp_results) && nrow(go_bp_results) > 0) {
  cat("正在简化 GO-BP 结果...\n")
  go_bp_simple <- simplify(go_bp_results, cutoff = 0.7, by = "p.adjust", select_fun = min)
} else { cat("警告：GO-BP 未返回任何显著结果。\n") }

go_cc_simple <- NULL
if (!is.null(go_cc_results) && nrow(go_cc_results) > 0) {
  cat("正在简化 GO-CC 结果...\n")
  go_cc_simple <- simplify(go_cc_results, cutoff = 0.7, by = "p.adjust", select_fun = min)
} else { cat("警告：GO-CC 未返回任何显著结果。\n") }

go_mf_simple <- NULL
if (!is.null(go_mf_results) && nrow(go_mf_results) > 0) {
  cat("正在简化 GO-MF 结果...\n")
  go_mf_simple <- simplify(go_mf_results, cutoff = 0.7, by = "p.adjust", select_fun = min)
} else { cat("警告：GO-MF 未返回任何显著结果。\n") }


# --- 5. KEGG 通路分析 ---
cat("\n--- 正在运行 KEGG 通路分析... ---\n")
kegg_results <- NULL
tryCatch({
  kegg_results <- enrichKEGG(gene          = genes_specific_entrez,
                             organism      = "hsa", 
                             universe      = genes_universe_entrez,
                             pAdjustMethod = "BH",
                             pvalueCutoff  = 0.01,
                             qvalueCutoff  = 0.05)
  cat("KEGG 分析完毕。\n")
  if (is.null(kegg_results) || nrow(kegg_results) == 0) {
     cat("警告：KEGG 分析未返回任何显著结果。\n")
  }
}, error = function(e) {
  cat(paste("KEGG 分析失败 (可能是网络问题):", e$message, "\n"))
})


# --- 6. 绘图与保存 (Nature 风格气泡图) ---

# 定义 "Nature 风格" 绘图函数 (美化 + 高亮)
plot_enrichment_dotplot <- function(enrich_result, title, show_n = 20) {
  
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    cat(paste("跳过绘图:", title, "(因为没有显著结果)\n"))
    return(NULL)
  }
  
  cat(paste("正在绘制图表:", title, "\n"))
  
  # 1. 生成基础图表对象
  p_base <- dotplot(enrich_result, 
                    showCategory = show_n, 
                    font.size = 10,
                    title = title)
  
  # 2. 定义您关心的关键词 (用于标红)
  keywords_to_highlight <- c(
    "cytokine", "chemokine", "inflammatory", "macrophage", 
    "leukocyte", "NF-kappa B", "JAK-STAT", "Phagosome", 
    "inflammasome", "NADPH oxidase", "chemotaxis", "immune", 
    "receptor binding"
  )
  
  # 3. 创建关键词搜索模式
  pattern <- paste(keywords_to_highlight, collapse = "|")
  
  # 4. 获取 Y 轴标签 (按 ggplot 绘图的顺序)
  plot_labels <- levels(p_base$data$Description)
  
  # 5. 根据关键词创建颜色向量 ("red" 或 "black")
  label_colors <- ifelse(grepl(pattern, plot_labels, ignore.case = TRUE), 
                           "red",  # 高亮颜色
                           "black") # 默认颜色
  
  # 6. 应用 "Nature" 风格 和 高亮颜色
  p_styled <- p_base + 
    theme_bw() + # 使用白底黑框
    theme(
      panel.grid.major = element_blank(), # 移除主网格线
      panel.grid.minor = element_blank(), # 移除次网格线
      panel.border = element_rect(color = "black", fill = NA, size = 1), # 加强边框
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # 标题居中加粗
      
      # *** 应用高亮颜色 ***
      axis.text.y = element_text(color = label_colors, face = "bold"), # Y轴标签颜色
      axis.text.x = element_text(face = "bold"),
      
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    ) +
    
    # 7. (可选) 替换为更"高级"的颜色梯度 (匹配您热图的 蓝-红)
    scale_color_gradient(low = "#a73336", high = "#333aab", name = "FDR") +
    
    # 8. (可选) 重命名图例
    labs(size = "Gene Count")
  
  return(p_styled)
}


# --- 6.1 绘制 GO 气泡图 (使用简化后的结果) ---
p_go_bp <- plot_enrichment_dotplot(go_bp_simple, "GO Biological Process (Top 20)", 20)
p_go_cc <- plot_enrichment_dotplot(go_cc_simple, "GO Cellular Component (Top 20)", 20)
p_go_mf <- plot_enrichment_dotplot(go_mf_simple, "GO Molecular Function (Top 20)", 20)

# (保存)
if(!is.null(p_go_bp)) ggsave(file.path(output_dir, "DotPlot_GO_BP_v2_Nature.png"), p_go_bp, width = 10, height = 10, dpi = 300)
if(!is.null(p_go_bp)) ggsave(file.path(output_dir, "DotPlot_GO_BP_v2_Nature.pdf"), p_go_bp, width = 10, height = 10)
if(!is.null(p_go_cc)) ggsave(file.path(output_dir, "DotPlot_GO_CC_v2_Nature.png"), p_go_cc, width = 10, height = 8, dpi = 300)
if(!is.null(p_go_cc)) ggsave(file.path(output_dir, "DotPlot_GO_CC_v2_Nature.pdf"), p_go_cc, width = 10, height = 8)
if(!is.null(p_go_mf)) ggsave(file.path(output_dir, "DotPlot_GO_MF_v2_Nature.png"), p_go_mf, width = 10, height = 8, dpi = 300)
if(!is.null(p_go_mf)) ggsave(file.path(output_dir, "DotPlot_GO_MF_v2_Nature.pdf"), p_go_mf, width = 10, height = 8)

# (显示)
if(!is.null(p_go_bp)) print(p_go_bp)
if(!is.null(p_go_cc)) print(p_go_cc)
if(!is.null(p_go_mf)) print(p_go_mf)

# (保存数据表)
if (!is.null(go_bp_simple)) write.csv(as.data.frame(go_bp_simple), file.path(output_dir, "Enrichment_Results_GO_BP_Simple.csv"), row.names = FALSE)
if (!is.null(go_cc_simple)) write.csv(as.data.frame(go_cc_simple), file.path(output_dir, "Enrichment_Results_GO_CC_Simple.csv"), row.names = FALSE)
if (!is.null(go_mf_simple)) write.csv(as.data.frame(go_mf_simple), file.path(output_dir, "Enrichment_Results_GO_MF_Simple.csv"), row.names = FALSE)


# --- 6.2 绘制 KEGG 气泡图 ---
if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
  
  kegg_results_readable <- setReadable(kegg_results, OrgDb = 'org.Hs.eg.db', keyType = "ENTREZID")
  
  # *** 调用新的绘图函数 ***
  p_kegg <- plot_enrichment_dotplot(kegg_results_readable, "KEGG Pathways (Top 20)", 20)
  
  # (保存)
  if(!is.null(p_kegg)) ggsave(file.path(output_dir, "DotPlot_KEGG_v2_Nature.png"), p_kegg, width = 10, height = 10, dpi = 300)
  if(!is.null(p_kegg)) ggsave(file.path(output_dir, "DotPlot_KEGG_v2_Nature.pdf"), p_kegg, width = 10, height = 10)
  
  if(!is.null(p_kegg)) print(p_kegg)
  
  # (保存数据表)
  write.csv(as.data.frame(kegg_results_readable), file.path(output_dir, "Enrichment_Results_KEGG.csv"), row.names = FALSE)
}

cat("\n--- 富集分析流程全部完成 (已应用美化和高亮) --- \n")
cat(paste("请检查输出目录:\n", output_dir, "\n"))
cat("新的图片文件已保存为 v2_Nature.png / .pdf \n")



# -------------------------------------------------------------------------
# 最终火山图绘制脚本 (V6 V4 版 - 标注 P-value Top 10 Up/Down)
# -------------------------------------------------------------------------

# --- 1. 环境准备 ---
rm(list = ls()) 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel)) 
suppressPackageStartupMessages(library(ggnewscale)) 

# --- 2. 定义文件路径 ---
input_file <- "G:/XY/29/结题报告-20251105-RNA-CQT2025092609-F001/src/4.DiffEnrich/diff.merged.annotated.xls"
output_dir <- "G:/XY/29/analysis"
# (新文件名)
output_png <- file.path(output_dir, "KD_vs_NC_volcano_plot_v4_Top_Pvalue.png")
output_pdf <- file.path(output_dir, "KD_vs_NC_volcano_plot_v4_Top_Pvalue.pdf")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# --- 3. 读取与预处理数据 ---
cat(paste("正在读取文件:", input_file, "\n"))
data <- read_tsv(input_file, col_types = cols(.default = "c")) 
data <- as.data.frame(data)
cols_to_numeric <- c("logFC.KD-vs-NC", "PValue.KD-vs-NC", "FDR.KD-vs-NC")
for (col in cols_to_numeric) {
  data[, col] <- as.numeric(data[, col])
}
data$log2FoldChange <- data$`logFC.KD-vs-NC`
data$pvalue_raw     <- data$`PValue.KD-vs-NC`
data$pvalue_adj     <- data$`FDR.KD-vs-NC`
data <- data[!is.na(data$log2FoldChange) & !is.na(data$pvalue_raw) & !is.na(data$pvalue_adj), ]
data$neg_log10_pvalue <- -log10(data$pvalue_raw)
infinite_rows <- is.infinite(data$neg_log10_pvalue)
if (any(infinite_rows)) {
  max_finite_p <- max(data$neg_log10_pvalue[!infinite_rows], na.rm = TRUE)
  data$neg_log10_pvalue[infinite_rows] <- max_finite_p * 1.1 
}
plot_data <- data

# --- 4. 定义阈值并拆分数据 ---
logfc_threshold <- 2  
fdr_threshold   <- 0.01 
data_ns <- subset(plot_data, abs(log2FoldChange) <= logfc_threshold | pvalue_adj >= fdr_threshold)
data_sig <- subset(plot_data, abs(log2FoldChange) > logfc_threshold & pvalue_adj < fdr_threshold)
data_up <- subset(data_sig, log2FoldChange > 0)
data_down <- subset(data_sig, log2FoldChange < 0)
n_up <- nrow(data_up)
n_down <- nrow(data_down)
cat(paste("筛选完毕：总计", (n_up + n_down), "个差异基因 (上调:", n_up, ", 下调:", n_down, ")\n"))

# --- 4.1 过滤极高P值的点 ---
y_limit_display <- 320
plot_data <- subset(plot_data, neg_log10_pvalue <= y_limit_display)
data_ns     <- subset(data_ns,   neg_log10_pvalue <= y_limit_display)
data_up     <- subset(data_up,   neg_log10_pvalue <= y_limit_display)
data_down   <- subset(data_down, neg_log10_pvalue <= y_limit_display)
cat(paste("已过滤掉 neg_log10(PValue) >", y_limit_display, "的点，不再用于绘图。\n"))


# --- 5. (新) 准备要标注的标签 (按 P-value 排序) ---
cat("正在准备 Top 10 (P-value) 的上调和下调基因标签...\n")

# A. (上边) 上调 Top 10 (按 P-value 排序)
label_up_pvalue <- data_up[order(data_up$neg_log10_pvalue, decreasing = TRUE), ]
data_label_up <- head(label_up_pvalue, 10)

# B. (上边) 下调 Top 10 (按 P-value 排序)
label_down_pvalue <- data_down[order(data_down$neg_log10_pvalue, decreasing = TRUE), ]
data_label_down <- head(label_down_pvalue, 10)

# C. 合并所有标签 (共 20 个)
data_labels <- rbind(data_label_up, data_label_down)

cat(paste("将标注 Top 10 (P-value) Up 和 Top 10 (P-value) Down 共", nrow(data_labels), "个基因。\n"))
# --- 标签准备完毕 ---


# --- 6. 开始绘图 (V6 模板) ---
cat("开始绘制火山图...\n")
p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  
  geom_point(data = data_ns, color = "grey", size = 1.2, alpha = 0.6, stroke = 0) +
  geom_point(data = data_down, aes(color = neg_log10_pvalue, size = neg_log10_pvalue), alpha = 0.8, stroke = 0) +
  scale_color_gradient(low = "#b5dbe6", high = "#333aab", name = "-log10(p-value) Down") +
  new_scale_color() +
  geom_point(data = data_up, aes(color = neg_log10_pvalue, size = neg_log10_pvalue), alpha = 0.8, stroke = 0) +
  scale_color_gradient(low = "#fe8264", high = "#a73336", name = "-log10(p-value) Up") +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "#999999") +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "#999999") + 
  scale_size_continuous(range = c(1, 4)) +
  scale_y_continuous(expand = c(0, 0)) +  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  ) +
  guides(size = "none") + 
  
  # 6.1 统计标注
  annotate("text",
           x = 18, y = 300, 
           label = paste("Up:", n_up),
           size = 5, color = "#a73336", fontface = "bold") +
  annotate("text",
           x = -18, y = 300, 
           label = paste("Down:", n_down),
           size = 5, color = "#333aab", fontface = "bold") +
  
  # 6.2 基因标签 (使用 P-value 排序的 data_labels)
  geom_text_repel(
    data = data_labels, 
    aes(label = SYMBOL),  
    color = "black",      
    size = 3.5,           
    fontface = "bold",    
    max.overlaps = Inf,   
    box.padding = 0.6,    
    point.padding = 0.2,  
    segment.color = 'grey50', 
    segment.size = 0.3      
  ) +
  
  xlab("Log2 Fold Change (KD vs NC)") +
  ylab("-log10(P-VALUE)") +
  ggtitle("Volcano Plot: KD vs NC (Top 10 Significant Labeled)") + # 标题
  
  coord_cartesian(xlim = c(-20, 20), ylim = c(0, 320))

# --- 7. 显示并保存图片 ---
cat("绘图完毕，正在保存...\n")
ggsave(output_png, plot = p, width = 10, height = 8, units = "in", dpi = 300)
ggsave(output_pdf, plot = p, width = 10, height = 8, units = "in")
cat(paste("成功保存图片到:\n", output_png, "\n", output_pdf, "\n"))

print(p)