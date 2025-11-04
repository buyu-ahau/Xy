# 1. 加载必要的库
# (如果未安装，请先运行 install.packages("rlang") )
library(ggplot2)
library(rlang) # 用于检查R环境中的对象

# 2. 定义路径
base_dir <- "G:/XY/22号/analysis/08_DESeq2_results"
new_folder_name <- "Autophagy_GSEA_Analysis" # 新文件夹的名称
new_folder_path <- file.path(base_dir, new_folder_name)

# 3. 创建新文件夹
#    showWarnings = FALSE 确保如果文件夹已存在也不会报错
dir.create(new_folder_path, showWarnings = FALSE, recursive = TRUE)
cat("已创建或确认文件夹:", new_folder_path, "\n\n")

# 4. 移动已经生成的 Excel 和 PNG 文件
#    (这些文件目前在 'base_dir' 中)
files_to_move <- c(
  "GSEA_autophagy_summary.xlsx", # 这是上上步保存的 Excel
  "GSEA_NC_style_KEGG.png",      # 这是上一步保存的 PNG
  "GSEA_NC_style_GO_BP.png",
  "GSEA_NC_style_Reactome.png"
)

cat("--- 正在移动文件... ---\n")
for (file in files_to_move) {
  old_path <- file.path(base_dir, file)
  new_path <- file.path(new_folder_path, file)
  
  if (file.exists(old_path)) {
    file.rename(from = old_path, to = new_path)
    cat("  已移动:", file, "\n")
  } else {
    cat("  警告: 未在 '08_DESeq2_results' 中找到文件", file, "，跳过移动。\n")
  }
}

# 5. 保存 PDF 版本的图表
#    (这假设 kegg_plot, go_plot, reactome_plot 仍在 R 环境中)
cat("\n--- 正在生成 PDF... ---\n")

# 检查绘图对象是否存在
plot_objects_exist <- all(
  rlang::env_has(globalenv(), "kegg_plot"),
  rlang::env_has(globalenv(), "go_plot"),
  rlang::env_has(globalenv(), "reactome_plot")
)

if (plot_objects_exist) {
  # 保存 KEGG PDF
  kegg_pdf_path <- file.path(new_folder_path, "GSEA_NC_style_KEGG.pdf")
  ggsave(kegg_pdf_path, kegg_plot, width = 8, height = 4, device = cairo_pdf)
  cat("  已保存: GSEA_NC_style_KEGG.pdf\n")
  
  # 保存 GO PDF
  go_pdf_path <- file.path(new_folder_path, "GSEA_NC_style_GO_BP.pdf")
  ggsave(go_pdf_path, go_plot, width = 10, height = 10, device = cairo_pdf)
  cat("  已保存: GSEA_NC_style_GO_BP.pdf\n")
  
  # 保存 Reactome PDF
  reactome_pdf_path <- file.path(new_folder_path, "GSEA_NC_style_Reactome.pdf")
  ggsave(reactome_pdf_path, reactome_plot, width = 10, height = 8, device = cairo_pdf)
  cat("  已保存: GSEA_NC_style_Reactome.pdf\n")
  
} else {
  cat("  --- 警告 ---\n")
  cat("  未在R环境中找到 'kegg_plot', 'go_plot' 或 'reactome_plot' 对象。\n")
  cat("  无法自动生成 PDF 文件。\n")
  cat("  请先重新运行上一步的绘图代码，然后再运行此代码的第 5 部分。\n")
}

# 6. 保存我们最后使用的 R 脚本 (作为"讨论"的记录)
#    (这是上一步的绘图代码)
our_script_code <- '
# --- GSEA 绘图脚本 (Autophagy Analysis) ---
#
# 1. 加载R包
library(ggplot2)
library(dplyr)
library(writexl) # 确保已加载

# 2. 检查数据
# (此脚本假设 `found_kegg`, `found_go`, `found_reactome` 
#  已通过 GSEA 搜索代码加载到 R 环境中)
if (!exists("found_kegg") || !exists("found_go") || !exists("found_reactome")) {
  stop("错误: 未找到 'found_kegg', 'found_go' 或 'found_reactome' 变量。")
}

# 3. 定义 "Nature 风格" 的主题
nature_style_theme <- function() {
  theme_bw() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black")
  )
}

# 4. 定义颜色方案 (浅灰 到 深红)
color_scale <- scale_fill_gradient(
  low = "grey80",
  high = "darkred",
  name = "-log10(FDR)"
)

# 5. 定义一个辅助函数来计算 -log10(FDR)
calculate_log10_fdr <- function(df) {
  non_zero_min_fdr <- min(df$`FDR q-val`[df$`FDR q-val` > 0], na.rm = TRUE)
  if (!is.finite(non_zero_min_fdr)) { non_zero_min_fdr <- 1e-10 }
  
  df %>%
    mutate(
      log10_FDR = -log10(ifelse(`FDR q-val` == 0, non_zero_min_fdr / 10, `FDR q-val`))
    )
}

# --- 图 1: 绘制 KEGG ---
df_kegg_plot <- calculate_log10_fdr(found_kegg)

kegg_plot <- ggplot(df_kegg_plot, aes(x = NES, y = reorder(NAME, NES))) +
  geom_col(aes(fill = log10_FDR), width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  color_scale +
  nature_style_theme() +
  labs(
    title = "GSEA: KEGG Pathways",
    x = "Normalized Enrichment Score (NES)",
    caption = "NES > 0: Upregulated in Treat | NES < 0: Downregulated in Treat"
  )


# --- 图 2: 绘制 GO BP ---
df_go_plot <- calculate_log10_fdr(found_go)

go_plot <- ggplot(df_go_plot, aes(x = NES, y = reorder(NAME, NES))) +
  geom_col(aes(fill = log10_FDR), width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  color_scale +
  nature_style_theme() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(
    title = "GSEA: GO Biological Process Pathways",
    x = "Normalized Enrichment Score (NES)",
    caption = "NES > 0: Upregulated in Treat | NES < 0: Downregulated in Treat"
  )


# --- 图 3: 绘制 Reactome ---
df_reactome_plot <- calculate_log10_fdr(found_reactome)

reactome_plot <- ggplot(df_reactome_plot, aes(x = NES, y = reorder(NAME, NES))) +
  geom_col(aes(fill = log10_FDR), width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
  color_scale +
  nature_style_theme() +
  labs(
    title = "GSEA: Reactome Pathways",
    x = "Normalized Enrichment Score (NES)",
    caption = "NES > 0: Upregulated in Treat | NES < 0: Downregulated in Treat"
  )

# 打印图像 (如果你在 RStudio 中运行)
# print(kegg_plot)
# print(go_plot)
# print(reactome_plot)
'

# 将脚本写入文件
script_save_path <- file.path(new_folder_path, "autophagy_plotting_script.R")
writeLines(our_script_code, script_save_path)
cat("\n--- 正在保存R脚本... ---\n")
cat("  已保存: 我们的R脚本 (autophagy_plotting_script.R)\n")

cat("\n--- 全部完成! ---\n")
cat("所有结果 (Excel, 3x PNG, 3x PDF) 和 R 脚本都已保存到:\n")
cat(new_folder_path, "\n")