# 1. 加载 Seurat 库
library(Seurat)

cat("[1/4] 正在加载 Seurat 对象...\n")

# 2. 指定文件路径
file_path <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/WT_E10.5_Rep1_seurat3_scTransform.rds"

# 3. 使用 readRDS() 函数读取文件
# (这一步仍然可能显示旧对象)
scrna_data_old <- readRDS(file_path)

cat("[2/4] 文件加载完毕。正在更新 Seurat 对象版本...\n")

# 4. 【关键修复步骤】更新 Seurat 对象
# 这是解决您遇到的 Error 所必需的
scrna_data <- UpdateSeuratObject(scrna_data_old)

cat("[3/4] 对象更新完毕。\n\n")

# 5. 查看 Seurat 对象的整体概览
cat("--- Seurat 对象整体概览 ---\n")
print(scrna_data)
cat("\n")

# 6. 查看 Assay (分析) 信息
cat("--- 可用的 Assays (例如 'RNA', 'SCT') ---\n")
print(Assays(scrna_data))
cat("\n")

# 7. 查看 metadata (元数据) 的格式
cat("--- Metadata 结构 (前 6 行) ---\n")
print(head(scrna_data@meta.data))
cat("\n")

# 8. 查看已有的降维信息 (例如 'pca', 'umap')
cat("--- 已计算的降维聚类 ---\n")
print(Reductions(scrna_data))
cat("\n")

cat("[4/4] 数据格式检查完毕。\n")


# 1. 加载 Seurat 库
library(Seurat)

cat("[1/4] 正在加载所有 4 个 scRNA-seq Seurat 对象...\n\n")

# 2. 定义文件路径
# (我根据您的目录结构和样本命名来推测, 尤其是 E11.5)
file_paths <- c(
  E105 = "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/WT_E10.5_Rep1_seurat3_scTransform.rds",
  E115 = "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/WT_E11.5_Rep2_seurat3_scTransform.rds",
  E125 = "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/WT_E12.5_Rep1_seurat3_scTransform.rds",
  E145 = "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/WT_E14.5_Rep1_seurat3_scTransform.rds"
)

# 3. 循环读取、更新对象，并存入列表
# 我们使用 lapply 来处理这个列表
scrna_list <- lapply(names(file_paths), function(sample_name) {
  
  cat("  正在处理:", sample_name, "\n")
  cat("    路径:", file_paths[sample_name], "\n")
  
  # 3.1 读取 RDS 文件
  data_old <- readRDS(file_paths[sample_name])
  
  # 3.2 更新 Seurat 对象 (修复版本不兼容问题)
  cat("    正在运行 UpdateSeuratObject...\n")
  data_new <- UpdateSeuratObject(data_old)
  
  # 3.3 (可选但推荐) 统一元数据
  # 我们添加一个 'samplename' 列, 名字就是我们命名的 E105, E115...
  data_new$samplename <- sample_name
  
  cat("    处理完毕:", ncol(data_new), "个细胞\n\n")
  return(data_new)
})

# 4. 为列表命名 (方便后续查看)
names(scrna_list) <- names(file_paths)

cat("[2/4] 所有 4 个 scRNA-seq 对象加载并更新完毕。\n")
cat("      已创建 'scrna_list'，包含以下对象:\n")
print(names(scrna_list))
cat("\n")

# -- 整合流程开始 --

cat("[3/4] 准备开始 SCT 整合流程的第一步: 选择整合特征...\n")

# 5. 选择用于整合的特征 (SCT 流程)
# 默认会选择 3000 个高变基因
# 这一步是为下一步寻找锚点做准备
features <- SelectIntegrationFeatures(object.list = scrna_list, nfeatures = 3000)

cat("[4/4] 整合特征 (Integration features) 选择完毕。\n")

cat("[1/4] 开始执行 PrepSCTIntegration...\n")
cat("      (此步骤将为 'scrna_list' 中的每个对象准备 SCT 整合)\n")

# 1. 准备 SCT 整合
# (注意：这一步会修改 scrna_list 中的对象)
scrna_list <- PrepSCTIntegration(
  object.list = scrna_list, 
  anchor.features = features, 
  verbose = TRUE # 显示详细进度
)

cat("\n[2/4] PrepSCTIntegration 步骤完成。\n")
cat("--------------------------------------------------\n")
cat("[3/4] 开始寻找整合锚点 (FindIntegrationAnchors)...\n")
cat("      (这是计算量最大的一步，可能需要较长时间)\n")

# 2. 寻找整合锚点 (Anchors)
# 我们使用 SCT 作为归一化方法 (normalization.method = 'SCT')
# 我们使用 "cca" (Canonical Correlation Analysis) 来寻找锚点，这是标准流程
integration_anchors <- FindIntegrationAnchors(
  object.list = scrna_list, 
  normalization.method = "SCT", 
  anchor.features = features,
  reduction = "cca", # 推荐使用 cca 
  verbose = TRUE # 显示详细进度
)

cat("\n[4/4] 寻找锚点 (FindIntegrationAnchors) 步骤完成！\n")
cat("      已成功创建 'integration_anchors' 对象。\n")

cat("[1/3] 开始整合数据 (IntegrateData)...\n")
cat("      (此步骤将使用找到的锚点来创建一个新的、已整合的 Seurat 对象)\n")

# 1. 整合数据
# 我们将 'integration_anchors' 传入 IntegrateData 函数
# 这一步会创建一个新的 Seurat 对象
scrna_integrated <- IntegrateData(
  anchorset = integration_anchors,
  normalization.method = "SCT", # 必须和上一步保持一致
  verbose = TRUE # 显示详细进度
)

cat("\n[2/3] 数据整合完毕！已创建 'scrna_integrated' 对象。\n")
cat("--------------------------------------------------\n")

# 2. (重要) 切换默认 Assay
# 默认情况下，Seurat 可能会激活 'RNA' assay，
# 但我们所有的下游分析(PCA, UMAP)都必须在新的 'integrated' assay 上进行
DefaultAssay(scrna_integrated) <- "integrated"

# 3. 查看整合后的对象概览
cat("[3/3] 查看整合后的 Seurat 对象概览：\n\n")
print(scrna_integrated)

cat("\n--- 整合流程全部完成！---\n")

cat("[1/3] 正在创建用于保存结果的新文件夹...\n")

# 1. 定义输出路径和新文件夹名称
output_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/Integration_Output"

# 2. 创建新文件夹
# recursive = TRUE 确保在需要时创建父目录
# showWarnings = FALSE 避免在文件夹已存在时显示警告
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("      新文件夹已创建(或已存在):", output_dir, "\n")

# 3. 定义最终保存的文件路径
save_file_path <- file.path(output_dir, "scrna_integrated.rds")

cat("[2/3] 正在保存整合后的 Seurat 对象...\n")
cat("      (这可能需要几分钟，因为对象很大: 28,215 个细胞)\n")

# 4. 保存对象
saveRDS(scrna_integrated, file = save_file_path)

cat("[3/3] 恭喜！整合后的对象已成功保存到:\n")
cat("     ", save_file_path, "\n")

cat("[0/8] 加载 R 库和上次保存的对象...\n")
library(Seurat)
library(stringr) # 用于基因名转换

# 定义文件路径 (与之前相同)
save_file_path <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/Integration_Output/scrna_integrated.rds"

# 加载我们刚刚保存的对象
scrna_integrated <- readRDS(save_file_path)

cat("      对象加载完毕。\n")

cat("[1/8] 【关键修复】正在转换细胞周期基因名为小鼠格式...\n")

# 1. 转换人类基因列表 (e.g. "MCM5") 为小鼠基因列表 (e.g. "Mcm5")
# 我们通过 stringr::str_to_title() 函数实现
s.genes.mouse <- str_to_title(tolower(cc.genes.updated.2019$s.genes))
g2m.genes.mouse <- str_to_title(tolower(cc.genes.updated.2019$g2m.genes))

cat("      基因名转换完毕。\n")

cat("[2/8] 正在 (正确地) 重新运行细胞周期评分...\n")

# 2. 重新运行 CellCycleScoring, 使用修正后的小鼠基因列表
scrna_integrated <- CellCycleScoring(
  scrna_integrated, 
  s.features = s.genes.mouse, 
  g2m.features = g2m.genes.mouse, 
  assay = 'RNA'
)

cat("      细胞周期评分 (修正版) 完成。\n")
cat("      (您这次不应该再看到 'features are not present' 的警告了)\n\n")

cat("[3/8] 正在 (正确地) 重新运行 PCA 并回归细胞周期...\n")

# 3. 重新运行 PCA (修正版)
scrna_integrated <- RunPCA(
  scrna_integrated, 
  assay = "integrated",
  features = VariableFeatures(object = scrna_integrated), 
  vars.to.regress = c("S.Score", "G2M.Score"), # 现在这些分数是准确的
  verbose = FALSE
)

cat("[4/8] 正在 (正确地) 重新寻找细胞近邻 (FindNeighbors)...\n")

# 4. 重新寻找近邻 (FindNeighbors)
scrna_integrated <- FindNeighbors(
  scrna_integrated, 
  reduction = "pca", 
  dims = 1:30
)

cat("[5/8] 正在 (正确地) 重新进行聚类 (FindClusters)...\n")

# 5. 重新聚类 (FindClusters)
scrna_integrated <- FindClusters(
  scrna_integrated, 
  resolution = 0.5
)

cat("[6/8] 正在 (正确地) 重新计算 UMAP 坐标...\n")

# 6. 重新运行 UMAP
scrna_integrated <- RunUMAP(
  scrna_integrated, 
  reduction = "pca", 
  dims = 1:30
)

cat("[7/8] 正在计算 tSNE 坐标 (根据您的要求)...\n")
# 7. 运行 tSNE
scrna_integrated <- RunTSNE(
  scrna_integrated,
  reduction = "pca",
  dims = 1:30
)

cat("--- 所有分析均已正确完成！---\n")

# 8. (重要) 保存这个最终的、正确的对象
cat("\n[8/8] 正在保存最终的、已校正的对象 (覆盖原文件)...\n")
saveRDS(scrna_integrated, file = save_file_path)
cat("对象已更新并保存。\n")


