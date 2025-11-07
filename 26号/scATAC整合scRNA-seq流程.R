# 创建新的 signacmm10 环境并安装所有包
mamba create -n signacmm10 -c conda-forge -c bioconda \
  r-base \
  bioconductor-signac \
  r-seurat \
  bioconductor-genomicranges \
  bioconductor-rtracklayer \
  bioconductor-rsamtools \
  bioconductor-jaspar2020 \
  bioconductor-tfbstools \
  bioconductor-bsgenome.mmusculus.ucsc.mm10 \
  bioconductor-motifmatchr \
  r-matrix \
  r-dplyr \
  r-ggplot2 \
  r-patchwork \
  r-viridis

# 激活环境
conda activate signacmm10
# ============================================
library(Signac)
library(Seurat)
library(GenomicRanges)
library(Matrix)
library(rtracklayer)
# 查看当前 genes 对象有哪些列
cat("当前 genes 的列名:\n")
print(colnames(mcols(genes)))
# 将 gene_type 重命名为 gene_biotype
mcols(genes)$gene_biotype <- mcols(genes)$gene_type

# 验证一下
cat("现在的列名中是否包含 gene_biotype:\n")
print("gene_biotype" %in% colnames(mcols(genes)))
# 添加基因组注释到 Seurat 对象
cat("添加基因组注释到 Seurat 对象...\n")
Annotation(atac) <- genes
cat("✓ 注释添加成功!\n")
# 计算核小体信号
cat("计算核小体信号 (Nucleosome Signal)...\n")
atac <- NucleosomeSignal(object = atac)
cat("✓ 完成!\n\n")

# 查看结果
cat("核小体信号统计:\n")
print(summary(atac$nucleosome_signal))

# 计算 TSS 富集分数
cat("计算 TSS 富集分数 (TSS Enrichment)...\n")
cat("这需要几分钟，请耐心等待...\n")
atac <- TSSEnrichment(object = atac, fast = FALSE)
cat("✓ 完成!\n\n")

# 查看结果
cat("TSS 富集分数统计:\n")
print(summary(atac$TSS.enrichment))

# 计算 reads in peaks 比例和黑名单比例
cat("计算额外的 QC 指标...\n")
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$passed_filters
cat("✓ 完成!\n\n")

# 查看所有 QC 指标
cat("=== 所有 QC 指标统计 ===\n\n")

cat("【1】Peak region fragments:\n")
print(summary(atac$peak_region_fragments))

cat("\n【2】Pct reads in peaks (%):\n")
print(summary(atac$pct_reads_in_peaks))

cat("\n【3】TSS enrichment:\n")
print(summary(atac$TSS.enrichment))

cat("\n【4】Nucleosome signal:\n")
print(summary(atac$nucleosome_signal))

cat("\n【5】Blacklist ratio:\n")
print(summary(atac$blacklist_ratio))

# 设置过滤阈值
cat("=== 过滤分析 ===\n\n")
cat("推荐阈值:\n")
cat("  ✓ peak_region_fragments: 1,000 - 20,000\n")
cat("  ✓ pct_reads_in_peaks: > 15%\n")
cat("  ✓ TSS.enrichment: > 2\n")
cat("  ✓ nucleosome_signal: < 4\n")
cat("  ✓ blacklist_ratio: < 0.05\n\n")

# 计算过滤
n_before <- ncol(atac)
keep_cells <- atac$peak_region_fragments > 1000 &
  atac$peak_region_fragments < 20000 &
  atac$pct_reads_in_peaks > 15 &
  atac$TSS.enrichment > 2 &
  atac$nucleosome_signal < 4 &
  atac$blacklist_ratio < 0.05

n_after <- sum(keep_cells)

cat("过滤统计:\n")
cat("  - 过滤前细胞数:", n_before, "\n")
cat("  - 过滤后细胞数:", n_after, "\n")
cat("  - 保留比例:", round(n_after/n_before*100, 2), "%\n")
cat("  - 移除细胞数:", n_before - n_after, "\n\n")

# 查看各指标不合格的细胞数
cat("不合格细胞详情:\n")
cat("  - fragments < 1000:", sum(atac$peak_region_fragments < 1000), "\n")
cat("  - fragments > 20000:", sum(atac$peak_region_fragments > 20000), "\n")
cat("  - pct_reads_in_peaks < 15%:", sum(atac$pct_reads_in_peaks < 15), "\n")
cat("  - TSS.enrichment < 2:", sum(atac$TSS.enrichment < 2), "\n")
cat("  - nucleosome_signal > 4:", sum(atac$nucleosome_signal > 4), "\n")
cat("  - blacklist_ratio > 0.05:", sum(atac$blacklist_ratio > 0.05), "\n")


# 执行过滤
cat("\n执行细胞过滤...\n")
atac_filtered <- subset(
  x = atac,
  subset = peak_region_fragments > 1000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    TSS.enrichment > 2 &
    nucleosome_signal < 4 &
    blacklist_ratio < 0.05
)

cat("✓ 过滤完成!\n")
cat("  过滤后细胞数:", ncol(atac_filtered), "\n")
cat("  过滤后 peaks:", nrow(atac_filtered), "\n\n")

# 保存两个对象
output_dir <- sample_dir
cat("保存对象...\n")

# 保存带 QC 的原始对象
saveRDS(atac, file.path(output_dir, "E145_atac_with_qc.rds"))
cat("  ✓ 原始对象（带 QC）:", file.path(output_dir, "E145_atac_with_qc.rds"), "\n")

# 保存过滤后的对象
saveRDS(atac_filtered, file.path(output_dir, "E145_atac_filtered.rds"))
cat("  ✓ 过滤后对象:", file.path(output_dir, "E145_atac_filtered.rds"), "\n")

cat("\n=== E14.5 样本 QC 完成! ===\n")





# 批量处理剩余 3 个样本
library(Signac)
library(Seurat)
library(GenomicRanges)
library(Matrix)
library(rtracklayer)

# 定义样本信息
samples <- data.frame(
  sample_id = c("E105", "E115", "E125"),
  sample_prefix = c("E10-5WT", "E11-5WT", "E12-5WT"),
  stringsAsFactors = FALSE
)

# 基础路径
base_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate"
gtf_file <- "/disk192/users_dir/buyu/2.参考基因组/3.mmusculus/gencode.vM10.annotation.gtf.gz"

# 读取 GTF（只读一次，所有样本共用）
cat("读取 GTF 文件（只需一次）...\n")
gtf <- rtracklayer::import(gtf_file)
genes <- gtf[gtf$type == "gene"]
mcols(genes)$gene_biotype <- mcols(genes)$gene_type
genome(genes) <- "mm10"
cat("✓ GTF 加载完成，基因数:", length(genes), "\n\n")

# 循环处理每个样本
for (i in 1:nrow(samples)) {
  sample_id <- samples$sample_id[i]
  sample_prefix <- samples$sample_prefix[i]
  
  cat("\n", paste(rep("=", 70), collapse=""), "\n")
  cat("开始处理样本", i, "/3:", sample_id, "(", sample_prefix, ")\n")
  cat(paste(rep("=", 70), collapse=""), "\n\n")
  
  # 设置路径
  sample_dir <- file.path(base_dir, "scATAC", paste0(sample_id, "_palate_scatac"))
  matrix_dir <- file.path(sample_dir, paste0(sample_prefix, "_filtered_peak_bc_matrix"))
  
  # 1. 读取数据
  cat("[1/9] 读取 peak 矩阵...\n")
  counts <- readMM(file.path(matrix_dir, "matrix.mtx"))
  peaks <- read.table(file.path(matrix_dir, "peaks.bed"), header = FALSE, sep = "\t")
  peak_names <- paste0(peaks[,1], ":", peaks[,2], "-", peaks[,3])
  barcodes <- read.table(file.path(matrix_dir, "barcodes.tsv"), header = FALSE, sep = "\t")
  barcode_names <- barcodes[,1]
  rownames(counts) <- peak_names
  colnames(counts) <- barcode_names
  cat("      Peak 矩阵:", nrow(counts), "peaks ×", ncol(counts), "cells\n\n")
  
  # 2. 读取 metadata
  cat("[2/9] 读取 metadata...\n")
  metadata <- read.csv(
    file.path(sample_dir, paste0(sample_prefix, "_singlecell.csv")),
    header = TRUE,
    row.names = 1
  )
  cat("      Metadata:", nrow(metadata), "rows\n\n")
  
  # 3. 创建 ChromatinAssay
  cat("[3/9] 创建 ChromatinAssay...\n")
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = file.path(sample_dir, paste0(sample_prefix, "_fragments.tsv.gz")),
    min.cells = 10,
    min.features = 200
  )
  cat("      过滤后:", nrow(chrom_assay), "peaks ×", ncol(chrom_assay), "cells\n\n")
  
  # 4. 创建 Seurat 对象
  cat("[4/9] 创建 Seurat 对象...\n")
  atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  atac$sample <- sample_id
  cat("      完成\n\n")
  
  # 5. 添加注释
  cat("[5/9] 添加基因组注释...\n")
  Annotation(atac) <- genes
  cat("      完成\n\n")
  
  # 6. 计算核小体信号
  cat("[6/9] 计算核小体信号...\n")
  atac <- NucleosomeSignal(object = atac)
  cat("      完成\n\n")
  
  # 7. 计算 TSS 富集
  cat("[7/9] 计算 TSS 富集分数（需要几分钟）...\n")
  atac <- TSSEnrichment(object = atac, fast = FALSE)
  cat("      完成\n\n")
  
  # 8. 计算额外指标
  cat("[8/9] 计算额外 QC 指标...\n")
  atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
  atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$passed_filters
  
  cat("\n      === QC 指标统计 ===\n")
  cat("      Peak fragments - Median:", round(median(atac$peak_region_fragments), 0), "\n")
  cat("      % Reads in peaks - Median:", round(median(atac$pct_reads_in_peaks), 2), "%\n")
  cat("      TSS enrichment - Median:", round(median(atac$TSS.enrichment), 2), "\n")
  cat("      Nucleosome signal - Median:", round(median(atac$nucleosome_signal), 2), "\n")
  cat("      Blacklist ratio - Median:", round(median(atac$blacklist_ratio), 4), "\n\n")
  
  # 9. 过滤细胞
  cat("[9/9] 过滤细胞...\n")
  n_before <- ncol(atac)
  atac_filtered <- subset(
    x = atac,
    subset = peak_region_fragments > 1000 &
      peak_region_fragments < 20000 &
      pct_reads_in_peaks > 15 &
      TSS.enrichment > 2 &
      nucleosome_signal < 4 &
      blacklist_ratio < 0.05
  )
  n_after <- ncol(atac_filtered)
  cat("      过滤前:", n_before, "cells\n")
  cat("      过滤后:", n_after, "cells\n")
  cat("      保留:", round(n_after/n_before*100, 1), "%\n\n")
  
  # 10. 保存
  cat("保存对象...\n")
  saveRDS(atac, file.path(sample_dir, paste0(sample_id, "_atac_with_qc.rds")))
  saveRDS(atac_filtered, file.path(sample_dir, paste0(sample_id, "_atac_filtered.rds")))
  cat("      ✓ 原始:", paste0(sample_id, "_atac_with_qc.rds\n"))
  cat("      ✓ 过滤:", paste0(sample_id, "_atac_filtered.rds\n"))
  
  cat("\n✓✓✓ 样本", sample_id, "处理完成！\n")
  cat(paste(rep("=", 70), collapse=""), "\n\n")
}

cat("\n", paste(rep("#", 70), collapse=""), "\n")
cat("### 所有样本处理完成！###\n")
cat(paste(rep("#", 70), collapse=""), "\n")


library(Signac)
library(Seurat)

cat("=== Signac 标准分析流程 ===\n\n")

# 读取合并对象（如果需要）
# combined <- readRDS("/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/combined_4timepoints_atac.rds")

cat("[1/7] TF-IDF 标准化...\n")
combined <- RunTFIDF(combined)
cat("  ✓ 完成\n\n")

cat("[2/7] 选择高变 features...\n")
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
cat("  ✓ 完成\n\n")

cat("[3/7] 奇异值分解（SVD/LSI）降维...\n")
cat("      (需要 2-5 分钟)\n")
combined <- RunSVD(combined)
cat("  ✓ 完成\n\n")

# 查看降维结果 - 检查第一个成分（通常与测序深度相关，需要去除）
cat("[4/7] 检查 LSI 成分与测序深度的相关性...\n")
DepthCor(combined)

cat("=== 创建输出文件夹 ===\n\n")

# 创建输出目录
output_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("  ✓ 创建文件夹:", output_dir, "\n\n")
} else {
  cat("  文件夹已存在:", output_dir, "\n\n")
}

cat("[5/7] 非线性降维（UMAP）...\n")
cat("      使用 LSI 2-30 维（排除第一维）\n")
combined <- RunUMAP(
  object = combined,
  reduction = 'lsi',
  dims = 2:30
)
cat("  ✓ UMAP 完成\n\n")

cat("[6/7] 构建 KNN 图和聚类...\n")
combined <- FindNeighbors(
  object = combined,
  reduction = 'lsi',
  dims = 2:30
)
combined <- FindClusters(
  object = combined,
  verbose = FALSE,
  algorithm = 3,
  resolution = 0.5
)
cat("  ✓ 聚类完成\n")
cat("  聚类数:", length(unique(combined$seurat_clusters)), "\n\n")

cat("[7/7] 可视化并保存图片...\n")

# 按样本和聚类查看
library(patchwork)

p1 <- DimPlot(combined, group.by = 'sample', pt.size = 0.1) + 
  ggtitle("按样本分布")

p2 <- DimPlot(combined, label = TRUE, pt.size = 0.1) + 
  ggtitle("聚类结果")

# 保存图片
pdf(file.path(output_dir, "01_UMAP_sample_and_clusters.pdf"), width = 14, height = 6)
print(p1 | p2)
dev.off()

png(file.path(output_dir, "01_UMAP_sample_and_clusters.png"), width = 1400, height = 600, res = 100)
print(p1 | p2)
dev.off()

cat("  ✓ 图片已保存\n")
cat("    PDF:", file.path(output_dir, "01_UMAP_sample_and_clusters.pdf"), "\n")
cat("    PNG:", file.path(output_dir, "01_UMAP_sample_and_clusters.png"), "\n\n")

cat("  查看各样本的聚类分布:\n")
print(table(combined$sample, combined$seurat_clusters))

cat("\n=== 基础分析完成！===\n")

cat("[7/7] 可视化并保存图片...\n")

# 加载必要的包
library(ggplot2)
library(patchwork)

# 按样本和聚类查看
p1 <- DimPlot(combined, group.by = 'sample', pt.size = 0.1) + 
  ggtitle("按样本分布")

p2 <- DimPlot(combined, label = TRUE, pt.size = 0.1) + 
  ggtitle("聚类结果")

# 保存图片
pdf(file.path(output_dir, "01_UMAP_sample_and_clusters.pdf"), width = 14, height = 6)
print(p1 | p2)
dev.off()

png(file.path(output_dir, "01_UMAP_sample_and_clusters.png"), width = 1400, height = 600, res = 100)
print(p1 | p2)
dev.off()

cat("  ✓ 图片已保存\n")
cat("    PDF:", file.path(output_dir, "01_UMAP_sample_and_clusters.pdf"), "\n")
cat("    PNG:", file.path(output_dir, "01_UMAP_sample_and_clusters.png"), "\n\n")

cat("  查看各样本的聚类分布:\n")
print(table(combined$sample, combined$seurat_clusters))

cat("\n各聚类的细胞数:\n")
print(table(combined$seurat_clusters))

cat("\n=== 基础分析完成！===\n")
cat("聚类总数: 19\n")
cat("总细胞数:", ncol(combined), "\n")

cat("=== 聚类统计信息 ===\n\n")

cat("各聚类的细胞数:\n")
cluster_counts <- table(combined$seurat_clusters)
print(cluster_counts)

cat("\n\n各样本在各聚类的分布:\n")
cross_table <- table(combined$sample, combined$seurat_clusters)
print(cross_table)

cat("\n\n各聚类的样本组成比例:\n")
prop_table <- prop.table(cross_table, margin = 2) * 100
print(round(prop_table, 1))

# 计算 gene activity 并加入对象（推荐）
library(Signac)
library(Seurat)

cat("=== 计算 Gene Activity 并加入对象 ===\n")

# 1) 计算 gene activity（可调 promoter upstream 长度，例如 2000 bp）
gene.activities <- GeneActivity(combined, extend.upstream = 2000, extend.downstream = 0)
# gene.activities 是一个稀疏矩阵：基因 × 细胞

# 2) 将其加入为一个新的 assay，命名为 "RNA"（方便与 scRNA 直接整合）
combined[["RNA"]] <- CreateAssayObject(counts = gene.activities)

# 3) 规范化 / 标准化 gene activity（与 scRNA 的处理一致）
combined <- NormalizeData(combined, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 1e4)
combined <- FindVariableFeatures(combined, assay = "RNA", selection.method = "vst", nfeatures = 3000)
combined <- ScaleData(combined, assay = "RNA")

cat("  ✓ Gene activity 已计算并加入为 assay 'RNA'\n")

# 保存对象以便后续使用
saveRDS(combined, file.path(output_dir, "combined_4timepoints_atac_withRNAactivity.rds"))
cat("  ✓ 已保存: combined_4timepoints_atac_withRNAactivity.rds\n")



# ============================================
# scRNA-seq 和 scATAC-seq 整合分析
# ============================================

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)

cat("=== scRNA 和 scATAC 整合分析 ===\n\n")

# 1. 加载数据
cat("[1/7] 加载数据...\n")
scRNA_path <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scRNA_rds/注释过RNA/integrated_annotation/scrna_integrated_annotation.rds"
scATAC_path <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/results/combined_4timepoints_atac_withRNAactivity.rds"

scRNA <- readRDS(scRNA_path)
scATAC <- readRDS(scATAC_path)

cat("  ✓ scRNA:", ncol(scRNA), "cells,", nrow(scRNA), "genes\n")
cat("  ✓ scATAC:", ncol(scATAC), "cells,", nrow(scATAC), "peaks\n\n")

# 2. 准备 scRNA 数据
cat("[2/7] 准备 scRNA 参考数据...\n")
DefaultAssay(scRNA) <- "RNA"

# 确保 scRNA 有必要的处理
if (!"data" %in% slotNames(scRNA@assays$RNA) || length(scRNA@assays$RNA@data) == 0) {
  cat("  标准化 scRNA...\n")
  scRNA <- NormalizeData(scRNA)
}

if (!"scale.data" %in% slotNames(scRNA@assays$RNA) || length(scRNA@assays$RNA@scale.data) == 0) {
  cat("  缩放 scRNA...\n")
  scRNA <- FindVariableFeatures(scRNA, nfeatures = 3000)
  scRNA <- ScaleData(scRNA)
}

if (!"pca" %in% names(scRNA@reductions)) {
  cat("  运行 PCA...\n")
  scRNA <- RunPCA(scRNA, npcs = 50)
}

cat("  ✓ scRNA 准备完成\n\n")

# 3. 准备 scATAC 数据
cat("[3/7] 准备 scATAC 数据...\n")
DefaultAssay(scATAC) <- "RNA"  # 使用 Gene Activity

cat("  ✓ scATAC Gene Activity 已就绪\n\n")

# 4. 寻找整合锚点
cat("[4/7] 寻找整合锚点（Transfer Anchors）...\n")
cat("  这可能需要 10-20 分钟...\n\n")

transfer.anchors <- FindTransferAnchors(
  reference = scRNA,
  query = scATAC,
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = "cca",
  dims = 1:30
)

cat("  ✓ 找到", length(transfer.anchors@anchors), "个锚点\n\n")

# 5. 转移细胞类型注释
cat("[5/7] 转移细胞类型标签...\n")

# 转移 cell_type 注释
celltype.predictions <- TransferData(
  anchorset = transfer.anchors,
  refdata = scRNA$cell_type,
  weight.reduction = scATAC[["lsi"]],
  dims = 2:30
)

scATAC <- AddMetaData(scATAC, metadata = celltype.predictions)

cat("  ✓ 细胞类型标签转移完成\n\n")

# 6. 查看转移结果
cat("=== 转移结果统计 ===\n\n")

cat("各细胞类型的细胞数:\n")
print(table(scATAC$predicted.id))

cat("\n\n预测得分分布:\n")
print(summary(scATAC$prediction.score.max))

cat("\n\n各时间点的细胞类型分布:\n")
print(table(scATAC$sample, scATAC$predicted.id))

# 7. 保存结果
cat("\n[6/7] 保存整合后的对象...\n")
output_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/results"
saveRDS(scATAC, file.path(output_dir, "scATAC_with_celltype_annotation.rds"))
cat("  ✓ 已保存:", file.path(output_dir, "scATAC_with_celltype_annotation.rds"), "\n\n")

cat("=== 整合完成！===\n")

# 可视化：按细胞类型着色的 UMAP
library(ggplot2)

p1 <- DimPlot(scATAC, group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("scATAC: Predicted Cell Types")

# 保存图片
pdf(file.path(output_dir, "06_celltype_annotation_UMAP.pdf"), width = 12, height = 8)
print(p1)
dev.off()

png(file.path(output_dir, "06_celltype_annotation_UMAP.png"), width = 1200, height = 800, res = 100)
print(p1)
dev.off()

cat("✓ 图片已保存\n")

# 可视化预测得分
p2 <- FeaturePlot(scATAC, features = "prediction.score.max", pt.size = 0.1) +
  ggtitle("Prediction Score") +
  scale_color_viridis_c()

# 保存图片
pdf(file.path(output_dir, "07_prediction_score.pdf"), width = 10, height = 8)
print(p2)
dev.off()

png(file.path(output_dir, "07_prediction_score.png"), width = 1000, height = 800, res = 100)
print(p2)
dev.off()

cat("✓ 预测得分图已保存\n")

# 对比原始聚类和预测细胞类型
library(patchwork)

p_cluster <- DimPlot(scATAC, group.by = "seurat_clusters", label = TRUE, pt.size = 0.1) +
  ggtitle("Original Clusters")

p_celltype <- DimPlot(scATAC, group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 0.1) +
  ggtitle("Predicted Cell Types")

# 并排显示
p_compare <- p_cluster | p_celltype

# 保存
pdf(file.path(output_dir, "08_cluster_vs_celltype.pdf"), width = 16, height = 7)
print(p_compare)
dev.off()

png(file.path(output_dir, "08_cluster_vs_celltype.png"), width = 1600, height = 700, res = 100)
print(p_compare)
dev.off()

cat("✓ 对比图已保存\n")




# 基于细胞类型的差异可及性分析
cat("=== 基于细胞类型的差异 Peak 分析 ===\n\n")

# 切换到 peaks assay
DefaultAssay(scATAC) <- "peaks"

# 设置细胞类型为 Idents
Idents(scATAC) <- "predicted.id"

cat("开始分析...\n")
cat("这可能需要 15-30 分钟\n\n")

# 寻找每个细胞类型的 marker peaks
da_peaks_celltype <- FindAllMarkers(
  object = scATAC,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments',
  min.pct = 0.05,
  logfc.threshold = 0.2
)

cat("✓ 分析完成！\n")
cat("发现", nrow(da_peaks_celltype), "个差异 peaks\n")


# 查看各细胞类型的 marker peaks 数量
cat("=== 各细胞类型的 marker peaks 数量 ===\n")
print(table(da_peaks_celltype$cluster))

cat("\n=== Top 10 marker peaks (按 log2FC 排序) ===\n")
print(head(da_peaks_celltype[order(-da_peaks_celltype$avg_log2FC), ], 10))


# 保存差异 peaks 结果
output_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/results"

# 保存完整结果
write.csv(da_peaks_celltype, 
          file.path(output_dir, "DA_peaks_by_celltype_all.csv"), 
          row.names = FALSE)

# 保存每个细胞类型的 top 100 peaks
library(dplyr)
top100_peaks <- da_peaks_celltype %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top100_peaks, 
          file.path(output_dir, "DA_peaks_by_celltype_top100.csv"), 
          row.names = FALSE)

cat("✓ 结果已保存\n")
cat("  完整结果:", file.path(output_dir, "DA_peaks_by_celltype_all.csv"), "\n")
cat("  Top100:", file.path(output_dir, "DA_peaks_by_celltype_top100.csv"), "\n")

> # 查看特定细胞类型的 top peaks 附近的基因
library(Signac)

# 以 Epithelial_cells 为例
epi_peaks <- da_peaks_celltype %>%
  filter(cluster == "Epithelial_cells") %>%
  arrange(desc(avg_log2FC)) %>%
  head(5)

cat("上皮细胞 Top 5 marker peaks:\n")
print(epi_peaks)
上皮细胞 Top 5 marker peaks:
                                 p_val avg_log2FC pct.1 pct.2     p_val_adj
chr1-10637014-10637641   6.398023e-145   6.006255 0.056 0.001 1.293891e-139
chr7-130340034-130341279  0.000000e+00   5.770024 0.150 0.003  0.000000e+00
chr9-77586957-77587577   8.229859e-134   5.553234 0.056 0.001 1.664349e-128
chr3-121897823-121898831 1.997638e-308   5.543650 0.126 0.003 4.039884e-303
chr3-38585454-38586073   1.751019e-139   5.515853 0.057 0.001 3.541138e-134
                                  cluster                     gene
chr1-10637014-10637641   Epithelial_cells   chr1-10637014-10637641
chr7-130340034-130341279 Epithelial_cells chr7-130340034-130341279
chr9-77586957-77587577   Epithelial_cells   chr9-77586957-77587577
chr3-121897823-121898831 Epithelial_cells chr3-121897823-121898831
chr3-38585454-38586073   Epithelial_cells   chr3-38585454-38586073
> 

# 加载包
library(Rsamtools)

cat("=== 使用 Rsamtools 重建索引 ===\n\n")

# 定义 fragment 文件路径
fragment_files <- c(
  "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/E105_palate_scatac/E10-5WT_fragments.tsv.gz",
  "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/E115_palate_scatac/E11-5WT_fragments.tsv.gz",
  "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/E125_palate_scatac/E12-5WT_fragments.tsv.gz",
  "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/E145_palate_scatac/E14-5WT_fragments.tsv.gz"
)

# 逐个重建索引
for (frag_path in fragment_files) {
  cat("正在为", basename(frag_path), "重建索引...\n")
  indexTabix(file = frag_path, format = "bed")
  cat("  ✓ 完成\n")
}

cat("\n✓ 所有索引重建完成！\n")


cat("=== 移除 NA 细胞 ===\n\n")

# 查看总细胞数
cat("原始细胞数:", ncol(scATAC), "\n")
cat("NA 细胞数:", sum(is.na(scATAC$predicted.id)), "\n")

# 创建去除 NA 的子集
scATAC_clean <- subset(scATAC, subset = predicted.id != "NA")

cat("清理后细胞数:", ncol(scATAC_clean), "\n")

# 检查是否还有 NA
cat("剩余 NA 数:", sum(is.na(scATAC_clean$predicted.id)), "\n")

cat("\n✓ 清理完成\n")

# 使用你已有的 GTF 文件
cat("=== 加载完整 GTF 注释 ===\n\n")

conda activate signac

library(rtracklayer)
library(GenomicRanges)

# 使用 GRCm39 GTF
gtf_path <- "/disk192/users_dir/buyu/2.参考基因组/3.mmusculus/Mus_musculus.GRCm39.115.chr.gtf.gz"
cat("加载 GTF:", gtf_path, "\n")

# 读取完整注释
full_annotation <- rtracklayer::import(gtf_path)

cat("原始注释记录数:", length(full_annotation), "\n")

# 查看包含哪些类型
cat("\n注释类型统计:\n")
print(table(full_annotation$type))

# 只保留需要的类型
annotation_filtered <- full_annotation[full_annotation$type %in% c("gene", "transcript", "exon")]

cat("\n过滤后记录数:", length(annotation_filtered), "\n")

# 检查并转换染色体命名
cat("染色体风格:", seqlevelsStyle(annotation_filtered)[1], "\n")

if (seqlevelsStyle(annotation_filtered)[1] != "UCSC") {
  cat("转换为 UCSC 风格...\n")
  seqlevelsStyle(annotation_filtered) <- "UCSC"
}

cat("转换后染色体:", paste(head(seqlevels(annotation_filtered), 10), collapse = ", "), "\n")

# 更新 Seurat 对象的注释
Annotation(scATAC_clean) <- annotation_filtered

cat("\n✓ 完整注释已添加到对象\n")
# 保存清理后的对象
cat("=== 保存清理后的 scATAC 对象 ===\n\n")

output_file <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/scATAC_clean_annotated.rds"

cat("保存到:", output_file, "\n")

saveRDS(scATAC_clean, file = output_file)

cat("✓ 保存完成\n")
cat("文件大小:", file.size(output_file) / 1024^2, "MB\n")

# ===================================================================
# 12. Motif Enrichment Analysis - 找到调控转录因子
# ===================================================================

cat("\n=== 12. Motif Enrichment Analysis ===\n\n")

# 加载所有必要的包
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)

# 设置输出目录
output_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/results"

# 确保使用 peaks assay
DefaultAssay(scATAC_clean) <- "peaks"

# 获取 JASPAR motif 数据库
cat("从 JASPAR2020 获取转录因子 motif...\n")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = FALSE)  # 10090 = 小鼠
)

cat("获得", length(pfm), "个转录因子 motif\n\n")

# 添加 motif 信息到对象（这步需要几分钟）
cat("添加 motif 到 Seurat 对象...\n")
cat("注意：这一步可能需要 5-10 分钟，请耐心等待...\n")

scATAC_clean <- AddMotifs(
  object = scATAC_clean,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

cat("\n✓ Motif 信息已添加！\n\n")

# 保存添加了 motif 的对象
cat("保存添加了 motif 信息的对象...\n")
saveRDS(scATAC_clean, 
        file = "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/scATAC_clean_with_motifs.rds")
cat("✓ 已保存\n")

# 先加载对象
cat("=== 加载对象 ===\n\n")

library(Signac)
library(Seurat)

cat("加载清理后的对象...\n")
scATAC_clean <- readRDS("/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/scATAC_clean_annotated.rds")

cat("✓ 对象加载完成\n")
cat("细胞数:", ncol(scATAC_clean), "\n")
cat("Peak 数:", nrow(scATAC_clean[["peaks"]]), "\n\n")

# 然后运行 motif 分析
cat("=== 开始 Motif 分析 ===\n\n")

library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(motifmatchr)

# 设置输出目录
output_dir <- "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/results"

# 确保使用 peaks assay
DefaultAssay(scATAC_clean) <- "peaks"

# 获取 JASPAR motif 数据库
cat("从 JASPAR2020 获取转录因子 motif...\n")
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 10090, all_versions = FALSE)
)

cat("获得", length(pfm), "个转录因子 motif\n\n")

# 添加 motif 信息
cat("添加 motif 到对象（需要 5-10 分钟）...\n")
scATAC_clean <- AddMotifs(
  object = scATAC_clean,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

cat("\n✓ Motif 信息已添加！\n\n")

# 保存
saveRDS(scATAC_clean, 
        file = "/disk192/users_dir/buyu/1.布宇/supplement/26/Mouse_palate/scATAC/scATAC_clean_with_motifs.rds")
cat("✓ 已保存\n")

# ===================================================================
# 13. 上皮细胞 Motif Enrichment
# ===================================================================

cat("\n=== 13. 上皮细胞特异性 Motif Enrichment ===\n\n")

# 设置 ident 为细胞类型
Idents(scATAC_clean) <- "predicted.id"

# 对上皮细胞 vs 其他细胞做 motif enrichment
cat("分析上皮细胞特异性转录因子...\n")
cat("注意：这一步也需要几分钟...\n\n")

enriched_motifs <- FindMarkers(
  object = scATAC_clean,
  ident.1 = "Epithelial_cells",
  only.pos = TRUE,
  test.use = 'LR',  # Logistic Regression
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

cat("✓ 找到", nrow(enriched_motifs), "个显著富集的 motif\n\n")

# 查看 top 结果
cat("Top 10 富集的转录因子 motif:\n")
print(head(enriched_motifs, 10))

# 保存结果
write.csv(enriched_motifs, 
          file = file.path(output_dir, "12_epithelial_enriched_motifs.csv"))

cat("\n✓ 结果已保存\n")


