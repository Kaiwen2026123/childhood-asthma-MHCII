################################################################################
## 02_single_cell.R
## 单细胞分析全流程: 质量控制 -> 聚类注释 -> scDEG -> 交集 -> 模块评分
##              -> pseudobulk (limma-voom) -> CellChat (pooled + per-donor)
## 对应论文 Methods 2.8-2.11, Results 3.5-3.6
## 注意: 对象路径与队列路径按实际环境修改
################################################################################
library(Seurat); library(dplyr); library(tidyr); library(ggplot2)
library(edgeR); library(limma); library(pROC); library(clusterProfiler)

sc <- readRDS("/home/share/pbmc/sc_clean_with_mono_subtypes.rds")
sc <- sc[!grepl("^ENS", rownames(sc), ignore.case = TRUE), ]   # 去除无注释转录本
## QC阈值 (实际流程, 见Methods 2.8): UMI >= 1000 + 基因数/线粒体严格过滤 + 双细胞去除
best_genes <- c("HLA-DPB1","HLA-DOA","CIITA","CCL5","HLA-DPA1","CD3G","CD2")
Idents(sc) <- "celltype_clean"

## ---- 2.1 scDEG 全局 + 细胞级 (Methods 2.9) ----
Idents(sc) <- "Group"
sc_deg_official <- FindMarkers(sc, ident.1 = "Asthma", ident.2 = "Healthy",
                               logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox")
sc_deg_final <- sc_deg_official[!grepl("^ENS", rownames(sc_deg_official)), ]
## 细胞级 (CD4+ T细胞群因组间极度不平衡57/1720排除)
cell_deg_list <- list()
for (ct in sort(unique(sc$celltype_clean))){
  if (ct == "CD4+ T cells") next
  sub <- subset(sc, subset = celltype_clean == ct)
  if (min(table(sub$Group)) < 50) next
  Idents(sub) <- "Group"
  d <- tryCatch(FindMarkers(sub, ident.1 = "Asthma", ident.2 = "Healthy",
                            logfc.threshold = 0, min.pct = 0.05, test.use = "wilcox"),
                error = function(e) NULL)
  if (!is.null(d)){ d$gene <- rownames(d); d$celltype <- ct
    cell_deg_list[[ct]] <- d }
}
cell_deg <- bind_rows(cell_deg_list)

## ---- 2.2 bulk∩sc 交集 (11候选与scDEG双向上/下调) ----
## (在bulk会话跑; 见01/03脚本)

## ---- 2.3 单细胞模块评分 (Fig 6G/H/I) ----
Idents(sc) <- "celltype_clean"
t2_genes_sc <- intersect(c("IL33","CCL17","CCL26","SERPINB2","FCER1A","IL4R","CPA3","HPGDS"),
                         rownames(sc))
panel_sc <- intersect(best_genes, rownames(sc))
sc <- AddModuleScore(sc, features = list(t2_genes_sc), name = "T2_score")
sc <- AddModuleScore(sc, features = list(panel_sc), name = "Panel_score")
## C2签名 (bulk top50, 见04脚本生成) 同法投影

## ---- 2.4 donor-level pseudobulk, limma-voom (Methods 2.9 修订版) ----
genes_check <- unique(c(best_genes, "HLA-DRB1","HLA-DQA1","HLA-DQB1","CTSC","DPP4"))
av_ct <- AggregateExpression(sc, assays = "RNA", features = intersect(genes_check, rownames(sc)),
                             group.by = c("celltype_clean","orig.ident"), slot = "counts")$RNA
key_ct <- c("Cytotoxic CD8+ T cells","NK cells","Classical Monocytes","Naive CD4+ T cells")
pb_de <- lapply(key_ct, function(ct2){
  cols <- grep(paste0(ct2, "_"), colnames(av_ct), fixed = TRUE, value = TRUE)
  samps <- sub(paste0(ct2, "_"), "", cols, fixed = TRUE)
  grp <- factor(sc$Group[match(samps, sc$orig.ident)], levels = c("Healthy","Asthma"))
  dge <- DGEList(av_ct[, cols, drop = FALSE], group = grp)
  dge <- dge[filterByExpr(dge, group = grp, min.count = 5), , keep.lib.size = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, model.matrix(~ grp))
  fit <- eBayes(lmFit(v, model.matrix(~ grp)))
  tt <- topTable(fit, coef = 2, n = Inf); tt$gene <- rownames(tt); tt$celltype <- ct2; tt
})
pb_de <- bind_rows(pb_de)
write.csv(pb_de, "TableS2_pseudobulk_voom.csv", row.names = FALSE)

## ---- 2.5 CellChat pooled (Methods 2.10) ----
library(CellChat)
set.seed(123)
sc_ds <- subset(sc, downsample = 2000)
mk_cc <- function(obj){
  cc <- createCellChat(obj, group.by = "celltype_clean", assay = "RNA")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc) %>% identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>% aggregateNet()
  cc
}
cc_asthma <- mk_cc(subset(sc_ds, subset = Group == "Asthma"))
cc_healthy <- mk_cc(subset(sc_ds, subset = Group == "Healthy"))
save(cc_asthma, cc_healthy, file = "CellChat_objects.RData")

## ---- 2.6 CellChat per-donor (Fig S16) ----
donor_flow <- sapply(sort(unique(sc_ds$orig.ident)), function(d){
  sub <- subset(sc_ds, subset = orig.ident == d)
  cc <- mk_cc(sub)
  sum(cc@netP$prob[, , "MHC-II"], na.rm = TRUE)
})
save(donor_flow, file = "CellChat_per_donor_MHC2flow.RData")

sessionInfo()
