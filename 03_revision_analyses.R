################################################################################
## 03_revision_analyses.R
## 审稿修订新增分析 (Response to Reviewers):
##   R1-C1: 训练集独立筛选重现 + 嵌套CV
##   R1-C2: 完整pseudobulk (见02脚本2.4节)
##   R1-C3: 聚类稳定性 (silhouette + bootstrap) + GSE27011外部重现 + 敏感性聚类
##   R1-C4: per-donor CellChat (见02脚本2.6节)
##   R2-C9: MHC-II双模式定量 (19/19, 10/10)
##   R2-C10: GSE40888未刺激亚组重验证 + GSE65204 IgE相关
################################################################################
library(pROC); library(dplyr); library(ConsensusClusterPlus); library(cluster)
library(GSVA); library(ggplot2)

load("ML_input.RData"); load("ML_result.RData")
load("combat_data.RData"); load("sig.RData")
load("WGCNA_result.RData"); load("subtype.RData")
trainX <- tr$X; trainy <- tr$y
ca_samples <- subtype_df$sample

## ===== R1-C1a: 训练集独立筛选 (预设规则0.65 -> 原7基因重现) =====
auc_tr <- sapply(colnames(trainX), function(g)
  as.numeric(auc(roc(trainy, trainX[[g]], quiet = TRUE))))
print(round(sort(auc_tr, decreasing = TRUE), 3))

## ===== R1-C1b: 嵌套交叉验证 =====
set.seed(123)
folds <- caret::createFolds(factor(ifelse(trainy == 1, "CA", "Ctrl")), k = 10)
nested <- sapply(folds, function(te){
  tr_i <- setdiff(seq_along(trainy), te)
  a <- sapply(colnames(trainX), function(g)
    as.numeric(auc(roc(trainy[tr_i], trainX[tr_i, g], quiet = TRUE))))
  g <- names(a)[a >= 0.65]
  if (length(g) < 3) g <- names(sort(a, decreasing = TRUE))[1:7]
  fit <- glm(trainy[tr_i] ~ ., trainX[tr_i, g, drop = FALSE], family = binomial)
  as.numeric(auc(roc(trainy[te], predict(fit, newdata = trainX[te, g, drop = FALSE],
                                         type = "response"), quiet = TRUE)))
})
cat("Nested CV AUC =", round(mean(nested), 3), "+/-", round(sd(nested), 3), "\n")
save(nested, file = "nested_CV_result.RData")

## ===== R1-C3a: 聚类稳定性 =====
sub_input <- combat_data[intersect_genes, ca_samples]
sil <- silhouette(as.numeric(factor(subtype_df$subtype)), dist(t(sub_input)))
cat("Silhouette =", round(mean(sil[, 3]), 3), "\n")
ari_fun <- function(x, y){
  tab <- table(x, y); n <- sum(tab)
  a <- sum(choose(tab, 2)); b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(n, 2) - a - b - c
  (a - (a+b)*(a+c)/choose(n, 2)) / ((a+b+a+c)/2 - (a+b)*(a+c)/choose(n, 2))
}
set.seed(123)
boot_ari <- replicate(100, {
  idx <- sample(ncol(sub_input), round(0.8 * ncol(sub_input)))
  cl <- pam(t(sub_input[, idx]), k = 2)$clustering
  ari_fun(cl, as.numeric(factor(subtype_df$subtype[idx])))
})
cat("Bootstrap ARI median =", round(median(boot_ari), 3), "\n")
save(sil, boot_ari, file = "cluster_stability.RData")

## ===== R1-C3b: GSE27011外部重现 =====
load("validation_cohorts.RData")
g11 <- intersect(intersect_genes, rownames(r27011$exp))
ca270 <- r27011$cl$sample[r27011$cl$group == "CA"]
set.seed(123)
cc270 <- ConsensusClusterPlus(as.matrix(r27011$exp[g11, ca270]), maxK = 4, reps = 100,
                              pItem = 0.8, clusterAlg = "pam", distance = "euclidean",
                              title = "CC_GSE27011", plot = "png")
cl270 <- cc270[[2]]$consensusClass
save(cl270, file = "GSE27011_clusters.RData")

## ===== R1-C3c: 敏感性聚类 (排除7 signature基因) =====
deg_noPanel <- setdiff(rownames(deg_strict), best_genes)
set.seed(123)
ccp_sens <- ConsensusClusterPlus(as.matrix(combat_data[deg_noPanel, ca_samples]),
                                 maxK = 4, reps = 100, pItem = 0.8,
                                 clusterAlg = "pam", distance = "euclidean",
                                 title = "CC_sensitivity", plot = "png")
cat("Sensitivity ARI =", round(ari_fun(ccp_sens[[2]]$consensusClass,
      as.numeric(factor(subtype_df$subtype))), 3), "\n")

## ===== R2-C9: MHC-II双模式定量 (单细胞水平, exploratory) =====
cell_deg <- read.csv("sc_deg_by_celltype.csv")
hla6 <- c("HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPB1","HLA-DPA1","CIITA")
sub <- cell_deg %>% filter(gene %in% hla6, p_val_adj < 0.05)
t_types <- c("Naive CD4+ T cells","CD4+ T cells","Th cells","Th17 cells",
             "CD8+ T cells","Cytotoxic CD8+ T cells","Regulatory T cell","Effector memory T cells")
myeloid <- c("Classical Monocytes","Non classical Monocytes","NK cells","Cytotoxic CD8+ T cells")
drdq_t <- sub %>% filter(gene %in% c("HLA-DRB1","HLA-DQA1","HLA-DQB1"), celltype %in% t_types)
dp_m   <- sub %>% filter(gene %in% c("HLA-DPB1","HLA-DPA1","CIITA"), celltype %in% myeloid)
cat("DR/DQ down in", sum(drdq_t$avg_log2FC < 0), "/", nrow(drdq_t),
    "median", round(median(drdq_t$avg_log2FC), 2), "\n")
cat("DP/CIITA up in", sum(dp_m$avg_log2FC > 0), "/", nrow(dp_m),
    "median", round(median(dp_m$avg_log2FC), 2), "\n")

## ===== R2-C10a: GSE40888未刺激亚组重验证 =====
con <- gzfile("GSE40888_series_matrix.txt.gz")
lines <- readLines(con, n = 3000); close(con)
gsm  <- gsub('"', "", strsplit(lines[grepl("^!Sample_geo_accession", lines)], "\t")[[1]])[-1]
stim <- gsub('"', "", strsplit(lines[grepl('^!Sample_characteristics_ch1\t"stimulation', lines)], "\t")[[1]][-1])
names(stim) <- gsm
keep <- names(stim)[grepl("no stimulation", stim)]
cl_m <- r40888$cl[r40888$cl$sample %in% keep, ]
g <- intersect(best_genes, rownames(r40888$exp))
X <- as.data.frame(t(r40888$exp[g, cl_m$sample, drop = FALSE])); X[] <- lapply(X, as.numeric)
y <- ifelse(cl_m$group == "CA", 1, 0)
fit_best <- glm(trainy ~ ., trainX[, best_genes, drop = FALSE], family = binomial)
p_m <- predict(fit_best, newdata = X[, g, drop = FALSE], type = "response")
cat("GSE40888 unstimulated AUC =",
    round(as.numeric(auc(roc(y, p_m, quiet = TRUE))), 3), "\n")

## ===== R2-C10b: GSE65204 IgE相关 (exploratory) =====
con <- gzfile("GSE65204_series_matrix.txt.gz")
lines <- readLines(con, n = 3000); close(con)
gsm <- gsub('"', "", strsplit(lines[grepl("^!Sample_geo_accession", lines)], "\t")[[1]])[-1]
ige <- as.numeric(gsub("lnige: ", "", gsub('"', "",
  strsplit(lines[grepl('^!Sample_characteristics_ch1\t"lnige', lines)], "\t")[[1]][-1])))
ast <- gsub("asthma: ", "", gsub('"', "",
  strsplit(lines[grepl('^!Sample_characteristics_ch1\t"asthma', lines)], "\t")[[1]][-1]))
names(ige) <- gsm
ca652 <- intersect(names(ige)[ast == "TRUE" & !is.na(ige)],
                   Group$sample[Group$group == "CA"])
ct <- cor.test(colMeans(combat_data[best_genes, ca652]), ige[ca652])
cat("IgE vs panel score: r =", round(ct$estimate, 3),
    "p =", format(ct$p.value, digits = 3), "n =", length(ca652), "\n")

sessionInfo()

## ===== 补充: T2 / MHC-II ssGSEA 评分 (Methods 2.7, Fig 8) =====
library(GSVA); library(org.Hs.eg.db); library(AnnotationDbi)
t2_genes <- intersect(c("IL33","CCL17","CCL26","SERPINB2","FCER1A","IL4R","CPA3","HPGDS"),
                      rownames(combat_data))
t2_mat <- gsva(gsvaParam(as.matrix(combat_data), list(T2 = t2_genes)), verbose = FALSE)
go_df <- AnnotationDbi::select(org.Hs.eg.db, keys = "GO:0042590",
                               keytype = "GOALL", columns = "SYMBOL")
mhc2_genes <- intersect(unique(go_df$SYMBOL), rownames(combat_data))
mhc2_mat <- gsva(gsvaParam(as.matrix(combat_data), list(MHC2 = mhc2_genes)),
                 verbose = FALSE)
panel_score <- colMeans(combat_data[best_genes, ca_samples])
save(t2_mat, mhc2_mat, file = "bulk_scores.RData")

## ===== 补充: CIITA 主角分析 (Results 3.4) =====
ciita <- as.numeric(combat_data["CIITA", ca_samples])
cor.test(ciita, panel_score)                       # r = 0.722
cor.test(ciita, as.numeric(mhc2_mat[1, ca_samples]))  # r = 0.680
