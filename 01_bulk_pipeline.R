################################################################################
## 复刻 Li et al. 2026 (Artif Cells Nanomed Biotechnol 54:358-373) —— bulk 全流程
## 疾病: 儿童哮喘(CA) vs Control
## Training : GSE65204 + GSE103166 (ComBat 去批次合并)
## Validation: GSE19187 (独立外部验证, 不参与去批次)
## 流程: 数据处理 -> 去批次 -> DEG -> WGCNA -> DEG∩module -> 机器学习组合筛选
##       -> 外部验证 -> 一致性聚类分型 -> CIBERSORT
################################################################################
rm(list = ls())
dir.create("/home/u25034101613/9.21分析", showWarnings = FALSE, recursive = TRUE)
setwd("/home/u25034101613/9.21分析")

## ======================= 00. 环境准备 ========================================
cran_pkgs  <- c("tidyverse","dplyr","stringr","data.table","pheatmap","RColorBrewer",
                "glmnet","caret","mboost","plsRglm","randomForest","xgboost",
                "pROC","scatterplot3d","factoextra","FactoMineR")
bioc_pkgs  <- c("GEOquery","limma","sva","clusterProfiler","enrichplot",
                "org.Hs.eg.db","ConsensusClusterPlus","GSVA")
## 首次运行取消下面两行注释:
# install.packages(setdiff(cran_pkgs, rownames(installed.packages())))
# BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())))
invisible(lapply(c(cran_pkgs, bioc_pkgs), library, character.only = TRUE))

## ======================= 01. 数据下载与处理 ==================================
## 所有GEO数据(series matrix + GPL注释)均由 getGEO() 自动下载到本目录
## 仅 CIBERSORT.R + LM22.txt 需官网注册下载 (见第06步)

process_gse <- function(expr, fdata, id_col, gene_col, pdata, phe_cols, group_fn, outfile){
  genes <- fdata[, c(id_col, gene_col)]
  expr  <- data.frame(expr)
  expr  <- merge(expr, genes, by.x = 0, by.y = 1)
  ## 重复探针取均值最大者
  rowM <- apply(expr, 1, function(x) mean(as.numeric(x), na.rm = TRUE))
  expr <- expr[order(rowM, decreasing = TRUE), ]
  expr <- expr[!duplicated(expr[, ncol(expr)]), ]
  expr <- na.omit(expr)
  expr <- expr[expr[, ncol(expr)] != "", ]
  rownames(expr) <- expr[, ncol(expr)]
  expr <- expr[, -c(1, ncol(expr))]
  cl <- pdata[, phe_cols]
  colnames(cl) <- c("sample", "group")
  cl <- group_fn(cl)
  expr <- expr[, cl$sample]
  save(expr, cl, file = outfile)
  list(expr = expr, cl = cl)
}

## 读取GEO平台注释文件 (.annot.gz 或 _platform.soft), 自动找基因列
read_gpl <- function(file){
  con <- if (endsWith(file, ".gz")) gzfile(file) else file(file)
  L <- readLines(con); close(con)
  i <- which(startsWith(L, "ID"))[1]   # 表格表头行(以ID开头)
  df <- read.delim(text = paste(L[i:length(L)], collapse = "\n"),
                   check.names = FALSE, stringsAsFactors = FALSE)
  gc <- grep("symbol|gene_assignment", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  message("使用基因列: ", gc)
  df
}

## ---- 1.1 GSE103166 ----
eSet  <- getGEO("GSE103166", destdir = ".", getGPL = TRUE)
expr0 <- exprs(eSet[[1]])
fd0   <- fData(eSet[[1]])
pd0   <- pData(eSet[[1]])
r103166 <- process_gse(expr0, fd0, "ID", "ORF", pd0, c(2, 33),
  function(cl){ cl <- cl[cl$group != "Convalescent", ]
                cl$group <- ifelse(cl$group == "Control", "Control", "CA"); cl },
  "GSE103166_exp_clin.Rdata")
GSE103166_exp  <- r103166$expr
GSE103166_clin <- r103166$cl

## ---- 1.2 GSE65204 (series matrix和GPL均自动下载到工作目录) ----
eSet  <- getGEO("GSE65204", destdir = ".", getGPL = FALSE)[[1]]
expr0 <- exprs(eSet)
pd0   <- pData(eSet)
fd0   <- read_gpl("GPL14550_platform.soft")
r65204 <- process_gse(expr0, fd0, "ID", "GENE_SYMBOL", pd0, c(2, 42),
  function(cl){ cl$group <- ifelse(cl$group == "FALSE", "Control", "CA"); cl },
  "GSE65204_exp_clin.Rdata")
GSE65204_exp  <- r65204$expr
GSE65204_clin <- r65204$cl

## ---- 1.3 GSE19187 (验证集; GPL6244已是log2, 不要再log2转换) ----
eSet  <- getGEO("GSE19187", destdir = ".", getGPL = FALSE)[[1]]
expr0 <- exprs(eSet)
print(range(expr0, na.rm = TRUE))   ## 应在2~14之间
pd0   <- pData(eSet)
fd0   <- read_gpl("GPL6244_platform.soft")
fd0$gene_assignment <- trimws(str_split_fixed(fd0$gene_assignment, "//", 3)[, 2])
r19187 <- process_gse(expr0, fd0, "ID", "gene_assignment", pd0, c(2, 35),
  function(cl){ cl <- cl[!cl$group %in% c("Rhinitis","Rhinitis with controlled asthma"), ]
                cl$group <- ifelse(cl$group == "Healthy", "Control", "CA"); cl },
  "GSE19187_exp_clin.Rdata")
GSE19187_exp  <- r19187$expr
GSE19187_clin <- r19187$cl

## ======================= 02. 合并训练集 + ComBat 去批次 ======================
sameid <- intersect(rownames(GSE65204_exp), rownames(GSE103166_exp))
bind_data <- cbind(GSE65204_exp[sameid, ], GSE103166_exp[sameid, ])
batch_info <- rbind(data.frame(ID = GSE65204_clin$sample,  batch = "GSE65204"),
                    data.frame(ID = GSE103166_clin$sample, batch = "GSE103166"))
batch <- ifelse(batch_info$batch == "GSE65204", 1, 2)
combat_data <- ComBat(dat = bind_data, batch = batch, mod = NULL,
                      par.prior = TRUE, prior.plots = FALSE)
combat_data <- normalizeBetweenArrays(combat_data)
Group <- rbind(data.frame(sample = GSE65204_clin$sample,  group = GSE65204_clin$group),
               data.frame(sample = GSE103166_clin$sample, group = GSE103166_clin$group))
combat_data <- combat_data[, Group$sample]
table(Group$group)
save(combat_data, Group, file = "combat_data.RData")

## 去批次前后箱线图 + PCA
pdf("去批次前箱线图.pdf", width = 10, height = 8)
boxplot(bind_data, outline = FALSE, notch = FALSE, las = 2, xaxt = "n"); dev.off()
pdf("去批次后箱线图.pdf", width = 10, height = 8)
boxplot(combat_data, outline = FALSE, notch = FALSE, las = 2, xaxt = "n"); dev.off()
library(FactoMineR); library(factoextra)
pca_check <- function(mat, col_vec, ttl){
  p <- PCA(as.data.frame(t(mat)), graph = FALSE)
  fviz_pca_ind(p, geom.ind = "point", col.ind = col_vec,
               palette = c("#3498DB","#E74C3C"), addEllipses = TRUE,
               legend.title = "Batch", title = ttl)
}
ggsave("去批次前_PCA.pdf", pca_check(bind_data,    batch_info$batch, "Before"), width = 8, height = 7)
ggsave("去批次后_PCA.pdf", pca_check(combat_data,  batch_info$batch, "After"),  width = 8, height = 7)

## ======================= 03. limma 差异分析 ==================================
design <- model.matrix(~ 0 + factor(Group$group))
rownames(design) <- Group$sample
colnames(design) <- levels(factor(Group$group))
fit <- eBayes(contrasts.fit(lmFit(combat_data, design),
                            makeContrasts(CA - Control, levels = design)))
deg <- topTable(fit, n = Inf)
deg <- deg %>% mutate(group = case_when(
  logFC >= 0.5  & adj.P.Val <= 0.05 ~ "up",
  logFC <= -0.5 & adj.P.Val <= 0.05 ~ "down",
  TRUE ~ "unchanged"))
deg <- na.omit(deg)
sig <- deg[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.5, ]
## 对齐原文(Fig3韦恩用): 严格阈值 |logFC|>1 & FDR<0.05
deg_strict <- deg[abs(deg$logFC) > 1 & deg$adj.P.Val < 0.05, ]
save(deg, sig, deg_strict, file = "sig.RData")
write.csv(sig, "sig.csv")
cat("DEG(>0.5):", nrow(sig), " | DEG(>1, 对齐原文):", nrow(deg_strict), "\n")
## 若 deg_strict 太少(<5), 可退回 deg[abs(logFC)>0.585 & adj.P.Val<0.05, ],
## 但必须在方法学里写明实际阈值。

## ---- 3.1 火山图 ----
deg1 <- deg; deg1$Genes <- rownames(deg1)
p_vol <- ggplot(deg1, aes(logFC, -log10(adj.P.Val), color = group)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("down" = "#1965B0", "unchanged" = "grey", "up" = "#BC3C28")) +
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4, lwd = 0.8) +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.8) +
  labs(x = "log2(Fold change)", y = "-log10(adj.P.Val)") + theme_bw()
ggsave("火山图.pdf", p_vol, width = 7, height = 7)

## ---- 3.2 差异基因热图 ----
hm <- combat_data[rownames(sig), ]
hm <- t(scale(t(hm))); hm[hm > 1] <- 1; hm[hm < -1] <- -1
ann_col <- data.frame(group = Group$group); rownames(ann_col) <- Group$sample
hm <- hm[, rownames(ann_col)]
pheatmap(hm, cluster_rows = TRUE, cluster_cols = FALSE, border_color = NA,
         show_colnames = FALSE, show_rownames = FALSE,
         annotation_col = ann_col,
         annotation_colors = list(group = c(CA = "#FC5256", Control = "#4682B4")),
         filename = "差异基因热图.pdf", width = 10, height = 8)

## ---- 3.3 GSEA (KEGG legacy, 可选; 需下载 c2.cp.kegg_legacy gmt) ----
genelist <- bitr(rownames(deg), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
genelist <- merge(genelist, deg, by.x = 1, by.y = 0)
gsea_input <- genelist$logFC; names(gsea_input) <- genelist$ENTREZID
gsea_input <- sort(gsea_input, decreasing = TRUE)
gmt <- read.gmt("c2.cp.kegg_legacy.v2026.1.Hs.symbols.gmt")
gmt <- gmt %>% inner_join(bitr(gmt$gene, "SYMBOL", "ENTREZID", "org.Hs.eg.db"),
                          by = c("gene" = "SYMBOL")) %>% dplyr::select(term, ENTREZID)
kk <- GSEA(gsea_input, TERM2GENE = gmt, pvalueCutoff = 1, eps = 0)
kkTab <- as.data.frame(kk); kkTab <- kkTab[kkTab$p.adjust < 0.05, ]
write.csv(kkTab, "GSEA_KEGG.csv")
for (tag in c("up","down")){
  sub <- if (tag == "up") kkTab[kkTab$NES > 0, ] else kkTab[kkTab$NES < 0, ]
  if (nrow(sub) >= 5){
    pdf(paste0("GSEA_", tag, ".pdf"), width = 6.5, height = 5)
    print(gseaplot2(kk, rownames(sub)[1:5], base_size = 8,
                    title = paste("Enriched in", tag)))
    dev.off()
  }
}

## ======================= 04. WGCNA + DEG∩module ==============================
library(WGCNA)
dat <- combat_data[, Group$sample]
gv <- apply(dat, 1, var)
dat_w <- t(dat[gv >= quantile(gv, 0.5), ])
grp <- Group$group[match(rownames(dat_w), Group$sample)]
trait <- data.frame(CA = ifelse(grp == "CA", 1, 0)); rownames(trait) <- rownames(dat_w)

pdf("WGCNA_sampleTree.pdf", width = 10, height = 6)
plot(hclust(dist(dat_w)), main = "Sample clustering to detect outliers"); dev.off()

powers <- c(1:20)
sft <- pickSoftThreshold(dat_w, powerVector = powers, networkType = "unsigned", verbose = 5)
beta <- sft$powerEstimate
if (is.na(beta)) beta <- min(sft$fitIndices[,1][sft$fitIndices[,2] >= 0.85])
cat("soft power =", beta, "\n")
pdf("WGCNA_softThreshold.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold", ylab = "Scale Free R^2", type = "n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers)
abline(h = 0.85, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold", ylab = "Mean Connectivity", type = "n")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers)
dev.off()

net <- blockwiseModules(dat_w, power = beta, networkType = "unsigned",
                        minModuleSize = 100, deepSplit = 2, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 3)
MEs <- net$MEs
mCor <- cor(MEs, trait, use = "p"); mP <- corPvalueStudent(mCor, nrow(dat_w))
pdf("WGCNA_moduleTrait.pdf", width = 5, height = 8)
tm <- paste0(signif(mCor,2), "\n(", signif(mP,1), ")"); dim(tm) <- dim(mCor)
labeledHeatmap(mCor, names(trait), names(MEs), names(MEs), colorLabels = FALSE,
               colors = blueWhiteRed(50), textMatrix = tm, setStdMargins = FALSE,
               cex.text = 0.7, zlim = c(-1,1), main = "Module-trait relationships")
dev.off()
mod_df <- data.frame(module = names(MEs), cor = mCor[,1], p = mP[,1]) %>%
  arrange(desc(abs(cor)))
print(mod_df)
best_mod <- as.numeric(gsub("^ME", "", mod_df$module[1]))
module_genes <- colnames(dat_w)[net$colors == best_mod]
## module membership vs gene significance (对标原文 Fig3E)
gs <- as.numeric(cor(dat_w[, module_genes], trait$CA, use = "p"))
mm <- as.numeric(cor(dat_w[, module_genes], MEs[, paste0("ME", best_mod)], use = "p"))
pdf("WGCNA_MM_GS.pdf", width = 6, height = 6)
plot(mm, gs, pch = 16, col = "turquoise", cex = 0.5,
     xlab = "Module Membership", ylab = "Gene significance for CA",
     main = paste0("cor = ", round(cor(mm, gs), 2))); dev.off()

intersect_genes <- intersect(rownames(deg_strict), module_genes)
cat("DEG∩module =", length(intersect_genes), "\n"); print(intersect_genes)
save(deg_strict, net, module_genes, intersect_genes, file = "WGCNA_result.RData")
## 对标原文 Fig3F 韦恩
pdf("Fig3F_venn.pdf", width = 6, height = 6)
venn::venn(list(DEGs = rownames(deg_strict), WGCNA = module_genes), zcol = c("blue","gold"))
dev.off()

## ======================= 05. 机器学习组合筛选 + 外部验证 =====================
## 5.1 训练/验证矩阵
mk <- function(exp_mat, cl, genes){
  g <- intersect(genes, rownames(exp_mat))
  list(X = as.data.frame(t(exp_mat[g, cl$sample])),
       y = ifelse(cl$group == "CA", 1, 0), genes = g)
}
tr <- mk(combat_data, Group, intersect_genes)
va <- mk(GSE19187_exp, GSE19187_clin, intersect_genes)
trainX <- tr$X; trainy <- tr$y; valX <- va$X; valy <- va$y
save(tr, va, file = "ML_input.RData")

## 5.2 算法筛选器 (输入data.frame X与0/1 y, 返回基因名)
sel_glmnet <- function(X, y, alpha){
  m <- cv.glmnet(as.matrix(X), y, family = "binomial", alpha = alpha, nfolds = 5)
  cf <- as.matrix(coef(m, s = "lambda.min"))
  setdiff(rownames(cf)[cf[,1] != 0], "(Intercept)")
}
sel_step <- function(X, y, direction){
  s <- step(glm(y ~ ., data.frame(y = y, X), family = binomial),
            direction = direction, trace = 0)
  setdiff(names(coef(s)), "(Intercept)")
}
sel_glmboost <- function(X, y){
  m <- glmboost(y ~ ., data.frame(y = y, X), family = Binomial(),
                control = boost_control(mstop = 300))
  names(coef(m))[coef(m) != 0]
}
sel_plsr <- function(X, y){
  res <- tryCatch({
    r <- plsRglm(y ~ ., data.frame(y = y, X), nt = 5, pvals.expli = TRUE)
    names(r$ppvals.expli)[r$ppvals.expli < 0.05]
  }, error = function(e) colnames(X))
  if (length(res) < 3) res <- colnames(X)
  res
}
sel_svm <- function(X, y){
  ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 3)
  sizes <- c(2:min(15, ncol(X)))
  r <- rfe(X, factor(ifelse(y == 1, "CA", "Ctrl")), sizes = sizes,
           rfeControl = ctrl, method = "svmRadial")
  predictors(r)
}
sel_xgb <- function(X, y){
  m <- xgboost(xgb.DMatrix(as.matrix(X), label = y), nrounds = 100,
               max_depth = 3, eta = 0.3, objective = "binary:logistic", verbose = 0)
  imp <- xgb.importance(model = m)
  imp$Feature[1:min(10, nrow(imp))]
}
sel_rf <- function(X, y){
  m <- randomForest(factor(ifelse(y == 1, "CA", "Ctrl")) ~ ., data = X, ntree = 500)
  imp <- importance(m)[, "MeanDecreaseGini"]
  names(sort(imp, decreasing = TRUE))[1:min(10, length(imp))]
}
sel_all <- function(X, y) colnames(X)

run_algo <- function(name, X, y){
  if (grepl("^Enet", name)){
    sel_glmnet(X, y, as.numeric(str_match(name, "alpha=([0-9.]+)")[,2]))
  } else switch(name,
    "Lasso" = sel_glmnet(X, y, 1),
    "Ridge" = sel_glmnet(X, y, 0),
    "Stepglm[forward]"  = sel_step(X, y, "forward"),
    "Stepglm[backward]" = sel_step(X, y, "backward"),
    "Stepglm[both]"     = sel_step(X, y, "both"),
    "glmBoost" = sel_glmboost(X, y),
    "plsRglm"  = sel_plsr(X, y),
    "SVM"      = sel_svm(X, y),
    "XGBoost"  = sel_xgb(X, y),
    "RF"       = sel_rf(X, y),
    "LDA"      = sel_all(X, y),
    "NaiveBayes" = sel_all(X, y),
    stop("unknown: ", name))
}

## 5.3 组合清单 (单算法+两两组合; 原文107种取决于其算法清单, 框架一致即可)
algos <- c("Stepglm[forward]","Stepglm[backward]","Stepglm[both]",
           "Enet[alpha=0.1]","Enet[alpha=0.3]","Enet[alpha=0.5]","Enet[alpha=0.7]","Enet[alpha=0.9]",
           "glmBoost","plsRglm","SVM","LDA","XGBoost","NaiveBayes","RF","Lasso","Ridge")
combo_list <- c(as.list(algos), combn(algos, 2, simplify = FALSE))

eval_genes <- function(genes){
  genes <- intersect(genes, colnames(trainX))
  if (length(genes) < 3) return(NULL)   ## 对齐原文"只保留>3个变量的模型"
  fit  <- glm(trainy ~ ., trainX[, genes, drop = FALSE], family = binomial)
  c(AUC_train    = as.numeric(auc(trainy, predict(fit, type = "response"))),
    AUC_GSE19187 = as.numeric(auc(valy, predict(fit, newdata = valX[, genes, drop = FALSE],
                                                type = "response"))),
    n_genes = length(genes))
}
results <- list()
for (i in seq_along(combo_list)){
  cb <- combo_list[[i]]
  if (i %% 10 == 0) cat(i, "/", length(combo_list), "\n")
  g <- tryCatch({
    g1 <- run_algo(cb[1], trainX, trainy)
    if (length(cb) == 2) run_algo(cb[2], trainX[, g1, drop = FALSE], trainy) else g1
  }, error = function(e) NULL)
  m <- eval_genes(g)
  if (!is.null(m)) results[[i]] <- data.frame(model = paste(cb, collapse = "+"), t(m))
}
ml_res <- bind_rows(results) %>%
  mutate(AUC_avg = (AUC_train + AUC_GSE19187) / 2) %>% arrange(desc(AUC_avg))
write.csv(ml_res, "ML_combo_results.csv", row.names = FALSE)

## 5.4 复刻 Fig4A 排名热图
plot_df <- ml_res %>% head(40) %>%
  pivot_longer(cols = starts_with("AUC_") & !ends_with("avg"),
               names_to = "cohort", values_to = "auc")
p_rank <- ggplot(plot_df, aes(cohort, reorder(model, auc))) +
  geom_tile(aes(fill = auc), color = "white") +
  geom_text(aes(label = round(auc, 3)), size = 2.5) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0.8) +
  labs(y = NULL, fill = "AUC") + theme_bw() +
  theme(axis.text.y = element_text(size = 7))
ggsave("Fig4A_ML_ranking.pdf", p_rank, width = 6, height = 12)

## 5.5 最佳模型 -> 诊断基因
best_model <- ml_res$model[1]
best_genes <- local({
  idx <- which(sapply(combo_list, function(x) paste(x, collapse = "+")) == best_model)
  cb <- combo_list[[idx]]
  g1 <- run_algo(cb[1], trainX, trainy)
  if (length(cb) == 2) run_algo(cb[2], trainX[, g1, drop = FALSE], trainy) else g1
})
best_genes <- intersect(best_genes, colnames(trainX))
print(best_model); print(best_genes)
save(best_model, best_genes, ml_res, file = "ML_result.RData")

## 5.6 复刻 Fig4B/4C: 单基因箱线图 + 单基因ROC
expr_all <- rbind(
  data.frame(trainX[, best_genes], group = factor(ifelse(trainy == 1, "CA", "Control")), cohort = "Training"),
  data.frame(valX[, best_genes],   group = factor(ifelse(valy   == 1, "CA", "Control")), cohort = "GSE19187"))
box_df <- expr_all %>% pivot_longer(cols = all_of(best_genes), names_to = "gene", values_to = "expr")
p_box <- ggplot(box_df, aes(group, expr, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  facet_grid(cohort ~ gene, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4682B4", CA = "#FC5256")) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Fig4B_diagGenes_boxplot.pdf", p_box, width = 2.2 * length(best_genes), height = 6)

## ggroc版本, 每基因单独一个PDF
for (g in best_genes){
  r_tr <- roc(trainy, trainX[[g]], quiet = TRUE)
  r_va <- roc(valy,   valX[[g]],   quiet = TRUE)
  p_roc <- ggroc(list(Training = r_tr, GSE19187 = r_va), size = 1) +
    scale_color_manual(values = c(Training = "#FC5256", GSE19187 = "#4682B4")) +
    ggtitle(paste0(g, "\nTrain AUC = ", round(auc(r_tr), 3),
                   " | GSE19187 AUC = ", round(auc(r_va), 3))) +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 11))
  ggsave(paste0("Fig4C_ROC_", g, ".pdf"), p_roc, width = 6, height = 6)
}

## 组合模型整体ROC (训练+验证)
fit_best <- glm(trainy ~ ., trainX[, best_genes, drop = FALSE], family = binomial)
p_tr <- predict(fit_best, type = "response")
p_va <- predict(fit_best, newdata = valX[, best_genes, drop = FALSE], type = "response")
r_tr_m <- roc(trainy, p_tr, quiet = TRUE)
r_va_m <- roc(valy, p_va, quiet = TRUE)
p_model <- ggroc(list(Training = r_tr_m, GSE19187 = r_va_m), size = 1) +
  scale_color_manual(values = c(Training = "#FC5256", GSE19187 = "#4682B4")) +
  ggtitle(paste0("Diagnostic model (", length(best_genes), " genes)\nTrain AUC = ",
                 round(auc(r_tr_m), 3), " | GSE19187 AUC = ", round(auc(r_va_m), 3))) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
ggsave("Fig4_model_ROC.pdf", p_model, width = 6, height = 6)

## ======================= 06. 一致性聚类分型 + 免疫浸润 =====================
library(ConsensusClusterPlus)
ca_samples <- Group$sample[Group$group == "CA"]
sub_input <- combat_data[intersect_genes, ca_samples]
set.seed(123)
ccp <- ConsensusClusterPlus(as.matrix(sub_input), maxK = 9, reps = 100,
                            pItem = 0.8, pFeature = 1, clusterAlg = "pam",
                            distance = "euclidean", title = "ConsensusCluster",
                            plot = "png")
## 看 ConsensusCluster/consensus.pdf 和 CDF 曲线定 k; 原文 k=2
k <- 2
subtypes <- ccp[[k]]$consensusClass
subtype_df <- data.frame(sample = ca_samples, subtype = paste0("C", subtypes))
table(subtype_df$subtype)

## 6.1 3D PCA (对标原文 Fig5C)
pca <- prcomp(t(sub_input))
pdf("Fig5C_PCA_3D.pdf", width = 7, height = 7)
scatterplot3d(pca$x[,1], pca$x[,2], pca$x[,3],
              color = ifelse(subtypes == 1, "#8E44AD", "#D4AC0D"), pch = 16,
              xlab = "PC1", ylab = "PC2", zlab = "PC3")
legend("topright", c("C1","C2"), col = c("#8E44AD","#D4AC0D"), pch = 16)
dev.off()
save(subtype_df, file = "subtype.RData")

## 6.2 CIBERSORT (需官网下载 CIBERSORT.R + LM22.txt)
source("CIBERSORT.R")
ciber <- CIBERSORT("LM22.txt", as.matrix(combat_data[, ca_samples]), perm = 100)
ciber_frac <- as.data.frame(ciber[, 1:22])
ciber_frac$subtype <- subtype_df$subtype[match(rownames(ciber_frac), subtype_df$sample)]
library(ggpubr)
ciber_df <- ciber_frac %>% pivot_longer(cols = -subtype, names_to = "cell", values_to = "frac")
p_ciber <- ggplot(ciber_df, aes(cell, frac, fill = subtype)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.9)) +
  stat_compare_means(aes(group = subtype), label = "p.format", size = 2.5) +
  scale_fill_manual(values = c(C1 = "#1965B0", C2 = "#C52B20")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("Fig5D_CIBERSORT.pdf", p_ciber, width = 14, height = 6)

## 6.3 分型间DEG (为后续单细胞 scDEG 交集做准备)
design2 <- model.matrix(~ 0 + subtype_df$subtype)
colnames(design2) <- c("C1","C2")
fit2 <- eBayes(contrasts.fit(lmFit(combat_data[, ca_samples], design2),
                             makeContrasts(C1 - C2, levels = design2)))
subtype_deg <- topTable(fit2, n = Inf)
write.csv(subtype_deg, "subtype_deg_C1_vs_C2.csv")
save(ciber_frac, subtype_deg, file = "subtype_analysis.RData")

cat("\n========== bulk 流程完成 ==========\n")
cat("产出文件: sig.RData / WGCNA_result.RData / ML_combo_results.csv / ML_result.RData / subtype.RData\n")
