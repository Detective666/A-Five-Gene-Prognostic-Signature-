# ============================================================
# TNBC five-gene prognostic signature: complete analysis code
# ============================================================

# 0. Packages and paths ----
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(DESeq2)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(survival)
  library(survminer)
  library(glmnet)
  library(timeROC)
  library(IOBR)
  library(GSVA)
})

set.seed(123)
options(timeout = 600)

dir.create("data/TCGA", recursive = TRUE, showWarnings = FALSE)
dir.create("data/GSE58812", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

map_to_symbol <- function(mat, ensembl_ids = rownames(mat), aggregate_fun = c("sum", "mean")) {
  aggregate_fun <- match.arg(aggregate_fun)
  ids <- sub("\\..*$", "", ensembl_ids)
  symbols <- mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  keep <- !is.na(symbols) & symbols != ""
  x <- mat[keep, , drop = FALSE]
  rownames(x) <- symbols[keep]
  if (aggregate_fun == "sum") return(rowsum(x, rownames(x), reorder = FALSE))
  as.matrix(aggregate(as.data.frame(x), by = list(Gene = rownames(x)), FUN = mean) |> tibble::column_to_rownames("Gene"))
}

fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# ============================================================
# 1. TCGA-BRCA download and TNBC cohort ----
# ============================================================
query_clin <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "bcr xml"
)
GDCdownload(query_clin, method = "api", files.per.chunk = 30)
clinical_df <- GDCprepare_clinic(query_clin, clinical.info = "patient")

tnbc_barcodes <- clinical_df %>%
  filter(
    breast_carcinoma_estrogen_receptor_status == "Negative",
    breast_carcinoma_progesterone_receptor_status == "Negative",
    lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
  ) %>%
  distinct(bcr_patient_barcode) %>%
  pull(bcr_patient_barcode)

query_normal <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)
normal_barcodes <- unique(getResults(query_normal)$cases)

query_final <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = unique(c(tnbc_barcodes, normal_barcodes))
)
query_final$results[[1]] <- getResults(query_final) %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal"))

GDCdownload(query_final, files.per.chunk = 10, directory = "data/TCGA")
data_se <- GDCprepare(query_final, directory = "data/TCGA")
saveRDS(data_se, "data/TCGA/data_se.rds")

count_mat <- map_to_symbol(assay(data_se, "unstranded"), aggregate_fun = "sum")
expr_tpm <- map_to_symbol(as.matrix(assay(data_se, "tpm_unstrand")), aggregate_fun = "sum")
expr_log2 <- log2(expr_tpm + 1)

clin_df <- as.data.frame(colData(data_se))
clin_df$sample_id <- rownames(clin_df)
clin_df$patient_id <- substr(clin_df$sample_id, 1, 12)
clin_df$group <- ifelse(clin_df$sample_type == "Primary Tumor", "Tumor",
                        ifelse(clin_df$sample_type == "Solid Tissue Normal", "Normal", NA))
clin_df <- clin_df %>% filter(group %in% c("Tumor", "Normal"))
expr_tpm <- expr_tpm[, clin_df$sample_id, drop = FALSE]
expr_log2 <- expr_log2[, clin_df$sample_id, drop = FALSE]
count_mat <- count_mat[, clin_df$sample_id, drop = FALSE]

saveRDS(count_mat, "data/TCGA/count_mat_symbol.rds")
saveRDS(expr_tpm, "data/TCGA/tpm_symbol.rds")
saveRDS(clin_df, "data/TCGA/clinical_sample_level.rds")

# ============================================================
# 2. Fig. 1: differential expression ----
# ============================================================
coldata <- clin_df
rownames(coldata) <- coldata$sample_id
coldata$group <- factor(coldata$group, levels = c("Normal", "Tumor"))

dds <- DESeqDataSetFromMatrix(countData = round(count_mat), colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds <- DESeq(dds)
res_raw <- results(dds, contrast = c("group", "Tumor", "Normal"))
res_shrink <- lfcShrink(dds, coef = "group_Tumor_vs_Normal", type = "normal")
res <- as.data.frame(res_shrink)
res$Gene <- rownames(res)

sig_genes <- res %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
write.csv(sig_genes, "results/tables/DEGs_TNBC_vs_Normal.csv", row.names = FALSE)

volcano_df <- res %>%
  mutate(
    change = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    logP = -log10(padj)
  )

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, logP, color = change)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = c(Down = "#3B6FB6", NS = "grey75", Up = "#C73E3A")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic(base_size = 12) +
  labs(x = "log2 Fold Change", y = "-log10(adjusted P)", title = "TNBC vs Normal DEGs") +
  theme(legend.title = element_blank())
ggsave("results/figures/Fig1A_volcano.pdf", p_volcano, width = 7, height = 5.5)

vsd <- vst(dds, blind = FALSE)
mat_vst <- assay(vsd)
top_genes <- sig_genes %>% arrange(padj, desc(abs(log2FoldChange))) %>% slice_head(n = 50) %>% pull(Gene)
top_genes <- intersect(top_genes, rownames(mat_vst))
ann_col <- data.frame(Group = coldata$group, row.names = rownames(coldata))
ord <- order(ann_col$Group)
pdf("results/figures/Fig1B_DEG_heatmap.pdf", width = 8.5, height = 7)
pheatmap(mat_vst[top_genes, ord, drop = FALSE], scale = "row",
         annotation_col = ann_col[ord, , drop = FALSE], show_colnames = FALSE,
         cluster_cols = FALSE, fontsize_row = 6)
dev.off()

# ============================================================
# 3. Fig. 2: GO, KEGG and Hallmark GSEA ----
# ============================================================
gene_list <- sig_genes$Gene

ego <- enrichGO(
  gene = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
  pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE
)
entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
ekegg <- enrichKEGG(gene = entrez, organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.05)

p_go <- dotplot(ego, showCategory = 20, font.size = 8, label_format = 45) + ggtitle("GO biological process")
p_kegg <- dotplot(ekegg, showCategory = 20, font.size = 8, label_format = 45) + ggtitle("KEGG pathways")
ggsave("results/figures/Fig2A_GO.pdf", p_go, width = 7.5, height = 6.5)
ggsave("results/figures/Fig2B_KEGG.pdf", p_kegg, width = 7.5, height = 6.5)

res_gsea <- as.data.frame(res_raw)
res_gsea$Gene <- rownames(res_gsea)
res_gsea <- res_gsea %>% filter(!is.na(stat), Gene != "")
gene_rank <- res_gsea$stat
names(gene_rank) <- res_gsea$Gene
gene_rank <- sort(gene_rank, decreasing = TRUE)

hallmark_gmt <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>% distinct()

gsea_hallmark <- GSEA(
  geneList = gene_rank, TERM2GENE = hallmark_gmt, minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE
)
gsea_sig <- as.data.frame(gsea_hallmark) %>% filter(p.adjust < 0.05)
write.csv(gsea_sig, "results/tables/Hallmark_GSEA.csv", row.names = FALSE)

n_show <- min(20, nrow(gsea_sig))
plot_gsea <- gsea_sig %>%
  slice_max(order_by = abs(NES), n = n_show) %>%
  mutate(
    Direction = ifelse(NES > 0, "Enriched in TNBC", "Enriched in Normal"),
    Pathway = str_replace_all(str_remove(ID, "^HALLMARK_"), "_", " ")
  ) %>%
  arrange(NES) %>% mutate(Pathway = factor(Pathway, levels = Pathway))

p_gsea <- ggplot(plot_gsea, aes(NES, Pathway, fill = Direction)) +
  geom_col(width = 0.75) + geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Enriched in Normal" = "#E67C73", "Enriched in TNBC" = "#27A7A4")) +
  theme_classic(base_size = 11) + labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme(legend.title = element_blank(), legend.position = "top")
ggsave("results/figures/Fig2C_Hallmark_GSEA_NES.pdf", p_gsea, width = 7.5, height = 6)

# ============================================================
# 4. Fig. 3: xCell2 tumor-normal immune landscape ----
# ============================================================
# This section is intentionally written for xCell2 because the revised manuscript states xCell2.
# The final submitted figure and manuscript must be generated with the same algorithm.
if (!requireNamespace("xCell2", quietly = TRUE)) {
  stop("Package 'xCell2' is required for Fig. 3. Install it before running this section.")
}
ref_file <- "data/xCell2/TMECompendium.xCell2Ref.rds"
if (exists("TMECompendium.xCell2Ref")) {
  xcell2_ref <- TMECompendium.xCell2Ref
} else if (file.exists(ref_file)) {
  xcell2_ref <- readRDS(ref_file)
} else {
  stop("Provide the TMECompendium.xCell2Ref object or save it as data/xCell2/TMECompendium.xCell2Ref.rds.")
}
xcell2_scores <- xCell2::xCell2Analysis(
  mix = expr_log2,
  xcell2object = xcell2_ref
)

xcell2_scores <- as.matrix(xcell2_scores)
if (ncol(xcell2_scores) == nrow(clin_df) && all(colnames(xcell2_scores) %in% clin_df$sample_id)) {
  xcell2_long <- as.data.frame(t(xcell2_scores)) %>% tibble::rownames_to_column("sample_id")
} else if (nrow(xcell2_scores) == nrow(clin_df) && all(rownames(xcell2_scores) %in% clin_df$sample_id)) {
  xcell2_long <- as.data.frame(xcell2_scores) %>% tibble::rownames_to_column("sample_id")
} else {
  stop("Check xCell2 score orientation and sample identifiers.")
}

xcell2_long <- xcell2_long %>% left_join(clin_df %>% select(sample_id, group), by = "sample_id")
cell_cols <- setdiff(names(xcell2_long), c("sample_id", "group"))

xcell_stats <- lapply(cell_cols, function(cc) {
  x <- xcell2_long[[cc]]
  g <- xcell2_long$group
  wt <- wilcox.test(x ~ g)
  data.frame(CellType = cc, Pvalue = wt$p.value,
             Tumor_median = median(x[g == "Tumor"], na.rm = TRUE),
             Normal_median = median(x[g == "Normal"], na.rm = TRUE))
}) %>% bind_rows() %>% mutate(FDR = p.adjust(Pvalue, "BH"))
write.csv(xcell_stats, "results/tables/Fig3_xCell2_TNBC_vs_Normal.csv", row.names = FALSE)

xcell_mat <- t(as.matrix(xcell2_long[, cell_cols, drop = FALSE]))
colnames(xcell_mat) <- xcell2_long$sample_id
xcell_z <- t(scale(t(xcell_mat)))
xcell_z[is.na(xcell_z)] <- 0
ann3 <- data.frame(Group = xcell2_long$group, row.names = xcell2_long$sample_id)
ord3 <- order(ann3$Group)
pdf("results/figures/Fig3A_xCell2_heatmap.pdf", width = 10, height = 9)
pheatmap(xcell_z[, ord3, drop = FALSE], annotation_col = ann3[ord3, , drop = FALSE],
         cluster_cols = FALSE, show_colnames = FALSE, fontsize_row = 6)
dev.off()

plot_xcell <- xcell2_long %>%
  pivot_longer(all_of(cell_cols), names_to = "CellType", values_to = "Score")
p_xcell <- ggplot(plot_xcell, aes(group, Score, fill = group)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ CellType, scales = "free_y") +
  theme_classic(base_size = 9) + labs(x = NULL, y = "xCell2 enrichment score") +
  theme(legend.position = "none")
ggsave("results/figures/Fig3B-D_xCell2_boxplots.pdf", p_xcell, width = 10, height = 9)

# ============================================================
# 5. Survival cohort and five-gene model ----
# ============================================================
get_numeric <- function(df, candidates) {
  out <- rep(NA_real_, nrow(df))
  for (cc in intersect(candidates, names(df))) {
    out <- dplyr::coalesce(out, suppressWarnings(as.numeric(as.character(df[[cc]]))))
  }
  out
}

surv_df <- clin_df
surv_df$OS_event <- ifelse(surv_df$vital_status == "Dead", 1,
                            ifelse(surv_df$vital_status == "Alive", 0, NA))
dtd <- get_numeric(surv_df, c("days_to_death", "paper_days_to_death"))
dlfu <- get_numeric(surv_df, c("days_to_last_follow_up", "paper_days_to_last_followup", "paper_days_to_last_follow_up"))
surv_df$OS_time_days <- ifelse(surv_df$OS_event == 1, dtd, dlfu)
surv_df$OS_time_days <- ifelse(is.na(surv_df$OS_time_days), dplyr::coalesce(dtd, dlfu), surv_df$OS_time_days)
surv_df$OS_time_months <- surv_df$OS_time_days / 30.44

clin_tumor <- surv_df %>%
  filter(group == "Tumor", !is.na(OS_event), !is.na(OS_time_months), OS_time_months > 0) %>%
  arrange(patient_id, sample_id) %>% distinct(patient_id, .keep_all = TRUE) %>%
  filter(sample_id %in% colnames(expr_log2))

train_df <- clin_tumor
expr_surv_tumor <- expr_log2[, train_df$sample_id, drop = FALSE]
write.csv(train_df, "results/tables/TCGA_survival_cohort.csv", row.names = FALSE)

candidate_genes <- intersect(sig_genes$Gene, rownames(expr_surv_tumor))
run_univ_gene <- function(gene) {
  d <- data.frame(time = train_df$OS_time_months, event = train_df$OS_event,
                  expr = as.numeric(expr_surv_tumor[gene, train_df$sample_id])) %>% na.omit()
  if (nrow(d) < 10 || sd(d$expr) == 0) return(NULL)
  fit <- tryCatch(coxph(Surv(time, event) ~ expr, data = d), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)
  data.frame(Gene = gene, HR = s$conf.int[1,"exp(coef)"],
             Lower95 = s$conf.int[1,"lower .95"], Upper95 = s$conf.int[1,"upper .95"],
             Pvalue = s$coefficients[1,"Pr(>|z|)"])
}

univ_cox_df <- lapply(candidate_genes, run_univ_gene) %>% bind_rows() %>%
  mutate(FDR = p.adjust(Pvalue, "BH")) %>% arrange(Pvalue)
write.csv(univ_cox_df, "results/tables/univariate_Cox_DEGs.csv", row.names = FALSE)

lasso_genes <- univ_cox_df %>% filter(Pvalue < 0.01) %>% pull(Gene) %>% unique()
x <- t(expr_surv_tumor[lasso_genes, train_df$sample_id, drop = FALSE])
y <- Surv(train_df$OS_time_months, train_df$OS_event)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, standardize = TRUE)

coef_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))
coef_min <- as.matrix(coef(cvfit, s = "lambda.min"))
if (sum(coef_1se != 0) >= 2) {
  coef_use <- coef_1se
  lambda_used <- "lambda.1se"
} else {
  coef_use <- coef_min
  lambda_used <- "lambda.min"
}
coef_df <- data.frame(Gene = rownames(coef_use), Coefficient = as.numeric(coef_use)) %>%
  filter(Coefficient != 0) %>%
  left_join(univ_cox_df %>% select(Gene, Univ_P = Pvalue), by = "Gene") %>%
  arrange(Univ_P)

coef_final <- coef_df %>% slice_head(n = 5) %>% select(Gene, Coefficient, Univ_P)
write.csv(coef_df, "results/tables/LASSO_nonzero_genes.csv", row.names = FALSE)
write.csv(coef_final, "results/tables/final_five_gene_coefficients.csv", row.names = FALSE)
stopifnot(identical(coef_final$Gene, c("SPRR4", "P2RX6", "ARHGAP6", "IVL", "HES2")))

pdf("results/figures/Fig4A_LASSO_CV.pdf", width = 7, height = 6)
plot(cvfit)
dev.off()

# Standardize the five genes using TCGA development-cohort means and SDs.
train_ref <- list(
  mean = rowMeans(expr_surv_tumor[coef_final$Gene, train_df$sample_id, drop = FALSE], na.rm = TRUE),
  sd = apply(expr_surv_tumor[coef_final$Gene, train_df$sample_id, drop = FALSE], 1, sd, na.rm = TRUE)
)
train_ref$sd[is.na(train_ref$sd) | train_ref$sd == 0] <- 1

calc_score <- function(expr_mat, coef_table, ref) {
  genes <- coef_table$Gene
  stopifnot(all(genes %in% rownames(expr_mat)))
  x <- t(expr_mat[genes, , drop = FALSE])
  x <- sweep(x, 2, ref$mean[genes], "-")
  x <- sweep(x, 2, ref$sd[genes], "/")
  as.numeric(x %*% coef_table$Coefficient)
}
train_df$risk_score <- calc_score(expr_surv_tumor[, train_df$sample_id, drop = FALSE], coef_final, train_ref)
cutoff <- median(train_df$risk_score, na.rm = TRUE)
train_df$risk_group <- factor(ifelse(train_df$risk_score >= cutoff, "High", "Low"), levels = c("Low", "High"))
rownames(train_df) <- train_df$sample_id
write.csv(train_df, "results/tables/TCGA_risk_scores.csv", row.names = FALSE)

fit_km <- survfit(Surv(OS_time_months, OS_event) ~ risk_group, data = train_df)
p_km <- ggsurvplot(fit_km, data = train_df, pval = TRUE, risk.table = TRUE,
                   xlab = "Time (months)", ylab = "Overall survival probability",
                   legend.title = "Risk group")
pdf("results/figures/Fig4B_KM_TCGA.pdf", width = 7, height = 7)
print(p_km)
dev.off()

# ============================================================
# 6. GSE58812 external validation and Fig. 4C-E ----
# ============================================================
phe_GSE58812 <- readRDS("data/GSE58812/phe_GSE58812.rds")
expr_GSE58812 <- as.matrix(readRDS("data/GSE58812/count_GSE58812.rds"))

if (!all(coef_final$Gene %in% rownames(expr_GSE58812)) && all(coef_final$Gene %in% colnames(expr_GSE58812))) {
  expr_GSE58812 <- t(expr_GSE58812)
}
stopifnot(all(coef_final$Gene %in% rownames(expr_GSE58812)))

test_df <- phe_GSE58812
if (!"sample_id" %in% names(test_df)) {
  id_candidates <- intersect(c("geo_accession", "sample", "ID", "id"), names(test_df))
  if (length(id_candidates) == 0) stop("Add a sample_id column to phe_GSE58812.")
  test_df$sample_id <- as.character(test_df[[id_candidates[1]]])
}
if (!all(c("OS_time_months", "OS_event") %in% names(test_df))) {
  stop("phe_GSE58812 must contain OS_time_months and OS_event, as used in the original analysis.")
}

common_test <- intersect(test_df$sample_id, colnames(expr_GSE58812))
test_df <- test_df %>% filter(sample_id %in% common_test)
expr_test <- expr_GSE58812[, test_df$sample_id, drop = FALSE]
test_df$risk_score <- calc_score(expr_test, coef_final, train_ref)
test_df$risk_group <- factor(ifelse(test_df$risk_score >= cutoff, "High", "Low"), levels = c("Low", "High"))
rownames(test_df) <- test_df$sample_id
write.csv(test_df, "results/tables/GSE58812_risk_scores.csv", row.names = FALSE)

roc_times <- c(12, 36, 60)
make_time_roc <- function(df, cohort) {
  tt <- roc_times[roc_times < max(df$OS_time_months, na.rm = TRUE)]
  obj <- timeROC(T = df$OS_time_months, delta = df$OS_event, marker = df$risk_score,
                 cause = 1, times = tt, iid = TRUE)
  ci <- confint(obj, level = 0.95)$CI_AUC
  if (max(ci, na.rm = TRUE) > 1.5) ci <- ci / 100
  tab <- data.frame(Cohort = cohort, Time_months = tt, AUC = as.numeric(obj$AUC),
                    Lower95 = pmax(0, ci[, "2.5%"]), Upper95 = pmin(1, ci[, "97.5%"]))
  list(object = obj, table = tab)
}

roc_train <- make_time_roc(train_df, "TCGA")
roc_test <- make_time_roc(test_df, "GSE58812")
auc_summary <- bind_rows(roc_train$table, roc_test$table)
write.csv(auc_summary, "results/tables/time_dependent_AUCs.csv", row.names = FALSE)

plot_time_roc <- function(obj, tab, file, title) {
  pdf(file, width = 7, height = 6)
  plot(obj, time = tab$Time_months[1], col = 1, title = FALSE)
  if (nrow(tab) > 1) for (i in 2:nrow(tab)) plot(obj, time = tab$Time_months[i], add = TRUE, col = i)
  legend("bottomright", legend = paste0(tab$Time_months, " months AUC = ", sprintf("%.3f", tab$AUC)),
         col = seq_len(nrow(tab)), lwd = 2, bty = "n")
  title(title)
  dev.off()
}
plot_time_roc(roc_train$object, roc_train$table, "results/figures/Fig4C_ROC_TCGA.pdf", "TCGA")
plot_time_roc(roc_test$object, roc_test$table, "results/figures/Fig4D_ROC_GSE58812.pdf", "GSE58812")

ord <- order(train_df$risk_score)
mat5 <- expr_surv_tumor[coef_final$Gene, train_df$sample_id[ord], drop = FALSE]
ann5 <- data.frame(Risk = train_df$risk_group[ord], row.names = train_df$sample_id[ord])
pdf("results/figures/Fig4E_five_gene_heatmap.pdf", width = 10, height = 4)
pheatmap(mat5, scale = "row", annotation_col = ann5, cluster_cols = FALSE,
         show_colnames = FALSE, cluster_rows = FALSE)
dev.off()

# Bootstrap optimism-corrected C-index
boot_dat <- train_df %>% select(OS_time_months, OS_event, risk_score) %>% na.omit()
fit0 <- coxph(Surv(OS_time_months, OS_event) ~ risk_score, data = boot_dat)
lp0 <- predict(fit0, type = "lp")
c_app <- concordance(Surv(OS_time_months, OS_event) ~ lp0, data = boot_dat, reverse = TRUE)$concordance

B <- 1000
optimism <- rep(NA_real_, B)
for (b in seq_len(B)) {
  idx <- sample(seq_len(nrow(boot_dat)), replace = TRUE)
  db <- boot_dat[idx, , drop = FALSE]
  fitb <- tryCatch(coxph(Surv(OS_time_months, OS_event) ~ risk_score, data = db), error = function(e) NULL)
  if (is.null(fitb)) next
  lp_b <- predict(fitb, newdata = db, type = "lp")
  lp_o <- predict(fitb, newdata = boot_dat, type = "lp")
  cb <- concordance(Surv(OS_time_months, OS_event) ~ lp_b, data = db, reverse = TRUE)$concordance
  co <- concordance(Surv(OS_time_months, OS_event) ~ lp_o, data = boot_dat, reverse = TRUE)$concordance
  optimism[b] <- cb - co
}
bootstrap_summary <- data.frame(
  Apparent_C_index = c_app,
  Mean_optimism = mean(optimism, na.rm = TRUE),
  Optimism_corrected_C_index = c_app - mean(optimism, na.rm = TRUE),
  Bootstrap_iterations = B
)
write.csv(bootstrap_summary, "results/tables/Supplementary_Table_S1_bootstrap_Cindex.csv", row.names = FALSE)

get_event_summary <- function(df, cohort) {
  lapply(roc_times, function(tt) data.frame(
    Cohort = cohort, Time_months = tt, Total_N = nrow(df),
    Events_by_time = sum(df$OS_event == 1 & df$OS_time_months <= tt, na.rm = TRUE),
    Followed_beyond_time = sum(df$OS_time_months >= tt, na.rm = TRUE),
    Censored_before_time = sum(df$OS_event == 0 & df$OS_time_months < tt, na.rm = TRUE)
  )) %>% bind_rows()
}
event_summary <- bind_rows(get_event_summary(train_df, "TCGA"), get_event_summary(test_df, "GSE58812"))
write.csv(event_summary, "results/tables/Supplementary_Table_S1_event_summary.csv", row.names = FALSE)

# ============================================================
# 7. Clinical characteristics and Cox analyses ----
# ============================================================
clinical_analysis <- train_df

age <- get_numeric(clinical_analysis, c("age_at_diagnosis", "age_at_initial_pathologic_diagnosis"))
clinical_analysis$Age <- ifelse(age > 200, age / 365.25, age)
clinical_analysis$T_stage <- case_when(
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_t)), "^T1") ~ "T1",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_t)), "^T2") ~ "T2",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_t)), "^T3") ~ "T3",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_t)), "^T4") ~ "T4",
  TRUE ~ NA_character_
)
clinical_analysis$N_stage <- case_when(
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_n)), "^N0") ~ "N0",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_n)), "^N1") ~ "N1",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_n)), "^N2") ~ "N2",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_n)), "^N3") ~ "N3",
  TRUE ~ NA_character_
)
clinical_analysis$M_stage <- case_when(
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_m)), "^M0") ~ "M0",
  str_detect(toupper(as.character(clinical_analysis$ajcc_pathologic_m)), "^M1") ~ "M1",
  TRUE ~ NA_character_
)
stage_raw <- toupper(as.character(clinical_analysis$ajcc_pathologic_stage))
clinical_analysis$Stage_group <- case_when(
  str_detect(stage_raw, "STAGE I($|A|B)") ~ "Stage I",
  str_detect(stage_raw, "STAGE II($|A|B|C)") ~ "Stage II",
  str_detect(stage_raw, "STAGE III($|A|B|C)") ~ "Stage III",
  str_detect(stage_raw, "STAGE IV") ~ "Stage IV",
  TRUE ~ NA_character_
)

extract_treatment_status <- function(x, pattern) {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0 || !"treatment_type" %in% names(x)) return("Unknown")
  hit <- str_detect(as.character(x$treatment_type), regex(pattern, ignore_case = TRUE))
  if (!any(hit, na.rm = TRUE)) return("Unknown")
  if ("treatment_or_therapy" %in% names(x)) {
    z <- tolower(as.character(x$treatment_or_therapy[hit]))
    if (any(z == "yes", na.rm = TRUE)) return("Yes")
    if (length(z) > 0 && all(z == "no", na.rm = TRUE)) return("No")
  }
  "Yes"
}

if ("treatments" %in% names(clinical_analysis)) {
  clinical_analysis$Chemotherapy <- vapply(clinical_analysis$treatments, extract_treatment_status,
                                           pattern = "Chemotherapy|Pharmaceutical", FUN.VALUE = character(1))
  clinical_analysis$Radiotherapy <- vapply(clinical_analysis$treatments, extract_treatment_status,
                                           pattern = "Radiation", FUN.VALUE = character(1))
  clinical_analysis$Surgery <- vapply(clinical_analysis$treatments, extract_treatment_status,
                                      pattern = "Surgery", FUN.VALUE = character(1))
} else {
  clinical_analysis$Chemotherapy <- clinical_analysis$Radiotherapy <- clinical_analysis$Surgery <- "Unknown"
}

# Table 1: use the full group denominator for displayed percentages.
fmt_med_iqr <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, .25), quantile(x, .75))
}
fmt_n_pct_full <- function(x, level) sprintf("%d (%.1f%%)", sum(x == level, na.rm = TRUE), 100 * sum(x == level, na.rm = TRUE) / length(x))
fmt_missing <- function(x) sprintf("%d (%.1f%%)", sum(is.na(x)), 100 * sum(is.na(x)) / length(x))
get_cat_p <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  tab <- table(x[keep], group[keep])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  ex <- suppressWarnings(tryCatch(chisq.test(tab)$expected, error = function(e) NULL))
  if (is.null(ex) || any(ex < 5)) fisher.test(tab, simulate.p.value = TRUE, B = 10000)$p.value else chisq.test(tab, correct = FALSE)$p.value
}

low_df <- clinical_analysis %>% filter(risk_group == "Low")
high_df <- clinical_analysis %>% filter(risk_group == "High")

vars_levels <- list(
  "Pathological T category" = c("T1", "T2", "T3", "T4"),
  "Pathological N category" = c("N0", "N1", "N2", "N3"),
  "Pathological M category" = c("M0", "M1"),
  "Pathological stage" = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  "Chemotherapy" = c("Yes", "No", "Unknown"),
  "Radiotherapy" = c("Yes", "No", "Unknown"),
  "Surgery" = c("Yes", "No", "Unknown")
)
var_cols <- c("T_stage", "N_stage", "M_stage", "Stage_group", "Chemotherapy", "Radiotherapy", "Surgery")

tab1 <- data.frame(Characteristic = "Age, years",
                   Total = fmt_med_iqr(clinical_analysis$Age), Low = fmt_med_iqr(low_df$Age), High = fmt_med_iqr(high_df$Age),
                   P_value = fmt_p(wilcox.test(Age ~ risk_group, data = clinical_analysis)$p.value))
for (i in seq_along(vars_levels)) {
  nm <- names(vars_levels)[i]; vv <- var_cols[i]; levs <- vars_levels[[i]]
  p <- if (vv %in% c("Chemotherapy","Radiotherapy","Surgery")) get_cat_p(ifelse(clinical_analysis[[vv]] == "Unknown", NA, clinical_analysis[[vv]]), clinical_analysis$risk_group) else get_cat_p(clinical_analysis[[vv]], clinical_analysis$risk_group)
  tab1 <- bind_rows(tab1, data.frame(Characteristic = nm, Total = "", Low = "", High = "", P_value = fmt_p(p)))
  for (lv in levs) {
    tab1 <- bind_rows(tab1, data.frame(Characteristic = paste0("  ", lv),
      Total = fmt_n_pct_full(clinical_analysis[[vv]], lv), Low = fmt_n_pct_full(low_df[[vv]], lv), High = fmt_n_pct_full(high_df[[vv]], lv), P_value = ""))
  }
  if (!vv %in% c("Chemotherapy","Radiotherapy","Surgery")) {
    tab1 <- bind_rows(tab1, data.frame(Characteristic = "  Missing", Total = fmt_missing(clinical_analysis[[vv]]),
                                      Low = fmt_missing(low_df[[vv]]), High = fmt_missing(high_df[[vv]]), P_value = ""))
  }
}
colnames(tab1)[2:4] <- c(paste0("Total (n = ", nrow(clinical_analysis), ")"),
                          paste0("Low-risk (n = ", nrow(low_df), ")"),
                          paste0("High-risk (n = ", nrow(high_df), ")"))
write.csv(tab1, "results/tables/Table1_clinicopathological_characteristics.csv", row.names = FALSE)

# Expanded Cox analyses
cox_dat <- clinical_analysis %>% mutate(
  T_group = factor(ifelse(T_stage %in% c("T3","T4"), "T3-4", ifelse(T_stage %in% c("T1","T2"), "T1-2", NA)), levels = c("T1-2","T3-4")),
  N_group = factor(ifelse(N_stage == "N0", "N0", ifelse(N_stage %in% c("N1","N2","N3"), "N+", NA)), levels = c("N0","N+")),
  M_group = factor(M_stage, levels = c("M0","M1")),
  Stage_binary = factor(ifelse(Stage_group %in% c("Stage I","Stage II"), "I-II", ifelse(Stage_group %in% c("Stage III","Stage IV"), "III-IV", NA)), levels = c("I-II","III-IV")),
  Radiotherapy2 = factor(ifelse(Radiotherapy %in% c("Yes","No"), Radiotherapy, NA), levels = c("No","Yes")),
  Surgery2 = factor(ifelse(Surgery %in% c("Yes","No"), Surgery, NA), levels = c("No","Yes"))
)

cox_one <- function(formula, label) {
  fit <- coxph(formula, data = cox_dat)
  s <- summary(fit)
  data.frame(Variable = label, N = fit$n, Events = fit$nevent,
             HR = s$conf.int[1,"exp(coef)"], Lower95 = s$conf.int[1,"lower .95"],
             Upper95 = s$conf.int[1,"upper .95"], Pvalue = s$coefficients[1,"Pr(>|z|)"])
}

expanded_univ <- bind_rows(
  cox_one(Surv(OS_time_months, OS_event) ~ risk_score, "Risk score"),
  cox_one(Surv(OS_time_months, OS_event) ~ Age, "Age"),
  cox_one(Surv(OS_time_months, OS_event) ~ T_group, "T3-4 vs T1-2"),
  cox_one(Surv(OS_time_months, OS_event) ~ N_group, "N+ vs N0"),
  cox_one(Surv(OS_time_months, OS_event) ~ M_group, "M1 vs M0"),
  cox_one(Surv(OS_time_months, OS_event) ~ Stage_binary, "Stage III-IV vs I-II"),
  cox_one(Surv(OS_time_months, OS_event) ~ Radiotherapy2, "Radiotherapy: Yes vs No"),
  cox_one(Surv(OS_time_months, OS_event) ~ Surgery2, "Surgery: Yes vs No")
)
write.csv(expanded_univ, "results/tables/Supplementary_Table_S2_expanded_univariable_Cox.csv", row.names = FALSE)

primary_fit <- coxph(Surv(OS_time_months, OS_event) ~ risk_score + Age + Stage_binary, data = cox_dat)
primary_summary <- broom::tidy(primary_fit, exponentiate = TRUE, conf.int = TRUE)
write.csv(primary_summary, "results/tables/Table2_primary_multivariable_Cox.csv", row.names = FALSE)

sensitivity_formulas <- list(
  Stage = Surv(OS_time_months, OS_event) ~ risk_score + Age + Stage_binary,
  T = Surv(OS_time_months, OS_event) ~ risk_score + Age + T_group,
  N = Surv(OS_time_months, OS_event) ~ risk_score + Age + N_group,
  Radiotherapy = Surv(OS_time_months, OS_event) ~ risk_score + Age + Radiotherapy2
)
sens <- lapply(names(sensitivity_formulas), function(nm) {
  fit <- coxph(sensitivity_formulas[[nm]], data = cox_dat)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% mutate(Model = nm, N = fit$n, Events = fit$nevent)
}) %>% bind_rows()
write.csv(sens, "results/tables/Supplementary_Table_S3_multivariable_sensitivity_Cox.csv", row.names = FALSE)

# ============================================================
# 8. Fig. 5: CIBERSORT, ESTIMATE, GSVA and gene-pathway links ----
# ============================================================
expr_iobr <- expr_surv_tumor[, train_df$sample_id, drop = FALSE]
iobr_res <- list(
  cibersort = deconvo_tme(eset = expr_iobr, method = "cibersort"),
  estimate = deconvo_tme(eset = expr_iobr, method = "estimate")
)

cib <- as.data.frame(iobr_res$cibersort)
est <- as.data.frame(iobr_res$estimate)
cib$risk_group <- train_df$risk_group[match(cib$ID, train_df$sample_id)]
est$risk_group <- train_df$risk_group[match(est$ID, train_df$sample_id)]

immune_cols <- intersect(c(
  "B_cells_naive_CIBERSORT", "T_cells_CD8_CIBERSORT", "T_cells_regulatory_(Tregs)_CIBERSORT",
  "NK_cells_activated_CIBERSORT", "Macrophages_M1_CIBERSORT", "Macrophages_M2_CIBERSORT"
), names(cib))

cib_plot <- cib %>% select(ID, risk_group, all_of(immune_cols)) %>%
  pivot_longer(all_of(immune_cols), names_to = "CellType", values_to = "Fraction") %>%
  mutate(CellType = str_remove(CellType, "_CIBERSORT$"))
p5a <- ggplot(cib_plot, aes(risk_group, Fraction, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.4) + facet_wrap(~ CellType, scales = "free_y") +
  theme_classic(base_size = 10) + labs(x = NULL, y = "Estimated fraction") + theme(legend.position = "none")
ggsave("results/figures/Fig5A_CIBERSORT.pdf", p5a, width = 9, height = 6)

estimate_cols <- intersect(c("StromalScore_ESTIMATE", "ImmuneScore_ESTIMATE", "StromalScore", "ImmuneScore"), names(est))
est_plot <- est %>% select(ID, risk_group, all_of(estimate_cols)) %>%
  pivot_longer(all_of(estimate_cols), names_to = "ScoreType", values_to = "Score")
p5b <- ggplot(est_plot, aes(risk_group, Score, fill = risk_group)) +
  geom_boxplot() + facet_wrap(~ ScoreType, scales = "free_y") + theme_classic(base_size = 11) +
  labs(x = NULL, y = "ESTIMATE score") + theme(legend.position = "none")
ggsave("results/figures/Fig5B_ESTIMATE.pdf", p5b, width = 6.5, height = 4.5)

cor_dat <- train_df %>% select(sample_id, risk_score) %>% left_join(cib, by = c("sample_id" = "ID"))
cor_res <- lapply(immune_cols, function(cc) {
  ct <- cor.test(cor_dat$risk_score, cor_dat[[cc]], method = "spearman", exact = FALSE)
  data.frame(Cell = str_remove(cc, "_CIBERSORT$"), rho = as.numeric(ct$estimate), Pvalue = ct$p.value)
}) %>% bind_rows() %>% mutate(FDR = p.adjust(Pvalue, "BH"))
p5c <- ggplot(cor_res, aes(rho, reorder(Cell, rho))) +
  geom_segment(aes(x = 0, xend = rho, yend = reorder(Cell, rho))) + geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") + theme_classic(base_size = 11) +
  labs(x = "Spearman rho", y = NULL)
ggsave("results/figures/Fig5C_risk_immune_correlation.pdf", p5c, width = 7, height = 4.5)

hallmark_list <- split(hallmark_gmt$gene_symbol, hallmark_gmt$gs_name)
ssgsea_par <- ssgseaParam(exprData = as.matrix(expr_surv_tumor), geneSets = hallmark_list)
gsva_hallmark <- gsva(ssgsea_par, verbose = FALSE)
gsva_df <- as.data.frame(t(gsva_hallmark)) %>% tibble::rownames_to_column("sample_id") %>%
  left_join(train_df %>% select(sample_id, risk_score, risk_group), by = "sample_id")

focus_pathways <- intersect(c(
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS", "HALLMARK_APICAL_JUNCTION", "HALLMARK_P53_PATHWAY"
), names(gsva_df))

gsva_plot <- gsva_df %>% select(sample_id, risk_group, all_of(focus_pathways)) %>%
  pivot_longer(all_of(focus_pathways), names_to = "Pathway", values_to = "Score")
p5d <- ggplot(gsva_plot, aes(risk_group, Score, fill = risk_group)) +
  geom_boxplot() + facet_wrap(~ Pathway, scales = "free_y") + theme_classic(base_size = 9) +
  labs(x = NULL, y = "ssGSEA score") + theme(legend.position = "none")
ggsave("results/figures/Fig5D_GSVA.pdf", p5d, width = 9, height = 6)

five_genes <- coef_final$Gene
pathways_cor <- intersect(c(
  "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE"
), names(gsva_df))

gene_expr_df <- as.data.frame(t(expr_surv_tumor[five_genes, train_df$sample_id, drop = FALSE])) %>%
  tibble::rownames_to_column("sample_id")
link_df <- gene_expr_df %>% left_join(gsva_df %>% select(sample_id, all_of(pathways_cor)), by = "sample_id")

gene_pathway_cor <- expand_grid(Gene = five_genes, Pathway = pathways_cor) %>%
  rowwise() %>%
  mutate(
    N = sum(complete.cases(link_df[[Gene]], link_df[[Pathway]])),
    test = list(cor.test(link_df[[Gene]], link_df[[Pathway]], method = "spearman", exact = FALSE)),
    Rho = as.numeric(test$estimate), Pvalue = test$p.value
  ) %>%
  ungroup() %>% select(-test) %>% mutate(FDR = p.adjust(Pvalue, "BH"),
    Sig = case_when(FDR < .001 ~ "***", FDR < .01 ~ "**", FDR < .05 ~ "*", TRUE ~ ""))
write.csv(gene_pathway_cor, "results/tables/Fig5E_gene_pathway_correlations.csv", row.names = FALSE)

pathway_labels <- c(
  HALLMARK_TGF_BETA_SIGNALING = "TGF-beta signaling",
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = "EMT",
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_ANGIOGENESIS = "Angiogenesis",
  HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma response",
  HALLMARK_INFLAMMATORY_RESPONSE = "Inflammatory response"
)
gene_pathway_cor$PathwayLabel <- pathway_labels[gene_pathway_cor$Pathway]
p5e <- ggplot(gene_pathway_cor, aes(PathwayLabel, factor(Gene, levels = five_genes), fill = Rho)) +
  geom_tile(color = "white") + geom_text(aes(label = paste0(sprintf("%.2f", Rho), Sig)), size = 3.3) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C73E3A", midpoint = 0, limits = c(-1,1)) +
  theme_classic(base_size = 10) + labs(x = NULL, y = NULL, fill = "rho") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
ggsave("results/figures/Fig5E_five_gene_pathway_heatmap.pdf", p5e, width = 8.5, height = 4.8)

