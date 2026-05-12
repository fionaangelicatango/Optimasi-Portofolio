# =========================================================
# K-Medoids (PAM) + Penentuan k Optimal
# Data R_2024: tahun 2010-2024 (baris = waktu, kolom = saham)
# =========================================================

suppressPackageStartupMessages({
  library(cluster)
  library(fpc)
  library(KneeArrower)
})

set.seed(1)

# ---------------------------------------------------------
# 1) Standardisasi
# ---------------------------------------------------------
R_2024 <- column_to_rownames(R_2024, var = "Tanggal")
R_2024 <- R_2024 %>% 
  select(-ASRI, -BHIT, -BBKP, -BKSL, -BNBR, -BUMI, -DOID, -ELTY, -INTA, -LPKR, -MIRA, -MLPL, -MPPA, -SMCB, -SSIA, -TPIA, -UNSP)
std <- scale(R_2024)

# ---------------------------------------------------------
# 2) Transpose
# ---------------------------------------------------------
X <- t(std)  # baris = saham, kolom = waktu
rownames(X) <- colnames(R_2024)

# ---------------------------------------------------------
# 3) Jarak antar saham
# ---------------------------------------------------------
d <- dist(X, method = "euclidean")
D <- as.matrix(d)

k_max   <- min(50, nrow(X) - 1)
k_range <- 2:k_max

# =========================================================
# A) Silhouette
# =========================================================
sil_score <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  fit <- pam(d, k = k, diss = TRUE)
  sil_score[i] <- fit$silinfo$avg.width
}

k_opt_sil <- k_range[which.max(sil_score)]
clust_sil <- pam(d, k = k_opt_sil, diss = TRUE)
hasil_cluster_sil <- data.frame(
  Saham   = rownames(X),
  Cluster = factor(clust_sil$clustering)
)

# =========================================================
# B) Elbow
# =========================================================
wss_pam <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  fit <- pam(d, k = k, diss = TRUE)
  wss_pam[i] <- fit$objective[1]
}

knee_pt <- findCutoff(x = k_range, y = wss_pam, method = "curvature")
k_opt_elb <- if (length(knee_pt$x) > 1) min(knee_pt$x) else knee_pt$x
clust_elb <- pam(d, k = k_opt_elb, diss = TRUE)
hasil_cluster_elb <- data.frame(
  Saham   = rownames(X),
  Cluster = factor(clust_elb$clustering)
)

# =========================================================
# C) Gap Statistic
# =========================================================
pam_clust_gap <- function(x, k) {
  dd <- dist(x, method = "euclidean")
  fit <- pam(dd, k = k, diss = TRUE)
  list(cluster = fit$clustering)
}

gap <- clusGap(
  x = as.matrix(X),
  FUNcluster = pam_clust_gap,
  K.max = k_max,
  B = 200
)

k_opt_gap <- maxSE(gap$Tab[,"gap"], gap$Tab[,"SE.sim"], method="Tibs2001SEmax")
clust_gap <- pam(d, k = k_opt_gap, diss = TRUE)
hasil_cluster_gap <- data.frame(
  Saham   = rownames(X),
  Cluster = factor(clust_gap$clustering)
)

# =========================================================
# D) Calinski–Harabasz (CH)
# =========================================================
ch_vals <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  fit <- pam(d, k = k, diss = TRUE)
  st  <- cluster.stats(d = d, clustering = fit$clustering)
  ch_vals[i] <- st$ch
}

k_opt_CH <- k_range[which.max(ch_vals)]
clust_CH <- pam(d, k = k_opt_CH, diss = TRUE)
hasil_cluster_CH <- data.frame(
  Saham   = rownames(X),
  Cluster = factor(clust_CH$clustering)
)

# =========================================================
# E) Davies–Bouldin (DB)
# =========================================================
db_medoids <- function(D, clusters, medoid_ids) {
  cl_ids <- sort(unique(clusters))
  k <- length(cl_ids)
  
  # Si: rata-rata jarak anggota cluster i ke medoid i
  S <- numeric(k)
  for (idx in seq_along(cl_ids)) {
    cl <- cl_ids[idx]
    members <- which(clusters == cl)
    m <- medoid_ids[idx]
    S[idx] <- if (length(members) <= 1) 0 else mean(D[members, m], na.rm = TRUE)
  }
  
  # jarak antar medoid
  M <- matrix(NA_real_, k, k)
  for (i in 1:k) for (j in 1:k) {
    if (i != j) M[i, j] <- D[medoid_ids[i], medoid_ids[j]]
  }
  
  # Rij = (Si+Sj)/Mij ; DB = mean_i max_j Rij
  R <- matrix(NA_real_, k, k)
  for (i in 1:k) for (j in 1:k) {
    if (i != j && !is.na(M[i, j]) && M[i, j] > 0) {
      R[i, j] <- (S[i] + S[j]) / M[i, j]
    }
  }
  
  mean(apply(R, 1, max, na.rm = TRUE), na.rm = TRUE)
}

db_vals <- numeric(length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  fit <- pam(d, k = k, diss = TRUE)
  db_vals[i] <- db_medoids(D, fit$clustering, fit$id.med)
}

k_opt_DB <- k_range[which.min(db_vals)]
clust_DB <- pam(d, k = k_opt_DB, diss = TRUE)
hasil_cluster_DB <- data.frame(
  Saham   = rownames(X),
  Cluster = factor(clust_DB$clustering)
)

# =========================================================
# F) Stability (Bootstrap Jaccard) - fpc::clusterboot
#    (lebih besar lebih baik)
# =========================================================
stab_score <- rep(NA_real_, length(k_range))

for (i in seq_along(k_range)) {
  k <- k_range[i]
  
  cb <- clusterboot(
    as.matrix(X),
    B = 50,
    bootmethod   = "boot",
    clustermethod = pamkCBI,
    k = k
  )
  
  stab_score[i] <- mean(cb$bootmean, na.rm = TRUE)
}

k_opt_stab <- k_range[which.max(stab_score)]
clust_stab <- pam(d, k = k_opt_stab, diss = TRUE)
hasil_cluster_stab <- data.frame(
  Saham   = rownames(X),
  Cluster = factor(clust_stab$clustering)
)

# =========================================================
# Ringkasan kandidat k
# =========================================================
k_candidates <- c(
  Silhouette = k_opt_sil,
  Elbow      = k_opt_elb,
  Gap        = k_opt_gap,
  CH         = k_opt_CH,
  DB         = k_opt_DB,
  Stability  = k_opt_stab
)
k_candidates

# =========================================================
# Hasil clustering final (k=2)
# =========================================================
hasil_cluster_final <- hasil_cluster_stab
medoids_final <- clust_stab$medoids

hasil_cluster_final
medoids_final

library(tidyverse)
saham_stab1 <- hasil_cluster_stab$Saham[hasil_cluster_stab$Cluster == 1]
stab_cluster1 <- R_2024 %>%
  rownames_to_column(var = "Tanggal") %>%
  select(Tanggal, any_of(saham_stab1))

saham_stab2 <- hasil_cluster_stab$Saham[hasil_cluster_stab$Cluster == 2]
stab_cluster2 <- R_2024 %>%
  rownames_to_column(var = "Tanggal") %>%
  select(Tanggal, any_of(saham_stab2))

#write_xlsx(hasil_cluster_final, "Hasil Re-Cluster.xlsx")

# =========================================================
# Karakteristik Klaster
# =========================================================
library(dplyr)
library(ggplot2)

# kalau kolom "Tanggal" ada sebagai kolom, buang dulu:
cluster1 <- stab_cluster1[ , !(names(stab_cluster1) %in% "Tanggal")]

mean_cluster1 <- sapply(cluster1, mean, na.rm = TRUE)
sd_cluster1   <- sapply(cluster1, sd,   na.rm = TRUE)
p5_cluster1   <- sapply(cluster1, quantile, probs = 0.05, na.rm = TRUE, type = 7)

hasil_cluster1 <- data.frame(
  saham = names(mean_cluster1),
  mean  = as.numeric(mean_cluster1),
  sd    = as.numeric(sd_cluster1),
  p5    = as.numeric(p5_cluster1),
  row.names = NULL
)

cluster2 <- stab_cluster2[ , !(names(stab_cluster2) %in% "Tanggal")]

mean_cluster2 <- sapply(cluster2, mean, na.rm = TRUE)
sd_cluster2   <- sapply(cluster2, sd,   na.rm = TRUE)
p5_cluster2   <- sapply(cluster2, quantile, probs = 0.05, na.rm = TRUE, type = 7)

hasil_cluster2 <- data.frame(
  saham = names(mean_cluster2),
  mean  = as.numeric(mean_cluster2),
  sd    = as.numeric(sd_cluster2),
  p5    = as.numeric(p5_cluster2),
  row.names = NULL
)

hasil_cluster <- bind_rows(
  hasil_cluster1 %>% mutate(cluster = "cluster 1"),
  hasil_cluster2 %>% mutate(cluster = "cluster 2")
)

# =========================================================
# PLOT 1: Mean vs Standar Deviasi
# =========================================================
ggplot(hasil_cluster, aes(mean, sd, color = factor(cluster))) +
  geom_point(size = 2, alpha = 0.9) +
  labs(x = "Mean Return", y = "Standar Deviasi", color = "Cluster") +
  theme_minimal() +
  theme(legend.position = "bottom")

# =========================================================
# PLOT 2: Mean vs Percentile 5%
# =========================================================
ggplot(hasil_cluster, aes(x = mean, y = p5, color = factor(cluster))) +
  geom_point(size = 2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Mean Return", y = "Persentil 5%", color = "Cluster") +
  theme_minimal() +
  theme(legend.position = "bottom")