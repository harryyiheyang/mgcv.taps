files <- list.files(file.path("example", "re"), pattern = "_re_pvalues\\.rds$",
                    full.names = TRUE)
rows <- list()
for (f in files) {
  fam <- sub("_re_pvalues\\.rds$", "", basename(f))
  pv <- readRDS(f)
  for (s in names(pv)) {
    p <- pv[[s]]
    rows[[paste(fam, s)]] <- data.frame(
      family = fam, scenario = s, n_sim = length(p),
      reject_0.01 = mean(p < 0.01), reject_0.05 = mean(p < 0.05),
      reject_0.10 = mean(p < 0.10), median_p = median(p),
      ks_p = suppressWarnings(ks.test(p, "punif")$p.value)
    )
  }
}
res <- do.call(rbind, rows)
res <- res[order(res$family, res$scenario), ]
print(res, row.names = FALSE, digits = 3)
write.csv(res, file.path("example", "re", "re_summary.csv"), row.names = FALSE)
