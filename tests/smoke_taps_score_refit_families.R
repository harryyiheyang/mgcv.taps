library(mgcv.taps)

nsim <- 20L
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) nsim <- as.integer(args[1L])
if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
  stop("nsim must be a positive integer.")
}

example_files <- c(
  gaussian = "example/gaussian_example.R",
  binary_probit = "example/binary_example.R",
  binomial_counts = "example/Binomial_example.R",
  negative_binomial = "example/nb_example.R",
  scaled_t = "example/scat_example.R",
  beta_regression = "example/sbeta_example.R",
  tweedie = "example/tweedie_example.R",
  ordered_categorical = "example/ocat_example.R",
  zero_inflated_poisson = "example/ziP_example.R",
  cox_ph = "example/cox_example.R"
)

base_score <- mgcv.taps::taps_score_test
results <- vector("list", length(example_files))
names(results) <- names(example_files)
warn0 <- getOption("warn")

for (family_name in names(example_files)) {
  code <- readLines(example_files[[family_name]], warn = FALSE)
  code <- code[!grepl("^devtools::load_all", code)]
  code <- sub("^nsim <- 100L$", paste0("nsim <- ", nsim, "L"), code)
  code <- gsub(
    "c(null = 100L, alternative = 100L)",
    "c(null = nsim, alternative = nsim)", code, fixed = TRUE
  )

  run <- new.env(parent = globalenv())
  run$base_score <- base_score
  run$family_name <- family_name
  run$nsim_run <- nsim
  run$call_id <- 0L
  run$rows <- vector("list", 2L * nsim)
  run$saveRDS <- function(...) invisible(NULL)
  run$write.csv <- function(...) invisible(NULL)
  run$taps_score_test <- eval(quote(function(fit, ...) {
    call_id <<- call_id + 1L
    time0 <- system.time(ans0 <- base_score(fit, ..., refit = FALSE))[["elapsed"]]
    time1 <- system.time(ans1 <- base_score(fit, ..., refit = TRUE))[["elapsed"]]
    rows[[call_id]] <<- data.frame(
      family = family_name,
      scenario = if (call_id <= nsim_run) "null" else "alternative",
      simulation = (call_id - 1L) %% nsim_run + 1L,
      pvalue_no_refit = ans0$smooth.pvalue,
      pvalue_refit = ans1$smooth.pvalue,
      time_no_refit = time0,
      time_refit = time1
    )
    ans0
  }), envir = run)

  eval(parse(text = code, keep.source = TRUE), envir = run)
  options(warn = warn0)
  if (run$call_id != 2L * nsim) {
    stop(family_name, " did not produce the expected number of score tests.")
  }
  results[[family_name]] <- do.call(rbind, run$rows)
  cat(family_name, "completed\n")
}

pairs <- do.call(rbind, results)
row.names(pairs) <- NULL
pairs$difference <- pairs$pvalue_refit - pairs$pvalue_no_refit
pairs$absolute_difference <- abs(pairs$difference)
if (any(!is.finite(as.matrix(pairs[, 4:9])))) {
  stop("The family smoke test produced non-finite results.")
}

keys <- unique(pairs[c("family", "scenario")])
summary_rows <- vector("list", nrow(keys))
for (i in seq_len(nrow(keys))) {
  ii <- pairs$family == keys$family[i] & pairs$scenario == keys$scenario[i]
  x <- pairs[ii, ]
  summary_rows[[i]] <- data.frame(
    family = keys$family[i], scenario = keys$scenario[i], n = nrow(x),
    mean_no_refit = mean(x$pvalue_no_refit),
    mean_refit = mean(x$pvalue_refit),
    mean_difference = mean(x$difference),
    median_absolute_difference = median(x$absolute_difference),
    max_absolute_difference = max(x$absolute_difference),
    correlation = cor(x$pvalue_no_refit, x$pvalue_refit),
    mean_time_no_refit = mean(x$time_no_refit),
    mean_time_refit = mean(x$time_refit)
  )
}
summary <- do.call(rbind, summary_rows)

write.csv(
  pairs, "tests/smoke_taps_score_refit_family_pairs.csv", row.names = FALSE
)
write.csv(
  summary, "tests/smoke_taps_score_refit_family_summary.csv", row.names = FALSE
)
print(summary)
