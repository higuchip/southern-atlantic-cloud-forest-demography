dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

build_family <- function(coefficients, variance, gradient_label) {
  stopifnot(
    nrow(coefficients) == 6,
    nrow(variance) == 6,
    identical(as.character(coefficients$Modelo), as.character(variance$Modelo))
  )

  responses <- c(
    "Mortality", "Recruitment", "Basal-area loss", "Basal-area gain",
    "Net abundance change", "Net basal-area change"
  )

  data.frame(
    gradient = gradient_label,
    response = responses,
    beta = coefficients$Coef,
    ci_lower = coefficients$IC_inf,
    ci_upper = coefficients$IC_sup,
    p = coefficients$P_val,
    q = p.adjust(coefficients$P_val, method = "BH"),
    r2_fixed_pct = variance$R2_marginal_pct,
    r2_area_increment_pct = variance$R2_aleatorio_pct,
    stringsAsFactors = FALSE
  )
}

fertility <- build_family(
  res_coefs_fert,
  res_var_fert_export,
  "Soil fertility (PC1)"
)

canopy_associated <- build_family(
  res_coefs_cd,
  res_var_cd_export,
  "Canopy-associated gradient (PC2)"
)

table1 <- rbind(fertility, canopy_associated)
write.csv(
  table1,
  "results/tables/Table1_model_results_full_precision.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

formatted <- transform(
  table1,
  beta = sprintf("%.4f", beta),
  ci = sprintf("%.4f to %.4f", ci_lower, ci_upper),
  p = sprintf("%.3f", p),
  q = sprintf("%.3f", q),
  r2_fixed_pct = sprintf("%.1f", r2_fixed_pct),
  r2_area_increment_pct = sprintf("%.1f", r2_area_increment_pct)
)[c(
  "gradient", "response", "beta", "ci", "p", "q",
  "r2_fixed_pct", "r2_area_increment_pct"
)]

write.csv(
  formatted,
  "results/tables/Table1_model_results_formatted.csv",
  row.names = FALSE,
  quote = TRUE,
  fileEncoding = "UTF-8"
)

header <- paste(
  "| Gradient | Response | beta | 95% CI | p | q |",
  "R2 fixed (%) | R2 area increment (%) |"
)
separator <- "|---|---|---:|---|---:|---:|---:|---:|"
rows <- apply(formatted, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
writeLines(
  c(header, separator, rows),
  "results/tables/Table1_model_results.md",
  useBytes = TRUE
)

cat("Table 1 files written to results/tables\n")
