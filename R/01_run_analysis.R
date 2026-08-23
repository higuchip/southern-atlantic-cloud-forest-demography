# =============================================================================
# ANALISE PRINCIPAL
# Entradas: data/tree_census.csv e data/environmental.csv
# Resultados: objetos usados diretamente pelos scripts 02 e 03.
# =============================================================================

set.seed(20260822)

# --- Pacotes ---
library(dplyr)
library(readr)
library(stringr)
library(FactoMineR)
library(lme4)
library(lmerTest)
library(performance)
library(broom.mixed)

cat("\n============================================================\n")
cat("REGENERANDO RESULTADOS NUMERICOS\n")
cat("============================================================\n\n")

# #############################################################################
# PARTE A: PCA AMBIENTAL
# #############################################################################

cat(">>> PARTE A: PCA de Solos <<<\n")

# Ler dados ambientais completos (arquivo correto da Driely)
dados_solos <- read.csv2("data/environmental.csv",
                         header = TRUE,
                         encoding = "UTF-8",
                         dec = ",")
# Corrigir BOM na primeira coluna
names(dados_solos)[1] <- "Area"

stopifnot(
  nrow(dados_solos) == 80,
  n_distinct(dados_solos$Area) == 4,
  !anyDuplicated(dados_solos[c("Area", "Parcela")])
)

area_original <- dados_solos$Area
dados_solos_num <- dados_solos %>% select(-Area, -Parcela, -Id_unico)

# Remover NAs
indices_completos <- complete.cases(dados_solos_num)
dados_solos_num <- dados_solos_num[indices_completos, ]
area_final <- area_original[indices_completos]

cat("Solos: ", nrow(dados_solos_num), " amostras, ", ncol(dados_solos_num), " variaveis\n")

# PCA
dados_scaled <- scale(dados_solos_num)
pca_result <- PCA(
  dados_scaled,
  ncp = min(ncol(dados_scaled), nrow(dados_scaled) - 1),
  graph = FALSE
)
eigenvalues <- data.frame(
  eigenvalue = pca_result$eig[, 1],
  variance.percent = pca_result$eig[, 2],
  cumulative.variance.percent = pca_result$eig[, 3],
  row.names = rownames(pca_result$eig)
)

# #############################################################################
# PARTE B: DINAMICA FLORESTAL
# #############################################################################

cat("\n>>> PARTE B: Dinamica Florestal <<<\n")

# Ler a exportacao analitica sem coordenadas precisas.
dados <- readr::read_csv(
  "data/tree_census.csv",
  locale = readr::locale(encoding = "UTF-8"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  name_repair = "minimal"
) %>% as.data.frame(check.names = FALSE)

stopifnot(
  nrow(dados) == 1609,
  nrow(distinct(dados, AREA_SIGLA, P)) == 80,
  n_distinct(dados$AREA_SIGLA) == 4
)

cat("Dados carregados:", nrow(dados), "linhas\n")

# Funcao para processar multiplos fustes
processar_multiplos_fustes <- function(cap_string) {
  if (is.na(cap_string) || cap_string == "" || cap_string == "0") return(0)
  cap_string <- gsub(",", ".", cap_string)
  fustes <- strsplit(cap_string, "\\+")[[1]]
  fustes_num <- as.numeric(fustes)
  fustes_num <- fustes_num[!is.na(fustes_num)]
  if (length(fustes_num) == 0) return(0)
  dap_fustes <- fustes_num / pi
  dap_equivalente <- sqrt(sum(dap_fustes^2))
  return(dap_equivalente)
}

# Converter CAP em DAP
dados$DAP1 <- sapply(dados$CAP1, processar_multiplos_fustes)
dados$DAP2 <- sapply(dados$CAP2, processar_multiplos_fustes)
dados$DAP3 <- sapply(dados$CAP3, processar_multiplos_fustes)

dados$DAP1[is.na(dados$DAP1)] <- 0
dados$DAP2[is.na(dados$DAP2)] <- 0
dados$DAP3[is.na(dados$DAP3)] <- 0

dados$Parcela_ID <- paste(dados$AREA_SIGLA, dados$P, sep = "_")

# Funcao para calcular taxas por parcela
# `periodo_anos` pode ser escalar (ex.: 4) ou vetor nomeado por sigla de area
# (ex.: c(PNSJ1=4, PNSJ2=4, PNSJ3=4, PNSJ4=3)) para intervalos distintos por area.
calcular_taxas_parcela <- function(data, dap1_col, dap2_col, periodo_anos) {
  df_trabalho <- data.frame(
    Parcela = data$Parcela_ID,
    Especie = data$SPP,
    Area = data$AREA_SIGLA,
    DAP1 = data[[dap1_col]],
    DAP2 = data[[dap2_col]]
  )

  if (length(periodo_anos) == 1) {
    df_trabalho$intervalo <- periodo_anos
  } else {
    df_trabalho$intervalo <- unname(periodo_anos[as.character(df_trabalho$Area)])
    if (any(is.na(df_trabalho$intervalo))) {
      stop("Vetor 'periodo_anos' nao cobre todas as areas presentes nos dados.")
    }
  }

  df_trabalho <- df_trabalho[!(df_trabalho$DAP1 == 0 & df_trabalho$DAP2 == 0), ]

  taxas_parcela <- df_trabalho %>%
    group_by(Parcela) %>%
    summarise(
      intervalo = first(intervalo),
      n_sobreviventes = sum(DAP1 > 0 & DAP2 > 0),
      n_mortas = sum(DAP1 > 0 & DAP2 == 0),
      n_recrutas = sum(DAP1 == 0 & DAP2 > 0),
      n0 = n_sobreviventes + n_mortas,
      n1 = n_sobreviventes + n_recrutas,
      ab0 = sum((pi * DAP1^2) / 40000),
      ab1 = sum((pi * DAP2^2) / 40000),
      ab_mortas = sum((pi * DAP1^2) / 40000 * (DAP1 > 0 & DAP2 == 0)),
      ab_recrutas = sum((pi * DAP2^2) / 40000 * (DAP1 == 0 & DAP2 > 0)),
      ab_ganho_sob = sum(((pi * DAP2^2) / 40000 - (pi * DAP1^2) / 40000) *
                           (DAP1 > 0 & DAP2 > 0) * ((DAP2 - DAP1) > 0)),
      ab_perda_sob = sum(((pi * DAP2^2) / 40000 - (pi * DAP1^2) / 40000) *
                           (DAP1 > 0 & DAP2 > 0) * ((DAP2 - DAP1) < 0)),
      .groups = "drop"
    ) %>%
    mutate(
      taxa_mortalidade = ifelse(n0 > 0, (1 - ((n0 - n_mortas) / n0)^(1/intervalo)) * 100, 0),
      taxa_recrutamento = ifelse(n1 > 0, (1 - (1 - n_recrutas/n1)^(1/intervalo)) * 100, 0),
      taxa_mudanca_liq = ifelse(n0 > 0, ((n1/n0)^(1/intervalo) - 1) * 100, 0),
      taxa_rotatividade = (taxa_mortalidade + taxa_recrutamento) / 2,
      ab_ganho_total = ab_ganho_sob + ab_recrutas,
      ab_perda_total = abs(ab_perda_sob) + ab_mortas,
      taxa_perda_ab = ifelse(ab0 > 0, (1 - ((ab0 - ab_perda_total) / ab0)^(1/intervalo)) * 100, 0),
      taxa_ganho_ab = ifelse(ab1 > 0, (1 - (1 - ab_ganho_total/ab1)^(1/intervalo)) * 100, 0),
      taxa_mudanca_ab = ifelse(ab0 > 0, ((ab1/ab0)^(1/intervalo) - 1) * 100, 0),
      taxa_rotatividade_ab = (taxa_perda_ab + taxa_ganho_ab) / 2
    )
  return(taxas_parcela)
}

# Periodos
# Periodo 1: PNSJ1-3 = 4 anos (2016-2020); PNSJ4 = 3 anos (2017-2020)
taxas_periodo1 <- calcular_taxas_parcela(
  dados, "DAP1", "DAP2",
  c(PNSJ1 = 4, PNSJ2 = 4, PNSJ3 = 4, PNSJ4 = 3)
)
taxas_periodo1$Periodo <- "2016/17-2020"

# Periodo 2: 4 anos para todas as areas
taxas_periodo2 <- calcular_taxas_parcela(dados, "DAP2", "DAP3", 4)
taxas_periodo2$Periodo <- "2020-2024"

# Medias dos dois periodos por parcela
taxas_medias_parcela <- taxas_periodo1 %>%
  select(Parcela,
         mort_p1 = taxa_mortalidade, recr_p1 = taxa_recrutamento,
         mud_liq_p1 = taxa_mudanca_liq, rot_p1 = taxa_rotatividade,
         perda_ab_p1 = taxa_perda_ab, ganho_ab_p1 = taxa_ganho_ab,
         mud_ab_p1 = taxa_mudanca_ab, rot_ab_p1 = taxa_rotatividade_ab) %>%
  full_join(
    taxas_periodo2 %>%
      select(Parcela,
             mort_p2 = taxa_mortalidade, recr_p2 = taxa_recrutamento,
             mud_liq_p2 = taxa_mudanca_liq, rot_p2 = taxa_rotatividade,
             perda_ab_p2 = taxa_perda_ab, ganho_ab_p2 = taxa_ganho_ab,
             mud_ab_p2 = taxa_mudanca_ab, rot_ab_p2 = taxa_rotatividade_ab),
    by = "Parcela"
  ) %>%
  mutate(
    taxa_mortalidade_media = (mort_p1 + mort_p2) / 2,
    taxa_recrutamento_media = (recr_p1 + recr_p2) / 2,
    taxa_mudanca_liq_media = (mud_liq_p1 + mud_liq_p2) / 2,
    taxa_rotatividade_media = (rot_p1 + rot_p2) / 2,
    taxa_perda_ab_media = (perda_ab_p1 + perda_ab_p2) / 2,
    taxa_ganho_ab_media = (ganho_ab_p1 + ganho_ab_p2) / 2,
    taxa_mudanca_ab_media = (mud_ab_p1 + mud_ab_p2) / 2,
    taxa_rotatividade_ab_media = (rot_ab_p1 + rot_ab_p2) / 2,
    Area = str_extract(Parcela, "^[^_]+")
  )

# Tabela final com medias arredondadas
tabela_final <- taxas_medias_parcela %>%
  select(
    Parcela, Area,
    Mortalidade_P1 = mort_p1, Recrutamento_P1 = recr_p1,
    Mudanca_Liq_P1 = mud_liq_p1, Rotatividade_P1 = rot_p1,
    Mortalidade_P2 = mort_p2, Recrutamento_P2 = recr_p2,
    Mudanca_Liq_P2 = mud_liq_p2, Rotatividade_P2 = rot_p2,
    Mortalidade_Media = taxa_mortalidade_media,
    Recrutamento_Media = taxa_recrutamento_media,
    Mudanca_Liq_Media = taxa_mudanca_liq_media,
    Rotatividade_Media = taxa_rotatividade_media,
    Perda_AB_Media = taxa_perda_ab_media,
    Ganho_AB_Media = taxa_ganho_ab_media,
    Mudanca_AB_Media = taxa_mudanca_ab_media,
    Rotatividade_AB_Media = taxa_rotatividade_ab_media
  ) %>%
  arrange(Area, Parcela) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

# #############################################################################
# PARTE C: INDICE DE FERTILIDADE (necessario para modelagem)
# #############################################################################

cat("\n>>> PARTE C: Indice de Fertilidade <<<\n")

# Extrair loadings PC1
loadings_pc1 <- pca_result$var$coord[, 1]

# Determinar direcao do indice (variaveis do arquivo AMBIENTAIS_URUBICI_COMPLETO_DRIELLY)
interpretacao_pc1 <- data.frame(
  Variavel = names(loadings_pc1),
  Loading = loadings_pc1,
  Contribuicao_pct = pca_result$var$contrib[, 1]
) %>%
  mutate(
    Sinal_Fertilidade = case_when(
      Variavel %in% c("Ca", "Mg", "K", "CTC", "SatBase") & Loading > 0 ~ "+",
      Variavel %in% c("Ca", "Mg", "K", "CTC", "SatBase") & Loading < 0 ~ "-",
      Variavel == "Al" & Loading > 0 ~ "-",
      Variavel == "Al" & Loading < 0 ~ "+",
      Variavel == "pH" & Loading > 0 ~ "+",
      Variavel == "pH" & Loading < 0 ~ "-",
      Variavel == "P" & Loading > 0 ~ "+",
      Variavel == "P" & Loading < 0 ~ "-",
      Variavel == "MO" & Loading > 0 ~ "+",
      Variavel == "MO" & Loading < 0 ~ "-",
      TRUE ~ "0"
    )
  )

variaveis_importantes <- interpretacao_pc1 %>%
  filter(Contribuicao_pct > 100/nrow(interpretacao_pc1))

fertilidade_positiva <- sum(variaveis_importantes$Sinal_Fertilidade == "+")
fertilidade_negativa <- sum(variaveis_importantes$Sinal_Fertilidade == "-")
fator_multiplicacao <- ifelse(fertilidade_positiva > fertilidade_negativa, 1, -1)

scores_pc1 <- pca_result$ind$coord[,1]
indice_fertilidade <- scores_pc1 * fator_multiplicacao
indice_fertilidade_0_100 <- scales::rescale(indice_fertilidade, to = c(0, 100))

dados_fertilidade <- data.frame(
  Area = area_final,
  Parcela = dados_solos$Parcela[indices_completos],
  Id_unico = dados_solos$Id_unico[indices_completos],
  PC1_original = round(scores_pc1, 3),
  Indice_Fertilidade = round(indice_fertilidade, 3),
  Indice_Fertilidade_0_100 = round(indice_fertilidade_0_100, 1),
  Classe_Fertilidade = cut(indice_fertilidade_0_100,
                           breaks = c(0, 25, 50, 75, 100),
                           labels = c("Baixa", "Media-Baixa", "Media-Alta", "Alta"),
                           include.lowest = TRUE)
)

cat("  Fertilidade calculada para", nrow(dados_fertilidade), "parcelas\n")

modelo_anova_fert <- aov(Indice_Fertilidade_0_100 ~ factor(Area), data = dados_fertilidade)
anova_summary_fert <- summary(modelo_anova_fert)
anova_pvalue_fert <- anova_summary_fert[[1]][["Pr(>F)"]][1]

tukey_fert <- TukeyHSD(modelo_anova_fert)

anova_fert_export <- data.frame(
  axis = "Soil fertility (PC1)",
  df_between = anova_summary_fert[[1]][["Df"]][1],
  df_within = anova_summary_fert[[1]][["Df"]][2],
  F = anova_summary_fert[[1]][["F value"]][1],
  p = anova_pvalue_fert
)
tukey_fert_matrix <- tukey_fert$`factor(Area)`
tukey_fert_export <- data.frame(
  axis = "Soil fertility (PC1)",
  comparison = rownames(tukey_fert_matrix),
  difference = tukey_fert_matrix[, "diff"],
  lower = tukey_fert_matrix[, "lwr"],
  upper = tukey_fert_matrix[, "upr"],
  p_adjusted = tukey_fert_matrix[, "p adj"],
  row.names = NULL
)


# #############################################################################
# PARTE C2: MODELOS MISTOS PARA O GRADIENTE DE FERTILIDADE
# #############################################################################

cat("\n>>> PARTE C2: Modelos mistos - fertilidade <<<\n")

# Preparar join fertilidade x taxas
criar_id_fert <- function(a, p) {
  clean <- str_replace_all(p, "-", "_")
  paste0(a, "_", str_extract(clean, "\\d+"), str_extract(clean, "[DE]"))
}

part_fert <- dados_fertilidade %>%
  mutate(Parcela_ID = criar_id_fert(Area, Parcela)) %>%
  select(Parcela_ID, Area, Indice_Fertilidade_0_100)

part_taxas_fert <- taxas_medias_parcela %>%
  mutate(Parcela_ID = str_replace(Parcela, "^PNSJ", "")) %>%
  select(Parcela_ID, starts_with("taxa_"))

model_data_fert <- inner_join(part_fert, part_taxas_fert, by = "Parcela_ID") %>% na.omit()
cat("  Dataset integrado:", nrow(model_data_fert), "parcelas\n")

# Ajustar modelos lineares mistos
vars_resp_f <- c("taxa_mortalidade_media", "taxa_recrutamento_media",
                 "taxa_perda_ab_media", "taxa_ganho_ab_media",
                 "taxa_mudanca_liq_media", "taxa_mudanca_ab_media")

nomes_respostas_f <- c("Mortalidade", "Recrutamento", "Perda AB",
                       "Ganho AB", "Mudanca Ind.", "Mudanca AB")

res_coefs_fert <- data.frame()
res_var_fert <- data.frame()

for(i in 1:length(vars_resp_f)) {
  var <- vars_resp_f[i]
  nome <- nomes_respostas_f[i]

  mod <- lmer(as.formula(paste0(var, " ~ Indice_Fertilidade_0_100 + (1|Area)")),
              data = model_data_fert, REML = TRUE)

  tidied <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "Indice_Fertilidade_0_100")

  sig <- ifelse(tidied$p.value < 0.05, "Sim", "Nao")
  r2 <- performance::r2(mod)

  res_coefs_fert <- rbind(res_coefs_fert, data.frame(
    Modelo = nome, Coef = tidied$estimate, Erro = tidied$std.error,
    IC_inf = tidied$conf.low, IC_sup = tidied$conf.high,
    P_val = tidied$p.value, Significativo = sig
  ))

  res_var_fert <- rbind(res_var_fert, data.frame(
    Modelo = nome,
    Fixo = r2$R2_marginal,
    Aleatorio = r2$R2_conditional - r2$R2_marginal,
    Residual = 1 - r2$R2_conditional
  ))
}

# Decomposicao da variancia (R2 marginal / aleatorio / residual)
res_var_fert_export <- res_var_fert %>%
  mutate(
    R2_marginal_pct = round(Fixo * 100, 2),
    R2_aleatorio_pct = round(Aleatorio * 100, 2),
    R2_residual_pct = round(Residual * 100, 2)
  )


# #############################################################################
# PARTE D: GRADIENTE ASSOCIADO AO DOSSEL
# #############################################################################

cat("\n>>> PARTE D: Gradiente associado ao dossel <<<\n")

# Extrair Scores do PC2
scores_cd <- pca_result$ind$coord[,2]
indice_cd <- scales::rescale(scores_cd, to = c(0, 100))

idx_originais <- as.numeric(rownames(pca_result$ind$coord))

dados_cd <- data.frame(
  Area = area_final,
  Parcela = dados_solos$Parcela[idx_originais],
  Score_PC2 = scores_cd,
  Indice_CD = indice_cd
)

# ANOVA + Tukey
anova_cd <- aov(Indice_CD ~ factor(Area), data = dados_cd)
anova_summary_cd <- summary(anova_cd)
anova_pvalue_cd <- anova_summary_cd[[1]][["Pr(>F)"]][1]
tukey_cd <- TukeyHSD(anova_cd)

anova_cd_export <- data.frame(
  axis = "Canopy-associated gradient (PC2)",
  df_between = anova_summary_cd[[1]][["Df"]][1],
  df_within = anova_summary_cd[[1]][["Df"]][2],
  F = anova_summary_cd[[1]][["F value"]][1],
  p = anova_pvalue_cd
)
tukey_cd_matrix <- tukey_cd$`factor(Area)`
tukey_cd_export <- data.frame(
  axis = "Canopy-associated gradient (PC2)",
  comparison = rownames(tukey_cd_matrix),
  difference = tukey_cd_matrix[, "diff"],
  lower = tukey_cd_matrix[, "lwr"],
  upper = tukey_cd_matrix[, "upr"],
  p_adjusted = tukey_cd_matrix[, "p adj"],
  row.names = NULL
)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(
  rbind(anova_fert_export, anova_cd_export),
  "results/tables/environmental_axis_anova.csv",
  row.names = FALSE
)
write.csv(
  rbind(tukey_fert_export, tukey_cd_export),
  "results/tables/environmental_axis_tukey.csv",
  row.names = FALSE
)


# #############################################################################
# PARTE E: MODELOS MISTOS PARA O GRADIENTE ASSOCIADO AO DOSSEL
# #############################################################################

cat("\n>>> PARTE E: Modelos mistos - gradiente associado ao dossel <<<\n")

# Preparar join
criar_id_cd <- function(a, p) {
  clean <- str_replace_all(p, "-", "_")
  paste0(a, "_", str_extract(clean, "\\d+"), str_extract(clean, "[DE]"))
}

part_cd <- dados_cd %>%
  mutate(Parcela_ID = criar_id_cd(Area, Parcela)) %>%
  select(Parcela_ID, Area, Indice_CD)

part_taxas_cd <- taxas_medias_parcela %>%
  mutate(Parcela_ID = str_replace(Parcela, "^PNSJ", "")) %>%
  select(Parcela_ID, starts_with("taxa_"))

model_data_cd <- inner_join(part_cd, part_taxas_cd, by = "Parcela_ID") %>% na.omit()
cat("  Dataset integrado:", nrow(model_data_cd), "parcelas\n")

# Ajustar modelos lineares mistos
vars_resp_cd <- c("taxa_mortalidade_media", "taxa_recrutamento_media",
                  "taxa_perda_ab_media", "taxa_ganho_ab_media",
                  "taxa_mudanca_liq_media", "taxa_mudanca_ab_media")

nomes_respostas_cd <- c("Mortalidade", "Recrutamento", "Perda AB",
                        "Ganho AB", "Mudanca Ind.", "Mudanca AB")

res_coefs_cd <- data.frame()
res_var_cd <- data.frame()

for(i in 1:length(vars_resp_cd)) {
  var <- vars_resp_cd[i]
  nome <- nomes_respostas_cd[i]

  mod <- lmer(as.formula(paste0(var, " ~ Indice_CD + (1|Area)")),
              data = model_data_cd, REML = TRUE)

  tidied <- broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE) %>%
    filter(term == "Indice_CD")

  sig <- ifelse(tidied$p.value < 0.05, "Sim", "Nao")
  r2 <- performance::r2(mod)

  res_coefs_cd <- rbind(res_coefs_cd, data.frame(
    Modelo = nome, Coef = tidied$estimate, Erro = tidied$std.error,
    IC_inf = tidied$conf.low, IC_sup = tidied$conf.high,
    P_val = tidied$p.value, Significativo = sig
  ))

  res_var_cd <- rbind(res_var_cd, data.frame(
    Modelo = nome,
    Fixo = r2$R2_marginal,
    Aleatorio = r2$R2_conditional - r2$R2_marginal,
    Residual = 1 - r2$R2_conditional
  ))
}

# Decomposicao da variancia (R2 marginal / aleatorio / residual)
res_var_cd_export <- res_var_cd %>%
  mutate(
    R2_marginal_pct = round(Fixo * 100, 2),
    R2_aleatorio_pct = round(Aleatorio * 100, 2),
    R2_residual_pct = round(Residual * 100, 2)
  )

cat("\n============================================================\n")
cat("ANALISE CONCLUIDA COM SUCESSO!\n")
cat("============================================================\n")
cat("Resultados mantidos em memoria para as tabelas e figuras finais.\n")
cat("============================================================\n")
