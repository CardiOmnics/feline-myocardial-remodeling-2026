# ==============================================================================
# Quantitative Histopathology Analysis
# ==============================================================================
#
# This script includes all analyses for the feline heart failure study manuscript.
#
# Analyses included:
#   1. AI pixel classifier validation against human observers (fibrosis & adipocyte)
#   2. Myocardial fibrosis quantification (Masson's trichrome staining)
#   3. Adipocyte size and percentage quantification (H&E staining)
#   4. Intramural vessel morphometry (wall-to-lumen ratio, wall thickness, WCSA)
#   5. Cardiomyocyte nucleus shape analysis (mixed-effects models)
#   6. Assembly of per-cat phenotype summary file
#   7. Phenotype correlation analysis (Pearson/Spearman with bootstrap CIs)
#   8. Nanopore transcriptome correlation and GO pathway enrichment (GSEA)
#
# Input:  QuPath measurements, metadata files, Nanopore featureCounts tables,
#         human vs AI validation data
# Output: CSV statistics files, PDF figures, merged phenotype summary file
#
# Cats: 38 cats total. Cat 8 excluded from all analyses because of cardiac neoplasia.
# Disease groups (3-group): "No", "Yes_ATE" (HCM with ATE),
#                           "Yes_Congestive_HF" (HCM with congestive heart failure)
# Known cardiac death (2-group): "Yes", "No"
#
# ==============================================================================


# ==============================================================================
# SECTION 1: LIBRARIES
# ==============================================================================

# --- Genomics / pathway analysis (load before dplyr to avoid masking) ---
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(STRINGdb)
library(fgsea)

# --- Statistical testing ---
library(car)
library(rstatix)
library(lme4)
library(lmerTest)
library(nlme)
library(emmeans)
library(nortest)
library(FSA)
library(DHARMa)
library(pROC)
library(broom)

# --- Visualization ---
library(ggplot2)
library(gridExtra)
library(ggtext)
library(ggcorrplot)
library(scales)
library(ggpubr)
library(ggridges)
library(forcats)

# --- Core data manipulation (load dplyr LAST so select/filter/rename win) ---
library(tidyr)
library(reshape2)
library(data.table)
library(readxl)
library(openxlsx)
library(stringr)
library(purrr)
library(dplyr)

# Resolve common namespace conflicts: ensure dplyr versions win
conflictRules("dplyr", exclude = NULL)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate
summarise <- dplyr::summarise

set.seed(123)


# ==============================================================================
# SECTION 2: CONFIGURATION
# ==============================================================================
# Set the working directory to the project root before running this script.
# Expected folder structure:
#   data/metadata/    - Metadata_HE.txt, Metadata_MT.txt, Clinical_data_table_1.txt
#   data/qupath/      - QuPath measurement CSVs
#   data/nanopore/    - featureCounts .tabular files
#   data/validation/  - Human vs AI comparison CSV
#   results/          - Output directory for CSVs and PDFs

prefix <- format(Sys.Date(), "%d.%m.%Y")

# --- Relative paths (for GitHub) ---
# output          <- "results/"
# metadata_path   <- "data/metadata/"
# qupath_path     <- "data/qupath/"
# nanopore_path   <- "data/nanopore/"
# validation_path <- "data/validation/"

# --- Absolute paths  ---
base_path <- "specify path"
output          <- file.path(base_path, "Github/results/")
metadata_path   <- file.path(base_path, "Metadata/")
qupath_path     <- file.path(base_path, "QuPath measurements/")
nanopore_path   <- file.path(base_path, "Nanopore/featureCounts/")
validation_path <- file.path(base_path, "QuPath measurements/")

# Create output directory if it doesn't exist
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Cats to include (cat 8 excluded from all analyses)
INCLUDE_CATS <- c(1:7, 9:38)

# Adipocyte analyses also exclude cat 9 (insufficient tissue for adipocyte quantification)
INCLUDE_CATS_ADIPO <- c(1:7, 10:38)

# Fibrosis analyses exclude additional cats with missing tissue sections
INCLUDE_CATS_FIBROSIS <- c(1:7, 9:25, 27:30, 32:36, 38)

# Helper: standardize ATE labels across all data sources
standardize_ATE <- function(ate_values) {
  dplyr::recode(as.character(ate_values),
    "Yes:ATE"           = "Yes_ATE",
    "Yes:Congestive HF" = "Yes_Congestive_HF",
    "Yes:CHF"           = "Yes_Congestive_HF",
    "HCM with ATE"      = "Yes_ATE",
    "HCM only"          = "Yes_Congestive_HF",
    "No ATE"            = "No",
    "Yes_ATE"           = "Yes_ATE",
    "Yes_Congestive_HF" = "Yes_Congestive_HF",
    "No"                = "No",
    .default            = as.character(ate_values)
  )
}

# Helper: standardize Known cardiac death labels (0/1 -> "No"/"Yes")
standardize_KCD <- function(kcd_values) {
  ifelse(as.numeric(kcd_values) == 1, "Yes", "No")
}


# ==============================================================================
# SECTION 3: SHARED HELPER FUNCTIONS
# ==============================================================================

# Pairwise ATE 3-group comparison (Shapiro-Wilk -> Levene -> t-test/Wilcoxon)
# Used for fibrosis and nuclei analyses.
# Groups: "Yes_Congestive_HF", "Yes_ATE", "No"
# P-value correction: Benjamini-Hochberg
analyze_pairwise <- function(data, params, output, prefix, filename) {

  all_results <- list()

  for (param in params) {

    group_HCMonly     <- na.omit(data[[param]][data$ATE == "Yes_Congestive_HF"])
    group_HCMwithATE  <- na.omit(data[[param]][data$ATE == "Yes_ATE"])
    group_NoATE       <- na.omit(data[[param]][data$ATE == "No"])

    # Check normality (Shapiro-Wilk requires n >= 3)
    safe_shapiro <- function(x) if (length(x) >= 3) shapiro.test(x)$p.value else NA
    shapiro_HCMonly     <- safe_shapiro(group_HCMonly)
    shapiro_HCMwithATE  <- safe_shapiro(group_HCMwithATE)
    shapiro_NoATE       <- safe_shapiro(group_NoATE)
    shapiro_pvals <- c(shapiro_HCMonly, shapiro_HCMwithATE, shapiro_NoATE)
    normal_all <- all(shapiro_pvals > 0.05, na.rm = TRUE) & !any(is.na(shapiro_pvals))

    # Check homogeneity of variances (var.test requires n >= 2 in each group)
    safe_vartest <- function(x, y) {
      if (length(x) >= 2 & length(y) >= 2) var.test(x, y)$p.value else NA
    }
    levene_HCMonly_HCMwithATE <- safe_vartest(group_HCMonly, group_HCMwithATE)
    levene_HCMonly_NoATE      <- safe_vartest(group_HCMonly, group_NoATE)
    levene_HCMwithATE_NoATE   <- safe_vartest(group_HCMwithATE, group_NoATE)
    levene_pvals <- c(levene_HCMonly_HCMwithATE, levene_HCMonly_NoATE, levene_HCMwithATE_NoATE)
    levene_all <- all(levene_pvals > 0.05, na.rm = TRUE) & !any(is.na(levene_pvals))

    use_parametric <- normal_all & levene_all

    comparisons <- list(
      "HCMonly vs HCMwithATE" = list(group_HCMonly, group_HCMwithATE),
      "HCMonly vs NoATE"      = list(group_HCMonly, group_NoATE),
      "HCMwithATE vs NoATE"   = list(group_HCMwithATE, group_NoATE)
    )

    results_list <- list()

    for (comp in names(comparisons)) {
      g1 <- comparisons[[comp]][[1]]
      g2 <- comparisons[[comp]][[2]]

      if (length(g1) > 1 & length(g2) > 1) {
        if (use_parametric) {
          test_result <- t.test(g1, g2, var.equal = TRUE)
          method <- "Unpaired t-Test"
          ci_lower <- round(test_result$conf.int[1], 3)
          ci_upper <- round(test_result$conf.int[2], 3)
        } else {
          test_result <- wilcox.test(g1, g2)
          method <- "Wilcoxon Rank-Sum Test"
          ci_lower <- NA
          ci_upper <- NA
        }

        levene_value <- switch(comp,
          "HCMonly vs HCMwithATE" = levene_HCMonly_HCMwithATE,
          "HCMonly vs NoATE"      = levene_HCMonly_NoATE,
          "HCMwithATE vs NoATE"   = levene_HCMwithATE_NoATE,
          NA
        )

        results_list[[comp]] <- data.frame(
          Comparison = comp,
          Parameter = param,
          Method = method,
          Mean_Group_1 = round(mean(g1, na.rm = TRUE), 3),
          Mean_Group_2 = round(mean(g2, na.rm = TRUE), 3),
          Median_Group_1 = round(median(g1, na.rm = TRUE), 3),
          Median_Group_2 = round(median(g2, na.rm = TRUE), 3),
          Test.Statistic = round(test_result$statistic, 3),
          n_Group_1 = length(g1),
          n_Group_2 = length(g2),
          Normality_p_Group_1 = if (length(g1) >= 3) signif(shapiro.test(g1)$p.value, 3) else NA,
          Normality_p_Group_2 = if (length(g2) >= 3) signif(shapiro.test(g2)$p.value, 3) else NA,
          Levene_p = signif(levene_value, 3),
          p.value = signif(test_result$p.value, 3),
          CI_Lower_95 = ci_lower,
          CI_Upper_95 = ci_upper
        )
      }
    }

    if (length(results_list) > 0) {
      param_results <- do.call(rbind, results_list)
      param_results$p.value <- as.numeric(param_results$p.value)
      valid_pvalues <- na.omit(param_results$p.value)

      if (length(valid_pvalues) > 0) {
        param_results$Adjusted_p_BH <- NA
        param_results$Adjusted_p_BH[!is.na(param_results$p.value)] <- p.adjust(valid_pvalues, method = "BH")
      } else {
        param_results$Adjusted_p_BH <- NA
      }

      all_results[[param]] <- param_results
    }
  }

  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    write.csv(combined_results,
              file = paste0(output, prefix, "_", filename),
              row.names = FALSE)
    print(combined_results)
  } else {
    message("No valid comparisons were found.")
  }
}

# Two-group comparison: Known cardiac death (Yes/No)
# Shapiro-Wilk -> Levene -> t-test/Wilcoxon
analyze_two_groups <- function(data, params, group_var, group_levels, output, prefix, filename) {
  all_results <- list()

  for (param in params) {
    group_0 <- na.omit(data[[param]][data[[group_var]] == group_levels[1]])
    group_1 <- na.omit(data[[param]][data[[group_var]] == group_levels[2]])

    shapiro_0 <- if (length(group_0) >= 3) shapiro.test(group_0) else list(p.value = NA)
    shapiro_1 <- if (length(group_1) >= 3) shapiro.test(group_1) else list(p.value = NA)
    levene_test <- leveneTest(as.formula(paste(param, "~", group_var)), data = data)

    normal_both <- !is.na(shapiro_0$p.value) & !is.na(shapiro_1$p.value) &
      shapiro_0$p.value > 0.05 & shapiro_1$p.value > 0.05
    homoscedastic <- levene_test$`Pr(>F)`[1] > 0.05

    if (normal_both & homoscedastic) {
      test_result <- t.test(as.formula(paste(param, "~", group_var)),
                            data = data, var.equal = TRUE, na.action = na.omit)
      method <- "Unpaired t-Test"
    } else {
      test_result <- wilcox.test(as.formula(paste(param, "~", group_var)),
                                  data = data, na.action = na.omit)
      method <- "Wilcoxon Rank-Sum Test (Mann-Whitney U)"
    }

    all_results[[param]] <- data.frame(
      Parameter = param, Method = method,
      Mean_No = round(mean(group_0, na.rm = TRUE), 3),
      Mean_Yes = round(mean(group_1, na.rm = TRUE), 3),
      Median_No = round(median(group_0, na.rm = TRUE), 3),
      Median_Yes = round(median(group_1, na.rm = TRUE), 3),
      Test.Statistic = round(test_result$statistic, 3),
      p.value = signif(test_result$p.value, 3),
      CI_Lower_95 = ifelse(method == "Unpaired t-Test", round(test_result$conf.int[1], 3), NA),
      CI_Upper_95 = ifelse(method == "Unpaired t-Test", round(test_result$conf.int[2], 3), NA),
      Normality_p_No = round(shapiro_0$p.value, 3),
      Normality_p_Yes = round(shapiro_1$p.value, 3),
      Variance_Homogeneity_p = round(levene_test$`Pr(>F)`[1], 3)
    )
  }

  combined_results <- do.call(rbind, all_results)
  write.csv(combined_results,
            file = file.path(output, paste0(prefix, "_", filename)),
            row.names = FALSE)
  print(combined_results)
}


# ==============================================================================
# SECTION 4: CLASSIFIER VALIDATION (AI vs Human Observers)
# ==============================================================================
# Compares AI pixel classifier performance against lab observers and a
# pathologist for fibrosis and adipocyte quantification.

# Load validation data
validation_data <- read.csv(
  file.path(validation_path, "Fibrosis_quantification_ Human vs AI_270225.csv"),
  dec = ","
)
rownames(validation_data) <- validation_data$Observer
validation_data$Observer <- NULL

# Split into fibrosis and adipocyte datasets
fibrosis_wide  <- validation_data[, 1:6]
adipocyte_wide <- validation_data[, 7:12]
fibrosis_wide$observer  <- rownames(validation_data)
adipocyte_wide$observer <- rownames(validation_data)

# --- Helper functions ---

make_long <- function(df_wide) {
  df_wide %>%
    pivot_longer(cols = -observer, names_to = "variable", values_to = "value") %>%
    mutate(
      observer = as.character(observer),
      variable = as.character(variable),
      value = as.numeric(value)
    )
}

analyze_model <- function(long_df, title_text, y_label, pdf_name, csv_name,
                          ai_name = "AI", pathologist_name = "Pathologist") {

  lab_df <- long_df %>% filter(observer != ai_name, observer != pathologist_name)
  ai_df  <- long_df %>% filter(observer == ai_name)
  pathologist_df <- long_df %>% filter(observer == pathologist_name)

  # Summary per measurement area from lab observers
  lab_summary <- lab_df %>%
    group_by(variable) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value   = sd(value, na.rm = TRUE),
      min_value  = min(value, na.rm = TRUE),
      max_value  = max(value, na.rm = TRUE),
      .groups = "drop"
    )

  # Sort by mean lab score
  variable_order <- lab_summary %>% arrange(mean_value) %>% pull(variable)

  lab_summary <- lab_summary %>%
    mutate(variable = factor(variable, levels = variable_order)) %>%
    arrange(variable) %>%
    mutate(variable_numeric = seq_len(n()))

  lab_df <- lab_df %>% mutate(variable = factor(variable, levels = variable_order))
  ai_df  <- ai_df %>%
    mutate(variable = factor(variable, levels = variable_order)) %>%
    arrange(variable) %>%
    mutate(variable_numeric = seq_len(n()))
  pathologist_df <- pathologist_df %>%
    mutate(variable = factor(variable, levels = variable_order)) %>%
    arrange(variable) %>%
    mutate(variable_numeric = seq_len(n()))

  # Human observer variability relative to lab mean
  human_errors <- lab_df %>%
    left_join(lab_summary %>% dplyr::select(variable, mean_value), by = "variable") %>%
    group_by(observer) %>%
    summarise(
      MAE  = mean(abs(value - mean_value), na.rm = TRUE),
      RMSE = sqrt(mean((value - mean_value)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  avg_human_mae  <- mean(human_errors$MAE, na.rm = TRUE)
  avg_human_rmse <- mean(human_errors$RMSE, na.rm = TRUE)

  # AI performance relative to lab mean
  ai_scores <- ai_df %>%
    dplyr::select(variable, ai_value = value) %>%
    left_join(lab_summary %>% dplyr::select(variable, mean_value), by = "variable") %>%
    mutate(diff = ai_value - mean_value)

  ai_mae  <- mean(abs(ai_scores$diff), na.rm = TRUE)
  ai_rmse <- sqrt(mean(ai_scores$diff^2, na.rm = TRUE))
  ai_bias <- mean(ai_scores$diff, na.rm = TRUE)

  # Summary statistics table
  stat_results <- data.frame(
    Section = title_text,
    Comparison = c("Lab observer mean", "Pixel classifier",
                   "Lab observer mean", "Pixel classifier",
                   "Pixel classifier vs lab observer mean",
                   "Pixel classifier vs lab observer mean",
                   "Pixel classifier vs lab observer mean",
                   "Lab observers vs lab observer mean",
                   "Lab observers vs lab observer mean"),
    Metric = c("Mean percentage", "Mean percentage",
               "Median percentage", "Median percentage",
               "MAE", "RMSE", "Bias",
               "Average observer MAE", "Average observer RMSE"),
    Value = round(c(
      mean(ai_scores$mean_value, na.rm = TRUE),
      mean(ai_scores$ai_value, na.rm = TRUE),
      median(ai_scores$mean_value, na.rm = TRUE),
      median(ai_scores$ai_value, na.rm = TRUE),
      ai_mae, ai_rmse, ai_bias,
      avg_human_mae, avg_human_rmse
    ), 3)
  )

  write.csv(stat_results,
            file = file.path(output, paste0(prefix, "_", csv_name)),
            row.names = FALSE)

  # Plot
  annotation_label <- paste0("MAE = ", round(ai_mae, 2),
                             "\nRMSE = ", round(ai_rmse, 2),
                             "\nBias = ", round(ai_bias, 2))

  p <- ggplot() +
    geom_ribbon(data = lab_summary,
                aes(x = variable_numeric, ymin = min_value, ymax = max_value),
                fill = "orange", alpha = 0.2) +
    geom_line(data = lab_summary,
              aes(x = variable_numeric, y = mean_value, color = "Lab observers"),
              linewidth = 1) +
    geom_line(data = ai_df,
              aes(x = variable_numeric, y = value, color = "Pixel classifier"),
              linewidth = 1.2) +
    geom_line(data = pathologist_df,
              aes(x = variable_numeric, y = value, color = "Pathologist"),
              linewidth = 1.2) +
    annotate("text", x = 1, y = max(lab_summary$max_value, na.rm = TRUE) * 0.95,
             label = annotation_label, hjust = 0, size = 4, color = "black") +
    scale_x_continuous(breaks = lab_summary$variable_numeric,
                       labels = as.character(lab_summary$variable)) +
    scale_color_manual(values = c("Lab observers" = "orange",
                                  "Pixel classifier" = "blue",
                                  "Pathologist" = "red")) +
    labs(title = title_text,
         x = "Measurement area (sorted by mean lab observer percentage)",
         y = y_label, color = "Legend") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

  pdf(file.path(output, paste0(prefix, "_", pdf_name)), width = 6, height = 5)
  print(p)
  dev.off()

  print(p)
  print(stat_results)
  print(human_errors)

  invisible(list(plot = p, stat_results = stat_results,
                 human_errors = human_errors, lab_summary = lab_summary,
                 ai_scores = ai_scores))
}

# --- Run fibrosis validation ---
fibrosis_long <- make_long(fibrosis_wide)
fibrosis_results <- analyze_model(
  long_df    = fibrosis_long,
  title_text = "Fibrosis quantification: Lab observers vs pixel classifier",
  y_label    = "Fibrosis percentage",
  pdf_name   = "Fibrosis_model_performance.pdf",
  csv_name   = "Fibrosis_model_performance.csv"
)

# --- Run adipocyte validation ---
adipocyte_long <- make_long(adipocyte_wide)
adipocyte_results <- analyze_model(
  long_df    = adipocyte_long,
  title_text = "Adipocyte quantification: Lab observers vs pixel classifier",
  y_label    = "Adipocyte percentage",
  pdf_name   = "Adipocyte_model_performance.pdf",
  csv_name   = "Adipocyte_model_performance.csv"
)


# ==============================================================================
# SECTION 5: FIBROSIS ANALYSIS
# ==============================================================================
# Quantifies myocardial and total fibrosis from Masson's trichrome staining.
# Compares disease groups (Known cardiac death, ATE status).

# --- Load data ---
metadata_MT <- read.delim(file.path(metadata_path, "Metadata_MT.txt"))
measurements_without_borders <- read.csv(file.path(qupath_path, "measurements_fibrosis_without_borders_281024.csv"))
measurements_whole_heart     <- read.csv(file.path(qupath_path, "measurements_fibrosis_whole_heart_281024.csv"))
clinical_data <- read.delim(file.path(metadata_path, "260225_Clinical_data_table_1.txt"))
clinical_data$ATE <- standardize_ATE(clinical_data$ATE)

# --- Calculate fibrosis percentages ---
metadata_MT <- metadata_MT %>%
  mutate(image_id = sub("_.*", "", Image))

# Whole heart fibrosis (total)
measurements_whole_heart <- measurements_whole_heart %>%
  mutate(image_id = sub("\\.ndpi", "", Image)) %>%
  left_join(metadata_MT %>% dplyr::select(image_id, Cat_number, Tissue_part), by = "image_id")

fibrosis_total <- measurements_whole_heart %>%
  group_by(Cat_number) %>%
  summarise(
    Total_Fibrosis_Area = sum(Cardiomyocyte_fibrosis_210824..Fibrosis.area.µm.2, na.rm = TRUE),
    Total_Cardiomyocyte_Area = sum(Cardiomyocyte_fibrosis_210824..Cardiomyocyte.area.µm.2, na.rm = TRUE)
  ) %>%
  mutate(Mean_Fibrosis_Percentage = Total_Fibrosis_Area /
           (Total_Fibrosis_Area + Total_Cardiomyocyte_Area) * 100)

# Myocardial fibrosis (without epicardial/endocardial borders)
measurements_without_borders <- measurements_without_borders %>%
  mutate(image_id = sub("\\.ndpi", "", Image)) %>%
  left_join(metadata_MT %>% dplyr::select(image_id, Cat_number, Tissue_part), by = "image_id")

fibrosis_myocardial <- measurements_without_borders %>%
  group_by(Cat_number) %>%
  summarise(
    Total_Fibrosis_Area = sum(Cardiomyocyte_fibrosis_210824..Fibrosis.area.µm.2, na.rm = TRUE),
    Total_Cardiomyocyte_Area = sum(Cardiomyocyte_fibrosis_210824..Cardiomyocyte.area.µm.2, na.rm = TRUE)
  ) %>%
  mutate(Mean_Myocardial_Fibrosis_Percentage = Total_Fibrosis_Area /
           (Total_Fibrosis_Area + Total_Cardiomyocyte_Area) * 100)

# Combine and filter cats
measurements_mean <- fibrosis_total %>%
  left_join(
    fibrosis_myocardial %>%
      dplyr::select(Cat_number, Mean_Myocardial_Fibrosis_Percentage),
    by = "Cat_number"
  )

measurements_mean <- measurements_mean[measurements_mean$Cat_number %in% INCLUDE_CATS_FIBROSIS, ]

# Calculate non-myocardial fibrosis
measurements_mean$Mean_NonMyocardial_Fibrosis_Percentage <-
  measurements_mean$Mean_Fibrosis_Percentage - measurements_mean$Mean_Myocardial_Fibrosis_Percentage

# --- Add clinical data for cardiac death and ATE status ---
clinical_data$Known_cardiac_death <- standardize_KCD(clinical_data$Known_cardiac_death)

measurements_mean$Cat_number <- as.numeric(measurements_mean$Cat_number)
measurements_mean <- measurements_mean %>%
  left_join(
    clinical_data %>% dplyr::select(Cat_number, Known_cardiac_death),
    by = "Cat_number"
  )

# --- Bar plots ---

# Total vs myocardial fibrosis bar plot
cat1_40 <- measurements_mean[measurements_mean$Cat_number %in% c(1:38), ]

pdf(file.path(output, paste0(prefix, "_Mean_fibrosis_percentage.pdf")), width = 5, height = 3.5)
p <- ggplot(cat1_40, aes(x = reorder(Cat_number, Mean_Myocardial_Fibrosis_Percentage))) +
  geom_col(aes(y = Mean_Fibrosis_Percentage, fill = "Mean_Fibrosis_Percentage"), position = "dodge") +
  geom_col(aes(y = Mean_Myocardial_Fibrosis_Percentage, fill = "Mean_Myocardial_Fibrosis_Percentage"), position = "dodge") +
  xlab("Cat number") + ylab("Fibrosis (%)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
        axis.line = element_line(size = 0.2, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.text = element_text(size = 8), legend.title = element_blank()) +
  scale_fill_manual(values = c("Mean_Fibrosis_Percentage" = "lightgrey",
                                "Mean_Myocardial_Fibrosis_Percentage" = "grey10"))
print(p)
dev.off()

# Bar plot colored by cardiac death status
cat1_40$Cat_label <- factor(cat1_40$Cat_number)

pdf(file.path(output, paste0(prefix, "_Mean_fibrosis_percentage_cardiac_death.pdf")), width = 5, height = 3.5)
p <- ggplot(cat1_40, aes(x = reorder(Cat_label, Mean_Myocardial_Fibrosis_Percentage))) +
  geom_col(aes(y = Mean_Fibrosis_Percentage, fill = "Mean_Fibrosis_Percentage"), position = "dodge") +
  geom_col(aes(y = Mean_Myocardial_Fibrosis_Percentage, fill = "Mean_Myocardial_Fibrosis_Percentage"), position = "dodge") +
  xlab("Cat number") + ylab("Fibrosis (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(size = 6, angle = 90, hjust = 1, face = "bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
    axis.line = element_line(size = 0.2, colour = "black"),
    legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
    legend.text = element_text(size = 8), legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("Mean_Fibrosis_Percentage" = "grey70",
                                "Mean_Myocardial_Fibrosis_Percentage" = "grey10")) +
  scale_x_discrete(labels = function(x) {
    ifelse(cat1_40$Known_cardiac_death[match(x, cat1_40$Cat_label)] == "Yes",
           paste0("<span style='color: red;'>", x, "</span>"),
           paste0("<span style='color: black;'>", x, "</span>"))
  })
print(p)
dev.off()

# Bar plot colored by ATE status
clinical_data$Cat_number <- as.numeric(clinical_data$Cat_number)
cat1_40$Cat_number <- as.numeric(cat1_40$Cat_number)
cat1_40 <- cat1_40 %>%
  left_join(clinical_data %>% dplyr::select(Cat_number, ATE), by = "Cat_number")

pdf(file.path(output, paste0(prefix, "_Mean_fibrosis_percentage_ATE.pdf")), width = 5, height = 3.5)
p <- ggplot(cat1_40, aes(x = reorder(Cat_label, Mean_Myocardial_Fibrosis_Percentage))) +
  geom_col(aes(y = Mean_Fibrosis_Percentage, fill = "Mean_Fibrosis_Percentage"), position = "dodge") +
  geom_col(aes(y = Mean_Myocardial_Fibrosis_Percentage, fill = "Mean_Myocardial_Fibrosis_Percentage"), position = "dodge") +
  xlab("Cat number") + ylab("Fibrosis (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_markdown(size = 6, angle = 90, hjust = 1, face = "bold"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
    axis.line = element_line(size = 0.2, colour = "black"),
    legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
    legend.text = element_text(size = 8), legend.title = element_blank()
  ) +
  scale_fill_manual(values = c("Mean_Fibrosis_Percentage" = "grey70",
                                "Mean_Myocardial_Fibrosis_Percentage" = "grey10")) +
  scale_x_discrete(labels = function(x) {
    ate_value <- cat1_40$ATE[match(x, cat1_40$Cat_label)]
    color <- ifelse(ate_value == "No", "black",
                    ifelse(ate_value == "Yes_ATE", "red", "green"))
    paste0("<span style='color: ", color, ";'>", x, "</span>")
  })
print(p)
dev.off()

# --- Statistics: Known cardiac death (2-group comparison) ---
measurements_mean_pheno <- measurements_mean
measurements_mean_pheno$Known_cardiac_death <- as.factor(measurements_mean_pheno$Known_cardiac_death)

fibrosis_params <- c("Mean_Myocardial_Fibrosis_Percentage",
                     "Mean_Fibrosis_Percentage",
                     "Mean_NonMyocardial_Fibrosis_Percentage")

analyze_two_groups(data = measurements_mean_pheno, params = fibrosis_params,
                   group_var = "Known_cardiac_death", group_levels = c("No", "Yes"),
                   output = output, prefix = prefix,
                   filename = "All_Fibrosis_Perc_Results.csv")

# --- Statistics: ATE (3-group pairwise comparison) ---
measurements_mean_pheno2 <- measurements_mean
clinical_data$Cat_number <- as.numeric(clinical_data$Cat_number)
measurements_mean_pheno2$Cat_number <- as.numeric(measurements_mean_pheno2$Cat_number)
measurements_mean_pheno2 <- measurements_mean_pheno2 %>%
  left_join(clinical_data %>% dplyr::select(Cat_number, ATE), by = "Cat_number")

for (parameter in fibrosis_params) {
  analyze_pairwise(
    data = measurements_mean_pheno2,
    params = parameter,
    output = output,
    prefix = prefix,
    filename = paste(parameter, "Fibrosis_Pairwise_Stats.csv", sep = "_")
  )
}

# --- Boxplots for ATE groups ---
fibrosis_y_labels <- c(
  "Mean_Fibrosis_Percentage" = "Total fibrosis %",
  "Mean_Myocardial_Fibrosis_Percentage" = "Myocardial fibrosis %",
  "Mean_NonMyocardial_Fibrosis_Percentage" = "Non-myocardial fibrosis %"
)

color_mapping_ate <- c("No" = "blue",
                       "Yes_ATE" = "red",
                       "Yes_Congestive_HF" = "red")

boxplot_list <- list()
for (parameter in fibrosis_params) {
  y_label <- fibrosis_y_labels[[parameter]]
  p <- ggplot(measurements_mean_pheno2, aes(x = factor(ATE), y = .data[[parameter]], colour = ATE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2) +
    scale_color_manual(values = color_mapping_ate) +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    labs(x = "Known cardiac death", y = y_label) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 14, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 16, face = "bold", color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  boxplot_list[[parameter]] <- p
}

pdf(file.path(output, paste0(prefix, "_Fibrosis_measurements_ATE.pdf")), width = 8, height = 5)
grid.arrange(grobs = boxplot_list, nrow = 1, ncol = 3)
dev.off()

# --- Export fibrosis per cat ---
write.csv(measurements_mean,
          file = paste0(output, prefix, "_Fibrosis_percentage_per_cat.csv"),
          row.names = FALSE)


# ==============================================================================
# SECTION 6: ADIPOCYTE ANALYSIS
# ==============================================================================
# Quantifies adipocyte size and percentage from H&E staining.
# Compares left vs right ventricle and disease groups.

# --- Load data ---
metadata_HE <- read.delim(file.path(metadata_path, "Metadata_HE.txt"))
Adipo_size <- read.csv(file.path(qupath_path, "Adipocyte_size_271124.csv"))
Adipo_perc <- read.csv(file.path(qupath_path, "Adipocyte_percentage_141124.csv"))

# --- Prepare metadata ---
metadata_HE <- metadata_HE %>%
  mutate(image_id = sub("\\.ndpi", "", image_id))

# --- Adipocyte size ---
Adipo_size <- Adipo_size %>%
  mutate(image_id = sub("\\.ndpi", "", Image)) %>%
  left_join(metadata_HE %>% dplyr::select(image_id, Cat_number, Tissue_part), by = "image_id")

# Remove cats 8 and 9
Adipo_size <- Adipo_size[Adipo_size$Cat_number %in% INCLUDE_CATS_ADIPO, ]

# Filter by circularity and area thresholds
Adipo_size_filtered <- Adipo_size %>%
  filter(Circularity > 0.35, Area.µm.2 < 3000, Area.µm.2 > 70)

Adipo_size_mean <- Adipo_size_filtered %>%
  group_by(Cat_number) %>%
  summarise(Mean_Adipocyte_size = mean(Area.µm.2, na.rm = TRUE))

Adipo_size_mean <- Adipo_size_mean[Adipo_size_mean$Cat_number %in% INCLUDE_CATS_ADIPO, ]

# --- Adipocyte percentage ---
Adipo_perc <- Adipo_perc %>%
  mutate(image_id = sub("\\.ndpi", "", Image)) %>%
  left_join(metadata_HE %>% dplyr::select(image_id, Cat_number, Tissue_part), by = "image_id")

# Combine area columns from two batches
Adipo_perc <- Adipo_perc %>%
  mutate(
    Tissue_area = ifelse(!is.na(Tissue_detection_quantification_130924..Tissue.area.µm.2),
                         Tissue_detection_quantification_130924..Tissue.area.µm.2,
                         Tissue_detection_quantification_131124..Tissue.area.µm.2),
    Adipocytes_area = ifelse(!is.na(Adipocyte_130924..Adipocytes.area.µm.2),
                             Adipocyte_130924..Adipocytes.area.µm.2,
                             Adipocyte_131124..Adipocytes.area.µm.2),
    Percentage = Adipocytes_area / (Tissue_area + Adipocytes_area) * 100
  )

Adipo_perc_mean <- Adipo_perc %>%
  group_by(Cat_number) %>%
  summarise(
    Total_Tissue_area = sum(Tissue_area, na.rm = TRUE),
    Total_Adipocytes_area = sum(Adipocytes_area, na.rm = TRUE)
  ) %>%
  mutate(Mean_Adipocyte_perc = Total_Adipocytes_area / (Total_Tissue_area + Total_Adipocytes_area) * 100)

Adipo_perc_mean <- Adipo_perc_mean[Adipo_perc_mean$Cat_number %in% INCLUDE_CATS_ADIPO, ]

# --- Remove size measurements from cats with very low fat (<3%) ---
Adipo_size_mean <- Adipo_size_mean %>%
  left_join(Adipo_perc_mean %>% dplyr::rename(percentage = Mean_Adipocyte_perc), by = "Cat_number")
Adipo_size_mean$Mean_Adipocyte_size[Adipo_size_mean$percentage < 3] <- NA
Adipo_size_mean <- Adipo_size_mean[, 1:2]

# --- Left vs Right ventricle comparison ---
paired_data <- Adipo_perc %>%
  filter(Tissue_part %in% c("R1", "R2", "L1", "L2")) %>%
  group_by(Cat_number, Tissue_part) %>%
  summarize(Mean_Percentage = mean(Percentage, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = Cat_number, names_from = Tissue_part, values_from = Mean_Percentage) %>%
  mutate(
    Right = rowMeans(dplyr::select(., R1, R2), na.rm = TRUE),
    Left  = rowMeans(dplyr::select(., L1, L2), na.rm = TRUE)
  ) %>%
  filter(!is.na(Right) & !is.na(Left))

wilcox_test_lr <- wilcox.test(paired_data$Right, paired_data$Left,
                              paired = TRUE, na.action = na.omit)

lr_test_results <- data.frame(
  Test_Type = "Wilcoxon Signed-Rank Test (WSRT)",
  Median_Right = round(median(paired_data$Right, na.rm = TRUE), 3),
  Median_Left = round(median(paired_data$Left, na.rm = TRUE), 3),
  Test_Statistic = wilcox_test_lr$statistic,
  p.value = signif(wilcox_test_lr$p.value, 3),
  n = nrow(paired_data)
)
print(lr_test_results)
write.csv(lr_test_results,
          file = paste0(output, prefix, "_Left_vs_Right_Adipocyte_Test_Results.csv"),
          row.names = FALSE)

# --- Statistics: ATE pairwise comparisons (BH correction) ---

# Percentage analysis
Adipo_perc_mean_pheno <- Adipo_perc_mean
clinical_data$Cat_number <- as.numeric(clinical_data$Cat_number)
Adipo_perc_mean_pheno$Cat_number <- as.numeric(Adipo_perc_mean_pheno$Cat_number)
Adipo_perc_mean_pheno <- Adipo_perc_mean_pheno %>%
  left_join(clinical_data %>% dplyr::select(Cat_number, ATE), by = "Cat_number")

analyze_pairwise(data = Adipo_perc_mean_pheno, params = c("Mean_Adipocyte_perc"),
                 output = output, prefix = prefix,
                 filename = "Pairwise_AdipocytePerc_Results.csv")

# Size analysis
Adipo_size_mean_pheno <- Adipo_size_mean
Adipo_size_mean_pheno$Cat_number <- as.numeric(Adipo_size_mean_pheno$Cat_number)
Adipo_size_mean_pheno <- Adipo_size_mean_pheno %>%
  left_join(clinical_data %>% dplyr::select(Cat_number, ATE), by = "Cat_number")

analyze_pairwise(data = Adipo_size_mean_pheno, params = c("Mean_Adipocyte_size"),
                 output = output, prefix = prefix,
                 filename = "Pairwise_AdipocyteSize_Results.csv")

# --- Export per-cat summaries ---
write.csv(Adipo_perc_mean,
          file = paste0(output, prefix, "_Adipocyte_percentage_per_cat.csv"),
          row.names = FALSE)
write.csv(Adipo_size_mean,
          file = paste0(output, prefix, "_Adipocyte_size_per_cat.csv"),
          row.names = FALSE)


# ==============================================================================
# SECTION 7: VESSEL ANALYSIS
# ==============================================================================
# Calculates vessel morphometry (OD, LD, WLR, WT, WCSA) and compares
# ATE groups stratified by vessel size.

# --- Load and calculate vessel metrics ---
Vessels <- read.csv(file.path(qupath_path, "Vessels_100125.csv"))
Vessels$Diameter_Type <- ifelse(Vessels$Parent == "Root object (Image)", "Outer Diameter", "Inner Diameter")
Vessels$Vessel_ID <- cumsum(Vessels$Diameter_Type == "Outer Diameter")

Vessels <- Vessels %>%
  mutate(
    Diameter = 2 * sqrt(Area.µm.2 / pi),
    OD = ifelse(Diameter_Type == "Outer Diameter", Diameter, NA),
    LD = ifelse(Diameter_Type == "Inner Diameter", Diameter, NA)
  )

Vessels_results <- Vessels %>%
  group_by(Vessel_ID, Image) %>%
  summarise(
    OD = max(OD, na.rm = TRUE),
    LD = max(LD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    WLR  = (OD - LD) / LD,
    WT   = (OD - LD) / 2,
    WCSA = pi / 4 * (OD^2 - LD^2)
  )

write.csv(Vessels_results,
          file = paste0(output, prefix, "_Vessel_Calculations.csv"),
          row.names = FALSE)

# --- Join with metadata ---
Vessels_results <- Vessels_results %>%
  mutate(image_id = sub("\\.ndpi", "", Image)) %>%
  left_join(metadata_MT %>% dplyr::select(image_id, Cat_number, Tissue_part), by = "image_id")

Vessels_results <- Vessels_results[Vessels_results$Cat_number %in% INCLUDE_CATS, ]

# Mean per cat
Vessel_measurements_mean <- Vessels_results %>%
  group_by(Cat_number) %>%
  summarise(Mean_WLR = mean(WLR, na.rm = TRUE),
            Mean_WT = mean(WT, na.rm = TRUE),
            Mean_WCSA = mean(WCSA, na.rm = TRUE))

write.csv(Vessel_measurements_mean,
          file = file.path(output, paste0(prefix, "_Vessel_measurements_per_cat.csv")),
          row.names = FALSE)

# --- Numbers for results section ---
Vessel_stats <- Vessels_results %>%
  group_by(Cat_number) %>%
  summarise(
    Vessels_detected = n(),
    Mean_WLR = mean(WLR, na.rm = TRUE), Variance_WLR = var(WLR, na.rm = TRUE),
    Mean_WT = mean(WT, na.rm = TRUE), Variance_WT = var(WT, na.rm = TRUE),
    Mean_WCSA = mean(WCSA, na.rm = TRUE), Variance_WCSA = var(WCSA, na.rm = TRUE)
  )

Overall_stats <- Vessel_stats %>%
  summarise(
    Min_vessels_per_cat = min(Vessels_detected),
    Max_vessels_per_cat = max(Vessels_detected),
    Mean_vessels_per_cat = mean(Vessels_detected),
    SD_vessels_per_cat = sd(Vessels_detected)
  )
print(Vessel_stats)
print(Overall_stats)

# --- Merge with clinical data for group comparisons ---
Vessels_results_pheno <- merge(
  Vessels_results,
  clinical_data[, c("Cat_number", "Known_cardiac_death", "ATE", "Age_years")],
  by = "Cat_number", all = FALSE
)
Vessels_results_pheno$Known_cardiac_death <- as.factor(Vessels_results_pheno$Known_cardiac_death)

cat("Overall mean WLR:", mean(Vessels_results_pheno$WLR), "\n")

# --- Distribution check: QQ plots and Gamma GLMM ---
qq_raw <- ggplot(Vessels_results_pheno, aes(sample = WLR)) +
  stat_qq() + stat_qq_line() + theme_minimal() +
  labs(title = "(A) QQ Plot of Raw WLR Data", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

glmm_model <- glmer(WLR ~ ATE + (1 | Cat_number),
                     data = Vessels_results_pheno, family = Gamma(link = "log"))
Vessels_results_pheno$residuals <- residuals(glmm_model)

qq_residuals <- ggplot(Vessels_results_pheno, aes(sample = residuals)) +
  stat_qq() + stat_qq_line() + theme_minimal() +
  labs(title = "(B) QQ Plot of Model Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

boxplot_resid <- ggplot(Vessels_results_pheno, aes(x = ATE, y = residuals, fill = ATE)) +
  geom_boxplot(outlier.shape = NA) + theme_minimal() +
  labs(title = "(C) Homogeneity of Variance in Residuals", x = "Disease Group", y = "Residuals") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

final_plot <- ggarrange(qq_raw, qq_residuals, boxplot_resid,
                        ncol = 2, nrow = 2, labels = c("A", "B", "C"), heights = c(1, 0.6))
pdf(file.path(output, paste0(prefix, "_Supplemental_QQ_Residuals.pdf")), width = 11, height = 7)
print(final_plot)
dev.off()

# --- Size stratification ---
Vessels_results_pheno$Size_Group <- cut(Vessels_results_pheno$OD,
    breaks = c(-Inf, 100, 250, Inf),
    labels = c("Small", "Medium", "Large"))

# Gamma GLMMs per size group
Vessels_results_pheno$ATE <- factor(Vessels_results_pheno$ATE)
Vessels_results_pheno$ATE <- relevel(Vessels_results_pheno$ATE, ref = "No")

glmm_small  <- glmer(WLR ~ ATE + (1 | Cat_number),
                      data = Vessels_results_pheno[Vessels_results_pheno$Size_Group == "Small", ],
                      family = Gamma(link = "log"))
glmm_medium <- glmer(WLR ~ ATE + (1 | Cat_number),
                      data = Vessels_results_pheno[Vessels_results_pheno$Size_Group == "Medium", ],
                      family = Gamma(link = "log"))
glmm_large  <- glmer(WLR ~ ATE + (1 | Cat_number),
                      data = Vessels_results_pheno[Vessels_results_pheno$Size_Group == "Large", ],
                      family = Gamma(link = "log"))

summary(glmm_small)
summary(glmm_medium)
summary(glmm_large)

emmeans(glmm_small, pairwise ~ ATE, adjust = "none")
emmeans(glmm_medium, pairwise ~ ATE, adjust = "none")
emmeans(glmm_large, pairwise ~ ATE, adjust = "none")

# --- Pairwise models with Kruskal-Wallis + LMM contrasts ---
analyze_ATE_groups_vessel <- function(data, size_group) {
  df <- data %>% filter(Size_Group == size_group)

  kruskal_result <- kruskal.test(WLR ~ ATE, data = df)
  kruskal_p <- kruskal_result$p.value
  n_vessels <- nrow(df)
  n_cats <- length(unique(df$Cat_number))

  results <- data.frame(
    Comparison = size_group, Parameter = "WLR", Method = "Kruskal-Wallis",
    Mean_Group_1 = NA, Mean_Group_2 = NA, Median_Group_1 = NA, Median_Group_2 = NA,
    Test.Statistic = kruskal_result$statistic, Raw_p_value = kruskal_p,
    Adjusted_p_value = NA, n_Vessels = n_vessels, n_Cats = n_cats,
    Model_Used = NA, Estimate = NA, SE = NA, df = NA, t_ratio = NA
  )

  model_list <- list()

  if (kruskal_p < 0.05) {
    pairwise_comparisons <- list(
      c("No", "Yes_Congestive_HF"),
      c("No", "Yes_ATE"),
      c("Yes_Congestive_HF", "Yes_ATE")
    )

    posthoc_results_list <- list()
    all_p_values <- c()

    for (contrast in pairwise_comparisons) {
      df_pairwise <- df %>% filter(ATE %in% contrast)

      means   <- aggregate(WLR ~ ATE, data = df_pairwise, mean)
      medians <- aggregate(WLR ~ ATE, data = df_pairwise, median)

      df_pairwise$logWLR <- log(df_pairwise$WLR)
      glmm_gauss <- lmer(logWLR ~ ATE + (1 | Cat_number), data = df_pairwise)
      model_list[[paste(contrast, collapse = " vs ")]] <- glmm_gauss

      em_results <- summary(emmeans(glmm_gauss, pairwise ~ ATE, adjust = "none"))$contrasts

      posthoc_results <- as.data.frame(em_results) %>%
        dplyr::rename(Contrast = contrast, Estimate = estimate, SE = SE,
                      df = df, t_ratio = t.ratio, Raw_p_value = p.value) %>%
        mutate(Comparison = size_group, Parameter = "WLR", Method = "Gaussian LMM",
               Mean_Group_1 = means$WLR[1], Mean_Group_2 = means$WLR[2],
               Median_Group_1 = medians$WLR[1], Median_Group_2 = medians$WLR[2],
               n_Vessels = nrow(df_pairwise),
               n_Cats = length(unique(df_pairwise$Cat_number)),
               Model_Used = "Log-Gaussian LMM") %>%
        dplyr::select(Comparison, Parameter, Method, Mean_Group_1, Mean_Group_2,
                      Median_Group_1, Median_Group_2, Contrast,
                      Test.Statistic = t_ratio, Raw_p_value,
                      n_Vessels, n_Cats, Model_Used, Estimate, SE, df)

      all_p_values <- c(all_p_values, posthoc_results$Raw_p_value)
      posthoc_results_list[[paste(contrast, collapse = " vs ")]] <- posthoc_results
    }

    adjusted_p_values <- p.adjust(all_p_values, method = "BH")
    index <- 1
    for (contrast_name in names(posthoc_results_list)) {
      n_rows <- nrow(posthoc_results_list[[contrast_name]])
      posthoc_results_list[[contrast_name]]$Adjusted_p_value <- adjusted_p_values[index:(index + n_rows - 1)]
      index <- index + n_rows
    }

    posthoc_results_df <- bind_rows(posthoc_results_list)
    results <- bind_rows(results, posthoc_results_df)
  }

  return(list(results = results, models = model_list))
}

# Run for each size group
analysis_results <- lapply(c("Small", "Medium", "Large"), function(size) {
  analyze_ATE_groups_vessel(Vessels_results_pheno, size)
})

final_vessel_results <- bind_rows(lapply(analysis_results, function(res) res$results))
final_models <- lapply(analysis_results, function(res) res$models)

print(final_vessel_results)
write.csv(final_vessel_results, file.path(output, "ATE_group_comparisons.csv"), row.names = FALSE)

# --- DHARMa residual diagnostics ---
save_DHARMa_plot <- function(model, comparison_name) {
  simulationOutput <- simulateResiduals(fittedModel = model, plot = FALSE)
  pdf(file.path(output, paste0(comparison_name, "_DHARMa.pdf")), width = 10, height = 6)
  plot(simulationOutput)
  dev.off()
}

for (i in seq_along(final_models)) {
  for (contrast_name in names(final_models[[i]])) {
    save_DHARMa_plot(final_models[[i]][[contrast_name]],
                     paste(c("Small", "Medium", "Large")[i], contrast_name, sep = "_"))
  }
}

# --- Boxplots by vessel size group ---
vessel_boxplot <- function(data, size_group, decimals = 2) {
  data_filtered <- data %>% filter(Size_Group == size_group)
  data_filtered$ATE <- factor(data_filtered$ATE, levels = c("No", "Yes_ATE", "Yes_Congestive_HF"))

  is_large <- size_group == "Large"
  y_limit <- ifelse(is_large, 1.8, NA)
  max_WLR <- max(data_filtered$WLR, na.rm = TRUE)

  p <- ggplot(data_filtered, aes(x = ATE, y = WLR, colour = ATE)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.4) +
    scale_color_manual(values = color_mapping_ate) +
    scale_y_continuous(labels = scales::number_format(accuracy = 10^(-decimals))) +
    labs(title = paste(size_group, "Vessels"), x = "Known cardiac death", y = "Wall-to-lumen ratio") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 14, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 16, face = "bold", color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))

  if (is_large) {
    p <- p +
      coord_cartesian(ylim = c(0, y_limit)) +
      annotate("text", x = 2.9, y = 1.7, label = "Outlier at 3.48", size = 3, color = "black") +
      annotate("segment", x = 2.9, xend = 2.9, y = 1.8, yend = max_WLR,
               arrow = arrow(length = unit(0.2, "cm")))
  }
  return(p)
}

boxplot_list_vessels <- lapply(c("Small", "Medium", "Large"), function(sg) {
  vessel_boxplot(Vessels_results_pheno, sg)
})

pdf(file.path(output, paste0(prefix, "_Vessel_measurements_ATE.pdf")), width = 8, height = 5)
grid.arrange(grobs = boxplot_list_vessels, nrow = 1, ncol = 3)
dev.off()

# --- Multi-parameter Kruskal-Wallis ---
analyze_ATE_groups_kruskal <- function(data, size_group, parameters) {
  df <- data %>% filter(Size_Group == size_group)
  n_vessels <- nrow(df)
  n_cats <- length(unique(df$Cat_number))

  results_list <- list()
  for (param in parameters) {
    shapiro_p <- shapiro.test(pull(df, param))$p.value
    kruskal_result <- kruskal.test(pull(df, param) ~ df$ATE)

    means   <- df %>% group_by(ATE) %>% summarize(mean_value = mean(!!sym(param), na.rm = TRUE), .groups = "drop")
    medians <- df %>% group_by(ATE) %>% summarize(median_value = median(!!sym(param), na.rm = TRUE), .groups = "drop")

    results_list[[param]] <- data.frame(
      Comparison = size_group, Parameter = param, Method = "Kruskal-Wallis",
      Mean_HCM_only = means$mean_value[1], Mean_HCM_with_ATE = means$mean_value[2],
      Mean_No_ATE = means$mean_value[3],
      Median_HCM_only = medians$median_value[1], Median_HCM_with_ATE = medians$median_value[2],
      Median_No_ATE = medians$median_value[3],
      Normality_p_value = shapiro_p,
      Test.Statistic = kruskal_result$statistic, Raw_p_value = kruskal_result$p.value,
      n_Vessels = n_vessels, n_Cats = n_cats
    )
  }
  bind_rows(results_list)
}

vessel_params <- c("OD", "LD", "WT", "WLR", "WCSA")
final_kruskal <- bind_rows(
  analyze_ATE_groups_kruskal(Vessels_results_pheno, "Small", vessel_params),
  analyze_ATE_groups_kruskal(Vessels_results_pheno, "Medium", vessel_params),
  analyze_ATE_groups_kruskal(Vessels_results_pheno, "Large", vessel_params)
)
print(final_kruskal)
write.csv(final_kruskal,
          file = paste0(output, prefix, "_ATE_group_comparisons_Kruskal.csv"),
          row.names = FALSE)


# ==============================================================================
# SECTION 8: NUCLEI ANALYSIS
# ==============================================================================
# Analyzes cardiomyocyte nucleus shape from H&E staining using
# mixed-effects models. Includes small hematoxylin-dense nuclei analysis.

# --- Load data ---
nuclei <- read.csv(file.path(qupath_path, "Nuclei_HE_300725.csv"))

# --- Prepare data ---
nuclei <- nuclei %>% mutate(image_id = sub("\\.ndpi", "", Image))
nuclei_clean <- nuclei %>%
  left_join(metadata_HE %>% dplyr::select(image_id, Cat_number, Tissue_part), by = "image_id")
nuclei_clean <- nuclei_clean[nuclei_clean$Cat_number %in% INCLUDE_CATS, ]

# Add clinical data for ATE and Known_cardiac_death
nuclei_clean <- nuclei_clean %>%
  left_join(
    clinical_data %>% dplyr::select(Cat_number, ATE, Known_cardiac_death),
    by = "Cat_number"
  )

# Convert to data.table, filter, clean
nuclei_dt <- as.data.table(nuclei_clean)
nuclei_dt <- nuclei_dt[Tissue_part %in% c("L1", "L2", "WH")]
nuclei_dt <- nuclei_dt[complete.cases(nuclei_dt$Length.µm)]
nuclei_dt[, Known_cardiac_death := factor(Known_cardiac_death, levels = c("No", "Yes"))]
nuclei_dt[, ATE := factor(ATE)]

# --- Additional shape metrics ---
nuclei_dt <- nuclei_dt %>%
  mutate(
    Aspect_Ratio = Max.diameter.µm / Min.diameter.µm,
    Roundness = 4 * Area.µm.2 / (pi * Max.diameter.µm^2),
    Form_Factor = Area.µm.2 / Perimeter.µm,
    Eccentricity = sqrt(1 - (Min.diameter.µm^2 / Max.diameter.µm^2))
  )

# --- Split by orientation ---
transversal  <- nuclei_dt[nuclei_dt$Parent == "Annotation (Transversal)", ]
longitudinal <- nuclei_dt[nuclei_dt$Parent == "Annotation (Longitudinal)", ]

# --- Model selection and fitting functions ---
fit_best_model <- function(data, parameter, direction, output, prefix) {
  print(paste("Fitting models for:", parameter, "in", direction, "data"))

  model_gaussian <- lmer(as.formula(paste(parameter, "~ ATE + (1|Cat_number)")), data = data)
  model_log <- lmer(as.formula(paste0("log(", parameter, " + 1) ~ ATE + (1|Cat_number)")), data = data)
  model_gls <- gls(as.formula(paste(parameter, "~ ATE")), data = data, weights = varExp())
  model_gamma <- tryCatch(
    glmer(as.formula(paste(parameter, "~ ATE + (1|Cat_number)")), family = Gamma(link = "log"), data = data),
    error = function(e) NULL
  )

  aic_values <- data.frame(
    Model = c("Gaussian", "Log", "GLS", "Gamma"),
    AIC = c(AIC(model_gaussian), AIC(model_log), AIC(model_gls),
            if (!is.null(model_gamma)) AIC(model_gamma) else NA)
  )

  best_model <- aic_values$Model[which.min(aic_values$AIC)]
  best_model_object <- switch(best_model,
                              "Gaussian" = model_gaussian, "Log" = model_log,
                              "GLS" = model_gls, "Gamma" = model_gamma)

  sim_best <- simulateResiduals(fittedModel = best_model_object, plot = FALSE)
  plot_file <- file.path(output, paste0(prefix, "_", direction, "_", parameter, "_", best_model, "_DHARMa.pdf"))
  pdf(plot_file, width = 10, height = 6)
  plot(sim_best, main = paste("DHARMa Plot - Best Model:", best_model, "for", parameter))
  dev.off()

  print(aic_values)
  return(best_model)
}

run_final_model <- function(data, parameter, direction, model_choice, output, prefix) {
  print(paste("Running final model for:", parameter, "in", direction, "data using", model_choice))

  final_model <- switch(model_choice,
    "Gaussian" = lmer(as.formula(paste(parameter, "~ ATE + (1|Cat_number)")), data = data),
    "Log"      = lmer(as.formula(paste0("log(", parameter, " + 1) ~ ATE + (1|Cat_number)")), data = data),
    "GLS"      = gls(as.formula(paste(parameter, "~ ATE")), data = data, weights = varExp()),
    "Gamma"    = glmer(as.formula(paste(parameter, "~ ATE + (1|Cat_number)")), family = Gamma(link = "log"), data = data)
  )

  group_stats <- data %>%
    group_by(ATE) %>%
    summarise(Mean = mean(get(parameter), na.rm = TRUE),
              Median = median(get(parameter), na.rm = TRUE),
              n_nuclei = n(), n_cats = length(unique(Cat_number))) %>%
    arrange(ATE)

  test_statistic <- NA
  model_p_value <- NA
  model_aic <- AIC(final_model)

  if (model_choice %in% c("Gaussian", "Log")) {
    anova_results <- tryCatch(anova(final_model, type = 3), error = function(e) NULL)
    if (!is.null(anova_results) && "Pr(>F)" %in% colnames(anova_results)) {
      test_statistic <- anova_results$`F value`[1]
      model_p_value <- anova_results$`Pr(>F)`[1]
    }
  } else if (model_choice == "Gamma") {
    anova_results <- tryCatch(car::Anova(final_model, type = 3), error = function(e) NULL)
    if (!is.null(anova_results) && "Pr(>Chisq)" %in% colnames(anova_results)) {
      test_statistic <- anova_results$Chisq[1]
      model_p_value <- anova_results$`Pr(>Chisq)`[1]
    }
  } else if (model_choice == "GLS") {
    anova_results <- tryCatch(anova(final_model), error = function(e) NULL)
    if (!is.null(anova_results) && "p-value" %in% colnames(anova_results)) {
      test_statistic <- anova_results$`F value`[1]
      model_p_value <- anova_results$`p-value`[1]
    }
  }

  model_summary <- data.frame(
    Direction = direction, contrast = "Model Summary", Parameter = parameter,
    Method = model_choice, Mean_Group_1 = NA, Mean_Group_2 = NA,
    Median_Group_1 = NA, Median_Group_2 = NA,
    Test_Statistic = test_statistic, Raw_p_value = model_p_value,
    Adjusted_p_value = NA, n_nuclei_Group_1 = NA, n_nuclei_Group_2 = NA,
    n_cats_Group_1 = NA, n_cats_Group_2 = NA,
    estimate = NA, SE = NA, df = NA, AIC = model_aic
  )

  posthoc <- emmeans(final_model, pairwise ~ ATE, adjust = "BH")
  posthoc_results <- as.data.frame(summary(posthoc$contrasts))

  posthoc_results <- posthoc_results %>%
    mutate(Adjusted_p_value = p.adjust(p.value, method = "BH")) %>%
    rowwise() %>%
    mutate(
      Group_1 = strsplit(contrast, " - ")[[1]][1],
      Group_2 = strsplit(contrast, " - ")[[1]][2],
      Mean_Group_1 = group_stats$Mean[group_stats$ATE == Group_1],
      Mean_Group_2 = group_stats$Mean[group_stats$ATE == Group_2],
      Median_Group_1 = group_stats$Median[group_stats$ATE == Group_1],
      Median_Group_2 = group_stats$Median[group_stats$ATE == Group_2],
      n_nuclei_Group_1 = group_stats$n_nuclei[group_stats$ATE == Group_1],
      n_nuclei_Group_2 = group_stats$n_nuclei[group_stats$ATE == Group_2],
      n_cats_Group_1 = group_stats$n_cats[group_stats$ATE == Group_1],
      n_cats_Group_2 = group_stats$n_cats[group_stats$ATE == Group_2]
    ) %>%
    ungroup()

  posthoc_results <- posthoc_results %>%
    mutate(Direction = direction, Parameter = parameter, Method = model_choice,
           Test_Statistic = NA, AIC = NA) %>%
    dplyr::select(Direction, contrast, Parameter, Method,
                  Mean_Group_1, Mean_Group_2, Median_Group_1, Median_Group_2,
                  Test_Statistic, p.value, Adjusted_p_value,
                  n_nuclei_Group_1, n_nuclei_Group_2,
                  n_cats_Group_1, n_cats_Group_2,
                  estimate, SE, df, AIC) %>%
    dplyr::rename(Raw_p_value = p.value)

  results_table <- bind_rows(model_summary, posthoc_results)

  file_path <- paste0(output, prefix, "_final_model_results_", direction, "_", parameter, ".csv")
  write.csv(results_table, file = file_path, row.names = FALSE)
  return(results_table)
}

# --- Run models for all parameters and orientations ---
nuclei_parameters <- c("Area.µm.2", "Solidity", "Eccentricity")
datasets <- list(transversal = transversal, longitudinal = longitudinal)
all_nuclei_results <- list()

for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  for (param in nuclei_parameters) {
    message("Running for: ", dataset_name, " | Parameter: ", param)
    selected_model <- fit_best_model(dataset, param, dataset_name, output = output, prefix = prefix)
    result <- run_final_model(dataset, param, dataset_name, selected_model, output = output, prefix = prefix)
    all_nuclei_results[[paste(dataset_name, param, sep = "_")]] <- list(model = selected_model, result = result)
  }
}

# Combine all results
nuclei_results_df <- purrr::map_dfr(names(all_nuclei_results), function(key) {
  entry <- all_nuclei_results[[key]]
  parts <- strsplit(key, "_")[[1]]
  data.frame(Dataset = parts[1], Parameter = parts[2], Model = entry$model,
             entry$result, stringsAsFactors = FALSE)
})
write.csv(nuclei_results_df,
          file = file.path(output, paste0(prefix, "_Final_Model_Results.csv")),
          row.names = FALSE)

# --- Boxplots: longitudinal nuclei ---
agg_data_long <- longitudinal[, lapply(.SD, mean, na.rm = TRUE),
                              by = .(Cat_number, ATE), .SDcols = nuclei_parameters]
agg_data_long$ATE <- as.factor(agg_data_long$ATE)
agg_data_long[, (nuclei_parameters) := lapply(.SD, as.numeric), .SDcols = nuclei_parameters]

dotplot_list_long <- list()
for (param in nuclei_parameters) {
  param_sym <- sym(param)
  p <- ggplot(agg_data_long, aes(x = ATE, y = !!param_sym, color = ATE)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
    scale_color_manual(values = color_mapping_ate) +
    labs(x = "Known cardiac death", y = param) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 14, color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))
  dotplot_list_long[[param]] <- p
}

pdf(file.path(output, paste0(prefix, "_longitudinal_Nucleus_Measurements.pdf")), width = 13, height = 4.5)
grid.arrange(grobs = dotplot_list_long, nrow = 1, ncol = 3)
dev.off()

# --- Boxplots: all nuclei (transversal + longitudinal) ---
agg_data_all <- nuclei_dt[, lapply(.SD, mean, na.rm = TRUE),
                          by = .(Cat_number, ATE), .SDcols = nuclei_parameters]
agg_data_all$Group <- factor(agg_data_all$ATE, levels = c("No", "Yes_ATE", "Yes_Congestive_HF"))

dotplot_list_all <- list()
for (param in nuclei_parameters) {
  param_sym <- sym(param)
  p <- ggplot(agg_data_all, aes(x = Group, y = !!param_sym, color = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
    scale_color_manual(values = color_mapping_ate) +
    labs(x = "Known cardiac death", y = param) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(color = "black"))
  dotplot_list_all[[param]] <- p
}

pdf(file.path(output, paste0(prefix, "_Transversal_Nucleus_Measurements.pdf")), width = 13, height = 4)
grid.arrange(grobs = dotplot_list_all, nrow = 1, ncol = 3)
dev.off()

# Numbers for paper
nuclei_summary <- nuclei_dt %>%
  group_by(Cat_number) %>%
  summarise(n_nuclei = n()) %>%
  summarise(Mean_Nuclei_Per_Cat = mean(n_nuclei), SD_Nuclei_Per_Cat = sd(n_nuclei))
print(nuclei_summary)

# --- Small hematoxylin-dense nuclei ---
longitudinal <- longitudinal %>%
  group_by(Cat_number) %>%
  mutate(Hematoxylin_cutoff = quantile(Hematoxylin..Median, 0.5, na.rm = TRUE)) %>%
  ungroup()

candidates <- longitudinal %>%
  filter(Area.µm.2 > 0, Area.µm.2 < 35, Hematoxylin..Median > Hematoxylin_cutoff)

candidate_counts <- candidates %>%
  group_by(Cat_number) %>%
  summarise(N_small_round_dark = n(), .groups = "drop")

total_counts <- longitudinal %>%
  group_by(Cat_number, ATE) %>%
  summarise(Total_nuclei = n(), .groups = "drop")

density_summary <- total_counts %>%
  left_join(candidate_counts, by = "Cat_number") %>%
  mutate(N_small_round_dark = tidyr::replace_na(N_small_round_dark, 0),
         Percentage = N_small_round_dark / Total_nuclei * 100)

density_summary$ATE <- factor(density_summary$ATE, levels = c("No", "Yes_ATE", "Yes_Congestive_HF"))

# Statistical test
normality_pvals <- density_summary %>%
  group_by(ATE) %>%
  summarise(p_value = shapiro.test(N_small_round_dark)$p.value)
print(normality_pvals)

all_normal <- all(normality_pvals$p_value > 0.05)
if (all_normal) {
  levene_p <- leveneTest(N_small_round_dark ~ ATE, data = density_summary)$`Pr(>F)`[1]
  homoscedastic <- levene_p > 0.05
} else {
  homoscedastic <- FALSE
}

if (all_normal && homoscedastic) {
  test_result_nuclei <- aov(N_small_round_dark ~ ATE, data = density_summary)
  print(summary(test_result_nuclei))
  cat("\nTest used: ANOVA\n")
  print(TukeyHSD(test_result_nuclei))
} else {
  test_result_nuclei <- kruskal.test(N_small_round_dark ~ ATE, data = density_summary)
  print(test_result_nuclei)
  cat("\nTest used: Kruskal-Wallis\n")
  print(dunnTest(N_small_round_dark ~ ATE, data = density_summary, method = "bh"))
}

# Plot
p_nuclei <- ggplot(density_summary, aes(x = ATE, y = N_small_round_dark, color = ATE)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  scale_color_manual(values = color_mapping_ate) +
  theme_minimal() +
  labs(x = "Known cardiac death",
       y = "Small, hematoxylin-dense nuclei",
       title = "Nucleus count by group") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 16, face = "bold", color = "black"),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1))

pdf(file.path(output, paste0(prefix, "_Small_hematoxylin_dense_nuclei.pdf")), width = 4.5, height = 5)
print(p_nuclei)
dev.off()

# Pairwise stats for small nuclei
analyze_pairwise(
  data = density_summary,
  params = c("N_small_round_dark"),
  output = output,
  prefix = prefix,
  filename = "SmallHematoxylinDenseNuclei_Pairwise_Stats.csv"
)

# --- Export nucleus measurements per cat ---
nuclei_dt[, Orientation := ifelse(grepl("Longitudinal", Parent), "Longitudinal",
                                  ifelse(grepl("Transversal", Parent), "Transversal", NA))]
nuclei_dt$Orientation <- as.factor(nuclei_dt$Orientation)

export_params <- c("Area.µm.2", "Circularity", "Solidity", "Max.diameter.µm",
                   "Min.diameter.µm", "Perimeter.µm", "Eccentricity")

average_per_cat_orientation <- nuclei_dt[, lapply(.SD, mean, na.rm = TRUE),
                                         by = .(Cat_number, Orientation),
                                         .SDcols = export_params]

average_per_cat_orientation <- dcast(average_per_cat_orientation,
                                     Cat_number ~ Orientation,
                                     value.var = export_params, fill = NA)

colnames(average_per_cat_orientation) <- c("Cat_number",
                                           paste0("Longitudinal_", export_params),
                                           paste0("Transversal_", export_params))

average_per_cat_orientation <- merge(
  average_per_cat_orientation,
  density_summary[, c("Cat_number", "N_small_round_dark")],
  by = "Cat_number", all.x = TRUE
)
names(average_per_cat_orientation)[names(average_per_cat_orientation) == "N_small_round_dark"] <-
  "N_small_hematoxylin_dense_nuclei"

write.csv(average_per_cat_orientation,
          file = file.path(output, paste0(prefix, "_Nucleus_Measurements.csv")),
          row.names = FALSE)


# ==============================================================================
# SECTION 9: PHENOTYPE FILE ASSEMBLY
# ==============================================================================
# Merges per-cat summary data from all measurement types into a single
# phenotype file used for correlation and GSEA analyses.

Adipo_perc_assembled <- read.csv(paste0(output, prefix, "_Adipocyte_percentage_per_cat.csv"))
Adipo_size_assembled <- read.csv(paste0(output, prefix, "_Adipocyte_size_per_cat.csv"))
Fibrosis_assembled   <- read.csv(paste0(output, prefix, "_Fibrosis_percentage_per_cat.csv"))
Fibrosis_assembled$Cat_number <- as.integer(Fibrosis_assembled$Cat_number)
Vessels_assembled    <- read.csv(paste0(output, prefix, "_Vessel_Calculations.csv"))
# Note: use Vessel_measurements_mean computed above for per-cat means
Clinical_assembled   <- read.delim(file.path(metadata_path, "260225_Clinical_data_table_1.txt"))
Nuclei_assembled     <- read.csv(paste0(output, prefix, "_Nucleus_Measurements.csv"))

data_list <- list(Adipo_perc_assembled, Adipo_size_assembled, Fibrosis_assembled,
                  Vessel_measurements_mean, Clinical_assembled, Nuclei_assembled)

merged_data <- Reduce(function(x, y) merge(x, y, by = "Cat_number", all = TRUE), data_list)

# Resolve duplicate columns (.x/.y pairs)
cols_x <- grep("\\.x$", names(merged_data), value = TRUE)
cols_y <- gsub("\\.x$", ".y", cols_x)

for (i in seq_along(cols_x)) {
  if (cols_y[i] %in% names(merged_data)) {
    na_x <- sum(is.na(merged_data[[cols_x[i]]]))
    na_y <- sum(is.na(merged_data[[cols_y[i]]]))
    if (na_x > na_y) {
      merged_data <- merged_data %>% dplyr::select(-all_of(cols_x[i]))
      names(merged_data)[names(merged_data) == cols_y[i]] <- gsub("\\.y$", "", cols_y[i])
    } else {
      merged_data <- merged_data %>% dplyr::select(-all_of(cols_y[i]))
      names(merged_data)[names(merged_data) == cols_x[i]] <- gsub("\\.x$", "", cols_x[i])
    }
  }
}

# Remove cat 8
merged_data <- merged_data[merged_data$Cat_number %in% INCLUDE_CATS, ]
write.csv(merged_data, file = paste0(output, prefix, "_Phenotype_file.csv"))

# Age calculations per group
merged_data$Age_years <- as.numeric(merged_data$Age_years)
merged_data$ATE <- standardize_ATE(merged_data$ATE)
cat("Mean age Yes_Congestive_HF:",
    mean(merged_data$Age_years[merged_data$ATE == "Yes_Congestive_HF"],
         na.rm = TRUE), "\n")
cat("Mean age Yes_ATE:",
    mean(merged_data$Age_years[merged_data$ATE == "Yes_ATE"],
         na.rm = TRUE), "\n")
cat("Mean age No:",
    mean(merged_data$Age_years[merged_data$ATE == "No"],
         na.rm = TRUE), "\n")


# ==============================================================================
# SECTION 10: PHENOTYPE CORRELATIONS
# ==============================================================================
# Computes pairwise correlations between all phenotype variables using
# Pearson (if normally distributed) or Spearman (otherwise), with
# Benjamini-Hochberg correction and bootstrap confidence intervals.

# --- Load phenotype data ---
pheno_corr <- read.csv(paste0(output, prefix, "_Phenotype_file.csv"))

pheno_corr[pheno_corr == "Unknown"] <- NA
pheno_corr$Age_years <- as.numeric(pheno_corr$Age_years)
pheno_corr$ATE <- standardize_ATE(pheno_corr$ATE)
pheno_corr$Known_cardiac_death <- ifelse(
  standardize_KCD(pheno_corr$Known_cardiac_death) == "Yes", 1, 0)
pheno_corr$HFwithATE <- ifelse(pheno_corr$ATE == "Yes_ATE", 1, 0)
pheno_corr$HFnoATE   <- ifelse(pheno_corr$ATE == "Yes_Congestive_HF", 1, 0)
pheno_corr$X <- NULL
pheno_corr$Female <- NULL

pheno_corr <- data.frame(lapply(pheno_corr, function(x) {
  x[is.infinite(x)] <- NA
  return(x)
}))

# All columns for correlation statistics (CSV)
selected_columns <- c("Age_years", "Male", "Body_Condition", "Known_cardiac_death",
                       "Mean_NonMyocardial_Fibrosis_Percentage", "Mean_Adipocyte_perc",
                       "Mean_Adipocyte_size", "Mean_Myocardial_Fibrosis_Percentage",
                       "Mean_WLR", "Longitudinal_Eccentricity", "Longitudinal_Area.µm.2",
                       "Longitudinal_Solidity", "Transversal_Eccentricity",
                       "Transversal_Area.µm.2", "Transversal_Solidity",
                       "HFwithATE", "HFnoATE", "N_small_hematoxylin_dense_nuclei")

# Columns for heatmap (excludes adipocyte variables)
heatmap_columns <- setdiff(selected_columns, c("Mean_Adipocyte_perc", "Mean_Adipocyte_size"))
data_corr <- pheno_corr[, selected_columns]
colnames(data_corr) <- gsub("^Mean_", "", colnames(data_corr))
data_corr <- data.frame(lapply(data_corr, as.numeric))

cat("Missing values in selected columns:\n")
print(colSums(is.na(data_corr)))

# Normality tests
normality_results <- apply(data_corr, 2, function(col)
  ifelse(sum(!is.na(col)) > 3, shapiro.test(col)$p.value, NA))
print(normality_results)

# --- Correlation functions ---
cor_matrix_mixed <- function(data, normality_results, conf.level = 0.95) {
  n <- ncol(data)
  cor_matrix <- matrix(NA, n, n, dimnames = list(colnames(data), colnames(data)))
  p_matrix   <- matrix(NA, n, n, dimnames = list(colnames(data), colnames(data)))

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (sum(!is.na(data[, i])) > 3 && sum(!is.na(data[, j])) > 3) {
        method <- if (normality_results[i] > 0.05 && normality_results[j] > 0.05) "pearson" else "spearman"
        test <- cor.test(data[, i], data[, j], use = "pairwise.complete.obs",
                         method = method, conf.level = conf.level)
        cor_matrix[i, j] <- cor_matrix[j, i] <- test$estimate
        p_matrix[i, j]   <- p_matrix[j, i]   <- test$p.value
      }
    }
  }
  diag(cor_matrix) <- NA
  diag(p_matrix) <- NA
  list(cor_matrix = cor_matrix, p_matrix = p_matrix)
}

bootstrap_ci <- function(x, y, method = "spearman", conf.level = 0.95, n_boot = 1000) {
  cor_boot <- replicate(n_boot, {
    idx <- sample(seq_along(x), replace = TRUE)
    if (sum(!is.na(x[idx]) & !is.na(y[idx])) > 3) cor(x[idx], y[idx], method = method) else NA
  })
  cor_boot <- cor_boot[!is.na(cor_boot)]
  if (length(cor_boot) > 0) {
    quantile(cor_boot, c((1 - conf.level) / 2, 1 - (1 - conf.level) / 2), na.rm = TRUE)
  } else {
    c(NA, NA)
  }
}

ggcorrplot_lower_triangle <- function(cor_matrix, ci_matrix, p_matrix, sig.level = 0.05) {
  vars <- colnames(cor_matrix)
  n <- length(vars)

  # Build data frame from lower triangle indices directly
  idx <- which(lower.tri(cor_matrix), arr.ind = TRUE)
  corr_data <- data.frame(
    Var1        = factor(vars[idx[, 1]], levels = vars),
    Var2        = factor(vars[idx[, 2]], levels = vars),
    value       = cor_matrix[idx],
    p_value     = p_matrix[idx],
    ci_label    = ci_matrix[idx],
    stringsAsFactors = FALSE
  )
  corr_data$significant <- ifelse(corr_data$p_value < sig.level, "bold", "plain")
  corr_data$cor_label   <- sprintf("%.2f", corr_data$value)

  ggplot(corr_data, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = cor_label, fontface = significant), size = 3.5, na.rm = TRUE) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10),
          axis.text.y = element_text(vjust = 1, hjust = 1, size = 10),
          axis.title = element_blank(), panel.grid = element_blank(),
          axis.line.y = element_line(color = "black"),
          axis.ticks.y = element_blank(),
          axis.text.y.left = element_text(hjust = 1, margin = margin(r = -10))) +
    scale_y_discrete(position = "right")
}

extract_statistics <- function(cor_matrix, raw_p_matrix, adjusted_p_matrix, ci_matrix, methods_matrix) {
  n <- ncol(cor_matrix)
  stats_list <- list()
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (!is.na(cor_matrix[i, j])) {
        stats_list[[length(stats_list) + 1]] <- data.frame(
          Variable1 = colnames(cor_matrix)[i],
          Variable2 = colnames(cor_matrix)[j],
          Correlation = cor_matrix[i, j],
          Raw_P_Value = raw_p_matrix[i, j],
          Adjusted_P_Value = adjusted_p_matrix[i, j],
          Test_Method = methods_matrix[i, j],
          Confidence_Interval = ci_matrix[i, j]
        )
      }
    }
  }
  do.call(rbind, stats_list)
}

# --- Compute correlations ---
mixed_results <- cor_matrix_mixed(data_corr, normality_results)
cor_matrix <- mixed_results$cor_matrix
p_matrix   <- mixed_results$p_matrix

# BH correction
p_matrix_corrected <- p.adjust(p_matrix[lower.tri(p_matrix)], method = "BH")
p_matrix[lower.tri(p_matrix)] <- p_matrix_corrected
p_matrix[upper.tri(p_matrix)] <- t(p_matrix)[upper.tri(p_matrix)]

# Bootstrap CIs
ci_matrix <- matrix(NA, ncol(data_corr), ncol(data_corr),
                    dimnames = list(colnames(data_corr), colnames(data_corr)))

for (i in 1:(ncol(data_corr) - 1)) {
  for (j in (i + 1):ncol(data_corr)) {
    if (sum(!is.na(data_corr[, i])) > 3 && sum(!is.na(data_corr[, j])) > 3) {
      method <- if (normality_results[i] > 0.05 && normality_results[j] > 0.05) "pearson" else "spearman"
      ci <- bootstrap_ci(data_corr[, i], data_corr[, j], method = method)
      ci_matrix[i, j] <- ci_matrix[j, i] <- sprintf("[%.2f, %.2f]", ci[1], ci[2])
    }
  }
}

# Heatmap (without adipocyte variables)
heatmap_vars <- gsub("^Mean_", "", heatmap_columns)
cor_matrix_hm <- cor_matrix[heatmap_vars, heatmap_vars]
ci_matrix_hm  <- ci_matrix[heatmap_vars, heatmap_vars]
p_matrix_hm   <- p_matrix[heatmap_vars, heatmap_vars]

pdf(file.path(output, paste0(prefix, "_Phenotype_Correlations_LowerTriangle.pdf")),
    width = 11, height = 8)
ggcorrplot_lower_triangle(cor_matrix_hm, ci_matrix_hm, p_matrix_hm)
dev.off()

# --- Export statistics ---
methods_matrix <- matrix(NA, ncol(data_corr), ncol(data_corr),
                         dimnames = list(colnames(data_corr), colnames(data_corr)))

mixed_results2 <- cor_matrix_mixed(data_corr, normality_results)
cor_matrix2    <- mixed_results2$cor_matrix
raw_p_matrix   <- mixed_results2$p_matrix

for (i in 1:(ncol(data_corr) - 1)) {
  for (j in (i + 1):ncol(data_corr)) {
    if (!is.na(cor_matrix2[i, j])) {
      methods_matrix[i, j] <- methods_matrix[j, i] <-
        ifelse(normality_results[i] > 0.05 && normality_results[j] > 0.05, "Pearson", "Spearman")
    }
  }
}

adjusted_p_matrix2 <- raw_p_matrix
p_corrected2 <- p.adjust(raw_p_matrix[lower.tri(raw_p_matrix)], method = "BH")
adjusted_p_matrix2[lower.tri(adjusted_p_matrix2)] <- p_corrected2
adjusted_p_matrix2[upper.tri(adjusted_p_matrix2)] <- t(adjusted_p_matrix2)[upper.tri(adjusted_p_matrix2)]

correlation_statistics <- extract_statistics(cor_matrix2, raw_p_matrix, adjusted_p_matrix2,
                                             ci_matrix, methods_matrix)
write.csv(correlation_statistics,
          file = file.path(output, paste0(prefix, "_Phenotype_Correlations_with_CI_and_Methods.csv")),
          row.names = FALSE)


# ==============================================================================
# SECTION 11: NANOPORE CORRELATIONS & GSEA
# ==============================================================================
# Correlates Nanopore gene expression (featureCounts) with phenotype variables
# using Spearman correlation, then runs GO Biological Process GSEA.

# --- Load and align count data ---
nanopore_files <- list.files(nanopore_path, pattern = "\\.tabular$", full.names = TRUE)

process_file <- function(file) {
  data <- read.delim(file, header = TRUE)
  catID <- sub("catID_", "", regmatches(file, regexpr("catID_\\d+", file)))
  colnames(data)[2] <- catID
  setNames(data[2], catID)
}

count_data <- do.call(cbind, lapply(nanopore_files, process_file))
rownames(count_data) <- read.delim(nanopore_files[1], header = TRUE)[[1]]
write.csv(count_data, file = paste0(output, prefix, "_count_data.csv"))

# Align with phenotype data
nanopore_phenotype <- read.csv(paste0(output, prefix, "_Phenotype_file.csv"))
nanopore_phenotype$Cat_number <- as.character(nanopore_phenotype$Cat_number)

matching_rows <- nanopore_phenotype %>%
  filter(Cat_number %in% colnames(count_data))

count_data_aligned <- count_data[, matching_rows$Cat_number]
write.csv(matching_rows, file = paste0(output, prefix, "_Phenotype_data_nanopore_cats.csv"))

# --- Correlation + GSEA loop ---
matching_rows$Age_years <- as.numeric(matching_rows$Age_years)

phenotype_columns <- c("Mean_Adipocyte_perc", "Mean_Myocardial_Fibrosis_Percentage",
                       "Mean_WLR", "Longitudinal_Area.µm.2", "Longitudinal_Solidity",
                       "Longitudinal_Eccentricity", "Transversal_Area.µm.2",
                       "Transversal_Solidity", "Transversal_Eccentricity",
                       "Age_years", "N_small_hematoxylin_dense_nuclei")

significant_gene_list <- list()
top_pathways_list <- list()

for (phenotype_column in phenotype_columns) {
  cat("\n--- Running analysis for:", phenotype_column, "---\n")

  # Spearman correlation per gene
  cor_results <- apply(count_data_aligned, 1, function(gene_counts) {
    cor_test <- cor.test(gene_counts, matching_rows[[phenotype_column]],
                         use = "complete.obs", method = "spearman")
    c(correlation = cor_test$estimate, pvalue = cor_test$p.value)
  })

  cor_results_df <- as.data.frame(t(cor_results))
  colnames(cor_results_df) <- c("Correlation", "Pvalue")
  cor_results_df$Gene <- rownames(count_data_aligned)
  cor_results_df <- cor_results_df[, c("Gene", "Correlation", "Pvalue")]
  cor_results_df <- cor_results_df %>% arrange(desc(abs(Correlation)))

  write.csv(cor_results_df,
            file = paste0(output, prefix, "_Correlation_", phenotype_column, ".csv"),
            row.names = FALSE)

  # Significant genes (p < 0.05)
  significant_correlations <- cor_results_df %>% filter(Pvalue < 0.05)
  significant_correlations$Correlation <- as.vector(significant_correlations$Correlation)

  if (nrow(significant_correlations) > 0) {
    significant_correlations$Phenotype <- phenotype_column
    significant_gene_list[[phenotype_column]] <- significant_correlations
  }

  # GSEA
  gene_list <- cor_results_df$Correlation
  names(gene_list) <- cor_results_df$Gene
  gene_list <- sort(gene_list, decreasing = TRUE)

  set.seed(123)
  gsea_results <- gseGO(
    geneList = gene_list, ont = "BP", OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL", minGSSize = 10, maxGSSize = 500,
    pvalueCutoff = 0.05, verbose = TRUE
  )

  gsea_df <- as.data.frame(gsea_results)
  write.csv(gsea_df,
            file = paste0(output, prefix, "_GSEA_Association_with_", phenotype_column, ".csv"),
            row.names = FALSE)

  # Simplify
  simplified_gsea <- clusterProfiler::simplify(gsea_results, cutoff = 0.5, by = "p.adjust",
                                               select_fun = min, measure = "Wang")
  simplified_df <- as.data.frame(simplified_gsea)
  write.csv(simplified_df,
            file = paste0(output, prefix, "_Simplified_GSEA_Association_with_", phenotype_column, ".csv"),
            row.names = FALSE)

  # Top 10 pathways
  top10 <- simplified_df %>%
    arrange(p.adjust) %>%
    head(10) %>%
    mutate(Phenotype = phenotype_column) %>%
    dplyr::select(Phenotype, ID, Description, setSize, enrichmentScore,
                  NES, pvalue, p.adjust, qvalue, rank, leading_edge, core_enrichment)

  top_pathways_list[[phenotype_column]] <- top10
}

# Combined reports
final_top_pathways <- bind_rows(top_pathways_list)
write.csv(final_top_pathways,
          file = paste0(output, prefix, "_Top10_Pathways_per_Phenotype.csv"),
          row.names = FALSE)

final_significant_genes <- bind_rows(significant_gene_list)
write.csv(final_significant_genes,
          file = paste0(output, prefix, "_All_Significant_Genes_Pvalue_lt_0.05.csv"),
          row.names = FALSE)

# --- Top 3 GO BP pathways bubble plot ---
selected_phenotypes <- c("Mean_Myocardial_Fibrosis_Percentage", "Mean_WLR",
                         "Longitudinal_Area.µm.2", "Longitudinal_Solidity",
                         "Longitudinal_Eccentricity", "Transversal_Area.µm.2",
                         "Transversal_Solidity", "Transversal_Eccentricity",
                         "Age_years", "N_small_hematoxylin_dense_nuclei")

top_pathways_plot <- final_top_pathways %>%
  filter(!is.na(NES), !is.na(p.adjust), !is.na(Description)) %>%
  filter(Phenotype %in% selected_phenotypes) %>%
  group_by(Phenotype) %>%
  arrange(p.adjust) %>%
  slice_head(n = 3) %>%
  ungroup()

top_pathways_plot$Description_short <- str_trunc(top_pathways_plot$Description, 80)

p_gsea <- ggplot(top_pathways_plot,
                 aes(x = Phenotype,
                     y = fct_reorder(Description_short, NES),
                     size = -log10(p.adjust),
                     color = NES)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "Top 3 enriched GO pathways per phenotype parameter",
       x = "Phenotype", y = "GO Pathway",
       color = "NES", size = "-log10(p.adjust)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

pdf(file.path(output, paste0(prefix, "_Top3_GO_BP.pdf")), width = 12, height = 10)
print(p_gsea)
dev.off()

# ==============================================================================
# SECTION 12: COMBINED EXCEL EXPORT
# ==============================================================================
# Creates a single Excel file with 3 sheets:
#   Sheet 1: All statistical results combined
#   Sheet 2: Phenotype correlations
#   Sheet 3: Top 10 GSEA pathways per phenotype

# --- Sheet 1: Collect all statistics CSVs ---
stats_files <- list(
  Fibrosis_CardiacDeath = paste0(output, prefix,
    "_All_Fibrosis_Perc_Results.csv"),
  Fibrosis_ATE_Myocardial = paste0(output, prefix,
    "_Mean_Myocardial_Fibrosis_Percentage",
    "_Fibrosis_Pairwise_Stats.csv"),
  Fibrosis_ATE_Total = paste0(output, prefix,
    "_Mean_Fibrosis_Percentage",
    "_Fibrosis_Pairwise_Stats.csv"),
  Fibrosis_ATE_NonMyocardial = paste0(output, prefix,
    "_Mean_NonMyocardial_Fibrosis_Percentage",
    "_Fibrosis_Pairwise_Stats.csv"),
  Adipocyte_Perc_ATE = paste0(output, prefix,
    "_Pairwise_AdipocytePerc_Results.csv"),
  Adipocyte_Size_ATE = paste0(output, prefix,
    "_Pairwise_AdipocyteSize_Results.csv"),
  Adipocyte_LR_Ventricle = paste0(output, prefix,
    "_Left_vs_Right_Adipocyte_Test_Results.csv"),
  Vessel_ATE_Pairwise = file.path(output,
    "ATE_group_comparisons.csv"),
  Vessel_ATE_Kruskal = paste0(output, prefix,
    "_ATE_group_comparisons_Kruskal.csv"),
  Nuclei_Models = file.path(output,
    paste0(prefix, "_Final_Model_Results.csv")),
  Nuclei_SmallHematoxylinDense = paste0(output, prefix,
    "_SmallHematoxylinDenseNuclei_Pairwise_Stats.csv")
)

# Read all stats and add a Source column
all_stats <- lapply(names(stats_files), function(name) {
  f <- stats_files[[name]]
  if (file.exists(f)) {
    df <- read.csv(f)
    df$Source <- name
    df
  } else {
    message("File not found: ", f)
    NULL
  }
})
all_stats_combined <- bind_rows(all_stats)

# --- Sheet 2: Phenotype correlations ---
corr_file <- file.path(output, paste0(prefix,
  "_Phenotype_Correlations_with_CI_and_Methods.csv"))
corr_stats <- if (file.exists(corr_file)) {
  read.csv(corr_file)
} else {
  correlation_statistics
}

# --- Sheet 3: Top 10 pathways ---
pathways_file <- paste0(output, prefix,
  "_Top10_Pathways_per_Phenotype.csv")
pathways_data <- if (file.exists(pathways_file)) {
  read.csv(pathways_file)
} else {
  final_top_pathways
}

# --- Write Excel ---
wb <- createWorkbook()
addWorksheet(wb, "All_Statistics")
addWorksheet(wb, "Phenotype_Correlations")
addWorksheet(wb, "Top10_Pathways")

writeData(wb, "All_Statistics", all_stats_combined)
writeData(wb, "Phenotype_Correlations", corr_stats)
writeData(wb, "Top10_Pathways", pathways_data)

excel_file <- file.path(output,
  paste0(prefix, "_Summary_Results.xlsx"))
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat("Excel summary saved to:", excel_file, "\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
