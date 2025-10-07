
# ------------------------------
# Full PLS modeling script  CLR ABUNDANCES
# ------------------------------

# Version random splitting 2

library(tidyverse)
library(tidymodels)
library(plsmod)     # for mixOmics PLS/sPLS
library(mixOmics)   # sPLS backend
library(future)

# ----------------FUNCTIONS-----------------

# --- Patch for plsmod mixOmics bug ("incorrect number of dimensions")
fix_multi_numeric_preds <- dget("./custom-functions/fix_multi_numeric_preds.R")
# Register the patch in the main session
assignInNamespace("multi_numeric_preds", fix_multi_numeric_preds, ns = "plsmod")

## UPDATED VERSION 17 SEP 2025 to introduce randomness into the test window splitting (and get mean values for repeated test fitting)
# To mitigate sensitivity to a specific test-window placement, we introduced stochastic tie-breaking into the variance-balanced window selection.
# When multiple candidate windows exhibited variance close to the median, one was selected at random (with fixed seeds for reproducibility).
# This ensured that evaluation results were not overly dependent on a single deterministic choice of test windows,
# thereby providing a more robust assessment of model generalizability across temporal contexts.


choose_random_windows <- function(df, outcome, assess, n_windows = 2, seed = NULL) {
  # assess = number of consecutive samples per test window
  # n_windows = how many test windows to select

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(df)
  starts <- 1:(n - assess + 1)

  selected <- list()
  used_idx <- integer()

  while (length(selected) < n_windows && length(starts) > 0) {
    pick <- sample(starts, 1)   # random start
    this_window <- pick:(pick + assess - 1)
    selected <- append(selected, list(this_window))
    used_idx <- c(used_idx, this_window)

    # prevent overlaps (embargo 1 either side)
    embargo_range <- (min(this_window) - assess):(max(this_window) + assess)
    starts <- setdiff(starts, embargo_range)
  }

  return(selected)
}

# ---- Make an h-v block (purged) rset on TRAIN ONLY ----
# assess = length of each validation window (in rows)
# h      = embargo size on each side of the validation window (in rows)
# reps   = number of random windows (repeats)
make_hv_block_rset <- function(df, assess = 12, h = 2, reps = 50) {
  n <- nrow(df)
  stopifnot(assess > 0, assess < n)

  starts <- 1:(n - assess + 1)
  starts <- sample(starts, size = reps, replace = (reps > length(starts)))

  splits <- lapply(seq_along(starts), function(i) {
    s <- starts[i]
    assess_idx  <- s:(s + assess - 1)

# Purge/embargo around assessment window
    embargo_lo <- max(1, s - h)
    embargo_hi <- min(n, s + assess - 1 + h)
    analysis_idx <- setdiff(seq_len(n), embargo_lo:embargo_hi)

   rsample::make_splits(
      list(analysis = analysis_idx, assessment = assess_idx),
      data = df
    )
  })
    rsample::manual_rset(splits, ids = paste0("hv_", seq_along(splits)))
}


# -------------- SETTINGS ------------------

data_files <- list(
  "species_allclr"      = "./data/mdata.species.clr_all.csv"  # <<--- chosen 12 Sep 2025
  # Add more as needed
)


outcomes <- c("A", "T", "DN", "T2", "C",
              "C2", "H", "H2",
              "V", "V2")  # Change these to your actual response variable names

#outcomes <- c("A_roll6", "T_roll6", "DNA_g_kgVS_roll6", "T2_roll6", "C_roll6",
#              "C2_roll6", "H_roll6", "H2_roll6",
#              "V_roll6", "V2_roll6")  # Change these to your actual response variable names

#outcomes <- c(
#   'h2__hydrogenotrophic_methanogens_roll6', "su__sugar_degraders_fermenters_roll6", "aa__amino_acid_degraders_roll6",
#   "c4__butyrate_and_valerate_degraders_roll6", "fa__long_chain_fatty_acid_degraders_roll6", "pro__propionate_degraders_roll6",
#   "ac__acetoclastic_methanogens_roll6"
#)

#outcomes <- c(
#  'h2__hydrogenotrophic_methanogens', "su__sugar_degraders_fermenters", "aa__amino_acid_degraders",
#  "c4__butyrate_and_valerate_degraders", "fa__long_chain_fatty_acid_degraders", "pro__propionate_degraders",
#  "ac__acetoclastic_methanogens"
#)

id_col    <- "Digester"
time_col  <- "Sample_Date"

plan(multisession, workers = parallel::detectCores() - 1)

# Parameters
n_repeats   <- 10   # how many different test-window selections
assess_len  <- 12   # length of test/CV window
h_embargo   <- 2    # purge size
reps_cv     <- 20   # inner h–v block CV reps

# -------------- OUTPUT STORAGE ------------------
results <- list()

# -------------- MAIN LOOP ------------------

for (data_name in names(data_files)) {
  dataset_path <- data_files[[data_name]]
  df <- read.csv(dataset_path)

  library(zoo)
  # for numerical responses only (not for OTU abundance predictors) create a 6-day rolling average
  # A 6-day sample window approximates the sludge turnover in the digesters (2 weeks).
# df <- df %>%
#    arrange(Digester, Sample_Date) %>%   # sort within each group
#    group_by(Digester) %>%               # group by digester (or other variable)
#    mutate(across(
#      where(is.numeric) & !all_of("Sample_Date") & !contains("_otu"),  # if numeric but not for otu-labeled columns apply the rollapply function
#      ~ rollapply(.x, width = 6, FUN = mean,                  # take the mean of 6 samples in each selected column
#                  align = "right", fill = NA, na.rm = TRUE),  # right aligned = no future leakage future leakage
#      .names = "{.col}_roll6"
#    )) %>%
#    ungroup() %>%
#    dplyr::select(Digester, Sample_Date, ends_with("_roll6"),contains("_otu")) # Keep only group/date + rolled columns

  # Auto-detect OTU cols for this dataset
  otu_cols <- grep("_otu$|FLASV", colnames(df), value = TRUE)

  if (length(otu_cols) == 0) {
    stop(paste("No OTU columns detected in dataset", data_name))
  }

  # Sort by digester and time
  df <- df %>%
    arrange(.data[[id_col]], .data[[time_col]])

  for (outcome in outcomes) {

    all_repeats <- list()

    df.sel <- df %>%
      dplyr::select(all_of(outcome), id_col, time_col,
                    contains("_otu") | contains("FLASV")) %>%
      arrange(.data[[time_col]]) %>%
      drop_na()

    n_samples <- nrow(df.sel)
    # Verbose update
    message("Processing dataset: ", data_name,
            " | Outcome: ", outcome,
            " | n = ", nrow(df.sel))


    # -------------------------
    # Outer: balanced test windows
    # -------------------------

    for (r in seq_len(n_repeats)) {
      message("Dataset: ", data_name, " | Outcome: ", outcome, " | Repeat: ", r)

    test_windows <- choose_random_windows(
      df.sel,
      outcome  = outcome,
      assess   = assess_len,
      n_windows = 2,   # number of test windows
      seed = r
    )

    test_idx <- sort(unique(unlist(test_windows)))

    test_dat  <- df.sel[test_idx, ]
    train_dat <- df.sel[-test_idx, ]

    # -------------------------
    # Inner: h–v block CV on training
    # -------------------------
    resamps <- make_hv_block_rset(
      df  = train_dat %>% arrange(.data[[time_col]]),
      assess = assess_len, h = h_embargo, reps = reps_cv
    )

    # -------------------------
    # Recipe
    # -------------------------
      frm <- as.formula(paste0("`", outcome, "` ~ ."))
      rec <- recipe(frm, data = train_dat) %>%
        update_role(all_of(id_col),   new_role = "id") %>%
        update_role(all_of(time_col), new_role = "time") %>%
        #update_role(.block_id,        new_role = "cv_block") %>%
        #step_log(any_of(!!otu_cols), offset = 1) %>%
        step_zv(all_predictors()) %>%
        step_normalize(all_predictors())

      # -------------------------
      # Model spec
      # -------------------------
      pls_spec <- parsnip::pls(mode = "regression",
                               num_comp       = tune(),
                               predictor_prop = tune()
      ) %>%
        set_engine("mixOmics")

      wf <- workflow() %>%
        add_recipe(rec) %>%
        add_model(pls_spec)

      grid <- grid_regular(
        num_comp(range = c(2L, 20L)),
        predictor_prop(range = c(0.6, 1.0)),
        levels = 8
      )

      tuned <- tune_grid(
        wf,
        resamples = resamps,
        grid      = grid,
        metrics   = metric_set(rmse, rsq, mae, mape),
        control   = control_grid(save_pred = TRUE,
                                 allow_par = TRUE,
                                 parallel_over = "resamples")
      )

      # Skip if tuning failed
      if (nrow(collect_metrics(tuned)) == 0) {
        warning(paste("Tuning failed for", data_name, outcome, "block", b_out))
        next
      }

      best <- select_best(tuned, metric = "rmse")

      # Refit on training
      final_wf <- finalize_workflow(wf, best)
      final_fit <- fit(final_wf, data = train_dat)

      # ---- CV average (already averaged by collect_metrics)
      cv_metrics <- tuned %>%
        collect_metrics() %>%
        dplyr::filter(.metric %in% c("rmse","rsq","mae", "mape"),
                      num_comp == best$num_comp,
                      predictor_prop == best$predictor_prop) %>%
        transmute(
          dataset  = "cv(train)",
          .metric  = .metric,
          estimate = mean,
          repeat_id = r
        )


      # ---- Refit training metrics ----
      train_metrics <- metric_set(rmse, rsq, mae, mape)(
        predict(final_fit, new_data = train_dat) %>%
          bind_cols(truth = train_dat[[outcome]]),
        truth = truth, estimate = .pred
      ) %>%
        transmute(
        dataset   = "train(refit)",
        .metric   = .metric,
        estimate  = .estimate,
        repeat_id = r
      )


      # ---- Test(outer)
      test_metrics <- metric_set(rmse, rsq, mae,mape)(
        predict(final_fit, new_data = test_dat) %>%
          bind_cols(truth = test_dat[[outcome]]),
        truth = truth, estimate = .pred
      ) %>%
        transmute(
        dataset   = "test(outer)",
        .metric   = .metric,
        estimate  = .estimate,
        repeat_id = r
      )

      # Combine / Store
      all_repeats[[r]] <- bind_rows(cv_metrics, train_metrics, test_metrics) %>%
        mutate(Dataset = data_name,
               Outcome = outcome,
               n_samples = n_samples)
    }

    # Save averaged metrics across repeats
    results[[paste(data_name, outcome, sep = "::")]] <-
      bind_rows(all_repeats) %>%
      group_by(Dataset, Outcome, dataset, .metric, n_samples) %>%
      summarise(mean = mean(estimate, na.rm = TRUE),
                sd   = sd(estimate, na.rm = TRUE),
                .groups = "drop")
  }
}

plan(sequential)

# -------------- SUMMARISE RESULTS ------------------

# ---- Final summary table ----


summary_tbl <- bind_rows(results) %>%
  group_by(Dataset, Outcome) %>%
  mutate(n_repeats = n_repeats) %>%
  ungroup() %>%
  pivot_wider(
    names_from = c(dataset, .metric),
    values_from = c(mean, sd),
    names_sep = "_"
  )

summary_tbl2 <- summary_tbl

# Collapse mean ± sd into one column per metric
  # Collapse mean ± sd for each metric column, then drop sd_ columns
  mean_cols <- grep("^mean_", names(summary_tbl2), value = TRUE)
  sd_cols   <- sub("^mean_", "sd_", mean_cols)

  summary_tbl2[mean_cols] <- Map(function(m, s) {
    m <- signif(m, 2)   # 1 significant figure
    s <- signif(s, 2)
    paste0(m, " +/- ", s)
  }, summary_tbl2[mean_cols], summary_tbl2[sd_cols])

summary_tbl2 <- summary_tbl2 %>% dplyr::select(-all_of(sd_cols))

# --- Desired column order ---
dataset_order <- c("train(refit)", "cv(train)", "test(outer)")
metric_order  <- c("rmse", "mae", "mape", "rsq" )

col_order <- c(
  "Dataset", "Outcome", "n_repeats",
  unlist(lapply(dataset_order, function(d) {
    paste0("mean_", d, "_", metric_order)
  }))
)

summary_tbl2 <- summary_tbl2 %>%
  dplyr::select(all_of(col_order))

#readr::write_csv(summary_tbl, "pls_model_summary_roll6.csv")

# Scatterplot: CV vs Test R²
plot_df <- summary_tbl %>%
  dplyr::select(Dataset, Outcome,
                cv_rsq   = `mean_cv(train)_rsq`,
                test_rsq = `mean_test(outer)_rsq`) %>%
  tidyr::drop_na()

p <- ggplot(plot_df, aes(x = cv_rsq, y = test_rsq, label = Outcome)) +
  geom_point(size = 3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    #title = "Optimism gap between CV and Test performance",
    x = "CV R²",
    y = "Test R²"
  ) +
  coord_equal(xlim = c(0,1), ylim = c(0,1))

p

# Convert summary_tbl into long format for plotting
# Example for R² only
plot_df_err <- summary_tbl %>%
  dplyr::select(Dataset, Outcome,
                mean_train = `mean_train(refit)_rsq`,  sd_train = `sd_train(refit)_rsq`,
                mean_cv    = `mean_cv(train)_rsq`,     sd_cv    = `sd_cv(train)_rsq`,
                mean_test  = `mean_test(outer)_rsq`,   sd_test  = `sd_test(outer)_rsq`) %>%
  tidyr::pivot_longer(
    cols = -c(Dataset, Outcome),
    names_to = c(".value", "dataset"),
    names_sep = "_"
  )

p2 <- ggplot(plot_df_err, aes(x = Outcome, y = mean, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                position = position_dodge(width = 0.8), width = 0.3) +
  theme_minimal(base_size = 14) +
  labs(
    title = "clr_hvblock_random_aquasim_roll6",
    x = "Outcome",
    y = expression(R^2),
    fill = "Data"
  ) +
  theme(axis.text.x = element_text(angle = 75, hjust = 1))

p2
