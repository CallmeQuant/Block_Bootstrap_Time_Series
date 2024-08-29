compare_param_estimation <- function(coefs, model_type, order, n = 120, 
                                     block_length = NULL, reps = 999,
                                     overwrite_fig = TRUE,
                                     fig_path = NULL) {
  set.seed(123)
  if (is.null(block_length)){block_length = floor(sqrt(n))}
  # Simulate process
  if (model_type == "AR") {
    coef <- coefs$AR[1:order]
    y <- arima.sim(n = n, model = list(ar = coef))
  } else if (model_type == "MA") {
    coef <- coefs$MA[1:order]
    y <- arima.sim(n = n, model = list(ma = coef))
  } else {
    stop("Invalid model type. Use 'AR' or 'MA'.")
  }
  
  # Plot the simulated time series
  plot(y, type = "l", main = paste("Simulated", 
                                   model_type, "(", order, ") Process"))
  
  # Estimation function
  estimate_params <- function(data) {
    if (model_type == "AR") {
      model <- arima(data, order = c(order, 0, 0))
    } else {
      model <- arima(data, order = c(0, 0, order))
    }
    return(c(coef(model)[1:order]))
  }
  
  
  bootf <- function (data) {
    if (model_type == "AR") {
      model <- arima(data, order = c(order, 0, 0))
    } else {
      model <- arima(data, order = c(0, 0, order))
    }
    return(coef(model)[1:order])}
  
  # Perform bootstraps
  mbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed")
  mbb_result <- bootstrap_and_store(mbb_boot, "MBB")
  
  sb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "geom")
  sb_result <- bootstrap_and_store(sb_boot, "SB")
  
  cbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed", 
                     endcorr = TRUE)
  cbb_result <- bootstrap_and_store(cbb_boot, "CBB")
  
  nbb_results <- nbb_param_estimate(y, block_length, reps, estimate_params)
  nbb_result <- as_tibble(nbb_results) |> mutate(type = "NBB")
  
  # Combine all results
  combined_result <- bind_rows(mbb_result, sb_result, cbb_result, nbb_result)
  
  # Summarize mean, 95% interval
  result_summary <- combined_result |>
    group_by(type) |>
    summarise(across(starts_with("param"), 
                     list(lower = ~quantile(., 0.025, na.rm = TRUE),
                          upper = ~quantile(., 0.975, na.rm = TRUE),
                          mean = ~mean(., na.rm = TRUE))))
  
  # Create and arrange plots
  plots <- lapply(1:order, function(i) {
    create_plot(combined_result, result_summary, 
                paste("Parameter", i), paste0("param", i))
  })
  ncol <- ifelse(order == 1, 1, ceiling(sqrt(order)))
  nrow <- ceiling(order / ncol)
  main_title <- paste(model_type, order, "Parameter Estimation")
  
  # Adjust layout for odd number of plots >= 3
  if (order == 3) {
    layout_matrix = rbind(c(1,2), 3)
    
    combined_plot <- grid.arrange(
      grobs = plots,
      ncol = ncol,
      nrow = nrow,
      layout_matrix = layout_matrix,
      top = textGrob(main_title, gp = gpar(fontsize = 12, fontface = "bold")))
  } else{
    combined_plot <- grid.arrange(
      grobs = plots,
      ncol = ncol,
      nrow = nrow,
      top = textGrob(main_title, gp = gpar(fontsize = 12, fontface = "bold")))
  }
  
  # Save plot
  # Create fig_path to be relative to the project root directory
  if (is.null(fig_path)) {
    fig_path <- here("figures") 
  }
  dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
  file_name <- file.path(fig_path, paste0(model_type, order, "param.png"))
  if (!file.exists(file_name) || 
      overwrite_fig == TRUE){
    cat(paste("Saving plot", model_type, "_", order, "\n"))
    ggsave(filename = paste(model_type, order, "param.png", sep = ""), 
           plot = combined_plot, 
           path = fig_path, 
           width = 10, height = 8)}
  else{cat(paste("Plot already exists. Skipping save.\n"))}
  
  # Summary for table result 
  table_result <- table_summary_result(combined_result,
                                       result_summary, 
                                       coef)
  return(list(summary = result_summary, result = combined_result, 
              table = table_result))
}

compare_param_estimation_calibrate <- function(coefs, model_type, order, 
                                               n_values = c(50, 100, 500, 1000), 
                                               block_length = NULL, 
                                               sim_reps = 1000, 
                                               boot_reps = 999,
                                               coverage_steps = seq(0, 1, 0.1),
                                               overwrite_fig=TRUE) {
  set.seed(123)
  
  results <- list()
  
  for (n in n_values) {
    if (is.null(block_length)) {
      block_length <- floor(sqrt(n))
    }
    
    coverage_counts <- matrix(0, nrow = length(coverage_steps), ncol = 4)
    interval_lengths <- matrix(0, nrow = length(coverage_steps), ncol = 4)
    colnames(coverage_counts) <- colnames(interval_lengths) <- 
      c("MBB", "SB", "CBB", "NBB")
    
    for (i in 1:sim_reps) {
      # Simulate process
      if (model_type == "AR") {
        coef <- coefs$AR[1:order]
        y <- arima.sim(n = n, model = list(ar = coef))
      } else if (model_type == "MA") {
        coef <- coefs$MA[1:order]
        y <- arima.sim(n = n, model = list(ma = coef))
      }
      
      # Estimate function
      estimate_params <- function(data) {
        if (model_type == "AR") {
          model <- arima(data, order = c(order, 0, 0))
        } else {
          model <- arima(data, order = c(0, 0, order))
        }
        return(coef(model)[1:order])
      }
      
      # Perform bootstraps
      mbb_boot <- tsboot(y, estimate_params, R = boot_reps, 
                         l = block_length, sim = "fixed")
      sb_boot <- tsboot(y, estimate_params, R = boot_reps, 
                        l = block_length, sim = "geom")
      cbb_boot <- tsboot(y, estimate_params, R = boot_reps, 
                         l = block_length, sim = "fixed", endcorr = TRUE)
      nbb_boot <- nbb_param_estimate(y, block_length, 
                                     boot_reps, estimate_params)
      
      # Extract bootstrap results
      boot_results <- list(
        MBB = mbb_boot$t,
        SB = sb_boot$t,
        CBB = cbb_boot$t,
        NBB = nbb_boot  
      )
      
      # Calculate intervals and check coverage
      for (c_idx in seq_along(coverage_steps)) {
        alpha <- coverage_steps[c_idx]
        for (j in 1:4) {
          intervals <- apply(boot_results[[j]], 2, quantile, 
                             probs = c(alpha/2, 1-alpha/2))
          interval_lengths[c_idx, j] <- interval_lengths[c_idx, j] + 
            mean(diff(intervals))
          coverage_counts[c_idx, j] <- coverage_counts[c_idx, j] + 
            all(coef >= intervals[1,] & coef <= intervals[2,])
        }
      }
    }
    
    # Calculate empirical frequencies and median lengths
    empirical_freq <- coverage_counts / sim_reps
    median_lengths <- interval_lengths / sim_reps
    
    # Calculate KPI
    kpi <- empirical_freq / sqrt(median_lengths + 1)
    
    # Store results
    results[[as.character(n)]] <- list(
      empirical_freq = empirical_freq,
      median_lengths = median_lengths,
      kpi = kpi
    )
  }
  
  # summary dataframe
  summary_df <- do.call(rbind, lapply(names(results), function(n) {
    data.frame(
      n = as.numeric(n),
      coverage = rep(coverage_steps, 4),
      method = rep(c("MBB", "SB", "CBB", "NBB"), each = length(coverage_steps)),
      empirical_freq = c(results[[n]]$empirical_freq),
      median_length = c(results[[n]]$median_lengths),
      kpi = c(results[[n]]$kpi)
    )
  }))
  summary_df |> 
    group_by(method)|> 
    mutate(empirical_freq = sort(empirical_freq)) -> summary_df
  
  # Create calibration plot
  calibration_plot <- ggplot(summary_df, aes(x = empirical_freq, y = coverage, 
                                             color = method)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    facet_wrap(~ n, scales = "free") +
    labs(y = "Theoretical Coverage", x = "Empirical Frequency", 
         title = "Calibration Plot by Sample Size and Method") +
    mytheme
  
  # Save plot
  # Create fig_path to be relative to the project root directory
  if (is.null(fig_path)) {
    fig_path <- here("figures") 
  }
  dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
  file_name <- file.path(fig_path, paste0(model_type, order, "calibration.png"))
  if (!file.exists(file_name) || 
      overwrite_fig == TRUE){
    cat(paste("Saving plot", model_type, "_", order, "\n"))
    ggsave(filename = paste(model_type, order, "calibration.png", sep = ""), 
           plot = calibration_plot, 
           path = fig_path, 
           width = 10, height = 8)}
  else{cat(paste("Plot already exists. Skipping save.\n"))}
  
  # Summary for table result 
  return(list(summary = summary_df, calibration_plot = calibration_plot))
}

compare_arma_estimation <- function(coefs_ar, coefs_ma, order_ar, order_ma, n = 120, 
                                    block_length = NULL, reps = 999, 
                                    overwrite_fig = TRUE) {
  set.seed(123)
  if (is.null(block_length)){ block_length = floor(sqrt(n)) }
  
  # Simulate ARMA process
  y <- arima.sim(n = n, model = list(ar = coefs_ar[1:order_ar], 
                                     ma = coefs_ma[1:order_ma]))
  
  # Plot the simulated time series
  plot(y, type = "l", main = paste("Simulated ARMA(", 
                                   order_ar, ",", order_ma, ") Process"))
  
  # Estimation function
  estimate_params <- function(data) {
    model <- arima(data, order = c(order_ar, 0, order_ma))
    return(c(coef(model)[1:order_ar], 
             coef(model)[(order_ar+1):(order_ar+order_ma)]))
  }
  
  bootf <- function (data) {
    model <- arima(data, order = c(order_ar, 0, order_ma))
    return(coef(model)[1:(order_ar + order_ma)])
  }
  
  # Perform bootstraps
  mbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed")
  mbb_result <- bootstrap_and_store(mbb_boot, "MBB")
  
  sb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "geom")
  sb_result <- bootstrap_and_store(sb_boot, "SB")
  
  cbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed", 
                     endcorr = TRUE)
  cbb_result <- bootstrap_and_store(cbb_boot, "CBB")
  
  nbb_results <- nbb_param_estimate(y, block_length, reps, estimate_params)
  nbb_result <- as_tibble(nbb_results) |> mutate(type = "NBB")
  
  # Combine all results
  combined_result <- bind_rows(mbb_result, sb_result, cbb_result, nbb_result)
  
  # Summarize mean, 95% interval
  result_summary <- combined_result |>
    group_by(type) |>
    summarise(across(starts_with("param"), 
                     list(lower = ~quantile(., 0.025, na.rm = TRUE),
                          upper = ~quantile(., 0.975, na.rm = TRUE),
                          mean = ~mean(., na.rm = TRUE))))
  
  # Create and arrange plots
  plots <- lapply(1:(order_ar + order_ma), function(i) {
    create_plot(combined_result, result_summary, 
                paste("Parameter", i), paste0("param", i))
  })
  ncol <- ifelse(order_ar + order_ma == 1, 1, 
                 ceiling(sqrt(order_ar + order_ma)))
  nrow <- ceiling((order_ar + order_ma) / ncol)
  main_title <- paste("ARMA", order_ar, order_ma, "Parameter Estimation")
  
  # Adjust layout for odd number of plots >= 3
  if (order_ar + order_ma == 3) {
    layout_matrix = rbind(c(1, 2), 3)
    
    combined_plot <- grid.arrange(
      grobs = plots,
      ncol = ncol,
      nrow = nrow,
      layout_matrix = layout_matrix,
      top = textGrob(main_title, gp = gpar(fontsize = 12, fontface = "bold")))
  } else {
    combined_plot <- grid.arrange(
      grobs = plots,
      ncol = ncol,
      nrow = nrow,
      top = textGrob(main_title, gp = gpar(fontsize = 12, fontface = "bold")))
  }
  
  # Save plot
  # Create fig_path to be relative to the project root directory
  if (is.null(fig_path)) {
    fig_path <- here("figures") 
  }
  dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
  file_name <- file.path(fig_path, paste0("ARMA", order_ar, 
                                          order_ma, "param.png"))
  if (!file.exists(file_name) || 
      overwrite_fig == TRUE){
    cat(paste("Saving plot ARMA_", order_ar, "_", order_ma, "\n"))
    ggsave(filename = paste("ARMA", order_ar, order_ma, "param.png", sep = ""), 
           plot = combined_plot, 
           path = fig_path, 
           width = 10, height = 8)
  } else {
    cat(paste("Plot already exists. Skipping save.\n"))
  }
  
  # Summary for table result 
  table_result <- table_summary_result(combined_result,
                                       result_summary, 
                                       c(coefs_ar[1:order_ar], coefs_ma[1:order_ma]))
  return(list(summary = result_summary, result = combined_result, 
              table = table_result))
}

# Function to run ACF estimation experiments
compare_acf_estimation <- function(coefs, model_type, order, n = 120, 
                                   block_length = NULL, reps = 999, lags = 5,
                                   overwrite_fig = TRUE) {
  set.seed(123)
  if (is.null(block_length)) { block_length = floor(sqrt(n)) }
  
  # Simulate process
  if (model_type == "AR") {
    coef <-coefs$AR[1:order]
    y <- arima.sim(n = n, model = list(ar = coef))
  } else if (model_type == "MA") {
    coef <- coefs$MA[1:order]
    y <- arima.sim(n = n, model = list(ma = coef))
  } else {
    stop("Invalid model type. Use 'AR' or 'MA'.")
  }
  
  # Plot the simulated time series
  plot(y, type = "l", main = paste("Simulated", 
                                   model_type, "(", order, ") Process"))
  
  # Estimate ACF (used by NBB )
  estimate_acf <- function(data) {
    acf_values <- acf(data, plot = FALSE, lag.max = lags)$acf
    return(acf_values[2:(lags + 1)])
  }
  
  bootf <- function(data) {
    acf_values <- acf(data, plot = FALSE, lag.max = lags)$acf
    return(acf_values[2:(lags + 1)])
  }
  
  # Perform bootstraps
  mbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed")
  mbb_result <- bootstrap_and_store(mbb_boot, "MBB")
  
  sb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "geom")
  sb_result <- bootstrap_and_store(sb_boot, "SB")
  
  cbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed", 
                     endcorr = TRUE)
  cbb_result <- bootstrap_and_store(cbb_boot, "CBB")
  
  nbb_results <- nbb_param_estimate(y, block_length, reps, estimate_acf)
  nbb_result <- as_tibble(nbb_results) |> mutate(type = "NBB")
  
  # Combine all results
  combined_result <- bind_rows(mbb_result, sb_result, cbb_result, nbb_result)
  
  acf_plot <- create_acf_plot(y, combined_result, lags, model_type, order)
  
  # Save plot
  # Create fig_path to be relative to the project root directory
  if (is.null(fig_path)) {
    fig_path <- here("figures") 
  }
  dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
  file_name <- file.path(fig_path, paste0(model_type, order, "acf.png"))
  if (!file.exists(file_name) || 
      overwrite_fig == TRUE){
    cat(paste("Saving plot", model_type, "_", order, "\n"))
    ggsave(filename = paste(model_type, order, "acf.png", sep = ""), 
           plot = acf_plot, 
           path = fig_path, 
           width = 10, height = 8)}
  else{cat(paste("Plot already exists. Skipping save.\n"))}
  
  # Summary table 
  true_params <- estimate_acf(y)
  table_result <- table_summary_result_acf(combined_result, true_params)
  
  return(list(result = combined_result, table_se = table_result$SE,
              table_bias = table_result$Bias))
}

# Main function for running ARMA ACF estimation experiments
compare_arma_acf_estimation <- function(coefs_ar, coefs_ma, order_ar, order_ma, 
                                        n = 120, 
                                        block_length = NULL, 
                                        reps = 999, 
                                        lags = 5, 
                                        overwrite_fig = TRUE) {
  set.seed(123)
  if (is.null(block_length)) { block_length = floor(sqrt(n)) }
  
  # Simulate ARMA process
  y <- arima.sim(n = n, model = list(ar = coefs_ar[1:order_ar], 
                                     ma = coefs_ma[1:order_ma]))
  
  # Plot the simulated time series
  plot(y, type = "l", main = paste("Simulated ARMA(", 
                                   order_ar, ",", order_ma, ") Process"))
  
  # Define function to estimate ACF (used by NBB )
  estimate_acf <- function(data) {
    acf_values <- acf(data, plot = FALSE, lag.max = lags)$acf
    return(acf_values[2:(lags + 1)])
  }
  
  bootf <- function(data) {
    acf_values <- acf(data, plot = FALSE, lag.max = lags)$acf
    return(acf_values[2:(lags + 1)])
  }
  
  # Perform bootstraps
  mbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed")
  mbb_result <- bootstrap_and_store(mbb_boot, "MBB")
  
  sb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "geom")
  sb_result <- bootstrap_and_store(sb_boot, "SB")
  
  cbb_boot <- tsboot(y, bootf, R = reps, l = block_length, sim = "fixed", 
                     endcorr = TRUE)
  cbb_result <- bootstrap_and_store(cbb_boot, "CBB")
  
  nbb_results <- nbb_param_estimate(y, block_length, reps, estimate_acf)
  nbb_result <- as_tibble(nbb_results) |> mutate(type = "NBB")
  
  # Combine all results
  combined_result <- bind_rows(mbb_result, sb_result, cbb_result, nbb_result)
  
  # Create ACF plot
  acf_plot <- create_acf_plot(y, combined_result, lags, paste("ARMA", 
                                                              order_ar,
                                                              order_ma),
                              order = paste(order_ar, order_ma, sep = ",")
  )
  
  # Save plot
  # Create fig_path to be relative to the project root directory
  if (is.null(fig_path)) {
    fig_path <- here("figures") 
  }
  dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)
  file_name <- file.path(fig_path, paste0("ARMA", order_ar, 
                                          order_ma, "acf.png"))
  if (!file.exists(file_name) || 
      overwrite_fig == TRUE){
    cat(paste("Saving plot ARMA_", order_ar, "_", order_ma, "\n"))
    ggsave(filename = paste("ARMA", order_ar, order_ma, "acf.png", sep = ""), 
           plot = acf_plot, 
           path = fig_path, 
           width = 10, height = 8)
  } else {
    cat(paste("Plot already exists. Skipping save.\n"))
  }
  
  # Summary table 
  true_params <- estimate_acf(y)
  table_result <- table_summary_result_acf(combined_result, 
                                           true_params)
  
  return(list(result = combined_result, table_se = table_result$SE,
              table_bias = table_result$Bias))
}