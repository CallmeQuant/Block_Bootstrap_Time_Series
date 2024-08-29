non_overlapping_block_bootstrap <- function(y, block_length, reps) {
  n <- length(y)
  num_blocks <- ceiling(n / block_length)
  
  bootstrap_samples <- matrix(NA, nrow = reps, ncol = n)
  
  for (i in 1:reps) {
    block_indices <- sample(1:num_blocks, num_blocks, replace = TRUE)
    bootstrap_sample <- numeric(n)
    
    for (j in 1:num_blocks) {
      start_idx <- (block_indices[j] - 1) * block_length + 1
      end_idx <- min(start_idx + block_length - 1, n)
      
      # Calculate where to place the block in the bootstrap sample
      sample_idx <- (j - 1) * block_length + 1
      sample_end <- min(sample_idx + (end_idx - start_idx), n)
      
      # Length of the block to be copied
      block_length_actual <- end_idx - start_idx + 1
      
      if (sample_end > n) {
        sample_end <- n
      }
      
      # only fill up to the remaining space in the bootstrap sample
      if (sample_idx <= n) {
        fill_length <- min(block_length_actual, n - sample_idx + 1)
        bootstrap_sample[sample_idx:(sample_idx + fill_length - 1)] <- 
          y[start_idx:(start_idx + fill_length - 1)]
      }
    }
    
    bootstrap_samples[i, ] <- bootstrap_sample
  }
  
  return(bootstrap_samples)
}



# NBB with parameter estimation
nbb_param_estimate <- function(y, block_length, reps, estimate_func, ...) {
  nbb_samples <- non_overlapping_block_bootstrap(y, block_length, reps)
  params <- matrix(NA, nrow = reps, ncol = length(estimate_func(y, ...)))
  
  for (i in 1:reps) {
    # Remove NAs and estimate parameters for each bootstrap sample
    sample_ts <- na.omit(nbb_samples[i,])
    params[i,] <- estimate_func(sample_ts, ...)
  }
  
  colnames(params) <- paste0("param", seq_len(ncol(params)))
  return(params)
}

# Function to perform bootstrap, 
# extract values for all parameters, 
# and store them in the result tibble
bootstrap_and_store <- function(boot_object, method_name) {
  colnames(boot_object$t) <- paste0("param", seq_len(ncol(boot_object$t)))
  
  boot_values_tibble <- as_tibble(boot_object$t)
  
  boot_values_tibble <- boot_values_tibble |>
    mutate(type = method_name)
  
  return(boot_values_tibble)
}

# Function to create distribution of bootstrap samples
create_plot <- function(data, summary, param_name, param_column) {
  ggplot(data, aes(x = .data[[param_column]], fill = type)) +
    geom_density(alpha = 0.5) +
    geom_vline(data = summary, 
               aes(xintercept = .data[[paste0(param_column, "_mean")]], 
                   color = type),
               linetype = "solid", linewidth = 1) +
    geom_vline(data = summary, 
               aes(xintercept = .data[[paste0(param_column, "_lower")]], 
                   color = type),
               linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = summary, 
               aes(xintercept = .data[[paste0(param_column, "_upper")]], 
                   color = type),
               linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ type, scales = "free") +
    labs(title = "",
         x = param_name,
         y = "Density") +
    mytheme +
    theme(legend.position = "none")
}

# Function to create summary table for parameter estimation
table_summary_result <- function(result, result_summary, true_params) {
  order <- length(true_params)
  methods <- c("NBB", "MBB", "SB", "CBB")
  df <- data.frame(matrix(nrow = 4 * order, ncol = length(methods)))
  colnames(df) <- methods
  rownames(df) <- c(
    paste0("Param", rep(1:order, 4), "_", 
           rep(c("Mean", "CI", "SE", "Bias"), each = order))
  )
  
  for (method in methods) {
    # Extract mean estimates
    means <- result_summary |>
      filter(type == method) |>
      select(ends_with("_mean")) |>
      unlist()
    df[paste0("Param", 1:order, "_Mean"), method] <- round(means, 3)
    
    # Extract confidence intervals
    lower <- result_summary |>
      filter(type == method) |>
      select(ends_with("_lower")) |>
      unlist()
    upper <- result_summary |>
      filter(type == method) |>
      select(ends_with("_upper")) |>
      unlist()
    df[paste0("Param", 1:order, "_CI"), method] <- 
      paste0("(", round(lower, 3), ", ", round(upper, 3), ")")
    
    # Calculate standard errors
    se <- result |>
      filter(type == method) |>
      select(starts_with("param")) |>
      summarise(across(everything(), sd)) |>
      unlist()
    df[paste0("Param", 1:order, "_SE"), method] <- round(se, 3)
    
    # Calculate bias
    bias <- means - true_params
    df[paste0("Param", 1:order, "_Bias"), method] <- round(bias, 3)
  }
  
  return(df)
}

# Function to create summary table for ACF estimation
table_summary_result_acf <- function(result, true_params) {
  methods <- c("NBB", "MBB", "SB", "CBB")
  params <- paste0("param", 1:length(true_params))
  
  # SE table
  se_table <- result |>
    group_by(type) |>
    summarise(across(all_of(params), ~sd(.x, na.rm = TRUE))) |>
    rename_with(~str_replace(., "param", "Lag "), starts_with("param")) |>
    rename(Method = type) |>
    mutate(across(starts_with("Lag "), ~round(.x, 3)))
  
  # Bias table
  bias_table <- result |>
    group_by(type) |>
    summarise(across(all_of(params), ~mean(.x, na.rm = TRUE))) |>
    mutate(across(
      all_of(params), 
      ~.x - true_params[which(params == cur_column())])) |>
    rename_with(~str_replace(., "param", "Lag "), starts_with("param")) |>
    rename(Method = type) |>
    mutate(across(starts_with("Lag "), ~round(.x, 3)))
  
  se_table <- se_table |> slice(match(methods, Method))
  bias_table <- bias_table |> slice(match(methods, Method))
  
  return(list(SE = se_table, Bias = bias_table))
}

# Function to create ACF plot with bootstrap estimates
create_acf_plot <- function(original_series, bootstrap_results, max_lag, 
                            model_type, order) {
  # Compute original ACF
  original_acf <- acf(original_series, plot = FALSE, lag.max = max_lag)
  acf_data <- data.frame(
    Lag = 1:max_lag,
    ACF = original_acf$acf[2:(max_lag + 1)]
  )
  
  # Prepare bootstrap ACF estimates
  bootstrap_acf_data <- bootstrap_results |>
    group_by(type) |>
    summarise(across(starts_with("param"), mean, na.rm = TRUE)) |>
    pivot_longer(cols = starts_with("param"), 
                 names_to = "Lag", values_to = "BootstrapACF")
  
  # Convert Lag to numeric
  bootstrap_acf_data$Lag <- as.numeric(gsub("param", "", 
                                            bootstrap_acf_data$Lag))
  
  # Create plot
  acf_plot <- ggplot() +
    geom_segment(data = acf_data, aes(x = Lag, xend = Lag, y = 0, yend = ACF), 
                 color = "blue") +
    geom_point(data = acf_data, aes(x = Lag, y = ACF), 
               color = "blue", size = 2) +
    geom_line(data = bootstrap_acf_data, aes(x = Lag, y = BootstrapACF, 
                                             color = type, group = type)) +
    geom_point(data = bootstrap_acf_data, aes(x = Lag, y = BootstrapACF, 
                                              color = type)) +
    labs(title = paste("ACF Plot with Bootstrap Estimates (", model_type, 
                       order, ")", sep = ""),
         x = "Lag", y = "ACF") +
    mytheme +
    theme(legend.position = "bottom") +
    scale_color_discrete(name = "Bootstrap Method") +
    geom_hline(yintercept = 0)
  
  return(acf_plot)
}

# MBB version from forecast library
MBB <- function(x, window_size) {
  
  bx <- array(0, (floor(length(x)/window_size)+2)*window_size)
  for (i in 1:(floor(length(x)/window_size)+2)) {
    c <- sample(1:(length(x)-window_size+1),1)
    bx[((i-1)*window_size+1):(i*window_size)] <- x[c:(c+window_size-1)]
  }
  start_from <- sample(0:(window_size-1),1) + 1
  bx[start_from:(start_from+length(x)-1)]
}

# Smooth version of STL-MBB method
smo.bootstrap <- function(data, num, block_size=NULL, alpha = 0.05)
{
  freq <- frequency(data)
  if(is.null(block_size))
  {
    block_size <- ifelse(freq > 1, 2*freq, min(8, floor(length(data)/ 2)))
  }
  
  data_boot <- list()
  data_boot[[1]] <- data # the first series is the original one
  
  if (num>1) {
    # Box-Cox transformation
    lambda <- BoxCox.lambda(data, lower=0.01, upper=1)
    data_bc <- BoxCox(data, lambda)
    
    if (freq>1) {
      # STL decomposition
      data_stl <- stl(ts(data_bc, frequency=freq), "per", 
                      robust = TRUE)$time.series
      seasonal <- data_stl[,1]
      trend <- data_stl[,2]
      # Smoothing the remainder part by Holt-Winters exponential smoothing
      remainder <- ts(c(0, HoltWinters(data_stl[,3], 
                                       alpha = alpha, 
                                       beta = FALSE, 
                                       gamma = FALSE)$fitted[,1]), freq = freq)
      
    } else {
      # Loess
      trend <- 1:length(data)
      suppressWarnings(data_loess <- loess(data_bc ~ trend, 
                                           span=6/length(data), degree=1))
      seasonal <- rep(0, length(data))
      trend <- data_loess$fitted
      remainder <- data_loess$residuals
    }
    
    # Bootstrap remainder using MBB then add trend seasonal part and 
    # inverse box-cox
    for (i in 2:num) {
      data_boot[[i]] <- InvBoxCox(trend + seasonal + 
                                    MBB(remainder, block_size), lambda)
    }
  }
  data_boot
}

# Function to create time series plot
create_ts_plot <- function(data, title) {
  ggplot(data.frame(x = 1:length(data), y = data), aes(x = x, y = y)) +
    geom_line() +
    labs(title = title, x = "Time", y = "Value") +
    mytheme
}

# Function to create ACF plot
create_acf_plot <- function(data, title) {
  acf_data <- acf(data, plot = FALSE)
  acf_df <- with(acf_data, data.frame(lag, acf))
  
  ggplot(acf_df, aes(x = lag, y = acf)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_segment(aes(xend = lag, yend = 0)) +
    labs(title = title, x = "Lag", y = "ACF") +
    mytheme
}