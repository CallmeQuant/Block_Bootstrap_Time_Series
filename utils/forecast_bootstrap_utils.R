wide_to_long <- function(data) {
  start_date <- date("2011-01-29")
  
  data |>
    select(id, starts_with("d_")) |>
    pivot_longer(starts_with("d_"), names_to = "dates", values_to = "sales") |>
    mutate(dates = as.integer(str_remove(dates, "d_"))) |>
    mutate(dates = start_date + dates - 1) |>
    mutate(dates = ymd(dates)) |>
    mutate(id = str_remove(id, "_validation"))
}

split_data <- function(data, split_duration) {
  if (is.numeric(split_duration)) {
    split_period <- split_duration
  } else if (is.character(split_duration)) {
    split_duration <- tolower(split_duration)
    # Parse duration string
    if (grepl("week", split_duration)) {
      split_period <- weeks(as.numeric(gsub("[^0-9]", "", split_duration)))
    } else if (grepl("day", split_duration)) {
      split_period <- days(as.numeric(gsub("[^0-9]", "", split_duration)))
    } else {
      stop("Unknown duration format. Please use 'weeks' or 'days'.")
    }
  } else {
    stop("split_duration must be numeric or character.")
  }
  
  data <- data |> arrange(dates)
  
  end_train <- max(data$dates) - split_period
  
  train_data <- data |> filter(dates <= end_train)
  test_data <- data |> filter(dates > end_train)
  
  return(list(train = train_data, test = test_data))
}

bootstrap_model <- function(model, model_name, forecast_horizon) {
  bootstrapped_forecast <- forecast(model, h = forecast_horizon, 
                                    bootstrap = TRUE, simulate = TRUE)
  
  forecast_data <- bootstrapped_forecast |>
    as.data.frame() |>
    mutate(.model = model_name)
  
  forecast_data <- as_tibble(forecast_data)
  
  return(forecast_data)
}

bootstrap_all_models <- function(fit, forecast_horizon) {
  model_names <- names(fit)
  models <- fit |> as.list()
  
  # Store results
  results <- list()
  
  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    model <- models[[i]]
    
    result <- bootstrap_model(model, model_name, forecast_horizon)
    
    # Store result in the list
    results[[model_name]] <- result
  }
  
  combined_results <- bind_rows(results)
  
  return(combined_results)
}

merge_forecast_bootstrap <- function(forecast_df, bootstrap_df) {
  # Ensure column names are consistent
  forecast_df <- forecast_df |>
    rename(
      forecast_dates = dates,
      forecast_model = .model,
      forecast_sales = sales,
      forecast_mean = .mean
    )
  
  bootstrap_df <- bootstrap_df |>
    rename(
      bootstrap_dates = dates,
      bootstrap_model = .model,
      bootstrap_sales = sales,
      bootstrap_mean = .mean
    )
  
  # Merge the data frames
  combined_results <- forecast_df |>
    left_join(
      bootstrap_df,
      by = c("forecast_dates" = "bootstrap_dates", 
             "forecast_model" = "bootstrap_model"),
      suffix = c("_forecast", "_bootstrap")
    )
  
  return(combined_results)
}

tibble_to_fable <- function(forecast_df, pred_interval = c(0.025, 0.975)) {
  # Grouped by .model and date
  forecast_df <- forecast_df |>
    group_by(.model, dates)
  
  # Calculate mean and prediction intervals
  summarized_results <- forecast_df |>
    reframe(
      .mean = mean(sales),
      .lower = quantile(sales, pred_interval[1]),
      .upper = quantile(sales, pred_interval[2]),
      sales = list(sales)
    )
  
  # Create a fable structure
  fable_result <- summarized_results |>
    mutate(
      sales = distributional::dist_sample(sales),
      .distribution = sales
    ) |>
    as_fable(index = dates, key = .model, 
             distribution = .distribution, response = "sales")
  
  return(fable_result)
}

plot_combined_forecasts <- function(forecast_data, 
                                    train_data, 
                                    test_data, 
                                    model_name, 
                                    pred_interval = 95) {
  # Filter data for the specified model
  model_data <- forecast_data |>
    filter(.model == model_name)
  
  # Combine train and test data
  full_data <- bind_rows(
    train_data |> mutate(type = "Train"),
    test_data |> mutate(type = "Test")
  )
  color_forecast = c(paste("Bootstrap", pred_interval, "% PI"),
                     paste("Regular", pred_interval, "% PI"))
  # Create the plot
  ggplot() +
    # Plot train data
    geom_line(data = full_data |> 
                filter(type == "Train", 
                       dates <= "2016-04-24" & dates >= "2016-03-24"),
              aes(x = dates, y = sales, color = "Train"), size = 0.7) +
    
    # Plot test data
    geom_line(data = full_data |> filter(type == "Test"), 
              aes(x = dates, y = sales, color = "Test"), size = 0.7) +
    
    # Plot bootstrap forecast
    geom_line(data = model_data, 
              aes(x = dates, y = .mean_boot, color = "Bootstrap Forecast"), 
              size = 0.7) +
    geom_ribbon(data = model_data, 
                aes(x = dates, ymin = .lower_boot, ymax = .upper_boot, 
                    fill = color_forecast[1]), 
                alpha = 0.2) +
    
    # Plot regular forecast
    geom_line(data = model_data, 
              aes(x = dates, y = .mean_regular, color = "Regular Forecast"), 
              size = 0.7) +
    geom_ribbon(data = model_data, 
                aes(x = dates, ymin = .lower_regular, ymax = .upper_regular, 
                    fill = color_forecast[2]), 
                alpha = 0.2) +
    
    # Customize the plot
    labs(title = paste("Forecast Comparison for", model_name),
         y = "Sales",
         x = "Date") +
    scale_color_manual(name = "Data",
                       values = c("Train" = "black", 
                                  "Test" = "darkgreen",
                                  "Bootstrap Forecast" = "blue", 
                                  "Regular Forecast" = "red")) +
    scale_fill_manual(name = "Prediction Interval",
                      values = setNames(c("blue", "red"), color_forecast)) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.margin = margin()) +
    guides(color = guide_legend(order = 1),
           fill = guide_legend(order = 2))
}