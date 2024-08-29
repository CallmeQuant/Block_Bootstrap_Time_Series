source("./utils/bootstrap_utils.R")
source("./utils/post_processing_utils.R")
source("./utils/forecast_bootstrap_utils.R")

# List of libraries 
libraries <- c("feasts", "tidyverse", "tsibble", "fable", "bsts",
               "lubridate", "timetk", "timeDate", "here","dplyr",
               "readr", "vroom", "tibble", "tidyr", "gridExtra",
               "kableExtra", "grid", "parallel", "foreach",
               "doParallel", "stats", "boot", "meboot")

check_and_load_libraries(libraries)

# Change data path if needed 
data_path <- here("data_report_bootstrap")

# Preventing clash from Dplyr and MASS package for select function

select <- dplyr::select

train_df <- vroom(here(data_path, "sales_train_validation.csv"),
                  delim = ",",
                  col_types = cols())

# Retrieving special holidays
holidays <- c(USThanksgivingDay(2011:2016), USChristmasDay(2011:2016))@Data
holidays <- as.Date(holidays)

start_date = date("2011-01-29")

temp <- train_df |>
  group_by(cat_id, state_id) |>
  summarise(across(starts_with("d_"), sum), .groups = "drop") |>
  select(ends_with("id"), starts_with("d_")) |>
  pivot_longer(starts_with("d_"), names_to = "dates", values_to = "sales") |>
  mutate(dates = as.integer(str_remove(dates, "d_"))) |>
  mutate(dates = start_date + dates - 1) |>
  mutate(dates = ymd(dates)) |>
  mutate(dummy = if_else(dates %in% holidays, 1, 0)) |>
  mutate(dow = lubridate::wday(dates, label = TRUE))
# filter(!str_detect(as.character(dates), '..-12-25'))
# replace_na(list(sales =  0))

temp_ts <- temp |>
  as_tsibble(key = c(cat_id, state_id), index = dates)

food.CA <- temp_ts |>
  filter(cat_id == "FOODS" & state_id == "CA") |>
  select(dates, sales, dummy, dow) 

# Check for missing dates
full_dates <- seq.Date(min(food.CA$dates), max(food.CA$dates), by = "day")
missing_dates <- setdiff(full_dates, food.CA$dates)
missing_dates


splitter <- split_data(food.CA, "2 week")
train_data <- splitter$train
test_data <- splitter$test

fit_food_CA_week <- train_data |>
  model(
    stlf = decomposition_model(
      STL(sales ~ season(window = 'periodic') + trend(window = 15), 
          robust = TRUE),
      ARIMA(season_adjust ~ PDQ(0,0,0)),
      SNAIVE(season_year)
    ),
    # additive_ets = ETS(sales ~ error("A") + trend("A") + season("A")),
    # multiplicative_ets = ETS(sales ~ error("M") + trend("A") + season("M")),
    arima = ARIMA(sales))

# Normal approximation forecast
food_CA_fc_week <- fit_food_CA_week |>
  forecast(h = "2 week") |> as_tibble()

# Bootstrap forecast
food_CA_boot_fc_week <- bootstrap_all_models(fit_food_CA_week, "2 week")

# Change frome tibble to fable (necessary to use forecast library)
food_CA_boot_fc_week <- tibble_to_fable(food_CA_boot_fc_week)
food_CA_fc_week <- tibble_to_fable(food_CA_fc_week)

# Merged forecasts
merged_forecasts <- full_join(
  food_CA_boot_fc_week,
  food_CA_fc_week,
  by = c(".model", "dates"),
  suffix = c("_boot", "_regular")
)

plot_combined_forecasts(merged_forecasts, train_data, test_data, "stlf")

merged_forecasts

# Evaluate forecasts
fc_for_measure <- fit_food_CA_week |>
  forecast(h = "2 week")

acc_regular <- accuracy(fc_for_measure, test_data, 
                        measures = list(rmse=RMSE, 
                                        mae=MAE, 
                                        crps=CRPS))

fc_boot <- forecast(fit_food_CA_week, h = "2 week", 
                    bootstrap = TRUE, simulate = TRUE)

acc_boot <- accuracy(fc_boot, test_data, measures = list(rmse=RMSE, 
                                                         mae=MAE, 
                                                         crps=CRPS))
acc_regular <- acc_regular|>mutate(.type = "regular")
acc_boot <- acc_boot|>mutate(.type = "boot")

acc_combined <- bind_rows(acc_regular, acc_boot) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))