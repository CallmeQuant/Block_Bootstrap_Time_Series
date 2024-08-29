source("./exp_func/compare_calibrate_parameter.R")
source("./utils/bootstrap_utils.R")
source("./utils/post_processing_utils.R")

# List of libraries 
libraries <- c("feasts", "tidyverse", "tsibble", "fable", "bsts",
               "lubridate", "timetk", "timeDate", "here","dplyr",
               "readr", "vroom", "tibble", "tidyr", "gridExtra",
               "kableExtra", "grid", "parallel", "foreach",
               "doParallel", "stats", "boot", "meboot")

check_and_load_libraries(libraries)

######################################
# Example AR/MA parameter estimation #
######################################
coefs <- list(AR = c(0.5, -0.3, 0.2),
              MA = c(-0.3, 0.4, 0.5))
# AR(1)
ar1_results <- compare_param_estimation(coefs,"AR", 1)
# MA(1)
ma1_results <- compare_param_estimation(coefs,"MA", 1)

# ARMA(1,1)
arma11_results <- compare_arma_estimation(coefs_ar = c(-0.9), 
                                          coefs_ma = c(0.5), 
                                          order_ar = 1, order_ma = 1)

################################
# Example AR/MA ACF estimation #
################################
coefs <- list(AR = c(0.5, -0.3, 0.2),
              MA = c(-0.3, 0.4, 0.5))
# AR(1)
ar1_results_acf <- compare_acf_estimation(coefs,"AR", 1)

# MA(1)
ma1_results_acf <- compare_acf_estimation(coefs,"MA", 1)

# ARMA(1,1)
arma11_results_acf <- compare_arma_acf_estimation(coefs_ar = c(-0.9), 
                                                  coefs_ma = c(0.5), 
                                                  order_ar = 1, order_ma = 1)

###########################
# Example AR/MA calibrate #
###########################

coefs <- list(AR = c(0.5, -0.3), MA = c(-0.3, 0.4))
n_values <- c(100)
ar1_sim_results <- compare_param_estimation_calibrate(coefs, model_type = "AR", 
                                                      order = 1, 
                                                      n_values=n_values,
                                                      sim_reps = 1)

ma1_sim_results <- compare_param_estimation_calibrate(coefs, model_type = "MA", 
                                                      order = 1, 
                                                      n_values=n_values,
                                                      sim_reps = 1)
#########################
# Example SM+MBB method #
#########################
library(Mcomp)
data(M3)
period <- 1 # yearly period 
# period <- 12 # monthly period
yearly_M3 <- subset(M3, "yearly")
# monthly_M3 <- subset(M3, "monthly")
data_ts <- as.numeric(yearly_M3$N0001$x)
# data_ts <- as.numeric(monthly_M3$N1499$x)
data_boot_mbb <- bld.mbb.bootstrap(ts(data_ts, freq = period), 100)
data_plot <- tibble(
  Value = unlist(data_boot_mbb),
  ID = rep(1:100, each = length(data_ts)),
  Time = rep(1:length(data_ts), 100)
)

ggplot(data_plot) +
  geom_line(aes(x = Time, y = Value, group = ID, color = "Bootstrap Samples"), 
            alpha = 0.5) +
  geom_line(data = filter(data_plot, ID == 1), 
            aes(x = Time, y = Value, color = "Original Series"),
            alpha = 0.9, size = 0.8) +
  scale_color_manual(values = c("Bootstrap Samples" = "steelblue", 
                                "Original Series" = "darkorange")) +
  labs(color = "Legend") +
  mytheme

#####################################
# Example time series decomposition #
#####################################
df <- read_csv("https://raw.githubusercontent.com/philipperemy/
               keras-tcn/master/tasks/monthly-milk-production-pounds-p.csv")
timeindex <- seq(as.Date("1962-01-01"), as.Date("1975-12-01"), by = "1 month")
df <- df |>
  rename(milksale = milk_production_pounds) |>
  mutate(month = yearmonth(timeindex))


ts_df <- tsibble(df)

ts_df |> model(STL(milksale ~ trend(window = 28) + 
                     season(window = "periodic"),
                   robust = TRUE)) |> 
  components() |>
  autoplot() -> milk_decomp_plot

milk_decomp_plot + mytheme + theme(axis.title.x = element_blank(),
                                   axis.title.y = element_blank(),
                                   plot.title = element_blank())
ts_df