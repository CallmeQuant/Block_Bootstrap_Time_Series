# Block Bootstrap Time Series Analysis Project

## Project Overview

This project is focused on implementing and evaluating four block bootstrap methods for time series analysis including
+ Nonoverlapping block bootstrap (NBB)
+ Moving block bootstrap (MBB)
+ Stationary block bootstrap (SBB)
+ Circular block bootstrap (CBB)

The main objective is to estimate parameters of traditional time series models and assess their accuracy through standard error, bias, and confidence intervals.
Additionally, the project involves using bootstrap techniques to calculate forecast means and prediction intervals, comparing them with traditional normal approximation methods on real-world data.

## Project Structure

The project is organized into the following directories and files:

### `data_report_bootstrap/`
- **sales_train_evaluation.csv**: This file contains the real-world time series data used for analysis and forecasting.

### `exp_func/`
- **compare_calibrate_parameter.R**: This script compares and calibrates different parameters for the bootstrap methods and the traditional time series models.

### `utils/`
- **bootstrap_utils.R**: Contains utility functions related to implementing the block bootstrap algorithm and other summary tables creation.
- **forecast_bootstrap_utils.R**: Functions specifically for handling bootstrap forecasting tasks.
- **post_processing_utils.R**: Includes functions for post-processing the results, such as converting dataframes to latex tables or plotting.
- **run_experiment_bootstrap.R**: Runs experiments using bootstrap methods for parameter estimation and interval calculation.
- **run_experiment_forecast.R**: Executes the forecasting tasks, comparing bootstrap prediction intervals with normal approximation methods.
