# Function to check, install, and load libraries
check_and_load_libraries <- function(libs) {
  for (lib in libs) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      message(paste("Installing", lib))
      install.packages(lib)
    }
    library(lib, character.only = TRUE)
  }
}

df_to_latex <- function(df){
  df |>
    kable(format = 'latex', booktabs = TRUE)
}

write_tables_to_latex <- function(results_list, file_path) {
  # Open the file connection
  con <- file(file_path, "w")
  
  # Loop through each result object
  for (model in names(results_list)) {
    result <- results_list[[model]]
    
    # Write identifier for SE table
    cat(paste0("% ", model, " Standard Error Table\n"), file = con)
    
    # Convert SE table to latex
    se_latex <- df_to_latex(result$table_se)
    cat(se_latex, file = con, sep = "\n")
    cat("\n\n", file = con)  # Add some space between tables
    
    # Write identifier for Bias table
    cat(paste0("% ", model, " Bias Table\n"), file = con)
    
    # Convert Bias table to latex
    bias_latex <- df_to_latex(result$table_bias)
    cat(bias_latex, file = con, sep = "\n")
    cat("\n\n", file = con)  # Add some space between models
  }
  
  # Close the file connection
  close(con)
}

library(ggplot2)
mytheme <- theme_classic() +
  theme(
    panel.background = element_blank(),
    strip.background = element_rect(colour = NA, fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "plain", size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 9.5),
    text = element_text(family = "sans")
  )

