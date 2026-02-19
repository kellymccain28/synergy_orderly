# Load required library
library(tidyverse)

# Set the path to the outputs folder
outputs_path <- "R:/Kelly/synergy_orderly/src/sim_cohort_generic/outputs"

# Get all folders in the outputs directory
folders <- list.dirs(outputs_path, full.names = TRUE, recursive = FALSE)

# Filter for folders that match the pattern outputs_YYYY-MM-DD*
folder_names <- basename(folders)
output_folders <- folders[grepl("^outputs_\\d{4}-\\d{2}-\\d{2}", folder_names)]

# Initialize empty vectors to store results
folder_col <- character()
notes_col <- character()

# Loop through each folder and read sim_notes.txt
for (folder in output_folders) {
  folder_name <- basename(folder)
  notes_file <- file.path(folder, "sim_notes.txt")
  
  # Check if sim_notes.txt exists
  if (file.exists(notes_file)) {
    # Read the contents of sim_notes.txt
    notes_content <- paste(readLines(notes_file, warn = FALSE), collapse = "\n")
    
    # Add to results
    folder_col <- c(folder_col, folder_name)
    notes_col <- c(notes_col, notes_content)
  } else {
    # If file doesn't exist, add NA or a note
    folder_col <- c(folder_col, folder_name)
    notes_col <- c(notes_col, NA)
    warning(paste("sim_notes.txt not found in", folder_name))
  }
}

# Create data frame
results_df <- data.frame(
  folder_name = folder_col,
  sim_notes = notes_col,
  stringsAsFactors = FALSE
)

# Extract dates from folder names and filter for Dec 1, 2024 onwards
results_df <- results_df %>%
  mutate(date = as.Date(str_extract(folder_name, "\\d{4}-\\d{2}-\\d{2}"))) %>%
  filter(date >= as.Date("2024-12-01")) %>%
  arrange(desc(date), desc(folder_name)) %>%
  dplyr::select(-date)

# Write to CSV
output_csv <- file.path(outputs_path, "sim_notes_summary.csv")
write.csv(results_df, output_csv, row.names = FALSE)

# Print summary
cat("Extraction complete!\n")
cat("Total folders processed:", length(output_folders), "\n")
cat("Folders with sim_notes.txt:", sum(!is.na(results_df$sim_notes)), "\n")
cat("Output saved to:", output_csv, "\n")

# Display first few rows
print(head(results_df))