# Function to read in outputs from the simulations over a grid in sim_cohort_grid

read_in_outputs <- function(output_dir = 'simulation_outputs') {
  
  # Get all simulation directories
  sim_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  sim_dirs <- sim_dirs[grepl("^parameter_set_", basename(sim_dirs))]
  
  # Read in the parameters and baseline params 
  baseline_inputs <- readRDS(file.path(output_dir, "base_inputs.rds"))
  
  # Read all data types
  results <- list()
  
  # Infection records
  results$infection_records <- map_dfr(sim_dirs, function(sim_dir) {
    file_path <- file.path(sim_dir, "infection_records.rds")
    if (file.exists(file_path)) {
      data <- readRDS(file_path)
      data$sim_id <- basename(sim_dir)
      return(data)
    }
    return(NULL)
  })
  
  # Child metadata
  results$child_metadata <- map_dfr(sim_dirs, function(sim_dir) {
    file_path <- file.path(sim_dir, "child_metadata.rds")
    if (file.exists(file_path)) {
      data <- readRDS(file_path)
      data$sim_id <- basename(sim_dir)
      return(data)
    }
    return(NULL)
  })
  
  # Parasitemia data
  # results$parasitemia_data <- map_dfr(sim_dirs, function(sim_dir) {
  #   file_path <- file.path(sim_dir, "parasitemia_data.rds")
  #   if (file.exists(file_path)) {
  #     data <- readRDS(file_path)
  #     data$sim_id <- basename(sim_dir)
  #     return(data)
  #   }
  #   return(NULL)
  # })
  
  # Parameters (as a lookup table)
  results$parameters <- map_dfr(sim_dirs, function(sim_dir) {
    file_path <- file.path(sim_dir, "parameters.rds")
    if (file.exists(file_path)) {
      params <- readRDS(file_path)
      params$sim_id <- basename(sim_dir)
      return(as.data.frame(params))
    }
    return(NULL)
  })
  
  # # All expected efficacies ((this is for compare task))
  # results$model_outputs <- map_dfr(sim_dirs, function(sim_dir) {
  #   file_path <- file.path(sim_dir, "parameters.rds")
  #   if (file.exists(file_path)) {
  #     params <- readRDS(file_path)
  #     params$sim_id <- basename(sim_dir)
  #     return(as.data.frame(params))
  #   }
  #   return(NULL)
  # })
  
  return(list(baseline_inputs, results))
}


# read_simulation_results <- function(output_dir = , file_name, add_params = TRUE) {
#   
#   # Get all simulation directories  
#   sim_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
#   sim_dirs <- sim_dirs[grepl("^sim_", basename(sim_dirs))]
#   
#   # Read specified files
#   all_data <- map_dfr(sim_dirs, function(sim_dir) {
#     
#     target_file <- file.path(sim_dir, file_name)
#     
#     if (file.exists(target_file)) {
#       # Read the data
#       data <- readRDS(target_file)
#       
#       # Add simulation identifier
#       data$sim_id <- basename(sim_dir)
#       
#       # Optionally add parameters
#       if (add_params) {
#         param_file <- file.path(sim_dir, "parameters.rds")
#         if (file.exists(param_file)) {
#           params <- readRDS(param_file)
#           for (param_name in names(params)) {
#             data[[param_name]] <- params[[param_name]]
#           }
#         }
#       }
#       
#       return(data)
#     } else {
#       warning("No ", file_name, " found in ", sim_dir)
#       return(NULL)
#     }
#   })
#   
#   return(all_data)
# }