fit_prob_bite_spline <- function(params_row,  # Add params_row as argument
                                 metadata_df, 
                                 base_inputs, 
                                 output_dir,
                                 incidence_trial,
                                 setup,
                                 country_to_run) {
  
  optimization_state <- list(
    counter = 0,
    best_rmse = Inf,
    best_coefs = NULL,
    starting_coefs_logit = setup$starting_coefs, # logit scale
    history = list(),
    start_time = Sys.time()
  )
  
  save_checkpoint <- function() {
    saveRDS(optimization_state, file = file.path(output_dir, "optim_checkpoint.rds"))
  }
  
  X <- setup$X
  starting_coefs <- setup$starting_coefs # these are on logit scale
  ns_basis = setup$ns_basis
  
  # Add small random noise to help exploration
  set.seed(123)  # for reproducibility
  starting_coefs <- starting_coefs + rnorm(length(starting_coefs), 0, 0.03)
  
  # Scale the design matrix to improve conditioning
  X_scale <- apply(X, 2, sd)
  X_scale[1] <- 1  # don't scale intercept
  X_scaled <- sweep(X, 2, X_scale, "/")
  
  # Adjust starting coefficients for scaling
  starting_coefs_scaled <- starting_coefs * X_scale
  
  objectivefunc <- function(coefs_scaled){
    # Unscale coefficients
    coefs <- coefs_scaled / X_scale
    
    # Print parameters to see if they're changing
    cat("Params received:", round(coefs_scaled[1:5], 2), "...\n")
    
    # Update counter
    optimization_state$counter <<- optimization_state$counter + 1
    iter <- optimization_state$counter
    
    # Transform to probability space
    linear_pred <- as.numeric(X %*% coefs) # same as predict() with model obj; gives log odds 
    bite_prob <- plogis(linear_pred) # converts log odds to probabilities 
    # cap at 0.5 maximum
    bite_prob <- pmin(bite_prob, 0.3)
    
    cat("  Probability summary: min=", min(bite_prob), 
        " mean=", mean(bite_prob), 
        " max=", max(bite_prob), "\n")
    cat("  % zeros:", 100 * mean(bite_prob < 1e-5), "%\n")
    
    if(max(bite_prob) < 1e-2) {
      cat("  ⚠️ All probabilities near zero, returning large error\n")
      return(1e10)
    }
    
    params_row_wprob <- params_row
    params_row_wprob$p_bite <- list(bite_prob)
    
    #Run model and transform data to get output of rmse
    o <- run_cohort_simulation(params_row_wprob, # this should have max smc kill rate, lambda, kappa, lag, simid, and pbite
                               metadata_df,
                               base_inputs,
                               output_dir = output_dir,
                               return_parasitemia = FALSE,
                               save_outputs = FALSE)
    message('finished simulation')
    
    if(is.null(o) || is.null(o$infection_records) || nrow(o$infection_records) == 0) {
      cat("  ⚠️ No infections detected - skipping processing, returning large error\n")
      return(1e10)
    }
    
    cat(sprintf("  ✅ %d infections detected\n", nrow(o$infection_records)))
    
    tryCatch({
      o$infection_records$sim_id <- params_row$sim_id
      
      # remove any infections that occurred within 7 days 
      infs <- o$infection_records %>% 
        group_by(rid) %>%
        arrange(rid, detection_day) %>%
        mutate(previous_detday = lag(detection_day),
               diff = detection_day - previous_detday) %>%
        filter(diff > 7 | is.na(diff)) %>% select(-diff, -previous_detday) %>%
        ungroup()
      message('removed infections within 7 days')
      
      infs_formatted <- format_model_output(model_data = infs, 
                                            cohort = country_to_run, 
                                            start_cohort = as.Date('2017-04-01'),
                                            simulation = params_row_wprob$sim_id)
      message('formatted infection df')
      
      inci <- get_incidence(df_children = metadata_df,
                            casedata = infs_formatted) %>%
        mutate(sim_id = params_row_wprob$sim_id)
      message('got incidence')
      
      metrics <- compare_incidence(incidence_model = inci, 
                                   incidence_trial = incidence_trial,
                                   output_dir = output_dir)
      rmse <- metrics$rmse
      # Ensure rmse is finite
      if(!is.finite(rmse)) {
        cat("  ⚠️ Non-finite RMSE, returning large value\n")
        return(1e10)
      }
      
    }, error = function(e) {
      cat("  Processing error:", e$message, "\n")
      return(1e10)  # Return from objectivefunc on error
    })
    
    # If we got an error and returned early, we need to check
    if(!exists("rmse")) return(1e10)
    
    # Update best if improved
    if(rmse < optimization_state$best_rmse && rmse < 1e9) {
      optimization_state$best_rmse <<- rmse
      optimization_state$best_coefs <<- coefs_scaled
      
      saveRDS(list(
        coefs = coefs_scaled,
        bite_prob = bite_prob,
        rmse = rmse,
        iter = iter,
        inci = inci,
        params_row = params_row_wprob
      ), file = file.path(output_dir, "best_so_far.rds"))

      
      cat(sprintf(" new best RMSE = %.4f\n", rmse))
    }
    
    # Store in history (limit size to prevent memory issues)
    if(iter <= 50) {  # Keep last 50
      optimization_state$history[[iter]] <<- list(
        coefs = coefs_scaled,
        rmse = rmse,
        time = Sys.time()
      )
    }
    
    
    # Save checkpoint every 5 iterations
    if(iter %% 5 == 0) {
      save_checkpoint()
      cat("Checkpoint saved at iteration", iter, "\n")
    }
    
    # Print progress
    elapsed <- difftime(Sys.time(), optimization_state$start_time, units = "hours")
    cat(sprintf("Iter %d: RMSE = %.4f (best: %.4f), Elapsed: %.1f hrs\n", 
                iter, rmse, optimization_state$best_rmse, elapsed))
    
    return(rmse)
  }
  
  optimresults <- nloptr::nloptr(x0 = starting_coefs_scaled,
                                 eval_f = objectivefunc,
                                 opts = list(algorithm = "NLOPT_LN_SBPLX", 
                                             maxeval = 35,
                                             ftol_rel = 1e-3,
                                             print_level = 2))
  message('got optim results')
  # Get final probabilities
  final_prob <- plogis(as.numeric(X %*% optimresults$solution))
  
  return(list(
    coefficients = optimresults$solution,
    probabilities = final_prob,
    X = X,
    ns_basis = ns_basis,
    optimresults = optimresults,
    best_rmse = optimization_state$best_rmse,
    history = optimization_state$history
  ))
}


resume_optimization <- function(output_dir, max_additional_evals = 30) {
  # Load checkpoint
  checkpoint_file <- file.path(output_dir, "optim_checkpoint.rds")
  if(!file.exists(checkpoint_file)) {
    stop("No checkpoint found in ", output_dir)
  }
  
  state <- readRDS(checkpoint_file)
  cat("Resuming from iteration", state$counter, "\n")
  cat("Best RMSE so far:", state$best_rmse, "\n")
  
  # You'll need to recreate the objective function with the saved state
  # This is trickier - you might want to save the entire environment
  
  return(state)
}
# If optimization crashes, just run:
# result <- resume_optimization()