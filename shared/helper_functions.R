# Helper functions to run model, format data, plot data, run simulation


# Wrapper function to run, format, and get time to threshold 
run_process_model <- function(n_particles = 1L,
                              n_threads = 1L,
                              PEV_on,
                              SMC_on,
                              t_inf_vax,
                              tt , # sequence of timesteps (2 days for each 1 step)
                              VB = 1e6,
                              num_bites = 1,
                              det_mode = FALSE,
                              vmin = 0,
                              infection_start_day, # external time that infection begins
                              SMC_time, # vector with same length as smc_kill_vec
                              SMC_kill_vec, # per-parasite kill rate per 2-day timestep
                              alpha_ab = 1.32, # default values from White 2013 
                              beta_ab = 6.62, # default values from White 2013
                              tboost1 = 364,
                              tboost2 = 729
                              ){ 
  
  # Run model  out <- run_model(args)
  out <- run_model(n_particles = n_particles,
                   n_threads = n_threads,
                   PEV_on = PEV_on,
                   SMC_on = SMC_on,
                   t_inf_vax = t_inf_vax,
                   tt = tt,
                   VB = VB,
                   vmin = vmin,
                   num_bites = num_bites,
                   det_mode = det_mode,
                   infection_start_day = infection_start_day, # external time that infection begins
                   SMC_time = SMC_time,
                   SMC_kill_vec = SMC_kill_vec, 
                   alpha_ab = alpha_ab, # default values from White 2013 
                   beta_ab = beta_ab, # default values from White 2013
                   tboost1 = tboost1,
                   tboost2 = tboost2)
  
  # Format data 
  out_formatted <- format_data(out,
                               tt = tt,
                               infection_start_day = infection_start_day, 
                               n_particles = n_particles)
  
  # Find first day threshold is reached 
  threshold_day <- get_ttoinf(out_formatted) %>% pull(time_withinhost2)
  
  return(list(trajectory = out_formatted, 
              threshold_day = threshold_day))
}

#' @param n_particles number of particles to create
#' @param n_threads number of threads to use in parallelisable calculations
#' @param PEV_on 0 for no PEV, 1 if PEV is delivered
#' @param SMC_on 0 for no SMC, 1 if SMC is delivered
#' @param t_inf_vax time of infection relative to vaccination (i.e. 20 means the vaccine was given 20 days prior to infectious bite); 
# should be >=0 if vaccination will protect from infection; (this influences ab titre at time of infection)
#' @param tt sequence of timesteps in 2-day increments 
#' @param VB volume of blood in microL; 1e6 is for a small child (1L); adult male is 5L
#' @param det_mode logical, if model should be run in deterministic mode
run_model <- function(n_particles = 1L,
                      n_threads = 1L,
                      PEV_on,
                      SMC_on,
                      t_inf_vax,
                      tt,
                      VB,
                      vmin = 0,
                      num_bites = 1,
                      det_mode = FALSE,
                      infection_start_day = 0, # external time that infection begins 
                      SMC_time, 
                      SMC_kill_vec,
                      alpha_ab = 1.32, # default values from White 2013 
                      beta_ab = 6.62, # default values from White 2013
                      tboost1 = 364, # timesteps after 3rd dose that the first booster is delivered
                      tboost2 = 729# timesteps after 1st booster that the second booster is delivered 
){
  
  # Get antibody curve over 3 years of cohort study 
  ts <- seq(0,365*5)
  phases <- ifelse(ts < tboost1, 1,
                   ifelse(ts < tboost2, 2, 3))
  ab <- antibody_titre(t = ts,
                       phase = phases,
                       peak1 = c(621,0.35) ,
                       peak2 = c(277, 0.35),
                       peak3 = c(277,0.35),
                       duration1 = c(45,16),
                       duration2 = c(591,245),
                       rho1 = c(2.37832, 1.00813),
                       rho2 = c(1.034, 1.027),
                       rho3 = c(1.034, 1.027),
                       t_boost1 = tboost1,
                       t_boost2 = tboost2)
  ab_user <- if(PEV_on == 1 & t_inf_vax > 0) ab[t_inf_vax] else 0 # on day 0, the ab haven't started yet 
  
  # Set parameters
  pars <- list(PEV_on = PEV_on,
               SMC_on = SMC_on,
               ab_user = ab_user,
               VB = VB,
               vmin = vmin,
               num_bites = num_bites,
               tt = max(tt),
               infection_start_day = infection_start_day, # external time that infection begins 
               SMC_time = unlist(SMC_time),
               SMC_kill_vec = unlist(SMC_kill_vec),
               alpha_ab = alpha_ab, # default values from White 2013 
               beta_ab = beta_ab # default values from White 2013
  )
  
  sys <- dust_system_create(gen_bs, 
                            pars, 
                            n_particles = n_particles,
                            n_threads = n_threads,
                            # seed = 12345L,
                            deterministic = det_mode)
  # Set initial values using initial() equations in model
  dust_system_set_state_initial(sys)
  # Get starting states 
  init_states <- dust_unpack_state(sys, dust_system_state(sys))
  # Calculate the median value of merozoites initiating infection across all particles
  # median_mero_init <- median(init_states$PB)
  # Get number of infections blocked 
  # sum(init_states$PB == 0)
  
  # Run model 
  bs_model <- dust_system_simulate(sys, tt)
  
  out <- dust_unpack_state(sys, bs_model)
  out$mero_init_out <- init_states$mero_init_out
  
  return(out)
}

format_data <- function(out, tt, infection_start_day, n_particles){
  # tt <- seq(1, ts, by = 1)
  
  if(n_particles == 1){
    df <- as.data.frame(out[names(out) != "buffer"]) # removing the buffer column which is not necessary 
    df <- df %>% select(-starts_with('p_history'))
    df$time <- tt
    df$run <- "run1"
    df <- df %>%
      rename(parasites = PB, 
             innate_imm = sc,
             genadaptive_imm = sm,
             varspecific_imm = sv) %>%
      mutate(time_withinhost = time, 
             time_withinhost2 = time_withinhost*2,
             time = time * 2 + infection_start_day)
  }
  
  if(n_particles > 1){
    # pull out only the blood stage 
    bl <- as.data.frame(t(out$PB))
    colnames(bl) <- paste0('run',seq(1:n_particles))
    bl$time <- tt
    df_long <- bl %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'parasites') #%>%
    # use threshold for the maximum value of parasites (limit of detection)
    # mutate(parasites = ifelse(parasites > 5e7, 5e7, 
    #                           ifelse(parasites < 1, 0, parasites)))
    
    sc <- as.data.frame(t(out$sc))
    colnames(sc) <- paste0('run',seq(1:n_particles))
    sc$time <- tt
    sc_long <- sc %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'innate_imm') 
    
    sm <- as.data.frame(t(out$sm))
    colnames(sm) <- paste0('run',seq(1:n_particles))
    sm$time <- tt
    sm_long <- sm %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'genadaptive_imm') 
    
    sv <- as.data.frame(t(out$sv))
    colnames(sv) <- paste0('run',seq(1:n_particles))
    sv$time <- tt
    sv_long <- sv %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'varspecific_imm') 
    
    gr <- as.data.frame(t(out$growth))
    colnames(gr) <- paste0('run',seq(1:n_particles))
    gr$time <- tt
    gr_long <- gr %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'growth') 
    
    m <- as.data.frame(t(out$m))
    colnames(m) <- paste0('run',seq(1:n_particles))
    m$time <- tt
    m_long <- m %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'm')
    
    smc <- as.data.frame(t(out$SMC_kill_rateout))
    colnames(smc) <- paste0('run',seq(1:n_particles))
    smc$time <- tt
    smc_long <- smc %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'smcrate')
    
    psmckill <- as.data.frame(t(out$prob_smckill))
    colnames(psmckill) <- paste0('run',seq(1:n_particles))
    psmckill$time <- tt
    psmckill_long <- psmckill %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'smc_prob')
    
    nsmckill <- as.data.frame(t(out$numkillSMC))
    colnames(nsmckill) <- paste0('run',seq(1:n_particles))
    nsmckill$time <- tt
    nsmckill_long <- nsmckill %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'nkillsmc')
    
    meros <- as.data.frame(t(out$mero_init_out))
    colnames(meros) <- paste0('run',seq(1:n_particles))
    meros$time <- 1
    meros_long <- meros %>% 
      pivot_longer(cols = starts_with('run'),
                   names_to = 'run',
                   values_to = 'mero_init_out')

    df <- df_long %>% 
      left_join(sc_long) %>% 
      left_join(sm_long) %>%
      left_join(sv_long) %>%
      left_join(gr_long) %>%
      left_join(m_long) %>%
      left_join(smc_long) %>%
      left_join(psmckill_long) %>%
      left_join(nsmckill_long) %>%
      left_join(meros_long) %>%
      # Fix time to be in outside days (*2) - not including inf start day as below because it doesn't match with ts since start of bs x labels 
      mutate(time_withinhost = time, # original timing in within-host model (2-day timesteps)
             time_withinhost2 = time_withinhost*2, # converted to 1-day timesteps 
             time = time * 2 + infection_start_day) # 1 day timesteps + day of start of BS infection 
      # Fix time to be dependent on when the infection was, relative to the start of the season / vaccination
      # mutate(time_ext = time*2 + t_inf_vax - 1)#,
    # infection_start = infection_start_day) # because the infection start day is in days and time is in timestpes, 
  }
  return(df)
}

make_plots <- function(df){
  # multiplier <- if(tstep == 1) 2 else if (tstep == 2) 1
  colorsplot <- c("#4123E8","#A6BF19")
  # 7 + time * multiplier
  df <- df %>%
    group_by(run) %>%
    mutate(cleared = ifelse(any(parasites < 1e-5), 'cleared', 'not cleared'),
           cleared = factor(cleared, levels = c('not cleared', 'cleared')),
           time = time, # timesteps are by 2 days, so we want to get the two day steps in plots
           detectable = ifelse(any(parasites >= threshold), 'detectable', 'not detectable')
    ) %>%
    ungroup() %>%
    group_by(time_withinhost2) %>%
    mutate(median_parasites = median(parasites))
  
  p <- ggplot(df) + 
    geom_line(aes(x = time_withinhost2, y = parasites, group = run, color = detectable), alpha = 0.6, linewidth = 0.6) + #, color = cleared
    geom_line(aes(x = time_withinhost2, y = median_parasites, color = detectable), linewidth = 0.6)+
    geom_hline(aes(yintercept = threshold), linetype = 2, color = 'darkred', linewidth = 1) + # this is the detection limit (followiung Challenger et al.)
    geom_hline(aes(yintercept = 1e-5), linetype = 2, color = 'darkgreen', linewidth = 1) + # this is the clearance threshold
    scale_y_log10(labels = scales::label_log(),
                  guide = "axis_logticks") +
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    scale_color_manual(values = colorsplot) +
    labs(x = 'Days since start of blood stage',
         y = 'PRBCs',
         caption = paste0('proportion of runs with no infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% 
                             count() %>% pull(n))/ n_particles)) + 
    theme_bw() + 
    labs() + 
    theme(legend.position = 'none')
  
  scplt <- ggplot(df )+#%>% filter(detectable == 'not detectable')) + 
    geom_line(aes(x = time_withinhost2, y = innate_imm, group = run, color = detectable), alpha = 0.7) + #, color = cleared
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'Innate immunity',
         caption = paste0('proportion of runs with no infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% 
                             count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  smplt <- ggplot(df)+#%>% filter(detectable == 'not detectable')) + 
    geom_line(aes(x = time_withinhost2, y = genadaptive_imm, group = run, color = detectable), alpha = 0.7) + #, color = cleared
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'General adaptive immunity',
         caption = paste0('proportion of runs with no infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  svplt <- ggplot(df)+#%>% filter(detectable == 'not detectable')) + 
    geom_line(aes(x = time_withinhost2, y = varspecific_imm, group = run, color = detectable), alpha = 0.7) + #, color = cleared
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'Var-specific immunity',
         caption = paste0('proportion of runs with no infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  grplt <- ggplot(df)+# %>% filter(detectable == 'not detectable')) + 
    geom_line(aes(x = time_withinhost2, y = growth, group = run, color = detectable), alpha = 0.7) + #
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'Growth rate (m * sc * sm * sv)',
         caption = paste0('proportion of bites that do not lead to infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  mplt <- ggplot(df)+#  %>% filter(cleared == 'not cleared')) + 
    geom_line(aes(x = time_withinhost2, y = m, group = run, color = detectable), alpha = 0.7) + #
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'm',
         caption = paste0('proportion of bites that do not lead to infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  smcplt <- ggplot(df) + 
    geom_line(aes(x = time_withinhost2, y = smcrate, group = run, color = detectable), alpha = 0.7, linewidth = 1) + #
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time), 14)),#(max(df$time)+7) * multiplier
      limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'SMC kill rate',
         caption = paste0('proportion of bites that do not lead to infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  psmckillplt <- ggplot(df) + 
    geom_line(aes(x = time_withinhost2, y = smc_prob, group = run, color = detectable), alpha = 0.7, linewidth = 1) + #
    scale_x_continuous(breaks = c(0, 7, seq(14, max(df$time_withinhost2), 28)),#(max(df$time)+7) * multiplier
                       limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'SMC per parasite/uL kill probability',
         caption = paste0('proportion of bites that do not lead to infection at time 0 = ', 
                          (df %>% filter(time == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    theme(legend.position = 'none')
  
  nsmckillplt <- ggplot(df ) + 
    geom_line(aes(x = time_withinhost2, y = nkillsmc, group = run, color = detectable), alpha = 0.7, linewidth = 0.6) + #
    scale_x_continuous(#breaks = c(0, 7, seq(14, max(df$time_withinhost2), 14)),#(max(df$time)+7) * multiplier
                       limits = c(0, max(df$time_withinhost2)))+#(max(df$time)+7)*multiplier
    labs(x = 'Days since start of blood stage',
         y = 'Parasites/uL killed by SMC',
         caption = paste0('proportion of bites that do not lead to infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == 1, parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    scale_y_log10(labels = scales::label_log(),
                  guide = "axis_logticks") +
    theme(legend.position = 'none')
  
  meroinitplt <- ggplot(df) + 
    geom_boxplot(aes(x = detectable, y = mero_init_out, color = detectable), alpha = 0.7) + #
    geom_jitter(aes(x = detectable, y = mero_init_out, color = detectable), alpha = 0.7) + #
    labs(x = 'Infection status',
         y = 'Merozoites initating infection',
         caption = paste0('proportion of bites that do not lead to infection at time 0 = ', 
                          (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)) + 
    scale_color_manual(values = colorsplot) +
    theme_bw() + 
    scale_y_log10(labels = scales::label_log(),
                  guide = "axis_logticks") +
    theme(legend.position = 'none')
  
  message('prop runs with no infection = ', (df %>% filter(time_withinhost2 == min(df$time_withinhost2), parasites == 0) %>% count() %>% pull(n))/ n_particles)
  
  return(list('pb'= p, 
              'sc'=scplt, 
              'sm'=smplt,
              'sv'=svplt,
              'growth'=grplt, 
              'data'=df, 
              'basicm'=mplt, 
              'smckill'=smcplt, 
              'probsmc'=psmckillplt, 
              'nsmc'=nsmckillplt, 
              'meroinit'=meroinitplt))
}

# Time to infection
# time between start of blood-stage and the day it reaches the threshold and can be detected
# time_withinhost2 is the external/cohort time value of the number of days since the beginning of the blood-stage
get_ttoinf <- function(df){
  df <- df %>%
    group_by(run) %>%
    mutate(threshold_reached = any(parasites >= threshold)) %>%
    filter(if (threshold_reached[1]) parasites >= threshold else time_withinhost2 == max(time_withinhost2)) %>%# filter(parasites >=10) %>%
    slice_min(time_withinhost2, n = 1) %>%
    # slice(1) %>%
    ungroup() %>%
    mutate(time_withinhost2 = ifelse(!threshold_reached, NA_real_, as.numeric(time_withinhost2)))
  
  return(df)
}

# Function to calculate time since SMC dose 
# where timings is the vector of days since April 1, 2017 that SMC was delivered for an individual child
# and days is a vector of days (0:end of cohort sim) 
# outputs at each day how long it has been since the last dose which can be used to calculate the kill rate due to SMC per day 
calc_time_since_dose <- function(timings, days) {
  suppressWarnings(sapply(days, function(d) {
    last_dose <- max(timings[timings <= d])
    ifelse(is.finite(last_dose), d - last_dose, NA)
  }))
}


#Calculate lagged vectors of probabilty of infectious bite for all unique lags
# This is for an input into the cohort sim
calc_lagged_vectors <- function(prob_data, lags, start_date = as.Date('2017-04-01'), 
                                end_date = '2020-04-01', burnints) {
  
  lag_list <- purrr::map(lags, function(lag_val) {
    # prob_lagged <- prob_data %>% 
    #   mutate(prob_lagged = dplyr::lag(prob_infectious_bite, n = lag_val),
    #          date_lagged = dplyr::lag(date, n = lag_val))
    
    if (lag_val >= 0) {
      prob_lagged <- prob_data %>% 
        mutate(
          prob_lagged = dplyr::lag(prob_infectious_bite, n = lag_val),
          date_lagged = dplyr::lag(date, n = lag_val)
        )
    } else {
      prob_lagged <- prob_data %>% 
        mutate(
          prob_lagged = dplyr::lead(prob_infectious_bite, n = abs(lag_val)),
          date_lagged = dplyr::lead(date, n = abs(lag_val))
        )
    }
    
    # Get start date minus burnin 
    start_date_pbite <- start_date - burnints
    prob_filtered <- prob_lagged[prob_lagged$date_lagged >= start_date_pbite & # &#
                                   prob_lagged$date_lagged < end_date & 
                                   !is.na(prob_lagged$date_lagged),]
    
    # # Now keep only every 2 days and sum the probability over those 2 days 
    # prob_filtered2 <- prob_filtered %>%
    #   mutate(twodaygroup = ceiling(rownumber(date))) %>%
    #   group_by(twodaygroup) %>%
    #   summarize(prob_infectious_bite = sum(prob_infectious_bite)) %>%
    #   left_join(prob_filtered)
    
    # instead of median, am now filtering to start date - burnin above
    # c(rep(median(prob_filtered$prob_lagged, na.rm = TRUE), burnints), # this is to have a probability of bite before the burnin
    #   prob_filtered$prob_lagged)
    return(prob_filtered)
  })
  
  names(lag_list) <- paste0("lag_", lags)
  return(lag_list)
}

# Script to calculate the SMC kill vector and the SMC kill timings and subset according to an infection start day 
#' @param smc_dose_days vector of days on which SMC is delivered per person, in vector format (e.g. list(c(1,30,60)))
#' @param ts default is a year in 2 day timesteps; this determines the length of the smc kill vector and is how many ts the model will run for
#' @param burnin length of burnin in 1 day timesteps; default is 0 which means that we need no padding (not run within cohort sim so the t means t) 
#' @param max_SMC_kill_rate scalar of the max kill rate in Weibull survival curve
#' @param lambda Weibull parameter for SMC kill curve
#' @param kappa Weibull parameter for SMC kill curve
#' @param infection_start_day within cohort, this is the day relative to the start of the sim (so t+burnin); for individual runs, this will subset the SMC kill vector from this number to the end of the vector
# output vector will be length (ts+burnin) / 2
get_smc_vectors <- function(smc_dose_days, 
                            ts = ceiling(365/2), 
                            burnin = 0, 
                            max_SMC_kill_rate, 
                            lambda,
                            kappa#, 
                            # infection_start_day
){
  
  time_since_smc <- lapply(list(smc_dose_days), 
                           calc_time_since_dose, 
                           days = 0:(ts*2 -1))
  smckillvec <- lapply(time_since_smc,
                       function(.){
                         kill <- max_SMC_kill_rate * exp(-(./ lambda)^kappa)  # calculate kill rate with hill function
                         kill[is.na(kill)] <- 0                                       # change NAs (when there is no SMC) to 0 
                         kill <- c(rep(0, burnin), kill)                              # pad front of vector with 0s for burnin period so that SMC begins after burnin
                         # Sum every two consecutive days 
                         n <- length(kill)
                         if (n %% 2 == 1) kill <- c(kill, kill[n])                         # pad with last value if odd length
                         colSums(matrix(kill, nrow = 2))
                       })
  return(unlist(smckillvec))
}
