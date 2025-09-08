# pak::pak("mrc-ide/site")
library(site)
library(countrycode)
library(tidyr)
library(terra)
library(tidyverse)
library(malariasimulation)
library(orderly2)

orderly_strict_mode()

orderly_resource(files = c('shapefiles/',
                           'chirps/'))

orderly_artefact(description = 'two datasets (Burkina Faso and Mali) with vectors of probability of infectious bite during the time of the trial',
                 files = c('prob_bite_BFA.rds',
                           'prob_bite_MLI.rds'))


# Get the site files
isos <- countrycode(sourcevar = c('Burkina Faso','Mali'),
                    origin = 'country.name',
                    destination = 'iso3c')

bf <- fetch_site(iso3c = isos[1], admin_level = 1)
hounderain <- bf$seasonality$fourier_prediction[bf$seasonality$fourier_prediction$name_1 == 'Haut-Bassins',]
# plot(hounderain$profile)

hounde <- subset_site(
  site = bf,
  site_filter = data.frame(
    country = "Burkina Faso",
    iso3c = "BFA",
    name_1 = "Haut-Bassins")
)

mali <- fetch_site(iso3c = isos[2], admin_level = 1)
# mali$sites
bougounirain <- mali$seasonality$fourier_prediction[mali$seasonality$fourier_prediction$name_1 == 'Sikasso',]
# plot(bougounirain$profile)

bougouni <- subset_site(
  site = mali,
  site_filter = data.frame(
    country = "Mali",
    iso3c = "MLI",
    name_1 = "Sikasso")
)

# umbrella
# devtools::install_github("mrc-ide/umbrella")
library(umbrella)

# urls_2017 <- get_urls(2017)
# urls_2018 <- get_urls(2018)
# urls_2019 <- get_urls(2019)
# urls_2020 <- get_urls(2020)
# urls <- c(urls_2017, urls_2018, urls_2019, urls_2020)
# # Create a vector of properly formatted dates
# start_date <- as.Date("2017-01-01")  # April 1st, 2017
# end_date <- as.Date("2020-12-31")    # December 31st, 2017
# dates <- seq(start_date, end_date, by = "day")
# date_strings <- format(dates, "%Y.%m.%d")
# date_strings <- date_strings[date_strings!='2020.02.29']

# Code below only needs to be run once -- it is to download the chirps files 
# daily_rainfall <- mapply(urls,
#                          FUN = download_raster,
#                          destination_file = paste0('chirps/chirps', date_strings, '.tif.gz'))
file_list <- list.files("chirps/", full.names = TRUE)

# Or read multiple files into a list
daily_rainfall <- terra::rast(file_list[-1462])

bf_sf <- read_sf("shapefiles/gadm41_BFA_1.shp") %>%
  dplyr::filter(NAME_1 == "Haut-Bassins")
mali_sf <- read_sf("shapefiles/gadm41_MLI_1.shp") %>%
  dplyr::filter(NAME_1 == "Sikasso")
# plot(bf_sf)
# plot(mali_sf)
extract_data_bf <- terra::extract(daily_rainfall, bf_sf, fun = mean)
# extract_data_bf
extract_data_mali <- terra::extract(daily_rainfall, mali_sf, fun = mean)
# extract_data_mali

# Tidy output 
raindata_bf <- extract_data_bf %>%
  pivot_longer(cols = starts_with('chirp'), 
               names_to = 'date',
               names_pattern = "chirps(.*)",
               values_to = 'rainfall') %>%
  select(-ID) %>%
  mutate(date = as.Date(gsub("\\.", "-", date))) %>%
  filter(date != as.Date('2020-02-29')) %>% # remove this date so the 4 years each have 365 days exactly 
  mutate(
    country = 'Mali',
    year = lubridate::year(date),
    day = row_number())
# ggplot(raindata_bf) + geom_point(aes(x = date, rainfall)) + scale_x_date(breaks = '1 month') + 
#   theme(axis.text.x = element_text(angle = 90))

raindata_mali <- extract_data_mali %>%
  pivot_longer(cols = starts_with('chirp'), 
               names_to = 'date',
               names_pattern = "chirps(.*)",
               values_to = 'rainfall') %>%
  select(-ID) %>%
  mutate(date = as.Date(gsub("\\.", "-", date))) %>%
  filter(date != as.Date('2020-02-29')) %>% # remove this date so the 4 years each have 365 days exactly 
  mutate(
    country = 'Mali',
    year = lubridate::year(date),
    day = row_number())
# plot(raindata_mali$date, raindata_mali$rainfall)

# Function to convert output of fitted fourier params to format that works with malsim
convert_coefs <- function(df, country, name, iso){
  coefs <- as.data.frame(t(df$coefficients)) %>%
    mutate(country = country,
           name_1 = name,
           iso3c = iso)
  return(coefs)
}

# Fit fourier parameters to BF rainfall 
fit1bf <- fit_fourier(rainfall = raindata_bf[raindata_bf$year == 2017,]$rainfall, 
                      t = raindata_bf[raindata_bf$year == 2017,]$day, 
                      floor = 0.1)
predictb_1 <- fourier_predict(coef = fit1bf$coefficients, t = 1:365, floor = 0.1)

fit2bf <- fit_fourier(rainfall = raindata_bf[raindata_bf$year == 2018,]$rainfall, 
                      t = raindata_bf[raindata_bf$year == 2018,]$day - 365, 
                      floor = 0.1)
predictb_2 <- fourier_predict(coef = fit2bf$coefficients, t = 1:365, floor = 0.1)

fit3bf <- fit_fourier(rainfall = raindata_bf[raindata_bf$year == 2019,]$rainfall, 
                      t = raindata_bf[raindata_bf$year == 2019,]$day - 730, 
                      floor = 0.1)
predictb_3 <- fourier_predict(coef = fit3bf$coefficients, t = 1:365, floor = 0.1)

fit4bf <- fit_fourier(rainfall = raindata_bf[raindata_bf$year == 2020,]$rainfall, 
                      t = raindata_bf[raindata_bf$year == 2020,]$day - 1095, 
                      floor = 0.1)
predictb_4 <- fourier_predict(coef = fit4bf$coefficients, t = 1:365, floor = 0.1)

bf_fits_p <- ggplot() + 
  geom_point(data = raindata_bf,
             aes(x = day, y = rainfall),
             alpha = 0.4, color = 'chartreuse3') + 
  geom_line(data = predictb_1,
            aes(x = t, y = profile),
            color = 'deeppink', linewidth = 1) + 
  geom_line(data = predictb_2,
            aes(x = t + 365, y = profile),
            color = 'dodgerblue', linetype = 2, linewidth = 1) + 
  geom_line(data = predictb_3,
            aes(x = t + 730, y = profile), 
            color = 'green4', linetype = 2, linewidth = 1) + 
  geom_line(data = predictb_4,
            aes(x = t + 1095, y = profile), 
            color = 'darkorchid3', linetype = 2, linewidth = 1) + 
  theme_bw()
ggsave('BFA_rainfallfits.png', bf_fits_p)

bffits <- list(fit1bf, fit2bf, fit3bf, fit4bf)
bfpredicts <- bind_rows(list(predictb_1, predictb_2, predictb_3, predictb_4))%>%
  mutate(day = 1:(365*4))
saveRDS(bfpredicts, 'bf_predictions.rds')

allcoefsbf <- lapply(bffits, convert_coefs, 
                     country = 'Burkina Faso',
                     name = 'Haut-Bassins',
                     iso = 'BFA')

# Fit fourier parameters to Mali rainfall 
fit1mali <- fit_fourier(rainfall = raindata_mali[raindata_mali$year == 2017,]$rainfall, 
                        t = raindata_mali[raindata_mali$year == 2017,]$day, 
                        floor = 0.1)
predictm_1 <- fourier_predict(coef = fit1mali$coefficients, t = 1:365, floor = 00.1)

fit2mali <- fit_fourier(rainfall = raindata_mali[raindata_mali$year == 2018,]$rainfall, 
                        t = raindata_mali[raindata_mali$year == 2018,]$day - 365, 
                        floor = 0.1)
predictm_2 <- fourier_predict(coef = fit2mali$coefficients, t = 1:365, floor = 0.1)

fit3mali <- fit_fourier(rainfall = raindata_mali[raindata_mali$year == 2019,]$rainfall, 
                        t = raindata_mali[raindata_mali$year == 2019,]$day - 730, 
                        floor = 0.1)
predictm_3 <- fourier_predict(coef = fit3mali$coefficients, t = 1:365, floor = 0.1)

fit4mali <- fit_fourier(rainfall = raindata_mali[raindata_mali$year == 2020,]$rainfall, 
                        t = raindata_mali[raindata_mali$year == 2020,]$day - 1095, 
                        floor = 0.1)
predictm_4 <- fourier_predict(coef = fit4mali$coefficients, t = 1:365, floor = 0.1)

mali_fits_p <- ggplot() + 
  geom_point(data = raindata_mali,
             aes(x = day, y = rainfall),
             alpha = 0.4, color = 'chartreuse3') + 
  geom_line(data = predictm_1,
            aes(x = t, y = profile),
            color = 'deeppink', linewidth = 1) + 
  geom_line(data = predictm_2,
            aes(x = t + 365, y = profile),
            color = 'dodgerblue', linetype = 2, linewidth = 1) + 
  geom_line(data = predictm_3,
            aes(x = t + 730, y = profile), 
            color = 'green4', linetype = 2, linewidth = 1) + 
  geom_line(data = predictm_4,
            aes(x = t + 1095, y = profile), 
            color = 'darkorchid3', linetype = 2, linewidth = 1) + 
  geom_vline(aes(xintercept = 365 * c(1,2,3,4))) + 
  theme_bw()
ggsave('MLI_rainfallfits.png', bf_fits_p)

malifits <- list(fit1mali, fit2mali, fit3mali, fit4mali)
malipredicts <- bind_rows(list(predictm_1, predictm_2, predictm_3, predictm_4)) %>%
  mutate(day = 1:(365*4))
saveRDS(malipredicts, 'mali_predictions.rds')

allcoefsmali <- lapply(malifits, convert_coefs, 
                       country = 'Mali',
                       name = 'Sikasso',
                       iso = 'MLI')

# Now, need to run malariasimilation or the odin model to get a curve of EIR or prevalence

# First, specify site parameters along with year-specific rainfall parameters for each of the years
bf_pars <- lapply(allcoefsbf, 
                  function(x) {
                    pars <- site_parameters(
                      interventions = bougouni$interventions,
                      demography = bougouni$demography,
                      vectors = bougouni$vectors$vector_species,
                      seasonality = x,
                      eir = bougouni$eir$eir,
                      overrides = list(
                        human_population = 10000
                      ))
                    
                    pars <- set_equilibrium(pars, 
                                            init_EIR = pars$init_EIR)
                    
                    return(pars)} )
mali_pars <- lapply(allcoefsmali, 
                    function(x) {
                      pars <- site_parameters(
                        interventions = bougouni$interventions,
                        demography = bougouni$demography,
                        vectors = bougouni$vectors$vector_species,
                        seasonality = x,
                        eir = bougouni$eir$eir,
                        overrides = list(
                          human_population = 10000
                        ))
                      
                      pars <- set_equilibrium(pars, 
                                              init_EIR = pars$init_EIR)} )

# now, run the simulation for each of these parameter sets 
outputs_bf <- lapply(bf_pars, 
                     function(x){
                       malariasimulation::run_simulation(
                         timesteps = x$timesteps,
                         parameters = x
                       )
                     })

outputs_mali <- lapply(mali_pars, 
                       function(x){
                         malariasimulation::run_simulation(
                           timesteps = x$timesteps,
                           parameters = x
                         )
                       })

# Process each of the outputs 
start_year = 2017
# for each of the fittet fourier parameter sets above, they are fitted to a single years worth of rainfall from that region 
# then, i run the simulation with site files for all of the years int eh site fiels (2000-2025) with those parameters 
# then i will need to filter out the years so that the year equals the rainfall year 
# create a date variable and filter to dates in the trial 
# create probabiltyof infectious bite variable to use for input to cohort simulation 

outputs_bf_all <- lapply(seq_along(outputs_bf), function(i) {
  outputs_bf[[i]] %>%
    mutate(rainfall_year = start_year + i - 1, # this is the fitted fourier parameters for that year
           year = floor(timestep / 365) + 2000,# this is the year of simulation within the site files 
           date = as.Date(timestep - (365*17) + 1, origin = '2017-01-01')) 
}) %>%
  bind_rows() %>%
  filter(year == rainfall_year# & # keep only the years for which the rainfall from that year was used to predict
           # (date >= '2017-04-01' & date < '2020-04-01') # filter to study period 
  ) %>%
  # calculate probability of bite per day per person 
  mutate(prob_infectious_bite = ifelse(!is.na(n_bitten), 
                                       n_bitten/(n_age_0_1824 + n_age_1825_5474 + n_age_5475_36499), 0))


outputs_mali_all <- lapply(seq_along(outputs_mali), function(i) {
  outputs_mali[[i]] %>%
    mutate(rainfall_year = start_year + i - 1, # this is the fitted fourier parameters for that year
           year = floor(timestep / 365) + 2000,# this is the year of simulation within the site files 
           date = as.Date(timestep - (365*17) + 1, origin = '2017-01-01')) 
}) %>%
  bind_rows() %>%
  filter(year == rainfall_year #&
           # (date >= '2017-04-01' & date < '2020-04-01')
         ) %>%
  # calculate probability of bite per day per person 
  mutate(prob_infectious_bite = ifelse(!is.na(n_bitten), 
                                       n_bitten/(n_age_0_1824 + n_age_1825_5474 + n_age_5475_36499), 0))




# Plot probability of bite for each of the sites 
# ggplot(outputs_bf_all) + 
#   geom_line(aes(x = timestep/365, y = prob_infectious_bite)) + 
#   labs(y = 'N bitten / N',
#        x = "Year") 
# 
# ggplot(outputs_mali_all) + 
#   geom_line(aes(x = timestep/365, y = prob_infectious_bite)) + 
#   labs(y = 'N bitten / N',
#        x = "Year") 

# output these proabilities to save for the cohort simulation 
saveRDS(outputs_bf_all %>% select(prob_infectious_bite, date, year, rainfall_year), 
        file = 'prob_bite_BFA.rds')
saveRDS(outputs_mali_all %>% select(prob_infectious_bite, date, year, rainfall_year), 
        file = 'prob_bite_MLI.rds')

# then can add in the interventions and see how the seasonality matches those cases 
