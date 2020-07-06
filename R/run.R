#' Run the vaccine model
#'
#' @param population Population vector (for each age group). Default = NULL,
#'   which will cause population to be sourced from \code{country}
#' @param country Character for country beign simulated. WIll be used to
#'   generate \code{population} and \code{contact_matrix_set} if
#'   unprovided. Either \code{country} or \code{population} and
#'   \code{contact_matrix_set} must be provided.
#' @param contact_matrix_set Contact matrices used in simulation. Default =
#'   NULL, which will generate this based on the \code{country}.
#' @param tt_contact_matrix Time change points for matrix change. Default = 0
#' @param R0 Basic Reproduction Number. Default = 3
#' @param tt_R0 Change time points for R0. Default = 0
#' @param beta_set Alternative parameterisation via beta rather than R0.
#'   Default = NULL, which causes beta to be estimated from R0
#' @param time_period Length of simulation. Default = 365
#' @param dt Time Step. Default = 0.1
#' @param replicates  Number of replicates. Default = 10
#' @param seeding_cases Initial number of cases seeding the epidemic
#' @param seed Random seed used for simulations. Deafult = runif(1, 0, 10000)
#' @param init Data.frame of initial conditions. Default = NULL
#' @param prob_hosp probability of hospitalisation by age.
#'   Default = c(0.001127564, 0.000960857, 0.001774408, 0.003628171,
#'   0.008100662, 0.015590734, 0.024597885, 0.035377529,
#'   0.04385549, 0.058495518, 0.08747709, 0.109730508,
#'   0.153943118, 0.177242143, 0.221362219, 0.267628264)
#' @param prob_severe Probability of developing severe symptoms by age.
#'   Default = c(3.73755e-05, 3.18497e-05, 5.88166e-05, 0.000120264,
#'   0.000268514, 0.000516788, 0.00081535, 0.001242525,
#'   0.001729275, 0.002880196, 0.00598205, 0.010821894,
#'   0.022736324, 0.035911156, 0.056362032, 0.081467057)
#' @param prob_non_severe_death_treatment Probability of death from non severe
#'   treated infection.
#'   Default = c(0.0125702, 0.0125702, 0.0125702, 0.0125702,
#'   0.0125702, 0.0125702, 0.0125702, 0.013361147,
#'   0.015104687, 0.019164124, 0.027477519, 0.041762108,
#'   0.068531658, 0.105302319, 0.149305732, 0.20349534)
#' @param prob_severe_death_treatment Probability of death from severe infection
#'   that is treated. Default = rep(0.5, 16)
#' @param prob_non_severe_death_no_treatment Probability of death in non severe
#'   hospital inections that aren't treated
#' @param prob_severe_death_no_treatment Probability of death from severe infection
#'   that is not treated. Default = rep(0.95, 16)
#' @param p_dist Preferentiality of age group receiving treatment relative to
#'   other age groups when demand exceeds healthcare capacity.
#' @param dur_E Mean duration of incubation period (days). Default = 4.6
#' @param dur_IMild Mean duration of mild infection (days). Default = 2.1
#' @param dur_ICase Mean duration from symptom onset to hospitil admission (days).
#'   Default = 4.5
#' @param dur_get_ox_survive Mean duration of oxygen given survive. Default = 5
#' @param dur_get_ox_die Mean duration of oxygen given death. Default = 5
#' @param dur_not_get_ox_survive Mean duration without oxygen given survive.
#'   Default = 5
#' @param dur_not_get_ox_die Mean duration without  oxygen given death.
#'  Default = 5
#' @param dur_get_mv_survive Mean duration of ventilation given survive.
#'   Default = 7.3
#' @param dur_get_mv_die Mean duration of ventilation given death. Default = 6
#' @param dur_not_get_mv_survive Mean duration without ventilation given
#'   survive. Default = 7.3
#' @param dur_not_get_mv_die Mean duration without ventilation given
#'   death. Default = 1
#' @param dur_rec Duration of recovery after coming off ventilation. Default = 2
#' @param hosp_bed_capacity General bed capacity. Can be single number of vector if capacity time-varies.
#' @param ICU_bed_capacity ICU bed capacity. Can be single number of vector if capacity time-varies.
#' @param tt_hosp_beds Times at which hospital bed capacity changes (Default = 0 = doesn't change)
#' @param tt_ICU_beds Times at which ICU bed capacity changes (Default = 0 = doesn't change)
#' @param seeding_cases Initial number of cases seeding the epidemic
#' @param dur_R Mean duration of naturally acquired immunity (days)
#' @param vaccination_target Index of age group targets for vaccination. Must be 0
#' (not vaccinated) or 1 (vaccinated) for each age group.
#' @param dur_V Mean duration of vaccine-derived immunity (days)
#' @param vaccine_efficacy_infection Efficacy of vaccine against infection (by age).
#' An efficacy of 1 will reduce FOI by 100 percent, an efficacy of 0.2 will reduce FOI by 20 percent etc.
#' @param vaccine_efficacy_disease Efficacy of vaccine against severe (requiring hospitilisation) disease (by age).
#' An efficacy of 1 will reduce the probability of hospitalisation by 100 percent,
#' an efficacy of 0.2 will reduce the probability of hospitalisation by 20 percent etc.
#' @param max_vaccine The maximum number of individuals who can be vaccinated per day.
#' @param tt_vaccine Time change points for vaccine capacity (\code{max_vaccine}).
#' @param dur_vaccine_delay Mean duration of period from vaccination to vaccine protection.
#' @param framework Model framework to run: stochastic or deterministic.
#'
#' @return Simulation output
#' @export
run <- function(

  # demography
  country = NULL,
  population = NULL,
  tt_contact_matrix = 0,
  contact_matrix_set = NULL,

  # transmission
  R0 = 3,
  tt_R0 = 0,
  beta_set = NULL,

  # initial state, duration, reps
  time_period = 365,
  dt = 0.1,
  replicates = 10,
  init = NULL,
  seed = stats::runif(1, 0, 100000000),

  # parameters
  # probabilities
  prob_hosp = probs$prob_hosp,
  prob_severe = probs$prob_severe,
  prob_non_severe_death_treatment = probs$prob_non_severe_death_treatment,
  prob_non_severe_death_no_treatment = probs$prob_non_severe_death_no_treatment,
  prob_severe_death_treatment = probs$prob_severe_death_treatment,
  prob_severe_death_no_treatment = probs$prob_severe_death_no_treatment,
  p_dist = probs$p_dist,

  # durations
  dur_E  = 4.6,
  dur_IMild = 2.1,
  dur_ICase = 4.5,

  dur_get_ox_survive = 9.5,
  dur_get_ox_die = 7.6,
  dur_not_get_ox_survive = 9.5*0.5,
  dur_not_get_ox_die = 7.6*0.5,

  dur_get_mv_survive = 11.3,
  dur_get_mv_die = 10.1,
  dur_not_get_mv_survive = 11.3*0.5,
  dur_not_get_mv_die = 1,

  dur_rec = 3.4,

  # vaccine
  dur_R = vaccine_pars$dur_R,
  vaccination_target = vaccine_pars$vaccination_target,
  dur_V = vaccine_pars$dur_V,
  vaccine_efficacy_infection = vaccine_pars$vaccine_efficacy_infection,
  vaccine_efficacy_disease = vaccine_pars$vaccine_efficacy_disease,
  max_vaccine = vaccine_pars$max_vaccine,
  tt_vaccine = vaccine_pars$tt_vaccine,
  dur_vaccine_delay = vaccine_pars$dur_vaccine_delay,

  # health system capacity
  hosp_bed_capacity = NULL,
  ICU_bed_capacity = NULL,
  tt_hosp_beds = 0,
  tt_ICU_beds = 0,

  seeding_cases = NULL,
  framework = "deterministic"
) {

  # Deal with framework shortcuts
  if(framework == "d"){
    framework <- "deterministic"
  }
  if(framework == "s"){
    framework <- "stochastic"
  }

  # Grab function arguments
  args <- as.list(environment())
  set.seed(seed)

  # create parameter list
  pars <- parameters(country = country,
                             population = population,
                             tt_contact_matrix = tt_contact_matrix ,
                             contact_matrix_set = contact_matrix_set,
                             R0 = R0,
                             tt_R0 = tt_R0 ,
                             beta_set = beta_set,
                             time_period = time_period,
                             dt = dt,
                             init = init,
                             seeding_cases = seeding_cases,
                             prob_hosp = prob_hosp,
                             prob_severe = prob_severe,
                             prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                             prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                             prob_severe_death_treatment = prob_severe_death_treatment,
                             prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                             p_dist = p_dist,
                             dur_E = dur_E,
                             dur_IMild = dur_IMild,
                             dur_ICase = dur_ICase,
                             dur_get_ox_survive = dur_get_ox_survive,
                             dur_get_ox_die = dur_get_ox_die,
                             dur_not_get_ox_survive = dur_not_get_ox_survive,
                             dur_not_get_ox_die = dur_not_get_ox_die,
                             dur_get_mv_survive = dur_get_mv_survive,
                             dur_get_mv_die = dur_get_mv_die,
                             dur_not_get_mv_survive = dur_not_get_mv_survive,
                             dur_not_get_mv_die = dur_not_get_mv_die,
                             dur_rec = dur_rec,
                             dur_R = dur_R,
                             hosp_bed_capacity = hosp_bed_capacity,
                             ICU_bed_capacity = ICU_bed_capacity,
                             tt_hosp_beds = tt_hosp_beds,
                             tt_ICU_beds = tt_ICU_beds,
                             vaccination_target = vaccination_target,
                             dur_V = dur_V,
                             vaccine_efficacy_infection = vaccine_efficacy_infection,
                             vaccine_efficacy_disease = vaccine_efficacy_disease,
                             max_vaccine = max_vaccine,
                             tt_vaccine = tt_vaccine ,
                             dur_vaccine_delay = dur_vaccine_delay,
                             framework = framework)

  # Set model type
  if(framework == "deterministic"){
    replicates <- 1
    mod_gen = vaccine_deterministic
  } else {
    mod_gen = vaccine_stochastic
  }

  # Running the Model
  mod <- mod_gen(user = pars, unused_user_action = "ignore")
  if(framework == "deterministic"){
    t <- seq(from = dt, to = time_period, by = dt)
  } else {
    t <- round(seq(from = 1, to = time_period/dt))
  }
  results <- mod$run(t, replicate = replicates)

  # coerce to array
  results <- array(results, dim = c(dim(results),1), dimnames = dimnames(results))

  # Summarise inputs
  parameters <- args
  parameters$population <- pars$population
  parameters$hosp_bed_capacity <- pars$hosp_beds
  parameters$ICU_bed_capacity <- pars$ICU_beds
  parameters$beta_set <- pars$beta_set
  parameters$seeding_cases <- pars$E1_0
  parameters$contact_matrix_set <- pars$contact_matrix_set

  out <- list(output = results, parameters = parameters, model = mod)
  out <- structure(out, class = "nimue_simulation")
  return(out)

}