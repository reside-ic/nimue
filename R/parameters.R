#' Return the default probabilities for modelling
#' @return list of default probabilities
default_probs <- function() {
  c(squire:::default_probs(), list(rel_infectiousness = rep(1, 17)))
}

probs <- default_probs()


#' Return the default vaccine parameters for modelling
#' @return list of default vaccine parameters
default_vaccine_pars <- function() {
  list(dur_R = Inf,
       dur_V = 365,
       vaccine_efficacy_infection = rep(0.95, 17),
       vaccine_efficacy_disease = rep(0.95, 17),
       max_vaccine = 1000,
       tt_vaccine = 0,
       dur_vaccine_delay = 14,
       vaccine_coverage_mat = matrix(0.8, ncol = 17, nrow = 1))
}

vaccine_pars <- default_vaccine_pars()

#' Vaccine parameters
#'
#' @details All durations are in days.
#'
#' @param ... nimue parameters
#'
#' @inheritParams get_nimue_parameters
#'
#' @export
parameters <- function(...) {

  list2env(list(...), environment())

  # Handle country population args
  cpm <- squire:::parse_country_population_mixing_matrix(country = country,
                                                         population = population,
                                                         contact_matrix_set = contact_matrix_set)
  country <- cpm$country
  population <- cpm$population
  contact_matrix_set <- cpm$contact_matrix_set

  # Standardise contact matrix set
  if(is.matrix(contact_matrix_set)){
    contact_matrix_set <- list(contact_matrix_set)
  }

  # populate contact matrix set if not provided
  if (length(contact_matrix_set) == 1) {
    baseline <- contact_matrix_set[[1]]
    contact_matrix_set <- vector("list", length(tt_contact_matrix))
    for(i in seq_along(tt_contact_matrix)) {
      contact_matrix_set[[i]] <- baseline
    }
  }


  # populate hospital and ICU bed capacity if not provided
  if (is.null(hosp_bed_capacity)) {
    if (!is.null(country)) {
      beds <- squire::get_healthcare_capacity(country)
      hosp_beds <- beds$hosp_beds
      hosp_bed_capacity <- rep(round(hosp_beds * sum(population)/1000), length(tt_hosp_beds))
    } else {
      hosp_bed_capacity <- round(5 * sum(population)/1000)
    }
  }
  if (is.null(ICU_bed_capacity)) {
    if (!is.null(country)) {
      beds <- squire::get_healthcare_capacity(country)
      ICU_beds <- beds$ICU_beds
      ICU_bed_capacity <- rep(round(ICU_beds * sum(population)/1000), length(tt_ICU_beds))
    } else {
      ICU_bed_capacity <- round(3 * hosp_bed_capacity/100)
    }
  }

  # Initial state and matrix formatting
  # ----------------------------------------------------------------------------

  # Initialise initial conditions
  mod_init <- init(population, seeding_cases, seeding_age_order)

  # Convert contact matrices to input matrices
  matrices_set <- squire:::matrix_set_explicit(contact_matrix_set, population)

  # If a vector is put in for matrix targeting
  if(is.vector(vaccine_coverage_mat)){
    vaccine_coverage_mat <- matrix(vaccine_coverage_mat, ncol = 17)
  }

  # Input checks
  # ----------------------------------------------------------------------------
  mc <- squire:::matrix_check(population[-1], contact_matrix_set)
  stopifnot(length(R0) == length(tt_R0))
  stopifnot(length(contact_matrix_set) == length(tt_contact_matrix))
  stopifnot(length(hosp_bed_capacity) == length(tt_hosp_beds))
  stopifnot(length(ICU_bed_capacity) == length(tt_ICU_beds))
  stopifnot(length(max_vaccine) == length(tt_vaccine))
  tc <- lapply(list(tt_R0/dt, tt_contact_matrix/dt), squire:::check_time_change, time_period/dt)
  tc2 <- lapply(list(tt_hosp_beds/dt, tt_ICU_beds/dt), squire:::check_time_change, time_period/dt)
  stopifnot(ncol(vaccine_coverage_mat) == 17)

  assert_pos(dt)
  assert_pos(dur_E)
  assert_pos(dur_IMild)
  assert_pos(dur_ICase)
  assert_pos(dur_get_ox_survive)
  assert_pos(dur_get_ox_die)
  assert_pos(dur_not_get_ox_survive)
  assert_pos(dur_not_get_ox_die)
  assert_pos(dur_get_mv_survive)
  assert_pos(dur_get_mv_die)
  assert_pos(dur_not_get_mv_survive)
  assert_pos(dur_not_get_mv_die)
  assert_pos(dur_R)
  assert_pos(dur_V)
  assert_pos(time_period)
  assert_pos(hosp_bed_capacity)
  assert_pos(ICU_bed_capacity)
  assert_pos(max_vaccine)
  assert_pos(dur_vaccine_delay)

  assert_length(prob_hosp, length(population))
  assert_length(prob_severe, length(population))
  assert_length(prob_non_severe_death_treatment, length(population))
  assert_length(prob_non_severe_death_no_treatment, length(population))
  assert_length(prob_severe_death_treatment, length(population))
  assert_length(prob_severe_death_no_treatment, length(population))
  assert_length(rel_infectiousness, length(population))
  assert_length(p_dist, length(population))

  assert_numeric(prob_hosp, length(population))
  assert_numeric(prob_severe, length(population))
  assert_numeric(prob_non_severe_death_treatment, length(population))
  assert_numeric(prob_non_severe_death_no_treatment, length(population))
  assert_numeric(prob_severe_death_treatment, length(population))
  assert_numeric(prob_severe_death_no_treatment, length(population))
  assert_numeric(rel_infectiousness, length(population))
  assert_numeric(p_dist, length(population))

  assert_leq(prob_hosp, 1)
  assert_leq(prob_severe, 1)
  assert_leq(prob_non_severe_death_treatment, 1)
  assert_leq(prob_non_severe_death_no_treatment, 1)
  assert_leq(prob_severe_death_treatment, 1)
  assert_leq(prob_severe_death_no_treatment, 1)
  assert_leq(rel_infectiousness, 1)
  assert_leq(p_dist, 1)

  assert_greq(prob_hosp, 0)
  assert_greq(prob_severe, 0)
  assert_greq(prob_non_severe_death_treatment, 0)
  assert_greq(prob_non_severe_death_no_treatment, 0)
  assert_greq(prob_severe_death_treatment, 0)
  assert_greq(prob_severe_death_no_treatment, 0)
  assert_greq(rel_infectiousness, 0)
  assert_greq(p_dist, 0)


  # Convert and Generate Parameters As Required
  # ----------------------------------------------------------------------------

  # durations
  gamma_E = 2 * 1/dur_E
  gamma_IMild = 1/dur_IMild
  gamma_ICase = 2 * 1/dur_ICase
  gamma_get_ox_survive = 2 * 1/dur_get_ox_survive
  gamma_get_ox_die = 2 * 1/dur_get_ox_die
  gamma_not_get_ox_survive = 2 * 1/dur_not_get_ox_survive
  gamma_not_get_ox_die = 2 * 1/dur_not_get_ox_die
  gamma_get_mv_survive = 2 * 1/dur_get_mv_survive
  gamma_get_mv_die = 2 * 1/dur_get_mv_die
  gamma_not_get_mv_survive = 2 * 1/dur_not_get_mv_survive
  gamma_not_get_mv_die = 2 * 1/dur_not_get_mv_die
  gamma_rec = 2 * 1/dur_rec
  gamma_R <- 2 * 1/dur_R
  gamma_V <- 2 * 1/dur_V
  gamma_vaccine_delay <- 2 * 1 / dur_vaccine_delay

  if (is.null(beta_set)) {
    baseline_matrix <- squire:::process_contact_matrix_scaled_age(contact_matrix_set[[1]], population)
    beta_set <- beta_est_infectiousness(dur_IMild = dur_IMild,
                                          dur_ICase = dur_ICase,
                                          prob_hosp = prob_hosp,
                                          mixing_matrix = baseline_matrix,
                                          rel_infectiousness = rel_infectiousness,
                                          R0 = R0)
  }

  # normalise to sum to 1
  p_dist <- matrix(rep(p_dist, 6), nrow = 17, ncol = 6)
  p_dist <- p_dist/mean(p_dist)

  # Format vaccine-specific parameters
  vaccine_efficacy_infection = 1 - vaccine_efficacy_infection
  prob_hosp_vaccine = (1 - vaccine_efficacy_disease) * prob_hosp

  # age X vaccine efficacy parameters
  prob_hosp <- matrix(c(prob_hosp, prob_hosp, prob_hosp,
                        prob_hosp_vaccine, prob_hosp_vaccine,
                        prob_hosp), nrow = 17, ncol = 6)
  vaccine_efficacy_infection <- matrix(c(rep(1, 17 * 3),
                                         vaccine_efficacy_infection, vaccine_efficacy_infection,
                                         rep(1, 17)), nrow = 17, ncol = 6)
  gamma_vaccine <- c(0, gamma_vaccine_delay, gamma_vaccine_delay, gamma_V, gamma_V, 0)


  # Collate Parameters Into List
  pars <- c(mod_init,
            list(N_age = length(population),
                 gamma_E = gamma_E,
                 gamma_IMild = gamma_IMild,
                 gamma_ICase = gamma_ICase,
                 gamma_get_ox_survive = gamma_get_ox_survive,
                 gamma_get_ox_die = gamma_get_ox_die,
                 gamma_not_get_ox_survive = gamma_not_get_ox_survive,
                 gamma_not_get_ox_die = gamma_not_get_ox_die,
                 gamma_get_mv_survive = gamma_get_mv_survive,
                 gamma_get_mv_die = gamma_get_mv_die,
                 gamma_not_get_mv_survive = gamma_not_get_mv_survive,
                 gamma_not_get_mv_die = gamma_not_get_mv_die,
                 gamma_rec = gamma_rec,
                 gamma_R = gamma_R,
                 prob_hosp = prob_hosp,
                 prob_severe = prob_severe,
                 prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                 prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                 prob_severe_death_treatment = prob_severe_death_treatment,
                 prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                 rel_infectiousness = rel_infectiousness,
                 p_dist = p_dist,
                 hosp_beds = hosp_bed_capacity,
                 ICU_beds = ICU_bed_capacity,
                 tt_hosp_beds = tt_hosp_beds,
                 tt_ICU_beds = tt_ICU_beds,
                 tt_matrix = tt_contact_matrix,
                 mix_mat_set = matrices_set,
                 tt_beta = tt_R0,
                 beta_set = beta_set,
                 dt = dt,
                 population = population,
                 contact_matrix_set = contact_matrix_set,
                 max_vaccine = max_vaccine,
                 tt_vaccine = tt_vaccine,
                 vaccine_efficacy_infection = vaccine_efficacy_infection,
                 vaccine_coverage_mat = vaccine_coverage_mat,
                 N_vaccine = 6,
                 N_prioritisation_steps = nrow(vaccine_coverage_mat),
                 gamma_vaccine = gamma_vaccine))

  class(pars) <- c("vaccine_parameters", "nimue_parameters")

  return(pars)
}

#' Estimate beta parameter for explicit model
#'
#' @param dur_IMild Duration of mild infectiousness (days)
#' @param dur_ICase Delay between symptom onset and requiring hospitalisation (days)
#' @param prob_hosp Probability of hospitilisation by ages
#' @param rel_infectiousness Relative infectiousness of age categories relative
#'   to maximum infectiousness age category
#' @param mixing_matrix Mixing matrix
#' @param R0 Basic reproduction number
#'
#' @return Beta parameter
#' @export
#'
# #' @examples
beta_est_infectiousness <- function(dur_IMild,
                                    dur_ICase,
                                    prob_hosp,
                                    rel_infectiousness,
                                    mixing_matrix,
                                    R0) {

  # assertions
  assert_single_pos(dur_ICase, zero_allowed = FALSE)
  assert_single_pos(dur_IMild, zero_allowed = FALSE)
  assert_numeric(prob_hosp)
  assert_numeric(rel_infectiousness)
  assert_same_length(prob_hosp, rel_infectiousness)
  assert_numeric(mixing_matrix)
  assert_square_matrix(mixing_matrix)
  assert_same_length(mixing_matrix[,1], prob_hosp)
  assert_pos(R0, zero_allowed = FALSE)

  if(sum(is.na(prob_hosp)) > 0) {
    stop("prob_hosp must not contain NAs")
  }

  if(sum(is.na(rel_infectiousness)) > 0) {
    stop("rel_infectiousness must not contain NAs")
  }

  if(sum(is.na(mixing_matrix)) > 0) {
    stop("mixing_matrix must not contain NAs")
  }

  relative_R0_by_age <- prob_hosp*dur_ICase + (1-prob_hosp)*dur_IMild
  adjusted_eigen <- Re(eigen(mixing_matrix*relative_R0_by_age*rel_infectiousness)$values[1])
  R0 / adjusted_eigen
}

#' Return the parameters for the vaccine model
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
#' @param prob_hosp probability of hospitalisation by age.
#'   Default = c(0.000744192, 0.000634166,0.001171109, 0.002394593, 0.005346437,
#'   0.010289885, 0.016234604, 0.023349169, 0.028944623, 0.038607042,
#'   0.057734879, 0.072422135, 0.101602458, 0.116979814, 0.146099064,
#'   0.176634654 ,0.180000000)
#' @param prob_severe Probability of developing severe symptoms by age.
#'   Default = c(0.05022296,	0.05022296,	0.05022296,	0.05022296,	0.05022296,
#'   0.05022296,	0.05022296,	0.053214942, 0.05974426,	0.074602879,
#'   0.103612417, 0.149427991, 0.223777304,	0.306985918,
#'   0.385779555, 0.461217861, 0.709444444)
#' @param prob_non_severe_death_treatment Probability of death from non severe
#'   treated infection.
#'   Default = c(0.0125702,	0.0125702,	0.0125702,	0.0125702,
#'   0.0125702,	0.0125702,	0.0125702,	0.013361147,
#'   0.015104687,	0.019164124,	0.027477519,	0.041762108,
#'   0.068531658,	0.105302319,	0.149305732,	0.20349534,	0.5804312)
#' @param prob_severe_death_treatment Probability of death from severe infection
#'   that is treated. Default = rep(0.5, 17)
#' @param prob_non_severe_death_no_treatment Probability of death in non severe
#'   hospital inections that aren't treated
#' @param prob_severe_death_no_treatment Probability of death from severe infection
#'   that is not treated. Default = rep(0.95, 17)
#' @param p_dist Preferentiality of age group receiving treatment relative to
#'   other age groups when demand exceeds healthcare capacity.
#' @param rel_infectiousness Relative infectiousness per age category relative
#'   to maximum infectiousness category. Default = rep(1, 17)
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
#' @param seeding_age_order Vector specifying the order in which seeds are allocated to ages.
#'   If NULL, seeds are distributed randomly within working ages. If specified, must be a vector
#'   of length 17 specifying the order seeds are allocated, e.g. 1:17 will allocate first seed
#'   to the youngest age group, then the second youngest and so on. Default = NULL
#' @param dur_R Mean duration of naturally acquired immunity (days)
#' @param dur_V Mean duration of vaccine-derived immunity (days)
#' @param vaccine_efficacy_infection Efficacy of vaccine against infection (by age).
#' An efficacy of 1 will reduce FOI by 100 percent, an efficacy of 0.2 will reduce FOI by 20 percent etc.
#' @param vaccine_efficacy_disease Efficacy of vaccine against severe (requiring hospitilisation) disease (by age).
#' An efficacy of 1 will reduce the probability of hospitalisation by 100 percent,
#' an efficacy of 0.2 will reduce the probability of hospitalisation by 20 percent etc.
#' @param max_vaccine The maximum number of individuals who can be vaccinated per day.
#' @param tt_vaccine Time change points for vaccine capacity (\code{max_vaccine}).
#' @param dur_vaccine_delay Mean duration of period from vaccination to vaccine protection.
#' @param vaccine_coverage_mat Vaccine coverage targets by age (columns) and priority (row)
#' @param rk Use rk4 solver with a fixed timestep
#'
#' @return A list of parameters for use in the model
#'
#' @export
get_nimue_parameters <- function(
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

  # onward infectiousness
  rel_infectiousness = probs$rel_infectiousness,

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
  dur_V = vaccine_pars$dur_V,
  vaccine_efficacy_infection = vaccine_pars$vaccine_efficacy_infection,
  vaccine_efficacy_disease = vaccine_pars$vaccine_efficacy_disease,
  max_vaccine = vaccine_pars$max_vaccine,
  tt_vaccine = vaccine_pars$tt_vaccine,
  dur_vaccine_delay = vaccine_pars$dur_vaccine_delay,
  vaccine_coverage_mat = vaccine_pars$vaccine_coverage_mat,

  # health system capacity
  hosp_bed_capacity = NULL,
  ICU_bed_capacity = NULL,
  tt_hosp_beds = 0,
  tt_ICU_beds = 0,

  seeding_cases = 20,
  seeding_age_order = NULL,
  rk = TRUE
  ) {
  as.list(environment())
}
