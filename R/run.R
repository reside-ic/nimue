#' Run the vaccine model
#'
#' @inheritParams get_nimue_parameters
#'
#' @param ... nimue parameters
#'
#' @return Simulation output
#' @export
run <- function(...) {

  args <- get_nimue_parameters(...)

  # put arguments in current environment
  list2env(args, environment())

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
                     seeding_cases = seeding_cases,
                     seeding_age_order = seeding_age_order,
                     prob_hosp = prob_hosp,
                     prob_severe = prob_severe,
                     prob_non_severe_death_treatment = prob_non_severe_death_treatment,
                     prob_non_severe_death_no_treatment = prob_non_severe_death_no_treatment,
                     prob_severe_death_treatment = prob_severe_death_treatment,
                     prob_severe_death_no_treatment = prob_severe_death_no_treatment,
                     p_dist = p_dist,
                     rel_infectiousness = rel_infectiousness,
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
                     dur_V = dur_V,
                     vaccine_efficacy_infection = vaccine_efficacy_infection,
                     vaccine_efficacy_disease = vaccine_efficacy_disease,
                     max_vaccine = max_vaccine,
                     tt_vaccine = tt_vaccine ,
                     dur_vaccine_delay = dur_vaccine_delay,
                     vaccine_coverage_mat = vaccine_coverage_mat)

  # Set model type
  replicates <- 1
  mod_gen = vaccine

  # Running the Model
  mod <- mod_gen(user = pars, unused_user_action = "ignore")

  # Daily output by default
  t <- round(seq(from = 1, to = time_period))

  if(rk){
    results <- mod$run(t, replicate = replicates, method = "rk4", hini = 0.05)
  } else {
    results <- mod$run(t, replicate = replicates)
  }

  # coerce to array
  results <- array(results, dim = c(dim(results), 1), dimnames = dimnames(results))

  # Summarise inputs
  parameters <- args
  parameters$population <- pars$population
  parameters$hosp_bed_capacity <- pars$hosp_beds
  parameters$ICU_bed_capacity <- pars$ICU_beds
  parameters$beta_set <- pars$beta_set
  parameters$seeding_cases <- pars$E1_0
  parameters$contact_matrix_set <- pars$contact_matrix_set

  out <- list(output = results, parameters = parameters, model = mod, odin_parameters = pars)
  out <- structure(out, class = "nimue_simulation")
  return(out)

}
