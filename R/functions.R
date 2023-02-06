#' Setup a default BioGeoBEARS run
#'
#' Equivalent to DEC model
#'
#' @param trfn Path to BioGeoBEARS tree file (in Newick format)
#' @param geogfn Path to BioGeoBEARS geography file (in phylip format)
#' @param max_range_size Maximum number of range states allowed in model
#' @param timesfn_fn Path to time slices file, optional
#' @param dispersal_multipliers_fn Path to dispersal multipliers file, optional
#' @param min_branchlength Minimum branchlength to treat tip as a direct
#' ancestor (no speciation event)
#' @param num_cores_to_use Number of cores to use during analysis
#' @param cluster_already_open Logical; is a cluster already open? Optional.
#' @param use_optimx Type of optimization to use; choose "optim", "optimx",
#' or "GenSA"
#'
#' @return List
#'
setup_bgb_default <- function(
  trfn,
  geogfn,
  max_range_size,
  timesfn = NULL,
  dispersal_multipliers_fn = NULL,
  min_branchlength =  0.000001,
  num_cores_to_use = 1,
  cluster_already_open = NULL,
  use_optimx = "GenSA") {

  # Basic settings
  # - default model is DEC
  bgb_obj <- BioGeoBEARS::define_BioGeoBEARS_run()
  bgb_obj$trfn <- trfn
  bgb_obj$geogfn <- geogfn
  bgb_obj$max_range_size <- max_range_size
  bgb_obj$min_branchlength <- min_branchlength # Min to treat tip as a direct ancestor (no speciation event) # nolint
  bgb_obj$include_null_range <- TRUE # set to FALSE for e.g. DEC* model, DEC*+J, etc. # nolint

  # (Optional) specify time slices file
  if (!is.null(timesfn)) bgb_obj$timesfn <- timesfn

  # (Optional) specify dispersal multipliers file
  if (!is.null(dispersal_multipliers_fn)) {
    bgb_obj$dispersal_multipliers_fn <- dispersal_multipliers_fn
  }

  # Speed options and multicore processing if desired
  bgb_obj$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check) #nolint
  bgb_obj$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params) #nolint
  bgb_obj$use_optimx <- use_optimx    # if FALSE, use optim() instead of optimx() #nolint
  bgb_obj$num_cores_to_use <- num_cores_to_use
  bgb_obj$force_sparse <- FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale #nolint
  bgb_obj$cluster_already_open <- cluster_already_open

  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work! #nolint
  # (It also runs some checks on these inputs for certain errors.)
  bgb_obj <- BioGeoBEARS::readfiles_BioGeoBEARS_run(
    bgb_obj)

  # Divide the tree up by timeperiods/strata for stratified analysis
  if (!is.null(timesfn)) {
    bgb_obj <- BioGeoBEARS::section_the_tree(
      inputs = bgb_obj,
      make_master_table = TRUE, plot_pieces = FALSE)
  }

  # Good default settings to get ancestral states
  bgb_obj$return_condlikes_table <- TRUE
  bgb_obj$calc_TTL_loglike_from_condlikes_table <- TRUE
  bgb_obj$calc_ancprobs <- TRUE    # get ancestral states from optim run #nolint

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(bgb_obj)

  bgb_obj
}


#' Setup a BioGeoBEARS run using the DEC model
#'
#' @param trfn Path to BioGeoBEARS tree file (in Newick format)
#' @param geogfn Path to BioGeoBEARS geography file (in phylip format)
#' @param max_range_size Maximum number of range states allowed in model
#' @param jump Logical; should "jump" parameter be used?
#' @param res_dec Output of DEC model; required if `jump` is TRUE
#' @param timesfn_fn Path to time slices file, optional
#' @param dispersal_multipliers_fn Path to dispersal multipliers file, optional
#' @param min_branchlength Minimum branchlength to treat tip as a direct
#' ancestor (no speciation event)
#' @param num_cores_to_use Number of cores to use during analysis
#' @param cluster_already_open Logical; is a cluster already open? Optional.
#' @param use_optimx Type of optimization to use; choose "optim", "optimx",
#' or "GenSA"
#'
#' @return List
#'
setup_bgb_dec <- function(
  trfn,
  geogfn,
  max_range_size,
  jump = FALSE,
  res_dec = NULL,
  timesfn = NULL,
  dispersal_multipliers_fn = NULL,
  min_branchlength =  0.000001,
  num_cores_to_use = 1,
  cluster_already_open = NULL,
  use_optimx = "GenSA") {

  # Setup default model (DEC)
  bgb_obj <- setup_bgb_default(
    trfn = trfn,
    geogfn = geogfn,
    max_range_size = max_range_size,
    timesfn = timesfn,
    dispersal_multipliers_fn = dispersal_multipliers_fn,
    min_branchlength = min_branchlength,
    num_cores_to_use = num_cores_to_use,
    cluster_already_open = cluster_already_open,
    use_optimx = use_optimx
  )

  # DEC+J # nolint
  if (jump == TRUE) {
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart <- res_dec$outputs@params_table["d", "est"]
    estart <- res_dec$outputs@params_table["e", "est"]
    jstart <- 0.0001

    # Input starting values for d, e
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "init"] <- dstart
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "est"] <- dstart
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "init"] <- estart
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "est"] <- estart
    # Add j as a free parameter
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart
  }

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(bgb_obj)

  bgb_obj

}

#' Setup a BioGeoBEARS run using the DIVALIKE model
#'
#' @param trfn Path to BioGeoBEARS tree file (in Newick format)
#' @param geogfn Path to BioGeoBEARS geography file (in phylip format)
#' @param max_range_size Maximum number of range states allowed in model
#' @param jump Logical; should "jump" parameter be used?
#' @param res_diva Output of DIVALIKE model; required if `jump` is TRUE
#' @param timesfn_fn Path to time slices file, optional
#' @param dispersal_multipliers_fn Path to dispersal multipliers file, optional
#' @param min_branchlength Minimum branchlength to treat tip as a direct
#' ancestor (no speciation event)
#' @param num_cores_to_use Number of cores to use during analysis
#' @param cluster_already_open Logical; is a cluster already open? Optional.
#' @param use_optimx Type of optimization to use; choose "optim", "optimx",
#' or "GenSA"
#'
#' @return List
#'
setup_bgb_diva <- function(
  trfn,
  geogfn,
  max_range_size,
  jump = FALSE,
  res_diva = NULL,
  timesfn = NULL,
  dispersal_multipliers_fn = NULL,
  min_branchlength = 0.000001,
  num_cores_to_use = 1,
  cluster_already_open = NULL,
  use_optimx = "GenSA") {

  # Setup default model (DEC)
  bgb_obj <- setup_bgb_default(
    trfn = trfn,
    geogfn = geogfn,
    max_range_size = max_range_size,
    timesfn = timesfn,
    dispersal_multipliers_fn = dispersal_multipliers_fn,
    min_branchlength = min_branchlength,
    num_cores_to_use = num_cores_to_use,
    cluster_already_open = cluster_already_open,
    use_optimx = use_optimx
  )

  if (jump == FALSE) {
    # Set up DIVALIKE model
    #
    # Remove subset-sympatry
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

    bgb_obj$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
    bgb_obj$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
    bgb_obj$BioGeoBEARS_model_object@params_table["y", "type"] <- "ysv*1/2"
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"

    # Allow classic, widespread vicariance; all events equiprobable
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01v", "init"] <- 0.5
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01v", "est"] <- 0.5
  } else if (jump == TRUE) {
    # Set up DIVALIKE+J model
    #
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart <- res_diva$outputs@params_table["d", "est"]
    estart <- res_diva$outputs@params_table["e", "est"]
    jstart <- 0.0001

    # Input starting values for d, e
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "init"] <- dstart
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "est"] <- dstart
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "init"] <- estart
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "est"] <- estart

    # Remove subset-sympatry
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

    bgb_obj$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
    bgb_obj$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
    bgb_obj$BioGeoBEARS_model_object@params_table["y", "type"] <- "ysv*1/2"
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"

    # Allow classic, widespread vicariance; all events equiprobable
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01v", "init"] <- 0.5
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01v", "est"] <- 0.5

    # Add jump dispersal/founder-event speciation
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart

    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in
    # DEC+J)
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "min"] <- 0.00001
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "max"] <- 1.99999
  }

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(bgb_obj)

  bgb_obj

}

#' Setup a BioGeoBEARS run using the BAYAREALIKE model
#'
#' @param trfn Path to BioGeoBEARS tree file (in Newick format)
#' @param geogfn Path to BioGeoBEARS geography file (in phylip format)
#' @param max_range_size Maximum number of range states allowed in model
#' @param jump Logical; should "jump" parameter be used?
#' @param res_bayarea Output of BAYAREALIKE model; required if `jump` is TRUE
#' @param timesfn_fn Path to time slices file, optional
#' @param dispersal_multipliers_fn Path to dispersal multipliers file, optional
#' @param min_branchlength Minimum branchlength to treat tip as a direct
#' ancestor (no speciation event)
#' @param num_cores_to_use Number of cores to use during analysis
#' @param cluster_already_open Logical; is a cluster already open? Optional.
#' @param use_optimx Type of optimization to use; choose "optim", "optimx",
#' or "GenSA"
#'
#' @return List
#'
setup_bgb_bayarea <- function(
  trfn,
  geogfn,
  max_range_size,
  jump = FALSE,
  res_bayarea = NULL,
  timesfn = NULL,
  dispersal_multipliers_fn = NULL,
  min_branchlength = 0.000001,
  num_cores_to_use = 1,
  cluster_already_open = NULL,
  use_optimx = "GenSA") {

  # Setup default model (DEC)
  bgb_obj <- setup_bgb_default(
    trfn = trfn,
    geogfn = geogfn,
    max_range_size = max_range_size,
    timesfn = timesfn,
    dispersal_multipliers_fn = dispersal_multipliers_fn,
    min_branchlength = min_branchlength,
    num_cores_to_use = num_cores_to_use,
    cluster_already_open = cluster_already_open,
    use_optimx = use_optimx
  )

  if (jump == FALSE) {

    # Set up BAYAREALIKE model
    # No subset sympatry
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

    # No vicariance
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "est"] <- 0.0

    # Adjust linkage between parameters
    bgb_obj$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
    bgb_obj$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/1"
    bgb_obj$BioGeoBEARS_model_object@params_table["y", "type"] <- "1-j"

    # Only sympatric/range-copying (y) events allowed, and with
    # exact copying (both descendants always the same size as the ancestor)
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01y", "init"] <- 0.9999
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01y", "est"] <- 0.9999

  } else if (jump == TRUE) {
    # Set up BAYAREALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart <- res_bayarea$outputs@params_table["d", "est"]
    estart <- res_bayarea$outputs@params_table["e", "est"]
    jstart <- 0.0001

    # Input starting values for d, e
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "init"] <- dstart
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "est"] <- dstart
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "init"] <- estart
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "est"] <- estart

    # No subset sympatry
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
    bgb_obj$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

    # No vicariance
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
    bgb_obj$BioGeoBEARS_model_object@params_table["v", "est"] <- 0.0

    # *DO* allow jump dispersal/founder-event speciation (set the starting value
    # close to 0)
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart

    # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in
    # DEC+J) or 2 (as in DIVALIKE+J)
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "max"] <- 0.99999

    # Adjust linkage between parameters
    bgb_obj$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
    bgb_obj$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/1"
    bgb_obj$BioGeoBEARS_model_object@params_table["y", "type"] <- "1-j"

    # Only sympatric/range-copying (y) events allowed, and with
    # exact copying (both descendants always the same size as the ancestor)
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01y", "init"] <- 0.9999
    bgb_obj$BioGeoBEARS_model_object@params_table["mx01y", "est"] <- 0.9999

    # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers,
    # usually Windows machines. I can't replicate this on my Mac machines, but
    # it is almost certainly just some precision under-run issue, when
    # optim/optimx tries some parameter value just below zero.  The "min" and
    # "max" options on each parameter are supposed to prevent this, but
    # apparently optim/optimx sometimes go slightly beyond these limits.
    # Anyway, if you get a crash, try raising "min" and lowering "max" slightly
    # for each parameter:
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "min"] <- 0.0000001
    bgb_obj$BioGeoBEARS_model_object@params_table["d", "max"] <- 4.9999999

    bgb_obj$BioGeoBEARS_model_object@params_table["e", "min"] <- 0.0000001
    bgb_obj$BioGeoBEARS_model_object@params_table["e", "max"] <- 4.9999999

    bgb_obj$BioGeoBEARS_model_object@params_table["j", "min"] <- 0.00001
    bgb_obj$BioGeoBEARS_model_object@params_table["j", "max"] <- 0.99999
  }

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(bgb_obj)

  bgb_obj

}

#' Run BioGeoBEARS::bears_optim_run() quietly
#'
#' Instead of printing all the results to the screen, which
#' it does by default
#'
#' @param bgb_obj List; a BioGeoBEARS object
#' @param ... Other arguments not used by this function, but meant for
#' tracking with drake.
#'
#' @return List
#'
run_bgb_quietly <- function(bgb_obj, ...) {

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(bgb_obj)

  # Suppress output printed to screen
  invisible(capture.output(results <- BioGeoBEARS::bears_optim_run(bgb_obj)))

  results

}

#' Run BioGeoBEARS::runBSM() quietly
#'
#' Instead of printing all the results to the screen, which
#' it does by default
#'
#' @param ... Arguments passed to BioGeoBEARS::runBSM()
#'
#' @return List
#'
run_bsm_quietly <- function(...) {

  # Suppress output printed to screen
  invisible(capture.output(results <- BioGeoBEARS::runBSM(...)))

  results

}

#' Compare BioGeoBEARS models using AIC
#'
#' @param dec_model Results of running BioGeoBEARS using the DEC model
#' @param dec_const_model Results of running BioGeoBEARS using the DEC
#' constrained model
#' @param decj_model Results of running BioGeoBEARS using the DECj model
#' @param decj_const_model Results of running BioGeoBEARS using the DECj
#' constrained model
#' @param diva_model Results of running BioGeoBEARS using the DIVALIKE model
#' @param diva_const_model Results of running BioGeoBEARS using the DIVALIKE
#' constrained model
#' @param divaj_model Results of running BioGeoBEARS using the DIVALIKE+J model
#' @param divaj_const_model Results of running BioGeoBEARS using the DIVALIKE+J
#' constrained model
#' @param bayarea_model Results of running BioGeoBEARS using the BAYAREALIKE
#' model
#' @param bayarea_const_model Results of running BioGeoBEARS using the
#' BAYAREALIKE constrained model
#' @param bayareaj_model Results of running BioGeoBEARS using the BAYAREALIKE+J
#' model
#' @param bayareaj_const_model Results of running BioGeoBEARS using the
#' BAYAREALIKE+J constrained model
#'
#' @return dataframe
get_bgb_stats <- function(
  dec_model,
  dec_const_model,
  decj_model,
  decj_const_model,
  diva_model,
  diva_const_model,
  divaj_model,
  divaj_const_model,
  bayarea_model,
  bayarea_const_model,
  bayareaj_model,
  bayareaj_const_model,
  include_model = FALSE) {

  models <- list(
      dec = dec_model,
      dec_const = dec_const_model,
      decj = decj_model,
      decj_const = decj_const_model,
      diva = diva_model,
      diva_const = diva_const_model,
      divaj = divaj_model,
      divaj_const = divaj_const_model,
      bayarea = bayarea_model,
      bayarea_const = bayarea_const_model,
      bayareaj = bayareaj_model,
      bayareaj_const = bayareaj_const_model
    )

tibble::tibble(mod_name = names(models), model = models) %>%
  mutate(
    params = map(
      models,
      ~extract_params_from_BioGeoBEARS_results_object(
        ., returnwhat = "table", addl_params = "j", paramsstr_digits = 4
      ))
  ) %>%
  unnest(params) %>%
  rename(lnl = LnL) %>%
  mutate(aic = 2 * numparams - 2 * lnl)
}

# Format the model AIC table for manuscript
format_model_table <- function(bgb_stats) {
  bgb_stats %>%
  # Arrange by constraint type, jump type, then model name
  mutate(
    constraint = if_else(
      str_detect(mod_name, "const"),
      "constrained",
      "unconstrained"),
    constraint = factor(constraint, levels = c("unconstrained", "constrained")),
    jump = if_else(
      str_detect(mod_name, "j"),
      "jump_yes",
      "jump_no"),
    jump = factor(jump, levels = c("jump_no", "jump_yes")),
    model_type = case_when(
      str_detect(mod_name, "dec") ~ "dec",
      str_detect(mod_name, "diva") ~ "diva",
      str_detect(mod_name, "bayarea") ~ "bayarea"
    ),
    model_type = factor(model_type, levels = c("dec", "diva", "bayarea"))
  ) %>%
  arrange(constraint, jump, model_type) %>%
  select(mod_name, constraint, lnl, numparams, d, e, j, aic)
}

# Get the best-scoring model from AIC table
get_best_model <- function(bgb_stats) {
  best_mod_row <-
  bgb_stats %>%
    slice_min(aic)
  message(paste("The best model by AIC was", best_mod_row$mod_name))
  best_mod_row$model[[1]]
}

#' Count the number of dispersal events to realms
#' 
#' Note that cladogenetic events do not include dispersal to new biomes,
#' so those don't get used in this function
#' 
#' See
#' https://github.com/idiv-biodiversity/2019_big_data_biogeography/blob/master/wed_DEC_bsm.Rmd
#'
#' @param bsm_results Results of stochastic character mapping to estimate range shifts over time 
#' @param max_age Maximum age to bin branches
#' @param bin_width Width of each bin
#' @param phy Phylogenetic tree
#'
#' @return Dataframe
#' 
count_normalized_dispersal_events <- function(bsm_results, bin_width, phy, max_age = NULL) {

  ana_events_tables <- bsm_results$RES_ana_events_tables

  # From Nick:
  # 
  # abs_event_time is absolute event time, this is usually what you want
  # event_time is IIRC the time along an individual branch
  # time_bp is the time-before-present of a node 
  # [JOEL: node subtending the event in the case of anagensis],
  # the same as you would see in the tree table from prt(tree).
  # The cladogenesis events always happen at the nodes.

  count_angen <- function (data) {
    data %>% 
      select(abs_event_time, event_type, event_txt, dispersal_to) %>% 
      mutate(from = str_split_fixed(event_txt, n = 2, pattern = "->")[,1 ])
  }

  # event_type is "e" for extinction, or "d" for dispersal
  angen <- map_df(ana_events_tables, count_angen, .id = "rep") %>%
    as_tibble()

  # Optionally filter by maximum age of dispersal event
  if(!is.null(max_age)) angen <- filter(angen, abs_event_time < max_age)

  # Count the number of dispersals per time bin
  disps <-
    angen %>%
    # FIXME: check if cladogenetic events need to be excluded
    filter(event_type == "d") %>%
    mutate(bin = cut_width(abs_event_time, width = bin_width, center = bin_width / 2)) %>%
    group_by(rep, dispersal_to, bin) %>%
    summarize(abs_count = n(),
              .groups = "drop") %>%
    group_by(dispersal_to, bin) %>%
    summarize(lower = quantile(abs_count, probs = 0.025),
              mean = mean(abs_count),
              upper = quantile(abs_count, probs = 0.975),
              sd = sd(abs_count),
              n = n(),
              .groups = "drop")

  # Bin branch lengths
  brls <- bin_brln(phy, bin_width, max_age)

  # Normalize the number of shifts by the total amount of branch length per time bin
  res <- left_join(disps, brls, by = "bin") %>% 
    mutate(lower = lower / tot_brl) %>% 
    mutate(mean = mean / tot_brl) %>% 
    mutate(upper = upper / tot_brl) %>%
    mutate(sd = sd / tot_brl)
}

bin_brln <- function(phy, bin_width, max_age = NULL) {
  # Get the total branch length per time bin (for normalizing)
  branchtimes <- ape::branching.times(phy)

  if(!is.null(max_age)) branchtimes <- branchtimes[branchtimes < max_age]

  # calculate maximum number (integer) of bins
  maxmilcatogories <-  ceiling(max(branchtimes) / bin_width)
  # make an empty matrix to fill with total branch length per bin
  milgaps_mod <- matrix(ncol = 2, nrow = maxmilcatogories, 0)
  colnames(milgaps_mod) <- c("time","tot_brl")
  # fill time (assuming bin width is a time in same units as tree, eg, million-years)
  milgaps_mod[,1] <- 1:maxmilcatogories*bin_width
  ## Sum times per bin
  v <- branchtimes
  for(i in 1:length(milgaps_mod[, 2])){
    s <- milgaps_mod[i, 1]        # time start
    e <- milgaps_mod[i, 1] - bin_width # time end
    l1 <- length(v[v > s]) * bin_width + min(length(v[v > s]), 1) * bin_width
    l2 <- sum(v[v > e & v<s] - e)
    
    if (l1 == 0){
      l2 = l2*2
    } # root
    
    milgaps_mod[i,2] <- l1 + l2
  }

  # Bin branch lengths
  milgaps_mod %>%
    as.data.frame %>% 
    as_tibble %>%
    # Shift time so it doesn't get cut exactly at the border
    mutate(time = time - 1) %>%
    mutate(bin = cut_width(time, width = bin_width, center = bin_width / 2))
}

# BSM ----

# Count cladogeneic events in a single rep
count_cladgen_single <- function(clado_event_table) {
  clado_event_table %>%
      select(
        node,
        event_type = clado_event_type,
        event_txt = clado_event_txt,
        dispersal_to = clado_dispersal_to,
        # since cladogenic events take place at nodes,
        # time_bp of the node is same as the abs event time
        abs_event_time = time_bp) %>%
    mutate(from = str_split_fixed(event_txt, n = 2, pattern = "->")[,1 ])
}

# Count cladogeneic events across reps
count_cladgen_across_reps <- function(clado_event_tables) {
  map_df(
  clado_event_tables,
  count_cladgen_single, .id = "rep") %>%
  as_tibble() %>%
  filter(event_txt != "")
}

# Count anagenic events in a single rep
count_angen_single <- function (ana_event_table) {
  ana_event_table %>% 
    select(
      node, abs_event_time, event_type, event_txt, dispersal_to) %>%
    mutate(from = str_split_fixed(event_txt, n = 2, pattern = "->")[,1 ])
}

# Count anagenic events across reps
count_angen_across_reps <- function(ana_events_tables) {
  map_df(ana_events_tables, count_angen_single, .id = "rep") %>%
  as_tibble()
}

# Count cladogenic and anagenic events across reps
count_bsm_across_reps <- function(bsm_results) {
  cladgen <- count_cladgen_across_reps(bsm_results$RES_clado_events_tables)
  angen <- count_angen_across_reps(bsm_results$RES_ana_events_tables)
  bind_rows(cladgen, angen) %>%
    rename(to = dispersal_to) %>%
    select(rep, node, abs_event_time, event_type, event_txt, from, to)
}

#' Count mean number of events per replicate
#' @param ... Bare names of events (variables) to count
count_mean_bsm <- function(events, ...) {
  events %>%
    # Add count of each event type per replicate
    group_by(..., rep) %>%
    summarize(count = dplyr::n(), .groups = "drop") %>%
    # This count won't include zeros since they don't appear in the raw data.
    # Add them
    complete(..., rep) %>%
    mutate(count = replace_na(count, 0)) %>%
    # Calculate mean number of events across replicates
    group_by(...) %>%
    summarize(
      mean = mean(count), sd = sd(count),
      .groups = "drop") %>%
    arrange(desc(mean))
}

#' Bin biogeographic events
#'
#' @param events Table of cladogenic or anagenic events from BSM
#' @param event_type_select Character vector; Type of event to bin;
#' should be a column name or names in `events`
#' @param bin_width Temporal width of bin
#'
bin_events <- function(
  events,
  event_type_select,
  bin_width) {

  assertthat::assert_that(
    all(event_type_select %in% events$event_type),
    msg = "One or more selected event types not in events table"
  )

  vars_1 <- unique(c("rep", vars))

# Count the number of dispersals per time bin
  events %>%
  filter(event_type %in% event_type_select) %>%
  mutate(
    bin = cut_width(
      abs_event_time,
      width = bin_width,
      center = bin_width / 2)
  ) %>%
  group_by(bin, rep) %>%
  summarize(abs_count = dplyr::n(),
            .groups = "drop") %>%
  # Account for missing observations
  complete(bin, rep) %>%
  mutate(abs_count = replace_na(abs_count, 0)) %>%
  group_by(bin) %>%
  summarize(lower = quantile(abs_count, probs = 0.025),
            mean = mean(abs_count),
            upper = quantile(abs_count, probs = 0.975),
            sd = sd(abs_count),
            n = dplyr::n(),
            .groups = "drop")
}

#' Normalize binned beogeographic events
#'
#' @param events_binned Binned events; output of bin_events()
#' @param blrn_binned Binned branchlengths; output of bin_brln()
normalize_events <- function(events_binned, blrn_binned) {
  left_join(events_binned, blrn_binned, by = "bin") %>%
    mutate(lower = lower / tot_brl) %>%
    mutate(mean = mean / tot_brl) %>%
    mutate(upper = upper / tot_brl) %>%
    mutate(sd = sd / tot_brl) %>%
  # If no events observed in that bin, convert to NA
  mutate(mean = case_when(
    mean == 0 & sd == 0 ~ NaN,
    TRUE ~ mean
  ))
}

#' Normalize dispersal events