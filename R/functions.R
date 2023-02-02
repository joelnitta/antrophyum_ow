
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

  # Basic settings
  # - default model is DEC
  BioGeoBEARS_run_object <- BioGeoBEARS::define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn <- trfn
  BioGeoBEARS_run_object$geogfn <- geogfn
  BioGeoBEARS_run_object$max_range_size <- max_range_size
  BioGeoBEARS_run_object$min_branchlength <- min_branchlength # Min to treat tip as a direct ancestor (no speciation event) # nolint
  BioGeoBEARS_run_object$include_null_range <- TRUE # set to FALSE for e.g. DEC* model, DEC*+J, etc. # nolint

  # (Optional) specify time slices file
  if (!is.null(timesfn)) BioGeoBEARS_run_object$timesfn <- timesfn

  # (Optional) specify dispersal multipliers file
  if (!is.null(dispersal_multipliers_fn)) {
    BioGeoBEARS_run_object$dispersal_multipliers_fn <- dispersal_multipliers_fn
  }

  # Speed options and multicore processing if desired
  BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check) #nolint
  BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params) #nolint
  BioGeoBEARS_run_object$use_optimx <- use_optimx    # if FALSE, use optim() instead of optimx() #nolint
  BioGeoBEARS_run_object$num_cores_to_use <- num_cores_to_use
  BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale #nolint
  BioGeoBEARS_run_object$cluster_already_open <- cluster_already_open

  # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work! #nolint
  # (It also runs some checks on these inputs for certain errors.)
  BioGeoBEARS_run_object <- BioGeoBEARS::readfiles_BioGeoBEARS_run(
    BioGeoBEARS_run_object)

  # Divide the tree up by timeperiods/strata for stratified analysis
  if (!is.null(timesfn)) {
    BioGeoBEARS_run_object <- BioGeoBEARS::section_the_tree(
      inputs = BioGeoBEARS_run_object,
      make_master_table = TRUE, plot_pieces = FALSE)
  }

  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table <- TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
  BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run #nolint

  # Set up DEC model
  # (nothing to do; defaults)

  # DEC+J
  if (jump == TRUE) {
    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart <- res_dec$outputs@params_table["d","est"]
    estart <- res_dec$outputs@params_table["e","est"]
    jstart <- 0.0001

    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] <- dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] <- estart
    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] <- jstart
  }

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  BioGeoBEARS_run_object

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

  # Basic settings
  # - default model is DEC
  BioGeoBEARS_run_object <- BioGeoBEARS::define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn <- trfn
  BioGeoBEARS_run_object$geogfn <- geogfn
  BioGeoBEARS_run_object$max_range_size <- max_range_size
  BioGeoBEARS_run_object$min_branchlength <- min_branchlength
  BioGeoBEARS_run_object$include_null_range <- TRUE

  # (Optional) specify time slices file
  if (!is.null(timesfn)) BioGeoBEARS_run_object$timesfn <- timesfn

  # (Optional) specify dispersal multipliers file
  if (!is.null(dispersal_multipliers_fn)) {
    BioGeoBEARS_run_object$dispersal_multipliers_fn <- dispersal_multipliers_fn
  }

  # Speed options and multicore processing
  BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check) #nolint
  BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params) #nolint
  BioGeoBEARS_run_object$use_optimx <- use_optimx    # if FALSE, use optim() instead of optimx() #nolint
  BioGeoBEARS_run_object$num_cores_to_use <- num_cores_to_use
  BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse = TRUE causes pathology & isn't much faster at this scale #nolint
  BioGeoBEARS_run_object$cluster_already_open <- cluster_already_open

  # Loads the dispersal multiplier matrix etc.
  BioGeoBEARS_run_object <- BioGeoBEARS::readfiles_BioGeoBEARS_run(
    BioGeoBEARS_run_object)

  # Divide the tree up by timeperiods/strata for stratified analysis
  if (!is.null(timesfn)) {
    BioGeoBEARS_run_object <- BioGeoBEARS::section_the_tree(
      inputs = BioGeoBEARS_run_object,
      make_master_table = TRUE, plot_pieces = FALSE)
  }

  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table <- TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
  BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run #nolint

  if (jump == FALSE) {
    # Set up DIVALIKE model
    #
    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <-"fixed" #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0 #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0 #nolint

    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "2-j" #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/2" #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "ysv*1/2" #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "ysv*1/2" #nolint

    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] <- "fixed" #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] <- 0.5 #nolint
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] <- 0.5 #nolint
  } else if (jump == TRUE) {
    # Set up DIVALIKE+J model
    #
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart <- res_diva$outputs@params_table["d", "est"]
    estart <- res_diva$outputs@params_table["e", "est"]
    jstart <- 0.0001

    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "init"] <- dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d", "est"] <- dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "init"] <- estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e", "est"] <- estart

    # Remove subset-sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y", "type"] <- "ysv*1/2"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"

    # Allow classic, widespread vicariance; all events equiprobable
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "init"] <- 0.5
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "est"] <- 0.5

    # Add jump dispersal/founder-event speciation
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] <- jstart

    # Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "min"] <- 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "max"] <- 1.99999
  }

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  BioGeoBEARS_run_object

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

  # Basic settings
  BioGeoBEARS_run_object <- BioGeoBEARS::define_BioGeoBEARS_run()
  BioGeoBEARS_run_object$trfn <- trfn
  BioGeoBEARS_run_object$geogfn <- geogfn
  BioGeoBEARS_run_object$max_range_size <- max_range_size
  BioGeoBEARS_run_object$min_branchlength <- min_branchlength
  BioGeoBEARS_run_object$include_null_range <- TRUE

  # (Optional) specify time slices file
  if (!is.null(timesfn)) BioGeoBEARS_run_object$timesfn <- timesfn

  # (Optional) specify dispersal multipliers file
  if (!is.null(dispersal_multipliers_fn)) {
    BioGeoBEARS_run_object$dispersal_multipliers_fn <- dispersal_multipliers_fn
  }

  # Speed options and multicore processing
  BioGeoBEARS_run_object$on_NaN_error <- -1e50    # returns very low lnL if parameters produce NaN error (underflow check) #nolint
  BioGeoBEARS_run_object$speedup <- TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params) #nolint
  BioGeoBEARS_run_object$use_optimx <- use_optimx    # if FALSE, use optim() instead of optimx() #nolint
  BioGeoBEARS_run_object$num_cores_to_use <- num_cores_to_use
  BioGeoBEARS_run_object$force_sparse <- FALSE    # force_sparse = TRUE causes pathology & isn't much faster at this scale #nolint
  BioGeoBEARS_run_object$cluster_already_open <- cluster_already_open

  # Loads the dispersal multiplier matrix etc.
  BioGeoBEARS_run_object <- BioGeoBEARS::readfiles_BioGeoBEARS_run(
    BioGeoBEARS_run_object)

  # Divide the tree up by timeperiods/strata for stratified analysis
  if (!is.null(timesfn)) {
    BioGeoBEARS_run_object <- BioGeoBEARS::section_the_tree(
      inputs = BioGeoBEARS_run_object,
      make_master_table = TRUE, plot_pieces = FALSE)
  }

  # Good default settings to get ancestral states
  BioGeoBEARS_run_object$return_condlikes_table <- TRUE
  BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
  BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run #nolint

  if (jump == FALSE) {

    # Set up BAYAREALIKE model
    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "init"] <- 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "est"] <- 0.0

    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y", "type"] <- "1-j"

    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "init"] <- 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y", "est"] <- 0.9999

  } else if (jump == TRUE) {
    # Set up BAYAREALIKE+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart <- resBAYAREALIKE$outputs@params_table["d","est"]
    estart <- resBAYAREALIKE$outputs@params_table["e","est"]
    jstart <- 0.0001

    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] <- dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] <- dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] <- estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] <- estart

    # No subset sympatry
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0

    # No vicariance
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] <- 0.0
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] <- 0.0

    # *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] <- jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] <- jstart

    # Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 0.99999

    # Adjust linkage between parameters
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/1"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "1-j"

    # Only sympatric/range-copying (y) events allowed, and with 
    # exact copying (both descendants always the same size as the ancestor)
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] <- 0.9999

    # NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers,
    # usually Windows machines. I can't replicate this on my Mac machines, but
    # it is almost certainly just some precision under-run issue, when
    # optim/optimx tries some parameter value just below zero.  The "min" and
    # "max" options on each parameter are supposed to prevent this, but
    # apparently optim/optimx sometimes go slightly beyond these limits.
    # Anyway, if you get a crash, try raising "min" and lowering "max" slightly
    # for each parameter:
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] <- 0.0000001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] <- 4.9999999

    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] <- 0.0000001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] <- 4.9999999

    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] <- 0.00001
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] <- 0.99999
  }

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  BioGeoBEARS_run_object

}

#' Run BioGeoBEARS::bears_optim_run() quietly
#' 
#' Instead of printing all the results to the screen, which
#' it does by default
#'
#' @param bgb_run_object List; a BioGeoBEARS object
#' @param ... Other arguments not used by this function, but meant for
#' tracking with drake.
#'
#' @return List
#' 
run_bgb_quietly = function (bgb_run_object, ...) {

  # Check that setup is valid
  BioGeoBEARS::check_BioGeoBEARS_run(bgb_run_object)

  # Suppress output printed to screen
  invisible(capture.output(results <- BioGeoBEARS::bears_optim_run(bgb_run_object)))

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
run_bsm_quietly = function (...) {

  # Suppress output printed to screen
  invisible(capture.output(results <- BioGeoBEARS::runBSM(...)))

  results

}

#' Compare BioGeoBEARS models using AIC
#'
#' @param resDEC Results of running BioGeoBEARS using the DEC model 
#' @param resDECj Results of running BioGeoBEARS using the DECj model 
#' @param resDIVALIKE Results of running BioGeoBEARS using the DIVALIKE model 
#' @param resDIVALIKEj Results of running BioGeoBEARS using the DIVALIKE+J model 
#' @param resBAYAREALIKE Results of running BioGeoBEARS using the BAYAREALIKE model 
#' @param resBAYAREALIKEj Results of running BioGeoBEARS using the BAYAREALIKE+J model 
#' @param tr Phylogenetic tree used in BioGeoBEARS analysis
#' @param return String, either "tests" or "aic"
#'
#' @return One of two dataframes, depending on `return`
#'   - if return = "tests": results comparing each pair of models (e.g, DEC vs DEC+J)
#'   - if return = "aic": AIC and corrected AIC comparing all models
#' 
get_bgb_stats <- function(resDEC,
                          resDECj,
                          resDIVALIKE,
                          resDIVALIKEj,
                          resBAYAREALIKE,
                          resBAYAREALIKEj,
                          tr,
                          return = c("tests", "aic")) {

  assert_that(return %in% c("tests", "aic"))

  # Set up empty tables to hold the statistical results
  restable = NULL
  teststable = NULL

  ### Statistics -- DEC vs. DEC+J ###

  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object(resDEC)
  LnL_1 = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object(resDECj)

  numparams1 = 3
  numparams2 = 2
  stats = BioGeoBEARS::AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

  # DEC, null model for Likelihood Ratio Test (LRT)
  res2 = BioGeoBEARS::extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DEC+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = BioGeoBEARS::extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

  # The null hypothesis for a Likelihood Ratio Test (LRT) is that two models
  # confer the same likelihood on the data. See: Brian O'Meara's webpage:
  # http://www.brianomeara.info/tutorials/aic
  # ...for an intro to LRT, AIC, and AICc
  tmp_tests = BioGeoBEARS::conditional_format_table(stats)

  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)

  ### Statistics -- DIVALIKE vs. DIVALIKE+J ###

  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
  LnL_1 = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

  numparams1 = 3
  numparams2 = 2
  stats = BioGeoBEARS::AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

  # DIVALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = BioGeoBEARS::extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = BioGeoBEARS::extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

  tmp_tests = conditional_format_table(stats)

  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)

  ### Statistics -- BAYAREALIKE vs. BAYAREALIKE+J ###

  # We have to extract the log-likelihood differently, depending on the 
  # version of optim/optimx
  LnL_2 = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
  LnL_1 = BioGeoBEARS::get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

  numparams1 = 3
  numparams2 = 2
  stats = BioGeoBEARS::AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

  # BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
  res2 = BioGeoBEARS::extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
  res1 = BioGeoBEARS::extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)

  tmp_tests = BioGeoBEARS::conditional_format_table(stats)

  restable = rbind(restable, res2, res1)
  teststable = rbind(teststable, tmp_tests)

  ### ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J ###

  teststable$alt = c("DEC+J", "DIVALIKE+J", "BAYAREALIKE+J")
  teststable$null = c("DEC", "DIVALIKE", "BAYAREALIKE")
  row.names(restable) = c("DEC", "DEC+J", "DIVALIKE", "DIVALIKE+J", "BAYAREALIKE", "BAYAREALIKE+J")
  restable = BioGeoBEARS::put_jcol_after_ecol(restable)

  ### Model weights of all six models ###

  # With AICs:
  AICtable = BioGeoBEARS::calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
  restable2 = cbind(restable, AICtable)
  restable_AIC_rellike = BioGeoBEARS::AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AIC")
  restable_AIC_rellike = BioGeoBEARS::put_jcol_after_ecol(restable_AIC_rellike)

  # With AICcs (factors in sample size)
  samplesize = length(tr$tip.label)
  AICtable = BioGeoBEARS::calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
  restable2 = cbind(restable, AICtable)
  restable_AICc_rellike = BioGeoBEARS::AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
  restable_AICc_rellike = BioGeoBEARS::put_jcol_after_ecol(restable_AICc_rellike)

  ### Format as tibbles for output ###

  teststable <- as_tibble(teststable) %>%
    mutate_if(is.list, purrr::flatten_chr) %>%
    mutate_at(vars(-alt, -null, -test, -tail), readr::parse_double)

  restable <-
    restable %>% 
    rownames_to_column("model") %>%
    as_tibble

  AIC_table <-
    restable_AIC_rellike %>% 
    rownames_to_column("model") %>%
    as_tibble

  AICc_table <-
    restable_AICc_rellike %>% 
    rownames_to_column("model") %>%
    as_tibble

  AIC_table <- inner_join(AIC_table, AICc_table, by = c("model", "LnL", "numparams", "d", "e", "j"))

  if (return == "tests") return(teststable)

  if (return == "aic") return(AIC_table)

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
