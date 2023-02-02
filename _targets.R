library(targets)
library(tarchetypes)
library(tidyverse)

source("R/functions.R")

tar_plan(
  tar_file(antro_trfn, "Antrophyum.newick"),
  tar_file(antro_geogfn, "Antrophyum_geodata.data"),
  tar_file(antro_multfn, "manual_dispersal_multipliers.txt"),
  antro_max_range = 6,
  # Set up BioGeoBEARS models
  # - DEC
  dec_settings = setup_bgb_dec(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range
  ),
  # - DEC constrained
  dec_const_settings = setup_bgb_dec(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range
  ),
  # - DEC-j
  decj_settings = setup_bgb_decj(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range,
    resDEC = dec_model
  ),
  # - DEC-j constrained
  decj_const_settings = setup_bgb_decj(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range,
    resDEC = dec_model
  ),
  # - DIVA
  diva_settings = setup_bgb_diva(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range
  ),
  # - DIVA constrained
  diva_const_settings = setup_bgb_diva(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range
  ),
  # - DIVA-j
  divaj_settings = setup_bgb_divaj(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range,
    resDIVALIKE = diva_model
  ),
  # - DIVA-j constrained
  divaj_const_settings = setup_bgb_divaj(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range,
    resDIVA = diva_model
  ),
  # - BAYAREA
  bayarea_settings = setup_bgb_bayarea(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range
  ),
  # - BAYAREA constrained
  bayarea_const_settings = setup_bgb_bayarea(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range
  ),
  # - BAYAREA-j
  bayareaj_settings = setup_bgb_bayareaj(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range,
    resBAYAREA = bayarea_model
  ),
  # - BAYAREA-j constrained
  bayareaj_const_settings = setup_bgb_bayareaj(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range,
    resBAYAREA = bayarea_model
  ),
  # Run BioGeoBEARS
  # - DEC
  dec_model = run_bgb_quietly(dec_settings),
  dec_const_model = run_bgb_quietly(dec_const_settings),
  decj_model = run_bgb_quietly(decj_settings),
  decj_const_model = run_bgb_quietly(decj_const_settings),
  # - DIVA
  diva_model = run_bgb_quietly(diva_settings),
  diva_const_model = run_bgb_quietly(diva_const_settings),
  divaj_model = run_bgb_quietly(divaj_settings),
  divaj_const_model = run_bgb_quietly(divaj_const_settings),
  # - BAYAREA
  bayarea_model = run_bgb_quietly(bayarea_settings),
  bayarea_const_model = run_bgb_quietly(bayarea_const_settings),
  bayareaj_model = run_bgb_quietly(bayareaj_settings),
  bayareaj_const_model = run_bgb_quietly(bayareaj_const_settings),
  # Do stochastic character mapping to estimate range shifts over time
  stochastic_mapping_inputs_list = 
    BioGeoBEARS::get_inputs_for_stochastic_mapping(res = dec_model),
  bsm_results = run_bsm_quietly(
    res = dec_model,
    stochastic_mapping_inputs_list = stochastic_mapping_inputs_list,
    maxnum_maps_to_try = 2000,
    nummaps_goal = 100,
    maxtries_per_branch = 40000,
    save_after_every_try = TRUE,
    savedir = "_targets/user/biogeobears/",
    seedval = 12345,
    wait_before_save = 0.01),
  
  # Count normalized dispersal events through time
  antro_phy = ape::read.tree(antro_trfn) #,
  # dispersal_events = count_normalized_dispersal_events(
  #   bsm_results = bsm_results,
  #   max_age = 20,
  #   bin_width = 2,
  #   tr = test_phy,
  #   realm_names = test_realms %>% pull(realm) %>% unique %>% sort
  # )

)