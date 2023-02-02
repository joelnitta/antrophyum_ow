library(targets)
library(tarchetypes)
library(tidyverse)
library(optimx)
library(GenSA)
library(FD)
library(snow)
library(parallel)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

source("R/functions.R")

tar_plan(
  # Specify data file paths ----
  tar_file(antro_trfn, "_targets/user/data/Antrophyum.newick"),
  tar_file(antro_geogfn, "_targets/user/data/Antrophyum_geodata.data"),
  tar_file(antro_multfn,
    "_targets/user/data/manual_dispersal_multipliers.txt"),
  # Set up BioGeoBEARS models ----
  antro_max_range = 6,
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
  decj_settings = setup_bgb_dec(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range,
    jump = TRUE,
    res_dec = dec_model
  ),
  # - DEC-j constrained
  decj_const_settings = setup_bgb_dec(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range,
    jump = TRUE,
    res_dec = dec_model
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
  divaj_settings = setup_bgb_diva(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range,
    jump = TRUE,
    res_diva = diva_model
  ),
  # - DIVA-j constrained
  divaj_const_settings = setup_bgb_diva(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range,
    jump = TRUE,
    res_diva = diva_model
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
  bayareaj_settings = setup_bgb_bayarea(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    max_range_size = antro_max_range,
    jump = TRUE,
    res_bayarea = bayarea_model
  ),
  # - BAYAREA-j constrained
  bayareaj_const_settings = setup_bgb_bayarea(
    trfn = antro_trfn,
    geogfn = antro_geogfn,
    dispersal_multipliers_fn = antro_multfn,
    max_range_size = antro_max_range,
    jump = TRUE,
    res_bayarea = bayarea_model
  ),
  # Run BioGeoBEARS ----
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
  # Compare models with AIC ----
  bgb_stats = get_bgb_stats(
    dec_model = dec_model,
    dec_const_model = dec_const_model,
    decj_model = decj_model,
    decj_const_model = decj_const_model,
    diva_model = diva_model,
    diva_const_model = diva_const_model,
    divaj_model = divaj_model,
    divaj_const_model = divaj_const_model,
    bayarea_model = bayarea_model,
    bayarea_const_model = bayarea_const_model,
    bayareaj_model  = bayareaj_model,
    bayareaj_const_model = bayareaj_const_model),
  # Format table for MS
  bgb_stats_pretty = format_model_table(bgb_stats),
  # Extract best-scoring model for BSM
  best_model = get_best_model(bgb_stats),
  # Run BSM ----
  stochastic_mapping_inputs_list =
    BioGeoBEARS::get_inputs_for_stochastic_mapping(res = best_model),
  bsm_results = run_bsm_quietly(
    res = best_model,
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