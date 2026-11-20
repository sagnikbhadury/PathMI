#!/usr/bin/env bash
# run_patient_chain.sh -- the downstream chain of the patient-level rerun (19 Nov 2026).
# Run from anywhere after pl_grid.R (patient mode) has finished. Each step is logged to
# TRIDENT/patient_level/chain.log with its exit status; the chain continues past a failure
# so independent steps still run, and a failed step can be rerun alone.
ROOT="/c/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe"
R="/c/Users/bhadury/AppData/Local/Programs/R/R-4.4.3/bin/Rscript.exe"
cd "$ROOT" || exit 1
LOG="TRIDENT/patient_level/chain.log"; OUTD="TRIDENT/patient_level/chain_out"; mkdir -p "$OUTD"
START_AT="${1:-1}"
step() { local n="$1"; shift; local name="$1"; shift
  if [ "$n" -lt "$START_AT" ]; then return; fi
  echo "$(date '+%H:%M:%S') START $n $name" >> "$LOG"
  "$@" > "$OUTD/$(printf '%02d' "$n")_$name.out" 2>&1; local rc=$?
  echo "$(date '+%H:%M:%S') END   $n $name rc=$rc" >> "$LOG"; }
step 1  incremental      "$R" TRIDENT/patient_level/pl_incremental.R patient
step 2  promote          "$R" TRIDENT/patient_level/pl_promote.R
step 3  promote_post     python TRIDENT/patient_level/pl_promote_post.py
step 4  clear_stale      bash -c 'rm -f TRIDENT/learning/*.csv TRIDENT/learning/*.rds TRIDENT/learning/*.log TRIDENT/record/*.csv TRIDENT/record/*.json; rm -rf TRIDENT/record/json; rm -f TRIDENT/usable/usable_infold_slide*.csv'
step 5  usable           env N_WORKERS=${NW:-20} "$R" TRIDENT/usable_infold.R
step 6  subsampling      env N_WORKERS=${NW:-20} "$R" TRIDENT/robustness_subsampling.R
step 7  prior_bootstrap  env N_WORKERS=${NW:-20} "$R" TRIDENT/robustness_prior_bootstrap.R
step 8  dimension        env N_WORKERS=${NW:-20} "$R" TRIDENT/robustness_dimension.R
step 9  burden           "$R" TRIDENT/sampling_burden_all_models.R
step 10 burden_boot      "$R" TRIDENT/burden_ordering_bootstrap.R
step 11 burden_table     "$R" TRIDENT/burden_table.R
step 12 tiles_crc        "$R" TRIDENT/tile_maps_hoptimus1.R crc
step 13 tiles_prad       "$R" TRIDENT/tile_maps_hoptimus1.R prad
step 14 fig_tiles        "$R" TRIDENT/figures_tile_maps.R
step 15 record_pass1     "$R" TRIDENT/recoverability_record.R
step 16 learning         env N_WORKERS=${NW:-20} "$R" TRIDENT/learning_curves.R
step 17 learning_summary "$R" TRIDENT/learning_summary.R
step 18 record_pass2     "$R" TRIDENT/recoverability_record.R
step 19 resolution       "$R" TRIDENT/resolution_units.R
step 20 resolution_table "$R" TRIDENT/resolution_table.R
step 21 usable_stats     "$R" TRIDENT/usable_infold_stats.R
step 22 usable_table     "$R" TRIDENT/usable_table.R
step 23 robust_tables    "$R" TRIDENT/robustness_tables.R
step 24 subsamp_table    "$R" TRIDENT/subsampling_table.R
step 25 numbers          "$R" TRIDENT/results_four_backbones_numbers.R
step 26 config_tests     "$R" TRIDENT/config_comparison_tests.R
step 27 fig_four         "$R" TRIDENT/figures_four_backbones.R
step 28 fig_results      "$R" TRIDENT/figures_results_four.R
step 29 fig_robust       "$R" TRIDENT/figures_robustness.R
step 30 fig_learning     "$R" TRIDENT/figures_learning_context.R
echo "$(date '+%H:%M:%S') CHAIN DONE" >> "$LOG"
