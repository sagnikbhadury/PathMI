#!/bin/bash
# 27 Nov 2026 final-edits rerun: back up affected outputs, then regenerate everything downstream of the
# max-16 check and the revised tolerance. Run from the project root.
set -u
R="/c/Users/bhadury/AppData/Local/Programs/R/R-4.4.3/bin/Rscript.exe"
BK="TRIDENT/archive_pre_final_edits_20261127"
OUTD="TRIDENT/patient_level/final_edits_out"
mkdir -p "$BK/record" "$BK/learning" "$BK/robustness" "$BK/manuscript_tables" "$BK/manuscript_figures" "$OUTD"
if [ ! -f "$BK/.done" ]; then
  cp -r TRIDENT/record/. "$BK/record/"
  cp TRIDENT/learning/learning_summary.csv TRIDENT/learning/n_evidence.csv "$BK/learning/"
  cp TRIDENT/robustness/resolution_units.csv "$BK/robustness/" 2>/dev/null
  cp manuscript/table_*.tex "$BK/manuscript_tables/"
  cp manuscript/figures/*.pdf "$BK/manuscript_figures/"
  touch "$BK/.done"
fi
step() { local n=$1; shift; echo "== step $n: $*"; "$@" > "$OUTD/$n.out" 2>&1; local rc=$?; echo "   exit $rc"; tail -3 "$OUTD/$n.out"; }
step 01_tolerance        "$R" TRIDENT/tolerance.R
step 02_learning_summary "$R" TRIDENT/learning_summary.R
step 03_record           "$R" TRIDENT/recoverability_record.R
step 04_status_recheck   "$R" TRIDENT/status_recheck.R
step 05_resolution       "$R" TRIDENT/resolution_units.R
step 06_resolution_table "$R" TRIDENT/resolution_table.R
step 07_fig_four         "$R" TRIDENT/figures_four_backbones.R
step 08_fig_results      "$R" TRIDENT/figures_results_four.R
step 09_fig_robust       "$R" TRIDENT/figures_robustness.R
step 10_fig_learning     "$R" TRIDENT/figures_learning_context.R
step 11_fig_tiles        "$R" TRIDENT/figures_tile_maps.R
step 12_fig_pipeline     "$R" TRIDENT/figures_pipeline_walkthrough.R
echo "== chain done"
