################################################################
# mirror_new_backbone_results.R
#
# Copies the finished H-optimus-1 and CONCH v1.5 result folders from
# the Turbo results area (M:/bhadury/mi_results/) into the project,
# following the convention used for UNI v2 and Virchow2: one copy
# under the cohort folder (CRC/ or PRAD/) and one under TRIDENT/,
# byte-identical, with TRIDENT/ being the copy the merged tables and
# pathmi_figures.R read.
#
#   M:/bhadury/mi_results/crc_kpca_subho_hoptimus1        -> CRC/,  TRIDENT/
#   M:/bhadury/mi_results/crc_kpca_subho_conch_v15        -> CRC/,  TRIDENT/
#   M:/bhadury/mi_results/TCGA_PRAD_kpca_subho_hoptimus1  -> PRAD/, TRIDENT/
#   M:/bhadury/mi_results/TCGA_PRAD_kpca_subho_conch_v15  -> PRAD/, TRIDENT/
#
# A folder is only mirrored once its pipeline has written the final
# rds (kpca_subho_results_*.rds) AND at least one figures_*.png, i.e.
# the run reached the last block. Unfinished runs are reported and
# skipped, so this can be rerun safely at any time.
#
# Usage (from anywhere): Rscript TRIDENT/mirror_new_backbone_results.R
################################################################

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root_dir  <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
M_RES     <- "M:/bhadury/mi_results"

jobs <- list(
  list(src = "crc_kpca_subho_hoptimus1",       cohort = "CRC"),
  list(src = "crc_kpca_subho_conch_v15",       cohort = "CRC"),
  list(src = "TCGA_PRAD_kpca_subho_hoptimus1", cohort = "PRAD"),
  list(src = "TCGA_PRAD_kpca_subho_conch_v15", cohort = "PRAD")
)

for (j in jobs) {
  src <- file.path(M_RES, j$src)
  if (!dir.exists(src)) { cat(sprintf("[skip] %s: not on M:/ yet\n", j$src)); next }
  files <- list.files(src, full.names = TRUE)
  done  <- any(grepl("^kpca_subho_results_.*[.]rds$", basename(files))) &&
           any(grepl("^figures_.*[.]png$",            basename(files)))
  if (!done) {
    cat(sprintf("[skip] %s: run not finished (%d files so far)\n", j$src, length(files)))
    next
  }
  for (dst_parent in c(j$cohort, "TRIDENT")) {
    dst <- file.path(root_dir, dst_parent, j$src)
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(files, dst, overwrite = TRUE)
    # Verify byte-identity of every copied file against M:/.
    same <- vapply(files, function(f)
      isTRUE(all.equal(file.size(f), file.size(file.path(dst, basename(f))))) &&
      identical(unname(tools::md5sum(f)), unname(tools::md5sum(file.path(dst, basename(f))))),
      logical(1))
    cat(sprintf("[ok]   %s -> %s/  (%d / %d files, %d identical)\n",
                j$src, dst_parent, sum(ok), length(files), sum(same)))
  }
}
