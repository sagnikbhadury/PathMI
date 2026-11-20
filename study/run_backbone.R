################################################################
# run_backbone.R
#
# Runs one cohort x backbone pass of the kernel-PCA SUBHO pipeline
# as a standalone Rscript process and copies the finished result
# folder from the Turbo results area (M:/bhadury/mi_results/) into
# TRIDENT/, next to the four existing result folders.
#
# Usage (from anywhere; paths are resolved from this file's location):
#   Rscript TRIDENT/run_backbone.R <cohort> <backbone>
#     cohort   : crc | prad
#     backbone : hoptimus1 | conch_v15
#
# Drives the new-backbone sibling scripts
# kPCA_SUBHO_{CRC,PRAD}_trident_hoptimus_conch.R, which read
# PATHMI_EMBED_TYPE, so nothing is edited per run. The original
# UNI v2 / Virchow2 scripts are not touched by this runner. Four
# runs can be launched in parallel as separate processes; each
# writes to its own EMBED_TYPE-suffixed results directory.
#
# Runtime reference (July 2026 runs, same settings): PRAD uni_v2
# took about 5.3 h for 16 decompositions x 2 targets x 201 SUBHO
# fits; CRC (one target) about half that.
#
# Use R 4.4.3: R 4.5.2 segfaults on this machine.
################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L)
  stop("usage: Rscript run_backbone.R <crc|prad> <hoptimus1|conch_v15>")

cohort   <- tolower(args[1])
backbone <- args[2]
stopifnot(cohort %in% c("crc", "prad"),
          backbone %in% c("hoptimus1", "conch_v15"))

# Locate this file, then the project root (its parent).
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
trident_dir <- normalizePath(dirname(this_file), winslash = "/")
root_dir    <- normalizePath(file.path(trident_dir, ".."), winslash = "/")

# The pipeline scripts source("PathMI_SUBHO_logistic.R") relative to
# the working directory, which must therefore be the project root.
setwd(root_dir)

script <- switch(cohort,
                 crc  = file.path(trident_dir, "kPCA_SUBHO_CRC_trident_hoptimus_conch.R"),
                 prad = file.path(trident_dir, "kPCA_SUBHO_PRAD_trident_hoptimus_conch.R"))

Sys.setenv(PATHMI_EMBED_TYPE = backbone)

t0 <- Sys.time()
cat(sprintf("[run_backbone] %s | cohort=%s backbone=%s\n  script=%s\n  cwd=%s\n",
            format(t0), cohort, backbone, script, getwd()))

# The pipeline script calls rm(list = ls()) at its top, which clears
# every variable defined above from the global environment. Stash
# what the copy step needs in options(), which rm() does not touch.
options(pathmi_runner = list(trident_dir = trident_dir, t0 = t0))

# RESULTS_DIR and EMBED_TYPE are left in the global environment by
# the pipeline and are read back below.
source(script, echo = FALSE)

st <- getOption("pathmi_runner")
t1 <- Sys.time()
cat(sprintf("\n[run_backbone] pipeline finished in %.1f h\n",
            as.numeric(difftime(t1, st$t0, units = "hours"))))

# Mirror the M:/ result folder into TRIDENT/<same folder name>.
src <- sub("/+$", "", RESULTS_DIR)
dst <- file.path(st$trident_dir, basename(src))
dir.create(dst, recursive = TRUE, showWarnings = FALSE)
files <- list.files(src, full.names = TRUE, recursive = FALSE)
ok <- file.copy(files, dst, overwrite = TRUE)
cat(sprintf("[run_backbone] copied %d / %d files -> %s\n", sum(ok), length(files), dst))
cat(sprintf("[run_backbone] done %s\n", format(Sys.time())))
