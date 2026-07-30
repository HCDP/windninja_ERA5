#!/usr/bin/env Rscript
# make_era5_hourly.R
#
# Split an ERA5 single-level GRIB into SINGLE-TIMESTEP WRF-mimic NetCDF files,
# one per hour, named by UTC valid time: wrfout_era5_YYYYMMDD_HHMM.nc
# (the file naming wn_era5_hi_1hr.py and the test run scripts expect).
#
# This is the driver used to produce the per-hour forecast files that WindNinja
# needs (it runs EVERY band in a supplied forecast file, so 1-band files are how
# an hourly run runs one hour). It wraps the proven converter era5_to_wrfout.R
# (same directory) - one conversion call per timestep in the GRIB.
#
# Usage:
#   Rscript make_era5_hourly.R --input <era5.grib> --outdir <dir> [--res m]
#
#   --res  target Mercator grid spacing in meters. Default 27000 (~ERA5 native
#          0.25 deg). For SMALL downloaded domains the converter requires at
#          least a 4x4 target grid - reduce res accordingly (e.g. a 1.25 deg
#          box needs res <= ~13000).
#
# Runs anywhere with R + terra + RNetCDF (developed/run on the local Windows
# box; the outputs are then copied to the server's input/er5_wrf/hourly/ tree).

suppressMessages(library(terra))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}

in_path <- get_arg("--input")
out_dir <- get_arg("--outdir")
res_m   <- get_arg("--res", "27000")
if (is.null(in_path) || is.null(out_dir)) {
  stop("usage: Rscript make_era5_hourly.R --input <era5.grib> --outdir <dir> [--res m]")
}

converter <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])),
                       "era5_to_wrfout.R")
if (!file.exists(converter)) stop("converter not found next to this script: ", converter)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setGDALconfig("GRIB_NORMALIZE_UNITS", "NO")
times <- sort(unique(time(rast(in_path))))
cat(sprintf("%s: %d hourly timesteps, %s to %s UTC\n",
            basename(in_path), length(times), format(times[1]), format(times[length(times)])))

rscript <- file.path(R.home("bin"), "Rscript")
n_ok <- 0; n_skip <- 0; n_fail <- 0
for (tt in times) {
  tt <- as.POSIXct(tt, origin = "1970-01-01", tz = "UTC")
  out_name <- format(tt, "wrfout_era5_%Y%m%d_%H%M.nc")
  out_path <- file.path(out_dir, out_name)
  if (file.exists(out_path)) { n_skip <- n_skip + 1; next }
  ts <- format(tt, "%Y-%m-%d %H:%M")
  status <- system2(rscript, c(converter,
                               "--input", shQuote(in_path),
                               "--output", shQuote(out_path),
                               "--start", shQuote(ts), "--stop", shQuote(ts),
                               "--res", res_m),
                    stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
    n_fail <- n_fail + 1
    cat("FAILED:", out_name, "\n")
    cat(tail(status, 4), sep = "\n")
  } else {
    n_ok <- n_ok + 1
  }
}
cat(sprintf("done: %d written, %d skipped (existing), %d failed -> %s\n",
            n_ok, n_skip, n_fail, out_dir))
if (n_fail > 0) quit(status = 1)
