#!/usr/bin/env Rscript
# plot_wind_compare.R
#
# Paired comparison of two WindNinja runs over the same DEM: for every hour
# present in BOTH output directories, a 3-panel PNG (run1 | run2 | run2-run1
# difference) with the speed panels on ONE shared scale across both runs and
# all hours, and the difference panel on a symmetric diverging scale (global
# max |diff|). Also assembles the panels into an animated GIF.
#
# Expects WindNinja CLI ascii outputs <dem>_<MM-DD-YYYY>_<HHMM>_<res>m_vel.asc
# in each directory (timestamps local to the run's time_zone).
#
# Usage:
#   Rscript plot_wind_compare.R --indir1 <dir> --indir2 <dir> --outdir <dir> [options]
#
#   --indir1/--indir2  the two run output directories                (required)
#   --outdir           output directory for PNGs + GIF               (required)
#   --label1/--label2  panel titles (default "Run 1"/"Run 2")
#   --title            figure super-title (default "WindNinja 10 m wind comparison")
#   --name             output basename (default "compare");
#                      PNGs <name>_<YYYYMMDD_HHMM>.png, GIF <name>_diff.gif
#   --delay            GIF frame delay in 1/100 s (default 200 = 2 s)
#   --tz               time zone label for titles (default "HST")
#   --coast            coastline vector file; reprojected/cropped automatically;
#                      default Hawaii coastline; skipped if missing; "none" disables
#
# Example (Iniki pure ERA5 vs ERA5/ERA5-Land blend):
#   Rscript plot_wind_compare.R \
#     --indir1 .../output/wn_output/iniki_1992_KA               --label1 "Pure ERA5" \
#     --indir2 .../output/wn_output/test_iniki_KA_19920911_blend --label2 "ERA5/ERA5-Land blend" \
#     --outdir .../plots/iniki_1992_compare --title "Hurricane Iniki - WindNinja 10 m wind, Kauai"
#
# Needs R with terra + magick (runs locally, not on the server).

suppressMessages({
  library(terra)
  library(magick)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}
in1 <- get_arg("--indir1"); in2 <- get_arg("--indir2"); out_dir <- get_arg("--outdir")
if (is.null(in1) || is.null(in2) || is.null(out_dir)) {
  stop("usage: Rscript plot_wind_compare.R --indir1 <dir> --indir2 <dir> --outdir <dir> [options]")
}
lab1 <- get_arg("--label1", "Run 1")
lab2 <- get_arg("--label2", "Run 2")
title_1 <- get_arg("--title", "WindNinja 10 m wind comparison")
run_name <- get_arg("--name", "compare")
delay_cs <- as.numeric(get_arg("--delay", "200"))
tz_lab   <- get_arg("--tz", "HST")
coast_in <- get_arg("--coast", "C:/Users/mpluc/gis_data/hi_coastline/Coastline.shp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stamp_of <- function(f) regmatches(basename(f), regexpr("\\d{2}-\\d{2}-\\d{4}_\\d{4}", basename(f)))
vel_in <- function(d) {
  fs <- list.files(d, pattern = "_vel\\.asc$", full.names = TRUE)
  setNames(fs, vapply(fs, stamp_of, character(1)))
}
v1 <- vel_in(in1); v2 <- vel_in(in2)
common <- sort(intersect(names(v1), names(v2)))
if (length(common) == 0) stop("no common timestamps between the two directories")
tt <- as.POSIXct(common, format = "%m-%d-%Y_%H%M", tz = "UTC")
common <- common[order(tt)]; tt <- sort(tt)
cat(sprintf("%d common hours (%s .. %s %s); run1-only: %d, run2-only: %d\n",
            length(common), format(tt[1], "%Y-%m-%d %H:%M"),
            format(tt[length(common)], "%Y-%m-%d %H:%M"), tz_lab,
            length(setdiff(names(v1), names(v2))), length(setdiff(names(v2), names(v1)))))

# ---- shared scales across both runs and all hours --------------------------
vmax <- 0; dmax <- 0
pairs <- list()
for (k in seq_along(common)) {
  s <- common[k]
  r1 <- rast(v1[[s]]); r2 <- rast(v2[[s]])
  stopifnot(ncol(r1) == ncol(r2), nrow(r1) == nrow(r2))
  pairs[[s]] <- list(r1 = r1, r2 = r2, dif = r2 - r1)
  vmax <- max(vmax, values(r1), values(r2))
  dmax <- max(dmax, abs(values(pairs[[s]]$dif)))
}
vmax <- ceiling(vmax); dmax <- ceiling(dmax * 10) / 10
cat(sprintf("shared speed scale 0..%d m/s; difference scale +/-%.1f m/s\n", vmax, dmax))

pal <- rev(rainbow(100, end = 0.8))       # house-style ramp for speed
pal_div <- hcl.colors(101, "Blue-Red 2")  # diverging for the difference, white at 0

coast <- NULL
if (tolower(coast_in) != "none" && file.exists(coast_in)) {
  r0 <- pairs[[1]]$r1
  coast <- tryCatch(crop(project(vect(coast_in), crs(r0)), ext(r0)),
                    error = function(e) NULL)
}

# ---- one 3-panel figure per hour -------------------------------------------
png_files <- character(0)
for (k in seq_along(common)) {
  s <- common[k]; p <- pairs[[s]]
  stamp <- paste(format(tt[k], "%Y-%m-%d %H:%M"), tz_lab)
  png_name <- file.path(out_dir, format(tt[k], paste0(run_name, "_%Y%m%d_%H%M.png")))
  png(png_name, width = 2100, height = 800, res = 140)
  par(mfrow = c(1, 3), mar = c(1, 1, 3, 5), oma = c(0, 0, 3.2, 0))
  plot(p$r1, col = pal, range = c(0, vmax), axes = FALSE, mar = NA,
       main = paste0(lab1, " (m/s)"), cex.main = 1.2, plg = list(title = "m/s"))
  if (!is.null(coast)) lines(coast, col = "white", lwd = 2)
  plot(p$r2, col = pal, range = c(0, vmax), axes = FALSE, mar = NA,
       main = paste0(lab2, " (m/s)"), cex.main = 1.2, plg = list(title = "m/s"))
  if (!is.null(coast)) lines(coast, col = "white", lwd = 2)
  plot(p$dif, col = pal_div, range = c(-dmax, dmax), axes = FALSE, mar = NA,
       main = sprintf("Difference: %s - %s (m/s)", lab2, lab1),
       cex.main = 1.2, plg = list(title = "m/s"))
  if (!is.null(coast)) lines(coast, col = "black", lwd = 2)
  sbar(10000, xy = "bottomleft", type = "bar", divs = 2, below = "m", cex = 0.8)
  mtext(sprintf("%s - %s (speed panels share 0..%d m/s)", title_1, stamp, vmax),
        outer = TRUE, cex = 1.1, font = 2)
  dev.off()
  png_files <- c(png_files, png_name)
  cat(sprintf("wrote %s (diff min %.1f, max %.1f, mean %.2f m/s)\n",
              basename(png_name), min(values(p$dif)), max(values(p$dif)),
              mean(values(p$dif))))
}

# ---- animate ---------------------------------------------------------------
gif_path <- file.path(out_dir, paste0(run_name, "_diff.gif"))
gif <- image_animate(image_join(image_read(png_files)), delay = delay_cs)
image_write(gif, gif_path)
cat("wrote", gif_path, "\n")
