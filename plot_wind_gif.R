#!/usr/bin/env Rscript
# plot_wind_gif.R
#
# Per-hour PNG maps of WindNinja 10 m wind speed (fill) + direction (arrows)
# from any directory of WindNinja ASCII outputs, all on a COMMON speed scale
# (global max across hours), then an animated GIF.
#
# Expects WindNinja CLI ascii output pairs in --indir:
#   <dem>_<MM-DD-YYYY>_<HHMM>_<res>m_vel.asc   (speed, m/s)
#   <dem>_<MM-DD-YYYY>_<HHMM>_<res>m_ang.asc   (met direction, deg FROM)
# Timestamps in the filenames are local (whatever time_zone the run used).
#
# Usage:
#   Rscript plot_wind_gif.R --indir <dir> --outdir <dir> [options]
#
#   --indir    directory with *_vel.asc / *_ang.asc pairs        (required)
#   --outdir   output directory for PNGs + GIF                   (required)
#   --title    plot title line 1 (default "WindNinja 10 m wind")
#   --name     basename for outputs (default: basename of indir);
#              PNGs: <name>_<YYYYMMDD_HHMM>.png, GIF: <name>_wind.gif
#   --delay    GIF frame delay in 1/100 s (default 100 = 1 s per hour)
#   --fact     arrow aggregation factor in cells (default 6; at 250 m
#              resolution that is one arrow per 1.5 km)
#   --tz       time zone label for the title timestamps (default "HST")
#   --coast    coastline vector file (shp/gpkg) drawn over the map; reprojected
#              to the raster CRS automatically. Default is the local Hawaii
#              coastline (C:/Users/mpluc/gis_data/hi_coastline/Coastline.shp);
#              silently skipped if the file does not exist. --coast none disables.
#
# Example (Hurricane Iniki, Kauai):
#   Rscript plot_wind_gif.R \
#     --indir  C:/Users/mpluc/claude_proj/er5_to_WN/output/wn_output/iniki_1992_KA \
#     --outdir C:/Users/mpluc/claude_proj/er5_to_WN/plots/iniki_1992 \
#     --title  "Hurricane Iniki - WindNinja 10 m wind, Kauai" --name iniki_1992_KA
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

in_dir  <- get_arg("--indir")
out_dir <- get_arg("--outdir")
if (is.null(in_dir) || is.null(out_dir)) {
  stop("usage: Rscript plot_wind_gif.R --indir <dir> --outdir <dir> [--title t] [--name n] [--delay 100] [--fact 15] [--tz HST]")
}
title_1 <- get_arg("--title", "WindNinja 10 m wind")
run_name <- get_arg("--name", basename(normalizePath(in_dir, mustWork = FALSE)))
delay_cs <- as.numeric(get_arg("--delay", "100"))
fact     <- as.integer(get_arg("--fact", "6"))
tz_lab   <- get_arg("--tz", "HST")
coast_in <- get_arg("--coast", "C:/Users/mpluc/gis_data/hi_coastline/Coastline.shp")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

vel_files <- sort(list.files(in_dir, pattern = "_vel\\.asc$", full.names = TRUE))
if (length(vel_files) == 0) stop("no *_vel.asc files found in ", in_dir)

# pair each vel with its ang; parse the local timestamp out of the filename
runs <- lapply(vel_files, function(vf) {
  af <- sub("_vel\\.asc$", "_ang.asc", vf)
  if (!file.exists(af)) stop("missing direction file for ", basename(vf))
  m  <- regmatches(basename(vf), regexpr("\\d{2}-\\d{2}-\\d{4}_\\d{4}", basename(vf)))
  if (length(m) == 0) stop("no MM-DD-YYYY_HHMM timestamp in filename: ", basename(vf))
  tt <- as.POSIXct(m, format = "%m-%d-%Y_%H%M", tz = "UTC") # tz irrelevant, label only
  list(vel = vf, ang = af, time = tt)
})
runs <- runs[order(sapply(runs, function(r) r$time))]

# ---- common scale: global speed max across all hours -----------------------
vmax <- max(sapply(runs, function(r) max(values(rast(r$vel)), na.rm = TRUE)))
vmax <- ceiling(vmax)
cat(sprintf("%d frames, global max speed -> shared scale 0..%d m/s\n",
            length(runs), vmax))

pal <- hcl.colors(64, "Viridis")          # sequential, perceptually uniform, CVD-safe
arrow_col <- "white"                       # neutral overlay, readable on dark high-speed fill

# optional coastline, reprojected to the raster CRS and cropped to the map
coast <- NULL
if (tolower(coast_in) != "none" && file.exists(coast_in)) {
  r0 <- rast(runs[[1]]$vel)
  coast <- tryCatch(crop(project(vect(coast_in), crs(r0)), ext(r0)),
                    error = function(e) { cat("coastline skipped:", conditionMessage(e), "\n"); NULL })
} else if (tolower(coast_in) != "none") {
  cat("coastline file not found, skipping:", coast_in, "\n")
}

# ---- one map per hour ------------------------------------------------------
png_files <- character(0)
for (r in runs) {
  spd <- rast(r$vel)
  ang <- rast(r$ang)

  # arrow field: aggregate direction via u/v components (never average angles)
  spd_c <- aggregate(spd, fact = fact, fun = "mean", na.rm = TRUE)
  u_c <- aggregate(-spd * sin(ang * pi / 180), fact = fact, fun = "mean", na.rm = TRUE)
  v_c <- aggregate(-spd * cos(ang * pi / 180), fact = fact, fun = "mean", na.rm = TRUE)

  xy  <- crds(spd_c, na.rm = FALSE)
  uu  <- values(u_c)[, 1]; vv <- values(v_c)[, 1]; ss <- values(spd_c)[, 1]
  ok  <- is.finite(uu) & is.finite(vv) & is.finite(ss) & ss > 0
  # arrow length proportional to speed on the SHARED scale
  max_len <- fact * res(spd)[1] * 0.9
  sc <- (ss[ok] / vmax) * max_len / sqrt(uu[ok]^2 + vv[ok]^2)

  stamp    <- paste(format(r$time, "%Y-%m-%d %H:%M"), tz_lab)
  png_name <- file.path(out_dir, format(r$time, paste0(run_name, "_%Y%m%d_%H%M.png")))

  png(png_name, width = 1100, height = 900, res = 130)
  par(mar = c(1.5, 1.5, 3.5, 1))
  plot(spd, col = pal, range = c(0, vmax), axes = FALSE, mar = NA,
       main = paste0(title_1, "\n", stamp,
                     "   (fill: speed m/s, arrows: flow direction)"),
       cex.main = 1.05, plg = list(title = "m/s"))
  arrows(xy[ok, 1], xy[ok, 2],
         xy[ok, 1] + uu[ok] * sc, xy[ok, 2] + vv[ok] * sc,
         length = 0.025, lwd = 0.8, col = arrow_col)
  if (!is.null(coast)) lines(coast, col = "black", lwd = 1.2)
  sbar(10000, xy = "bottomleft", type = "bar", divs = 2, below = "m", cex = 0.7)
  dev.off()

  png_files <- c(png_files, png_name)
  cat("wrote", basename(png_name), sprintf(" (frame max %.1f m/s)\n",
      max(values(spd), na.rm = TRUE)))
}

# ---- animate ---------------------------------------------------------------
gif_path <- file.path(out_dir, paste0(run_name, "_wind.gif"))
frames <- image_read(png_files)
gif <- image_animate(image_join(frames), delay = delay_cs)
image_write(gif, gif_path)
cat("wrote", gif_path, "\n")
