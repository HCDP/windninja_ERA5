#!/usr/bin/env Rscript
# regrid_era5_to_land.R
#
# 1) Nearest-neighbor resample every timestep of every ERA5 variable
#    (10u, 10v, 2t, tcc) onto the EXACT grid of an ERA5-Land GRIB (0.1 deg).
# 2) Blend: ERA5-Land values (10u/10v/2t) on their native land cells,
#    nearest-neighbor ERA5 filling the ocean/missing cells and supplying tcc
#    (tcc does not exist in ERA5-Land).
#
# Outputs two self-describing NetCDF files (variables as (time, lat, lon),
# CF-style coordinates, UTC time axis) plus a land_source mask in the blend:
#   era5_nn_0.1deg_<tag>.nc         resampled ERA5 only
#   era5land_blend_0.1deg_<tag>.nc  blended product
#
# Usage:
#   Rscript regrid_era5_to_land.R [--era5 f.grib] [--land f.grib]
#                                 [--outdir dir] [--tag 1992_09]
#
# Defaults are the Hurricane Iniki files.

suppressMessages({
  library(terra)
  library(RNetCDF)
})

proj_root <- "C:/Users/mpluc/claude_proj/er5_to_WN"
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}
era5_path <- get_arg("--era5", file.path(proj_root, "data/era5/1c0c774fe812285e8d59bc1dd2f53f70.grib"))
land_path <- get_arg("--land", file.path(proj_root, "data/era5/3ae6404a1d1aed4b54bb6854951deda5.grib"))
out_dir   <- get_arg("--outdir", file.path(proj_root, "data/era5_regrid/1992_09"))
tag       <- get_arg("--tag", "1992_09")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setGDALconfig("GRIB_NORMALIZE_UNITS", "NO")

classify <- function(nm) {
  if (grepl("u wind component", nm, ignore.case = TRUE)) return("u10")
  if (grepl("v wind component", nm, ignore.case = TRUE)) return("v10")
  if (grepl("2 metre temperature", nm, ignore.case = TRUE)) return("t2")
  if (grepl("cloud cover", nm, ignore.case = TRUE)) return("tcc")
  NA_character_
}

# ---- load both sources, classify bands ------------------------------------
era5 <- rast(era5_path)
land <- rast(land_path)
e_var <- vapply(names(era5), classify, character(1), USE.NAMES = FALSE)
l_var <- vapply(names(land), classify, character(1), USE.NAMES = FALSE)
if (anyNA(e_var)) stop("unrecognized band(s) in ERA5 file")
if (anyNA(l_var)) stop("unrecognized band(s) in ERA5-Land file")
e_time <- time(era5); l_time <- time(land)

times <- sort(unique(e_time))
nt <- length(times)
if (!setequal(times, unique(l_time))) stop("ERA5 and ERA5-Land timesteps differ")
cat(sprintf("%d hourly timesteps, %s to %s UTC\n",
            nt, format(times[1]), format(times[nt])))

nx <- ncol(land); ny <- nrow(land)
cat(sprintf("target grid (ERA5-Land): %d x %d @ %.3g deg\n", nx, ny, res(land)[1]))

vars_nn    <- c("u10", "v10", "t2", "tcc")
vars_blend <- c("u10", "v10", "t2")     # tcc has no ERA5-Land counterpart

# band pickers
pick <- function(r, vlab, want, tt) {
  k <- which(vlab == want & time(r) == tt)
  if (length(k) != 1) stop("expected 1 band for ", want, " at ", format(tt))
  r[[k]]
}

# ---- resample + blend, hour by hour ---------------------------------------
# arrays [lon, lat(north-first), time] to match the NetCDF layout below
mk_arr <- function() array(NA_real_, dim = c(nx, ny, nt))
nn <- setNames(lapply(vars_nn, function(v) mk_arr()), vars_nn)
bl <- setNames(lapply(vars_nn, function(v) mk_arr()), vars_nn)

as_mat <- function(r) t(matrix(values(r)[, 1], ncol = nx, byrow = TRUE))  # [lon, lat N-first]

land_mask <- NULL
for (i in seq_len(nt)) {
  tt <- times[i]
  for (v in vars_nn) {
    r_nn <- resample(pick(era5, e_var, v, tt), land, method = "near")
    m_nn <- as_mat(r_nn)
    nn[[v]][, , i] <- m_nn
    if (v %in% vars_blend) {
      m_land <- as_mat(pick(land, l_var, v, tt))
      bl[[v]][, , i] <- ifelse(is.na(m_land), m_nn, m_land)
      msk <- !is.na(m_land)
      if (is.null(land_mask)) land_mask <- msk
      else if (!identical(land_mask, msk)) cat("NOTE: ERA5-Land valid-cell pattern varies at", format(tt), "\n")
    } else {
      bl[[v]][, , i] <- m_nn   # tcc: ERA5 everywhere
    }
  }
}
cat(sprintf("ERA5-Land supplies %d of %d cells (%.0f%%); ERA5 fills the rest + tcc\n",
            sum(land_mask), nx * ny, 100 * sum(land_mask) / (nx * ny)))

# ---- verification ----------------------------------------------------------
stopifnot(!anyNA(nn$u10), !anyNA(nn$tcc), !anyNA(bl$u10), !anyNA(bl$t2), !anyNA(bl$tcc))
# NN correctness: sampled hours, output cell == source cell containing it
xy <- xyFromCell(land, seq_len(ncell(land)))
for (i in c(1, nt)) {
  src <- pick(era5, e_var, "t2", times[i])
  v_chk <- terra::extract(src, xy, method = "simple")[, 1]
  stopifnot(max(abs(as.vector(nn$t2[, , i]) - v_chk)) < 1e-4)
}
# blend correctness at one hour: land cells = ERA5-Land, others = NN
i <- 1
m_land <- as_mat(pick(land, l_var, "t2", times[i]))
stopifnot(all(abs(bl$t2[, , i][land_mask] - m_land[land_mask]) < 1e-6))
stopifnot(all(abs(bl$t2[, , i][!land_mask] - nn$t2[, , i][!land_mask]) < 1e-6))
cat("verification: no NAs; NN cells match source; blend = ERA5-Land on land, ERA5 elsewhere\n")

# ---- write NetCDF ----------------------------------------------------------
lon <- sort(unique(xy[, 1]))                 # ascending
lat <- sort(unique(xy[, 2]), decreasing = TRUE)  # north-first, matches array layout
time_num <- as.numeric(times)                # seconds since epoch, UTC

meta <- list(
  u10 = c("10 metre u wind component", "m s-1"),
  v10 = c("10 metre v wind component", "m s-1"),
  t2  = c("2 metre temperature", "K"),
  tcc = c("total cloud cover", "1"))

write_nc <- function(path, dat, title, extra_src_mask = NULL) {
  nc <- create.nc(path, format = "netcdf4")
  dim.def.nc(nc, "time", unlim = TRUE)
  dim.def.nc(nc, "latitude", length(lat))
  dim.def.nc(nc, "longitude", length(lon))
  var.def.nc(nc, "time", "NC_DOUBLE", "time")
  att.put.nc(nc, "time", "units", "NC_CHAR", "seconds since 1970-01-01 00:00:00")
  att.put.nc(nc, "time", "calendar", "NC_CHAR", "standard")
  var.def.nc(nc, "latitude", "NC_DOUBLE", "latitude")
  att.put.nc(nc, "latitude", "units", "NC_CHAR", "degrees_north")
  var.def.nc(nc, "longitude", "NC_DOUBLE", "longitude")
  att.put.nc(nc, "longitude", "units", "NC_CHAR", "degrees_east")
  for (v in names(dat)) {
    var.def.nc(nc, v, "NC_FLOAT", c("longitude", "latitude", "time"))
    att.put.nc(nc, v, "long_name", "NC_CHAR", meta[[v]][1])
    att.put.nc(nc, v, "units", "NC_CHAR", meta[[v]][2])
  }
  if (!is.null(extra_src_mask)) {
    var.def.nc(nc, "land_source", "NC_BYTE", c("longitude", "latitude"))
    att.put.nc(nc, "land_source", "long_name", "NC_CHAR",
               "1 = value from ERA5-Land (u10/v10/t2), 0 = nearest-neighbor ERA5; tcc is ERA5 everywhere")
  }
  att.put.nc(nc, "NC_GLOBAL", "title", "NC_CHAR", title)
  att.put.nc(nc, "NC_GLOBAL", "source_era5", "NC_CHAR", basename(era5_path))
  att.put.nc(nc, "NC_GLOBAL", "source_era5_land", "NC_CHAR", basename(land_path))
  att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR",
             "regrid_era5_to_land.R: nearest-neighbor ERA5 (0.25 deg) onto the ERA5-Land 0.1 deg grid")
  var.put.nc(nc, "time", time_num)
  var.put.nc(nc, "latitude", lat)
  var.put.nc(nc, "longitude", lon)
  for (v in names(dat)) var.put.nc(nc, v, dat[[v]])
  if (!is.null(extra_src_mask)) var.put.nc(nc, "land_source", ifelse(extra_src_mask, 1L, 0L))
  close.nc(nc)
  cat("wrote", path, sprintf("(%.0f KB)\n", file.size(path) / 1024))
}

write_nc(file.path(out_dir, paste0("era5_nn_0.1deg_", tag, ".nc")), nn,
         "ERA5 surface fields nearest-neighbor resampled to the ERA5-Land 0.1 deg grid")
write_nc(file.path(out_dir, paste0("era5land_blend_0.1deg_", tag, ".nc")), bl,
         "Blend: ERA5-Land (u10/v10/t2) on land cells; nearest-neighbor ERA5 elsewhere and for tcc",
         extra_src_mask = land_mask)

# ---- smoothed blend --------------------------------------------------------
# Bilinear smoothing of the blend ON THE SAME GRID: a straight bilinear
# resample onto an identical grid is a no-op (cell centers align), so the
# smoothing is done with the bilinear tent kernel - mathematically identical
# to bilinear-resampling onto a half-cell-shifted grid and back:
#   K = [1 2 1; 2 4 2; 1 2 1] / 16
# Edge cells (where the kernel would run off the grid) keep their blend
# values. Afterwards the ERA5-Land land pixels (u10/v10/t2) are restored to
# their ORIGINAL values, so smoothing only alters the ERA5-filled cells
# (tcc has no ERA5-Land source and stays smoothed everywhere).
K <- matrix(c(1, 2, 1, 2, 4, 2, 1, 2, 1), 3, 3) / 16
tmpl_r <- rast(land[[1]])

sm <- setNames(lapply(vars_nn, function(v) mk_arr()), vars_nn)
for (i in seq_len(nt)) {
  for (v in vars_nn) {
    r <- tmpl_r
    values(r) <- as.vector(bl[[v]][, , i])
    r_sm <- focal(r, w = K, fun = sum, na.rm = FALSE)
    m_sm <- as_mat(r_sm)
    m_bl <- bl[[v]][, , i]
    m_sm[is.na(m_sm)] <- m_bl[is.na(m_sm)]         # edges keep blend values
    if (v %in% vars_blend) m_sm[land_mask] <- m_bl[land_mask]  # restore ERA5-Land pixels
    sm[[v]][, , i] <- m_sm
  }
}

# verify: no NAs; land pixels exactly equal the blend (= ERA5-Land originals);
# interior ocean cells actually changed
stopifnot(!anyNA(sm$u10), !anyNA(sm$t2), !anyNA(sm$tcc))
stopifnot(all(abs(sm$t2[, , 1][land_mask] - bl$t2[, , 1][land_mask]) < 1e-9))
interior <- matrix(FALSE, nx, ny); interior[2:(nx - 1), 2:(ny - 1)] <- TRUE
chg <- interior & !land_mask
stopifnot(mean(abs(sm$t2[, , 1][chg] - bl$t2[, , 1][chg])) > 0)
cat("verification (smoothed): no NAs; land pixels = ERA5-Land originals; ocean cells smoothed\n")

write_nc(file.path(out_dir, paste0("era5land_blend_smooth_0.1deg_", tag, ".nc")), sm,
         "Smoothed blend: bilinear tent-kernel smoothing of the blend, then ERA5-Land land pixels (u10/v10/t2) restored to original values; tcc smoothed everywhere",
         extra_src_mask = land_mask)
cat("done\n")
