#!/usr/bin/env Rscript
# subset_blend_overlap.R
#
# Subset a blended NetCDF (from regrid_era5_to_land.R) to the OVERLAPPING
# datetime range of its two source GRIBs (ERA5 and ERA5-Land): only timesteps
# present in BOTH sources are kept. When the sources fully overlap (the usual
# case for paired downloads) the output equals the input; the script says so.
#
# Usage:
#   Rscript subset_blend_overlap.R [--smooth f.nc] [--era5 f.grib]
#                                  [--land f.grib] [--out f.nc]
#
# Defaults are the Hurricane Iniki files; --out defaults to the input name
# with "_overlap" inserted before .nc.

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
smooth_path <- get_arg("--smooth", file.path(proj_root, "data/era5_regrid/1992_09/era5land_blend_smooth_0.1deg_1992_09.nc"))
era5_path   <- get_arg("--era5",   file.path(proj_root, "data/era5/1c0c774fe812285e8d59bc1dd2f53f70.grib"))
land_path   <- get_arg("--land",   file.path(proj_root, "data/era5/3ae6404a1d1aed4b54bb6854951deda5.grib"))
out_path    <- get_arg("--out",    sub("\\.nc$", "_overlap.nc", smooth_path))

setGDALconfig("GRIB_NORMALIZE_UNITS", "NO")

# ---- source time axes and their intersection -------------------------------
t_era5 <- sort(unique(as.numeric(time(rast(era5_path)))))
t_land <- sort(unique(as.numeric(time(rast(land_path)))))
t_common <- sort(intersect(t_era5, t_land))
fmt <- function(x) format(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H:%M")
cat(sprintf("ERA5:      %d steps, %s .. %s UTC\n", length(t_era5), fmt(min(t_era5)), fmt(max(t_era5))))
cat(sprintf("ERA5-Land: %d steps, %s .. %s UTC\n", length(t_land), fmt(min(t_land)), fmt(max(t_land))))
if (length(t_common) == 0) stop("sources have NO overlapping timesteps")
cat(sprintf("overlap:   %d steps, %s .. %s UTC\n", length(t_common), fmt(min(t_common)), fmt(max(t_common))))

# ---- read the smoothed file and pick the overlapping timesteps -------------
nc_in <- open.nc(smooth_path)
t_nc <- var.get.nc(nc_in, "time")
keep <- which(t_nc %in% t_common)
cat(sprintf("smoothed file: %d steps; keeping %d, dropping %d\n",
            length(t_nc), length(keep), length(t_nc) - length(keep)))
if (length(keep) == 0) { close.nc(nc_in); stop("no smoothed timesteps fall in the overlap") }
noop <- length(keep) == length(t_nc)
if (noop) cat("NOTE: sources fully overlap -> output equals input (no-op subset)\n")

lat <- var.get.nc(nc_in, "latitude")
lon <- var.get.nc(nc_in, "longitude")
vars <- c("u10", "v10", "t2", "tcc")
dat <- lapply(vars, function(v) var.get.nc(nc_in, v)[, , keep, drop = FALSE])
names(dat) <- vars
get_att <- function(v, a) att.get.nc(nc_in, v, a)
long_names <- vapply(vars, function(v) get_att(v, "long_name"), character(1))
units_v    <- vapply(vars, function(v) get_att(v, "units"), character(1))
src_mask   <- tryCatch(var.get.nc(nc_in, "land_source"), error = function(e) NULL)
g_atts <- lapply(c("title", "source_era5", "source_era5_land", "history"),
                 function(a) tryCatch(att.get.nc(nc_in, "NC_GLOBAL", a), error = function(e) NA))
names(g_atts) <- c("title", "source_era5", "source_era5_land", "history")
close.nc(nc_in)

# ---- write the subset ------------------------------------------------------
nc <- create.nc(out_path, format = "netcdf4")
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
for (v in vars) {
  var.def.nc(nc, v, "NC_FLOAT", c("longitude", "latitude", "time"))
  att.put.nc(nc, v, "long_name", "NC_CHAR", long_names[[v]])
  att.put.nc(nc, v, "units", "NC_CHAR", units_v[[v]])
}
if (!is.null(src_mask)) {
  var.def.nc(nc, "land_source", "NC_BYTE", c("longitude", "latitude"))
  att.put.nc(nc, "land_source", "long_name", "NC_CHAR",
             "1 = value from ERA5-Land (u10/v10/t2), 0 = nearest-neighbor ERA5; tcc is ERA5 everywhere")
}
for (a in names(g_atts)) if (!is.na(g_atts[[a]])) att.put.nc(nc, "NC_GLOBAL", a, "NC_CHAR", g_atts[[a]])
att.put.nc(nc, "NC_GLOBAL", "subset", "NC_CHAR",
           sprintf("subset_blend_overlap.R: kept %d of %d timesteps = ERA5/ERA5-Land overlap %s .. %s UTC",
                   length(keep), length(t_nc), fmt(min(t_common)), fmt(max(t_common))))
var.put.nc(nc, "time", t_nc[keep])
var.put.nc(nc, "latitude", lat)
var.put.nc(nc, "longitude", lon)
for (v in vars) var.put.nc(nc, v, dat[[v]])
if (!is.null(src_mask)) var.put.nc(nc, "land_source", src_mask)
close.nc(nc)

# verify round trip
chk <- open.nc(out_path)
stopifnot(length(var.get.nc(chk, "time")) == length(keep),
          all(var.get.nc(chk, "time") == t_nc[keep]),
          all(abs(var.get.nc(chk, "t2") - dat$t2) < 1e-6))
close.nc(chk)
cat("verification: subset time axis and values round-trip OK\n")
cat("wrote", out_path, sprintf("(%.0f KB)\n", file.size(out_path) / 1024))
