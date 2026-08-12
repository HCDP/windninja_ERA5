#!/usr/bin/env Rscript
# blend_to_wrfout.R
#
# Convert the blended ERA5/ERA5-Land NetCDF (from regrid_era5_to_land.R /
# subset_blend_overlap.R; variables u10, v10, t2, tcc on the 0.1 deg grid)
# into SINGLE-TIMESTEP WRF-mimic NetCDF files for WindNinja - one file per
# hour, named by UTC valid time (wrfout_era5_YYYYMMDD_HHMM.nc), drop-in
# compatible with the existing run scripts (point WN_FORECAST_DIR at --outdir).
#
# Same verified wrfout schema as era5_to_wrfout.R (Mercator MAP_PROJ=3,
# bottom-up rows, QCLOUD = tcc*100 percent, T2 Kelvin, units "m s-1"/"K",
# TITLE containing WRF, mandatory globals). Differences:
#   - input is the blended NetCDF, not an ERA5 GRIB
#   - NEAREST-NEIGHBOR resampling onto the Mercator grid, preserving the
#     exact blended pixel values (including the restored ERA5-Land land cells)
#   - default resolution 11000 m (~the blend's native 0.1 deg)
#
# Usage:
#   Rscript blend_to_wrfout.R [--input blend.nc] [--outdir dir] [--res 11000]
#                             [--start "YYYY-MM-DD HH:MM"] [--stop "..."]  (UTC)

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
in_path <- get_arg("--input", file.path(proj_root, "data/era5_regrid/1992_09/era5land_blend_smooth_0.1deg_1992_09_overlap.nc"))
out_dir <- get_arg("--outdir", file.path(proj_root, "output/er5_wrf_hourly/1992_09_blend"))
res_m   <- as.numeric(get_arg("--res", "11000"))
start_s <- get_arg("--start"); stop_s <- get_arg("--stop")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

parse_utc <- function(s) as.POSIXct(s, tz = "UTC",
  tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d"))

# ---- read the blend --------------------------------------------------------
vars <- c("u10", "v10", "t2", "tcc")
stk <- lapply(vars, function(v) {
  r <- rast(in_path, subds = v)
  if (crs(r) == "") crs(r) <- "EPSG:4326"
  r
})
names(stk) <- vars
times <- as.POSIXct(time(stk$u10), tz = "UTC")
sel <- rep(TRUE, length(times))
if (!is.null(start_s)) sel <- sel & times >= parse_utc(start_s)
if (!is.null(stop_s))  sel <- sel & times <= parse_utc(stop_s)
idx <- which(sel)
if (length(idx) == 0) stop("no timesteps in requested window")
cat(sprintf("input: %s\n  %d timesteps, using %d (%s .. %s UTC)\n",
            basename(in_path), length(times), length(idx),
            format(times[idx[1]]), format(times[idx[length(idx)]])))

# ---- target Mercator grid --------------------------------------------------
r0 <- stk$u10[[1]]
e_src <- ext(r0)
lon0 <- round((xmin(e_src) + xmax(e_src)) / 2, 2)
crs_merc <- sprintf("+proj=merc +lon_0=%.2f +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", lon0)
e_prj <- project(e_src, from = crs(r0), to = crs_merc)
cx <- (xmin(e_prj) + xmax(e_prj)) / 2
cy <- (ymin(e_prj) + ymax(e_prj)) / 2
nx <- floor((xmax(e_prj) - xmin(e_prj)) / res_m)
ny <- floor((ymax(e_prj) - ymin(e_prj)) / res_m)
if (nx < 4 || ny < 4) stop("target grid too small; reduce --res")
e_tgt <- ext(cx - nx / 2 * res_m, cx + nx / 2 * res_m,
             cy - ny / 2 * res_m, cy + ny / 2 * res_m)
tgt <- rast(e_tgt, ncols = nx, nrows = ny, crs = crs_merc)
cen <- crds(project(vect(cbind(cx, cy), crs = crs_merc), "EPSG:4326"))
cen_lon <- cen[1, 1]; cen_lat <- cen[1, 2]
cat(sprintf("target grid: Mercator lon_0=%.2f, %d x %d @ %.0f m, center %.4f, %.4f, resample=near\n",
            lon0, nx, ny, res_m, cen_lat, cen_lon))

# XLAT/XLONG for the target grid (bottom-up storage, [x, y south-first])
xy <- xyFromCell(tgt, seq_len(ncell(tgt)))
ll <- crds(project(vect(xy, crs = crs_merc), "EPSG:4326"))
lat_m <- matrix(ll[, 2], nrow = ny, byrow = TRUE)
lon_m <- matrix(ll[, 1], nrow = ny, byrow = TRUE)
xlat_arr  <- t(lat_m)[, ny:1]
xlong_arr <- t(lon_m)[, ny:1]

to_wrf <- function(w) {              # SpatRaster layer -> [x, y south-first]
  m <- matrix(values(w)[, 1], nrow = ny, byrow = TRUE)
  t(m)[, ny:1]
}

# ---- one mimic file per timestep -------------------------------------------
n_ok <- 0
for (i in idx) {
  tt <- times[i]
  U <- to_wrf(project(stk$u10[[i]], tgt, method = "near"))
  V <- to_wrf(project(stk$v10[[i]], tgt, method = "near"))
  T2 <- to_wrf(project(stk$t2[[i]], tgt, method = "near"))
  QC <- to_wrf(project(stk$tcc[[i]], tgt, method = "near")) * 100
  QC <- pmin(pmax(QC, 0), 100)
  if (anyNA(U) || anyNA(T2) || anyNA(QC)) stop("NA cells after resampling at ", format(tt))
  if (min(T2) < 180 || max(T2) > 340) stop("T2 out of WindNinja range at ", format(tt))

  out_path <- file.path(out_dir, format(tt, "wrfout_era5_%Y%m%d_%H%M.nc"))
  nc <- create.nc(out_path, format = "offset64")
  dim.def.nc(nc, "Time", unlim = TRUE)
  dim.def.nc(nc, "DateStrLen", 19)
  dim.def.nc(nc, "west_east", nx)
  dim.def.nc(nc, "south_north", ny)
  var.def.nc(nc, "Times", "NC_CHAR", c("DateStrLen", "Time"))
  def_field <- function(name, desc, units) {
    var.def.nc(nc, name, "NC_FLOAT", c("west_east", "south_north", "Time"))
    att.put.nc(nc, name, "FieldType", "NC_INT", 104L)
    att.put.nc(nc, name, "MemoryOrder", "NC_CHAR", "XY ")
    att.put.nc(nc, name, "description", "NC_CHAR", desc)
    att.put.nc(nc, name, "units", "NC_CHAR", units)
    att.put.nc(nc, name, "stagger", "NC_CHAR", "")
  }
  def_field("U10", "U at 10 M (ERA5/ERA5-Land blend, earth-relative)", "m s-1")
  def_field("V10", "V at 10 M (ERA5/ERA5-Land blend, earth-relative)", "m s-1")
  def_field("T2", "TEMP at 2 M (ERA5/ERA5-Land blend)", "K")
  def_field("QCLOUD", "Total cloud cover as percent (ERA5 tcc*100); WindNinja ingests this band as percent 0-100", "%")
  var.def.nc(nc, "XLAT", "NC_FLOAT", c("west_east", "south_north"))
  att.put.nc(nc, "XLAT", "description", "NC_CHAR", "LATITUDE, SOUTH IS NEGATIVE")
  att.put.nc(nc, "XLAT", "units", "NC_CHAR", "degree_north")
  var.def.nc(nc, "XLONG", "NC_FLOAT", c("west_east", "south_north"))
  att.put.nc(nc, "XLONG", "description", "NC_CHAR", "LONGITUDE, WEST IS NEGATIVE")
  att.put.nc(nc, "XLONG", "units", "NC_CHAR", "degree_east")
  att.put.nc(nc, "NC_GLOBAL", "TITLE", "NC_CHAR",
             " OUTPUT FROM ERA5/ERA5-LAND BLEND MIMICKING WRF V3.2 MODEL")
  att.put.nc(nc, "NC_GLOBAL", "MAP_PROJ", "NC_INT", 3L)
  att.put.nc(nc, "NC_GLOBAL", "DX", "NC_FLOAT", res_m)
  att.put.nc(nc, "NC_GLOBAL", "DY", "NC_FLOAT", res_m)
  att.put.nc(nc, "NC_GLOBAL", "CEN_LAT", "NC_FLOAT", cen_lat)
  att.put.nc(nc, "NC_GLOBAL", "CEN_LON", "NC_FLOAT", cen_lon)
  att.put.nc(nc, "NC_GLOBAL", "MOAD_CEN_LAT", "NC_FLOAT", 0)
  att.put.nc(nc, "NC_GLOBAL", "STAND_LON", "NC_FLOAT", lon0)
  att.put.nc(nc, "NC_GLOBAL", "TRUELAT1", "NC_FLOAT", 0)
  att.put.nc(nc, "NC_GLOBAL", "TRUELAT2", "NC_FLOAT", 0)
  att.put.nc(nc, "NC_GLOBAL", "BOTTOM-TOP_GRID_DIMENSION", "NC_INT", 1L)
  att.put.nc(nc, "NC_GLOBAL", "START_DATE", "NC_CHAR", format(tt, "%Y-%m-%d_%H:%M:%S"))
  att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR",
             sprintf("blend_to_wrfout.R: %s -> WRF-mimic (nearest-neighbor, %g m)", basename(in_path), res_m))
  var.put.nc(nc, "Times", format(tt, "%Y-%m-%d_%H:%M:%S"))
  var.put.nc(nc, "U10", array(U, dim = c(nx, ny, 1)))
  var.put.nc(nc, "V10", array(V, dim = c(nx, ny, 1)))
  var.put.nc(nc, "T2", array(T2, dim = c(nx, ny, 1)))
  var.put.nc(nc, "QCLOUD", array(QC, dim = c(nx, ny, 1)))
  var.put.nc(nc, "XLAT", xlat_arr)
  var.put.nc(nc, "XLONG", xlong_arr)
  close.nc(nc)
  n_ok <- n_ok + 1
}
cat(sprintf("wrote %d hourly mimic files -> %s\n", n_ok, out_dir))

# ---- verify the first written file -----------------------------------------
f1 <- file.path(out_dir, format(times[idx[1]], "wrfout_era5_%Y%m%d_%H%M.nc"))
ncv <- open.nc(f1)
xl_s <- var.get.nc(ncv, "XLAT", start = c(1, 1), count = c(1, 1))
xl_n <- var.get.nc(ncv, "XLAT", start = c(1, ny), count = c(1, 1))
t_back <- var.get.nc(ncv, "Times")
q_rng <- range(var.get.nc(ncv, "QCLOUD"))
close.nc(ncv)
stopifnot(xl_s < xl_n)                       # bottom-up storage
stopifnot(nchar(t_back) == 19)               # Times format
stopifnot(q_rng[1] >= 0, q_rng[2] <= 100)    # QCLOUD percent range
g <- rast(sprintf('NETCDF:"%s":U10', f1))    # GDAL-openable
stopifnot(nlyr(g) == 1)
cat("verification: bottom-up storage, Times format, QCLOUD 0-100, GDAL subdataset - all OK\n")
