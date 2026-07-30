#!/usr/bin/env Rscript
# era5_to_wrfout.R
#
# Convert an ERA5 single-level GRIB file (10u, 10v, 2t, tcc) into a NetCDF file
# that mimics a WRF-ARW wrfout surface file, for use as a WindNinja
# wxModelInitialization forecast input.
#
# Design constraints (verified against WindNinja source, see
# era5_to_windninja_wrf_mimic.md):
#   - MAP_PROJ must be a projected CRS (lat/long = MAP_PROJ 6 is rejected).
#     Mercator (MAP_PROJ = 3) is used: grid north == true north, so ERA5's
#     earth-relative winds are exactly grid-relative (no rotation error).
#   - Georeferencing comes ONLY from CEN_LAT/CEN_LON + DX/DY + raster size;
#     the grid must be uniformly spaced in the projected CRS.
#   - T2 in Kelvin (range check 180-340), QCLOUD as percent 0-100 (tcc*100;
#     WindNinja divides by 100 internally), U10/V10 in m/s.
#   - char Times(Time, DateStrLen=19), strings "YYYY-MM-DD_HH:MM:SS" (UTC).
#   - TITLE global attribute must contain "WRF".
#   - No 3-D U/V/W variables (they would route the file to the WRF-3D reader).
#   - No NA/nodata cells anywhere (range checks scan every pixel).
#
# Usage:
#   Rscript era5_to_wrfout.R --input <era5.grib> --output <wrfout.nc> \
#       --start "2023-08-08 00:00" --stop "2023-08-08 23:00" [--res 27000]
#
#   --start / --stop  UTC datetimes (YYYY-MM-DD HH:MM); timesteps in [start, stop]
#                     inclusive are written. Defaults: full range in the file.
#   --res             target grid spacing in Mercator meters (default 27000,
#                     ~ERA5 native 0.25 deg).

suppressMessages({
  library(terra)
  library(RNetCDF)
})

# ---------------------------------------------------------------- CLI parsing
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 1 && i < length(args)) return(args[i + 1])
  default
}

in_path  <- get_arg("--input")
out_path <- get_arg("--output", "wrfout_era5.nc")
start_s  <- get_arg("--start")
stop_s   <- get_arg("--stop")
res_m    <- as.numeric(get_arg("--res", "27000"))

if (is.null(in_path)) stop("--input <era5.grib> is required")

parse_utc <- function(s, what) {
  t <- as.POSIXct(s, tz = "UTC",
                  tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
                                 "%Y-%m-%d_%H:%M:%S", "%Y-%m-%dT%H:%M:%S",
                                 "%Y-%m-%d"))
  if (is.na(t)) stop(sprintf("could not parse %s datetime '%s'", what, s))
  t
}

# ------------------------------------------------------------- read the GRIB
# Keep native units: without this GDAL converts temperature K -> C.
setGDALconfig("GRIB_NORMALIZE_UNITS", "NO")

r <- rast(in_path)
btime <- time(r)
bname <- names(r)

classify <- function(nm) {
  if (grepl("u wind component", nm, ignore.case = TRUE)) return("U10")
  if (grepl("v wind component", nm, ignore.case = TRUE)) return("V10")
  if (grepl("2 metre temperature", nm, ignore.case = TRUE)) return("T2")
  if (grepl("cloud cover", nm, ignore.case = TRUE)) return("QCLOUD")
  NA_character_
}
bvar <- vapply(bname, classify, character(1), USE.NAMES = FALSE)
if (anyNA(bvar)) {
  stop("unrecognized GRIB bands: ", paste(unique(bname[is.na(bvar)]), collapse = "; "))
}

t0 <- if (is.null(start_s)) min(btime) else parse_utc(start_s, "--start")
t1 <- if (is.null(stop_s))  max(btime) else parse_utc(stop_s,  "--stop")
if (t1 < t0) stop("--stop is before --start")

sel_times <- sort(unique(btime[btime >= t0 & btime <= t1]))
nt <- length(sel_times)
if (nt == 0) {
  stop(sprintf("no timesteps in [%s, %s]; file covers %s to %s UTC",
               format(t0), format(t1), format(min(btime)), format(max(btime))))
}

vars <- c("U10", "V10", "T2", "QCLOUD")
layer_idx <- lapply(vars, function(v) {
  idx <- vapply(sel_times, function(tt) {
    k <- which(bvar == v & btime == tt)
    if (length(k) != 1) stop(sprintf("expected exactly 1 band for %s at %s, found %d",
                                     v, format(tt), length(k)))
    k
  }, integer(1))
  idx
})
names(layer_idx) <- vars

cat(sprintf("Input: %s\n  %d bands, %s to %s UTC\nSelected %d timesteps: %s to %s UTC\n",
            in_path, nlyr(r), format(min(btime)), format(max(btime)),
            nt, format(sel_times[1]), format(sel_times[nt])))

# ------------------------------------------------- build target Mercator grid
# Shrink to the cell-center extent so bilinear warping never produces edge NAs.
e_src <- ext(r)
half  <- res(r) / 2
e_ctr <- ext(xmin(e_src) + half[1], xmax(e_src) - half[1],
             ymin(e_src) + half[2], ymax(e_src) - half[2])

lon0 <- round((xmin(e_src) + xmax(e_src)) / 2, 2)
crs_merc <- sprintf("+proj=merc +lon_0=%.2f +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", lon0)

e_prj <- project(e_ctr, from = crs(r), to = crs_merc)
cx <- (xmin(e_prj) + xmax(e_prj)) / 2
cy <- (ymin(e_prj) + ymax(e_prj)) / 2
nx <- floor((xmax(e_prj) - xmin(e_prj)) / res_m)
ny <- floor((ymax(e_prj) - ymin(e_prj)) / res_m)
if (nx < 4 || ny < 4) stop("target grid too small; reduce --res")

e_tgt <- ext(cx - nx / 2 * res_m, cx + nx / 2 * res_m,
             cy - ny / 2 * res_m, cy + ny / 2 * res_m)
tgt <- rast(e_tgt, ncols = nx, nrows = ny, crs = crs_merc)

# CEN_LAT/CEN_LON: the grid center back in lat/long (WindNinja reprojects this
# point and anchors the geotransform on it).
cen <- crds(project(vect(cbind(cx, cy), crs = crs_merc), "EPSG:4326"))
cen_lon <- cen[1, 1]; cen_lat <- cen[1, 2]

cat(sprintf("Target grid: Mercator lon_0=%.2f, %d x %d cells @ %.0f m, center %.4f, %.4f\n",
            lon0, nx, ny, res_m, cen_lat, cen_lon))

# ------------------------------------------------------------- warp each var
warped <- lapply(vars, function(v) {
  w <- project(r[[layer_idx[[v]]]], tgt, method = "bilinear")
  if (any(is.na(values(w)))) stop("NA cells after warping ", v, " - grid extends beyond source data")
  w
})
names(warped) <- vars

# Units fixups ---------------------------------------------------------------
# T2: must be Kelvin. If GDAL normalized to Celsius anyway, convert back.
t2max <- max(values(warped$T2))
if (t2max < 200) {
  cat("T2 appears to be Celsius (max", round(t2max, 1), ") - converting to Kelvin\n")
  warped$T2 <- warped$T2 + 273.15
}
# QCLOUD: ERA5 tcc is 0-1 -> percent 0-100 (WindNinja divides by 100 on ingest).
warped$QCLOUD <- clamp(warped$QCLOUD * 100, 0, 100)

# WindNinja range checks (fail here rather than inside WindNinja) ------------
rng <- function(x) range(values(x))
tr <- rng(warped$T2); qr <- rng(warped$QCLOUD)
ur <- rng(warped$U10); vr <- rng(warped$V10)
cat(sprintf("Ranges:  T2 %.1f..%.1f K | QCLOUD %.1f..%.1f %% | U10 %.1f..%.1f | V10 %.1f..%.1f m/s\n",
            tr[1], tr[2], qr[1], qr[2], ur[1], ur[2], vr[1], vr[2]))
if (tr[1] < 180 || tr[2] > 340) stop("T2 outside WindNinja's accepted range 180-340 K")
if (qr[1] < -0.0001 || qr[2] > 100) stop("QCLOUD outside WindNinja's accepted range 0-100")

# --------------------------------------------------------------- write NetCDF
# Array layout: RNetCDF wants dims fastest-varying first -> (west_east,
# south_north, Time). WRF convention stores row 0 = southernmost (bottom-up;
# GDAL's netCDF driver assumes bottom-up and flips on read). terra arrays are
# north-first, so transpose and reverse the y axis.
to_wrf_array <- function(w) {
  a <- as.array(w)                          # [row = north-first y, col = x, t]
  a <- aperm(a, c(2, 1, 3))                 # [x, y north-first, t]
  a[, dim(a)[2]:1, , drop = FALSE]          # [x, y south-first, t]
}

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
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
def_field("U10", "U at 10 M (ERA5 10m u-component, earth-relative)", "m s-1")
def_field("V10", "V at 10 M (ERA5 10m v-component, earth-relative)", "m s-1")
def_field("T2", "TEMP at 2 M (ERA5 2m temperature)", "K")
def_field("QCLOUD", "Total cloud cover as percent (ERA5 tcc*100); WindNinja ingests this band as percent 0-100", "%")

# XLAT/XLONG for human inspection (WindNinja never reads them; no
# 'coordinates' attributes are set so GDAL will not treat them as
# geolocation arrays).
var.def.nc(nc, "XLAT", "NC_FLOAT", c("west_east", "south_north"))
att.put.nc(nc, "XLAT", "description", "NC_CHAR", "LATITUDE, SOUTH IS NEGATIVE")
att.put.nc(nc, "XLAT", "units", "NC_CHAR", "degree_north")
var.def.nc(nc, "XLONG", "NC_FLOAT", c("west_east", "south_north"))
att.put.nc(nc, "XLONG", "description", "NC_CHAR", "LONGITUDE, WEST IS NEGATIVE")
att.put.nc(nc, "XLONG", "units", "NC_CHAR", "degree_east")

# Global attributes (all of these are mandatory for WindNinja's WRF reader)
att.put.nc(nc, "NC_GLOBAL", "TITLE", "NC_CHAR",
           " OUTPUT FROM ERA5 REANALYSIS MIMICKING WRF V3.2 MODEL")
att.put.nc(nc, "NC_GLOBAL", "MAP_PROJ", "NC_INT", 3L)               # Mercator
att.put.nc(nc, "NC_GLOBAL", "DX", "NC_FLOAT", res_m)
att.put.nc(nc, "NC_GLOBAL", "DY", "NC_FLOAT", res_m)
att.put.nc(nc, "NC_GLOBAL", "CEN_LAT", "NC_FLOAT", cen_lat)
att.put.nc(nc, "NC_GLOBAL", "CEN_LON", "NC_FLOAT", cen_lon)
# MOAD_CEN_LAT = 0 keeps WindNinja's Mercator WKT (latitude_of_origin) exactly
# equivalent to the +lat_ts=0 Mercator used for the resampling above.
att.put.nc(nc, "NC_GLOBAL", "MOAD_CEN_LAT", "NC_FLOAT", 0)
att.put.nc(nc, "NC_GLOBAL", "STAND_LON", "NC_FLOAT", lon0)
att.put.nc(nc, "NC_GLOBAL", "TRUELAT1", "NC_FLOAT", 0)              # unused for Mercator
att.put.nc(nc, "NC_GLOBAL", "TRUELAT2", "NC_FLOAT", 0)              # unused for Mercator
att.put.nc(nc, "NC_GLOBAL", "BOTTOM-TOP_GRID_DIMENSION", "NC_INT", 1L)
att.put.nc(nc, "NC_GLOBAL", "START_DATE", "NC_CHAR",
           format(sel_times[1], "%Y-%m-%d_%H:%M:%S"))
att.put.nc(nc, "NC_GLOBAL", "history", "NC_CHAR",
           sprintf("era5_to_wrfout.R: ERA5 %s -> WRF-mimic for WindNinja; %d steps %s..%s UTC",
                   basename(in_path), nt,
                   format(sel_times[1], "%Y-%m-%d %H:%M"),
                   format(sel_times[nt], "%Y-%m-%d %H:%M")))

var.put.nc(nc, "Times", format(sel_times, "%Y-%m-%d_%H:%M:%S"))
for (v in vars) var.put.nc(nc, v, to_wrf_array(warped[[v]]))

xy <- xyFromCell(tgt, seq_len(ncell(tgt)))
ll <- crds(project(vect(xy, crs = crs_merc), "EPSG:4326"))
lat_m <- matrix(ll[, 2], nrow = ny, byrow = TRUE)   # north-first rows
lon_m <- matrix(ll[, 1], nrow = ny, byrow = TRUE)
var.put.nc(nc, "XLAT",  t(lat_m)[, ny:1])           # -> [x, y south-first]
var.put.nc(nc, "XLONG", t(lon_m)[, ny:1])

close.nc(nc)
cat("Wrote", out_path, "\n")

# ------------------------------------------------------------------ verify
# 1) Raw storage must follow the WRF convention: south_north index 0 =
#    southernmost row (bottom-up), exactly like a genuine wrfout. GDAL treats
#    a real wrfout and this mimic identically (verified empirically: with no
#    coordinate variables the netCDF driver presents storage order for both),
#    so matching the real file's convention guarantees WindNinja handles this
#    file exactly as it handles real WRF output.
ncv <- open.nc(out_path)
xl_s <- var.get.nc(ncv, "XLAT", start = c(1, 1),  count = c(1, 1))
xl_n <- var.get.nc(ncv, "XLAT", start = c(1, ny), count = c(1, 1))
t_back <- var.get.nc(ncv, "T2", start = c(1, 1, 1), count = c(nx, ny, 1))
times_back <- var.get.nc(ncv, "Times")
close.nc(ncv)
if (xl_s >= xl_n) stop("Verification FAILED: storage is not bottom-up (south_north index 0 must be southernmost)")

# 2) Value round-trip: what we wrote must equal the warped grid (modulo the
#    known transpose/flip between terra's north-first layout and WRF storage).
m_orig <- to_wrf_array(warped$T2)[, , 1]
if (max(abs(t_back - m_orig)) > 0.01) stop("Verification FAILED: T2 round-trip values do not match")

# 3) Times strings exact format
if (!all(nchar(times_back) == 19) ||
    !all(times_back == format(sel_times, "%Y-%m-%d_%H:%M:%S"))) {
  stop("Verification FAILED: Times strings malformed")
}

# 4) GDAL can open every variable as a subdataset with the full band count
for (v in vars) {
  g <- rast(sprintf('NETCDF:"%s":%s', out_path, v))
  if (nlyr(g) != nt) stop("Verification FAILED: GDAL sees ", nlyr(g), " bands for ", v)
}

cat("Verification: bottom-up storage, value round-trip, Times format, GDAL subdatasets - all OK\n")
cat(sprintf("Done: %d timesteps, %d x %d grid, %s\n", nt, nx, ny, out_path))
