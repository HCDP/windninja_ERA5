==========================================================================
 era5_wrf_TEST - ERA5-driven WindNinja hourly wind workflow (TEST)
==========================================================================

WHAT THIS IS
------------
A test edition of the HCDP hourly wind mapping workflow
(HCDP/windninja_hrly_wind) in which WindNinja is initialized from ERA5
reanalysis instead of weather stations. ERA5 surface fields (10 m u/v
wind, 2 m temperature, total cloud cover) are converted to a NetCDF file
that mimics a WRF-ARW wrfout file, which WindNinja accepts through its
weather-model initialization path (wxModelInitialization) with the
diurnal slope-flow model ENABLED.

Station observations are still fetched every hour, but they are used
ONLY as independent cross-validation truth - never as model input. All
hourly and daily metadata text files report cross-validation statistics
(RMSE, RME, BIAS, R-SQUARE for wind speed and direction).

Everything is quarantined from the operational point-initialization
workflow: this code lives in /home/wn1/wn_codes/era5_wrf_TEST and writes
only to *_era5 data trees (wn_era5_input, wn_era5_out,
windProducts_era5). Operational products are never touched.

DEPLOYED LOCATION
-----------------
  Code:      /home/wn1/wn_codes/era5_wrf_TEST/          (this directory)
  Forecast:  /home/wn1/data/input/er5_wrf/hourly/wrfout_era5_<YYYYMMDD_HHMM>.nc
             ONE SINGLE-TIMESTEP file per hour, named by UTC valid time
             (made locally by era5_to_wrfout.R with --start/--stop; 48 files
             cover 2023-08-08 00:00 - 2023-08-09 23:00 UTC = all of
             2023-08-08 HST). Per-hour files are required because WindNinja
             runs EVERY band in a supplied forecast file - start/stop flags
             are ignored for wx-model runs - so a 1-band file is what makes
             an hourly run run one hour. wn_era5_hi_1hr.py picks the file
             for its hour automatically (HST + 10h -> UTC name).

HOW TO RUN
----------
One HST day (defaults to 2023-08-08, the Lahaina wind event):

  bash /home/wn1/wn_codes/era5_wrf_TEST/wn_era5_datetime_loop.sh 2023-08-08

That is the only command needed. It fans out all 24 HST hours in
parallel (xargs, 16 at a time), then builds the daily products and
daily validation.

WORKFLOW STEPS AND CODE
-----------------------
STEP 0 (done beforehand, on any machine with R + terra + RNetCDF):
  era5_to_wrfout.R  (kept in the project root, not this dir)
    Converts an ERA5 single-level GRIB (10u, 10v, 2t, tcc) to the
    WRF-mimic NetCDF. Takes --start/--stop UTC datetimes and --res.
    Key science: QCLOUD = tcc*100 (percent), T2 in Kelvin, data
    resampled to a uniform Mercator grid (MAP_PROJ=3) because WindNinja
    rejects lat/long grids and Mercator keeps ERA5's earth-relative
    winds rotation-free.

STEP 1: wn_era5_datetime_loop.sh <YYYY-MM-DD>
    Builds the 24 hourly datetime strings for the HST day and runs
    STEP 2 for each hour in parallel. Afterwards runs STEP 7.

STEP 2: wn_era5_workflow_all_co.sh <YYYY-MM-DD_HH:00:00>
    Per-hour dispatcher. Runs STEPS 3-6 in order for the hour.

STEP 3: multi_station_aggregator_era5.py <datetime>
    Fetches MADIS + HADS daily parse files from the public ikewai
    archive and extracts the target hour's observations per station.
    ERA5 edition also fetches the NEXT UTC day's files for HST hours
    >= 13 (the obs window crosses UTC midnight). Output:
    /home/wn1/data/input/stationFiles/<date>/<hr>/hrlyTable/
      <date>_<hr>_multi_station_all.csv

STEP 4: wn_era5_sta_file_mkr.R <datetime>
    Builds per-county validation station tables from step 3. ALL RAW
    stations kept: basic wind range limits only (speed 0-44.7 mps,
    direction 0-360); NO 2-SD outlier filter, NO temperature filter.
    Output: /home/wn1/data/input/wn_era5_input/<date>/<hr>/
      <date>_<hr>_<CO>_allData.csv        (CO = BI, MN, OA, KA)

STEP 5 (per county BI MN OA KA): wn_era5_hi_1hr.py <CO> <datetime>
    Runs WindNinja_cli with:
      initialization_method = wxModelInitialization
      forecast_filename     = the ERA5 WRF-mimic NetCDF
      diurnal_winds         = true
      time_zone             = Pacific/Honolulu (HST hour -> UTC band)
      mesh/ascii resolution = 250 m
    Output ASCII grids: /home/wn1/data/output/wn_era5_out/<date>/<hr>/<CO>/

STEP 6 (per county): wspd_dir_2_uv_mkr_era5.R <CO> <datetime>
    Resamples speed/direction to WGS84, derives u/v wind, masks to
    land, writes GeoTIFF products, and CROSS-VALIDATES the hour:
    model u/v sampled bilinearly at station locations (before masking),
    modeled speed/dir derived from the sampled vector, compared to the
    station observations. Writes per-hour validation pairs CSV and the
    hourly metadata txt containing all statistics.
    Output: /home/wn1/data/output/windProducts_era5/<date>/<hr>/<CO>/
      spd_dir_wind/  uv_wind/  stationData/  metadata/

STEP 7: statewide_hrly_wind_mkr_era5.R <datetime>
    Mosaics the 4 county maps into statewide hourly grids and writes
    statewide hourly metadata with POOLED (all-county) validation
    statistics for the hour.

STEP 8: dailyWind_mkr_era5.R <YYYY-MM-DD>   (run by STEP 1 at the end)
    Builds daily mean / median / max products per county + statewide
    mosaics, then DAILY CROSS-VALIDATION: for each station reporting
    >= 18 of 24 hours, the station's own daily statistic (computed from
    its hourly obs) is compared to the daily product grid sampled at
    the station location. Direction mirrors the grid construction:
    component-wise U/V statistic for mean/median; direction at the
    hour of maximum speed for max. All statistics written into every
    county and statewide daily metadata txt (attributes prefixed
    dailyMean*/dailyMedian*/dailyMax*), pairs CSVs saved under
    day/<CO>/validation/.

MIMIC-MAKING CODE (runs locally, not on the server)
---------------------------------------------------
  era5_to_wrfout.R      The ERA5 GRIB -> WRF-mimic NetCDF converter
                        (--input/--output/--start/--stop/--res). Needs R with
                        terra + RNetCDF.
  make_era5_hourly.R    Driver that splits a whole ERA5 GRIB into the
                        single-timestep hourly files WindNinja needs:
                        Rscript make_era5_hourly.R --input <grib> --outdir <dir> [--res m]
                        Files are named wrfout_era5_YYYYMMDD_HHMM.nc (UTC).
                        For small downloaded domains reduce --res so the
                        target grid is at least 4x4 (a 1.25 deg box needs
                        res <= ~13000).

ERA5/ERA5-LAND BLEND CHAIN (runs locally; needs R with terra + RNetCDF)
------------------------------------------------------------------------
  regrid_era5_to_land.R    ERA5 GRIB + ERA5-Land GRIB -> three NetCDFs on the
                           ERA5-Land 0.1 deg grid: (1) nearest-neighbor
                           resampled ERA5; (2) blend = ERA5-Land u10/v10/t2 on
                           its land cells, NN ERA5 elsewhere + tcc (ERA5-Land
                           has no tcc); (3) smoothed blend = bilinear
                           tent-kernel smoothing with the ERA5-Land land
                           pixels restored to original values. Each file
                           carries a land_source mask; all verified in-script.
  subset_blend_overlap.R   Subset a blend NetCDF to the datetime range present
                           in BOTH source GRIBs (no-op when they fully
                           overlap; the script says so).
  blend_to_wrfout.R        Blended NetCDF -> hourly single-timestep WRF-mimic
                           files (wrfout_era5_YYYYMMDD_HHMM.nc), NEAREST-
                           NEIGHBOR onto an 11 km Mercator grid (preserves the
                           exact land-cell values). Drop-in for the run
                           scripts via WN_FORECAST_DIR; --start/--stop UTC
                           window optional. Used by test_iniki_KA_19920911_blend.sh.

VISUALIZATION (runs locally; needs R with terra + magick)
---------------------------------------------------------
  plot_wind_gif.R    Per-hour PNG maps (speed fill + direction arrows, all
                     frames on a shared 0..globalmax m/s scale) from any dir
                     of WindNinja *_vel.asc/*_ang.asc outputs, then an
                     animated GIF (default 1 s per frame):
                     Rscript plot_wind_gif.R --indir <wn_out_dir> --outdir <dir> \
                       [--title t] [--name n] [--delay 100] [--fact 15] [--tz HST]
                     Direction is aggregated through U/V components (never
                     averaged as angles); arrow length scales with speed.

EVENT RUN SCRIPTS
-----------------
  test_single_hr_MN.sh       Maui, 2023-08-08 14:00 HST (Lahaina), one hour,
                             momentum solver ON.
  test_iniki_KA_19920911.sh  Kauai, 1992-09-11 11:00-17:00 HST (Hurricane
                             Iniki), 7 hourly runs, mass-only by default
                             (MOMENTUM=true env var enables the momentum
                             solver). Forecast files expected in
                             $WN_DATA_DIR/input/er5_wrf/hourly/1992_09/.

SHARED CODE
-----------
  wind_validation_fx.R
    Sourced by steps 6, 7 and 8. Defines the validation metrics
    (unit-tested):
      RMSE     root mean square error (speed mps; direction degrees on
               errors wrapped to +/-180)
      BIAS     mean error, model minus observed (direction: signed
               circular mean of wrapped errors)
      RME      combined U/V component bias = sqrt(biasU^2 + biasV^2), mps
      R-SQUARE speed: Pearson r^2; direction: squared
               Jammalamadaka-SenGupta circular correlation
    Calm observations (speed = 0) are excluded from direction metrics
    only. Also provides met<->UV conversions and name-based metadata
    lookup (metaGet).

OUTPUT TREE (created automatically on first run)
------------------------------------------------
/home/wn1/data/output/windProducts_era5/     [$WN_ERA5_PRODUCTS]
+-- <YYYYMMDD>/                              e.g. 20230808
    +-- 0000/ ... 2300/                      one dir per HST hour
    |   +-- BI|MN|OA|KA/                     one dir per county
    |   |   +-- spd_dir_wind/               *_vel_mps_wgs84.tif, *_dir_deg_wgs84.tif
    |   |   +-- uv_wind/                    *_uwind_wgs84.tif, *_vwind_wgs84.tif
    |   |   +-- stationData/                *_validation_sta.csv,
    |   |   |                               *_era5_validation_pairs.csv
    |   |   +-- metadata/                   <date>_<hr>_<CO>wind_metadata.txt
    |   |                                   (hourly validation stats)
    |   +-- statewide/
    |       +-- spd_dir_wind/  uv_wind/     hourly statewide mosaics
    |       +-- metadata/                   <date>_<hr>_statewide_wind_metadata.txt
    |                                       (pooled hourly validation stats)
    +-- day/
        +-- BI|MN|OA|KA/
        |   +-- spd_dir_wind/               daily mean/median/max vel + dir tifs
        |   +-- uv_wind/                    daily mean/median/max u/v tifs
        |   +-- max_hr/                     hour-of-daily-max raster (max only)
        |   +-- validation/                 <date>_<CO>_<func>_era5_validation_pairs.csv
        |   +-- metadata/                   <date>_<CO>_<func>_wind_metadata.txt
        |                                   (daily validation stats per statistic)
        +-- statewide/
            +-- spd_dir_wind/  uv_wind/  max_hr/   daily statewide mosaics
            +-- metadata/                  <date>_statewide_<func>_wind_metadata.txt
                                           (pooled daily validation stats)

Intermediate trees (also auto-created, inputs to the above):
/home/wn1/data/input/stationFiles/<date>/<hhmm>/   raw + hourly obs tables   [$WN_STATION_DIR]
/home/wn1/data/input/wn_era5_input/<date>/<hhmm>/  county validation tables  [$WN_ERA5_INPUT]
/home/wn1/data/output/wn_era5_out/<date>/<hhmm>/   raw WindNinja ASCII grids [$WN_ERA5_OUT]

The layout deliberately parallels the operational windProducts/ tree,
one directory over - nothing is ever written into windProducts/ itself.

WHERE THE VALIDATION NUMBERS END UP
-----------------------------------
  Hourly, per county:
    windProducts_era5/<date>/<hr>/<CO>/metadata/<date>_<hr>_<CO>wind_metadata.txt
  Hourly, statewide (pooled):
    windProducts_era5/<date>/<hr>/statewide/metadata/<date>_<hr>_statewide_wind_metadata.txt
  Daily, per county and statistic (mean/median/max):
    windProducts_era5/<date>/day/<CO>/metadata/<date>_<CO>_<func>_wind_metadata.txt
  Daily, statewide (pooled):
    windProducts_era5/<date>/day/statewide/metadata/<date>_statewide_<func>_wind_metadata.txt
  Raw per-station pairs (obs vs model):
    hourly: .../<hr>/<CO>/stationData/<date>_<hr>_<CO>_era5_validation_pairs.csv
    daily:  .../day/<CO>/validation/<date>_<CO>_<func>_era5_validation_pairs.csv

PATH OVERRIDES (env vars, all optional)
---------------------------------------
  WN_CODE_DIR        this code dir        (default /home/wn1/wn_codes/era5_wrf_TEST)
  WN_DATA_DIR        data root            (default /home/wn1/data/)
  WN_CLI             WindNinja_cli path   (default /home/wn1/wn_build/src/cli/WindNinja_cli)
  WN_FORECAST_DIR    hourly mimic files   (default $WN_DATA_DIR/input/er5_wrf/hourly/)
  WN_FORECAST_FILE   force one exact file (default <WN_FORECAST_DIR>/wrfout_era5_<UTC>.nc)
  WN_STATION_DIR     raw obs tree         (default /home/wn1/data/input/stationFiles/)
  WN_ERA5_INPUT      validation tables    (default /home/wn1/data/input/wn_era5_input)
  WN_ERA5_OUT        WindNinja output     (default /home/wn1/data/output/wn_era5_out)
  WN_ERA5_PRODUCTS   final products       (default /home/wn1/data/output/windProducts_era5)
  WN_MASK_DIR        county masks         (default /home/wn1/data/mask)

REQUIREMENTS ALREADY ON THE SERVER
----------------------------------
  WindNinja_cli build; python3 (pandas, numpy, requests, pytz); Rscript
  (raster, terra); DEMs /home/wn1/data/lndscp/{BI,MN,OA,KA}_utm_z4.tif;
  masks /home/wn1/data/mask/{bi,mn,oa,ka}_mask.tif; /home/wn1/data/runLog/.
  Output trees are created automatically.

PROVENANCE / MORE DETAIL
------------------------
  Adapted from HCDP/windninja_hrly_wind (pristine originals kept in the
  project's reference_orig/ dir). Full design rationale and the
  WindNinja source-code evidence trail: era5_to_windninja_wrf_mimic.md;
  session work summary: WORK_SUMMARY.md (both in the project root).
  Contact: Matthew Lucas (mplucas@hawaii.edu).
