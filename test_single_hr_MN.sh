#!/bin/bash
# Single-hour WindNinja test run: Maui (MN), 2023-08-08 14:00 HST (= 2023-08-09 00:00 UTC),
# initialized from a SINGLE-TIMESTEP ERA5 WRF-mimic file. WindNinja runs every band in
# a supplied forecast file (start/stop flags are ignored for wx-model runs on this
# build), so the 1-band hourly file is what makes this a 1-hour run.
# Output goes to /home/wn1/data/wn_test_outs. Same flags as the full workflow
# (wn_era5_hi_1hr.py) except output_path, num_threads=4 (nothing else running), and
# momentum_flag=true: conservation of mass AND momentum (NinjaFOAM solver). Requires
# a WindNinja build with OpenFOAM support and runs much longer than the mass-only
# solver (tens of minutes vs ~90 s). The full workflow stays mass-only.

OUT_DIR=/home/wn1/data/wn_test_outs
mkdir -p "$OUT_DIR"

/home/wn1/wn_build/src/cli/WindNinja_cli \
  --num_threads=4 \
  --elevation_file=/home/wn1/data/lndscp/MN_utm_z4.tif \
  --initialization_method=wxModelInitialization \
  --forecast_filename=/home/wn1/data/input/er5_wrf/hourly/wrfout_era5_20230809_0000.nc \
  --time_zone=Pacific/Honolulu \
  --diurnal_winds=true \
  --output_speed_units=mps \
  --output_wind_height=10.0 \
  --units_output_wind_height=m \
  --mesh_resolution=250 \
  --units_mesh_resolution=m \
  --ascii_out_resolution=250 \
  --units_ascii_out_resolution=m \
  --write_goog_output=false \
  --write_shapefile_output=false \
  --write_ascii_output=true \
  --write_farsite_atm=false \
  --output_path="$OUT_DIR" \
  --momentum_flag=true

echo "WindNinja exit code: $?"
ls -la "$OUT_DIR"
