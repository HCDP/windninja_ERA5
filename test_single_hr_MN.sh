#!/bin/bash
# Single-hour WindNinja test run: Maui (MN), 2023-08-08 14:00 HST (= 2023-08-09 00:00 UTC),
# initialized from the ERA5 WRF-mimic forecast file. Output goes to /home/wn1/data/wn_test_outs.
# Same flags as the full workflow (wn_era5_hi_1hr.py) except output_path, single hour,
# and num_threads raised to 4 since nothing else is running in parallel.

OUT_DIR=/home/wn1/data/wn_test_outs
mkdir -p "$OUT_DIR"

/home/wn1/wn_build/src/cli/WindNinja_cli \
  --num_threads=4 \
  --elevation_file=/home/wn1/data/lndscp/MN_utm_z4.tif \
  --initialization_method=wxModelInitialization \
  --forecast_filename=/home/wn1/data/input/er5_wrf/wrfout_era5_20230808_20230809.nc \
  --time_zone=Pacific/Honolulu \
  --diurnal_winds=true \
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
  --momentum_flag=false \
  --start_year=2023 --start_month=08 --start_day=08 --start_hour=14 --start_minute=00 \
  --stop_year=2023  --stop_month=08  --stop_day=08  --stop_hour=14  --stop_minute=00

echo "WindNinja exit code: $?"
ls -la "$OUT_DIR"
