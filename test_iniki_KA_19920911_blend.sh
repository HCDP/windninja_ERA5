#!/bin/bash
# Hurricane Iniki WindNinja run - ERA5/ERA5-Land BLEND edition: Kauai (KA),
# 1992-09-11 11:00-17:00 HST (= 1992-09-11 21:00 UTC through 1992-09-12 03:00 UTC),
# initialized from the blended single-timestep mimic files made by
# blend_to_wrfout.R (ERA5-Land 0.1 deg values on land cells, smoothed ERA5
# elsewhere; nearest-neighbor onto an 11 km Mercator grid).
#
# Identical to test_iniki_KA_19920911.sh except the forecast directory and the
# output directory, so the two runs form a controlled comparison of pure-ERA5
# vs blended forcing over the same event/DEM.
#
# Forecast files expected in: /home/wn1/data/input/er5_wrf/hourly/1992_09_blend/
# Output:                     /home/wn1/data/wn_test_outs/iniki_1992_KA_blend/
#
# MOMENTUM=true enables the mass+momentum (NinjaFOAM) solver - much slower
# per hour and requires libWindNinja.so + applyInit built (see README).
# Default is the mass-only solver.

MOMENTUM=${MOMENTUM:-false}
FC_DIR=${WN_FORECAST_DIR:-/home/wn1/data/input/er5_wrf/hourly/1992_09_blend}
OUT_DIR=/home/wn1/data/wn_test_outs/iniki_1992_KA_blend
mkdir -p "$OUT_DIR"

# 11:00-17:00 HST on 1992-09-11 -> UTC valid times (HST + 10h)
UTC_TIMES="19920911_2100 19920911_2200 19920911_2300 19920912_0000 19920912_0100 19920912_0200 19920912_0300"

for utc in $UTC_TIMES; do
  fc="$FC_DIR/wrfout_era5_${utc}.nc"
  if [ ! -f "$fc" ]; then
    echo "MISSING forecast file, skipping: $fc"
    continue
  fi
  echo "=== running $utc UTC ($fc) ==="
  /home/wn1/wn_build/src/cli/WindNinja_cli \
    --num_threads=4 \
    --elevation_file=/home/wn1/data/lndscp/KA_utm_z4.tif \
    --initialization_method=wxModelInitialization \
    --forecast_filename="$fc" \
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
    --momentum_flag=$MOMENTUM
  echo "WindNinja exit code for $utc: $?"
done

ls -la "$OUT_DIR"
