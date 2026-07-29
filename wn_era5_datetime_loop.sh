#run all 24 HST hours of a day in parallel, ERA5 wx-model edition
# ERA5 edition of wn_datetime_loop.sh (HCDP/windninja_hrly_wind).
# Input: HST day as YYYY-MM-DD (defaults to 2023-08-08, the Lahaina wind event).
# Each hour HH becomes a WindNinja run at YYYY-MM-DD_HH:00:00 Pacific/Honolulu,
# initialized from the ERA5 WRF-mimic file (which stores UTC; 2023-08-08 HST
# maps to 2023-08-08 10:00 .. 2023-08-09 09:00 UTC, inside the mimic's range).
# Requires on the run host: WindNinja_cli, the {BI,MN,OA,KA}_utm_z4.tif DEMs,
# and the mimic NetCDF at $WN_DATA_DIR/input/ (see wn_era5_hi_1hr.py).
CODE_DIR="${WN_CODE_DIR:-/home/wn1/wn_codes/era5_wrf_TEST}"
DATA_DIR="${WN_DATA_DIR:-/home/wn1/data}"

start_date=${1:-2023-08-08}
echo "$start_date"
dts=() #empty to store

#loop to make date with hours
for i in $(seq 0 23); do
  # Make hr in seq
  hour=$(printf "%02d" $i)
  echo "$hour"

  # Construct the new date-time string
  new_dt="${start_date}_${hour}:00:00"
  dts+=("$new_dt")
done

output_string=$(printf "%s\n" "${dts[@]}")
echo $output_string

MAX_PARALLEL=16 # Set the maximum number of parallel processes

echo "$output_string" | xargs -n 1 -P "$MAX_PARALLEL" bash "$CODE_DIR"/wn_era5_workflow_all_co.sh

echo "Parallel processing finished with control (using xargs)."

#run max/mean/median day r-code with daily cross-validation
Rscript "$CODE_DIR"/dailyWind_mkr_era5.R "$start_date" > "$DATA_DIR"/runLog/${start_date}_era5_day_stats_out.txt 2> "$DATA_DIR"/runLog/${start_date}_era5_day_stats_err.txt
wait
echo "All $start_date day wind stats + cross-validation finished!"
