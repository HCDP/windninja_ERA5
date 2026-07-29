# ERA5 edition of wn_workflow_all_co.sh (HCDP/windninja_hrly_wind).
# Station data IS fetched (multi_station_aggregator_era5.py + wn_era5_sta_file_mkr.R)
# but is used ONLY as independent cross-validation truth - the WindNinja runs are
# initialized from the ERA5 WRF-mimic forecast file (wxModelInitialization).
# Set WN_CODE_DIR / WN_DATA_DIR to relocate; defaults match the wn1 server layout.
CODE_DIR="${WN_CODE_DIR:-/home/wn1/wn_codes/era5_wrf_TEST}"
DATA_DIR="${WN_DATA_DIR:-/home/wn1/data}"

echo "DT input starting: $1"
#fetch + aggregate station obs for the hour (validation truth, not model input)
python3 "$CODE_DIR"/multi_station_aggregator_era5.py "$1" > "$DATA_DIR"/runLog/${1}_era5_ag_out.txt 2> "$DATA_DIR"/runLog/${1}_era5_ag_err.txt
wait
Rscript "$CODE_DIR"/wn_era5_sta_file_mkr.R "$1" > "$DATA_DIR"/runLog/${1}_era5_fm_out.txt 2> "$DATA_DIR"/runLog/${1}_era5_fm_err.txt
wait
for co in BI MN OA KA
do
	echo "$co $1 starting WN (ERA5 wx-model init)..."
	python3 "$CODE_DIR"/wn_era5_hi_1hr.py $co "$1" > "$DATA_DIR"/runLog/${1}_wn_era5_out.txt 2> "$DATA_DIR"/runLog/${1}_${co}_wn_era5_err.txt
	wait
	Rscript "$CODE_DIR"/wspd_dir_2_uv_mkr_era5.R $co "$1" > "$DATA_DIR"/runLog/${1}_era5_re_out.txt 2> "$DATA_DIR"/runLog/${1}_${co}_era5_re_err.txt
	wait
	echo "$co $1 finished!"
done
wait
Rscript "$CODE_DIR"/statewide_hrly_wind_mkr_era5.R "$1" > "$DATA_DIR"/runLog/${1}_era5_hihr_out.txt 2> "$DATA_DIR"/runLog/${1}_era5_hihr_err.txt
wait
echo "ALL $1 DT finished!"
