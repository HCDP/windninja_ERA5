import sys
import subprocess
from datetime import datetime
from datetime import timedelta
import csv
import os

# ERA5 edition of wn_point_hi_1hr.py (HCDP/windninja_hrly_wind).
# Differences from the point-initialization original:
#   - initialization_method = wxModelInitialization, driven by the ERA5->WRF-mimic
#     NetCDF produced by era5_to_wrfout.R (must cover the requested hour, in UTC).
#   - station-related flags removed (match_points, fetch_type, wx_station_filename,
#     write_wx_station_kml).
#   - diurnal_winds ENABLED: the mimic file carries real gridded T2 and cloud cover,
#     which is the reason to use the WRF-mimic route at all.
#   - output tree wn_out/ -> wn_era5_out/ so station-based products are not touched.
# Paths default to the wn1 server layout; override with env vars WN_DATA_DIR,
# WN_CLI and WN_FORECAST_FILE when running elsewhere.

def task(extent,st):

    # Setting the data directory
    #full_path = os.path.realpath(__file__)
    directory_path = os.environ.get("WN_DATA_DIR", "/home/wn1/data/")

    # the input files
    elevation_file = extent + "_utm_z4.tif"

    # the ERA5 WRF-mimic forecast file (UTC times inside; WindNinja matches the
    # HST start/stop below to the right band via --time_zone)
    forecast_file = os.environ.get(
        "WN_FORECAST_FILE",
        directory_path + "input/wrfout_era5_20230808_20230809.nc")

    #get time flags
    st_dt = datetime.strptime(st,'%Y-%m-%d_%H:%M:%S')
    st_year = st_dt.strftime('%Y')
    st_mon = st_dt.strftime('%m')
    st_day = st_dt.strftime('%d')
    st_hr = st_dt.strftime('%H')
    st_min = st_dt.strftime('%M')

    en_dt = st_dt
    en_year = st_dt.strftime('%Y')
    en_mon = st_dt.strftime('%m')
    en_day = st_dt.strftime('%d')
    en_hr = st_dt.strftime('%H')
    en_min = st_dt.strftime('%M')

    #make DT and extent outdirs
    outDate = st_year + st_mon + st_day
    outTime = st_hr + st_min
    out_dir_day = directory_path + "output/wn_era5_out/" + outDate
    out_dir_hr = directory_path + "output/wn_era5_out/" + outDate + "/" + outTime
    os.makedirs(out_dir_day, exist_ok=True)
    os.makedirs(out_dir_hr, exist_ok=True)
    out_dir = out_dir_hr + "/" + extent
    os.makedirs(out_dir, exist_ok=True)

    ## command to run by WindNinja_cli ##
    cmd = os.environ.get("WN_CLI", "/home/wn1/wn_build/src/cli/WindNinja_cli")

    cmd += " --num_threads=1"
    cmd += " --elevation_file=" + directory_path + "lndscp/" + elevation_file ##for "dem/" vegetation flag must have value; for "lndscp/" all veg and veg hieght is contained in the file
    cmd += " --initialization_method=wxModelInitialization"
    cmd += " --forecast_filename=" + forecast_file
    cmd += " --time_zone=Pacific/Honolulu"
    cmd += " --diurnal_winds=true" ##diurnal wind: ON - ERA5 mimic supplies gridded T2 + cloud cover
    cmd += " --output_wind_height=10.0"
    cmd += " --units_output_wind_height=m"
    ##cmd += " --vegetation=grass" ##ON when elevation_file is dem; OFF when elevation_file is lndscp
    cmd += " --mesh_resolution=250"
    cmd += " --units_mesh_resolution=m"
    cmd += " --ascii_out_resolution=250"
    cmd += " --units_ascii_out_resolution=m"
    cmd += " --write_goog_output=false"
    cmd += " --write_shapefile_output=false"
    cmd += " --write_ascii_output=true"
    cmd += " --write_farsite_atm=false"
    cmd += " --output_path="+out_dir
    cmd += " --momentum_flag=false"
    cmd += " --start_year="+st_year
    cmd += " --start_month="+st_mon
    cmd += " --start_day="+st_day
    cmd += " --start_hour="+st_hr
    cmd += " --start_minute="+st_min
    cmd += " --stop_year="+en_year
    cmd += " --stop_month="+en_mon
    cmd += " --stop_day="+en_day
    cmd += " --stop_hour="+en_hr
    cmd += " --stop_minute="+en_min

    # creating a process to run the command
    print(cmd)
    proc = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)

    print("After proc")
    # wait for the process to complete
    return_code = proc.wait()
    print("Return code ",return_code)

    # terminate the process
    proc.terminate()
    print("Process complete")


if __name__=="__main__":
    extent_name = sys.argv[1]
    st_date = sys.argv[2]

    task(extent_name,st_date)
