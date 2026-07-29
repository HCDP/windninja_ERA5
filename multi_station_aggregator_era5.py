"""
ERA5-workflow edition of multi_station_aggregator_with_meso.py (HCDP/windninja_hrly_wind).
Changes from the original:
    - For HST target hours >= 13, ALSO fetches the next UTC day's MADIS/HADS parse
      files: the observation window (HST hour + 10h .. +11h) crosses UTC midnight for
      hours 13-23 HST, and the daily parse files are UTC-dated. Both days' rows are
      concatenated before hour extraction so late-HST hours are not empty.
    - Output directory overridable with WN_STATION_DIR (default unchanged) so
      historical re-runs do not have to touch the operational tree.
    Everything else (sources, QA, output table format) is unchanged: output feeds
    the ERA5 cross-validation, not WindNinja input.
"""
#CMD LINE ARG: TARGET HOUR [SET IN HST, SCRIPT CONVERTS TO UTC]
#OUTPUT TABLE IN UTC TIMESTAMPS BUT COVERS SELECTED HST TARGET HOUR
import sys
import os
import pytz
import subprocess
import requests
import pandas as pd
import numpy as np
from functools import reduce
from datetime import datetime,timedelta
from io import StringIO

hcdp_api_token = os.environ.get('API_TOKEN')
aqs_base = 'https://ikeauth.its.hawaii.edu/files/v2/download/public/system/ikewai-annotated-data/HCDP/workflow_data/preliminary/data_aqs/data_outputs/'
master_link = 'https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/main/Hawaii_Master_Station_Meta.csv'
convert_link = "https://raw.githubusercontent.com/ikewai/hawaii_wx_station_mgmt_container/refs/heads/main/hi_mesonet_sta_status.csv"
master_df = pd.read_csv(master_link)

madis_path = aqs_base + 'madis/parse/'
hads_path = aqs_base + 'hads/parse/'

#make date-hr dir name
date_string = sys.argv[1]
original_format = "%Y-%m-%d_%H:%M:%S"
date_format = "%Y%m%d"
hr_format = "%H%M"
dt_object = datetime.strptime(date_string, original_format)
outTime = dt_object.strftime(hr_format)
outDate = dt_object.strftime(date_format)

#out dirs
outDir = os.environ.get('WN_STATION_DIR', '/home/wn1/data/input/stationFiles/')
outDirDay = os.path.join(outDir, outDate)
outDirHr = os.path.join(outDirDay, outTime)
os.makedirs(outDirDay, exist_ok=True)
os.makedirs(outDirHr, exist_ok=True)
rawOut = outDirHr + '/rawSource'
hrlyOut = outDirHr + '/hrlyTable'
os.makedirs(rawOut, exist_ok=True)
os.makedirs(hrlyOut, exist_ok=True)

STATIC_HEADER = ['SKN','Time','LAT','LON','Elev','County']
VAR_NAMES = ['Speed','Dir','Temperature']
VAR_TABLE = {'madis':['windSpeed','windDir','temperature'],'hads':['US','UD','TA'],'hi_mesonet':['WS_1_Avg','WDuv_1_Avg','Tair_1_Avg']}
SRC_LABELS = {'madis':{'stn':'stationId','time':'time','varname':'varname','value':'value','merge_id':'NWS.id'},
             'hads':{'stn':'staID','time':'obs_time','varname':'var','value':'value','merge_id':'NESDIS.id'},
             'hi_mesonet':{'stn':'SKN','time':'timestamp','varname':'variable','value':'value','merge_id':'SKN'}}
CONVERT_FACTORS = {'Speed':[(1/2.23693629),0],'Dir':[1,0],'Temperature':[(5/9),-32]}

def check_file(fname):
    if os.path.exists(fname):
        if os.stat(fname).st_size > 0:
            return True
        else:
            return False
    else:
        return False

def utc_dates_needed(targ_hour):
    """UTC calendar dates spanned by the observation window [HST+10h, HST+11h]."""
    st_utc = targ_hour + pd.Timedelta(hours=10)
    en_utc = st_utc + pd.Timedelta(hours=1)
    dates = sorted({st_utc.strftime(date_format), en_utc.strftime(date_format)})
    return dates

def aqs_extractor(targ_hour,var_table,src_labels):
    #extraction wrapper, calls source specific extractor
    src_list = ['madis','hads']
    nvars = len(var_table[src_list[0]])
    hr_code = targ_hour.strftime('%H')
    fetch_dates = utc_dates_needed(targ_hour)
    concat_src = pd.DataFrame()
    for src in src_list:
        print(src)
        #concatenate all fetched UTC-day files for this source
        day_frames = []
        for dc in fetch_dates:
            local_name = rawOut + '/' + '_'.join((dc,hr_code,src,'parsed'))+'.csv'
            if check_file(local_name):
                day_frames.append(pd.read_csv(local_name))
            else:
                print(local_name,"does not exist")
        if len(day_frames) == 0:
            continue
        src_df = pd.concat(day_frames,axis=0,ignore_index=True)
        stn_key = src_labels[src]['stn']
        time_key = src_labels[src]['time']
        var_key = src_labels[src]['varname']
        val_key = src_labels[src]['value']
        merge_key = src_labels[src]['merge_id']
        merged_var_all = []
        for i,varname in enumerate(VAR_NAMES):
            varname = var_table[src][i]
            var_df = src_df[src_df[var_key].isin([varname])]
            hour_hst = targ_hour
            hour_utc = targ_hour + pd.Timedelta(hours=10)
            st_time = hour_utc
            en_time = st_time + pd.Timedelta(hours=1)
            obs_time = pd.to_datetime(var_df[time_key])
            var_df.loc[:,time_key] = pd.to_datetime(var_df[time_key])
            var_hour = var_df.loc[obs_time[(obs_time>=st_time)&(obs_time<=en_time)].index].sort_values(by=[stn_key,time_key])
            uni_stns = var_hour[stn_key].unique()
            hour_midpoint = targ_hour + pd.Timedelta(hours=0.5)
            all_sel_stns = pd.DataFrame()
            for stn in uni_stns:
                stn_hour = var_hour[var_hour[stn_key]==stn]
                hour_midpoint_utc = targ_hour + pd.Timedelta(hours=0.5) + pd.Timedelta(hours=10)
                hour_time = obs_time.loc[stn_hour.index]
                diff = abs(hour_time-hour_midpoint_utc)
                sel_ind = diff[diff==diff.min()].index
                sel_val = stn_hour.loc[sel_ind]
                all_sel_stns = pd.concat([all_sel_stns,sel_val],axis=0)
            all_sel_stns = all_sel_stns.rename(columns={stn_key:merge_key})
            if all_sel_stns.empty:
                print(f"No data available for {varname} from {src} at {hr_code}:00.")
                continue
            else:
                merged_sel = all_sel_stns.merge(master_df,how='inner',on=merge_key)[['SKN',time_key,'LAT','LON','ELEV.m.','Island',val_key]]
                merged_sel = merged_sel.rename(columns={'Island':'County',time_key:'time',val_key:VAR_NAMES[i]})
                #if Hads convert to metric
                if src == 'hads':
                    a = CONVERT_FACTORS[VAR_NAMES[i]][0]
                    b = CONVERT_FACTORS[VAR_NAMES[i]][1]
                    merged_sel[VAR_NAMES[i]] = merged_sel[VAR_NAMES[i]].map(lambda x:a*x+b)
                merged_sel = merged_sel.replace({'County':['MA','MO','KO','LA']},{'County':'MN'}).sort_values(by='SKN',ignore_index=True)
                merged_var_all.append(merged_sel)
        if len(merged_var_all) > 0:
            src_var_final = reduce(lambda left,right: pd.merge(left,right,how='inner',on=['SKN','time','LAT','LON','ELEV.m.','County']),merged_var_all)
            concat_src = pd.concat([concat_src,src_var_final],axis=0)
    concat_src = concat_src.drop_duplicates(subset=['SKN'],keep='last')
    concat_src = concat_src.sort_values(by='SKN',ignore_index=True)
    return concat_src

if __name__=="__main__":
    #SET TIME IN HST
    if len(sys.argv) == 1:
        hst = pytz.timezone('HST')
        targ_hour = pd.to_datetime(datetime.now(hst)).tz_localize(None)
    else:
        targ_hour = sys.argv[1]
        targ_hour = pd.to_datetime(targ_hour,format="%Y-%m-%d_%H:%M:%S")

    #fetch MADIS and HADS files from gateway - every UTC day the obs window touches
    hr_code = targ_hour.strftime('%H')
    fetch_dates = utc_dates_needed(targ_hour)
    data_flag = 0
    n_attempts = 0
    for src in ['madis','hads']:
        for dc in fetch_dates:
            src_name = aqs_base + src + '/parse/' + '_'.join((dc,src,'parsed'))+'.csv'
            local_name = rawOut + '/' + '_'.join((dc,hr_code,src,'parsed'))+'.csv'
            n_attempts += 1
            cmd = ["wget",src_name,"-O",local_name]
            prc = subprocess.run(cmd)
            print(prc.returncode)
            if prc.returncode > 0:
                print("Issues encountered in downloading",local_name)
                data_flag += 1
            else:
                print("Success. Downloaded",local_name)

    if data_flag >= n_attempts:
        print("No data acquired")
        quit()
    output_df = aqs_extractor(targ_hour,VAR_TABLE,SRC_LABELS)
    if output_df.empty:
        print("Warning: No data extracted from source files.")
    output_time = targ_hour.strftime('%Y%m%d_%H%M')
    output_file = hrlyOut + '/' +'_'.join((output_time,'multi_station_all'))+'.csv'
    output_df.to_csv(output_file,index=False)
