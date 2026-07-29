# CONVENTION: directory/file names use the HST target hour; all in-file `date_time` values are UTC (Zulu, "...Z").
#
# ERA5-workflow edition of wn_r_file_mkr.R (HCDP/windninja_hrly_wind).
# The station tables produced here are CROSS-VALIDATION TRUTH, not WindNinja input
# (the ERA5 runs are initialized from the WRF-mimic forecast file). Differences:
#   - NO 2-SD U/V outlier filter: all raw stations kept (user choice - validation
#     uses every reporting station after basic range limits).
#   - Range limits applied to wind only (speed 0-44.7 mps, direction 0-360);
#     the temperature limit is NOT applied because temperature is not validated
#     and a bad thermometer should not delete a good anemometer.
#   - Writes only the per-county {datetime}_{co}_allData.csv tables; the per-station
#     CSVs and {co}_stations.csv WindNinja lists are not needed and are skipped.
#   - Output tree defaults to /home/wn1/data/input/wn_era5_input (NOT wn_input) so
#     the operational point-workflow inputs are never overwritten.

#functions to make ERA5 validation station files

wndfTZtrunc<-function(DF){ #formerly wndfTZ2hst; no TZ shift anymore
  DF$time<-strptime(DF$time, format="%Y-%m-%d %H:%M", tz="UTC")
  DF$time<-format(DF$time,"%Y-%m-%d %H:00:00")
  return(DF)
} # keeps time in UTC, truncates to top of hour

sknListMkr<-function(df,by="County"){
  unq<-unique(df[,by])
  sknByList<-list()
  for(i in unq){
    skn<-unique(df[df[,by]==i,"SKN"])
    sknByList[[i]]<-as.character(skn)
  }
  return(sknByList)
} #makes a list with a df per SKN

makeWNfiles<-function(DTin,
                      inDir=Sys.getenv("WN_STATION_DIR","/home/wn1/data/input/stationFiles"),
                      outDir=Sys.getenv("WN_ERA5_INPUT","/home/wn1/data/input/wn_era5_input")){
  Sys.setenv(TZ="HST") #DTin and all dir/file names below are HST; in-file data timestamps stay UTC

  #make date time dirs
  dataDTre<-as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone())
  stafileDT<-format(dataDTre,"%Y%m%d_%H%M") #HST-derived
  dirDate<-format(dataDTre,"%Y%m%d") #HST-derived
  dirHr<-format(dataDTre,"%H%M") #HST-derived

  #make dt indir
  indir<-paste(inDir,dirDate,dirHr,"hrlyTable",sep="/")

  #make and create dt outdir
  dir.create(outDir,showWarnings = F,recursive = T)
  dir.create(paste(outDir,dirDate,sep="/"),showWarnings = F)
  dir.create(paste(outDir,dirDate,dirHr,sep="/"),showWarnings = F)
  outdir<-paste(outDir,dirDate,dirHr,sep="/")

  #check dirs
  message(paste("current wd:",getwd()))
  message(paste("indir:",indir))
  message(paste("outdir:",outdir))

  #get station table
  setwd(indir)
  wnAlldf<-read.csv(paste0(stafileDT,"_multi_station_all.csv"),stringsAsFactors = FALSE)

  ##qaqc: basic wind range limits ONLY (no SD outlier filter, no temperature filter)
  message("qaqc basic wind range limits (raw station set, no outlier filter)...")
  wnAlldf<-wnAlldf[wnAlldf$Speed>=0 & wnAlldf$Speed<=44.7,] #range limits for wind speed: 0mph-100mph
  wnAlldf<-wnAlldf[wnAlldf$Dir>=0 & wnAlldf$Dir<=360,] #valid met direction only

  ##trunc to top of hour (stays UTC, no TZ shift)
  wnAlldf<-wndfTZtrunc(DF=wnAlldf)

  #seperate by co domain
  sknByCounty<-sknListMkr(df=wnAlldf)

  finalWNinput<-data.frame(
    Station_Name=as.character(wnAlldf$SKN),
    Coord_Sys=tolower("GEOGCS"),
    Datum="WGS84",
    Lat=wnAlldf$LAT,
    Lon=wnAlldf$LON,
    Height=10,
    Height_Units="meters",
    Speed=wnAlldf$Speed,
    Speed_Units="mps",
    Direction=wnAlldf$Dir,
    Temperature=wnAlldf$Temperature,
    Temperature_Units="C",
    Cloud_Cover=0,
    Radius_of_Influence=-1,
    Radius_of_Influence_Units="km",
    date_time=paste0(wnAlldf$time,"Z"), #real UTC/Zulu values; folder/file names are HST
    stringsAsFactors = F)

  #reorder by time (may not be necessary)
  finalWNinput<-finalWNinput[order(as.character(finalWNinput$date_time)),]

  #write per-county validation tables
  allCoCount<-data.frame()
  for(j in c("KA","OA","MN","BI")){
    CoSKN<-sknByCounty[[j]]
    subWNinput<-finalWNinput[finalWNinput$Station_Name %in% CoSKN,]
    write.csv(subWNinput,paste0(outdir,"/",stafileDT,"_",j,"_allData.csv"),row.names=F)
    coRow<-data.frame(dt=dataDTre,nsta=nrow(subWNinput),county=j)
    allCoCount<-rbind(allCoCount,coRow)
  }#end county loop

  #track co station output
  return(allCoCount)
}

#Get all command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(paste("making ERA5 validation station files for:",args))
out<-makeWNfiles(DTin=args)
print(out)
print(paste(args,"ERA5 validation station files made!"))
message(paste(args,"wn_era5_sta_file_mkr.R complete!"))
#PAU
