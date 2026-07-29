# ERA5-workflow edition of wspd_dir_2_uv_mkr.R (HCDP/windninja_hrly_wind).
# Makes UV / resampled speed+dir products from the ERA5 wx-model WindNinja output
# AND cross-validates the hourly grids against independent station observations.
# Differences from the original:
#   - indir/outdir/wnInput point to the ERA5 trees (wn_era5_out, windProducts_era5,
#     wn_era5_input); overridable with env vars.
#   - Station data is validation truth (raw stations, basic range limits only) -
#     stations were NEVER model inputs, so this is a true independent validation.
#   - Model u/v are sampled (bilinear) at station locations; modeled speed/direction
#     derived from the sampled vector. Validation pairs written to stationData dir.
#   - Hourly metadata gains RMSE / RME (combined U/V bias) / BIAS / R-SQUARE for
#     speed and (circular) direction, and the dataStatement describes the ERA5 run.

#custom functions
domainHI<-function(co){
  domain<-if(co=="BI"){"Hawaii"}else
    if(co=="MN"){"Maui"}else
      if(co=="OA"){"Oahu"}else
        if(co=="KA"){"Kauai"}else
          if(co=="HI"){"Hawaii"}else{"NA"}
  return(domain)
}

dataconvert<-function(x){
  return(ifelse(is.numeric(x),as.character(round(x,5)),as.character(x)))
}

metaMakerCo<-function(co, DTin, dataDir, staObs, valStats, staFile){
  require(raster)

  #DT for dir names
  DTinre<-as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone())
  DTday<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%Y%m%d')
  DThr<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%H00')
  fileDT<-format(DTinre,"%Y%m%d_%H%M")

  #set final wn data dir
  dataDir<-paste0(dataDir,"/",DTday,"/",DThr,"/",co)

  #find files
  setwd(dataDir)
  tifPaths<-list.files(pattern = "tif",recursive = T) #get tifs
  tifFilesUL<-unlist(strsplit(tifPaths, "/", fixed = FALSE))
  tifFiles<-tifFilesUL[grep("tif",tifFilesUL)] #subset only file names

  #station obs min max (validation truth, not model input)
  minSta<-ifelse(nrow(staObs)>0,round(min(staObs$Speed),2),NA)
  maxSta<-ifelse(nrow(staObs)>0,round(max(staObs$Speed),2),NA)

  coHrVel<-raster(tifPaths[grep("vel",tifPaths)])
  minMap<-round(cellStats(coHrVel,min),2)
  maxMap<-round(cellStats(coHrVel,max),2)

  #make meta values
  statement<-paste("These", format(DTinre, "%B %d %Y %H:%M %P"),
                   "HST hourly wind maps of",
                   domainHI(co),ifelse(co=="OA","(Honolulu City and County)",ifelse(co=="HI","statewide","county")),
                   "are a high spatial resolution (~250m) gridded wind modeled prediction of 10 meter wind speed, direction, and the conversion of these variables to u-wind and v-wind components. These maps were produced by a weather-model initialization (wxModelInitialization) run of the windninja (https://ninjastorm.firelab.org/windninja/) land surface wind model, initialized with ERA5 reanalysis surface fields (10 m u/v wind, 2 m temperature, total cloud cover) converted to a WRF-format forecast file, with the diurnal slope-flow model enabled, to account for land cover and topographic influence on wind speed and direction. NO station data was used as model input. All available raw wind station observations within the county domain ( n=",
                   nrow(staObs),
                   ", basic range limits only, no outlier filtering) were used exclusively as independent cross-validation truth. Minimum and maximum wind speed observed at validation stations were:",
                   minSta, "mps and", maxSta,
                   "mps respectively with modeled map minimum and maximum wind speeds being",
                   minMap, "mps and", maxMap,
                   "mps respectively.",
                   valStatement(valStats),
                   "Validation metric definitions: RMSE = root mean square error; BIAS = mean error (model minus observed); RME = combined U/V wind component bias sqrt(biasU^2+biasV^2); R-SQUARE = squared Pearson correlation for speed and squared Jammalamadaka-SenGupta circular correlation for direction; direction errors wrapped to +/-180 degrees and calm observations excluded from direction metrics. These wind map data are an approximation, use data with caution.")
  #Make meta df wide
  metaWide<-data.frame(
    dataStatement=statement,
    keywords=paste0(domainHI(co), ", Hawaiian Islands, wind, hourly wind, wind speed, wind direction, u-wind, v-wind, ERA5, cross-validation"),
    county=ifelse(co=="HI","KA, OA, MN, BI",co),
    dataStartDateTime=DTinre,
    dataEndDateTime=DTinre+((59*60)+59),
    dataDate=format(DTinre,"%Y-%m-%d"),
    dataHour=format(DTinre,"%H:%M HST"),
    dateProduced=Sys.Date(),
    dataVersionType="experimental",
    windInputFile="ERA5 WRF-mimic forecast (see modelInitialization)",
    windHeight="10",
    windHeightUnit="meters",
    windDirUnit="degree",
    windSpeedUnit="mps",
    windUVunit="mps",
    stationCount=nrow(staObs),
    staWindMinMPS=minSta,
    staWindMaxMPS=maxSta,
    gridWindMinMPS=minMap,
    gridWindMaxMPS=maxMap,
    windGridFiles=paste(tifFiles,collapse = ", "),
    GeoCoordUnits="Decimal Degrees",
    GeoCoordRefSystem="EPSG:4326",
    Xresolution=res(coHrVel)[1],
    Yresolution=res(coHrVel)[2],
    ExtentXmin=extent(coHrVel)[1],
    ExtentXmax=extent(coHrVel)[2],
    ExtentYmin=extent(coHrVel)[3],
    ExtentYmax=extent(coHrVel)[4],
    modelInitialization="wxModelInitialization: ERA5 reanalysis via WRF-mimic NetCDF (era5_to_wrfout.R), diurnal_winds=true",
    validationStationFile=staFile,
    validationStationSet="all raw stations, basic range limits only, no outlier filter; independent of model input",
    credits="All data produced by University of Hawaii at Manoa Water Resource Research Center (WRRC). Support provided by Change Hawaii EPSCoR funded by the National Science Foundation under EPSCoR Research Infrastructure Improvement Award #OIA-2149133",
    contacts="Matthew Lucas (mplucas@hawaii.edu), Keri Kodama (kodamak8@hawaii.edu), Sayed Bateni (smbateni@hawaii.edu), Ryan Longman (rlongman@hawaii.edu), Thomas Giambelluca (thomas@hawaii.edu)"
  )

  final_meta<-rbind(data.frame(attribute=as.character(names(metaWide)),value=as.character(lapply(metaWide[1,], dataconvert))))
  #append the cross-validation metric rows
  final_meta<-rbind(final_meta,valAttrRows(valStats))
  return(final_meta)
}

wsDir2UV <- function(dir, ws) {
  require(raster)

  dir[dir == 360] <- 0 #change 360 degrees North to 0
  dir_rad <- dir * (pi / 180) # Convert direction to radians for trig functions

  # Calculate U (East-West) component
  # Positive U is wind blowing TO the East
  u_wind <- ws * -sin(dir_rad)
  names(u_wind) <- "uwind"
  # Calculate V (North-South) component
  # Positive V is wind blowing TO the North
  v_wind <- ws * -cos(dir_rad)
  names(v_wind) <- "vwind"

  # Combine the two layers into a single SpatRaster and return it
  uvwind<-stack(u_wind,v_wind)
  return(uvwind)
}

wnSaveReUV<-function(DTin,mask,co,indir, wnInput, outdir){
  require(raster)

  maskCo<-raster(mask)
  message("mask loaded")

  #DT for dir names
  DTinre<-as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone())
  DTday<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%Y%m%d')
  DThr<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%H00')
  fileDT<-format(DTinre,"%Y%m%d_%H%M")
  message("date loaded")

  #define day/hr/co indir
  indir<-paste0(indir,"/",DTday,"/",DThr,"/",co)

  #get wn output data
  setwd(indir)
  velWN<-list.files(pattern=paste0("*_250m_vel.asc"))
  message(paste("wn vel found:",velWN))

  angWN<-list.files(pattern=paste0("*_250m_ang.asc"))
  message(paste("wn ang found:",angWN))

  #make rasters
  print(velWN)
  print(angWN)
  ws<-raster(velWN) #speed
  dir<-raster(angWN) #direction
  message("wn rasters loaded")

  #resample speed and dir
  dirRe<-projectRaster(dir, maskCo, method = 'ngb')
  wsRe<-projectRaster(ws, maskCo)
  wsDirRe<-stack(wsRe,dirRe)
  names(wsDirRe)<-c("spdMPS","dirDeg")
  message("dir and ws resampled")

  #make uv wind
  uvwind<-wsDir2UV(ws=wsDirRe[["spdMPS"]],dir=wsDirRe[["dirDeg"]])
  message("uvwind calc done")

  ##cross-validation: sample model u/v at station locations BEFORE masking
  ##(stations near the coast can sit in masked cells)
  staFileName<-paste0(fileDT,"_",co,"_allData.csv")
  staPathFull<-paste0(wnInput,"/",DTday,"/",DThr,"/",staFileName)
  if(file.exists(staPathFull)){
    staObs<-read.csv(staPathFull,stringsAsFactors=FALSE)
  }else{
    message(paste("WARNING: no validation station file at",staPathFull))
    staObs<-data.frame(Station_Name=character(0),Lat=numeric(0),Lon=numeric(0),
                       Speed=numeric(0),Direction=numeric(0))
  }
  if(nrow(staObs)>0){
    pts<-cbind(staObs$Lon,staObs$Lat)
    modU<-raster::extract(uvwind[["uwind"]],pts,method="bilinear")
    modV<-raster::extract(uvwind[["vwind"]],pts,method="bilinear")
    modWS<-uv2ws(modU,modV)
    modDir<-uv2dir(modU,modV)
    valPairs<-data.frame(Station_Name=staObs$Station_Name,
                         Lat=staObs$Lat,Lon=staObs$Lon,
                         obsSpeed=staObs$Speed,obsDir=staObs$Direction,
                         modSpeed=modWS,modDir=modDir,
                         modU=modU,modV=modV,
                         obsU=met2u(staObs$Speed,staObs$Direction),
                         obsV=met2v(staObs$Speed,staObs$Direction))
  }else{
    valPairs<-data.frame(Station_Name=character(0),Lat=numeric(0),Lon=numeric(0),
                         obsSpeed=numeric(0),obsDir=numeric(0),
                         modSpeed=numeric(0),modDir=numeric(0),
                         modU=numeric(0),modV=numeric(0),
                         obsU=numeric(0),obsV=numeric(0))
  }
  valStats<-windValStats(obsWS=valPairs$obsSpeed,obsDir=valPairs$obsDir,
                         modWS=valPairs$modSpeed,modDir=valPairs$modDir)
  message(paste("hourly cross-validation:",co,DTin,"n =",valStats$n))

  #mask only land surface
  wsDirRe<-mask(wsDirRe,maskCo)
  uvwind<-mask(uvwind,maskCo)
  message("new rasters data masked")

  ##make dirs & write files
  setwd(outdir)
  dir.create(DTday,showWarnings = F)
  dir.create(paste0(DTday,"/",DThr),showWarnings = F)
  dirHRco<-paste0(DTday,"/",DThr,"/",co)
  dir.create(dirHRco,showWarnings = F)
  dir.create(paste0(dirHRco,"/spd_dir_wind"),showWarnings = F)
  dir.create(paste0(dirHRco,"/uv_wind"),showWarnings = F)
  dir.create(paste0(dirHRco,"/metadata"),showWarnings = F)
  dir.create(paste0(dirHRco,"/stationData"),showWarnings = F)

  #speed & dir
  setwd(paste0(outdir,"/",dirHRco,"/spd_dir_wind"))
  writeRaster(wsDirRe[["dirDeg"]],paste0(fileDT,"_",co,"_dir_deg_wgs84.tif"),overwrite=TRUE)
  writeRaster(wsDirRe[["spdMPS"]],paste0(fileDT,"_",co,"_vel_mps_wgs84.tif"),overwrite=TRUE)
  message(paste(co, DTin, "ws & dir raster written"))

  #u and v wind
  setwd(paste0(outdir,"/",dirHRco,"/uv_wind"))
  writeRaster(uvwind[["uwind"]],paste0(fileDT,"_",co,"_uwind_wgs84.tif"),overwrite=TRUE)
  writeRaster(uvwind[["vwind"]],paste0(fileDT,"_",co,"_vwind_wgs84.tif"),overwrite=TRUE)
  message(paste(co, DTin, "u v raster written"))

  #save validation station copy + validation pairs to export
  setwd(paste0(outdir,"/",dirHRco,"/stationData"))
  if(nrow(staObs)>0){
    wnStaName<-paste0(fileDT,"_",co,"_validation_sta.csv")
    write.csv(staObs,wnStaName,row.names = F)
    message(paste(wnStaName,"validation station file written"))
  }
  valPairsName<-paste0(fileDT,"_",co,"_era5_validation_pairs.csv")
  write.csv(valPairs,valPairsName,row.names = F)
  message(paste(valPairsName,"validation pairs file written"))

  #make metadata for county hr run
  metaDF<-metaMakerCo(co=co,DTin=DTin, dataDir=outdir, staObs=staObs,
                      valStats=valStats, staFile=staFileName)
  setwd(paste0(outdir,"/",dirHRco,"/metadata"))
  metaName<-paste0(fileDT,"_",co,"wind_metadata.txt")
  message(getwd())
  write.table(metaDF,file=metaName,quote=F,sep = "\t",row.names = FALSE)
  message(paste(metaName, "metadata file written"))

  #final
  message(paste(co, DTin, "all final wind files written"))
  return(paste(fileDT,co,sep="_"))
}

args <- commandArgs(trailingOnly = TRUE)
co=args[1]
DTin<-args[2]

#dirs: ERA5 trees (env-overridable), never the operational point-workflow trees
codeDir<-Sys.getenv("WN_CODE_DIR","/home/wn1/wn_codes/era5_wrf_TEST")
source(file.path(codeDir,"wind_validation_fx.R"))

DTday<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%Y%m%d')
DThr<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%H00')
indir<-Sys.getenv("WN_ERA5_OUT","/home/wn1/data/output/wn_era5_out")
wnInput<-Sys.getenv("WN_ERA5_INPUT","/home/wn1/data/input/wn_era5_input")
outdir<-Sys.getenv("WN_ERA5_PRODUCTS","/home/wn1/data/output/windProducts_era5")
dir.create(outdir,showWarnings=F,recursive=T)
mask<-paste0(Sys.getenv("WN_MASK_DIR","/home/wn1/data/mask"),"/",tolower(co),"_mask.tif") #co mask file in wgs84 dd.ddd

print(paste("making ERA5 WN products + validation for:",co,DTin))

wnSaveReUV(DTin=DTin,mask=mask,co=co,indir=indir,wnInput=wnInput,outdir=outdir)

message(paste("ERA5 wind UV reproject + validation code ran:",co,DTin))

#pau
