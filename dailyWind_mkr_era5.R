###make day ws dirs, rasters and metadata per county and statewide - ERA5 wx-model edition
# ERA5 edition of dailyWind_mkr.R (HCDP/windninja_hrly_wind). Changes:
#   - filedir -> windProducts_era5; validation station tables from wn_era5_input
#     (env-overridable).
#   - DAILY CROSS-VALIDATION (station daily stat vs grid daily stat): for each
#     validation station the daily mean/median/max is computed from its own hourly
#     observations (stations reporting >= 18 of 24 hours) and compared to the
#     corresponding daily product grid sampled at the station location.
#     Direction daily stats mirror the grid construction: mean/median use the
#     component-wise (U,V) statistic converted back to direction; max uses the
#     direction at the hour of maximum speed.
#   - RMSE / RME (combined U/V bias) / BIAS / R-SQUARE for speed and (circular)
#     direction written into every county and statewide daily metadata file.
#   - county metadata values read BY ATTRIBUTE NAME in statewide meta (original
#     used row numbers, which breaks when attributes are added).

#define main dirs
filedir<-Sys.getenv("WN_ERA5_PRODUCTS","/home/wn1/data/output/windProducts_era5")
wnInput<-Sys.getenv("WN_ERA5_INPUT","/home/wn1/data/input/wn_era5_input")
codeDir<-Sys.getenv("WN_CODE_DIR","/home/wn1/wn_codes/era5_wrf_TEST")
source(file.path(codeDir,"wind_validation_fx.R"))

minHrsDaily<-18 #minimum reporting hours (of 24) for a station to enter daily validation

#get datetime in
args <- commandArgs(trailingOnly = TRUE)
DateIn=as.Date(args[1])

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

UV2wsDir <- function(u, v) {
  require(terra)

  # Calculate wind speed using Pythagorean theorem
  wind_speed <- sqrt(u^2 + v^2)
  names(wind_speed) <- "spdMPS"

  # Calculate mathematical direction in degrees using atan2
  math_dir_deg <- atan2(v, u) * (180 / pi)

  # Convert mathematical direction to meteorological direction and normalize
  wind_dir <- (270 - math_dir_deg) %% 360
  names(wind_dir) <- "dirDeg"

  # Combine the two layers into a single SpatRaster and return it
  wsdir<-c(wind_speed, wind_dir)
  return(wsdir)
}

##gather all 24 hourly validation station tables for a county/day
gatherObsDay<-function(DateIn,co,wnInput){
  DateFile<-format(as.Date(DateIn),'%Y%m%d')
  out<-data.frame()
  for(h in 0:23){
    hhmm<-sprintf("%02d00",h)
    f<-file.path(wnInput,DateFile,hhmm,paste0(DateFile,"_",hhmm,"_",co,"_allData.csv"))
    if(file.exists(f)){
      hrDF<-read.csv(f,stringsAsFactors=FALSE)
      if(nrow(hrDF)>0){
        hrDF$hour<-h
        out<-rbind(out,hrDF[,c("Station_Name","Lat","Lon","Speed","Direction","hour")])
      }
    }
  }
  return(out)
}

##per-station daily statistic from hourly obs, mirroring the grid construction
staDailyStat<-function(obsDay,func,minHrs=minHrsDaily){
  if(nrow(obsDay)==0){return(data.frame())}
  out<-data.frame()
  for(stn in unique(obsDay$Station_Name)){
    sd_<-obsDay[obsDay$Station_Name==stn,]
    sd_<-sd_[!duplicated(sd_$hour),]
    if(nrow(sd_)<minHrs){next}
    u<-met2u(sd_$Speed,sd_$Direction)
    v<-met2v(sd_$Speed,sd_$Direction)
    if(func=="max"){
      iMax<-which.max(sd_$Speed)
      wsStat<-sd_$Speed[iMax]
      dirStat<-sd_$Direction[iMax]
    }else{
      f<-get(func)
      wsStat<-f(sd_$Speed)
      dirStat<-uv2dir(f(u),f(v)) #component-wise stat -> direction, matches grid method
    }
    out<-rbind(out,data.frame(Station_Name=stn,Lat=sd_$Lat[1],Lon=sd_$Lon[1],
                              nHours=nrow(sd_),obsSpeed=wsStat,obsDir=dirStat))
  }
  return(out)
}

##daily cross-validation for one county/func; needs the day grids already written
valDailyCo<-function(DateIn,co,func,filedir,wnInput){
  require(terra)
  DateFile<-format(as.Date(DateIn),'%Y%m%d')
  dayDir<-paste0(filedir,"/",DateFile,"/day/",co)

  obsDay<-gatherObsDay(DateIn,co,wnInput)
  staDaily<-staDailyStat(obsDay,func)
  emptyPairs<-data.frame(Station_Name=character(0),Lat=numeric(0),Lon=numeric(0),
                         nHours=numeric(0),obsSpeed=numeric(0),obsDir=numeric(0),
                         modSpeed=numeric(0),modDir=numeric(0))
  if(nrow(staDaily)==0){
    return(list(stats=windValStats(numeric(0),numeric(0),numeric(0),numeric(0)),
                pairs=emptyPairs,nStations=0))
  }
  pts<-data.frame(lon=staDaily$Lon,lat=staDaily$Lat)

  velTif<-paste0(dayDir,"/spd_dir_wind/",paste(DateFile,co,func,"vel_mps_wgs84.tif",sep="_"))
  modWS<-terra::extract(rast(velTif),pts,method="bilinear")[,2]
  if(func=="max"){
    #direction-at-max grid is not smoothly interpolable - nearest cell
    dirTif<-paste0(dayDir,"/spd_dir_wind/",paste(DateFile,co,func,"dir_deg_wgs84.tif",sep="_"))
    modDir<-terra::extract(rast(dirTif),pts,method="simple")[,2]
  }else{
    uTif<-paste0(dayDir,"/uv_wind/",paste(DateFile,co,func,"uwind_wgs84.tif",sep="_"))
    vTif<-paste0(dayDir,"/uv_wind/",paste(DateFile,co,func,"vwind_wgs84.tif",sep="_"))
    modU<-terra::extract(rast(uTif),pts,method="bilinear")[,2]
    modV<-terra::extract(rast(vTif),pts,method="bilinear")[,2]
    modDir<-uv2dir(modU,modV)
  }

  pairs<-data.frame(staDaily,modSpeed=modWS,modDir=modDir)
  stats<-windValStats(obsWS=pairs$obsSpeed,obsDir=pairs$obsDir,
                      modWS=pairs$modSpeed,modDir=pairs$modDir)

  #write the pairs file next to the day metadata
  valDir<-paste0(dayDir,"/validation")
  dir.create(valDir,showWarnings = F,recursive = T)
  pairsName<-paste0(DateFile,"_",co,"_",func,"_era5_validation_pairs.csv")
  write.csv(pairs,paste0(valDir,"/",pairsName),row.names = F)
  message(paste(pairsName,"daily validation pairs written"))

  return(list(stats=stats,pairs=pairs,nStations=nrow(pairs)))
}

metaMakerCo<-function(co, DateIn, func, dataDir, valStats, nValStations){
  require(raster)

  #make date
  DateIn<-as.Date(DateIn)
  DateFile<-format(DateIn,'%Y%m%d')

  #get hourly files to count
  DateDir<-paste(filedir,DateFile,sep="/")
  tifFilesHrs<-Sys.glob(paste0(DateDir,"/*/",co,"/spd_dir_wind/*vel*.tif"))
  tifFilesDay<-Sys.glob(paste0(DateDir,"/day/",co,"/*/*.tif"))
  dayFunc<-tifFilesDay[grep(func,tifFilesDay)] #save day func ws files
  ULdayFunc<-unlist(strsplit(gsub(DateDir,"",dayFunc), "/", fixed = FALSE)) #subset only file name
  dayFuncFiles<-ULdayFunc[grep(".tif",ULdayFunc)]

  #get only hrly files
  tifFilesHrs<-tifFilesHrs[grep("day",tifFilesHrs,invert = T)]
  ULtifFilesHrs<-unlist(strsplit(gsub(DateDir,"",tifFilesHrs), "/", fixed = FALSE)) #subset only file name
  hrlyWStiffs<-ULtifFilesHrs[grep(".tif",ULtifFilesHrs)]

  #get daily vel tif
  dayWS<-dayFunc[grep("vel_mps",dayFunc)]
  dayRaster<-raster(dayWS)
  minMap<-round(cellStats(dayRaster,min),2)
  maxMap<-round(cellStats(dayRaster,max),2)

  statement<-paste("These", format(DateIn, "%B %d %Y"),
                   "daily Hawaii wind maps of",
                   domainHI(co),if(co=="OA"){"(Honolulu City and County)"},
                   "county are a high spatial resolution (~250m) gridded wind modeled prediction of 10 meter wind speed, direction, and the conversion of these variables to u-wind and v-wind components. Hourly maps were produced by a weather-model initialization (wxModelInitialization) run of the windninja (https://ninjastorm.firelab.org/windninja/) land surface wind model, initialized with ERA5 reanalysis surface fields (10 m u/v wind, 2 m temperature, total cloud cover) converted to a WRF-format forecast file, with the diurnal slope-flow model enabled, to account for land cover and topographic influence on wind speed and direction. NO station data was used as model input.",
                   "The", func, "daily (", DateIn,
                   ") wind maps were produced by calculating the per pixel",
                   ifelse(func=="max","maximum",ifelse(func=="mean","average",func)),
                   "wind values from all available hourly wind maps data (n=",
                   length(hrlyWStiffs),
                   ") within the target day and county domain. Daily cross-validation compares, for each of n=",
                   nValStations,
                   "raw validation stations (basic range limits only, no outlier filtering, reporting >=",
                   minHrsDaily,
                   "of 24 hours, never used as model input), the station's own daily",func,
                   "computed from its hourly observations against the daily",func,
                   "product grid sampled at the station location (direction for mean/median from component-wise U/V statistics; direction for max taken at the hour of maximum speed).",
                   valStatement(valStats),
                   "Validation metric definitions: RMSE = root mean square error; BIAS = mean error (model minus observed); RME = combined U/V wind component bias sqrt(biasU^2+biasV^2); R-SQUARE = squared Pearson correlation for speed and squared Jammalamadaka-SenGupta circular correlation for direction; direction errors wrapped to +/-180 degrees and calm observations excluded from direction metrics. Consult hourly wind map metadata for hourly cross-validation statistics. These wind map data are an approximation, use data with caution.")
  #Make meta df wide
  metaWide<-data.frame(
    dataStatement=statement,
    keywords=paste0(domainHI(co), ", Hawaiian Islands, daily ",func," wind, hourly wind, wind speed, wind direction, u-wind, v-wind, ERA5, cross-validation"),
    county=co,
    dataStartDateTime=paste(DateIn,"00:00:00 HST"),
    dataEndDateTime=paste(DateIn,"23:59:59 HST"),
    dataDate=format(DateIn,"%Y-%m-%d"),
    dateProduced=Sys.Date(),
    dataVersionType="experimental",
    windHeight="10",
    windHeightUnit="meters",
    windDirUnit="degree",
    windSpeedUnit="mps",
    windUVunit="mps",
    dailyStatistic=func,
    hourlyWindMapCount=length(hrlyWStiffs),
    hourlyWindMapIncluded=paste(substr(hrlyWStiffs,10,13),collapse=","), #get hrs of day used in stat
    gridWindMinMPS=minMap,
    gridWindMaxMPS=maxMap,
    windGridFiles=paste(dayFuncFiles,collapse = ", "),
    GeoCoordUnits="Decimal Degrees",
    GeoCoordRefSystem="EPSG:4326",
    Xresolution=res(dayRaster)[1],
    Yresolution=res(dayRaster)[2],
    ExtentXmin=extent(dayRaster)[1],
    ExtentXmax=extent(dayRaster)[2],
    ExtentYmin=extent(dayRaster)[3],
    ExtentYmax=extent(dayRaster)[4],
    modelInitialization="wxModelInitialization: ERA5 reanalysis via WRF-mimic NetCDF (era5_to_wrfout.R), diurnal_winds=true",
    validationStationSet="all raw stations, basic range limits only, no outlier filter; independent of model input",
    valStationCount=nValStations,
    valMinHoursPerStation=minHrsDaily,
    valPairing="station daily statistic (from hourly obs) vs daily product grid at station location",
    credits="All data produced by University of Hawaii at Manoa Water Resource Research Center (WRRC). Support provided by Change Hawaii EPSCoR funded by the National Science Foundation under EPSCoR Research Infrastructure Improvement Award #OIA-2149133",
    contacts="Matthew Lucas (mplucas@hawaii.edu), Keri Kodama (kodamak8@hawaii.edu), Sayed Bateni (smbateni@hawaii.edu), Ryan Longman (rlongman@hawaii.edu), Thomas Giambelluca (thomas@hawaii.edu)"
  )
  final_meta<-rbind(data.frame(attribute=as.character(names(metaWide)),value=as.character(lapply(metaWide[1,], dataconvert))))
  #append the daily cross-validation metric rows (prefixed with the statistic)
  final_meta<-rbind(final_meta,valAttrRows(valStats,prefix=paste0("daily",toupper(substr(func,1,1)),substr(func,2,nchar(func)))))
  return(final_meta)
}

statewideMeta<-function(func,DateIn,filedir,writeFile=TRUE){

  #make dirs and file
  DateIn<-as.Date(DateIn)
  DateFile<-format(DateIn,'%Y%m%d')
  DateDir<-paste(filedir,gsub("_","/",DateFile),sep="/")
  stateDir<-paste(DateDir,"day/statewide",sep="/")
  outDirMeta<-paste(stateDir,"metadata",sep="/")

  require(raster)

  #get a statewide tif files
  stateDir<-paste0(DateDir,"/day/statewide")
  tifPaths<-Sys.glob(paste0(stateDir,"/*/*",func,"*.tif"))
  tifFilesUL<-unlist(strsplit(tifPaths, "/", fixed = FALSE))
  tifFiles<-tifFilesUL[grep("tif",tifFilesUL)] #subset only file names
  stateVelTiff<-tifPaths[grep("vel_mps",tifPaths)]
  stateWSRater<-raster(stateVelTiff)
  minMap<-round(cellStats(stateWSRater,min),2)
  maxMap<-round(cellStats(stateWSRater,max),2)

  #get county metadata files and some data from them (BY NAME, not row number)
  metaFiles<-Sys.glob(paste0(DateDir,"/day/*/metadata/",DateFile,"*",func,"_wind_metadata.txt"))
  meta_lists <- lapply(metaFiles, read.delim)
  counties<-paste(sapply(meta_lists, function(df) metaGet(df,"county")),collapse = ", ")
  countyNames<-paste(sapply(meta_lists, function(df) domainHI(metaGet(df,"county"))),collapse = ", ")
  coWindMapCount<-paste(sapply(meta_lists, function(df) metaGet(df,"hourlyWindMapCount")),collapse = ", ")
  minMapCO<-paste(sapply(meta_lists, function(df) metaGet(df,"gridWindMinMPS")),collapse = ", ")
  maxMapCO<-paste(sapply(meta_lists, function(df) metaGet(df,"gridWindMaxMPS")),collapse = ", ")
  coValN<-paste(sapply(meta_lists, function(df) metaGet(df,"valStationCount")),collapse = ", ")

  #pooled statewide daily cross-validation from the county daily pair files
  pairFiles<-Sys.glob(paste0(DateDir,"/day/*/validation/",DateFile,"_*_",func,"_era5_validation_pairs.csv"))
  if(length(pairFiles)>0){
    allPairs<-do.call(rbind,lapply(pairFiles,read.csv,stringsAsFactors=FALSE))
  }else{
    allPairs<-data.frame(obsSpeed=numeric(0),obsDir=numeric(0),modSpeed=numeric(0),modDir=numeric(0))
  }
  valStats<-windValStats(obsWS=allPairs$obsSpeed,obsDir=allPairs$obsDir,
                         modWS=allPairs$modSpeed,modDir=allPairs$modDir)
  message(paste("statewide daily",func,"cross-validation:",DateIn,"n =",valStats$n))

  statement<-paste("These", format(DateIn, "%B %d %Y"),
                   "daily Hawaii wind maps of",
                   countyNames, paste0("counties ( ",counties," ),"),
                   "are a high spatial resolution (~250m) gridded wind modeled prediction of 10 meter wind speed, direction, and the conversion of these variables to u-wind and v-wind components. Hourly county maps were produced by a weather-model initialization (wxModelInitialization) run of the windninja (https://ninjastorm.firelab.org/windninja/) land surface wind model, initialized with ERA5 reanalysis surface fields converted to a WRF-format forecast file, with the diurnal slope-flow model enabled. NO station data was used as model input.",
                   "The",func, "daily (", DateIn,
                   ") wind maps were produced by calculating the per pixel",
                   ifelse(func=="max","maximum",ifelse(func=="mean","average",func)),
                   "wind values from all available county (",
                   counties,
                   ") hourly wind maps data (n=",
                   paste(coWindMapCount,collapse = ", "),
                   ") respectively within the target day and county domains. Daily cross-validation compares each raw validation station's own daily",func,
                   "(computed from its hourly observations; stations reporting >=",minHrsDaily,
                   "of 24 hours; per-county n=",coValN,
                   ") against the daily product grids sampled at the station locations, pooled statewide.",
                   valStatement(valStats),
                   "Validation metric definitions: RMSE = root mean square error; BIAS = mean error (model minus observed); RME = combined U/V wind component bias sqrt(biasU^2+biasV^2); R-SQUARE = squared Pearson correlation for speed and squared Jammalamadaka-SenGupta circular correlation for direction; direction errors wrapped to +/-180 degrees and calm observations excluded from direction metrics. Consult county daily metadata for per-county statistics. These wind map data are an approximation, use data with caution.")
  #Make meta df wide
  metaWide<-data.frame(
    dataStatement=statement,
    keywords=paste0(countyNames, ", Hawaiian Islands, wind, daily ",func," wind, wind speed, wind direction, u-wind, v-wind, ERA5, cross-validation"),
    county=counties,
    dataStartDateTime=paste(DateIn,"00:00:00 HST"),
    dataEndDateTime=paste(DateIn,"23:59:59 HST"),
    dataDate=format(DateIn,"%Y-%m-%d"),
    dateProduced=Sys.Date(),
    dataVersionType="experimental",
    windHeight="10",
    windHeightUnit="meters",
    windDirUnit="degree",
    windSpeedUnit="mps",
    windUVunit="mps",
    dailyStatistic=func,
    hourlyWindMapCount=paste(coWindMapCount,collapse = ", "),
    gridCountyWindMinMPS=paste(minMapCO,collapse = ", "),
    gridCountyWindMaxMPS=paste(maxMapCO,collapse = ", "),
    gridStateWindMinMPS=minMap,
    gridStateWindMaxMPS=maxMap,
    windGridFiles=paste(tifFiles,collapse = ", "),
    GeoCoordUnits="Decimal Degrees",
    GeoCoordRefSystem="EPSG:4326",
    Xresolution=res(stateWSRater)[1],
    Yresolution=res(stateWSRater)[2],
    ExtentXmin=extent(stateWSRater)[1],
    ExtentXmax=extent(stateWSRater)[2],
    ExtentYmin=extent(stateWSRater)[3],
    ExtentYmax=extent(stateWSRater)[4],
    modelInitialization="wxModelInitialization: ERA5 reanalysis via WRF-mimic NetCDF (era5_to_wrfout.R), diurnal_winds=true",
    validationStationSet="all raw stations, basic range limits only, no outlier filter; pooled across counties; independent of model input",
    valStationCountByCounty=coValN,
    valMinHoursPerStation=minHrsDaily,
    valPairing="station daily statistic (from hourly obs) vs daily product grid at station location",
    credits="All data produced by University of Hawaii at Manoa Water Resource Research Center (WRRC). Support provided by Change Hawaii EPSCoR funded by the National Science Foundation under EPSCoR Research Infrastructure Improvement Award #OIA-2149133",
    contacts="Matthew Lucas (mplucas@hawaii.edu), Keri Kodama (kodamak8@hawaii.edu), Sayed Bateni (smbateni@hawaii.edu), Ryan Longman (rlongman@hawaii.edu), Thomas Giambelluca (thomas@hawaii.edu)"
  )

  final_meta<-rbind(data.frame(attribute=as.character(names(metaWide)),value=as.character(lapply(metaWide[1,], dataconvert))))
  #append pooled statewide daily validation metric rows
  final_meta<-rbind(final_meta,valAttrRows(valStats,prefix=paste0("daily",toupper(substr(func,1,1)),substr(func,2,nchar(func)))))

  if(writeFile){
    dir.create(outDirMeta,showWarnings = F)
    setwd(outDirMeta)
    metaName<-paste0(DateFile,"_statewide_",func,"_wind_metadata.txt")
    write.table(final_meta,file=metaName,quote=F,sep = "\t",row.names = FALSE)
    message(paste(metaName, "statewide metadata file written to",outDirMeta))
    return(metaName)
  }else{
    return(final_meta)
  }
}

hrlyWind2Daily<-function(DateIn,filedir,domain,func,writeFile=TRUE){
  require(terra)
  DateIn<-as.Date(DateIn)
  DateFile<-format(DateIn,'%Y%m%d')
  DateDir<-paste(filedir,DateFile,sep="/")
  tifFilesHrs<-Sys.glob(paste(DateDir,"*",domain,"*","*.tif",sep="/"))
  tifFilesHrs<-tifFilesHrs[grep("/day/",tifFilesHrs,invert = T)]#remove day for re-runing

  #uwind
  uFilesHrs<-tifFilesHrs[grep("uwind_wgs84.tif",tifFilesHrs)]
  uRastHrs <- rast(uFilesHrs)
  dayUwindFunc <- app(uRastHrs, fun = func)
  names(dayUwindFunc)<-paste(DateFile,domain,func,"uwind_wgs84.tif",sep="_" )

  #vwind
  vFilesHrs<-tifFilesHrs[grep("vwind_wgs84.tif",tifFilesHrs)]
  vRastHrs <- rast(vFilesHrs)
  dayVwindFunc <- app(vRastHrs, fun = func)
  names(dayVwindFunc)<-paste(DateFile,domain,func,"vwind_wgs84.tif",sep="_" )

  #ws vel
  wsFilesHrs<-tifFilesHrs[grep("vel_mps_wgs84.tif",tifFilesHrs)]
  wsRastHrs <- rast(wsFilesHrs)
  dayVelWindFunc <- app(wsRastHrs, fun = func)
  names(dayVelWindFunc)<-paste(DateFile,domain,func,"vel_mps_wgs84.tif",sep="_" )

  if(func=="max"){
    dirFilesHrs<-tifFilesHrs[grep("dir_deg_wgs84.tif",tifFilesHrs)]
    dirRastHrs <- rast(dirFilesHrs)

    # Get values as matrices (each row = pixel, each column = hour)
    ws_vals <- values(wsRastHrs)
    dir_vals <- values(dirRastHrs)

    # Find the column index (hour) with max wind speed for each pixel
    max_indices <- max.col(ws_vals, ties.method = "first")

    # Extract corresponding directions using row-column indexing
    dir_at_max <- dir_vals[cbind(seq_len(nrow(dir_vals)), max_indices)]

    # Create output max speed direction raster with same properties as input
    dayDirWindFunc <- rast(wsRastHrs, nlyrs = 1) #template raster
    values(dayDirWindFunc) <- dir_at_max #fill with values
    names(dayDirWindFunc)<-paste(DateFile,domain,func,"dir_deg_wgs84.tif",sep="_" )

    #make max ws hr map
    dayMaxWindHr <- which.max(wsRastHrs)
    dayMaxWindHr <- dayMaxWindHr-1 #remove 1 hr to start with hour 0
    names(dayMaxWindHr)<-paste(DateFile,domain,func,"vel_hour_wgs84.tif",sep="_" )

    }else{
      #convert UV to dir
      wsDirFunc<-UV2wsDir(u=dayUwindFunc,v=dayVwindFunc)
      dayDirWindFunc<-wsDirFunc[["dirDeg"]]
      names(dayDirWindFunc)<-paste(DateFile,domain,func,"dir_deg_wgs84.tif",sep="_" )
  }
  if(writeFile){
      #Make daily out dirs
      outDir<-paste0(filedir,"/",DateFile,"/day")
      dir.create(outDir,showWarnings = F)
      outDir<-paste0(outDir,"/",domain)
      dir.create(outDir,showWarnings = F)
      varDirs<-c("spd_dir_wind","uv_wind","metadata","max_hr")
      outDirs<-paste(outDir,varDirs,sep="/")
      dir.create(outDirs[1],showWarnings = F)
      dir.create(outDirs[2],showWarnings = F)
      dir.create(outDirs[3],showWarnings = F)
      dir.create(outDirs[4],showWarnings = F)

      #write files
      writeRaster(dayDirWindFunc,paste(outDirs[1],names(dayDirWindFunc),sep="/"),overwrite=TRUE)
      message(paste(names(dayDirWindFunc),"written to",outDirs[1]))
      writeRaster(dayVelWindFunc,paste(outDirs[1],names(dayVelWindFunc),sep="/"),overwrite=TRUE)
      message(paste(names(dayVelWindFunc),"written to",outDirs[1]))
      writeRaster(dayUwindFunc,paste(outDirs[2],names(dayUwindFunc),sep="/"),overwrite=TRUE)
      message(paste(names(dayUwindFunc),"written to",outDirs[2]))
      writeRaster(dayVwindFunc,paste(outDirs[2],names(dayVwindFunc),sep="/"),overwrite=TRUE)
      message(paste(names(dayVwindFunc),"written to",outDirs[2]))

      if(func=="max"){
        writeRaster(dayMaxWindHr,paste(outDirs[4],names(dayMaxWindHr),sep="/"),overwrite=TRUE)
        message(paste(names(dayMaxWindHr),"written to",outDirs[4]))
      }

      #daily cross-validation (station daily stat vs grid daily stat)
      val<-valDailyCo(DateIn=DateIn,co=domain,func=func,filedir=filedir,wnInput=wnInput)

      #make metadata for county day run
      metaDF<-metaMakerCo(co=domain, DateIn=DateIn, func=func, dataDir=filedir,
                          valStats=val$stats, nValStations=val$nStations)
      setwd(outDirs[3])
      metaName<-paste0(DateFile,"_",domain,"_",func,"_wind_metadata.txt")
      write.table(metaDF,file=metaName,quote=F,sep = "\t",row.names = FALSE)
      message(paste(metaName,"written to",outDirs[3]))

      #final write message
      message(paste(domain, DateIn, "all final wind files written!"))
  }
  return(paste(DateIn,domain,sep="_"))
}#end make daily rasters

dayfiles2mosaic<-function(var,func,DateIn,filedir,writeFile=TRUE){
  require(terra)
  dateFile<-format(as.Date(DateIn),'%Y%m%d')
  dateDir<-paste0(filedir,"/",dateFile,"/day")
  funcFiles<-Sys.glob(paste0(dateDir,"/*/*/*",func,"*.tif"))
  var_unit<-ifelse(var=="dir","dir_deg",ifelse(var=="vel","vel_mps",ifelse(var=="hour","vel_hour",var)))
  varFiles<-funcFiles[grep(var_unit,funcFiles)]
  raster_list <- lapply(varFiles, rast)
  raster_collection <- sprc(raster_list)
  mosaiced_raster <- mosaic(raster_collection, fun = mean)
  names(mosaiced_raster)<-paste0(dateFile,"_statewide_",func,"_",var_unit,"_wgs84.tif")
  if(writeFile){
    outDir<-paste(dateDir,"statewide",sep="/")
    dir.create(outDir,showWarnings = F)
    varDir<-ifelse(var=="vel"|var=="dir","spd_dir_wind",ifelse(var=="vwind"|var=="uwind","uv_wind",if(var=="hour"){"max_hr"}))
    outDir<-paste(outDir,varDir,sep="/")
    dir.create(outDir,showWarnings = F)
    writeRaster(mosaiced_raster,paste(outDir,names(mosaiced_raster),sep="/"),overwrite=TRUE)
    message(paste(names(mosaiced_raster),"written to",outDir))
  }
  return(names(mosaiced_raster))
}

## MEAN DAILY WIND ##
#county run hrly to daily mean all variables
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="BI",func="mean",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="MN",func="mean",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="OA",func="mean",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="KA",func="mean",writeFile=T)

#per daily mean variable to statewide
dayfiles2mosaic(var="dir",func="mean",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="vel",func="mean",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="uwind",func="mean",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="vwind",func="mean",DateIn=DateIn,filedir=filedir)

#statewide mean daily metadata
statewideMeta(func="mean",DateIn=DateIn,filedir=filedir)

## MEDIAN DAILY WIND ##
#county run hrly to daily median all variables
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="BI",func="median",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="MN",func="median",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="OA",func="median",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="KA",func="median",writeFile=T)

#per daily median variable to statewide
dayfiles2mosaic(var="dir",func="median",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="vel",func="median",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="uwind",func="median",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="vwind",func="median",DateIn=DateIn,filedir=filedir)

#statewide median daily metadata
statewideMeta(func="median",DateIn=DateIn,filedir=filedir)

## MAX DAILY WIND ##
#county run hrly to daily max all variables
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="BI",func="max",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="MN",func="max",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="OA",func="max",writeFile=T)
hrlyWind2Daily(DateIn=DateIn,filedir=filedir,domain="KA",func="max",writeFile=T)

#per daily max variable to statewide
dayfiles2mosaic(var="dir",func="max",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="vel",func="max",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="uwind",func="max",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="vwind",func="max",DateIn=DateIn,filedir=filedir)
dayfiles2mosaic(var="hour",func="max",DateIn=DateIn,filedir=filedir)

#statewide max daily metadata
statewideMeta(func="max",DateIn=DateIn,filedir=filedir)

print(paste(DateIn,"ERA5 statewide daily wind + cross-validation complete!"))
#pau
