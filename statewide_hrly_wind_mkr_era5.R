###make hrly statewide ws and dir rasters and meta - ERA5 wx-model edition
# ERA5 edition of statewide_hrly_wind_mkr.R (HCDP/windninja_hrly_wind). Changes:
#   - filedir -> windProducts_era5 (env-overridable)
#   - county metadata values read BY ATTRIBUTE NAME (original used row numbers,
#     which breaks when attributes are added)
#   - statewide hourly metadata gains pooled cross-validation metrics computed
#     from the four county validation-pair files for the hour

#define main dir
filedir<-Sys.getenv("WN_ERA5_PRODUCTS","/home/wn1/data/output/windProducts_era5")
codeDir<-Sys.getenv("WN_CODE_DIR","/home/wn1/wn_codes/era5_wrf_TEST")
source(file.path(codeDir,"wind_validation_fx.R"))

#custom functions
dataconvert<-function(x){
  return(ifelse(is.numeric(x),as.character(round(x,5)),as.character(x)))
}

domainHI<-function(co,county=FALSE){
  domain<-if(co=="BI"){"Hawaii"}else
    if(co=="MN"){"Maui"}else
      if(co=="OA"){"Oahu"}else
        if(co=="KA"){"Kauai"}else
          if(co=="HI"){"Hawaii"}else{"NA"}
  if(county){
    domain<-ifelse(domain=="Oahu","Oahu (Honolulu City and County)",domain)
  }
  return(domain)
}

DTfiles2mosaic<-function(var,DTin,filedir,writeFile=TRUE){
  require(terra)
  DTFile<-format(as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone()),'%Y%m%d_%H%M')
  DTDir<-paste(filedir,gsub("_","/",DTFile),sep="/")
  tifFiles<-Sys.glob(paste0(DTDir,"/*/*/",DTFile,"*.tif"))
  varFiles<-tifFiles[grep(var,tifFiles)]
  raster_list <- lapply(varFiles, rast)
  raster_collection <- sprc(raster_list)
  mosaiced_raster <- mosaic(raster_collection, fun = mean)
  var_unit<-ifelse(var=="dir","dir_deg",ifelse(var=="vel","vel_mps",var))
  names(mosaiced_raster)<-paste0(DTFile,"_statewide_",var_unit,"_wgs84.tif")
  if(writeFile){
    outDir<-paste(DTDir,"statewide",sep="/")
    dir.create(outDir,showWarnings = F)
    varDir<-ifelse(var=="vel"|var=="dir","spd_dir_wind",if(var=="vwind"|var=="uwind"){"uv_wind"})
    outDir<-paste(outDir,varDir,sep="/")
    dir.create(outDir,showWarnings = F)
    writeRaster(mosaiced_raster,paste(outDir,names(mosaiced_raster),sep="/"),overwrite=TRUE)
    message(paste(names(mosaiced_raster),"written to",outDir))
  }
  return(names(mosaiced_raster))
}

statewideMeta<-function(DTin,filedir,writeFile=TRUE){

  #make dirs and file
  DTinre<-as.POSIXlt(DTin,format="%Y-%m-%d_%H:%M:%S",tz=Sys.timezone())
  DTFile<-format(DTinre,'%Y%m%d_%H%M')
  DTDir<-paste(filedir,gsub("_","/",DTFile),sep="/")
  stateDir<-paste(DTDir,"statewide",sep="/")
  outDirMeta<-paste(stateDir,"metadata",sep="/")

  require(raster)

  #get a statewide tif files
  stateDir<-paste0(DTDir,"/statewide")
  tifPaths<-Sys.glob(paste0(stateDir,"/*/*.tif"))
  tifFilesUL<-unlist(strsplit(tifPaths, "/", fixed = FALSE))
  tifFiles<-tifFilesUL[grep("tif",tifFilesUL)] #subset only file names
  stateRater<-raster(tifPaths[1])

  #get county metadata files and some data from them (BY NAME, not row number)
  metaFiles<-Sys.glob(paste0(DTDir,"/*/metadata/",DTFile,"*_wind_metadata.txt"))
  meta_lists <- lapply(metaFiles, read.delim)
  counties<-paste(sapply(meta_lists, function(df) metaGet(df,"county")),collapse = ", ")
  countyNames<-paste(sapply(meta_lists, function(df) domainHI(metaGet(df,"county"),county=T)),collapse = ", ")
  staCounts<-paste(sapply(meta_lists, function(df) metaGet(df,"stationCount")),collapse = ", ")
  minSta<-paste(sapply(meta_lists, function(df) metaGet(df,"staWindMinMPS")),collapse = ", ")
  maxSta<-paste(sapply(meta_lists, function(df) metaGet(df,"staWindMaxMPS")),collapse = ", ")
  minMap<-paste(sapply(meta_lists, function(df) metaGet(df,"gridWindMinMPS")),collapse = ", ")
  maxMap<-paste(sapply(meta_lists, function(df) metaGet(df,"gridWindMaxMPS")),collapse = ", ")

  #pooled statewide cross-validation from the county validation pair files
  pairFiles<-Sys.glob(paste0(DTDir,"/*/stationData/",DTFile,"*_era5_validation_pairs.csv"))
  if(length(pairFiles)>0){
    allPairs<-do.call(rbind,lapply(pairFiles,read.csv,stringsAsFactors=FALSE))
  }else{
    allPairs<-data.frame(obsSpeed=numeric(0),obsDir=numeric(0),modSpeed=numeric(0),modDir=numeric(0))
  }
  valStats<-windValStats(obsWS=allPairs$obsSpeed,obsDir=allPairs$obsDir,
                         modWS=allPairs$modSpeed,modDir=allPairs$modDir)
  message(paste("statewide hourly cross-validation:",DTin,"n =",valStats$n))

  statement<-paste("These", format(DTinre, "%B %d %Y %H:%M %P"),
                   "HST hourly Hawaii statewide wind maps of",
                   countyNames, paste0("counties ( ",counties," ),"),
                   "are a high spatial resolution (~250m) gridded wind modeled prediction of 10 meter wind speed, direction, and the conversion of these variables to u-wind and v-wind components. Each of these variables were combined into one map extent by combining all available county (",
                   counties,
                   ") maps into a single mosaiced extent. Each county map was produced by a weather-model initialization (wxModelInitialization) run of the windninja (https://ninjastorm.firelab.org/windninja/) land surface wind model, initialized with ERA5 reanalysis surface fields (10 m u/v wind, 2 m temperature, total cloud cover) converted to a WRF-format forecast file, with the diurnal slope-flow model enabled. NO station data was used as model input; all available raw wind station observations per county ( n=",
                   staCounts,
                   ", basic range limits only, no outlier filtering) were used exclusively as independent cross-validation truth. Minimum and maximum wind speed observed at validation stations for each county (",
                   counties,
                   ") were:",
                   minSta, "mps and", maxSta,
                   "mps respectively with modeled map minimum and maximum wind speeds for each county (",
                   counties,
                   ") being:",
                   minMap, "mps and", maxMap,
                   "mps respectively.",
                   valStatement(valStats),
                   "Validation metric definitions: RMSE = root mean square error; BIAS = mean error (model minus observed); RME = combined U/V wind component bias sqrt(biasU^2+biasV^2); R-SQUARE = squared Pearson correlation for speed and squared Jammalamadaka-SenGupta circular correlation for direction; direction errors wrapped to +/-180 degrees and calm observations excluded from direction metrics. These wind map data are an approximation, use data with caution.")
  #Make meta df wide
  metaWide<-data.frame(
    dataStatement=statement,
    keywords=paste0(countyNames, ", Hawaiian Islands, wind, hourly wind, wind speed, wind direction, u-wind, v-wind, ERA5, cross-validation"),
    county=counties,
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
    stationCount=staCounts,
    staWindMinMPS=minSta,
    staWindMaxMPS=maxSta,
    gridWindMinMPS=minMap,
    gridWindMaxMPS=maxMap,
    windGridFiles=paste(tifFiles,collapse = ", "),
    GeoCoordUnits="Decimal Degrees",
    GeoCoordRefSystem="EPSG:4326",
    Xresolution=res(stateRater)[1],
    Yresolution=res(stateRater)[2],
    ExtentXmin=extent(stateRater)[1],
    ExtentXmax=extent(stateRater)[2],
    ExtentYmin=extent(stateRater)[3],
    ExtentYmax=extent(stateRater)[4],
    modelInitialization="wxModelInitialization: ERA5 reanalysis via WRF-mimic NetCDF (era5_to_wrfout.R), diurnal_winds=true",
    validationStationSet="all raw stations, basic range limits only, no outlier filter; pooled across counties; independent of model input",
    credits="All data produced by University of Hawaii at Manoa Water Resource Research Center (WRRC). Support provided by Change Hawaii EPSCoR funded by the National Science Foundation under EPSCoR Research Infrastructure Improvement Award #OIA-2149133",
    contacts="Matthew Lucas (mplucas@hawaii.edu), Keri Kodama (kodamak8@hawaii.edu), Sayed Bateni (smbateni@hawaii.edu), Ryan Longman (rlongman@hawaii.edu), Thomas Giambelluca (thomas@hawaii.edu)"
  )

  final_meta<-rbind(data.frame(attribute=as.character(names(metaWide)),value=as.character(lapply(metaWide[1,], dataconvert))))
  #append pooled statewide validation metric rows
  final_meta<-rbind(final_meta,valAttrRows(valStats))

  if(writeFile){
    dir.create(outDirMeta,showWarnings = F)
    setwd(outDirMeta)
    metaName<-paste0(DTFile,"_statewide_wind_metadata.txt")
    write.table(final_meta,file=metaName,quote=F,sep = "\t",row.names = FALSE)
    message(paste(metaName, "statewide metadata file written to",outDirMeta))
  }
  return(final_meta)
}

#get datetime in
args <- commandArgs(trailingOnly = TRUE)
DTin=args[1]

#run statewide mosaic for all 4 vars
DTfiles2mosaic(var="dir",DTin=DTin, filedir=filedir)
DTfiles2mosaic(var="vel",DTin=DTin, filedir=filedir)
DTfiles2mosaic(var="uwind",DTin=DTin, filedir=filedir)
DTfiles2mosaic(var="vwind",DTin=DTin, filedir=filedir)
message("all statewide hrly wind tifs written")

#make statewide hrly metadata
statewideMeta(DTin,filedir)
message("statewide hrly wind meta written! Statewide process finished.")

#pau statewide hrly wind process
