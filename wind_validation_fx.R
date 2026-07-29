# Shared wind cross-validation functions for the ERA5 WindNinja workflow.
# Sourced by wspd_dir_2_uv_mkr_era5.R, statewide_hrly_wind_mkr_era5.R, dailyWind_mkr_era5.R.
#
# Metric definitions (per user specification):
#   RMSE      - root mean square error (speed: m/s; direction: degrees on wrapped errors)
#   BIAS      - mean error (speed: m/s; direction: signed circular mean of wrapped errors)
#   RME       - combined U/V bias: sqrt(biasU^2 + biasV^2), the magnitude of the mean
#               vector error, m/s (bias of U wind and V wind then combined)
#   R-SQUARE  - speed: ordinary Pearson r^2;
#               direction: squared Jammalamadaka-SenGupta circular correlation
# Direction errors are circular: differences wrapped to +/-180 degrees.
# Pairs with calm observed wind (speed <= calmThresh) are excluded from the
# direction metrics only (direction is undefined at calm); they remain in the
# speed and U/V metrics.

wrap180<-function(x){((x+180)%%360)-180}

#met convention: dir = direction wind blows FROM; +U = wind TO east, +V = wind TO north
met2u<-function(ws,dir){-ws*sin(dir*pi/180)}
met2v<-function(ws,dir){-ws*cos(dir*pi/180)}
uv2ws<-function(u,v){sqrt(u^2+v^2)}
uv2dir<-function(u,v){(270-atan2(v,u)*(180/pi))%%360}

circMeanDeg<-function(deg){ #signed circular mean in +/-180 degrees
  r<-deg*pi/180
  atan2(mean(sin(r)),mean(cos(r)))*180/pi
}

circCorJS<-function(a,b){ #Jammalamadaka-SenGupta circular correlation, inputs degrees
  a<-a*pi/180; b<-b*pi/180
  am<-atan2(mean(sin(a)),mean(cos(a)))
  bm<-atan2(mean(sin(b)),mean(cos(b)))
  num<-sum(sin(a-am)*sin(b-bm))
  den<-sqrt(sum(sin(a-am)^2)*sum(sin(b-bm)^2))
  if(!is.finite(den)||den==0){return(NA)}
  num/den
}

windValStats<-function(obsWS,obsDir,modWS,modDir,calmThresh=0){
  #assemble pairs, drop rows with any NA
  ok<-is.finite(obsWS)&is.finite(obsDir)&is.finite(modWS)&is.finite(modDir)
  obsWS<-obsWS[ok]; obsDir<-obsDir[ok]; modWS<-modWS[ok]; modDir<-modDir[ok]
  n<-length(obsWS)
  out<-list(n=n,nDir=NA,
            spdRMSE=NA,spdBIAS=NA,spdR2=NA,
            dirRMSE=NA,dirBIAS=NA,dirR2=NA,
            biasU=NA,biasV=NA,RME=NA)
  if(n==0){return(out)}

  #speed
  dS<-modWS-obsWS
  out$spdRMSE<-sqrt(mean(dS^2))
  out$spdBIAS<-mean(dS)
  if(n>=3 && sd(obsWS)>0 && sd(modWS)>0){
    out$spdR2<-cor(obsWS,modWS)^2
  }

  #U/V component bias -> combined RME
  dU<-met2u(modWS,modDir)-met2u(obsWS,obsDir)
  dV<-met2v(modWS,modDir)-met2v(obsWS,obsDir)
  out$biasU<-mean(dU)
  out$biasV<-mean(dV)
  out$RME<-sqrt(out$biasU^2+out$biasV^2)

  #direction (circular), calm obs excluded
  keep<-obsWS>calmThresh
  out$nDir<-sum(keep)
  if(out$nDir>0){
    dD<-wrap180(modDir[keep]-obsDir[keep])
    out$dirRMSE<-sqrt(mean(dD^2))
    out$dirBIAS<-circMeanDeg(dD)
    if(out$nDir>=3){
      r<-circCorJS(obsDir[keep],modDir[keep])
      out$dirR2<-if(is.na(r)){NA}else{r^2}
    }
  }
  return(out)
}

valRound<-function(x){ifelse(is.na(x),"NA",as.character(round(x,3)))}

#metadata rows (attribute/value) from a windValStats list; prefix e.g. "" or "dailyMean"
valAttrRows<-function(stats,prefix=""){
  p<-function(nm){paste0(prefix,nm)}
  data.frame(
    attribute=c(p("valPairCount"),p("valPairCountDir"),
                p("spdRMSE_mps"),p("spdBIAS_mps"),p("spdR2"),
                p("dirRMSE_deg"),p("dirBIAS_deg"),p("dirR2circ"),
                p("biasUwind_mps"),p("biasVwind_mps"),p("RMEuv_mps")),
    value=c(as.character(stats$n),as.character(stats$nDir),
            valRound(stats$spdRMSE),valRound(stats$spdBIAS),valRound(stats$spdR2),
            valRound(stats$dirRMSE),valRound(stats$dirBIAS),valRound(stats$dirR2),
            valRound(stats$biasU),valRound(stats$biasV),valRound(stats$RME)),
    stringsAsFactors=FALSE)
}

#one-line summary for dataStatement text
valStatement<-function(stats){
  if(stats$n==0){return("Cross-validation was attempted but no station observations were available for this product.")}
  paste0("Cross-validation against n=",stats$n," independent station observations ",
         "(raw stations, basic range limits only; stations are NOT model inputs): ",
         "wind speed RMSE=",valRound(stats$spdRMSE)," mps, BIAS=",valRound(stats$spdBIAS),
         " mps, R-SQUARE=",valRound(stats$spdR2),
         "; wind direction (circular, n=",stats$nDir," non-calm) RMSE=",valRound(stats$dirRMSE),
         " deg, BIAS=",valRound(stats$dirBIAS)," deg, R-SQUARE=",valRound(stats$dirR2),
         "; combined U/V bias RME=",valRound(stats$RME)," mps (U bias=",valRound(stats$biasU),
         ", V bias=",valRound(stats$biasV)," mps).")
}

#name-based metadata lookup (replaces fragile row-number indexing)
metaGet<-function(metaDF,attr){
  v<-metaDF$value[match(attr,metaDF$attribute)]
  ifelse(is.null(v),NA,v)
}
