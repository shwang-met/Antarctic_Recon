;===============================================================================
; This script interpolates anomalies from a list of stations using 
; the kriging weights
; S.-H. Wang, 2024
;===============================================================================

begin

; Parameters
  reanl   = "era5"
  varName = "t2m"
; month   = "01"

  yr1 = 1958
  yr2 = 2022

; Directory/file/variable names
  statDir     = "../data/station/"
  seriesFname = "stations_"+varName+"_"+yr1+"-"+yr2+"_ano_"+month+".nc"
  stationFname = "stations_list.txt"

  wgtsDir     = "../data/"+reanl+"/"
  wgtsFname   = "wgts_"+varName+"_1979-"+yr2+"_"+month+".nc"

  reconDir    = "../data/"+reanl+"/"
  reconFname  = "recon_"+varName+"_"+yr1+"-"+yr2+"_ano_"+month+".nc"

; Load station names & coordinates
  stations = readAsciiTable(statDir+stationFname,1,"string",1)
  stations := stations(:,0)
  staNames = str_get_field(stations, 1, " ")
  staLat   = stringtofloat(str_get_field(stations, 2, " "))
  staLon   = stringtofloat(str_get_field(stations, 3, " "))
  staLon   = where(staLon .lt. 0, staLon + 360., staLon)
  nsta     = dimsizes(staNames)
  print(staNames)

; Read in the kriging weights (dims: station, nx, ny)
  fw       := addfile(wgtsDir+wgtsFname,"r")
  wgts     = fw->wgts
  wgtsDims = dimsizes(wgts)
  nsta1    = wgtsDims(0)
  ny       = wgtsDims(1)
  nx       = wgtsDims(2)
  wgts@_FillValue = -999.9
; print(wgts)

; Read in the stations time series
  fs       := addfile(statDir+seriesFname,"r")

  series=new((/yr2-yr1+1,nsta/),"float")
  do nn = 0,nsta-1
     series(:,nn)   = fs->$staNames(nn)$
  end do
; print(series)
  seriesDims = dimsizes(series)
  ntime      = seriesDims(0)
  nsta2      = seriesDims(1)

; Missing count
  do iyr = yr1,yr2
     N = num(ismissing(series(iyr-yr1,:)))
     print(iyr+"   "+N)
     if (N.gt.0) then
        miss_id := ind(ismissing(series(iyr-yr1,:)))
        print(staNames(miss_id))
     end if
  end do

; Quick check: make sure the number of stations in wgts and series match
  if nsta1.ne.nsta2 then
      print("Error: Numbers of stations don't match")
      exit
  end if

; Expands the wgts and series arrays so that they have the same shape. This
; allows us to use a vectorized notation (faster) in the next step instead of
; iterating over all (i,j) grid points
  dummy  = new((/ntime,nsta1,ny,nx/),byte) ; the type doesn't matter here
  wgts   := conform(dummy, wgts, (/1,2,3/))
  series := conform(dummy, series, (/0,1/))

; Now, reconstruct the pressure anomaly field
  recon = new((/ntime, nx, ny/), float)
  recon = 0.
  do sta = 0, nsta1 - 1
     incr  = series(:,sta,:,:) * wgts(:,sta,:,:) ; incr = increment
     recon = recon + incr
  end do
; print(recon)

; Add dimension names to the recon array
  f          = addfile(statDir + seriesFname, "r")
  recon!0    = "time"
  recon&time = f->time
  recon!1    = "y"
  recon!2    = "x"
  recon@_FillValue = -999.9

; Write recontructed field to netcdf file
  fout = reconDir + reconFname
  system("rm -f " + fout)
  ncdf = addfile(fout, "c")
  ncdf->recon = recon
  print("Reconstructed data saved to " + fout)

end
