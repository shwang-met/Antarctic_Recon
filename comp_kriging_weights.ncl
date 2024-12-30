;===============================================================================
; This script calculates the kriging weights by using reanalsis data.
; The script calculates the kriging weights (λ) by solving equation A · λ = B.
; Below, the variable L is used for λ. 
; For details about the maths, see Appendix A of Nicolas and Bromwich (2014).
; S.-H. Wang, 2024
;===============================================================================

begin

; Parameters
  reanl   = "era5"
  varName = "t2m"
  yyyy1   = 1979
  yyyy2   = 2022

; month   = "12"
  detrend = True
  maxIter = 10

; Directory/file/variable names
  dataDir     = "../data/"+reanl+"/"
  eraiFname   = reanl+"_"+varName+"_"+yyyy1+"-"+yyyy2+"_ano_60km_"+month+".nc"
  seriesFname = reanl+"_"+varName+"_"+yyyy1+"-"+yyyy2+"_ano_series_"+month+".txt"
  wgtsDir     = "../data/"+reanl+"/"
  wgtsFname   = "wgts_"+varName+"_"+yyyy1+"-"+yyyy2+"_"+month+".nc"

; Read in ERA 3D data (time, lat, lon)
  f        = addfile(dataDir + eraiFname, "r")
  data     = f->$varName$
  dataDims = dimsizes(data)
  ntime    = dataDims(0)
  ny       = dataDims(1)
  nx       = dataDims(2)

; Read in interpolated ERA data at station locations
  ncols      = numAsciiCol(dataDir + seriesFname)
  series     = readAsciiTable(dataDir + seriesFname, ncols, "float", 1)
  series     := series(:,1:)
  seriesDims = dimsizes(series)
  ntime2     = seriesDims(0)
  nsta       = seriesDims(1)

; Also get the station names (not required but handy for debugging)
  header   = readAsciiHead(dataDir + seriesFname, 1)
  staNames = str_split(header, " ")
  staNames := staNames(1:)

; Quick check: make sure the time dimensions of series and data match
  if ntime.ne.ntime2 then
    print("Error: Time dimensions don't match")
    exit
  end if

; Create arrays A(2D), B(1D), and L(1D)
  A = new((/nsta+1, nsta+1/), float)
  B = new(nsta+1, float)
  L = new(nsta, float)

; Detrend ERA data and series before computing the correlations
  if detrend then
    return_info = False
    data   = dtrend_leftdim(data,   return_info)
    series = dtrend_leftdim(series, return_info)
  end if

; Compute the spatial "footprint" of each station
  correl   = new((/nsta, ny, nx/), float)
  do sta = 0, nsta-1
    correl(sta,:,:) = (/ escorc(series(:,sta), data(y|:,x|:,time|:)) /)
  end do
  r2 = correl^2
  correl_sign = correl / abs(correl) ; may cause error when/where correl = 0

; Initialize array A with the cross-correlations between stations
  do sta1 = 0, nsta-1
    do sta2 = 0, nsta-1
        A(sta1,sta2) = (/ escorc(series(:,sta1), series(:,sta2)) /)
    end do
  end do
  A(:,nsta)    = 1.
  A(nsta,:)    = 1.
  A(nsta,nsta) = 0.
  A = A^2 ; convert r to r^2

; Create 3D array to store L values once they are calculated
  wgts   = new((/nsta, ny, nx/), float)
  wgts!0 = "station"
  wgts!1 = "y"
  wgts!2 = "x"

; Loop over all gridpoints
  do i = 0, ny - 1
  do j = 0, nx - 1

;   print("i=" + i + " j=" + j)

    ; Initialize array B with r^2 coeffs for gridpoint (i,j)
    B(:nsta-1) = r2(:,i,j)
    B(nsta)    = 1.
    B_sign     = correl_sign(:,i,j)
    B = B^2 ; convert r to r^2

    ; Initialize array L
    L = 1.

    ; Iterations until L > 0 for all stations. Note that L < 0 is mathematically
    ; possible, but it doesn't make physical sense. The value of maxIter can be
    ; set in the Parameters at the beginning of the script. It generally takes
    ; 4 iterations to get rid of the negative L values.
    iter = 0
    do while(iter.lt.maxIter)

      ; This line is optional. It excludes stations with very small weights
      ; Here, I am using 0.01 as a cut-off value.
      L = where((L.gt.0) .and. (L.lt.0.01), 0., L)

      ; Only keep stations with L > 0.
      i_pos = array_append_record(ind(L.gt.0), nsta, 0)
      nsta_pos = dimsizes(i_pos) - 1

      ; Condition 1: skip 1st iteration (i.e for iter=0, include all stations)
      if iter.ge.1 then
        i_neg = ind(L.lt.0)
        nsta_neg = num(.not.ismissing(i_neg))
        ; Condition 2: still some L < 0
        if (nsta_neg.gt.0) .and. (nsta_pos.gt.1) then
          L = where(L.lt.0, 0., L)
        else
          ; Condition 3: all L > 0
          if (nsta_neg.eq.0) .and. (nsta_pos.gt.1) then
            delete([/i_pos, i_neg/])
            break
          else
            ; Condition 4: only one L > 0 left
            if nsta_pos.eq.1 then
              L = 0.
              L(i_pos(0)) = 1.
              delete([/i_pos, i_neg/])
              break
            end if
          end if
        end if
        delete(i_neg)
      end if

      ; Redefine A and B by keeping only stations with L > 0
      AA = A(i_pos,i_pos)
      BB = B(i_pos)

      ; Solve linear system (λ = A⁻¹ ∙ B)
      LL = solve_linsys(AA, BB)
      L(i_pos(:nsta_pos-1)) = LL(:nsta_pos-1)

      ; Delete variables that may change shape in next iteration
      delete([/i_pos, AA, BB, LL/])

      iter = iter + 1

    end do

    ; Before storing the L values, multiply by sign of correlation coeff.
    wgts(:,i,j) = where(L.gt.0, L * B_sign, 0.)

  end do
  end do

; Write weights to netcdf file
  fout = wgtsDir + wgtsFname
  system("rm -f " + fout)
  ncdf = addfile(fout, "c")
  ncdf->wgts = wgts
  print("Weights saved to " + fout)

end
