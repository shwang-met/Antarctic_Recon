begin

  diri    = "../data/reader/"

  yrStrt  = 1958
  yrLast  = 2022

  nmo    = 12
  ncol   = nmo+1
  nyrs   = yrLast-yrStrt+1

; nstat  = 16
  nstat  = 12

;*********************************

  yyyymm := yyyymm_time(yrStrt,yrLast,"integer")
  yrfrac := yyyymm_to_yyyyfrac(yyyymm, 0.0)
  print(yrfrac)

  yyyy   := yyyymm/100
  month  := yyyymm - (yyyy*100)
  day    := conform(yyyy, 1, -1)
  hour   := conform(yyyy, 0, -1)
  minute := conform(yyyy, 0, -1)
  sec    := conform(yyyy, 0, -1)
; tunits := "days since "+yyyy(0)+"-1-1 00:00:00"
  tunits = "days since 1951-1-1 00:00:0.0"
  time   := cd_inv_calendar(yyyy,month,day,hour,minute,sec,tunits,0)
  time!0 = "time"
  time&time = time

;*********************************

  stat_id = new(nstat,string)
; stat_id = (/"BYRD","FDAY","ORCA","CASE","HALL", \
;             "SCOT","DAVS","MAWS","SPOL","DUMT", \
;             "MIRY","SYOW","ESPZ","NOVO","VOST","BELG"/)
  stat_id = (/"FDAY","ESPZ","ORCA","NOVO","SYOW","MAWS", \
              "DAVS","MIRY","CASE","DUMT","SCOT","BELG"/)
; print(stat_id)

;*********************************
; Read years and data as 2D array
;*********************************

; fn   = "../data/station/stations_t2m_"+yrStrt+"-"+yrLast+".nc"
  fn   = "../data/station/stations_msl_"+yrStrt+"-"+yrLast+".nc"

  system ("/bin/rm -f "+fn)
  fo = addfile (fn,"c")
  fo->time        = time
; fo->yyyymm      = yyyymm

  do nn = 0,nstat-1

  print(stat_id(nn))

; fi    := stat_id(nn)+".temperature.txt"
  fi    := stat_id(nn)+".msl.txt"

  dummy := readAsciiTable(diri+fi,ncol,"float",1)  
; if the last year doesnot have full year, open following section
; if (yrLast.eq.2012) then
;    dummy := readAsciiTable(diri+fi,ncol,"float",(/1,6/))  
; else
;    dummy := readAsciiTable(diri+fi,ncol,"float",1)  
; end if
  dummy@_FillValue = -999.9
; print(dummy)
  
  year  := dummy(:,0)
  data  := dummy(:,1:nmo)
  print(year)
; print(data)

  data_1d      := ndtooned(data)
  data_1d!0    = "time"
  data_1d&time =  time

  delete_VarAtts(yrfrac,(/"long_name"/))
  yrfrac&time = yyyymm

; print(data_1d)
  printVarSummary(data_1d)
  printVarSummary(yrfrac)

  fo->$stat_id(nn)$ = data_1d

  end do

end
