pro set_mjdrange_vista,mjdrange,period=period,verbose=verbose,help=help
;+
; NAME:
;       SET_MJDRANGE_VISTA
;
; PURPOSE:
;  Set the UKIDSS MJD range by data release(DR), can also report the data range
;
;  To set the date from the start of the EDR to the end of DR4, use
;  two calls as follows.
;
;  set_mjdrange,mjdrange_edr,dr=edr
;  set_mjdrange,mjdrange_dr4,dr=dr4
;  mjdrange=[mjdrange_edr[0],mjdrange_dr4[1]]
;
; MODIFICATION HISTORY:
;       Dec, 2006: Original by Richard McMahon
;
; 
;  $Id: set_mjdrange.pro,v 1.1 2008/07/03 07:52:59 rgm Exp rgm $
;
;  $Log: set_mjdrange.pro,v $
;  Revision 1.1  2008/07/03 07:52:59  rgm
;  Initial revision
;
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2008/07/03 07:52:59 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.1 $ 
;-

splog, traceback(/verbose)
splog, 'period: ', period

mjdrange_all=[50000.0,100000.0] 
mdrange=[0.0,0.0]
; maybe should compute the mjd from the dates

; data/run release date limits
; year, month, day
;Period 87:  Apr 2011 - Sep 2011
p84_date_range=[[2009,10,1],[2010,2,28]]
;print, 'p84_date_range: ', p84_date_range

p85_date_range=[[2010,3,1],[2010,9,30]]
;print, 'p85_date_range: ', p85_date_range

p86_date_range=[[2010,10,1],[2011,3,31]]

p87_date_range=[[2011,4,1],[2011,9,30]]
p88_date_range=[[2011,10,1],[2012,3,31]]

p89_date_range=[[2012,4,1],[2012,9,30]]
p90_date_range=[[2012,10,1],[2013,3,31]]

p91_date_range=[[2013,4,1],[2013,9,31]]
p92_date_range=[[2013,10,1],[2014,3,31]]

p93_date_range=[[2014,4,1],[2014,9,31]]
p94_date_range=[[2014,10,1],[2015,3,31]]

p95_date_range=[[2015,4,1],[2015,9,31]]
p96_date_range=[[2015,10,1],[2016,3,31]]


; trawl through the data releases

; note order needed for JULDAY
; JULDAY(month,day,year,hour,min,sec)

mjdrange_p84=make_array(2, /double)
if strlowcase(period) eq 'p84' or strlowcase(period) eq 'dryrun' $
  or strpos(strlowcase(period), 'dr1') ge 0 then begin
  for i=0,1 do begin
    message,/inf,'P84 Dry Run  date range:'
    mjdrange_p84[i]=$
     JULDAY(p84_date_range[1,i], $
      p84_date_range[2,i], p84_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
     'MJD: ' + STRING(mjdrange_p84[0]) + STRING(mjdrange_p84[1]) 
  message,/inf, 'ISO date: ' + $
     MJD_ISODATE(mjdrange_p84[0]) + '  ' + MJD_ISODATE(mjdrange_p84[1])

endif


mjdrange_p85=make_array(2, /double)
if strlowcase(period) eq 'p85' or $
   strpos(strlowcase(period), 'dr1') ge 0 then begin
   message,/inf,'P85 date range:'

  for i=0,1 do begin
     mjdrange_p85[i]=$
     JULDAY(p85_date_range[1,i], $
      p85_date_range[2,i], p85_date_range[0,i])-2400000.5
  endfor 
 
  message,/inf,$
     'MJD: ' + STRING(mjdrange_p85[0]) + STRING(mjdrange_p85[1]) 
  message,/inf, 'ISO date: ' + $
     MJD_ISODATE(mjdrange_p85[0]) + '  ' + MJD_ISODATE(mjdrange_p85[1])

endif


mjdrange_p86=make_array(2, /double)
if strlowcase(period) eq 'p86' or $
   strpos(strlowcase(period), 'dr2') ge 0 then begin
 message,/inf,'P86 date range:' 

 for i=0,1 do begin
     mjdrange_p86[i]=$
     JULDAY(p86_date_range[1,i], $
      p86_date_range[2,i], p86_date_range[0,i])-2400000.5
 endfor

    message,/inf,$
     'MJD: ' + STRING(mjdrange_p86[0]) + STRING(mjdrange_p86[1]) 
    message,/inf, 'ISO date: ' + $
     MJD_ISODATE(mjdrange_p86[0]) + '  ' + MJD_ISODATE(mjdrange_p86[1])

endif

mjdrange_p87=make_array(2, /double)
if strlowcase(period) eq 'p87' or $
   strpos(strlowcase(period), 'dr2') ge 0 then begin
    message,/inf,'P87 date range:'

  for i=0,1 do begin
    mjdrange_p87[i]=$
     JULDAY(p87_date_range[1,i], $
      p87_date_range[2,i],p87_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p87[0]) + STRING(mjdrange_p87[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p87[0]) + '  ' + MJD_ISODATE(mjdrange_p87[1])

endif

mjdrange_p88=make_array(2, /double)
if strlowcase(period) eq 'p88' then begin
    message,/inf,'P88 date range:'
  for i=0,1 do begin

    mjdrange_p88[i]=$
     JULDAY(p88_date_range[1,i], $
      p88_date_range[2,i],p88_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p88[0]) + STRING(mjdrange_p88[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p88[0]) + '  ' + MJD_ISODATE(mjdrange_p88[1])

endif


mjdrange_p89=make_array(2, /double)
if strlowcase(period) eq 'p89' then begin
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p89[i]=$
     JULDAY(p89_date_range[1,i], $
      p89_date_range[2,i],p89_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p89[0]) + STRING(mjdrange_p89[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p89[0]) + '  ' + MJD_ISODATE(mjdrange_p89[1])

endif


mjdrange_p90=make_array(2, /double)
if strlowcase(period) eq 'p90' then begin
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p90[i]=$
     JULDAY(p90_date_range[1,i], $
      p90_date_range[2,i],p90_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p90[0]) + STRING(mjdrange_p90[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p90[0]) + '  ' + MJD_ISODATE(mjdrange_p90[1])

endif


mjdrange_p91=make_array(2, /double)
if strlowcase(period) eq 'p91' then begin
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p91[i]=$
     JULDAY(p91_date_range[1,i], $
      p91_date_range[2,i],p91_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p91[0]) + STRING(mjdrange_p91[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p91[0]) + '  ' + MJD_ISODATE(mjdrange_p91[1])

endif


mjdrange_p92=make_array(2, /double)
if strlowcase(period) eq 'p92' then begin
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p92[i]=$
     JULDAY(p92_date_range[1,i], $
      p92_date_range[2,i],p92_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p92[0]) + STRING(mjdrange_p92[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p92[0]) + '  ' + MJD_ISODATE(mjdrange_p92[1])

endif

mjdrange_p93=make_array(2, /double)
if strlowcase(period) eq 'p93' then begin
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p93[i]=$
     JULDAY(p93_date_range[1,i], $
      p93_date_range[2,i],p93_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p93[0]) + STRING(mjdrange_p93[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p93[0]) + '  ' + MJD_ISODATE(mjdrange_p93[1])

endif

mjdrange_p94=make_array(2, /double)
if strlowcase(period) eq 'p94' then begin
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p94[i]=$
     JULDAY(p94_date_range[1,i], $
      p94_date_range[2,i],p94_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p94[0]) + STRING(mjdrange_p94[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p94[0]) + '  ' + MJD_ISODATE(mjdrange_p94[1])

endif

mjdrange_p95=make_array(2, /double)
if strlowcase(period) eq 'p95' then begin
  splog, traceback(/verbose)
  for i=0,1 do begin
    message,/inf, period + 'date range:'
    mjdrange_p95[i]=$
     JULDAY(p95_date_range[1,i], $
      p95_date_range[2,i],p95_date_range[0,i])-2400000.5
  endfor

  message,/inf,$
   'MJD: ' + STRING(mjdrange_p95[0]) + STRING(mjdrange_p95[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p95[0]) + '  ' + MJD_ISODATE(mjdrange_p95[1])

endif


mjdrange_vhsdr1=make_array(2, /double)
if strlowcase(period) eq 'vhsdr1' or $
   strpos(strlowcase(period), 'dr1') ge 0 then begin

  message,/inf, period + ' date range:'

  mjdrange_vhsdr1=[mjdrange_p84[0],mjdrange_p85[1]]

  message,/inf,$
   'MJD: ' + STRING(mjdrange_vhsdr1[0]) + STRING(mjdrange_vhsdr1[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_vhsdr1[0]) + '  ' + MJD_ISODATE(mjdrange_vhsdr1[1])

endif

mjdrange_vhsdr2=make_array(2, /double)
if strlowcase(period) eq 'vhsdr2' or $
   strpos(strlowcase(period), 'dr1') ge 0 then begin

  message,/inf, period + ' date range:'

  mjdrange_vhsdr2=[mjdrange_p86[0],mjdrange_p87[1]]

  message,/inf,$
   'MJD: ' + STRING(mjdrange_vhsdr2[0]) + STRING(mjdrange_vhsdr2[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_vhsdr2[0]) + '  ' + MJD_ISODATE(mjdrange_vhsdr2[1])

endif




IF KEYWORD_SET(help) THEN BEGIN

  message,/inf,'ESO Period date ranges'
  message,/inf,'P87 date range; ' + $
   'MJD: ' + STRING(mjdrange_p87[0]) + STRING(mjdrange_p87[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange_p87[0]) + '  ' + MJD_ISODATE(mjdrange_p87[1])

ENDIF

;dqc_set_mjdrange,mjdrange,dr=dr,data=data
; default is to set to full range of dates
mjdrange=mjdrange_all
dr=''
if n_elements(period) gt 0 then begin
  if strlowcase(period) eq 'p84' then mjdrange=mjdrange_p84
  if strlowcase(period) eq 'p85' then mjdrange=mjdrange_p85
  if strlowcase(period) eq 'p86' then mjdrange=mjdrange_p86
  if strlowcase(period) eq 'p87' then mjdrange=mjdrange_p87
  if strlowcase(period) eq 'p88' then mjdrange=mjdrange_p88
  if strlowcase(period) eq 'p89' then mjdrange=mjdrange_p89
  if strlowcase(period) eq 'p90' then mjdrange=mjdrange_p90
  if strlowcase(period) eq 'p91' then mjdrange=mjdrange_p91
  if strlowcase(period) eq 'p92' then mjdrange=mjdrange_p92
  if strlowcase(period) eq 'p93' then mjdrange=mjdrange_p93
  if strlowcase(period) eq 'p94' then mjdrange=mjdrange_p94

  if strlowcase(period) eq 'vhsdr1' then mjdrange=mjdrange_vhsdr1

  if strlowcase(period) eq 'vhsdr2' then mjdrange=mjdrange_vhsdr2

  if strlowcase(period)  eq 'vhsdr1_dr2' or $
     strlowcase(period) eq 'vhsdr1dr2' then $
   mjdrange=[mjdrange_vhsdr1[0],mjdrange_vhsdr2[1]]

  if dr eq 'edr_dr1' then mjdrange=[mjdrange_edr[0],mjdrange_dr1[1]]
  if dr eq 'edr_dr2' then mjdrange=[mjdrange_edr[0],mjdrange_dr2[1]]
  if dr eq 'edr_dr5' then mjdrange=[mjdrange_edr[0],mjdrange_dr5[1]]

endif

message,/inf, period + ' date range:' 
message,/inf,$
   'MJD: ' + STRING(mjdrange[0]) + STRING(mjdrange[1]) 
message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange[0]) + '  ' + MJD_ISODATE(mjdrange[1])

END
