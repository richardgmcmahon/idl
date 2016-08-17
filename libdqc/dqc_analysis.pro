; $Id: dqc_analysis.pro,v 1.5 2007/03/14 09:17:20 rgm Exp rgm $	
pro dqc_analysis,   phase3=phase3, $
  data=data, chipid=chipid, $
  qcstatus=qcstatus, $
  wsa=wsa, casu=casu, vsa=vsa, $
  vircam=vircam, wfcam=wfcam, $
  survey=survey, sv=sv, dryrun=dryrun, $
  filename=filename, esoperiod=esoperiod, $
  mf=mf, $
  drelease=drelease, zone=zone, $
  gps=gps, yband=yband, kband=kband, $
  ngp=ngp, sgp=sgp, $
  verbose=verbose, debug=debug, $
  stripe82=stripe82,cadence=cadence,msbstats=msbstats, $
  vhs=vhs, viking=viking, $
  offsets=offsets, $
  check_radec=check_radec, $
  report=report, ab=ab,$
  pause=pause, $
  airmass_range=airmass_range, $
  mjd_range=mjd_range, $
  tagnames=tagnames, label=label, $
  linethickness=linethickness, $
  pawprint=pawprint, tile=tile, $
  batch=batch, plotpath=plotpath, filelabel=filelabel, $
  publication=publication, $ ; removes annotation
  polyfill=polyfill, $ ; plots sky coverage using polygons
  dqc_qcstatus=dqc_qcstatus, $
  ob_stats=ob_stats


;+
;
; NAME:
;     DQC_ANALYSIS
;
; PURPOSE:
; Analyze QC data for imaging image supports VISTA, UKIDSS LAS using
; either WSA or CASU DQC. In places it is quite messy due to the
; various input formats supported.
;
;  TODO:
;    make more modular e.g. read dqc in, filter dqc, various plots
;    modules 
;    add inter-band pointing offset diagram.
;
; see also 
;  dqc_check which checks for invalid dqc values
;  dqc_compare which compares casu dqc with wsa dqc
;
; Also plot a histogram in Dec of the TelDecs to identify
; the UKIDSS Dec Stripes. Could split these into RA directories
; maybe for future analysis. 
;
;
; Plot the matching stats for a search radius of 30" and 300" to
; see the main features
;
; need to distinguish the old casu db from the new one.
;
; could split the casu one into dr1, dr2 and dr3
;
; add the dates to the preamble:
;|    | UT Date  range  |       | MJD 
;|    |         |       |
; DR1
; DR2
; DR3
; 
;
;History
;  Original: Aug, Feb 2006 for UKIDSS EDR release
;   Mar, 2007; rgm added interband frame offset analysis   
;   Mar, 2007, added some radec field checks
;   Mar, 2012, publication option
;   Mar, 2012, added polyfill option for the footprint
;   Apr, 2014, fixing so that H band data is not required for VHS
;
;
;
;	$Id: dqc_analysis.pro,v 1.5 2007/03/14 09:17:20 rgm Exp rgm $	
;
;       $Log: dqc_analysis.pro,v $
;       Revision 1.5  2007/03/14 09:17:20  rgm
;       *** empty log message ***
;
;       Revision 1.4  2006/12/26 16:52:39  rgm
;       *** empty log message ***
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2007/03/14 09:17:20 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.5 $ 
;
;-

COMPILE_OPT IDL2

; needed for splog logger
common com_splog, loglun, fullprelog



if not keyword_set(verbose) then $
 message,/inf,'verbose = ' + 'false' 
if keyword_set(verbose) then $
 message,/inf,'verbose = ' + 'true' + string(verbose)


splog, traceback()

splog, get_routine()
get_datestamp, datestamp

routine=get_routine()

splog, get_routine, ' vhs; ',vhs

logfile='dqc_analysis' + $
 cgtimestamp(11, RANDOM_DIGITS=2, /UTC) + '.logfile'
logger, logfile=logfile, /head, hdr=hdr
splog, filename=logfile, /append, systime(/utc)
flush, loglun

htmlfile= 'dqc_analysis_' + label + $
   '_' + datestamp + '.html'
openw, ilun_htmlfile, htmlfile, /get_lun


splog, routine, ': viking; ',viking
splog, routine, ': tile; ',tile
splog, routine, ': pawprint; ',pawprint
splog, routine, ': wsa; ',wsa
splog, routine, ': vsa; ',vsa
splog, routine, ': casu; ',casu
splog, routine, ': des; ',des

!x.thick=1.0
!y.thick=1.0
!p.thick=2.0
;!p.thick=1.0

if not keyword_set(plotpath) then plotpath=""
if keyword_set(plotpath) then begin
  plotpath=plotpath + "/"
  test=file_test(plotpath,/directory)
  if test eq 0 then begin
      message, /inf, 'plotpath: ' + plotpath + ' does not exist'
      message, /inf, 'creating plotpath: ' + plotpath 
      file_mkdir, plotpath
      test=file_test(plotpath,/directory)
     if test eq 0 then message, /inf, 'WARNING: plotpath not created'
     if test ne 0 then message, /inf, 'plotpath created'
  endif
endif

if not keyword_set(verbose) then begin
 message,/inf,'verbose = ' + 'false' 
 verbose=0
endif
if keyword_set(verbose) then $
 message,/inf,'verbose = ' + 'true' + string(verbose)

if not keyword_set(wsa) then wsa=0
if(strpos(filename,'casu')) gt 0 then wsa=0
if(strpos(filename,'fs')) gt 0 then fs=1

message, /inf, traceback()
filename_head=trim_filename(filename,/head)
test=file_test(filename)
if test eq 0 then begin
  message, filename + ' does not exist'
endif
message, /inf, 'filename_head: ' + filename_head
if keyword_set(filelabel) then filelabel=filename_head

if not keyword_set(label) then label= ""
message, /inf, 'label: '+ label
if keyword_set(filelabel) then label = label + '_' + filelabel
message, /inf, 'label: ' + label
if keyword_set(label) then label= label 
message, /inf, 'label: ' + label

; setup the graphics
window,0,xsize=800,ysize=600

if keyword_set(publication) then window,0,xsize=1200,ysize=600

!p.multi=[0,1,1]
!x.style=1
!y.style=1
!p.background=fsc_color('white')
!p.color=fsc_color('black')


if n_elements(survey) gt 0 then begin
  print, 'n_elements(survey): ', n_elements(survey)
  print, 'Survey: ', survey
endif

plotfile_survey=''
if n_elements(survey) lt 1 then survey=''
if n_elements(survey) gt 0 then plotfile_survey=survey

MESSAGE,/INF,'Survey: '+ survey
MESSAGE,/INF,'plotfile_survey: '+ plotfile_survey

IF NOT KEYWORD_SET(ab) THEN ab=0

IF KEYWORD_SET(drelease) THEN dr=drelease

if n_elements(dr) le 0 then begin
  dr='all' 
  print,'using whole mjd range'
endif

if keyword_set(dr) then begin
  set_mjdrange,mjdrange,dr=dr, verbose=verbose
  print,'Using date range: ',mjd_iso(mjdrange[0]),' : ',mjd_iso(mjdrange[1])
endif

if keyword_set(esoperiod) then begin

  upto_period = -1
  from_period = -1

  message,/inf,traceback()
  print, 'n_elements(esoperiod): ', n_elements(esoperiod)
  splog,'ESO period: ', esoperiod

  nperiods=n_elements(esoperiod)
 
  if nperiods eq 1 then begin
    splog, traceback()
    print, 'esoperoid:', esoperiod
    ipos=strpos(strlowcase(esoperiod), 'p')
    print, 'ipos: ', ipos
    if ipos gt 0 then begin
      ipos=strpos(strlowcase(esoperiod), 't')
      if ipos ge 0 then begin
          upto_period=1
          print, 'upto_period:', upto_period
      endif
      ipos=strpos(strlowcase(esoperiod), 'f')
      if ipos ge 0 then begin
          from_period=1
          print, 'from_period:', from_period
      endif
    endif

    ipos=strpos(strlowcase(esoperiod), 'p')
    esoperiod_lower=strmid(esoperiod, ipos)
    message,/inf,'ESO period: ' + esoperiod

    ipos=strpos(strlowcase(esoperiod), 'p')
    print, 'ipos: ', ipos
    if ipos gt 0 then upto_period = 1
    esoperiod_lower=strmid(esoperiod, ipos)
    message,/inf,'ESO period: ' + esoperiod

    splog, traceback()

  endif

  pause, batch=batch

  if nperiods eq 2 then begin
    esoperiod_lower=esoperiod[0]
    esoperiod_upper=esoperiod[1]
    splog, 'esoperiod_lower: ', esoperiod_lower
    splog, 'esoperiod_upper: ', esoperiod_upper
  endif

  pause, batch=batch

  set_mjdrange_vista, mjdrange, period=esoperiod_lower, verbose=verbose
  message,/inf,'ESO Period date ranges: ' + esoperiod 

  if nperiods eq 2 then begin
    set_mjdrange_vista, mjdrange_upper, period=esoperiod_upper, verbose=verbose
    mjdrange[1]=mjdrange_upper[1]
  endif

  pause, batch=batch

  if upto_period gt 0 then begin
    print, 'Get dates upto Period:', esoperiod_upper
    set_mjdrange_vista, mjdrange_test, period='dryrun'
    mjdrange[0]=mjdrange_test[0]
  endif

  if from_period gt 0 then begin
    mjdrange[1]=60000.0
  endif

  if nperiods gt 0 then begin

  message,/inf, $
   'MJD: ' + STRING(mjdrange[0]) + STRING(mjdrange[1]) 
  message,/inf, 'ISO date: ' + $
   MJD_ISODATE(mjdrange[0]) + '  ' + MJD_ISODATE(mjdrange[1])

endif


pause, batch=batch

end

if keyword_set(debug) then pause



; VISTA SV date range
if keyword_set(sv) then begin
   mjdrange=[55118.5, 55138.5]
endif

if keyword_set(dryrun) then begin
   mjdrange=[55138.0,55238.0]
endif

print,'zone, rarange, decrange'
if n_elements(zone) le 0 then begin
  zone='allsky' 
  rarange=[0.0,360.0]
  decrange=[-90.0,90.0]
endif

if n_elements(zone) ge 1 and not keyword_set(vircam) then $
 set_las_zone,rarange,decrange,zone=zone

print,zone,rarange,decrange

If n_elements(filename) le 0 then begin
  if keyword_set(gps) then begin
    ipathname=3
    ifilename=9
  endif
endif


; read in the dqc data within mjd range and ra, dec range
message,/inf,traceback()
t0=systime(1)

filehead=trim_filename(filename,/head)

addtags=1
message, /inf, traceback()
message, /inf, 'verbose: ' + string(verbose)
MESSAGE,/INF,'VHS: '+ string(vhs)
message, /inf, 'Entering dqc_read_data'
dqc_read_data,data,ipathname,ifilename, $
 casu=casu,wsa=wsa,vsa=vsa, vircam=vircam, vhs=vhs, viking=viking, $
 merged=merged, chipid=chipid, $
 filename=filename, $
 mjdrange=mjdrange, $
 rarange=rarange,decrange=decrange, $
 verbose=verbose, debug=debug, $
 survey=survey, pawprint=pawprint, tile=tile, addtags=addtags
message, /inf, traceback()
MESSAGE,/INF,'VHS: '+ string(vhs)
MESSAGE,/INF,'Survey: '+ survey

; add a single mjd column to the fs file mjd, mjdobs_start, mjdobs_end
;dqc_add_mjd, data=data




; analysis by ESO period, obsid

if keyword_set(ob_stats) then begin
  ob_stats_surveys, data=data, verbose=verbose, survey=survey

  obsid=data.obsid
  iobsid=UNIQ(obsid, SORT(obsid))
  itest=where(obsid[iobsid] gt 0, count)
  n_obsid=count
  print, 'Number of unique OB ids: ', n_obsid
endif

splog, traceback()
splog, 'qcstatus: ',qcstatus
if keyword_set(qcstatus) or n_elements(qcstatus) gt 0 then begin

  qcstatus_test=qcstatus
  qcstatus=strtrim(data.qcstatus,2)
  itest=UNIQ(qcstatus, SORT(qcstatus))
  n_qcstatus=n_elements(itest)
  splog, 'Number of unique QC status values: ' + string(n_qcstatus, '(i6)')
  splog, qcstatus[itest]

  ntest=n_elements(qcstatus_test)
  for i=0, ntest-1 do begin
    splog, 'qcstatus_test[i]: ', qcstatus_test[i]
    itest=where(strpos(qcstatus, qcstatus_test[i]) ge 0, count)
    splog, count
    if i eq 0 then newdata=data[itest]
    if i gt 0 then newdata=[newdata,data[itest]]
  endfor

  data=newdata
  ob_stats_surveys, data=data, verbose=verbose, survey=survey

endif

splog, traceback()
splog, 'dqc_qcstatus: ', dqc_qcstatus
if keyword_set(dqc_qcstatus) then dqc_qcstatus, $
 data, verbose=verbose, survey=survey

;print, 'phase3 = ', phase3
if keyword_set(phase3) then print, 'keyword_set(phase3)'
if not keyword_set(phase3) then print, 'not keyword_set(phase3)'

; qc_analysis_exectime
itag=tag_indx(data,'yobsexectime')
if itag ge 0 then begin

exectime=[data.yobsexectime, data.jobsexectime, $
  data.hobsexectime,data.ksobsexectime]
itest=UNIQ(exectime, SORT(exectime))
ntest=n_elements(itest)
for i=0, ntest-1 do begin
  icount=where(exectime eq exectime[itest[i]], count)
  print, 'OB execution time: ',i+1, exectime[itest[i]], count
endfor


xdata=exectime
itest=where(xdata gt 0, count)
xdata=xdata[itest]
xtitle='OB execution times'
print,'X; min, max: ',min(xdata),max(xdata) 
bin=1.0
title=filename
charsize=1.4
plothist,xdata,bin=bin,charsize=charsize,$
   title=title, $
   xtitle=xtitle
pause, batch=batch

; analyze OBS Progam IDs and OB ids
; exclude the null/default vlaues 0f -9999

; concatenate the obsProgIds
obsprogid=[data.yobsProgID, data.jobsprogid, $
 data.hobsprogid, data.ksobsprogid]

obsid=[data.yobsID, data.jobsid, $
 data.hobsid, data.ksobsid]

mjdobs=[data.ymjdobs, data.jmjdobs, $
 data.hmjdobs, data.ksmjdobs]


iobsid=UNIQ(obsid, SORT(obsid))
itest=where(obsid[iobsid] gt 0, count)
n_obsid=count
print, 'Number of unique OB ids: ', n_obsid


itest=UNIQ(obsprogid, SORT(obsprogid))
ntest=n_elements(itest)
for i=0, ntest-1 do begin
  itest2=where(obsprogid eq obsprogid[itest[i]], count)
  itest3=UNIQ(obsid[itest2], SORT(obsid[itest2]))
  ;obsid_unique=
  iobsid=where(obsid[itest3] gt 0, n_obsid)

  print, 'Obs Prog ID: ',i+1, ': ', obsprogid[itest[i]], count, n_obsid
endfor



; cycle through each period and summarise the number of OBs
; and different execution times
; exec times for each observing period
itest=UNIQ(obsprogid, SORT(obsprogid))
obsprogids=obsprogid[itest]
n_obsprogids=n_elements(obsprogids)
; cycle through the runs/periods
for irun=0, n_obsprogids - 1 do begin
  print, 'Obs Prog ID: ', irun+1, ': ', obsprogids[irun]
  itest1=where(obsprogid eq obsprogids[irun], count)

  test_mjdobs=mjdobs[itest1]
  mjd_min=min(test_mjdobs>0) 
  mjd_max=max(test_mjdobs>0)
  print, 'MJD range: ', mjd_min, mjd_max
  print, 'Julian date range: ',mjd_iso(mjd_min), ' : ', mjd_iso(mjd_max)

  exectime_test=exectime[itest1]
  itest2=UNIQ(exectime_test, SORT(exectime_test))
  ntest2=n_elements(itest2)
  for j=0, ntest2-1 do begin
    icount=where(exectime_test eq exectime_test[itest2[j]], count)
    print, 'OB execution time: ',j+1, exectime_test[itest2[j]], count
  endfor
endfor


obsprogids=obsprogid[itest]

exectime=[data.yobsexectime,data.jobsexectime,data.hobsexectime,data.ksobsexectime]
itest=UNIQ(exectime, SORT(exectime))
ntest=n_elements(itest)
for i=0, ntest-1 do begin
  print, 'OB execution time: ',i+1, exectime[itest[i]]
endfor


; analyze SADT pattern
sadtpattern=[data.ysadtpattern,data.jsadtpattern,data.hsadtpattern,data.kssadtpattern]
itest=UNIQ(sadtpattern, SORT(sadtpattern))
ntest=n_elements(itest)
for i=0, ntest-1 do begin
  print, 'SADT Patterns: ',i+1, sadtpattern[itest[i]]
endfor

ndata=n_elements(data)
message,/inf, 'File read in: ' + filename
message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
message,/inf,'DQC file read in: ' + string(ndata) + ' records'
print,'merged: ',merged,'; ','casu:   ',casu,'; ','wsa:    ',wsa

endif

if keyword_set(survey) then begin
  splog, traceback(/verbose)
  message,/inf,'survey= ' + survey
  if strpos(survey,'null') lt 0 and strpos(survey,'ALL') lt 0 then begin
    ob_select_survey, data=data, survey=survey, verbose=verbose, $
     index=index
    data=data[index]
  endif
  ndata=n_elements(data)
  message, /inf, 'DQC file read in: ' + string(ndata) + ' records'
endif
if verbose then pause

message, /inf, traceback()
if verbose then structure_info, data

message, /inf, traceback()
if casu eq 1 then ob_stats_surveys, data=data, verbose=verbose, vst=vst, $
 survey=survey, $
 ntiles_y=ntiles_y, ntiles_j=ntiles_j, $
 ntiles_h=ntiles_h, ntiles_k=ntiles_k, $
 n_ob_unique=n_ob_unique, n_filename_unique=n_filename_unique

message, /inf, traceback()
if verbose then structure_info, data

if wsa lt 0 then wsa=0

if casu eq 1 or keyword_set(casu) then begin
  message,/inf,traceback() + '; casu = ' + string(casu)
  ra=data.cen_ra
  dec=data.cen_dec
endif

if casu ne 1 and not keyword_set(casu) then begin
  message,/inf,traceback() + '; casu = ' + string(casu)
  ra=data.ra
  dec=data.dec
endif


if verbose then structure_info, data

message, /inf, traceback()
print,'RA range: ',min(ra), max(ra)
print,'Dec range: ',min(dec), max(dec)

; filter the dqc 
;dqc_filte
if n_elements(airmass_range) gt 0 then begin
  message,/inf,'Filter airmass range'
  print,'airmass range: ',airmass_range
  itest=where(data.airmass ge airmass_range[0] and  $
   data.airmass le airmass_range[1], count)
  data=data[itest] 
  print,'airmass limit applied: ',n_elements(data),' records'
  airmass=data.airmass
  itest=where(airmass ge airmass_range[0] and  $
   airmass le airmass_range[1], count)
  data=data[itest] 
  print,'airmass limit applied: ',n_elements(data),' records'
endif

if casu gt 0 then begin
  mjd_min=min(data.mjdobs)
  mjd_max=max(data.mjdobs)
endif

if keyword_set(vsa) then begin
  itag=tag_indx(data,'mjdobs')
  if itag ge 0 then begin
    mjd_min=min(data.mjdobs)
    mjd_max=max(data.mjdobs)
  endif
endif

if wsa gt 0 or keyword_set(vsa) then begin
  dqc_get_mjdrange, data=data, mjd_min, mjd_max
endif

message,/inf,'MJD date range: '+STRING(mjd_min)+' : '+STRING(mjd_max)
message,/inf,'Julian date range: '+mjd_iso(mjd_min)+' : '+mjd_iso(mjd_max)
message, /inf, traceback()

if mjd_max lt 0 then mjdrange=[-999.0,-9.9]

if keyword_set(verbose) and not keyword_set(batch) then help,data,/str

if keyword_set(check_radec) then dqc_check_radec,data

if keyword_set(offsets) then begin
; analyze wsa framesets
  if(strpos(filename,'fs')) gt 0 then  dqc_fs_offsets,data,filename=filename
; analyze ungrouped casu dqc file 
  if(strpos(filename,'fs')) lt 0 then begin
    dqc_offsets,data,filename=filename
  endif
endif

if not keyword_set(wsa) AND not keyword_set(vsa) then begin

  message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
  itest=UNIQ(data.filtname, SORT(data.filtname))
  nfilters=n_elements(itest)
  for i=0, nfilters-1 do begin 
    print, i, data[itest].filtname
  endfor
  message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
  message, /inf,'Number of unique filters: ' + string(nfilters)

  filtername=data.filtname
  message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
  ifilter=UNIQ(filtername, SORT(filtername))
  nfilters=n_elements(ifilter)
  message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
  message, /inf,'Number of unique filters: ' + string(nfilters)

endif

; need to make this smarter e.g make a list of unique filternames
wavebands=['Y','J','H','K']
if wsa ge 1 then wavebands=['Y','J_1','H','K']

if keyword_set(viking) then wavebands=['Z','Y','J','H','KS']
if keyword_set(sv) then wavebands=['Z','Y','J','H','K']
if keyword_set(vhs) then wavebands=['Y','J','H','KS']

if strpos(strlowcase(survey),'des') ge 0 then wavebands=['J', 'H', 'KS']
if strpos(strlowcase(survey),'gps') ge 0 then wavebands=['J', 'KS']
if strpos(strlowcase(survey),'atlas') ge 0 then wavebands=['Y','J','H','KS']

message, /inf, traceback()
print,'wavebands: ',string(wavebands)
splog, 'wavebands: ', wavebands


if keyword_set(vhs) then begin
  ; generate wavebands from the data since H is not always present
  if not keyword_set(phase3) then begin
    wavebands=strtrim(filtername[ifilter[0]],2)
    ; append the filter names
    for i=1, nfilters-1 do begin
      wavebands=[wavebands, strtrim(filtername[ifilter[i]],2)] 
    endfor
  endif  

  ;
  ;wavebands=['Y', 'J', 'H', 'KS']

endif

wavebands=strtrim(wavebands,2)

splog, 'wavebands: ', wavebands
print,'wavebands: ',string(wavebands)
nbands=n_elements(wavebands)


; make array of unique filters/wavebands from the input file
itag=tag_indx(data,'filtname')
if itag le 0 then itag=tag_indx(data,'filtername')
itest=UNIQ(data.(itag), SORT(data.(itag)))
wavebands=data[itest].(itag)
wavebands=strtrim(wavebands,2)

print, 'vhs: ', vhs
if keyword_set(vhs) then begin
  waveband_ordered=['Y','J','H','Ks']
  nmax=n_elements(waveband_ordered)
  nbands=n_elements(wavebands)
  order=make_array(nbands,/int)
  for i=0, nbands-1 do begin
    for j=0, nmax-1 do begin
      if strupcase(wavebands[i]) eq strupcase(waveband_ordered[j]) then order[i]=j
    endfor
  endfor
  print, wavebands
  print, 'filter order: ', order
  isort=sort(order)
  wavebands=wavebands[isort]
  print, wavebands 
endif


splog, 'wavebands: ', wavebands
print,'wavebands: ',string(wavebands)
nbands=n_elements(wavebands)

wavebands=strtrim(wavebands,2)

;dqc_check,data,casu=casu,wsa=wsa,merged=merged

; check the ra, dec coverage
;check_ra,data.cen_ra,data.cen_dec,/plot

if casu eq 1 or keyword_set(casu) then begin
  ra=data.cen_ra
  dec=data.cen_dec
endif

if keyword_set(vsa) then begin
    itag=tag_indx(data,'ra')
    if itag ge 0 then begin
      data.ra=data.ra*!radeg
     ra=data.ra
    endif
    itag=tag_indx(data,'rabase')
    if itag ge 0 then ra=data.rabase

    itag=tag_indx(data,'dec')
    if itag ge 0 then begin
      data.dec=data.dec*!radeg
      dec=data.dec
    endif
    itag=tag_indx(data,'decbase')
    if itag ge 0 then dec=data.decbase
endif

if wsa eq 1 then begin
  if merged eq 1 then begin
    ra=data.ra*!radeg
    dec=data.dec*!radeg
  endif
  if merged ne 1 then begin
    ra=data.centralra
    dec=data.centraldec
    print,'n_elements(ra): ',n_elements(ra)
  endif
endif


MESSAGE,/inf,'RA range: ' + $
 STRING(min(ra))+' : '+STRING(max(ra))

MESSAGE,/inf,'Dec range: ' + $
 STRING(min(dec))+' : '+STRING(max(dec))

if keyword_set(cadence) and merged eq 1 then begin
  dqc_cadence,data
endif

if keyword_set(msbstats) and merged eq 1 then begin
  dqc_msbstats,data
endif

;!p.multi=[0,2,2]
!p.multi=[0,1,1]
title=trim_filename(filename)

print,'n_elements(data): ',n_elements(data)
print,'mjdrange: ',mjdrange


; do some ra, dec checks
; an option to turn these off could be added
print,'RA range: ',min(ra),max(ra)
print,'Dec range: ',min(dec),max(dec)

;ra=ra[where(dec gt -50.0 and dec lt 70.0,count)]
;dec=dec[where(dec gt -50.0 and dec lt 70.0,count)]

if keyword_set(yband) then begin
  if casu eq 1 then begin
    print,'Limiting by mjdrange and to Y band'
    index=where(strpos(data.filtname,'Y') ge 0 and $
     data.mjdobs ge mjdrange[0] and $
     data.mjdobs le mjdrange[1],count)
  endif
  if wsa eq 1 then begin
    index=where(strmid(data.filtername[0:0],0,1) eq 'Y',count)
  endif
  print,' Y band count: ',count,count/21.0
  ra=ra[index]
  dec=dec[index]
  mjd=data[index].mjdobs
  ntest=where(ra ge 100.0 and ra le 300.0,count)
  print,' Y band count[NGP]: ',count,count/21.0
  data=data[index]
endif

IF NOT keyword_set(yband) THEN BEGIN
  IF keyword_set(casu) THEN mjd=data.mjdobs
  IF keyword_set(wsa) THEN mjd=data.ymjdobs
ENDIF

IF keyword_set(casu) THEN mjd=data.mjdobs
IF keyword_set(wsa) THEN mjd=data.ymjdobs
IF keyword_set(vsa) THEN begin
  mjd=-1
  itag=tag_indx(data,'mjdobs')
  if itag ge 0 then mjd=data.(itag)

  itag=tag_indx(data,'y'+'mjdobs')
  if itag ge 0 then mjd=data.(itag)
  if itag ge 0 then begin
    mjd_min=min(data.(itag)>0) 
    mjd_min=max(data.(itag)>0) 
  endif

 itag=tag_indx(data,'j'+'mjdobs')
  if itag ge 0 then mjd=data.(itag)
  if itag ge 0 then begin
    mjd_min=min([mjd_min,data.(itag)>0]) 
    mjd_min=max([mjd_max,data.(itag)>0])
  endif


endif

itest=where(mjd gt 1,count)
mjd_min=min(mjd[itest])
mjd_max=max(mjd[itest])

print,'MJD range of data: ',mjd_min, mjd_max
print,'Date range of data: ',mjd_iso(mjd_min),' : ',mjd_iso(mjd_max)


print,'RA range: ',min(ra),max(ra),n_elements(ra)
print,'Dec range: ',min(dec),max(dec),n_elements(dec)
title=trim_filename(filename)
if keyword_set(esoperiod) then title = esoperiod + ': ' + title
if keyword_set(qcstatus) then title = $
 'QCstatus:' + strjoin(qcstatus_test) + ' ' + title

xrange=[0.0,24.0]

if keyword_set(publication) then yrange=[-10.0,45.0]
if keyword_set(vircam) then message,/inf,'vircam'
if keyword_set(vsa) then message,/inf,'vsa'
if keyword_set(ngp) then message,/inf,'ngp'
if keyword_set(vsa) or keyword_set(vircam) then yrange=[-90.0,20.0] 
if keyword_set(ngp) then begin
  xrange=[8.0,16.0]
  yrange=[-10.0,20.0]
endif


;dqc_skyplot, data, title=title, xrange=xrange, yrange=yrange
;pro dqc_skyplot, data, title=title, xrange=xrange, yrange=yrange

xtitle=TexToIDL('Right Ascension (hours)')
ytitle=TexToIDL('Declination (degrees)')

ndata=n_elements(data)

xdata=ra/15.0
ydata=dec

itest=where(ra/15.0 ge xrange[0] and ra/15.0 le xrange[1] and $ 
            dec ge yrange[0] and dec le yrange[1],count)

message, /inf,'ndata= ' + string(ndata)
message, /inf,'count= ' + string(count)

; reverse the RA axis
;xrange=[24.0,0.0]

if keyword_set(sgp) then message,/inf,'sgp'
message,/inf,'zone: ' + zone
if zone eq 'sgc' or keyword_set(sgp) then begin
  for i=0, ndata-1 do begin
   ;print, i, xdata[i]
   if xdata[i] gt 18.0 then xdata[i]=xdata[i]-24.0
  endfor

  xrange=[-3.0,5.0]
  yrange=[-5.0,25.0]
  itest=where(xdata ge xrange[0] and xdata le xrange[1] and $ 
            ydata ge yrange[0] and ydata le yrange[1],count)
endif 

if keyword_set(publication) then title=''
print, 'Plot RA range: ',min(xdata), max(xdata), n_elements(xdata)
print, 'Plot Dec range: ',min(ydata), max(ydata), n_elements(ydata)
plotsym,8,/fill

if keyword_set(polyfill) then nodata=1
if keyword_set(publication) then begin
  xthick=2.0
  ythick=2.0
  thick=2.0
  charthick=2.0
  charsize=2.0
  !p.font=1
  ;DEVICE, SET_FONT='Helvetica', /TT_FONT, SET_CHARACTER_SIZE=[12,12]
  DEVICE, SET_FONT='Courier', /TT_FONT, SET_CHARACTER_SIZE=[9,12]
  ;DEVICE, SET_FONT='Times', /TT_FONT, SET_CHARACTER_SIZE=[12,15]
endif

;title='!5Title of my plot'
;title='Title of my plot'



charsize=1.4
plot,xdata,ydata, psym=8,symsize=0.70,charsize=charsize,$
 title=title, xtitle=xtitle,ytitle=ytitle,$
 xrange=xrange,yrange=yrange,nodata=nodata, $
 thick=thick, xthick=xthick,ythick=ythick,charthick=charthick
if not keyword_set(publication) then plotid, /right
message,/inf,traceback()

if keyword_set(polyfill) then begin
   fieldsize_width=(2048*0.40)/3600.00
  plot_radec_polyfill, fieldsize_width=fieldsize_width, $
   xdata=xdata, ydata=ydata
  message,/inf,traceback()
endif


; plot galactic plane
if not keyword_set(publication) then  begin
plot_galcoords, /oplot, b=0, color=fsc_color('black'), linestyle=0
if keyword_set(vircam) then begin
  plot_galcoords, /oplot, b=-5, color=fsc_color('red'), linestyle=2
  plot_galcoords, /oplot, b=5, color=fsc_color('red'), linestyle=2
endif
plot_galcoords, /oplot, b=30, color=fsc_color('red'), linestyle=3
plot_galcoords, /oplot, b=-30, color=fsc_color('red'), linestyle=3
;plot_galcoords, /oplot, b=20, color=fsc_color('red'), linestyle=3
;plot_galcoords, /oplot, b=-20, color=fsc_color('red'), linestyle=3

legend='npoints= '+string(count,format='(i8)')
if keyword_set(vircam) then begin
  legend='npoints= '+string(count,format='(i8)')
endif
endif

if not keyword_set(publication) then $
al_legend, legend, /right,/bottom,charsize=1.4

string = 'Date range: '+string(mjd_iso(mjd_min))+$
  ' : '+string(mjd_iso(mjd_max))
string2= 'MJD range:  ' + string(mjd_min) + ' : '  + string(mjd_max)
string = [string, string2]
if not keyword_set(publication) then $
  al_legend, string, /top, /left, charsize=1.4


if not keyword_set(publication) then  begin
  if keyword_set(zone) then al_legend, zone, /left, /bottom, charsize=1.4
endif

message,/inf,traceback()
message, /inf, 'label: ' + label
plotfile = plotpath + 'dqc_analysis_skyplot_' + $
 label + '_' + filehead + '_' + plotfile_survey + datestamp + '.png'

pause,plotfile=plotfile,batch=batch

if keyword_set(vsa) then begin
  xdata=-1
  itag=tag_indx(data,'stdrms')
  if itag ge 0 then xdata=data.stdcrms  
  itag=tag_indx(data,'numrms')
  ydata=-1
  if itag ge 0 then ydata=data.numrms
endif

if not keyword_set(vsa) and not keyword_set(wsa) then begin
  xdata=data.wcsrms
  ydata=data.numwcsfit
endif


message,/inf,traceback() + '; wsa = ' + string(wsa)
IF not keyword_set(wsa) and not keyword_set(vsa) THEN begin

print,'WCS rms(arcsecs) range: ',min(xdata), max(xdata)
print,'WCS nfit range: ',min(ydata), max(ydata)

; setup the graphics
window,0,xsize=600,ysize=600

title=trim_filename(filename)
if keyword_set(esoperiod) then title = esoperiod + ': ' + title
xtitle='WCS rms(arcsec)'
ytitle='WCS calibration stars'
charsize=1.4
yrange=[5,2000]
if keyword_set(tile) or keyword_set(phase3) then  yrange=[500,10000]
plot, xdata, ydata, /xlog, /ylog, psym=1, $
   xrange=[0.05,1.0],yrange=yrange, $
   charsize=charsize, $
   title=title,xtitle=xtitle,ytitle=ytitle
ndata=n_elements(xdata)

al_legend,string,/top,/right,charsize=1.4
legend='n= ' + string(ndata,'(i7)')
al_legend,legend,/bottom,/right,charsize=1.4
plotid, /right

message,/inf,traceback()
message, /inf, 'label: ' + label

plotfile=plotpath+'dqc_analysis_wcs_' + $
 label + '_' + filehead + '_' + datestamp + '.png'
pause, plotfile=plotfile,batch=batch


endif

; now band by band analysis
message, /inf, 'band by band analysis'
message, /inf, traceback()
splog, 'wavebands: ', wavebands
;itest=UNIQ(filtername, SORT(filtername))

nx=2
ny=2
window,0,xsize=800,ysize=800
if nbands gt 4 then ny=3
!p.multi=[0,nx,ny]

; determine the names of the unique filters
if not keyword_set(fs) then begin
  if keyword_set(casu)  then itag=tag_indx(data,'filtname')
  if keyword_set(wsa) then itag=tag_indx(data,'filtername')
  itest=UNIQ(data.(itag), SORT(data.(itag)))
  ntest=n_elements(itest)
  for i = 0, ntest-1 do begin
    print, i, data[itest[i]].(itag)
  endfor
endif

for iband=0,nbands-1 do begin

  charsize=2.0
 
  if not keyword_set(fs) then begin
  
    index = $
     where(strpos(strupcase(data.(itag)),strupcase(wavebands[iband])) $
     ge 0,count)

    print,'count= ',count,' :',wavebands[iband],casu,wsa

    if not keyword_set(vsa) then begin
      xdata=data[index].wcsrms
      ydata=data[index].numwcsfit
    endif

    if keyword_set(vsa) then begin
      xdata=data[index].stdcrms
      ydata=data[index].numrms
    endif

  endif

  ; vsa/wsa multiband framesets
  if keyword_set(fs) then begin
    tag=wavebands[iband] + 'MULTIFRAMEID'
    itag=tag_indx(data,tag)
    itest=where(data.(itag) gt 0,count)
    message,/inf,'Number of valid: ' + tag + ' ' + string(count)  
    xtag=wavebands[iband] + 'stdcrms'
    message,/inf,'xtag: ' + xtag 
    ytag=wavebands[iband] + 'numrms'
    message,/inf,'ytag: ' + ytag 
    ixtag=tag_indx(data,xtag)
    message,/inf,'ixtag: ' + string(ixtag)
    iytag=tag_indx(data,ytag)
    message,/inf,'iytag: ' + string(iytag)
    xdata=data[itest].(ixtag)
    ydata=data[itest].(iytag)
  endif

  title=trim_filename(filename)
  if keyword_set(esoperiod) then title = esoperiod + ': ' + title
  if (iband gt 0) then title=''
  xtitle='WCS rms(arcsec)'
  ytitle='Number of WCS calibration stars'
  charsize=sqrt(ny)
  xrange=[0.05,1.0]

  yrange=[5,1000]
  if keyword_set(tile) or keyword_set(phase3) then  yrange=[500,10000]
  psym=3

  ndata=n_elements(xdata)
  print, 'ndata: ', ndata
  if ndata gt 1 then begin

    plot, xdata, ydata, /xlog, /ylog, psym=psym, $
     xrange=xrange,yrange=yrange, $
     charsize=charsize, $
     title=title,xtitle=xtitle,ytitle=ytitle

    print, 'X range: ', minmax(xdata)
    print, 'Y range: ', minmax(ydata)
    legend='n = ' + string(ndata,'(i6)')
    al_legend, legend, /bottom, /right, charsize=1.4

    al_legend, wavebands[iband], /right

  endif
  
endfor

plotfile=plotpath+'dqc_analysis_wcs_wavebands_' + $
 label + '_' + filehead + '_' + datestamp + '.png'
plotid, /right
message,/inf,traceback()
message, /inf, 'label: ' + label
pause,plotfile=plotfile,/force,batch=batch

; initialise a list of columns for the dqc plots

if not keyword_set(vsa) then begin
  tagnames=['maglimit','seeing','ellip','wcsrms', 'zeropoint', 'skybrightness']
endif

if keyword_set(vsa) then begin
  tagnames=['maglimit','stdcrms','seeingpix','seeingarcsecs',$
   'ellip','tablerows', $
   'photzpcat','apercor3']
endif

if keyword_set(vircam) and keyword_set(casu) then begin
  tagnames=['maglimit', 'seeing', 'stdcrms', 'numwcsfit', $
   'magzpt','ellipticity', $
   'airmass', 'skybrightness', 'skylevel', 'exptime']
endif

ntagnames=n_elements(tagnames)

!p.multi=[0,1,1]
title=trim_filename(filename)
if keyword_set(esoperiod) then $
 title = esoperiod + ' ' + plotfile_survey + ': ' + title
print,'wavebands: ', wavebands

pause, batch=batch

for itagname=0, ntagnames-1 do begin

  tagname=tagnames[itagname]
  splog, itagname,': ', tagname
  dqc_plot_cumulative, data=data, tagname=tagname, $
   mjdrange=mjdrange, fs=fs, tile=tile, $
   casu=casu, wsa=wsa, vsa=vsa, vircam=vircam, $
   ab=ab, zone=zone, viking=viking, $
   title=title, wavebands=wavebands, verbose=verbose

  message,/inf
  plotfile=plotpath+'dqc_analysis_' + tagname + '_' + $ 
   label + '_' + filehead + '_' + datestamp + '.png'
  if keyword_set(ab) and $
   (tagname eq 'maglimit' or tagname eq 'zeropoint' or $
   tagname eq 'skybrightness' or tagname eq 'zp' or $
   tagname eq 'magzpt' or tagname eq 'zp' or $
   tagname eq 'photzpcat' ) then begin
    plotfile=plotpath+'dqc_analysis_' + tagname + '_' + $
     label + '_' + filehead + '_' + datestamp + '_ab.png'
  endif
  message, /inf, traceback()
  message, /inf, 'label: ' + label
  pause, plotfile=plotfile,batch=batch

endfor


!p.multi=[0,1,1]


skip:

; plot some dqc
if casu eq 1 then dqc_multiplot, data,/casu,mjdrange=mjdrange,$
 filename=filename, wavebands=wavebands, survey=survey, vircam=vircam, $
 viking=viking, tile=tile, $
 batch=batch, plotpath=plotpath, label=label

if wsa eq 1 then dqc_multiplot, data,/wsa,/gps, batch=batch, $
 viking=viking, tile=tile

message, /inf

; check the date range
itag=tag_indx(data,'mjdobs')
if itag ge 0 then mjd=data.(itag)
itag=tag_indx(data,'j_1mjdobs')
if itag ge 0 then mjd=data.(itag)
itag=tag_indx(data,'jmjdobs')
if itag ge 0 then mjd=data.(itag)

itest=where(mjd lt 40000.0,count)
print,'Number of bad MJD values: ',count
if count gt 0 then begin
  mjd_junk=mjd[itest]
  ;print,mjd_junk
endif

; mjd plots

;dqc_analysis_mjd
!p.multi=[0,1,2]

title=filename
xdata=mjd

if keyword_set(vsa) then begin
  itag=tag_indx(data,'y'+'stdcrms')
  ydata=data.(itag)
  index=where(ydata gt 0, count)
  ydata=ydata[index]
endif

if not keyword_set(vsa) then begin
  wcsrms_flag=where(data.wcsrms gt 5.0,wcsrms_nflag)
  print,'wcsrms_nflag= ',wcsrms_nflag

  index=where(mjd gt 40000.0 and data.wcsrms lt 0.5 and data.wcsrms gt 0.0 and $
      mjd gt mjdrange[0] and mjd lt mjdrange[1])
  xdata=xdata[index]
  wcsrms=data.wcsrms
  ydata=wcsrms[index]
  splog, traceback()
endif

xtitle='mjd'
ytitle='wcsrms'
print,'yrange: ',min(ydata),max(ydata)
plot,xdata,ydata,psym=3.0,$
 title=title,xtitle=xtitle,ytitle=ytitle,charsize=1.4

string='Date range of data: '+string(mjd_iso(min(mjd[index])))+$
  ' : '+string(mjd_iso(max(mjd[index])))
al_legend,string,/top,/left

message, /inf, traceback()

if n_elements(index) gt 1 then xdata=mjd[index]

splog, 'minmax(xdata): ',minmax(xdata)

;if keyword_set(vsa) then begin

if keyword_set(vsa) or not keyword_set(vsa) then begin

  message, /inf, traceback()

  if not keyword_set(vsa) then numwcsfit=data.numwcsfit
  if keyword_set(vsa) then numwcsfit=data.numrms

  ydata=numwcsfit[index]

  splog, 'minmax(ydata): ',minmax(ydata)

  print,'yrange: ',min(ydata),max(ydata)
  ytitle='numwcsfit'
  plot,xdata,ydata,psym=3.0,$
   title=title,xtitle=xtitle,ytitle=ytitle,charsize=1.4

  message, /inf, traceback()
  message, /inf, 'label: ' + label
  plotfile=plotpath+'dqc_analysis_mjd_wcs_' + $
     label + '_' + filehead + '_' + datestamp + '.png'
  pause, plotfile=plotfile, batch=batch

  message,/inf, traceback()

  dqc_plot_mjd, data, casu=casu, wsa=wsa, $
   filename=filename, mjdrange=mjdrange, label=label, $
   wavebands=wavebands, survey=survey, batch=batch, plotpath=plotpath

  message,/inf, traceback()
  message, /inf, 'label: ' + label
  ;analyze the zeropoints
  !p.multi=[0,1,nbands]
  charsize=1.4

  for iband=0,nbands-1 do begin

  if casu gt 0 then $
   index=where(strpos(strupcase(data.filtname), strupcase(wavebands[iband])) ge 0,count)
  if wsa gt 0 then $
    index=where(strmid(strupcase(data.filtername[0:0]),0,1) eq strupcase(wavebands[iband]),count)
  print,'count= ',count,' :',wavebands[iband],casu,wsa

  xdata=data[index].mjdobs
  if wsa gt 0 then ydata=data[index].photzpcat
  if casu gt 0 then ydata=data[index].magzpt

  xtitle='mjdobs'
  ytitle='photozpcat'

  charsize=2.0
  xthick=1.0
  ythick=1.0
  thick=2.0

  psym=3
  waveband=wavebands[iband]
  plot,xdata,ydata,psym=psym,charsize=charsize,$
;  xrange=xrange,yrange=yrange,$
  xtitle=xtitle,ytitle=ytitle,$
  thick=thick,xthick=xthick,ythick=ythick
  al_legend, waveband


  endfor

  message, /inf, traceback()

  message, /inf, 'label: ' + label
  plotfile=plotpath+'dqc_analysis_mjd_photozpcat_' + $
     label + '_' + filehead + '_' + datestamp + '.png'
  pause, batch=batch, plotfile=plotfile
  message,/inf,traceback()

  ;analyze the xpixsize verus mjd
  for iband=0,nbands-1 do begin

  if casu gt 0 then $
   index=where(strpos(strupcase(data.filtname), strupcase(wavebands[iband])) ge 0,count)
  if wsa gt 0 then $
    index=where(strmid(strupcase(data.filtername[0:0]),0,1) eq strupcase(wavebands[iband]),count)
  print,'count= ',count,' :',wavebands[iband]

  xdata=data[index].mjdobs
;  xdata=juliandaynum

  if wsa gt 0 then ydata=data[index].xpixsize
  if casu gt 0 then ydata=abs(data[index].cd11*(3600.0))
  xtitle='mjdobs'
  ytitle='xpixsize'
  yrange=[0.15,0.45]  
  yrange=[0.340,0.342]  
  waveband=wavebands[iband]

  print, 'minmax(y): ',minmax(ydata)
  charsize=3.0
  plot,xdata,ydata,psym=psym,charsize=charsize,$
;  xrange=xrange,
   yrange=yrange,$
   xtitle=xtitle,ytitle=ytitle  
  al_legend,waveband

  endfor
  plotfile=plotpath+'dqc_analysis_mjd_xpixsize_' + $
     label + '_' + filehead + '_' + datestamp + '.png'
  pause, batch=batch, plotfile=plotfile
  message,/inf,traceback()



; various different zeropoints
; photzpext
; photzpcat
; 


if wsa gt 0 then begin

window,1,xsize=1000,ysize=750
!p.multi=[0,5,nbands]
charsize=1.7

for iband=0,nbands-1 do begin

  if casu gt 0 then $
   index=where(strpos(strupcase(data.filtname), strupcase(wavebands[iband])) ge 0,count)
  if wsa gt 0 then $
    index=where(strmid(strupcase(data.filtername[0:0]),0,1) eq strupcase(wavebands[iband]),count)
  print,'count= ',count,' :',wavebands[iband]

  avstellarell=data[index].avstellarell
  seeing=data[index].seeing*data[index].xPixSize
  skynoise=data[index].skynoise
  skylevel=data[index].skylevel
  photzpcat=data[index].photzpcat
  xpixsize=data[index].xpixsize
  exptime=data[index].exptime
  apercor3=data[index].apercor3

; 5sigma limit magnitude in aper3 aperture which has a diameter of 2arcsecs
  maglimit=photzpcat-2.5*alog10(5.0*skynoise*sqrt(1.2*3.141593)/(xpixsize*exptime))-apercor3

; sky brightness in magnitudes per square arc second from skylevel per pixel
  skybrightness=photzpcat-2.5*alog10(skylevel/(xpixsize*xpixsize*exptime))
;g
  if iband eq 0 then begin 
    yrange=[19.4,21.1]
   endif
  if iband eq 1 then begin 
    yrange=[18.6,20.3]
   endif
  if iband eq 2 then begin 
    yrange=[17.9,19.6]
   endif
  if iband eq 3 then begin 
    yrange=[17.2,18.9]
   endif

;
  xdata=maglimit
  xtitle='maglimit'
  print,'X; min, max: ',min(xdata),max(xdata) 
  bin=0.10
  plothist,xdata,bin=bin,charsize=charsize,$
   xtitle=xtitle

  xdata=seeing
  xtitle='seeing(")'
  print,'X; min, max: ',min(xdata),max(xdata) 
  bin=0.05
  plothist,xdata,bin=bin,charsize=charsize,$
   xtitle=xtitle

;
  xdata=skybrightness
  xtitle='skybrightness'
  print,'X; min, max: ',min(xdata),max(xdata) 
  bin=0.10
  plothist,xdata,bin=bin,charsize=charsize,$
   xtitle=xtitle

;
  xdata=avstellarell
  xtitle='avstellarell'
  print,'X; min, max: ',min(xdata),max(xdata) 
  bin=0.01
  plothist,xdata,bin=bin,charsize=charsize,$
   xtitle=xtitle

; 
  xdata=photzpcat
  xtitle='photzpcat'
  print,'X; min, max: ',min(xdata),max(xdata) 
  bin=0.025
  plothist,xdata,bin=bin,charsize=charsize,$
   xtitle=xtitle

endfor

plotfile='test1.png'
pause, batch=batch, plotfile=plotfile

endif

message, /inf, traceback()

if wsa gt 0 then begin
!p.multi=[0,2,2]
charsize=1.4
for iband=0,nbands-1 do begin

  if casu gt 0 then $
   index=where(strpos(strupcase(data.filtname), strupcase(wavebands[iband])) ge 0,count)
  if wsa gt 0 then $
    index=where(strmid(strupcase(data.filtername[0:0]),0,1) eq strupcase(wavebands[iband]),count)
  print,'count= ',count,' :',wavebands[iband]

  photzpcat=data[index].photzpcat
  photzpext=data[index].photzpext

  xdata=photzpcat
  ydata=photzpext
  
  xtitle='photzpcat'
  ytitle='photzpext'

  plot,xdata,ydata,psym=psym,charsize=charsize,$
;  xrange=xrange,yrange=yrange,$
  xtitle=xtitle,ytitle=ytitle  

endfor


endif

message, /inf, traceback()

if wsa gt 0 then begin
!p.multi=[0,1,1]
nrows=n_elements(data)
print,'nrows= ',nrows
; check that the filtername is valid
nfiltername_problem=0
irow=0
print,'filtername:',irow,data[irow].filtername[0:0]
for irow=0L,nrows-1 do begin
  if (strmid(data[irow].filtername[0:0],0,1) ne 'Y' and $
   strmid(data[irow].filtername[0:0],0,1) ne 'J' and $
   strmid(data[irow].filtername[0:0],0,1) ne 'H' and $
   strmid(data[irow].filtername[0:0],0,1) ne 'K') then begin
   print,'filtername:',strmid(data[irow].filtername[0:0],0,1)
    nfiltername_problem=nfiltername_problem+1
  endif
endfor
print,'nfiltername_problem= ',nfiltername_problem
endif

message, /inf, traceback()

;message,/inf

;goto, exit

;end

if casu eq 1 then dqc_analysis_plot, data, /casu, $
 wavebands=wavebands, batch=batch

if wsa eq 1 then dqc_analysis_plot, data, /wsa, $
 wavebands=wavebands, batch=batch

exit:

close,ilun_htmlfile

message, /inf, traceback()

end

message, /inf, traceback()
message, /inf, 'DQC filename completed: ' + filename

end
