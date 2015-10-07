
pro plot_despolygon, ra_despolygon, dec_despolygon, $
 des_polyfill=des_polyfill

  print, 'n_elements(ra_despolygon):  ', n_elements(ra_despolygon)
  print, 'n_elements(dec_despolygon): ', n_elements(dec_despolygon)

  if keyword_set(des_polyfill) then begin
    polyfill, ra_despolygon, dec_despolygon, col=fsc_color('blue'), noclip=0
  endif

  if not keyword_set(des_polyfill) then begin
    oplot, ra_despolygon, dec_despolygon, col=fsc_color('blue'), thick=2
  endif

  ndata=n_elements(ra_despolygon)
  xdata=ra_despolygon*0.0
  for i=0, ndata-1 do begin
    xdata[i]=24.0+ra_despolygon[i]
  endfor

  if keyword_set(des_polyfill) then begin
    polyfill, xdata, dec_despolygon, col=fsc_color('blue'), noclip=0
  endif

  if not keyword_set(des_polyfill) then begin
    oplot, xdata, dec_despolygon, col=fsc_color('blue'), thick=2
  endif

  ; charthick=2 is horrible
  xyouts, 7.0, -85.0, 'DES: Round13-poly (Aug 20, 2013)', $
   charsize=1.4, color=fsc_color('Blue'), charthick=1

end
