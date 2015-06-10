pro mkpng,filename
;+
;
;
;
;-
print,'Running mkpng: ',!d.name

pngfile=filename

TVLCT,r,g,b,/GET

;img = TVRD()
; use tvread to get around colour problems on mac 
;img= tvread(); retired in 2010; rgm migration 20150610
img= cgsnapshot()
Write_PNG,pngfile,img,r,g,b

;snapshot=TVRD(true=1)
;snapshot=TVRD()
;write_png,filename,snapshot,r,g,b
;write_png,filename,snapshot

print,'png created: ',filename

end


