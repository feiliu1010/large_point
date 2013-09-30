pro gaussian_fit_regular
flag=0U
nx=200
ny=200
curLon=112.5738889
curLat =35.4655555
distance = 24000

For year= 2005, 2010 do begin
	For month=1,12 do begin
	        nymd = year*10000L+month*100L+1*1L
	        Yr4 = string(year,format='(i4.4)')
	        mon2 = string(month, format='(i2.2)')
		file1_in = string('/home/liufei/r5/large_point/OMI_NO2_yc_0.0167deg/OMI_NO2_yc_0.0167deg/OMI_NO2_'+Yr4+Mon2+'.nc')
		file1_in=strcompress(file1_in, /REMOVE_ALL)
		fid=NCDF_OPEN(file1_in)
		varid=NCDF_VARID(fid,'TropVCD')
		NCDF_VARGET, fid, varid, VCD_temp
                varid=NCDF_VARID(fid,'LONGITUDE')
                NCDF_VARGET, fid, varid, X
                varid=NCDF_VARID(fid,'LATITUDE')
                NCDF_VARGET, fid, varid, Y
		NCDF_CLOSE, fid

		if flag then begin
                	VCD=[[VCD],[VCD_temp]]
                endif else begin
                        VCD = VCD_temp
                        flag =1U
                endelse
;		For I=0,nx-1 do begin
;			For J=0,ny-1 do begin
;		        	if  MAP_2POINTS(curLon,curLat,X[I],Y[J],/METERS) le distance then begin		
;					if flag then begin
;						VCD=[[VCD],[VCD_temp]]
;					endif else begin
;						VCD = VCD_temp
;						flag =1U
;					endelse
;				endif else begin
;                                       if flag then begin
;                                              VCD=[[VCD],[replicate(0,nx,ny)]]
;                                        endif else begin
;                                                VCD = replicate(0,nx,ny)
;                                                flag =1U
;                                        endelse
;				endelse
;			endfor
;		endfor
	endfor
endfor
print,size(VCD)
print,size(x)
print,size(y)
s=size(VCD)
nz=s[2]/ny
z=reform(VCD,nx,ny,nz)
z=z[88:112,88:112,*]
x=x[88:112]
y=y[88:112]
;transform map coordinates from longitude and latitude to Cartesian (x, y) coordinates
;xloc=fltarr(n_elements(x),n_elements(y))
;yloc=fltarr(n_elements(x),n_elements(y))
;For num=0,n_elements(y)-1 do begin
;xloc[*,num]=x
;endfor
;For num=0,n_elements(x)-1 do begin
;yloc[num,*]=y
;endfor

;mapStruct = MAP_PROJ_INIT( 'Lambert Conic', $
;                              STANDARD_PAR1=25, STANDARD_PAR2=40, $
;                              CENTER_LONGITUDE=110, CENTER_LATITUDE=36 )
;orin=MAP_PROJ_FORWARD(curLon, curLat, MAP_STRUCTURE = mapStruct)/1000
;coordinates= MAP_PROJ_FORWARD( xloc, yloc, MAP_STRUCTURE = mapStruct)/1000
;xloc=coordinates[0,*]-curLon
;yloc=coordinates[1,*]-curLat

A=fltarr(nz,7)
For num=0,nz-1 do begin
	surface,z[*,*,num]
	image = tvrd(true =1)
	write_jpeg,'guassian_fit.jpg',image,true=1
	fit=GAUSS2DFIT( z[*,*,num],A[num,*],x,y )
endfor

;plot,A[0,*],psym=2,$
;yrange=[0,max(y)],$
;title='background VCD',xtitle='year',ytitle='VCD/(10^15moles/cm2)'

;plot,A[1,*],psym=2,$
;yrange=[0,max(A[1,*])],$
;title='VCD',xtitle='year',ytitle='VCD/(10^15moles/cm2)'

end
