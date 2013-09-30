pro gaussian_fit_point_VCD

lat_point=23.13152
lon_point=113.262863

file_in_lat=string('/home/liufei/Data/Large_point/Steffen/latitude.nc')
fid=NCDF_OPEN(file_in_lat)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,lat
NCDF_CLOSE, fid

file_in_lon=string('/home/liufei/Data/Large_point/Steffen/longitude.nc')
fid=NCDF_OPEN(file_in_lon)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,lon
NCDF_CLOSE, fid

file_in=string('/home/liufei/Data/Large_point/Steffen/map_TVCD.nc')
fid=NCDF_OPEN(file_in)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,DATA
NCDF_CLOSE, fid

print,'size of DATA',size(DATA)

cellsize=0.167
cellnum_lat=7
cellnum_lon_up=6
cellnum_lon_down=12


loc_lat=where( (lat le lat_point+cellsize*cellnum_lat) and (lat ge lat_point-cellsize*cellnum_lat) )
loc_lon=where( (lon le lon_point+cellsize*cellnum_lon_down) and (lon ge lon_point-cellsize*cellnum_lon_up) )

;map_TVCD(lat,lon,sea,winddirection,windspeed)
For dir = 0,8 do begin
	TVCD=fltarr(N_Elements(loc_lon),N_Elements(loc_lat))
        For J= loc_lat[0] , loc_lat[N_Elements(loc_lat)-1] do begin
                For I=loc_lon[0] , loc_lon[N_Elements(loc_lon)-1] do begin
;;season
			season=3
			if (finite(DATA[J,I,season,dir,0]) eq 0) or (DATA[J,I,season,dir,0] lt 0) then begin
                        	TVCD[I-loc_lon[0],J-loc_lat[0]]=0
			endif else begin
				TVCD[I-loc_lon[0],J-loc_lat[0]]=DATA[J,I,season,dir,0] 
			endelse
		endfor
	endfor
endfor

undefine,season

namein=string(lat_point)+'N'+string(lon_point)+'E'
name=strcompress(namein,/REMOVE)

;write_jpeg,'test'+name+'LD_east & west.jpg',image,true=1


For dir = 0,8 do begin
	flag=0
	For J= loc_lat[0] , loc_lat[N_Elements(loc_lat)-1] do begin
		For I=loc_lon[0] , loc_lon[N_Elements(loc_lon)-1] do begin
;;season	
			season=where(DATA[J,I,3,dir,0] gt 0.0)
			if not array_equal(season,[-1]) then begin
				if flag then begin
					wind_temp=[wind_temp,mean( DATA[J,I,season,dir,0] )]
					x_lat=[x_lat,lat[J]]
					x_lon=[x_lon,lon[I]]
				endif else begin
					wind_temp=[mean( DATA[J,I,season,dir,0] )]
					x_lat=[lat[J]]
                 			x_lon=[lon[I]]
					flag=1
				endelse
			endif 
			undefine,season
		endfor
	endfor
	
	x_temp=fltarr(N_Elements(x_lat))
	For num=0,N_Elements(x_lat)-1 do begin
		if (x_lon[num] ge lon_point) and (x_lat[num] ge lat_point) then begin 
			x_temp[num]=Map_2Points(lon_point,lat_point,x_lon[num],x_lat[num],/METERS)/1000
		endif else begin
			x_temp[num]=-Map_2Points(lon_point,lat_point,x_lon[num],x_lat[num],/METERS)/1000
		endelse
	endfor
	
	case dir+1 of 
	1:wind_TVCD1=wind_temp/10
	2:wind_TVCD2=wind_temp/10
	3:wind_TVCD3=wind_temp/10
	4:wind_TVCD4=wind_temp/10
	5:wind_TVCD5='no_data'
	6:wind_TVCD6=wind_temp/10
	7:wind_TVCD7=wind_temp/10
	8:wind_TVCD8=wind_temp/10
	9:wind_TVCD9=wind_temp/10
	endcase
	
	case dir+1 of
	1:x1=x_temp
	2:x2=x_temp
	3:x3=x_temp
	4:x4=x_temp
	5:x5='no_data'
	6:x6=x_temp
	7:x7=x_temp
	8:x8=x_temp
	9:x9=x_temp
	endcase
	undefine,wind_temp
	undefine,x_temp
	
endfor

		

;wind=[[wind_TVCD1,wind_TVCD9],[wind_TVCD2,wind_TVCD8],[wind_TVCD3,wind_TVCD7],[wind_TVCD4,wind_TVCD6]]
;x=[[x1,x9],[x2,x8],[x3,x7],[x4,x6]]
PLOT,x1,wind_TVCD1,psym=2,$
	title='southeast & northwest',xtitle='x(km)',ytitle='NO2 TVCD(10^15 molec/cm2)'
OPLOT,x9,wind_TVCD9,psym=4
image = tvrd(true =1)
write_jpeg,name+'_southeast & northwest.jpg',image,true=1

PLOT,x2,wind_TVCD2,psym=2,$
        title='south & north',xtitle='x(km)',ytitle='NO2 TVCD(10^15 molec/cm2)'
OPLOT,x8,wind_TVCD8,psym=4
image = tvrd(true =1)
write_jpeg,name+'_south & north.jpg',image,true=1

PLOT,x3,wind_TVCD3,psym=2,$
        title='southwest & northeast',xtitle='x(km)',ytitle='NO2 TVCD(10^15 molec/cm2)'
OPLOT,x7,wind_TVCD7,psym=4
image = tvrd(true =1)
write_jpeg,name+'_southwest & northesat.jpg',image,true=1

PLOT,x4,wind_TVCD4,psym=2,$
        title='east & west',xtitle='x(km)',ytitle='NO2 TVCD(10^15 molec/cm2)'
OPLOT,x6,wind_TVCD6,psym=4
image = tvrd(true =1)
write_jpeg,name+'_east & west.jpg',image,true=1


end


