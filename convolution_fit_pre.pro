pro convolution_fit_pre
;PRD
pp_loc_lon=[113.262863]
pp_loc_lat=[23.13152]

file_in_lat=string('/home/liufei/Data/Large_point/Steffen/latitude.nc')
fid=NCDF_OPEN(file_in_lat)
varid=NCDF_VARID(fid,'map_TVCD')
NCDF_VARGET,fid,varid,lat
NCDF_CLOSE, fid
;reverse the latitude
lat=reverse(lat)

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

file_in=string('/home/liufei/Data/Large_point/Steffen/map_TVCDstd.nc')
fid=NCDF_OPEN(file_in)
varid=NCDF_VARID(fid,'map_TVCDstd')
NCDF_VARGET,fid,varid,DATAstd
NCDF_CLOSE, fid

file_in=string('/home/liufei/Data/Large_point/Steffen/map_V.nc')
fid=NCDF_OPEN(file_in)
varid=NCDF_VARID(fid,'map_V')
NCDF_VARGET,fid,varid,wind
NCDF_CLOSE, fid

print,'size of DATA',size(DATA)

For pp=0,N_elements(pp_loc_lon)-1 do begin
lat_point=pp_loc_lat[pp]
lon_point=pp_loc_lon[pp]

;transform TVCD coordinates from longitude and latitude to Cartesian (x, y) coordinates
xloc=fltarr(n_elements(lon),n_elements(lat))
yloc=fltarr(n_elements(lon),n_elements(lat))
For num=0,n_elements(lat)-1 do begin
        xloc[*,num]=lon
endfor
For num=0,n_elements(lon)-1 do begin
        yloc[num,*]=lat
endfor


mapStruct = MAP_PROJ_INIT( 'Lambert Conic', $
                            STANDARD_PAR1=25, STANDARD_PAR2=40, $
                            CENTER_LONGITUDE=lon_point, CENTER_LATITUDE=lat_point )
coordinates= MAP_PROJ_FORWARD( xloc, yloc, MAP_STRUCTURE = mapStruct)/1000
lon_new=reform(coordinates[0,*],[N_Elements(lon),N_Elements(lat)])
lat_new=reform(coordinates[1,*],[N_Elements(lon),N_Elements(lat)])

;x direction interval/km
a=300
;y direction interval/km
b=120
;find the point which is nearest to the origin point
vector=fltarr(N_Elements(lon),N_Elements(lat))
For J= 0 , N_Elements(lat)-1 do begin
        For I= 0 , N_Elements(lon)-1 do begin
                vector[I,J]= distance_measure([[0,0],[lon_new[I,J],lat_new[I,J]]])

                if vector[I,J] eq 0 then begin
                        print,I,J,0
                endif

        endfor
endfor
print,'min of distance vector',min(vector)
nearest=where(vector eq min(vector))
nearest_new=nearest[0]
near_row= Fix( (nearest_new+1)/N_Elements(lon) )
near_col=( (nearest_new+1) mod N_Elements(lon))-1
print,'near_row',near_row,'near_col',near_col

;E_average is the mean of fitted E in mol/s, sample_num is the number of fit with R>=0.9
E_average=0
sample_num=0
t_average=0
;map_TVCD(lat,lon,sea,winddirection,windspeed)
;map_V(lat,lon,sea,winddirection,windspeed)
For dir = 0,8 do begin
        TVCD=fltarr(N_Elements(lon),N_Elements(lat))
        V=fltarr(N_Elements(lon),N_Elements(lat))
        STD=fltarr(N_Elements(lon),N_Elements(lat))
        For J= 0 , N_Elements(lat)-1 do begin
                For I= 0 , N_Elements(lon)-1 do begin
;************season**************
                        season=1
                        if (finite(DATA[J,I,season,dir,0]) eq 0) or (DATA[J,I,season,dir,0] lt 0) then begin
                                TVCD[I,J]=0
                        endif else begin
                                TVCD[I,J]= DATA[J,I,season,dir,0]
                        endelse

                        if (finite(DATAstd[J,I,season,dir,0]) eq 0) or (DATAstd[J,I,season,dir,0] lt 0) then begin
                                STD[I,J]=0
                        endif else begin
                                STD[I,J]= DATAstd[J,I,season,dir,0]
                        endelse


                        if (finite(wind[J,I,season,dir,0]) eq 0) then begin
                                V[I,J]=-999
                        endif else begin
                                V[I,J]= wind[J,I,season,dir,0]
                        endelse

                endfor
        endfor

