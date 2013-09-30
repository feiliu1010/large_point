function gaus_convol,x,y,sig,nsigma=nsigma

on_error,2
c=check_math(0,1)
if n_params() ne 3 then message,'usage: result=gaus_convol(x,y,sig)'

if n_elements(nsigma) ne 1 then nsigma=2.5      ; approximate width in
                                                ; units of sigma
nsigma2= NSIGMA*2
n      = n_elements(x)
conv   = (max(double(x))-min(double(x)))/(n -1) ; conversion, units/point
n_pts  = ceil(nsigma2*sig/conv)                 ; number of points
n_pts  = (n_pts > 2) < (n-2)                    ; restrictions on n_pts
if (n_pts/2*2) eq n_pts then n_pts = n_pts +1   ; make odd number of points
xvar = (dindgen(n_pts)/(n_pts-1)-0.5d0)*n_pts   ; approx. - NSIGMA < x < +NSIGMA
gaus=exp(-.5*(xvar/(sig/conv))^2)               ; gaussian of width sig.
return,convol(y,gaus,/center)/total(gaus)       ; do convolution.
end

;fit a function of the form F(x) = E*(e(x) convolute G(x))+B  to sample pairs contained in array x(distance) and y(observed NO2 line densities).
;e(x)=exp(-(x-Xp)/x0) or 0 (x<Xp)
;G(x)=1/(sqrt(2*!DPI)*theta)*exp(-x^2/(2*theta^2))
;define a procedure to return F(x),given x. Note that A is an array containing the values x0,E,B,Xp and theta
function myfunct,X,P
	amp=[0.1D,1.D,1.D,0.1D,0.1D]
        ;P=[x0,E,B,Xp,theta]
        ex=exp(-(X-P[3]/amp[3])/(P[0]/amp[0]))
        ;define an array to truncate ex
        zero=fltarr(N_elements(X))
        zero[where(X ge P[3]/amp[3])]=1
        ex=ex*zero
        model= P[1]/amp[1]*GAUS_CONVOL(X,ex,P[4]/amp[4],nsigma=1.5)+P[2]/amp[2]
        return,model
end

pro gaussian_fit_point

;Locations: P1 (30.0.N, 122.1.E); P2 (29.2.N, 119.5.E); P3 (40.3.N, 111.3.E); P4 (38.8.N, 110.2.E); P5 (23.1.N, 109.8.E); P6 (30.5.N, 106.8.E)
;pp_loc_lon=[121.600103,121.145434,118.125528,116.9231,116.24179,114.94014,113.295,112.5738889,111.3544444]
;pp_loc_lat=[31.352795,30.628704,24.30395,35.326499,37.44928,40.660743,40.02583333,35.46555556,40.19444444]
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
;For pp=7,7 do begin
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
print,'vector_zero',where (vector eq 0) 
print,'min of vector',min(vector)
nearest=where(vector eq min(vector))
nearest_new=nearest[0]
near_row= Fix( (nearest_new+1)/N_Elements(lon) )
near_col=( (nearest_new+1) mod N_Elements(lon))-1 
print,'near_row',near_row,'near_col',near_col

;E_average is the mean of fitted E in mol/s, sample_num is the number of fit with R>=0.9
E_average=0
sample_num=0
t_average=0
wind_average=0
;map_TVCD(lat,lon,sea,winddirection,windspeed)
;map_V(lat,lon,sea,winddirection,windspeed)
For dir = 0,8 do begin
	TVCD=fltarr(N_Elements(lon),N_Elements(lat))
	V=fltarr(N_Elements(lon),N_Elements(lat))
	STD=fltarr(N_Elements(lon),N_Elements(lat))
        For J= 0 , N_Elements(lat)-1 do begin
                For I= 0 , N_Elements(lon)-1 do begin
;************season
			season=0
			;if (finite(DATA[J,I,season,dir,0]) eq 0) or (DATA[J,I,season,dir,0] lt 0) then begin
			if (finite(DATA[J,I,season,dir,0]) eq 0) then begin
				TVCD[I,J]=-999
			endif else begin
				TVCD[I,J]= DATA[J,I,season,dir,0]
			endelse

                        if (finite(DATAstd[J,I,season,dir,0]) eq 0) or (DATAstd[J,I,season,dir,0] lt 0) then begin
                                STD[I,J]=0
                        endif else begin
                                STD[I,J]= DATAstd[J,I,season,dir,0]
                        endelse


			if (finite(wind[J,I,season,dir,0]) eq 0) or (wind[J,I,season,dir,0] gt 15)then begin
                        	V[I,J]=-999
                        endif else begin
	                        V[I,J]= wind[J,I,season,dir,0]
                        endelse

		endfor
	endfor
	;reverse TVCD at latitude
	TVCD=reverse(TVCD,2)
	STD=reverse(STD,2)
	V=reverse(V,2)
	;rotate the coordinates of lon_new, lat_new
	case dir+1 of 
	    1:theta=135
	    2:theta=90
	    3:theta=45
            4:theta=0
	    5:theta=0
	    6:theta=0
	    7:theta=45
	    8:theta=90
	    9:theta=135
	endcase
        theta=theta/!radeg
        c=cos(theta)
        s=sin(theta)
	lon_new_rot=fltarr(N_Elements(lon),N_Elements(lat))
	lat_new_rot=fltarr(N_Elements(lon),N_Elements(lat))
        For J= 0 , N_Elements(lat)-1 do begin
                For I= 0 , N_Elements(lon)-1 do begin
			lon_new_rot[I,J]=c*lon_new[I,J]+s*lat_new[I,J]
			lat_new_rot[I,J]=c*lat_new[I,J]-s*lon_new[I,J]
		endfor
        endfor
	;find boundary of interval
	bou_y=where(vector le b)
	bou_x=where(vector le a)
	y_temp_row= Fix( (bou_y+1)/N_Elements(lon) )
	y_temp_col=( (bou_y+1) mod N_Elements(lon))-1
        x_temp_row= Fix( (bou_x+1)/N_Elements(lon) )
        x_temp_col=( (bou_x+1) mod N_Elements(lon))-1 

	;east-west
	if (dir eq 3) or (dir eq 5) then begin
		max_y1= max( y_temp_row )
		max_y2= min( y_temp_row )
		max_x1= max( x_temp_col )
		max_x2= min( x_temp_col )
		num_y= max_y1 - max_y2 +1
		num_x= max_x1 - max_x2 +1
	
		;Line Density:LD
		LD=fltarr(num_x)
		L_Num=fltarr(num_x)
		;Number of NAN data
		L_NAN=fltarr(num_x)
		;STD Line Density:SD
                SD=fltarr(num_x)
		;wind Line Density
		WD=fltarr(num_x)
		;pixel number with valid wind speed
		W_Num=fltarr(num_x)
		flag_axle=0
		For x=max_x2 , max_x1 do begin
			For y=max_y2 , max_y1 do begin
			    if y eq 0 then begin
				dy =abs(lat_new[x,y+1]-lat_new[x,y])
			    endif else if y eq N_elements(lat)-1 then begin
				dy=abs(lat_new[x,y]-lat_new[x,y-1])
			    endif else begin
				dy=abs( (lat_new[x,y+1]-lat_new[x,y-1]) )/2
			    endelse
				if TVCD[x,y] ne -999 then begin
                                        LD[x-max_x2]+=TVCD[x,y]*dy/10000
                                        L_num[x-max_x2]+=1
                                endif else begin
					L_NAN[x-max_x2]+=1
				endelse
				SD[x-max_x2]+=STD[x,y]*dy/10000
				if V[x,y] ne -999 then begin
					WD[x-max_x2]+=V[x,y]
					W_num[x-max_x2]+=1
				endif		
			endfor
			x_axle_r= lon_new[max_x2:max_x1,near_row]
		endfor
		;convert 2D x_axle_r to 1D x_axle
                x_axle= x_axle_r[0:n_elements(x_axle_r)-1]
		;print,'east-west',LD
	;south-north
	endif else if (dir eq 1) or (dir eq 7) then begin
                max_y1= max( y_temp_col )
                max_y2= min( y_temp_col )
                max_x1= max( x_temp_row )
                max_x2= min( x_temp_row )
                num_y= max_y1 - max_y2 +1
                num_x= max_x1 - max_x2 +1
                ;Line Density:LD
                LD=fltarr(num_x)
		L_num=fltarr(num_x)
		SD=fltarr(num_x)
                ;wind Line Density
                WD=fltarr(num_x)
                ;pixel number with valid wind speed
                W_num=fltarr(num_x)
                flag_axle=0
                For x=max_x2 , max_x1 do begin
                        For y=max_y2 , max_y1 do begin
			    if y eq 0 then begin
				dy=abs( lon_new[y+1,x]-lon_new[y,x] )
                            endif else if y eq N_elements(lon)-1 then begin
                                dy=abs( lon_new[y,x]-lon_new[y-1,x])
			    endif else begin
                                dy=abs( (lon_new[y+1,x]-lon_new[y-1,x]) )/2
			    endelse
				if TVCD[y,x] ne -999 then begin
	                                LD[x-max_x2]+=TVCD[y,x]*dy/10000
					L_num[x-max_x2]+=1
				endif else begin
					L_NAN[x-max_x2]+=1
				endelse
				SD[x-max_x2]+=STD[y,x]*dy/10000
				if V[y,x] ne -999 then begin
                                        WD[x-max_x2]+=V[y,x]
                                        W_num[x-max_x2]+=1
                                endif
                        endfor
                        x_axle_r= lat_new[near_col,max_x2:max_x1]
                endfor
		;convert 2D x_axle_r to 1D x_axle
		x_axle= x_axle_r[0:n_elements(x_axle_r)-1]
		;print,'south-north',LD
	;southwest-northeast
	endif else if (dir eq 2) or (dir eq 6) then begin
		num_y_positive=0
	        num_y_negative=0
        	x_positive=[near_col,near_row]
	        x_negative=[near_col,near_row]
		For count=1,min([N_Elements(lon),N_Elements(lat)]) do begin
			if (not array_equal(where(y_temp_col eq (near_col-count)),-1))$
			   and(not array_equal(where(y_temp_row eq (near_row+count)),-1)) then begin
				num_y_positive+=1
			endif
			if (not array_equal(where(y_temp_col eq (near_col+count)),-1))$
                           and(not array_equal(where( y_temp_row eq (near_row-count)),-1)) then begin
                                num_y_negative+=1
			endif
		

			if (not array_equal(where(x_temp_col eq near_col+count-1),-1))$
                           and(not array_equal(where(x_temp_row eq (near_row+count)),-1)) then begin
                                x_positive=[x_positive,[near_col+count-1,near_row+count]]
			endif
			if (not array_equal(where(x_temp_col eq (near_col+count)),-1))$
                           and(not array_equal(where( x_temp_row eq (near_row+count)),-1)) then begin
                                x_positive=[x_positive,[near_col+count,near_row+count]]
                        endif
                        if (not array_equal(where(x_temp_col eq (near_col-count)),-1))$
                           and(not array_equal(where(x_temp_row eq near_row-count+1),-1)) then begin
                               x_negative=[x_negative,[near_col-count,near_row-count-1]]
			endif
			if (not array_equal(where( x_temp_col eq (near_col-count)),-1))$
                           and(not array_equal(where(x_temp_row eq (near_row-count)),-1)) then begin
                               x_negative=[x_negative,[near_col-count,near_row-count]]
                        endif
		endfor
		;print,'num_y_positive',num_y_positive,'num_y_negative',num_y_negative
		x_positive=reform(x_positive,[2,N_Elements(x_positive)/2])
		x_negative=reform(x_negative,[2,N_Elements(x_negative)/2])
		x_negative_r=reverse(x_negative,2)
		;exclude last parameter of x_negative_r in case double-counting
		x_negative_new=x_negative_r[*,0:N_elements(x_negative)/2-2]
		x_sequence=[[x_negative_new],[x_positive]]
		;print,'size of x_sequence',size(x_sequence)
		;print,'x_sequence',x_sequence
		LD=fltarr(N_elements(x_sequence)/2)
		L_num=fltarr(N_elements(x_sequence)/2)
		L_NAN=fltarr(N_elements(x_sequence)/2)
		SD=fltarr(N_elements(x_sequence)/2)
                ;wind Line Density
                WD=fltarr(N_elements(x_sequence)/2)
                ;pixel number with valid wind speed
                W_num=fltarr(N_elements(x_sequence)/2)
 
		x_axle=fltarr(N_elements(x_sequence)/2)
		For x=0,N_elements(x_sequence)/2-1 do begin
			I=x_sequence[0,x]
			J=x_sequence[1,x]
			For y=0, num_y_positive do begin
			  ;pay attention to the database boundary
			  if ( I-y ge 0) and ( J+y le  N_elements(lat)-1) then begin
			    if (I-y eq 0) and (J+y ne 0) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y+1,J+y-1])
			    endif else if (I-y eq 0) and (J+y eq 0) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y+1,J+y+1])
			    endif else if (I-y ne 0) and (J+y eq 0) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y-1,J+y+1])
			    endif else if (I-y eq N_elements(lon)-1) and (J+y ne N_elements(lat)-1) then begin
				dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y-1,J+y+1])
                            endif else if (I-y eq N_elements(lon)-1) and (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y-1,J+y-1])
                            endif else if (I-y ne N_elements(lon)-1) and (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J+y]-lat_new_rot[I-y+1,J+y-1])
			    endif else begin				
				dy=abs( (lat_new_rot[I-y-1,J+y+1]-lat_new_rot[I-y+1,J+y-1]) )/2
			    endelse
				if TVCD[I-y,J+y] ne -999 then begin
					LD[x]+= TVCD[I-y,J+y]*dy/10000
					L_num[x]+=1
				endif else begin
					L_NAN[x]+=1
				endelse
				SD[x]+= STD[I-y,J+y]*dy/10000
                                if V[I-y,J+y] ne -999 then begin
                                        WD[x]+=V[I-y,J+y]
                                        W_num[x]+=1
                                endif   
			  endif
			endfor
			For y=1,num_y_negative do begin
                          ;pay attention to the database boundary
                          if (J-y ge 0) and (I+y le  N_elements(lon)-1) then begin
                            if (I+y eq 0) and (J-y ne 0) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y+1,J-y-1])
                            endif else if (I+y eq 0) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y+1,J-y+1])
                            endif else if (I+y ne 0) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y-1,J-y+1])
                            endif else if (I+y eq N_elements(lon)-1) and (J-y ne N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y-1,J-y+1])
                            endif else if (I+y eq N_elements(lon)-1) and (J-y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y-1,J-y-1])
                            endif else if (I+y ne N_elements(lon)-1) and (J-y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J-y]-lat_new_rot[I+y+1,J-y-1])
                            endif else begin
				dy=abs( (lat_new_rot[I+y-1,J-y+1]-lat_new_rot[I+y+1,J-y-1]) )/2
                            endelse
				if TVCD[I+y,J-y] ne -999 then begin
					LD[x]+= TVCD[I+y,J-y]*dy/10000
					L_num[x]+=1
				endif else begin
					L_NAN[x]+=1
				endelse
				SD[x]+=STD[I+y,J-y]*dy/10000
                                if V[I+y,J-y] ne -999 then begin
                                        WD[x]+=V[I+y,J-y]
                                        W_num[x]+=1
                                endif
  
			  endif				
			endfor
			x_axle[x]=lon_new_rot[x_sequence[0,x],x_sequence[1,x]]
		endfor
        ;southeast-northwest
        endif else if (dir eq 0) or (dir eq 8) then begin
                num_y_positive=0
                num_y_negative=0
                x_positive=[near_col,near_row]
                x_negative=[near_col,near_row]
                For count=1,min([N_Elements(lon),N_Elements(lat)]) do begin
                        if (not array_equal(where(y_temp_col eq (near_col-count)),-1))$
                           and(not array_equal(where(y_temp_row eq (near_row-count)),-1)) then begin
                                num_y_positive+=1
                        endif
                        if (not array_equal(where(y_temp_col eq (near_col+count)),-1))$
                           and(not array_equal(where( y_temp_row eq (near_row+count)),-1)) then begin
                                num_y_negative+=1
                        endif


                        if (not array_equal(where(x_temp_col eq near_col-count),-1))$
                           and(not array_equal(where(x_temp_row eq (near_row+count-1)),-1)) then begin
                                x_positive=[x_positive,[near_col-count,near_row+count-1]]
                        endif
                        if (not array_equal(where(x_temp_col eq (near_col-count)),-1))$
                           and(not array_equal(where( x_temp_row eq (near_row+count)),-1)) then begin
                                x_positive=[x_positive,[near_col-count,near_row+count]]
                        endif
                        if (not array_equal(where(x_temp_col eq (near_col+count-1)),-1))$
                           and(not array_equal(where(x_temp_row eq near_row-count),-1)) then begin
                               x_negative=[x_negative,[near_col+count-1,near_row-count]]
                        endif
                        if (not array_equal(where( x_temp_col eq (near_col+count)),-1))$
                           and(not array_equal(where(x_temp_row eq (near_row-count)),-1)) then begin
                               x_negative=[x_negative,[near_col+count,near_row-count]]
                        endif
                endfor
                ;print,'num_y_positive',num_y_positive,'num_y_negative',num_y_negative
                x_positive=reform(x_positive,[2,N_Elements(x_positive)/2])
                x_negative=reform(x_negative,[2,N_Elements(x_negative)/2])
                x_negative_r=reverse(x_negative,2)
                ;exclude last parameter of x_negative_r in case double-counting
                x_negative_new=x_negative_r[*,0:N_elements(x_negative)/2-2]
                x_sequence=[[x_negative_new],[x_positive]]
                ;print,'size of x_sequence',size(x_sequence)
                ;print,'x_sequence',x_sequence
                LD=fltarr(N_elements(x_sequence)/2)
		L_num=fltarr(N_elements(x_sequence)/2)
		L_NAN=fltarr(N_elements(x_sequence)/2)
		SD=fltarr(N_elements(x_sequence)/2)
                ;wind Line Density
                WD=fltarr(N_elements(x_sequence)/2)
                ;pixel number with valid wind speed
                W_num=fltarr(N_elements(x_sequence)/2)

                x_axle=fltarr(N_elements(x_sequence)/2)
                For x=0,N_elements(x_sequence)/2-1 do begin
                        I=x_sequence[0,x]
                        J=x_sequence[1,x]
                        For y=0, num_y_positive do begin
                          ;pay attention to the database boundary
                          if ( I-y ge 0) and ( J-y ge 0) then begin
                            if (I-y eq 0) and (J-y ne  N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y+1,J-y+1])
                            endif else if (I-y eq 0) and (J-y eq  N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y+1,J-y])
                            endif else if (I-y ne N_elements(lon)-1) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y+1,J-y+1])
                            endif else if (I-y eq N_elements(lon)-1) and (J-y eq 0) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y,J-y+1])
                            endif else if (I-y eq N_elements(lon)-1) or (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I-y,J-y]-lat_new_rot[I-y-1,J-y-1])
                            endif else begin
                                dy=abs( (lat_new_rot[I-y-1,J-y-1]-lat_new_rot[I-y+1,J-y+1]) )/2
                            endelse
			    if TVCD[I-y,J-y] ne -999 then begin
	                            LD[x]+= TVCD[I-y,J-y]*dy/10000
				    L_num[x]+=1
			    endif else begin
				    L_NAN[x]+=1
			    endelse
			    SD[x]+= STD[I-y,J-y]*dy/10000
                            if V[I-y,J-y] ne -999 then begin
                                WD[x]+=V[I-y,J-y]
			        W_num[x]+=1
  			    endif
			  endif	
                        endfor

                        For y=1,num_y_negative do begin
                          ;pay attention to the database boundary
                          if (J+y le  N_elements(lat)-1) and (I+y le  N_elements(lon)-1) then begin
                            if (I+y eq 0) and (J+y ne  N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y+1,J+y+1])
                            endif else if (I+y eq 0) and (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y+1,J+y])
                            endif else if (I+y ne  N_elements(lon)-1) and (J+y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y+1,J+y+1])
                            endif else if (I+y eq N_elements(lon)-1) and (J+y eq 0) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y,J+y+1])
                            endif else if (I+y eq N_elements(lon)-1) or (J+y eq N_elements(lat)-1) then begin
                                dy=abs(lat_new_rot[I+y,J+y]-lat_new_rot[I+y-1,J+y-1])
                            endif else begin
                                dy=abs( (lat_new_rot[I+y+1,J+y+1]-lat_new_rot[I+y-1,J+y-1]) )/2
                            endelse
			    if TVCD[I+y,J+y] ne -999 then begin
                                LD[x]+= TVCD[I+y,J+y]*dy/10000
				L_num[x]+=1
			    endif else begin
				L_NAN[x]+=1
			    endelse
				SD[x]+= STD[I+y,J+y]*dy/10000
                            if V[I+y,J+y] ne -999 then begin
                                WD[x]+=V[I+y,J+y]
                                W_num[x]+=1
                            endif
 
                          endif
                        endfor
                        x_axle[x]=lon_new_rot[x_sequence[0,x],x_sequence[1,x]]
                endfor
		;print,'southeast-northwest',LD	
	endif else begin
		LD=0
		x_axle=0
	endelse
	
        if (dir eq 3) or (dir eq 6) or (dir eq 7) or (dir eq 8) then begin
                x_axle=-x_axle
        endif


	if dir ne 4 then begin
		;print,'LD',LD
		;print,'x_axle',x_axle
		old_num=N_elements(LD)
		LD=LD*10.0D
		SD=SD*10.0D
		;sort x_axle ascend
		LD=LD[sort(x_axle)]
		L_num=L_num[sort(x_axle)]
		L_NAN=L_NAN[sort(x_axle)]
		unvalid=L_NAN/(L_num+L_NAN)	
		SD=SD[sort(x_axle)]
		WD=WD[sort(x_axle)]
		W_num=W_num[sort(x_axle)]
		x_axle=x_axle[sort(x_axle)]
		print,'L_num',total(L_num)
		print,'W_num',total(W_num)
		;distinguish upwind & downwind
		upwind=-200
		downwind=300
		;discard data in y direction has >=30% NAN
		;range=where( (x_axle le downwind) and  (x_axle gt upwind) and (unvalid le 0.3) )
		range=where(unvalid le 0.3)
		LD=LD[range]
		L_num=L_num[range]
		print,'L_num',total(L_num)
		SD=SD[range]
		WD=WD[range]
		W_num=W_num[range]
		print,'W_num',total(W_num)
		x_axle=x_axle[range]
		;integrate wind speed and deduce mean
		wind_mean=total(WD)/total(W_num)
		print,'wind_mean',wind_mean
		;max_LD=max(LD[where(x_axle gt -10)])
		max_LD=max(LD)
		min_LD=min(LD)
		distance=x_axle[where(LD eq max_LD)]
		print,'distance',distance
		;discard the situation when LD_max is at the position of x_axle_max
		if distance ne  max(x_axle) then begin 
    			ERR=make_array(N_elements(x_axle),value=0.05)
        		;ERR=SD/LD
			weights = make_array(N_elements(x_axle),value=1.0)
	       		;Provide an initial guess of the function's parameters.
			amp=[0.1D,1.D,1.D,0.1D,0.1D]
	        	start = [30.D*amp[0],max_LD*amp[1],min_LD*amp[2],distance*amp[3],20.D*amp[4]]
		        parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                		       limits:[0.D,0],relstep:0.D,mpside:0}, 5)

		        parinfo[0].limited[0] = 1
	        	parinfo[0].limits[0]  =0.D*amp[0]
		        parinfo[1].limited[*] = 1
                        parinfo[1].limits[0]  = 0.D
                        parinfo[1].limits[1]  = 15.D
			parinfo[2].limited[*] = 1
		        parinfo[2].limits[0]  =0.D*amp[2]
		       	parinfo[2].limits[1]  =(min_LD+0.5D)*amp[2]
		        parinfo[3].limited[*] = 1
		        parinfo[3].limits[0]  =(distance-10.D)*amp[3]
		        parinfo[3].limits[1]  =(distance+10.D)*amp[3]
		       	parinfo[4].limited[0] = 1
		        parinfo[4].limits[0]  = 0.D*amp[4]

		        parinfo[*].RELSTEP=0.01
		        parinfo[*].mpside=2
		       	parinfo[*].value = [30.D*amp[0],max_LD*amp[1],min_LD*amp[2],distance*amp[3],20.D*amp[4]]
		        print,'start',parinfo[*].value
		        P = MPFITFUN('myfunct',x_axle,LD,ERR,start,WEIGHTS=weights,PARINFO=parinfo,PERROR=perror,yfit=yfit,/QUIET)
			;print,'Para_ERROR',PERROR
			print,'P',dir+1,':',P
			r=CORRELATE(LD,yfit)
			print,'r',r

			space=findgen(600)*0.5-100
			F=spline(x_axle,yfit,space)
			
			endif else begin
				print,'dir',dir+1,'LD_max is at the position of x_axle_max'
				F= make_array(N_elements(x_axle),value=0)
			endelse

		namein=string(lat_point)+'N'+string(lon_point)+'E'+'season'+string(season+1)+'dir'+string(dir+1)
		name=strcompress(namein,/REMOVE)
		PLOT,x_axle,LD,psym=2,$
		        title='wind direction'+string(dir+1),xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)',$
			xrange=[-100,200]
		OPLOT,space,F
		;lifetime,unit:h
		t=(P[1]/amp[1])/wind_mean
		;t=3.5
		;wind_mean=(P[1]/amp[1])/t
		;emit rate
		E=(P[1]/amp[1]*wind_mean*100)/6.023
		if r ge 0.9 then begin
			E_average+=E
			sample_num+=1
			t_average+=t
			wind_average+=wind_mean
		endif
		xyouts,-80,max_LD/1.1,$
			'E(mol/s):'+string(E)+'!Clifetime:'+string(t)+'!CR:'+string(r)+'!Ce-folding distance downwind:'+string(P[0])+'!Cwindspeed:'+string(wind_mean),$
			charsize=1.5, charthick=1.8
		image = tvrd(true =1)
		write_jpeg,name+'LD.jpg',image,true=1
	
       endif else begin
		F=0
       endelse
	
	case dir+1 of
	1:LD1=LD
	2:LD2=LD 
	3:LD3=LD
	4:LD4=LD
	5:LD5=LD
	6:LD6=LD
	7:LD7=LD
	8:LD8=LD
	9:LD9=LD
	endcase

	case dir+1 of 
	1:x_axle1=x_axle
	2:x_axle2=x_axle
	3:x_axle3=x_axle
	4:x_axle4=-x_axle
	5:x_axle5=x_axle
	6:x_axle6=x_axle
	7:x_axle7=-x_axle
	8:x_axle8=-x_axle
	9:x_axle9=-x_axle
	endcase

        case dir+1 of
        1:space1=space
        2:space2=space
        3:space3=space
        4:space4=space
        5:space5=space
        6:space6=space
        7:space7=space
        8:space8=space
        9:space9=space
        endcase


	case dir+1 of
        1:F1=F
        2:F2=F
        3:F3=F
        4:F4=F
        5:F5=F
        6:F6=F
        7:F7=F
        8:F8=F
        9:F9=F
        endcase

endfor

E_average=E_average/sample_num
t_average=t_average/sample_num
wind_average=wind_average/sample_num
print,'E_average:',E_average,'lifetime:',t_average,'sample_num:',sample_num,'wind_average:',wind_average


namein=string(lat_point)+'N'+string(lon_point)+'E'
name=strcompress(namein,/REMOVE)
PLOT,x_axle1,LD1,psym=2,$
        title='southeast & northwest',xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)',$
	xrange=[-200,200],yrange=[0,max([LD1,F1,LD9,F9])+0.5]
OPLOT,space1,F1
OPLOT,-x_axle9,LD9,psym=4
OPLOT,-space9,F9
image = tvrd(true =1)
write_jpeg,name+'LD_southeast & northwest.jpg',image,true=1


PLOT,x_axle2,LD2,psym=2,$
        title='south & north',xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)',$
	xrange=[-200,200],yrange=[0,max([LD2,F2,LD8,F8])+0.5]
OPLOT,space2,F2
OPLOT,-x_axle8,LD8,psym=4
OPLOT,-space8,F8
image = tvrd(true =1)
write_jpeg,name+'LD_south & north.jpg',image,true=1

PLOT,x_axle3,LD3,psym=2,$
        title='southwest & northeast',xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)',$
	xrange=[-200,200],yrange=[0,max([LD3,F3,LD7,F7])+0.5]
OPLOT,space3,F3
OPLOT,-x_axle7,LD7,psym=4
OPLOT,-space7,F7
image = tvrd(true =1)
write_jpeg,name+'LD_southwest & northeast.jpg',image,true=1

PLOT,-x_axle4,LD4,psym=2,$
        title='east & west',xtitle='x(km)',ytitle='NO2 Line Density(10^23 molec/cm)',$
	xrange=[-200,200],yrange=[0,max([LD4,F4,LD6,F6])+0.5]
OPLOT,-space4,F4
OPLOT,x_axle6,LD6,psym=4
OPLOT,space6,F6
image = tvrd(true =1)
write_jpeg,name+'LD_east & west.jpg',image,true=1

endfor

end


