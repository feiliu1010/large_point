pro test
data=fltarr(3,2647)
filename ='/home/liufei/r5/large_point/data_test.asc'
openr,lun,filename,/get_lun
readf,lun,data
free_lun,lun

PRINT,size(data)
xloc=data[0,*]
yloc=data[1,*]
z=data[2,*]

curLon =-111.6
curLat =57
mapStruct = MAP_PROJ_INIT( 'Lambert Conic', $
                              STANDARD_PAR1=25, STANDARD_PAR2=40, $
                              CENTER_LONGITUDE=110, CENTER_LATITUDE=36 )
orin=MAP_PROJ_FORWARD(curLon, curLat, MAP_STRUCTURE = mapStruct)/1000
coordinates= MAP_PROJ_FORWARD( xloc, yloc, MAP_STRUCTURE = mapStruct)/1000
print,'orin',orin
xloc=coordinates[0,*]-orin[0]
yloc=coordinates[1,*]-orin[1]

data_test=fltarr(3,2647)
data_test[0,*]=xloc
data_test[1,*]=yloc
data_test[2,*]=z
outfile = '/home/liufei/r5/large_point/data_test2.asc'
openw,lun,outfile,/get_lun
printf,lun,data_test
close,lun
free_lun,lun
;z2=krig2d(z,xloc,yloc,expon=[10,0.],nx=40,ny=40)
;z2=GRID_TPS(xloc,yloc,z)

;generage net
cell_x=20
cell_y=10
nx=ceil( (max(xloc)-min(xloc))/cell_x )
ny=ceil( (max(yloc)-min(yloc))/cell_y )
net_start=[min(xloc),min(yloc)]
print,'nx',nx,'ny',ny,'net_start',net_start
x_new=fltarr(nx)
y_new=fltarr(ny)
x_corner=fltarr(nx)
y_corner=fltarr(ny)
For I=0,nx-1 do begin
        x_new[I] = net_start[0]+I*cell_x
	x_corner[I]=net_start[0]-0.5*cell_x+I*cell_x
endfor

For J=0,ny-1 do begin
        y_new[J] = net_start[1]+J*cell_y
	y_corner[J]=net_start[1]-0.5*cell_y+J*cell_y
endfor


z_new=fltarr(nx,ny)
z_new_num=fltarr(nx,ny)
Iloc=VALUE_LOCATE (x_corner, xloc)
Jloc=VALUE_LOCATE (y_corner, yloc)
For num=0,2647-1 do begin
	z_new[Iloc[num],Jloc[num]]+=z[num]
	z_new_num[Iloc[num],Jloc[num]]+=1
endfor

no_data=0
For I=0,nx-1 do begin
	For J=0,ny-1 do begin
		if  z_new_num[I,J] ne 0.0 then begin
			 z_new[I,J]= z_new[I,J]/z_new_num[I,J]
;		endif
		endif else begin
			no_data+=1
;			 z_new[I,J]= !Values.F_NaN	

		endelse
	endfor
endfor
print,'no_data',no_data

data_test=z_new
outfile = '/home/liufei/r5/large_point/z_new.asc'
openw,lun,outfile,/get_lun
printf,lun,data_test
close,lun
free_lun,lun

surface,z_new,charsize=4

fit=GAUSS2DFIT( Z_new,A, x_new,y_new )
print,A[0],A[1],A[2],A[3]

end
