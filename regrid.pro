pro regrid
; Elapsed time
; For one month, 36km, whole china       : 100.11265 seconds (CPU 100%
; For one month, 4km,  Beijing           : 43.394784 seconds (CPU 100%)
; For year 2005-2011:  0.05deg, 200*200  : 11260.986 seconds (CPU 40-50%)

starttime = systime(/seconds)
;-----------------------------------------------------------------------
; Setting up directories and common environment variable
;-----------------------------------------------------------------------
; Set up directories
  DataDir = '/z6/satellite/OMI/no2/DOMINO_S_v2/'  ; note: end with '/'
  OutDir = 'results/'                                ; note: end with '/'

; Set up date to process
  STDATE = 200501    ; start date: YYYYMM
  EDDATE = 201112    ;   end date: YYYYMM

; remove large pixels
  Nlp = 10           ; the number of edge pixels to be removed for each side
;  Nlp = 0           ; Nlp = 0, take into account all pixels
  L_Month = 'Y'      ; Switch: generate monthly mean from original swath data?
  L_Year  = 'Y'      ; Switch: generate annual mean From monthly mean?
  Month_st = 6       ; Start month for annual mean calculation
  Month_ed = 8       ; end month for annual mean calculation

; Set up output coordinates & Domain
; Lambert Conic Projection

  L_PROJ = 'N'      ; Switch: use projection coordinates?

;  mapStruct = MAP_PROJ_INIT( 'Lambert Conformal Conic ', /GCTP, $
;                              DATUM='GRS 1980/WGS 84', $
  mapStruct = MAP_PROJ_INIT( 'Lambert Conic', $
                              STANDARD_PAR1=25, STANDARD_PAR2=40, $
                              CENTER_LONGITUDE=110, CENTER_LATITUDE=36 )

; output Domain
; 36km, whole china
;  Res = 36000      ; grid resolution
;  Row = 113        ; Number of Rows
;  Col = 169        ; Number of Columns
;  x0 = -3042000    ; Origin X coordinate(bottom right corner)
;  y0 = -2034000    ; Origin Y coordinate(bottom right corner)

; 4km, Beijing
;  Res = 4000
;  Row = 63
;  Col = 90
;  x0 = 378000
;  y0 = 378000
;  Row = 200
;  Col = 300
;  x0 = 0
;  y0 = 0

; 0.05deg, Beijing
;  Res = 0.05
;  Row = 200
;  Col = 200
;  x0 = 115
;  y0 = 39

; 0.05deg, Yangcheng(35.46N, 112.58E)
; Res = 0.05
 Res = 0.016666666667
 Row = 200
 Col = 200
 x0 = 121.924444 - Res*Col/2
 y0 = 42.790278 - Res*Row/2

;-----------------------------------------------------------------------
; Begin here
;-----------------------------------------------------------------------

;  Res = long(Res)    ; important!!!  long,int
;  xCorner = LonArr(Col+1)
;  xCenter = LonArr(Col)
;  yCorner = LonArr(Row+1)
;  yCenter = LonArr(Row)
;  for i = 1,Col do begin
;    xCorner[i-1] = long(x0 + Res*(i-1))
;    xCenter[i-1] = long(xCorner[i-1] + Res/2)
;  endfor
;  xCorner[Col] = long(xCorner[Col-1] + Res)
;  for j = 1,Row do begin
;    yCorner[j-1] = long(y0 + Res*(j-1))
;    yCenter[j-1] = long(yCorner[j-1] + Res/2)
;  endfor
;  yCorner[Row] = long(yCorner[Row-1] + Res)

  Res = float(Res)    ; important!!!  long,int
  xCorner = fltArr(Col+1)
  xCenter = fltArr(Col)
  yCorner = fltArr(Row+1)
  yCenter = fltArr(Row)
  for i = 1,Col do begin
    xCorner[i-1] = x0 + Res*(i-1)
    xCenter[i-1] = xCorner[i-1] + Res/2
  endfor
  xCorner[Col] = xCorner[Col-1] + Res
  for j = 1,Row do begin
    yCorner[j-1] = y0 + Res*(j-1)
    yCenter[j-1] = yCorner[j-1] + Res/2
  endfor
  yCorner[Row] = yCorner[Row-1] + Res

;-----------------------------------------------------------------------
; Looping by months
;-----------------------------------------------------------------------
n_orbits = 0
VC_temp = fltarr(Col,Row)
Sam_num_temp = lonarr(Col,Row)

for  JDATE = STDATE,EDDATE do begin
  Year = Fix(JDATE/100)
  Month = JDATE - Year*100L
if ( Month gt 12 ) then begin
  Year = Year + 1
  Month = 1
  JDATE = Year*100L + Month
endif
if ( JDATE le EDDATE) then begin
VC = fltarr(Col,Row)
Area = fltarr(Col,Row)
Sample_num = lonarr(Col,Row)
  Yr4 = string(Year,format='(i4.4)')
  Mon2 = string(Month,format='(i2.2)')
  InDir = DataDir + Yr4 + '/' + Mon2 + '/'
  Out_txt_file = OutDir + 'OMI_NO2_' + Yr4 + Mon2 + '.asc'
  Out_nc_file  = OutDir + 'OMI_NO2_' + Yr4 + Mon2 + '.nc'
  Out_txt_file2 = OutDir + 'OMI_NO2_' + Yr4 + '_' + string(Month_st,format='(i2.2)') + '-' + string(Month_ed,format='(i2.2)') + '.asc'
  Out_nc_file2  = OutDir + 'OMI_NO2_' + Yr4 + '_' + string(Month_st,format='(i2.2)') + '-' + string(Month_ed,format='(i2.2)') + '.nc'
  Out_txt_file3 = OutDir + 'SampleNumber_' + Yr4 + '_' + string(Month_st,format='(i2.2)') + '-' + string(Month_ed,format='(i2.2)') + '.asc'

if ( L_Month eq 'Y' ) then begin
;;;;;  Monthly,InDir,Out_txt_file,Nlp,mapStruct,xCorner,yCorner,n_orbits,Row,Col,x0,y0,Res,VC,Area,Sample_num

;;;;;endfor

;;;;;endtime = systime(/seconds)
;;;;;print,'Elapsed time (seconds) : ', endtime - starttime
;;;;;end

;;;;;pro Monthly,InDir,Out_txt_file,Nlp,mapStruct,xCorner,yCorner,n_orbits,Row,Col,x0,y0,Res,VC,Area,Sample_num

;-----------------------------------------------------------------------
; Looping by HDF files
;-----------------------------------------------------------------------
  LsCmd = 'ls ' + InDir +'*.he5'
  spawn,LsCmd,LsFile
if ( LsFile[0] ne '') then begin
  n_files = n_elements(LsFile)
;;;n_files = 4

  if ( EOS_EXISTS() eq 0 ) then Message, 'HDF not supported'

for  F =0,n_files-1 do begin
  Infile = LsFile[F]
  fid = H5F_OPEN(Infile)
    if ( fid lt 0 ) then Message, 'Error opening file!'
; Read data from HDF file
    LAT_NAME = '/HDFEOS/SWATHS/DominoNO2/Geolocation Fields/Latitude'
    lat_id = H5D_OPEN( fid, LAT_NAME)
    lat = H5D_READ(lat_id)
    H5D_CLOSE, lat_id

    LON_NAME = '/HDFEOS/SWATHS/DominoNO2/Geolocation Fields/Longitude'
    lon_id = H5D_OPEN( fid, LON_NAME)
    lon = H5D_READ(lon_id)
    H5D_CLOSE, lon_id

    SZA_NAME = '/HDFEOS/SWATHS/DominoNO2/Geolocation Fields/SolarZenithAngle'
    sza_id = H5D_OPEN( fid, SZA_NAME)
    sza = H5D_READ(sza_id)
    H5D_CLOSE, sza_id

    TroVCD_NAME = '/HDFEOS/SWATHS/DominoNO2/Data Fields/TroposphericVerticalColumn'
    TroVCD_id = H5D_OPEN( fid, TroVCD_NAME)
    TroVCD = H5D_READ(TroVCD_id)
    H5D_CLOSE, TroVCD_id

    CRD_NAME = '/HDFEOS/SWATHS/DominoNO2/Data Fields/CloudRadianceFraction'
    crd_id = H5D_OPEN( fid, CRD_NAME)
    crd = H5D_READ(crd_id)
    H5D_CLOSE, crd_id

    TIME_NAME = '/HDFEOS/SWATHS/DominoNO2/Geolocation Fields/Time'
    mtime_id = H5D_OPEN( fid, TIME_NAME)
    mtime = H5D_READ(mtime_id)
    H5D_CLOSE, mtime_id

  H5F_CLOSE, fid

; get the date
  JD_ref = JULDAY(1,1,1993,0,0,0)
  JD = JD_ref + mtime[0]/3600/24
  CALDAT,JD,MM,DD,YY,HH,MIN,SEC
  NYMD = YY*10000L + MM*100L + DD
print,NYMD,HH,MIN,SEC

; Remove row anomaly. URL: disc.sci.gsfc.nasa.gov/Aura/data-holdings/OMI/index.shtml#info
; Anomaly 1: Since June 25th, 2007, cross-track scenes 53-54 (0-based)
; Anomaly 2: Since May 11th, 2008, cross-track scenes 37-44 (0-based)
; Anomaly 3: Since January 24th, 2009, cross-track scenes 27-44 (0-based)
; Anomaly 4: Since July 5th, 2011, cross-track scenes 42-45 (0-based)
; Anomaly 5: Since Aug, 2011, cross-track scenes 41-45 (0-based)

  if ( NYMD ge 20070625 ) then begin
    TroVCD[53:54,*] = -999
  if ( NYMD ge 20080511 ) then begin
    TroVCD[37:44,*] = -999
  if ( NYMD ge 20090124 ) then begin
    TroVCD[27:44,*] = -999
  if ( NYMD ge 20110705 ) then begin
    TroVCD[42:45,*] = -999
  if ( NYMD ge 20110801 ) then begin
    TroVCD[41:45,*] = -999
  endif
  endif
  endif
  endif
  endif

; remove large pixels
  PixSize = size(lat,/DIMENSIONS)
  Across_track = PixSize[0]
  lat_new = lat( Nlp:Across_track-Nlp-1 ,*)
  lon_new = lon( Nlp:Across_track-Nlp-1 ,*)
  TroVCD_new = TroVCD( Nlp:Across_track-Nlp-1 ,*)
  sza_new = sza( Nlp:Across_track-Nlp-1 ,*)
  crd_new = crd( Nlp:Across_track-Nlp-1 ,*)

if(L_PROJ eq 'Y') then begin
; transforms map coordinates from longitude/latitude to (x, y) coordinates
  XY = MAP_PROJ_FORWARD( lon_new, lat_new, MAP_STRUCTURE = mapStruct )
endif else begin
  XY2 = [[lon_new(*)],[lat_new(*)]]
  XY = TRANSPOSE(XY2)
endelse

; Filter out solarZenithAngle, cloud fraction etc.
  Ind = where( XY[0,*] gt min(xCorner) and  XY[0,*] lt max(xCorner) and $
               XY[1,*] gt min(yCorner) and  XY[1,*] lt max(yCorner) and $
               TroVCD_new gt 0.0       and  sza_new lt 70.0 and $
               crd_new ge 0            and  crd_new/100.0 lt 50.0   )

  Sum_L = n_elements(Ind)

  if ( Ind[0] ge 0 ) then begin
    Ind_dims = size(lon_new,/DIMENSIONS)
    Ind2 = ARRAY_INDICES(Ind_dims,Ind,/DIMENSIONS)
    n_orbits = n_orbits + 1

  fid = H5F_OPEN(Infile)
    LATCORNER_NAME = '/HDFEOS/SWATHS/DominoNO2/Geolocation Fields/LatitudeCornerpoints'
    latcorner_id = H5D_OPEN( fid, LATCORNER_NAME)
    latcorner = H5D_READ(latcorner_id)
    H5D_CLOSE, latcorner_id

    LONCORNER_NAME = '/HDFEOS/SWATHS/DominoNO2/Geolocation Fields/LongitudeCornerpoints'
    loncorner_id = H5D_OPEN( fid, LONCORNER_NAME)
    loncorner = H5D_READ(loncorner_id)
    H5D_CLOSE, loncorner_id
  H5F_CLOSE, fid

    latcorner_new = latcorner( Nlp:Across_track-Nlp-1,*,*)
    loncorner_new = loncorner( Nlp:Across_track-Nlp-1,*,*)

if(L_PROJ eq 'Y') then begin
;    corner_xy = MAP_PROJ_FORWARD( loncorner_new(Ind2[0,*],Ind2[1,*],*), latcorner_new(Ind2[0,*],Ind2[1,*],*), MAP_STRUCTURE = mapStruct )
    corner_xy = MAP_PROJ_FORWARD( loncorner_new, latcorner_new, MAP_STRUCTURE = mapStruct )
endif else begin
    corner_xy2 = [[loncorner_new(*)],[latcorner_new(*)]]
    corner_xy = TRANSPOSE(corner_xy2)
endelse
help,loncorner_new
help,corner_xy

    n_cxy = n_elements(corner_xy)
    n_cxy2 = n_cxy/2L/4L
    corner_xy1 = REFORM(corner_xy,2,n_cxy2,4)
    corner_x = corner_xy1(0,Ind,*)
    corner_y = corner_xy1(1,Ind,*)

ttt=0
nsample = 100L
    for m=1,Sum_L do begin
      xi1 = max([Long(min((corner_x(0,m-1,*)-x0)/Res)),0])
      yi1 = max([Long(min((corner_y(0,m-1,*)-y0)/Res)),0])
;      xi2 = min([Long(max((corner_x(0,m-1,*)-x0)/Res)),Col-1])
;      yi2 = min([Long(max((corner_y(0,m-1,*)-y0)/Res)),Row-1])
      xi2 = min([Long(max((corner_x(0,m-1,*)-x0)/Res+1)),Col])
      yi2 = min([Long(max((corner_y(0,m-1,*)-y0)/Res+1)),Row])

      if( xi1 lt xi2 and yi1 lt yi2 ) then begin
        lon_xi = [corner_x(0,m-1,0),corner_x(0,m-1,1),corner_x(0,m-1,3),corner_x(0,m-1,2)]
        lat_yi = [corner_y(0,m-1,0),corner_y(0,m-1,1),corner_y(0,m-1,3),corner_y(0,m-1,2)]
        lon_xi2 = long(((lon_xi(*)-x0)/Res-xi1)*nsample)
        lat_yi2 = long(((lat_yi(*)-y0)/Res-yi1)*nsample)

        nx = xi2-xi1
        ny = yi2-yi1
        nxx = long((xi2-xi1)*nsample)
        nyy = long((yi2-yi1)*nsample)
        mark = fltarr(nxx,nyy)
        Ind_grid = POLYFILLV(lon_xi2,lat_yi2,nxx,nyy)

        if(Ind_grid[0] ge 0) then begin
        mark(Ind_grid) = 1.0
        frac = REBIN(mark,nx,ny)

        Area[xi1:xi2-1,yi1:yi2-1] = Area[xi1:xi2-1,yi1:yi2-1] + frac(*,*)
        VC[xi1:xi2-1,yi1:yi2-1] = VC[xi1:xi2-1,yi1:yi2-1] + frac(*,*)*TroVCD_new(Ind2[0,m-1],Ind2[1,m-1])
        Sample_num[xi1:xi2-1,yi1:yi2-1] = Sample_num[xi1:xi2-1,yi1:yi2-1] + 1L

ttt=ttt+1
help,ttt
        endif

      endif

    endfor

;xv=[0,36000,36000,0,0]
;yv=[0,0,36000,36000,0]
;area=poly_area(xv,yv)

;Verts = [[0,0,0],[0,2,0],[2,2,0],[2,0,0]]
;Conn = [4,0,1,2,3]
;Verts1 = [[1,1,0],[1,3,0],[3,3,0],[3,1,0]]
;Conn1 = [4,0,1,2,3]
;bbb=MESH_MERGE(Verts, Conn, Verts1, Conn1 )


  endif


endfor
endif

;-----------------------------------------------------------------------
; Write .asc/.txt output file
;-----------------------------------------------------------------------

  Ind_vc = where( Area gt 0.0)
  VC[Ind_vc] = VC[Ind_vc]/Area[Ind_vc]
  VC = REVERSE(VC, 2)    ; for ArcGIS plot

  openw, lun, Out_txt_file ,/GET_LUN
  printf, lun, 'ncols ', Col
  printf, lun, 'nrows ', Row
  printf, lun, 'xllcorner ', x0
  printf, lun, 'yllcorner ', y0
  printf, lun, 'cellsize ', Res
  printf, lun, 'NODATA_Value ', -999.0
  printf, lun, VC
  close, lun
  Free_LUN, lun

;-----------------------------------------------------------------------
; Write netCDF file
;-----------------------------------------------------------------------

  FileId = NCDF_Create( Out_nc_file, /Clobber )
  NCDF_Control, FileID, /NoFill
    xID   = NCDF_DimDef( FileID, 'X', Col  )
    yID   = NCDF_DimDef( FileID, 'Y', Row  )
    LonID = NCDF_VarDef( FileID, 'LONGITUDE',    [xID],     /Float )
    LatID = NCDF_VarDef( FileID, 'LATITUDE',     [yID],     /Float )
    VCDID = NCDF_VarDef( FileID, 'TropVCD',      [xID,yID], /Float )
    SamID = NCDF_VarDef( FileID, 'SampleNumber', [xID,yID], /Long  )
  NCDF_Attput, FileID, /Global, 'Title', 'NetCDF file created by GAMAP'
  NCDF_Control, FileID, /EnDef
  NCDF_VarPut, FileID, LonID, xCenter,   Count=[ Col ]
  NCDF_VarPut, FileID, LatID, yCenter,   Count=[ Row ]
  NCDF_VarPut, FileID, VCDID, VC,        Count=[ Col, Row ]
  NCDF_VarPut, FileID, SamID, Sample_num,Count=[ Col, Row ]
  NCDF_Close, FileID
endif
if(L_Year eq 'Y') then begin

  if((Month ge Month_st) and (Month le Month_ed)) then begin
    fIdnc1 = NCDF_Open(Out_nc_file) 
    VCDID1 = NCDF_VarId( fIdnc1, 'TropVCD')
    NCDF_VarGet, fIdnc1, VCDID1, VC_temp2
    SamID1 = NCDF_VarId( fIdnc1, 'SampleNumber')
    NCDF_VarGet, fIdnc1, SamID1, Sam_num_temp2
    NCDF_Close, fIdnc1

    VC_temp = VC_temp + VC_temp2 * Sam_num_temp2
    Sam_num_temp = Sam_num_temp + Sam_num_temp2

  if(Month eq Month_ed) then begin
    Ind_vc2 = where( Sam_num_temp gt 0)
    VC_temp[Ind_vc2] = VC_temp[Ind_vc2]/Sam_num_temp[Ind_vc2]

  openw, lun2, Out_txt_file2 ,/GET_LUN
  printf, lun2, 'ncols ', Col
  printf, lun2, 'nrows ', Row
  printf, lun2, 'xllcorner ', x0
  printf, lun2, 'yllcorner ', y0
  printf, lun2, 'cellsize ', Res
  printf, lun2, 'NODATA_Value ', -999.0
  printf, lun2, VC_temp
  close, lun2
  Free_LUN, lun2

  openw, lun3, Out_txt_file3 ,/GET_LUN
  printf, lun3, 'ncols ', Col
  printf, lun3, 'nrows ', Row
  printf, lun3, 'xllcorner ', x0
  printf, lun3, 'yllcorner ', y0
  printf, lun3, 'cellsize ', Res
  printf, lun3, 'NODATA_Value ', -999.0
  printf, lun3, Sam_num_temp
  close, lun3
  Free_LUN, lun3

    VC_temp(*) = 0.0
    Sam_num_temp(*) = 0

  endif

  endif

endif

endif
endfor

endtime = systime(/seconds)
print,'Elapsed time (seconds) : ', endtime - starttime

end
