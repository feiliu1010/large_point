clear;

%ncid=netcdf.create('map_TVCD.nc','NC_NOCLOBBER');
%ncid=netcdf.create('map_TVCDstd.nc','NC_NOCLOBBER');
ncid=netcdf.create('map_U.nc','NC_NOCLOBBER');
%ncid=netcdf.create('latitude.nc','NC_NOCLOBBER');
%ncid=netcdf.create('longitude.nc','NC_NOCLOBBER');
dimid1=netcdf.defDim(ncid,'lat',160);
dimid2=netcdf.defDim(ncid,'lon',260);
dimid3=netcdf.defDim(ncid,'sea',4);
dimid4=netcdf.defDim(ncid,'windd',9);
dimid5=netcdf.defDim(ncid,'winds',3);
%varid=netcdf.defVar(ncid,'map_TVCD','double',[dimid1,dimid2,dimid3,dimid4,dimid5]);
%varid=netcdf.defVar(ncid,'map_TVCDstd','double',[dimid1,dimid2,dimid3,dimid4,dimid5]);
varid=netcdf.defVar(ncid,'map_U','double',[dimid1,dimid2,dimid3,dimid4,dimid5]);
%varid=netcdf.defVar(ncid,'map_TVCD','double',dimid2);
netcdf.endDef(ncid);

target=load('/home/liufei/Data/Large_point/Steffen/wind_dependent_mean_v20120905_Far_East.mat');
%data=target.map_TVCD;
%data=target.map_TVCDstd;
data=target.map_TVCD;


%target=load('/home/liufei/Data/Large_point/Steffen/lat_lon_grid_Far_East.mat');
%data=target.lats;
%data=target.longs;
netcdf.putVar(ncid,varid,data);

netcdf.close(ncid);


