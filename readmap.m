clear;

data=load('/home/liufei/Data/Large_point/Steffen/wind_dependent_mean_v20120905_Far_East.mat');
TVCD=data.map_TVCD;
location=load('/home/liufei/Data/Large_point/Steffen/lat_lon_grid_Far_East.mat');
longs=location.longs;
lats=location.lats;

long_min=112;
long_max=113;
lat_min=35;
lat_max=36;
index_longs=find(longs >long_min & longs <long_max);
index_lats=find(lats>lat_min & lats<lat_max);
cellsize = 0.176;
ncols=size(index_longs,2)
nrows=size(index_lats,2)

wind1=TVCD(index_lats,index_longs,3,1,2);
wind2=TVCD(index_lats,index_longs,3,2,2);
wind3=TVCD(index_lats,index_longs,3,3,2);
wind4=TVCD(index_lats,index_longs,3,4,2);
wind5=TVCD(index_lats,index_longs,3,5,2);
wind6=TVCD(index_lats,index_longs,3,6,2);
wind7=TVCD(index_lats,index_longs,3,7,2);
wind8=TVCD(index_lats,index_longs,3,8,2);
wind9=TVCD(index_lats,index_longs,3,9,2);


for i = 1:9
    resFile = ['wind_' num2str(i) '.asc'];
    fid = fopen(resFile,'w');
    fprintf(fid,['ncols ' num2str(ncols) '\r\nnrows ' num2str(nrows) '\r\nxllcorner '  num2str(longs(index_longs(1))) '\r\nyllcorner '  num2str(lats(index_lats(nrows))) '\r\ncellsize ' num2str(cellsize) '\r\nNODATA_value  NaN\r\n']);
    fclose(fid);
    eval(['dlmwrite(resFile,' ['wind' num2str(i)] ',''-append'',''delimiter'','' '',''precision'',6);']);
end
