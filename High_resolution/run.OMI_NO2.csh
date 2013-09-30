#! /bin/csh -f
########################################################################
#set pp  = ( NMG_HNYM SD_DTHD GS_LZXG GS_HNPL SX_GHJJ SC_SCBS HN_DTXY LN_LNHD JL_CCRD )
#set xpp = ( 119.7733333 120.218745 103.625 106.7884889 110.165717 104.772803 114.141098 123.802634 125.142073 )
#set ypp = ( 48.54722222 36.040682 36.12916667 35.5017 38.736311 31.801032 32.108536 42.344963 43.778622 )
#set pp  = ( HN_HNHY HN_DTXT HN_GDYY GD_YXHB GD_HNST ZJ_GDBL SD_HDQD SD_HNWH NMG_HNYM )
#set xpp = ( 111.486254 113.002807 112.310953 111.669332 116.735822 121.4833333 120.33639 122.204958 119.7341667 )
#set ypp = ( 27.622488 27.824939 28.5935 21.53923 23.326046 29.93333333 36.119446 37.449551 49.23777778 )
#set pp  = ( SCLZ HNHNYY HNJGS HNYH HNWH )
#set xpp = ( 105.288548 113.163922 115.019557 121.148958 122.204958 )
#set ypp = ( 28.772136 29.444989 27.047897 28.115131 37.449551 )
set out = '/home/liufei/Data/Large_point/Regrid/results'
set pp  = ( chengdu P2 )
set xpp = ( 104.062347 121.6)
set ypp = ( 30.665872 31.3528)
set s10 = "OutDir = '/home/liufei/Data/Large_point/Regrid/results/'"
set s20 = "xpp = 119.7733333"
set s30 = "ypp = 48.54722222"

set num = $#pp
@ i = 1
while ($i <= $num)
  set OutDir = $out/${pp[$i]}/
  if (! -d $OutDir) then
    mkdir -p "$OutDir"
    echo "mkdir: $OutDir"
  endif

  set s1  = "OutDir = '$OutDir'"
  set s2  = "xpp = ${xpp[$i]}"
  set s3  = "ypp = ${ypp[$i]}"
  echo $s1
  echo $s2
  echo $s3
  sed -e '1s/OMI_NO2/OMI_NO2_temp/g' OMI_NO2.pro > OMI_NO2_temp.pro
  sed -i "s:${s10}:${s1}:g" OMI_NO2_temp.pro
  sed -i "s:${s20}:${s2}:g" OMI_NO2_temp.pro
  sed -i "s:${s30}:${s3}:g" OMI_NO2_temp.pro
idl < run.idl
@ i++
end
