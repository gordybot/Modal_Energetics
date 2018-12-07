#!/bin/csh -x
#PBS -N r256_260_1z_archv.2016
#PBS -j oe
#PBS -o r256_260_1z_archv.2016.log
#PBS -l select=5:ncpus=24:mpiprocs=24
#PBS -l walltime=20:00:00
#PBS -W umask=027
#PBS -A ONRDC38645449
#PBS -q standard  
#
#!/bin/csh
#
set echo
set time=1
#
# --- interpolate to 3-d z-levels from a HYCOM mean archive file.
# --- configured for K input layers and a small number of output levels.
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM).
# --- or use archv2ncdf3z for netCDF output.
#
# --- output is HYCOM .a files.
#
setenv R GLBc0.04
setenv K 41
setenv Y 2016
#setenv E 1111 
#
#setenv X `echo ${E} | awk '{printf("%04.1f", $1*0.1)}'`
setenv W ~wallcraf/hycom
#
setenv NUM 221
setenv EXPT expt_22.1
setenv DD regional.depth
setenv GG regional.grid 
setenv SD $WORKDIR/hycom/${R}/${EXPT}/data/input3zc/
#
#start day and end day
setenv DF 259   #1
setenv DL 260   #32
#
# interpolation 0,1,2 and z,zi
setenv IN 1z
#
cd ${SD}
#
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/cp -u ${SD}${DD}.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/cp -u ${SD}${DD}.b regional.depth.b
endif
# 
# copy something that is in the same directory is silly
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/cp -u ${SD}${GG}.a regional.grid.a
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/cp -u ${SD}${GG}.b regional.grid.b
endif
#
#
# loop over days and hours
#
 @ D = ${DF}
 while ($D <= ${DL})

##Turn eg 1 to 001 in the next line
#
setenv DA `echo $D | awk '{printf("%03d", $1)}'`

foreach h ( 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23)
#foreach h ( 00 01 )

if ( -e ${NUM}_archv.${Y}_${DA}_${h}.a ) then #exist?
echo ${NUM}_archv.${Y}_${DA}_${h}.a

#
# --- y,m select the archive files.
#
setenv FOR042A archv.${Y}_${IN}_${DA}_${h}_T.a
setenv FOR042  archv.${Y}_${IN}_${DA}_${h}_T.b
setenv FOR043A archv.${Y}_${IN}_${DA}_${h}_S.a
setenv FOR043  archv.${Y}_${IN}_${DA}_${h}_S.b
setenv FOR044A archv.${Y}_${IN}_${DA}_${h}_r.a
setenv FOR044  archv.${Y}_${IN}_${DA}_${h}_r.b
setenv FOR045A archv.${Y}_${IN}_${DA}_${h}_i.a
setenv FOR045  archv.${Y}_${IN}_${DA}_${h}_i.b
/bin/rm $FOR042  $FOR043  $FOR044  $FOR045
/bin/rm $FOR042A $FOR043A $FOR044A $FOR045A
#cat <<E-o-D >! archv2data3z.IN
cat <<E-o-D >! blkdat.input
${NUM}_archv.${Y}_${DA}_${h}.a
HYCOM
 000	'iexpt ' = experiment number x10 (000=from archive file)
   3	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
9000	'idm   ' = longitudinal array size
7055	'jdm   ' = latitudinal  array size
${K}	'kdm   ' = number of layers
   0.0	'thbase' = reference density (sigma units)
   0	'smooth' = smooth fields before plotting (0=F,1=T)
   1	'baclin' = extract baroclinic velocity (0=total,1=baroclinic)
   1	'xyward' = output original unrotated velocities (0=no,1=yes)
   1	'iorign' = i-origin of plotted subregion
   1	'jorign' = j-origin of plotted subregion
   0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
   0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
   2	'inbot ' = read in bot and/or zbot (1=bot,2=zbot,3=bot&zbot)
   0.0	'zbot  ' = depth above bottom for bottom value (none for zbot<=0.0)
   1	'itype ' = interpolation type (0=sample,1=linear,2=parabolic)
 240	'kz    ' = number of depths to sample
  12.5  'z     ' = sample depth   1
  37.5  'z     ' = sample depth   2
  62.5  'z     ' = sample depth   3
  87.5  'z     ' = sample depth   4
 112.5  'z     ' = sample depth   5
 137.5  'z     ' = sample depth   6
 162.5  'z     ' = sample depth   7
 187.5  'z     ' = sample depth   8
 212.5  'z     ' = sample depth   9
 237.5  'z     ' = sample depth  10
 262.5  'z     ' = sample depth  11
 287.5  'z     ' = sample depth  12
 312.5  'z     ' = sample depth  13
 337.5  'z     ' = sample depth  14
 362.5  'z     ' = sample depth  15
 387.5  'z     ' = sample depth  16
 412.5  'z     ' = sample depth  17
 437.5  'z     ' = sample depth  18
 462.5  'z     ' = sample depth  19
 487.5  'z     ' = sample depth  20
 512.5  'z     ' = sample depth  21
 537.5  'z     ' = sample depth  22
 562.5  'z     ' = sample depth  23
 587.5  'z     ' = sample depth  24
 612.5  'z     ' = sample depth  25
 637.5  'z     ' = sample depth  26
 662.5  'z     ' = sample depth  27
 687.5  'z     ' = sample depth  28
 712.5  'z     ' = sample depth  29
 737.5  'z     ' = sample depth  30
 762.5  'z     ' = sample depth  31
 787.5  'z     ' = sample depth  32
 812.5  'z     ' = sample depth  33
 837.5  'z     ' = sample depth  34
 862.5  'z     ' = sample depth  35
 887.5  'z     ' = sample depth  36
 912.5  'z     ' = sample depth  37
 937.5  'z     ' = sample depth  38
 962.5  'z     ' = sample depth  39
 987.5  'z     ' = sample depth  40
1012.5  'z     ' = sample depth  41
1037.5  'z     ' = sample depth  42
1062.5  'z     ' = sample depth  43
1087.5  'z     ' = sample depth  44
1112.5  'z     ' = sample depth  45
1137.5  'z     ' = sample depth  46
1162.5  'z     ' = sample depth  47
1187.5  'z     ' = sample depth  48
1212.5  'z     ' = sample depth  49
1237.5  'z     ' = sample depth  50
1262.5  'z     ' = sample depth  51
1287.5  'z     ' = sample depth  52
1312.5  'z     ' = sample depth  53
1337.5  'z     ' = sample depth  54
1362.5  'z     ' = sample depth  55
1387.5  'z     ' = sample depth  56
1412.5  'z     ' = sample depth  57
1437.5  'z     ' = sample depth  58
1462.5  'z     ' = sample depth  59
1487.5  'z     ' = sample depth  60
1512.5  'z     ' = sample depth  61
1537.5  'z     ' = sample depth  62
1562.5  'z     ' = sample depth  63
1587.5  'z     ' = sample depth  64
1612.5  'z     ' = sample depth  65
1637.5  'z     ' = sample depth  66
1662.5  'z     ' = sample depth  67
1687.5  'z     ' = sample depth  68
1712.5  'z     ' = sample depth  69
1737.5  'z     ' = sample depth  70
1762.5  'z     ' = sample depth  71
1787.5  'z     ' = sample depth  72
1812.5  'z     ' = sample depth  73
1837.5  'z     ' = sample depth  74
1862.5  'z     ' = sample depth  75
1887.5  'z     ' = sample depth  76
1912.5  'z     ' = sample depth  77
1937.5  'z     ' = sample depth  78
1962.5  'z     ' = sample depth  79
1987.5  'z     ' = sample depth  80
2012.5  'z     ' = sample depth  81
2037.5  'z     ' = sample depth  82
2062.5  'z     ' = sample depth  83
2087.5  'z     ' = sample depth  84
2112.5  'z     ' = sample depth  85
2137.5  'z     ' = sample depth  86
2162.5  'z     ' = sample depth  87
2187.5  'z     ' = sample depth  88
2212.5  'z     ' = sample depth  89
2237.5  'z     ' = sample depth  90
2262.5  'z     ' = sample depth  91
2287.5  'z     ' = sample depth  92
2312.5  'z     ' = sample depth  93
2337.5  'z     ' = sample depth  94
2362.5  'z     ' = sample depth  95
2387.5  'z     ' = sample depth  96
2412.5  'z     ' = sample depth  97
2437.5  'z     ' = sample depth  98
2462.5  'z     ' = sample depth  99
2487.5  'z     ' = sample depth 100
2512.5  'z     ' = sample depth 101
2537.5  'z     ' = sample depth 102
2562.5  'z     ' = sample depth 103
2587.5  'z     ' = sample depth 104
2612.5  'z     ' = sample depth 105
2637.5  'z     ' = sample depth 106
2662.5  'z     ' = sample depth 107
2687.5  'z     ' = sample depth 108
2712.5  'z     ' = sample depth 109
2737.5  'z     ' = sample depth 110
2762.5  'z     ' = sample depth 111
2787.5  'z     ' = sample depth 112
2812.5  'z     ' = sample depth 113
2837.5  'z     ' = sample depth 114
2862.5  'z     ' = sample depth 115
2887.5  'z     ' = sample depth 116
2912.5  'z     ' = sample depth 117
2937.5  'z     ' = sample depth 118
2962.5  'z     ' = sample depth 119
2987.5  'z     ' = sample depth 120
3012.5  'z     ' = sample depth 121
3037.5  'z     ' = sample depth 122
3062.5  'z     ' = sample depth 123
3087.5  'z     ' = sample depth 124
3112.5  'z     ' = sample depth 125
3137.5  'z     ' = sample depth 126
3162.5  'z     ' = sample depth 127
3187.5  'z     ' = sample depth 128
3212.5  'z     ' = sample depth 129
3237.5  'z     ' = sample depth 130
3262.5  'z     ' = sample depth 131
3287.5  'z     ' = sample depth 132
3312.5  'z     ' = sample depth 133
3337.5  'z     ' = sample depth 134
3362.5  'z     ' = sample depth 135
3387.5  'z     ' = sample depth 136
3412.5  'z     ' = sample depth 137
3437.5  'z     ' = sample depth 138
3462.5  'z     ' = sample depth 139
3487.5  'z     ' = sample depth 140
3512.5  'z     ' = sample depth 141
3537.5  'z     ' = sample depth 142
3562.5  'z     ' = sample depth 143
3587.5  'z     ' = sample depth 144
3612.5  'z     ' = sample depth 145
3637.5  'z     ' = sample depth 146
3662.5  'z     ' = sample depth 147
3687.5  'z     ' = sample depth 148
3712.5  'z     ' = sample depth 149
3737.5  'z     ' = sample depth 150
3762.5  'z     ' = sample depth 151
3787.5  'z     ' = sample depth 152
3812.5  'z     ' = sample depth 153
3837.5  'z     ' = sample depth 154
3862.5  'z     ' = sample depth 155
3887.5  'z     ' = sample depth 156
3912.5  'z     ' = sample depth 157
3937.5  'z     ' = sample depth 158
3962.5  'z     ' = sample depth 159
3987.5  'z     ' = sample depth 160
4012.5  'z     ' = sample depth 161
4037.5  'z     ' = sample depth 162
4062.5  'z     ' = sample depth 163
4087.5  'z     ' = sample depth 164
4112.5  'z     ' = sample depth 165
4137.5  'z     ' = sample depth 166
4162.5  'z     ' = sample depth 167
4187.5  'z     ' = sample depth 168
4212.5  'z     ' = sample depth 169
4237.5  'z     ' = sample depth 170
4262.5  'z     ' = sample depth 171
4287.5  'z     ' = sample depth 172
4312.5  'z     ' = sample depth 173
4337.5  'z     ' = sample depth 174
4362.5  'z     ' = sample depth 175
4387.5  'z     ' = sample depth 176
4412.5  'z     ' = sample depth 177
4437.5  'z     ' = sample depth 178
4462.5  'z     ' = sample depth 179
4487.5  'z     ' = sample depth 180
4512.5  'z     ' = sample depth 181
4537.5  'z     ' = sample depth 182
4562.5  'z     ' = sample depth 183
4587.5  'z     ' = sample depth 184
4612.5  'z     ' = sample depth 185
4637.5  'z     ' = sample depth 186
4662.5  'z     ' = sample depth 187
4687.5  'z     ' = sample depth 188
4712.5  'z     ' = sample depth 189
4737.5  'z     ' = sample depth 190
4762.5  'z     ' = sample depth 191
4787.5  'z     ' = sample depth 192
4812.5  'z     ' = sample depth 193
4837.5  'z     ' = sample depth 194
4862.5  'z     ' = sample depth 195
4887.5  'z     ' = sample depth 196
4912.5  'z     ' = sample depth 197
4937.5  'z     ' = sample depth 198
4962.5  'z     ' = sample depth 199
4987.5  'z     ' = sample depth 200
5012.5  'z     ' = sample depth 201
5037.5  'z     ' = sample depth 202
5062.5  'z     ' = sample depth 203
5087.5  'z     ' = sample depth 204
5112.5  'z     ' = sample depth 205
5137.5  'z     ' = sample depth 206
5162.5  'z     ' = sample depth 207
5187.5  'z     ' = sample depth 208
5212.5  'z     ' = sample depth 209
5237.5  'z     ' = sample depth 210
5262.5  'z     ' = sample depth 211
5287.5  'z     ' = sample depth 212
5312.5  'z     ' = sample depth 213
5337.5  'z     ' = sample depth 214
5362.5  'z     ' = sample depth 215
5387.5  'z     ' = sample depth 216
5412.5  'z     ' = sample depth 217
5437.5  'z     ' = sample depth 218
5462.5  'z     ' = sample depth 219
5487.5  'z     ' = sample depth 220
5512.5  'z     ' = sample depth 221
5537.5  'z     ' = sample depth 222
5562.5  'z     ' = sample depth 223
5587.5  'z     ' = sample depth 224
5612.5  'z     ' = sample depth 225
5637.5  'z     ' = sample depth 226
5662.5  'z     ' = sample depth 227
5687.5  'z     ' = sample depth 228
5712.5  'z     ' = sample depth 229
5737.5  'z     ' = sample depth 230
5762.5  'z     ' = sample depth 231
5787.5  'z     ' = sample depth 232
5812.5  'z     ' = sample depth 233
5837.5  'z     ' = sample depth 234
5862.5  'z     ' = sample depth 235
5887.5  'z     ' = sample depth 236
5912.5  'z     ' = sample depth 237
5937.5  'z     ' = sample depth 238
5962.5  'z     ' = sample depth 239
5987.5  'z     ' = sample depth 240
   0	'botio ' = bathymetry  I/O unit (0 no I/O)
   0	'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
   0.25	'tmljmq' = equiv. temp. jump across mixed-layer (degC,  0 no I/O)
   0.25 'tmlnav' = NAVO rp33 T  jump across mixed-layer (degC,  0 no I/O)
   0	'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
   0	'wvlio ' = w-velocity  I/O unit (0 no I/O)
   0	'uvlio ' = u-velocity  I/O unit (0 no I/O)
   0	'vvlio ' = v-velocity  I/O unit (0 no I/O)
   0	'splio ' = speed       I/O unit (0 no I/O)
  42	'temio ' = temperature I/O unit (0 no I/O)
  43	'salio ' = salinity    I/O unit (0 no I/O)
   0	'tthio ' = density     I/O unit (0 no I/O)
   0	'keio  ' = kinetic egy I/O unit (0 no I/O)
E-o-D
#cat archv2data3z.IN
cat blkdat.input
unset echo
#if (-e ${W}/ALLcnl/archive/src_2.2.35/archv2data3z) then
#   aprun -n 1 -m 31g ${W}/ALLcnl/archive/src_2.2.35/archv2data3z < archv2data3z.IN
#   aprun -n 1 -m 120g ${W}/ALLcnl/archive/src_2.2.35/archv2data3z < archv2data3z.IN
#   aprun -n 60 -N 12 -S 6 ${W}/ALLmpi/archive/src/archv2data3z < archv2data3z.IN
#   aprun -n 5 -N 1  ${W}/ALLmpi/archive/src/archv2data3z < archv2data3z.IN 
#   aprun -n 5 -N 1  ${W}/ALLmpi/archive/src/archv2data3z
   aprun -n 60 -N 12 -S 6  ${W}/ALLmpi/archive/src/archv2data3z
#else
#   ${W}/ALL/archive/src_2.2.35/archv2data3z < archv2data3z.IN
#endif

#
endif #if exist

end #h
 @ D = $D + 1
#
end #while

