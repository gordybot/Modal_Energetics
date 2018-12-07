#!/bin/csh -x
#PBS -N r245_250_2zi_archv.2016
#PBS -j oe
#PBS -o r245_250_2zi_archv.2016.log
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
setenv SD $WORKDIR/hycom/${R}/${EXPT}/data/input3z/
#
#start day and end day
setenv DF 245   #1
setenv DL 250   #32
#
# interpolation 0,1,2 and z,zi
setenv IN 2zi
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
#foreach h ( 00 01)

if ( -e ${NUM}_archv.${Y}_${DA}_${h}.a ) then #exist?
echo ${NUM}_archv.${Y}_${DA}_${h}.a

#
# --- y,m select the archive files.
#
setenv FOR042A archv.${Y}_${IN}_${DA}_${h}_u.a
setenv FOR042  archv.${Y}_${IN}_${DA}_${h}_u.b
setenv FOR043A archv.${Y}_${IN}_${DA}_${h}_v.a
setenv FOR043  archv.${Y}_${IN}_${DA}_${h}_v.b
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
   2	'itype ' = interpolation type (0=sample,1=linear,2=parabolic)
 240	'kzi   ' = number of depths to sample
   0.0  'zi    ' = sample cell interface    1
  25.0  'zi    ' = sample cell interface    2
  50.0  'zi    ' = sample cell interface    3
  75.0  'zi    ' = sample cell interface    4
 100.0  'zi    ' = sample cell interface    5
 125.0  'zi    ' = sample cell interface    6
 150.0  'zi    ' = sample cell interface    7
 175.0  'zi    ' = sample cell interface    8
 200.0  'zi    ' = sample cell interface    9
 225.0  'zi    ' = sample cell interface   10
 250.0  'zi    ' = sample cell interface   11
 275.0  'zi    ' = sample cell interface   12
 300.0  'zi    ' = sample cell interface   13
 325.0  'zi    ' = sample cell interface   14
 350.0  'zi    ' = sample cell interface   15
 375.0  'zi    ' = sample cell interface   16
 400.0  'zi    ' = sample cell interface   17
 425.0  'zi    ' = sample cell interface   18
 450.0  'zi    ' = sample cell interface   19
 475.0  'zi    ' = sample cell interface   20
 500.0  'zi    ' = sample cell interface   21
 525.0  'zi    ' = sample cell interface   22
 550.0  'zi    ' = sample cell interface   23
 575.0  'zi    ' = sample cell interface   24
 600.0  'zi    ' = sample cell interface   25
 625.0  'zi    ' = sample cell interface   26
 650.0  'zi    ' = sample cell interface   27
 675.0  'zi    ' = sample cell interface   28
 700.0  'zi    ' = sample cell interface   29
 725.0  'zi    ' = sample cell interface   30
 750.0  'zi    ' = sample cell interface   31
 775.0  'zi    ' = sample cell interface   32
 800.0  'zi    ' = sample cell interface   33
 825.0  'zi    ' = sample cell interface   34
 850.0  'zi    ' = sample cell interface   35
 875.0  'zi    ' = sample cell interface   36
 900.0  'zi    ' = sample cell interface   37
 925.0  'zi    ' = sample cell interface   38
 950.0  'zi    ' = sample cell interface   39
 975.0  'zi    ' = sample cell interface   40
1000.0  'zi    ' = sample cell interface   41
1025.0  'zi    ' = sample cell interface   42
1050.0  'zi    ' = sample cell interface   43
1075.0  'zi    ' = sample cell interface   44
1100.0  'zi    ' = sample cell interface   45
1125.0  'zi    ' = sample cell interface   46
1150.0  'zi    ' = sample cell interface   47
1175.0  'zi    ' = sample cell interface   48
1200.0  'zi    ' = sample cell interface   49
1225.0  'zi    ' = sample cell interface   50
1250.0  'zi    ' = sample cell interface   51
1275.0  'zi    ' = sample cell interface   52
1300.0  'zi    ' = sample cell interface   53
1325.0  'zi    ' = sample cell interface   54
1350.0  'zi    ' = sample cell interface   55
1375.0  'zi    ' = sample cell interface   56
1400.0  'zi    ' = sample cell interface   57
1425.0  'zi    ' = sample cell interface   58
1450.0  'zi    ' = sample cell interface   59
1475.0  'zi    ' = sample cell interface   60
1500.0  'zi    ' = sample cell interface   61
1525.0  'zi    ' = sample cell interface   62
1550.0  'zi    ' = sample cell interface   63
1575.0  'zi    ' = sample cell interface   64
1600.0  'zi    ' = sample cell interface   65
1625.0  'zi    ' = sample cell interface   66
1650.0  'zi    ' = sample cell interface   67
1675.0  'zi    ' = sample cell interface   68
1700.0  'zi    ' = sample cell interface   69
1725.0  'zi    ' = sample cell interface   70
1750.0  'zi    ' = sample cell interface   71
1775.0  'zi    ' = sample cell interface   72
1800.0  'zi    ' = sample cell interface   73
1825.0  'zi    ' = sample cell interface   74
1850.0  'zi    ' = sample cell interface   75
1875.0  'zi    ' = sample cell interface   76
1900.0  'zi    ' = sample cell interface   77
1925.0  'zi    ' = sample cell interface   78
1950.0  'zi    ' = sample cell interface   79
1975.0  'zi    ' = sample cell interface   80
2000.0  'zi    ' = sample cell interface   81
2025.0  'zi    ' = sample cell interface   82
2050.0  'zi    ' = sample cell interface   83
2075.0  'zi    ' = sample cell interface   84
2100.0  'zi    ' = sample cell interface   85
2125.0  'zi    ' = sample cell interface   86
2150.0  'zi    ' = sample cell interface   87
2175.0  'zi    ' = sample cell interface   88
2200.0  'zi    ' = sample cell interface   89
2225.0  'zi    ' = sample cell interface   90
2250.0  'zi    ' = sample cell interface   91
2275.0  'zi    ' = sample cell interface   92
2300.0  'zi    ' = sample cell interface   93
2325.0  'zi    ' = sample cell interface   94
2350.0  'zi    ' = sample cell interface   95
2375.0  'zi    ' = sample cell interface   96
2400.0  'zi    ' = sample cell interface   97
2425.0  'zi    ' = sample cell interface   98
2450.0  'zi    ' = sample cell interface   99
2475.0  'zi    ' = sample cell interface  100
2500.0  'zi    ' = sample cell interface  101
2525.0  'zi    ' = sample cell interface  102
2550.0  'zi    ' = sample cell interface  103
2575.0  'zi    ' = sample cell interface  104
2600.0  'zi    ' = sample cell interface  105
2625.0  'zi    ' = sample cell interface  106
2650.0  'zi    ' = sample cell interface  107
2675.0  'zi    ' = sample cell interface  108
2700.0  'zi    ' = sample cell interface  109
2725.0  'zi    ' = sample cell interface  110
2750.0  'zi    ' = sample cell interface  111
2775.0  'zi    ' = sample cell interface  112
2800.0  'zi    ' = sample cell interface  113
2825.0  'zi    ' = sample cell interface  114
2850.0  'zi    ' = sample cell interface  115
2875.0  'zi    ' = sample cell interface  116
2900.0  'zi    ' = sample cell interface  117
2925.0  'zi    ' = sample cell interface  118
2950.0  'zi    ' = sample cell interface  119
2975.0  'zi    ' = sample cell interface  120
3000.0  'zi    ' = sample cell interface  121
3025.0  'zi    ' = sample cell interface  122
3050.0  'zi    ' = sample cell interface  123
3075.0  'zi    ' = sample cell interface  124
3100.0  'zi    ' = sample cell interface  125
3125.0  'zi    ' = sample cell interface  126
3150.0  'zi    ' = sample cell interface  127
3175.0  'zi    ' = sample cell interface  128
3200.0  'zi    ' = sample cell interface  129
3225.0  'zi    ' = sample cell interface  130
3250.0  'zi    ' = sample cell interface  131
3275.0  'zi    ' = sample cell interface  132
3300.0  'zi    ' = sample cell interface  133
3325.0  'zi    ' = sample cell interface  134
3350.0  'zi    ' = sample cell interface  135
3375.0  'zi    ' = sample cell interface  136
3400.0  'zi    ' = sample cell interface  137
3425.0  'zi    ' = sample cell interface  138
3450.0  'zi    ' = sample cell interface  139
3475.0  'zi    ' = sample cell interface  140
3500.0  'zi    ' = sample cell interface  141
3525.0  'zi    ' = sample cell interface  142
3550.0  'zi    ' = sample cell interface  143
3575.0  'zi    ' = sample cell interface  144
3600.0  'zi    ' = sample cell interface  145
3625.0  'zi    ' = sample cell interface  146
3650.0  'zi    ' = sample cell interface  147
3675.0  'zi    ' = sample cell interface  148
3700.0  'zi    ' = sample cell interface  149
3725.0  'zi    ' = sample cell interface  150
3750.0  'zi    ' = sample cell interface  151
3775.0  'zi    ' = sample cell interface  152
3800.0  'zi    ' = sample cell interface  153
3825.0  'zi    ' = sample cell interface  154
3850.0  'zi    ' = sample cell interface  155
3875.0  'zi    ' = sample cell interface  156
3900.0  'zi    ' = sample cell interface  157
3925.0  'zi    ' = sample cell interface  158
3950.0  'zi    ' = sample cell interface  159
3975.0  'zi    ' = sample cell interface  160
4000.0  'zi    ' = sample cell interface  161
4025.0  'zi    ' = sample cell interface  162
4050.0  'zi    ' = sample cell interface  163
4075.0  'zi    ' = sample cell interface  164
4100.0  'zi    ' = sample cell interface  165
4125.0  'zi    ' = sample cell interface  166
4150.0  'zi    ' = sample cell interface  167
4175.0  'zi    ' = sample cell interface  168
4200.0  'zi    ' = sample cell interface  169
4225.0  'zi    ' = sample cell interface  170
4250.0  'zi    ' = sample cell interface  171
4275.0  'zi    ' = sample cell interface  172
4300.0  'zi    ' = sample cell interface  173
4325.0  'zi    ' = sample cell interface  174
4350.0  'zi    ' = sample cell interface  175
4375.0  'zi    ' = sample cell interface  176
4400.0  'zi    ' = sample cell interface  177
4425.0  'zi    ' = sample cell interface  178
4450.0  'zi    ' = sample cell interface  179
4475.0  'zi    ' = sample cell interface  180
4500.0  'zi    ' = sample cell interface  181
4525.0  'zi    ' = sample cell interface  182
4550.0  'zi    ' = sample cell interface  183
4575.0  'zi    ' = sample cell interface  184
4600.0  'zi    ' = sample cell interface  185
4625.0  'zi    ' = sample cell interface  186
4650.0  'zi    ' = sample cell interface  187
4675.0  'zi    ' = sample cell interface  188
4700.0  'zi    ' = sample cell interface  189
4725.0  'zi    ' = sample cell interface  190
4750.0  'zi    ' = sample cell interface  191
4775.0  'zi    ' = sample cell interface  192
4800.0  'zi    ' = sample cell interface  193
4825.0  'zi    ' = sample cell interface  194
4850.0  'zi    ' = sample cell interface  195
4875.0  'zi    ' = sample cell interface  196
4900.0  'zi    ' = sample cell interface  197
4925.0  'zi    ' = sample cell interface  198
4950.0  'zi    ' = sample cell interface  199
4975.0  'zi    ' = sample cell interface  200
5000.0  'zi    ' = sample cell interface  201
5025.0  'zi    ' = sample cell interface  202
5050.0  'zi    ' = sample cell interface  203
5075.0  'zi    ' = sample cell interface  204
5100.0  'zi    ' = sample cell interface  205
5125.0  'zi    ' = sample cell interface  206
5150.0  'zi    ' = sample cell interface  207
5175.0  'zi    ' = sample cell interface  208
5200.0  'zi    ' = sample cell interface  209
5225.0  'zi    ' = sample cell interface  210
5250.0  'zi    ' = sample cell interface  211
5275.0  'zi    ' = sample cell interface  212
5300.0  'zi    ' = sample cell interface  213
5325.0  'zi    ' = sample cell interface  214
5350.0  'zi    ' = sample cell interface  215
5375.0  'zi    ' = sample cell interface  216
5400.0  'zi    ' = sample cell interface  217
5425.0  'zi    ' = sample cell interface  218
5450.0  'zi    ' = sample cell interface  219
5475.0  'zi    ' = sample cell interface  220
5500.0  'zi    ' = sample cell interface  221
5525.0  'zi    ' = sample cell interface  222
5550.0  'zi    ' = sample cell interface  223
5575.0  'zi    ' = sample cell interface  224
5600.0  'zi    ' = sample cell interface  225
5625.0  'zi    ' = sample cell interface  226
5650.0  'zi    ' = sample cell interface  227
5675.0  'zi    ' = sample cell interface  228
5700.0  'zi    ' = sample cell interface  229
5725.0  'zi    ' = sample cell interface  230
5750.0  'zi    ' = sample cell interface  231
5775.0  'zi    ' = sample cell interface  232
5800.0  'zi    ' = sample cell interface  233
5825.0  'zi    ' = sample cell interface  234
5850.0  'zi    ' = sample cell interface  235
5875.0  'zi    ' = sample cell interface  236
5900.0  'zi    ' = sample cell interface  237
5925.0  'zi    ' = sample cell interface  238
5950.0  'zi    ' = sample cell interface  239
5975.0  'zi    ' = sample cell interface  240
6000.0  'zi    ' = sample cell interface  241
   0	'botio ' = bathymetry  I/O unit (0 no I/O)
   0	'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
   0.25	'tmljmq' = equiv. temp. jump across mixed-layer (degC,  0 no I/O)
   0.25 'tmlnav' = NAVO rp33 T  jump across mixed-layer (degC,  0 no I/O)
   0	'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
   0	'wvlio ' = w-velocity  I/O unit (0 no I/O)
  42	'uvlio ' = u-velocity  I/O unit (0 no I/O)
  43	'vvlio ' = v-velocity  I/O unit (0 no I/O)
   0	'splio ' = speed       I/O unit (0 no I/O)
   0	'temio ' = temperature I/O unit (0 no I/O)
   0	'salio ' = salinity    I/O unit (0 no I/O)
  44	'tthio ' = density     I/O unit (0 no I/O)
   0	'keio  ' = kinetic egy I/O unit (0 no I/O)
E-o-D
cat blkdat.input
unset echo
#aprun -n 5 -N 1  ${W}/ALLmpi/archive/src/archv2data3z
aprun -n 60 -N 12 -S 6  ${W}/ALLmpi/archive/src/archv2data3z

#
endif #if exist

end #h
 @ D = $D + 1
#
end #while

