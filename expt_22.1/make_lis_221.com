#! /bin/sh
#make_lis.com, nrl, mcb, 2016/06/25 
#create lis file with *.a file names
#max range is 1-32, 01-00

# filename in
FIN='archv.2016'
RNM='221_'

is=245 #days
ie=261 #long
js=1   #hours
je=23  #long

#directory where files are stored
SD=$WORKDIR'/hycom/GLBc0.04/expt_22.1/input'
SD2=$WORKDIR'/hycom/GLBc0.04/expt_22.1/input'

#last day  minus 1
ee=$(( ie - 1 ))
#echo $ee

#first day  plus 1
gg=$(( is + 1 ))
#echo $gg

# difference in days
ff=$((ie-is))
#echo $ff

# first do less than one day ---------------------------
#echo $FIN.lis

# save regional.depth.a regional.grid.a
echo "$SD2/regional.depth.a" > "$RNM$FIN.lis"
echo "$SD2/regional.grid.a" >> "$RNM$FIN.lis"

# temp fix to test it
echo "$SD2/tidal.rh.a" >> "$RNM$FIN.lis"
echo "$SD2/cb.a" >> "$RNM$FIN.lis"

hh=0
if [ $ff -eq 0 ]; then

for i in $(seq $is $is); do
    for j in $(seq   $js $je); do
        hh=$((hh+1))
        echo "$SD/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done
fi

# if more days -------------------------------------------

if [ $ff -ge 1 ]; then

# first day
for i in $(seq $is $is); do
    for j in $(seq   $js 23); do
        hh=$((hh+1))
        echo "$SD/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done


# then do whole days first
# -w does the padding :-)

for i in $(seq $gg  $ee); do
    for j in $(seq   0 23); do
        hh=$((hh+1))
        echo "$SD/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done

# last day
for i in $(seq $ie $ie); do
    for j in $(seq   0 $je); do
        hh=$((hh+1))
        echo "$SD/${RNM}${FIN}_$(printf %03d $i)_$(printf %02d $j).a" >> "$RNM$FIN.lis"
   done
done
fi

echo "number of time steps: ${hh}"
