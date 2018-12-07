#! /bin/sh
#reads *.a files from *.list file and extracts the data in the *.a files
#and appends it to a binary output file

#|: control operator to separate a sequence commands in a pipeline
#egrep: match pattern in file
#tail -7:  list last 7 lines
#wc -l: count lines
#>> extract_all.in: opens file extract_all.in and appends text
#> extract_all.in: opens file extract_all.in and (over)writes text
#awk prints text read in file test
#>& redirects both the stdout and the stderr of command to filename
#time: time a simple command or give resource usage

# run number
EXPT='06'
num='1'
RNNM=$EXPT$num

#range i:1-30, j:1-18
jblks=16
jblke=18
iblks=1
iblke=30

# path where the simulation is done
PTHNM=$WORKDIR'/hycom/GLBc0.08/expt_'$EXPT'.'$num'/'
echo $PTHNM


i=$iblks
j=$jblks

echo "along x: do ${i} to ${iblke}"
echo "along y: do ${j} to ${jblke}"
echo "$jblks    'jblks'  " > ${PTHNM}pmode_${RNNM}_${j}_${i}.in
echo "$jblke    'jblke'  " >> ${PTHNM}pmode_${RNNM}_${j}_${i}.in
echo "$iblks    'iblks'  " >> ${PTHNM}pmode_${RNNM}_${j}_${i}.in
echo "$iblke    'iblke'  " >> ${PTHNM}pmode_${RNNM}_${j}_${i}.in
echo "'$RNNM'   'run num'  " >> ${PTHNM}pmode_${RNNM}_${j}_${i}.in

# copy screen output to pbs file
echo "aprun -n 1 ./pmodes_pbot_61.x  < ./pmode_${RNNM}_${j}_${i}.in > ./output_pmode_${RNNM}_${j}_${i}.out"
