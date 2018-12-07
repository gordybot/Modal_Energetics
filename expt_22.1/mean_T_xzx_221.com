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
EXPT='22'
num='1'
RNNM=$EXPT$num
INTID='1z'
LISN=$RNNM'_archv.2016_'$INTID'_TS.lis'
echo $LISN

#number of time steps
NT=360
#NT=1

#range i:1-60, j:1-35
jblks=1
jblke=35
iblks=1
iblke=60
stp=5

# path where the simulation is done
PTHNM=$WORKDIR'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
PTHHO=$HOME'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
echo $PTHNM

for j in $(seq $jblks $stp $jblke ); do
  echo $j

  i=$iblks
  je=$((j+stp-1))

  echo "along x: do ${i} to ${iblke}"
  echo "along y: do ${j} to ${je}"
  egrep "\.a" $LISN | tail -10000 | wc -l > test
  echo "$j    'jblks'  " > ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  echo "$je    'jblke'  " >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  echo "$iblks    'iblks'  " >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  echo "$iblke    'iblke'  " >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  echo "'$INTID'   'interp ID'  " >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  echo "$NT       'num times'  " >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in
  egrep "\.a" $LISN | tail -10000 >> ${PTHNM}mean_T_${RNNM}_${j}_${i}.in

# make all the pbs files
  echo "#! /bin/sh"                       > ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -N Tj${j}_${i}"              >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -o Tj${j}_${i}.out"          >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -e Tj${j}_${i}.err"          >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -M gordon.stephenson@usm.edu" >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -m bea"                      >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -A ONRDC38645449"            >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -l walltime=30:00:00"        >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -q standard"                 >> ${PTHHO}pbsT_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbsT_j${j}_${i}
  echo "# Change to the specified directory" >> ${PTHHO}pbsT_j${j}_${i}
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsT_j${j}_${i}
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbsT_j${j}_${i}
  echo "aprun -n 1 ./mean_T_c004_v1.x  < ./mean_T_${RNNM}_${j}_${i}.in > ./mean_T_${RNNM}_${j}_${i}.out" >> ${PTHHO}pbsT_j${j}_${i}

  qsub ${PTHHO}pbsT_j${j}_${i}

done




