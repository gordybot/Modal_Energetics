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
#INTID='2zi'
LISN=$RNNM'_archv.2016_'$INTID'_TS.lis'
echo $LISN

#range i:1-35, j:1-60
jblks=1
jblke=35
iblks=1
iblke=60
stp=5

# path where the simulation is done
PTHNM=$WORKDIR'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
PTHHO=$HOME'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
echo $PTHNM
echo $PTHHO

for j in $(seq $jblks $stp $jblke ); do
  echo $j

  i=$iblks
  je=$((j+stp-1))

  echo "along x: do ${i} to ${iblke}"
  echo "along y: do ${j} to ${je}"
  egrep "\.a" $LISN | tail -10000 | wc -l > test
  echo "$j    'jblks'  " > ${PTHNM}mean_N_${RNNM}_${j}_${i}.in
  echo "$je    'jblke'  " >> ${PTHNM}mean_N_${RNNM}_${j}_${i}.in
  echo "$iblks    'iblks'  " >> ${PTHNM}mean_N_${RNNM}_${j}_${i}.in
  echo "$iblke    'iblke'  " >> ${PTHNM}mean_N_${RNNM}_${j}_${i}.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}mean_N_${RNNM}_${j}_${i}.in
  echo "'$INTID'   'interp ID'  " >> ${PTHNM}mean_N_${RNNM}_${j}_${i}.in

# copy screen output to pbs file
  echo "aprun -n 1 ./mean_N_c004_v1.x  < ./mean_N_${RNNM}_${j}_${i}.in > ./mean_N_${RNNM}_${j}_${i}.out"

# make all the pbs files
  echo "#! /bin/sh"                       > ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -N Nj${j}_${i}"             >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -o Nj${j}_${i}.out"         >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -e Nj${j}_${i}.err"         >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -M gordon.stephenson@usm.edu" >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -m bea"                      >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -A ONRDC38645449"            >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -l walltime=30:30:00"        >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -q standard"                 >> ${PTHHO}pbsN_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbsN_j${j}_${i}
  echo "# Change to the specified directory" >> ${PTHHO}pbsN_j${j}_${i}
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsN_j${j}_${i}
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbsN_j${j}_${i}
  echo "aprun -n 1 ./mean_N_c004_v1.x  < ./mean_N_${RNNM}_${j}_${i}.in > ./mean_N_${RNNM}_${j}_${i}.out" >> ${PTHHO}pbsN_j${j}_${i}

  qsub ${PTHHO}pbsN_j${j}_${i}
  
done

