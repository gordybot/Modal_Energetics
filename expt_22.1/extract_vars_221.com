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
LISN=$RNNM'_archv.2016.lis'
echo $LISN

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
echo $PTHHO

for j in $(seq $jblks $stp $jblke ); do
  echo $j

  i=$iblks
  je=$((j+stp-1))

  echo "along x: do ${i} to ${iblke}"
  echo "along y: do ${j} to ${je}"
  egrep "\.a" $LISN | tail -1000 | wc -l > test
  echo "$j    'jblks'  " > ${PTHNM}extract_${RNNM}_${j}_${i}.in
  echo "$je    'jblke'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}.in
  echo "$iblks    'iblks'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}.in
  echo "$iblke    'iblke'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}extract_${RNNM}_${j}_${i}.in
  awk '{print $0"     \x27numfls\x27"}' test >> ${PTHNM}extract_${RNNM}_${j}_${i}.in
  egrep "\.a" $LISN | tail -1000 >> ${PTHNM}extract_${RNNM}_${j}_${i}.in

# copy screen output to pbs file
  echo "aprun -n 1 ./extract_vars_c004_v12.x  < ./extract_${RNNM}_${j}_${i}.in > ./extract_${RNNM}_${j}_${i}.out"

# make all the pbs files
  echo "#! /bin/sh"                       > ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -N Xj${j}_${i}"             >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -o Xj${j}_${i}.out"         >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -e Xj${j}_${i}.err"         >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -M gordon.stephenson@usm.edu" >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -m bea"                      >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -A NRLSS03755018"            >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -l walltime=20:00:00"        >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -q standard"                    >> ${PTHHO}pbsX_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbsX_j${j}_${i}
  echo "# Change to the specified directory" >> ${PTHHO}pbsX_j${j}_${i}
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsX_j${j}_${i}
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbsX_j${j}_${i}
  echo "aprun -n 1 ./extract_vars_c004_v12.x  < ./extract_${RNNM}_${j}_${i}.in > ./extract_${RNNM}_${j}_${i}.out" >> ${PTHHO}pbsX_j${j}_${i}

  qsub ${PTHHO}pbsX_j${j}_${i}

done

