#! /bin/sh

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
  echo "$j    'jblks'  " > ${PTHNM}neig_${RNNM}_${j}_${i}.in
  echo "$je    'jblke'  " >> ${PTHNM}neig_${RNNM}_${j}_${i}.in
  echo "$iblks    'iblks'  " >> ${PTHNM}neig_${RNNM}_${j}_${i}.in
  echo "$iblke    'iblke'  " >> ${PTHNM}neig_${RNNM}_${j}_${i}.in
  echo "'$RNNM'   'run num'  " >> ${PTHNM}neig_${RNNM}_${j}_${i}.in
  echo "'$INTID'   'interp ID'  " >> ${PTHNM}neig_${RNNM}_${j}_${i}.in

# copy screen output to pbs file
  echo "aprun -n 1 ./NEIG_c004_v1.x  < ./neig_${RNNM}_${j}_${i}.in > ./neig_${RNNM}_${j}_${i}.out"

# make all the pbs files
  echo "#! /bin/sh"                       > ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -N EIGj${j}_${i}"             >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -o EIGj${j}_${i}.out"         >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -e EIGj${j}_${i}.err"         >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -M gordon.stephenson@usm.edu" >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -m bea"                      >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -A ONRDC38645449"            >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -l walltime=09:00:00"        >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -q standard"                 >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "# Change to the specified directory" >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbsEIG_j${j}_${i}
  echo "aprun -n 1 ./NEIG_c004_v1.x  < ./neig_${RNNM}_${j}_${i}.in > ./neig_${RNNM}_${j}_${i}.out" >> ${PTHHO}pbsEIG_j${j}_${i}

  qsub ${PTHHO}pbsEIG_j${j}_${i}

done

