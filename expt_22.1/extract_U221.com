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
#LISN=$RNNM'_archv.2012.lis'
LISN='E221_ufiles'
echo $LISN

#range i:1-60, j:1-35
jblks=1
jblke=35
iblks=1
iblke=60
stp=3

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
echo "$j    'jblks'  " > ${PTHNM}X_u_${RNNM}_${j}_${i}.in
echo "$je    'jblke'  " >> ${PTHNM}X_u_${RNNM}_${j}_${i}.in
echo "$iblks    'iblks'  " >> ${PTHNM}X_u_${RNNM}_${j}_${i}.in
echo "$iblke    'iblke'  " >> ${PTHNM}X_u_${RNNM}_${j}_${i}.in
echo "'$RNNM'   'run num'  " >> ${PTHNM}X_u_${RNNM}_${j}_${i}.in
echo "'2zi'     'intID'   " >> ${PTHNM}X_u_${RNNM}_${j}_${i}.in
echo "383       'num ufiles'" >>${PTHNM}X_u_${RNNM}_${j}_${i}.in
cat $LISN >>${PTHNM}X_u_${RNNM}_${j}_${i}.in

# copy screen output to pbs file
echo "aprun -n 1 ./extract_u_221.x  < ./X_u_${RNNM}_${j}_${i}.in > ./output_${RNNM}_${j}_${i}.out"

# make all the pbs files
  echo "#! /bin/sh"                       > ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -N Xuj${j}_${i}"             >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -o Xuj${j}_${i}.out"         >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -e Xuj${j}_${i}.err"         >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -M gordon.stephenson@usm.edu" >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -m bea"                      >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -A ONRDC38645449"            >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -S /bin/bash"                >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -l walltime=06:00:00"        >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -q standard"                    >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "# Change to the specified directory" >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "cd ${WORKDIR}/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbs_XU_j${j}_${i}
  echo "aprun -n 1 ./extract_u_221.x  < ./X_u_${RNNM}_${j}_${i}.in > ./X_u_${RNNM}_${j}_${i}.out" >> ${PTHHO}pbs_XU_j${j}_${i}

  qsub ${PTHHO}pbs_XU_j${j}_${i}

done

