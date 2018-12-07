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

# path where the simulation is done
PTHNM=$WORKDIR'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
PTHHO=$HOME'/hycom/GLBc0.04/expt_'$EXPT'.'$num'/'
echo $PTHNM
echo $PTHHO

# start a selected process
echo "'$RNNM'   'run num'  "   > ${PTHNM}comb_eigen_${RNNM}.in
echo "'$INTID'  'interp ID'  " >> ${PTHNM}comb_eigen_${RNNM}.in

echo $RNNM

echo "aprun -n 1 ./combine_eigen_c008_v1.x   < ./comb_eigen_${RNNM}.in > ./comb_eigen_${RNNM}.out"

# make all the pbs files
echo "#! /bin/sh"                       > ${PTHHO}pbscomb_eig_all
echo "#PBS -N CE_all"             >> ${PTHHO}pbscomb_eig_all
echo "#PBS -o CE_all.out"         >> ${PTHHO}pbscomb_eig_all
echo "#PBS -e CE_all.err"         >> ${PTHHO}pbscomb_eig_all
echo "#PBS -M gordon.stephenson@usm.edu" >> ${PTHHO}pbscomb_eig_all
echo "#PBS -m bea"                      >> ${PTHHO}pbscomb_eig_all
echo "#PBS -A ONRDC38645449"            >> ${PTHHO}pbscomb_eig_all
echo "#PBS -S /bin/bash"                >> ${PTHHO}pbscomb_eig_all
echo "#PBS -l walltime=00:30:00"        >> ${PTHHO}pbscomb_eig_all
echo "#PBS -q debug"                    >> ${PTHHO}pbscomb_eig_all
echo "#PBS -l select=1:ncpus=1"         >> ${PTHHO}pbscomb_eig_all
echo "# Change to the specified directory" >> ${PTHHO}pbscomb_eig_all
echo "cd $WORKDIR/hycom/GLBc0.04/expt_${EXPT}.${num}" >> ${PTHHO}pbscomb_eig_all
echo "# Execute the serial executable on 1 core"  >> ${PTHHO}pbscomb_eig_all
echo "aprun -n 1 ./combine_eigen_c004_v1.x  < ./comb_eigen_${RNNM}.in > ./comb_eigen_${RNNM}.out" >> ${PTHHO}pbscomb_eig_all

# why not submit the job too :-)
qsub ${PTHHO}pbscomb_eig_all

