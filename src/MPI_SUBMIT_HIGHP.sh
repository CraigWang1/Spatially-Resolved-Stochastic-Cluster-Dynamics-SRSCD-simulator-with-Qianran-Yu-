###### MPI_SUBMIT.sh START ######################
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
# Edit the line below as needed:
#$ -l h_rt=300:00:00,h_data=5G,arch=intel-gold*,highp
#$ -pe dc* 1 # notice you may also try w/ "-pe shared 4"
# Add multiple cores/nodes as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo "Job $JOB_ID will run on:   "
cat $PE_HOSTFILE
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
module unload intel
module load mpich #change if you want intelmpi or openmpi, I found mpich works best 
module li
echo " "

# substitute the <NAME OF YOUR EXECUTABLE> to run below:
echo '/usr/bin/time -v mpirun -n $NSLOTS ./scdexe >> output.$JOB_ID'
/usr/bin/time -v `which mpirun` -n $NSLOTS ./scdexe >> output.$JOB_ID

# echo job info on joblog:
echo " "
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
###### MPI_SUBMIT.sh STOP ######################

