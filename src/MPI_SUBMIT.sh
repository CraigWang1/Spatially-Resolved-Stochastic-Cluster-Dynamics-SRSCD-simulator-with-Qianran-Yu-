###### MPI_SUBMIT.sh START ######################
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
# Edit the line below as needed:
#$ -l h_rt=24:00:00,h_data=2G
# Add multiple cores/nodes as needed:
#$ -pe dc* 4   # notice you may also try w/ "-pe shared 4"
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
module load intel   # change if you use any other version of intelmpi or openmpi
module li
echo " "

# substitute the <NAME OF YOUR EXECUTABLE> to run below:
echo '/usr/bin/time -v mpirun -n $NSLOTS ./scdexe >> output.$JOB_ID'
/usr/bin/time -v mpirun -n $NSLOTS ./scdexe >> output.$JOB_ID

# echo job info on joblog:
echo " "
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
###### MPI_SUBMIT.sh STOP ######################

