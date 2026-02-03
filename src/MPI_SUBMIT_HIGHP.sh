###### MPI_SUBMIT.sh START ######################
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
# Edit the line below as needed:
#$ -l h_rt=336:00:00,h_data=5G,arch=intel-gold*,highp
#$ -pe dc* 1 # notice you may also try w/ "-pe shared 4"
# Add multiple cores/nodes as needed:
# Email address to notify
#$ -M $USER@mail
# Notify when
#$ -m bea

# get the folder name the script is run from
FOLDER_NAME=$(basename "$PWD")

# create a timestamp: eg. 2025-01-13_14-52-30
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")

# combined destination folder
DEST="${HOME}/simulation_results/${FOLDER_NAME}_${TIMESTAMP}"

# make the folder
mkdir -p "$DEST"

# Copy setup configuration for logging
cp restart.txt sink.txt constants.h "$DEST"/

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

# log output files
cp species*.txt sink*.txt Desorbed.txt "$DEST"/

###### MPI_SUBMIT.sh STOP ######################

