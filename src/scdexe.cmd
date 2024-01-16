#!/bin/csh -f
#  scdexe.cmd
#
#  UGE job for scdexe built Mon Jan 15 17:22:18 PST 2024
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog
#$ -o /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Threaded:     4-way threaded
#  Resources requested
#$ -pe shared 4
#$ -l h_data=1024M,h_rt=100:00:00,highp
# #$ -q jmarian_*.q
#  Name of application for log
#$ -v QQAPP=openmp
#  Email address to notify
#$ -M cwang@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
# Initialization for openmp threaded execution
#
  unalias *
  set qqversion = 
  set qqapp     = "openmp threads"
  set qqmtasks  = 4
  set qqidir    = /u/scratch/c/cwang/parallel/10_vacancies/src
  set qqjob     = scdexe
  set qqodir    = /u/scratch/c/cwang/parallel/10_vacancies/src
  cd     /u/scratch/c/cwang/parallel/10_vacancies/src
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for scdexe built Mon Jan 15 17:22:18 PST 2024"
  echo ""
  echo "  scdexe directory:"
  echo "    "/u/scratch/c/cwang/parallel/10_vacancies/src
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  'scratch' directory (on each node):"
  echo "    $qqscratch"
  echo "  scdexe 4-way threaded job configuration:"
#
  echo ""
  echo "scdexe started on:   "` hostname -s `
  echo "scdexe started at:   "` date `
  echo ""
#
# Run the user program
#

  source /u/local/Modules/default/init/modules.csh
  module load intel

  echo scdexe "" \>\& scdexe.output.$JOB_ID
  echo ""
  setenv OMP_NUM_THREADS 4

  /usr/bin/time /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe  >& /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.output.$JOB_ID

  echo ""
  echo "scdexe finished at:  "` date `

#
# Cleanup after openmp threaded execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/openmp.log.threaded
 if (`wc -l /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
        head -50 /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmp.log.threaded
        echo " "  >> /u/local/apps/queue.logs/openmp.log.threaded
        tail -10 /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmp.log.threaded
  else
        cat /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID >> /u/local/apps/queue.logs/openmp.log.threaded
  endif
#  cat            /u/scratch/c/cwang/parallel/10_vacancies/src/scdexe.joblog.$JOB_ID           >> /u/local/apps/queue.logs/openmp.log.threaded
  exit (0)
