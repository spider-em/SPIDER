#!/bin/tcsh
#PBS -q normal
#PBS -j oe

# USAGE: qsub -l nodes=1:ppn=4 -v DATEXT=$DATEXT,workDir=`pwd`,LOGFILE='log.ctffind4' master-ctffind4-qsub.sh

source /packages/run/modules-2.0/init/csh 
source /storage/brno2/home/tapu/local/spider.cshrc
module add xmipp-3.0.1
module add mpich3
module list

set DEFAULTLOG='log.ctffind4'
if ! $?LOGFILE then
    echo "LOGFILE not defined, using '$DEFAULTLOG'... "
    set LOGFILE=$DEFAULTLOG
endif

echo "Working directory $workDir"
echo "Data extension $DATEXT, #CPUs $PBS_NUM_PPN"

cd ${workDir}
set RESULTS=`Nextresults results.spi | awk -F . '{print $NF}'`
spider spi/$DATEXT @ctffind4 $RESULTS ncpus=$PBS_NUM_PPN
