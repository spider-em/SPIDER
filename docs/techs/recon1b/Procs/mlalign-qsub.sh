#!/bin/tcsh
#PBS -q normal
#PBS -j oe
#PBS -m e

# USAGE: qsub -v DATEXT='dat',results=1,scratchlink='tmpalign',workDir=`pwd` -l nodes=1:ppn=4:xeon,mem=2500mb,scratch=1152mb mlalign-qsub.sh

source /packages/run/modules-2.0/init/csh 
source /storage/brno2/home/tapu/local/spider.cshrc
module add mpich3
module list

echo "Working directory $workDir"
echo "Data extension $DATEXT, #CPUs $PBS_NUM_PPN"

cd ${workDir}
rm -f $scratchlink
ln -sv $SCRATCHDIR $scratchlink
spider spi/$DATEXT @framealign $results ncpus=$PBS_NUM_PPN 
rm -f $scratchlink

