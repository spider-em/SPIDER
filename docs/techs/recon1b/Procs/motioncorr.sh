#!/bin/bash
#PBS -N frameAlignOFgpu
#PBS -q gpu
#PBS -l gpu=1
#PBS -l walltime=23:59:00
#PBS -l scratch=96gb:ssd
#PBS -l mem=6gb
#PBS -j oe
#PBS -m e

# USAGE: 
# qsub -l nodes=1:ppn=2 -v workingDir=`pwd`,outputDir='Dosefgpu',refFrame=1,buf=50,DATEXT='dat' dosefgpu-qsub.sh
# (Micrographs are assumed to be in current directory and ending in '.mrcs')

# Adapted from frame_alignment_xmipp.bsh (Jirka Novacek)

#------------------------------------------------------------------------------
# Description
# - Falcon II frames alignment
# - makes "global" alignment according to Li et al. 2013 (Nat. Methods)
# - input files are stored as MRC stack in current directory, ending in '.mrcs'
# - Outputs single MRC file in the outputDir
# - scratch (ideally SSD) used for faster file access
# - cuda 5.0 needed 
#------------------------------------------------------------------------------

# Temporary filenames

# unprocessed micrographs will be stored locally here temporarily
tempUnalignedDir="Unaligned"

# completed-frames directory, relative to workingDir - where the already processed frames are moved temporarily
doneFramesDir="CompletedFrames"

# set utils path dosefgpu_driftcorr
export PATH=$PATH:/storage/brno2/home/tapu/bin/

# cuda 5.0 needed for global alignment (dosefgpu_driftcorr)
module add cuda-5.0

# load XMIPP to convert to SPIDER format
module add xmipp-3.0.1

#set first or last frame as reference
if [ "$refFrame" == 0 ]; then
    refFrame="nst"
else
    refFrame="ned"
fi

#Clean SCRATCHDIR after exit
trap 'clean_scratch' TERM EXIT

echo "Job started on `hostname` at `date`"

# create output directories in working folder
cd $workingDir

if [ ! -d $doneFramesDir ]; then
    mkdir -v $doneFramesDir/
fi

if [ ! -d $outputDir ]; then
    mkdir -v $outputDir/
fi

# create directories on scratch disk
if [ ! -d $SCRATCHDIR/$tempUnalignedDir ]; then
    mkdir -v $SCRATCHDIR/$tempUnalignedDir/
fi

if [ ! -d $SCRATCHDIR/$outputDir ]; then
    mkdir -v $SCRATCHDIR/$outputDir/
fi

# move files to fast scratch - process max $buf per batch
while [ $( ls | grep  '\.mrcs$' | head -n ${buf} | wc -l) -ge 1 ]; do
    
    nrOfBufferedImagesToProcess=$( ls | grep  '\.mrcs$' | head -n $buf | wc -l)
    
    echo "Copying $nrOfBufferedImagesToProcess micrographs to scratch buffer"
    date
    
    # Copy micrographs to scratch directory (in background)
    cp $( ls | grep  '\.mrcs$' | head -n $buf) $SCRATCHDIR/$tempUnalignedDir/ &

    cd $SCRATCHDIR/$tempUnalignedDir/
        while [ "$nrOfBufferedImagesToProcess" -ge 1 ]; do
            # wait until certain amount of files has been copied in the background to the scratch disk 
            sleep 300
            ls -ltr *.mrcs   # I'm curious how long it takes to copy
            
            for i in $(ls *.mrcs); do

                micStem=$(echo $i | awk '{printf "%s", substr($1,1,length($1)-5)}')

                # globally align (create drift corrected stack) frames according to Li et al. 2013 (Nat. Methods)
                dosefgpu_driftcorr $i -$refFrame 0 -ssc 0 -slg 0 -dsp 0 -fcs tmp$micStem.mrc
                
                # convert to SPIDER format
                xmipp_image_convert -i tmp$micStem.mrc -o $workingDir/$outputDir/$micStem.$DATEXT

                # move the source frames to $doneFramesDir (temporarily)
                mv $workingDir/${i} $workingDir/$doneFramesDir/
                    
                # remove scratch version of micrograph
                rm -v $SCRATCHDIR/$tempUnalignedDir/${i}
                    
                # decrement #remaining micrographs
                ((nrOfBufferedImagesToProcess--))
            done
        done
            
    cd $workingDir

done

rm tmp$micStem.mrc

# move micrographs back to original location
mv $workingDir/$doneFramesDir/* $workingDir/

echo "Job finished on `hostname` at `date`"
