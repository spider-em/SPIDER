#!/bin/csh
#
# SOURCE: /usr16/software/spider/docs/techs/recon1a/Utils/spr2tar.csh
#
# PURPOSE: Puts batch files listed in: ../spider/docs/techs/recon1a/mr1.html
#          into a zipped tar file:  spiproject.year.month.date.tar.gz
#
# NOTE:    mkdirs.py reads: mr1.html, and loads anything with a link to 
#          into the appropriate reconstruction project directory.
#
# CHANGES:   BB Oct-29-04
#            BR Jul-26-05
#            AL Mar-10-09  Added echo of file sources
#            AL Mar-10-09  Migrate to usr8
#            AL Aug-05-14  Altered for no defocus groups
#            AL Apr-25-16  Removed refine.html 
#            AL Dec-28-16  Save uncompressed for use in git
#            AL May-30-17  usr16

setenv ROOTDIR /usr16/software/spider/docs/techs/recon1a
setenv DOCSDIR $ROOTDIR/Docs
setenv PROJDIR myproject

echo ; echo " Read mr1.html  to list reconstruction procedure files -----"

$ROOTDIR/Utils/mkdirs.py -n $PROJDIR

# Create project Docs directory to hold important html files & some useful info.
mkdir -p $PROJDIR/Docs

echo ; echo " Kludge to copy directory specific 'DOT' .* files ----------------"
# Copy .**  files to project directories
cp -vu $ROOTDIR/Procs/Averages/.??*       $PROJDIR/Averages 
cp -vu $ROOTDIR/Procs/Micrographs/.??*    $PROJDIR/Micrographs 
cp -vu $ROOTDIR/Procs/Particles/.??*      $PROJDIR/Particles 
cp -vu $ROOTDIR/Procs/Power_Spectra/.??*  $PROJDIR/Power_Spectra 

echo ; echo " Copy other useful files to project "Docs" directory -------"

# Copy useful files to project "Docs" directory
cp -vu $DOCSDIR/mr1.html      $PROJDIR/Docs
cp -vu $DOCSDIR/mrstyle2.css  $PROJDIR/Docs

echo ; echo " Run SPIDER session and write version to info file -------- "
# Create file "info". "info" contains infomation about the SPIDER executable  
#    and the zipped file that is created by this script 
echo " Zipped Batch File Info.  : " >  $PROJDIR/Docs/info
echo " " >> $PROJDIR/Docs/info
echo "spiproject.`date +%y%m%d`.tar.gz" >> $PROJDIR/Docs/info
echo " " >> $PROJDIR/Docs/info
echo "SPIDER Version Info.  : " >> $PROJDIR/Docs/info
echo " " >> $PROJDIR/Docs/info

# Run a SPIDER session to write the output to info file
spider zyx/cba en >>  $PROJDIR/Docs/info 

# Remove unnecessary files created by SPIDER 
\rm -f LOG.zyx  results.zyx.*

echo ; echo " Tar files to desired location ---------------------"
tar  cvf $ROOTDIR/"spiproject.`date +%y%m%d`.tar"  myproject

echo ; echo " Zip file to desired location ---------------------"
gzip -cf $ROOTDIR/"spiproject.`date +%y%m%d`.tar" > $ROOTDIR/"spiproject.`date +%y%m%d`.tar.gz"

# Delete the temp. dir "myproject"
\rm -rf $PROJDIR

# Link the zipped file for access from the WEB page
echo ; echo " Link the zipped files for easy access ---------------------"
\rm     $ROOTDIR/spiproject.tar 
\rm     $ROOTDIR/spiproject.tar.gz
ln -s  "$ROOTDIR/spiproject.`date +%y%m%d`.tar.gz" $ROOTDIR/spiproject.tar.gz
ln -s  "$ROOTDIR/spiproject.`date +%y%m%d`.tar"    $ROOTDIR/spiproject.tar

echo ; echo " List final tar archive locations -------------------------"
ls -l  $ROOTDIR/spiproject.tar
ls -l  $ROOTDIR/spiproject.tar.gz
ls -l "$ROOTDIR/spiproject.`date +%y%m%d`.tar"
ls -l "$ROOTDIR/spiproject.`date +%y%m%d`.tar.gz"

echo " "
