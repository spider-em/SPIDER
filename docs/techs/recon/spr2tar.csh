#!/bin/csh
#
# Puts batch files on ../spider/docs/techs/recon/mr.html
# into a zipped tar file, spiproject.year.month.date.tar.gz
#
# mkdirs.py reads mr.html, and loads anything with a link to the newprogs
# directory.
#
# BB Oct-29-04
# BR Jul-26-05
# AL Mar-10-09  Added echo of file sources
# AL Mar-10-09  Migrate to usr8
# AL Dec-28-16  Save uncompressed for use in git

setenv SRCDIR /usr16/software/spider/docs/techs/recon

echo " "; echo " Read mr.html  to get list of desired reconstruction procedure files -----"

./mkdirs.py -n myproject

# Create Docs directory to hold important html files & some useful info.
mkdir -p myproject/Docs

echo " " ; echo " Copy other useful files to "Docs" directory -----"

# Copy useful files to "Docs" directory
cp -vu $SRCDIR/mr.html      myproject/Docs/.
cp -vu $SRCDIR/mrstyle.css  myproject/Docs/.
cp -vu $SRCDIR/mrspire.html myproject/Docs/.
cp -vu $SRCDIR/refine.html  myproject/Docs/.

echo " " ; echo " Run a SPIDER session and write version to info file ----- "
# Create file "info". "info" contains infomation about the spider executable  
#    and the zipped file that is created by this script 
echo "Zipped Batch File Info.  : " >  myproject/Docs/info
echo " " >> myproject/Docs/info
echo "spiproject.`date +%y%m%d`.tar.gz" >> myproject/Docs/info
echo " " >> myproject/Docs/info
echo "SPIDER Version Info.  : " >> myproject/Docs/info
echo " " >> myproject/Docs/info

# Run a SPIDER session to write the output to info file
spider zyx/cba en >>  myproject/Docs/info 

# Remove unnecessary files created by SPIDER 
\rm -f LOG.zyx  results.zyx.*

echo " " ; echo " Tar and  zip file to desired location -----"
# Tar, zip and move the file to desired location
tar cvf "spiproject.`date +%y%m%d`.tar" myproject

gzip -cf "spiproject.`date +%y%m%d`.tar" > "spiproject.`date +%y%m%d`.tar.gz"

\cp "spiproject.`date +%y%m%d`.tar.gz" $SRCDIR/Attic/batcharch/

# Delete the temp. dir "myproject"
#\rm -rf myproject

# Link the zipped file for access from the WEB page
unlink $SRCDIR/spiproject.tar.gz
ln -s "$SRCDIR/spiproject.`date +%y%m%d`.tar.gz" $SRCDIR/spiproject.tar.gz

echo " " ; echo " List finalized tar archive location -----"
ls -l "$SRCDIR/Attic/batcharch/spiproject.`date +%y%m%d`.tar.gz"
ls -l "$SRCDIR/spiproject.tar.gz"
ls -l "spiproject.`date +%y%m%d`.tar"

echo " "
