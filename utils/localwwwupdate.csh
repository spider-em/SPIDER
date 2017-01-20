#!/bin/csh
#
# SOURCE: spider/utils/localwwwupdate.csh
#
# PURPOSE: Updates local SPIDER HTML documents from ASCII manual pages
#          Add headers & trailers to HTML manuals and documents
#
# CHANGES: new                             ArDean Leith  May 94
#          removed src.txt file usage      ArDean Leith  Aug 05
#          for usr8                        ArDean Leith  Aug 10
#          exa images                      ArDean Leith  Jan 13


set spider_root = /usr8/spider

set LOGFILE     = "$spider_root/utils/localwwwupdate.LOG"

# Create a LOG for update. ------------------------------------------
\rm -f $LOGFILE
echo " Update performed :" >  $LOGFILE
date                       >> $LOGFILE
echo "  "                  >> $LOGFILE
echo " Creating HTML manual chapters from text format"

# Convert manual chapters from ASCII text to HTML -------------------

$spider_root/utils/old2raw.perl

echo " Created HTML manual chapters from text format"
echo " Created HTML manual chapters from text format" >> $LOGFILE
echo " ------------------------";

# Add HTML headers, etc. to docs, man, exa HTML pages ---------------

$spider_root/utils/raw2docs.perl

echo " Added html headers to docs & man pages." 
echo " Added html headers, etc." >> $LOGFILE
echo " ------------------------";

# Copy css files from rawdocs & man directories ---------------------
#      -z   Compress file data during the transfer
#      -u   Skip files that are newer
#      -E   Preserve executability
#      -t   Preserve modification times
#       ./man/niceman.css   ./user.css    ./ex.css    ./buttons.css     /exa/ex.css

rsync -zuEt $spider_root/man/*.css          $spider_root/docs/man
rsync -zuEt $spider_root/rawdocs/*.css      $spider_root/docs
rsync -zuEt $spider_root/rawdocs/exa/*.css  $spider_root/docs/exa

echo " ------------------------";

echo " Finished  " ; echo " " 

exit 1
