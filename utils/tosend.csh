#! /bin/csh
#
# SOURCE:  /usr16/software/spider/utils/tosend.csh
#
# PURPOSE: Update local copy of SPIDER distribution        
#          The distribution copy is currently in:  /usr16/software/send/spider/...
#
# NEXT:    The combined SPIDER/Web distribution can then be placed in a compressed
#            tar archive and put in the download directory on the
#            external SPIDER website.
#
# AUTHOR:  ArDean Leith                May 1994
#          FFTW added                  Nov 2005
#          Removed onecpp.perl call    Mar 2008
#          Nextresults                 Feb 2009
#          Rewrite for Linux & rsync   Mar 2009
#          Rewrite for /usr8           Jul 2010
#          Rewrite for /usr16          Apr 2017
#          spire plus tools            Nov 2018
#
echo 
echo "Did you run: /usr16/software/spider/docs/techs/recon/Utils/spr2tar.csh  first?"
echo "Did you run: /usr16/software/spider/docs/techs/recon1a/Utils/spr2tar.csh  first?"
echo "Did you run: /usr16/software/spider/docs/techs/recon1b/Utils/spr2tar.csh  first?"
echo "Did you run: /usr16/software/spider/utils/create-spire-dist.csh           first?"
echo 

# Set some variables for input locations
set spiroot      = /usr16/software/spider    

set srcdir       = $spiroot/src
set bindir       = $spiroot/bin
set mandir       = $spiroot/man
set procdir      = $spiroot/proc
set docdir       = $spiroot/docs
set tipsdir      = $spiroot/docs/tips
set fftwdir      = $spiroot/fftw
set spiredir     = $spiroot/spire 
set spiredistdir = $spiroot/spire-dist 
   
set jwebdir      = /usr16/software/web/jweb 

# Set some variables for output locations
set destroot     = /usr16/software/send     

set srcdest      = $destroot/spider/src
set docdest      = $destroot/spider/docs 
set mandest      = $destroot/spider/man
set procdest     = $destroot/spider/proc
set bindest      = $destroot/spider/bin
set fftwdest     = $destroot/spider/fftw
set spiredest    = $destroot/spider/spire
set jwebdest     = $destroot/web/jweb

# Make necessary dir
mkdir -p $destroot  
mkdir -p $destroot/spider     $srcdest $docdest $mandest 
mkdir -p $procdest $fftwdest  
mkdir -p $jwebdest

# Set rsync = compressed, update, preserve executability, 
#             preserve time, follow Symlinks
set sendit  = 'rsync -zuEt --out-format="%n%L" --delete '

set excludes = "--exclude="RCS" --exclude="Attic" --exclude="dev" "  

# ------------------ Copy source files --------------------------------

echo  'Copying src files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

$sendit -d  $srcdir/*.f      \
            $srcdir/*.INC    \
            $srcdir/Make*    \
            $srcdir/*.inc    \
            $srcdir/makeall  \
            $srcdir/Nextversion   $srcdest

# Replace Makebody.inc with distribution version
\cp -p $srcdir/Makebody.inc.send  $srcdest/Makebody.inc

$sendit -d  $excludes $srcdir/Makefile_samples/*  $srcdest/Makefile_samples
$sendit -dv $excludes $srcdir/ifort_mods/*        $srcdest/ifort_mods

# --------------------- Copy man files -------------------------------

echo 'Copying man files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

$sendit -d  $mandir/*.man    \
            $mandir/*.also   $mandest

# --------------------- Copy bin files -------------------------------

echo 'Copying bin/Nextresults  file. xxxxxxxxxxxxxxxxxxxxxx'

$sendit -d $bindir/Nextresults   $bindest

# ------------------- Copy external tip files ------------------------

echo "Copying external docs/tips files xxxxxxxxxxxxxxxxxxxx" 

$sendit -d $tipsdir                    $docdest
$sendit -d $tipsdir/index_send.html    $docdest/tips/index.html 
$sendit -d $tipsdir/utilities.html      \
           $tipsdir/spiprogramming.html \
           $tipsdir/timebprp.spi        \
           $tipsdir/timing.html        $docdest/tips            
 
# --------------------- Copy proc files ------------------------------

echo 'Copying proc files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
$sendit -d $excludes       \
           $procdir/*.dat  \
           $procdir/*.bat  \
           $procdir/*.spi  \
           $procdir/*.py   \
           $procdir/*.img  \
           $procdir/*.tom  \
           $procdir/*.perl      $procdest

# --------------------- Copy Spire and Python tools files ------------

echo 'Copying Spire and Python tools files. xxxxxxxxxxxxxx'

$sendit -r $spiredistdir/*    $spiredest

# --------------------- Copy JWeb files -----------------------------

echo 'Copying JWeb files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

$sendit -r $excludes                    \
           $jwebdir/src                 \
           $jwebdir/linux.tar.gz    $jwebdest

# --------------------- Copy  html doc & tech files -----------------

echo 'Copying html doc & tech files. xxxxxxxxxxxxxxxxxxxxxx'

$sendit -rL  $excludes --exclude="tips"                       \
           --exclude="tar_archive/rct2*"                      \
           --exclude="spiproject.1*"                          \
           --exclude="spiproject.tar.gz"                      \
           --exclude="techs/lgstr/data"                  \
           --exclude="techs/lgstr/output"                \
           --exclude="exa/images/bp3fpart*dat"                \
           --exclude="techs/recon1a/natproc_data_mics.tar.gz" \
           $docdir/*   $docdest

# \rm -rf /usr16/software/send/spider/docs/techs/lgstr/data/*
# \rm -rf /usr16/software/send/spider/docs/techs/lgstr/output/*

# ------------------------------------------------------------------

echo ' '
echo 'SPIDER successfully copied to: distribution dir '
echo ' '
echo 'Check for extra  tar archives with du -a | grep '\.tar' '
echo 'Update FFTW files with: send-fftw.sh '
echo 'Update executables with makeall in src dir (on gyan)'
echo 'Update Web  with:       /usr16/software/web/utils/tosend.sh '
echo 'touch /usr16/software/send/spider/bin/CONTAINS_SPIDER_RELEASE_26.03 '
echo 'Archive and compress the distribution in: /usr16/software/send '
echo 'set wwwdir = spider-stage:/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider '
echo 'scp -p /usr16/software/send/spiderweb.25.02.tar.gz  $wwwdir/download '
echo 'Edit:  /usr16/software/spider/docs/spi-download.html '
echo 'scp -p /usr16/software/spider/docs/spi-download.html  $wwwdir/docs '
echo 'Update external web pages using: /usr16/software/spider/utils/wwwupdate.csh '
echo ' '


exit 0
rm -v spider/docs/techs/recon*/spiproject.1*
