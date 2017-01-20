#! /bin/csh
#
# SOURCE:  /usr8/spider/utils/tosend.csh
#
# PURPOSE: Update local copy of SPIDER distribution        
#          The distribution copy is currently in:  /usr8/send/spider/...
#
# NEXT:    The combined SPIDER/WB distribution can then be placed in a compressed
#            tar archive and put in the download directory on the
#            external SPIDER website.
#
# AUTHOR:  ArDean Leith                May 1994
#          FFTW added                  Nov 2005
#          Removed onecpp.perl call    Mar 2008
#          Nextresults                 Feb 2009
#          Rewrite for Linux & rsync   Mar 2009
#          Rewrite for /usr8           Jul 2010

echo 
echo "Did you run: /usr8/spider/docs/techs/recon/spr2tar.csh          first?"
echo "Did you run: /usr8/spider/docs/techs/recon1a/Utils/spr2tar.csh  first?"
echo "Did you run: /usr8/spider/docs/techs/recon1b/Utils/spr2tar.csh  first?"
echo "Did you run: /usr8/spider/utils/create-tools-dist.csh           first?"
echo 

# Set some variables for input locations
set spiroot    = /usr8/spider    

set srcdir     = $spiroot/src
set bindir     = $spiroot/bin
set mandir     = $spiroot/man
set procdir    = $spiroot/proc
set docdir     = $spiroot/docs
set tipsdir    = $spiroot/docs/tips
set pubsubdir  = $spiroot/pubsub
set fftwdir    = $spiroot/fftw
set toolsdir   = $spiroot/tools
set spiredir   = $spiroot/spire 
set spiredistdir = $spiroot/spire-dist 
   
set jwebdir    = /usr8/web/jweb 

# Set some variables for output locations
set destroot   = /usr8/send     

set srcdest    = $destroot/spider/src
set docdest    = $destroot/spider/docs 
set mandest    = $destroot/spider/man
set procdest   = $destroot/spider/proc
set bindest    = $destroot/spider/bin
set fftwdest   = $destroot/spider/fftw
set pubsubdest = $destroot/spider/pubsub
set spiredest  = $destroot/spider/spire
set toolsdest  = $destroot/spider/tools
set jwebdest   = $destroot/web/jweb

# Make necessary dir
mkdir -p $destroot/spider    $srcdest $docdest $mandest 
mkdir -p $procdest $fftwdest $pubsubdest
mkdir -p $destroot $jwebdest

# Set rsync = verbose, compressed, update, 
#             preserve executability, preserve time, follow Symlinks
set sendit  = 'rsync -zuEt --out-format="%n%L" --delete '

set excludes = "--exclude="RCS" --exclude="Attic" --exclude="dev" "  


# ------------------ Copy source files --------------------------------

echo  'Copying src files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

$sendit -d  $srcdir/*.f   \
            $srcdir/*.INC \
            $srcdir/Make* \
            $srcdir/*.inc \
            $srcdir/Nextversion   $srcdest

# Replace Makebody.inc with distribution version
$sendit -d $srcdir/Makebody.inc.send  $srcdest/Makebody.inc

$sendit -d $excludes $srcdir/Makefile_samples/*  $srcdest/Makefile_samples
 
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
 
# --------------------- Copy PubSub* files ---------------------------

echo 'Copying pubsub files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

$sendit -d $excludes  $pubsubdir/*   $pubsubdest

# --------------------- Copy proc files ---------------------------

echo 'Copying proc files. xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
$sendit -d $excludes       \
           $procdir/*.dat  \
           $procdir/*.bat  \
           $procdir/*.spi  \
           $procdir/*.py   \
           $procdir/*.img  \
           $procdir/*.tom  \
           $procdir/*.perl      $procdest

# --------------------- Copy Python tools files ----------------------

echo 'Copying Python tools files  xxxxxxxxxxxxxxxxxxxxxxxxx'
$sendit -r $excludes \
           $toolsdir/readme*      \
           $toolsdir/install.html \
           $toolsdir/tools.tar.gz \
           $toolsdir/docs             $toolsdest

# --------------------- Copy Spire files ---------------------------

echo 'Copying Spire distribution files. xxxxxxxxxxxxxxxxxxx'

$sendit -r $spiredir/readme*                      $spiredest
$sendit -r $excludes $spiredir/doc                $spiredest
$sendit -r $spiredistdir/spire_linux-1.5.5.tar.gz $spiredest

# --------------------- Copy JWeb files ---------------------------

echo 'Copying JWeb Linux, & Windows files. xxxxxxxxxxxxxxxx'

$sendit -r $excludes       \
           $jwebdir/linux  \
           $jwebdir/win    \
           $jwebdir/src    $jwebdest

# --------------------- Copy  html doc & tech files --------------------------

echo 'Copying html doc & tech files. xxxxxxxxxxxxxxxxxxxxxx'

$sendit -r  $excludes --exclude="tips"                    \
           --exclude="techs/lgstr/tomo/data"              \
           --exclude="techs/lgstr/tomo/output"            \
           --exclude="exa/images/bp3fpart*dat"            \
           --exclude="techs/recon1a/natproc_data_mics.tar.gz" \
           $docdir/*   $docdest

# ----------------------------------------------------------------

echo ' '
echo 'SPIDER successfully copied to: distribution dir '
echo ' '
echo 'Check for extra  tar archives '
echo 'Update FFTW files with: send-fftw.sh '
echo 'Update executables with make in src dir '
echo 'Update Web  with:       /usr8/web/utils/tosend.sh '
echo 'touch /usr8/send/spider/bin/CONTAINS_SPIDER_RELEASE_24.00 '
echo 'Archive and compress the distribution in: /usr8/send '
echo 'set wwwdir = spider-stage:/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider '
echo 'scp -p /usr8/send/spiderweb.23.02.tar.gz  $wwwdir/download '
echo 'Edit: /usr8/spider/docs/spi-download.html '
echo 'scp -p /usr8/spider/docs/spi-download.html  $wwwdir/docs '
echo 'Update external web pages using: /usr8/spider/utils/wwwupdate.sh '
echo ' '

#echo 'Update executables from OSX spider '

exit 0
