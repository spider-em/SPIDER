#!/usr/bin/csh
#
# SOURCE:   /usr8/spider/utils/wwwupdate.csh 
#
# PURPOSE: Copies all website SPIDER documents to: CSS's desired location
#          CSS will then copy from that location into external 
#          Wadsworth website once a day
#
# USAGE:   wwwupdate.csh 

# CHANGES: new                                ArDean Leith  May 94
#          python tools                       ArDean Leith  Apr 14
#          tmpdir                             ArDean Leith  Jul 15

# Root location for our local source files
set spider_root  = /usr8/spider    
set tmpdir       = $spider_root/utils/jnk_raw_www_docs        
set spider_send  = /usr8/send/

set spiredir     = /usr8/spider/spire/    
set docsdir      = $spider_root/docs/
set srcdir       = $spider_root/src/ 
set procsdir     = $spider_root/proc/ 
set pubsubdir    = $spider_root/pubsub/ 
set tipsdir      = $spider_root/docs/tips/ 
set toolsdir     = $spider_root/tools-src/ 


# Root location for Wadsworth wwwinternal output files 
# OLD wwwhostdir = '/net/info/usr3/WWW/wwwinternal/spider_doc/'
# OLD wwwhostdir = 'nnewton:/usr3/WWW/wwwinternal/spider_doc/'    
set wwwhost       = spider-stage    
set wwwdir        = /export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider    
set wwwhostdir    = $wwwhost':'$wwwdir    

set wwwdocsdir    = $wwwhostdir/docs

# For Wadsworth WWW. Set rsync = verbose, compressed, update, 
#        preserve executability, preserve time, follow Symlinks
#set sendit = 'rsync -vzuEtL' pre 4/27/2012 for Linux target
set sendit = 'rsync -vzuptL --exclude="RCS" --exclude="Attic" '

pushd .
cd $spider_root/utils

set LOGFILE = $spider_root/utils/the-wwwupdate.LOG

# ----------- Create LOG for the update. -----------------------------

\rm -f $LOGFILE
echo " Update performed:" >  $LOGFILE
date                      >> $LOGFILE
echo "  "                 >> $LOGFILE

# Convert SPIDER manual chapters from text format --------------------
echo " Preparing raw manual chapters.  " 
$spider_root/utils/old2raw.perl 
echo " Converted raw manual chapters to html.  " 
echo " Converted raw manual chapters to html.  "  >> $LOGFILE

# Add headers, etc. to HTML, place in: $tmpdir  ----------------------
$spider_root/utils/raw2docs.perl -wadsworth
echo " Added html headers, trailers, etc. "
echo " Added html headers, trailers, etc. " >> $LOGFILE
echo " Created operations list.   "
echo " Created operations list. " >> $LOGFILE
 
# Copy headerized SPIDER man files  ----------------------------------
echo " Copying headerized SPIDER man files xxxxxxxxxxxxxxxxxxxxxxx" 
ssh $wwwhost mkdir -p $wwwdir/docs/man

$sendit     $tmpdir/man/*               $wwwdocsdir/man/

echo " Copied headerized SPIDER man files " >> $LOGFILE
echo " $sendit     $tmpdir/man/*               $wwwdocsdir/man/"

# Copy associated manual images --------------------------------------
$sendit  $spider_root/man/*.jpg  $wwwdocsdir/man
echo " Copied SPIDER man image files " >> $LOGFILE

# Copy headerized SPIDER doc files  ----------------------------------
echo " Copying headerized SPIDER doc files xxxxxxxxxxxxxxxxxxxxxxx" 
$sendit    $tmpdir/*                    $wwwdocsdir

echo " Copied headerized SPIDER doc files " >> $LOGFILE

# Copy exa, icons, buttons, img & spidui SPIDER doc files ------------
echo " Copying exa, icons, buttons, img & spidui SPIDER files xxxxxxxxxxxxxxxxxxx" 
$sendit    $docsdir/img/*               $wwwdocsdir/img/
$sendit    $docsdir/icons/*             $wwwdocsdir/icons/
$sendit    $docsdir/buttons/*           $wwwdocsdir/buttons/
$sendit    $docsdir/spidui/*            $wwwdocsdir/spidui/
$sendit    $docsdir/spidui/pics/*       $wwwdocsdir/spidui/pics/
$sendit    $docsdir/exa/*               $wwwdocsdir/exa/
$sendit    $docsdir/exa/images/*        $wwwdocsdir/exa/images/

echo " Copied exa, icons, buttons, img & spidui SPIDER files " >> $LOGFILE

# Copy non-headerized SPIDER doc files ------------------------------
echo " Copying non-headerized SPIDER doc files xxxxxxxxxxxxxxxxxxx" 
$sendit    $docsdir/spider.html        \
           $docsdir/spider_license.html $docsdir/spider_avail.html  \
           $docsdir/spi-register.html   $docsdir/spi-download.html  \
           $docsdir/agl.html            $docsdir/*.pdf              \
           $docsdir/*.gif               $docsdir/*.jpg              \
           $docsdir/*.spi               $docsdir/*.css              \
           $docsdir/spider78.html                                   $wwwdocsdir/

echo " Copied non-headerized SPIDER doc files " >> $LOGFILE

# Copy SPIDER source files ------------------------------------------
echo " Copying SPIDER src files xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" 
$sendit -r --exclude="*.o"   --exclude="*.mod"   --exclude="jnk*"   \
           --exclude="*.a"   --exclude="ifort/"  --exclude="gfort/" \
           --exclude="RCS"   --exclude="Attic" \
           $srcdir/*  $wwwhostdir/src
echo " Copied SPIDER source files " >> $LOGFILE

# Copy SPIDER python tools files ------------------------------------------
echo " Copying  $toolsdir   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" 
ssh  $wwwhost mkdir -p $wwwhostdir/tools
$sendit    $toolsdir/*.html   $toolsdir/*.py $wwwhostdir/tools
$sendit    $toolsdir/bin/*py  $wwwhostdir/tools/bin
$sendit -r $toolsdir/docs     $wwwhostdir/tools
echo " Copied python tools files" >> $LOGFILE


# Copy SPIDER pubsub files ------------------------------------------
echo  " Copying pubsub & procs files xxxxxxxxxxxxxxxxxxxxxxxxxxx" 
$sendit  $pubsubdir/*    $wwwhostdir/pubsub  
$sendit  $procsdir/*     $wwwhostdir/proc  
echo  " Copied  proc & pubsub  files "  >> $LOGFILE

# Copy spire files --------------------------------------------------
echo  " Copying spire files (except *.tar) xxxxxxxxxxxxxxxxxxxxx" 
$sendit -r --exclude="*.tar" --exclude="Attic"  --exclude="RCS" \
           --exclude="lib" \
           $spiredir/doc  $spiredir/readme.html $spiredir/tosend $wwwhostdir/spire  
echo  " Copied  spire  files "  >> $LOGFILE


# Copy external tip files -------------------------------------------
echo  " Copying external docs/tips files xxxxxxxxxxxxxxxxxxxxxx" 
ssh $wwwhost mkdir -p $wwwdir/docs/tips
$sendit    $tipsdir/index_send.html                     $wwwdocsdir/tips/index.html  
$sendit    $tipsdir/utilities.html  $tipsdir/spiprogramming.html   $wwwdocsdir/tips/ 
$sendit    $tipsdir/timing.html     $tipsdir/timebprp.spi          $wwwdocsdir/tips/ 
echo " Copied  external tips files "  >> $LOGFILE


# Copy external techs files -----------------------------------------
echo " Copying external docs/techs/* xxxxxxxxxxxxxxxxxxxxxxxxx" 

$sendit -r --exclude="bzvol.dat"  --exclude="dev/" \
           --exclude="backup/"    --exclude="Attic"  --exclude="RCS" \
           $docsdir/techs    $wwwdocsdir
echo " Copied external SPIDER tech files " >> $LOGFILE

echo " Copying various external tar files xxxxxxxxxxxxxxxxxxxxxxxxx" 

$sendit -v $docsdir/techs/supclass/tar_archive/*  $wwwdocsdir/techs/supclass/tar_archive 
$sendit -v $docsdir/techs/verify/tar_archive/*    $wwwdocsdir/techs/verify/tar_archive 
$sendit -v $docsdir/techs/recon/batcharch/*       $wwwdocsdir/techs/recon/batcharch 

#echo " Copying SPIDER distribution file xxxxxxxxxxxxxxxxxxxxxxxxx" 
#$sendit -v $spider_send/spiderweb.*tar.gz        $wwwhostdir/download

echo " FINISHED"    >> $LOGFILE

echo "  "
echo "  "
echo "  "
echo " You may need to check: spider/download/spider! "
echo " Use: ssh spider-stage"
exit 1



#############  Unused below here ##

# cp /usr8/send/spiderweb.tar.gz        $wwwhostdir/download/spiderweb.15.10.tar.gz
# scp -p spi-download.html $wwwhost':/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider/docs
# " scp -p reference_based.tar.gz nnewton:/usr3/WWW/wwwinternal/spider_doc/spider/download"
#  Legacy: sterecon stuff is in /~leith/www/sterecon
# scp -p spiderweb.21.00.tar.gz*   http-rev:/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider/download/spiderweb.21.00.tar.gz
