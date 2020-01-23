#!/usr/bin/csh
#
# SOURCE:  /usr16/software/spider/utils/wwwupdate.csh 
#
# PURPOSE: Copies all website SPIDER documents to: CSS's desired location
#          CSS will then copy from that location into external 
#          Wadsworth website once a day
#
# USAGE:   wwwupdate.csh 
#
# CHANGES: new                                ArDean Leith  May 94
#          python tools                       ArDean Leith  Apr 14
#          tmpdir                             ArDean Leith  Jul 15
#          usr16                              ArDean Leith  Apr 17
#
# Root location for our local source files
set spider_root  = /usr16/software/spider    
set tmpdir       = /usr16/software/spider/utils/jnk_raw_www_docs        
set spider_send  = /usr16/software/send
set spiredir     = /usr16/software/spider/spire/
    
set docsdir      = $spider_root/docs/
set srcdir       = $spider_root/src/ 
set procsdir     = $spider_root/proc/ 
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
set sendit   = 'rsync -vzuptL  --exclude="RCS" --exclude="Attic" '

alias senditr   'rsync -vzuptLr --delete --exclude="RCS" --exclude="Attic" --exclude="dev" '

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

# Copy headerized SPIDER doc files  ----------------------------------
echo " Copying headerized SPIDER doc files xxxxxxxxxxxxxxxxxxxxxxx" 
senditr    $tmpdir/*                    $wwwdocsdir
echo " Copied headerized SPIDER doc files " >> $LOGFILE

# Copy exa, icons, buttons, img & spidui SPIDER doc files ------------
echo " Copying exa, icons, buttons, img & spidui SPIDER files xxxxxxxxxxxxxxxxxxx" 
$sendit    $docsdir/img/*               $wwwdocsdir/img/
$sendit    $docsdir/buttons/*           $wwwdocsdir/buttons/

$sendit -r    $docsdir/spipylib            $wwwdocsdir 
$sendit -r    $docsdir/exa                 $wwwdocsdir

echo " Copied exa, icons, buttons, img & spidui SPIDER files " >> $LOGFILE

# Copy non-headerized SPIDER doc files ------------------------------
echo " Copying non-headerized SPIDER doc files xxxxxxxxxxxxxxxxxxx" 
$sendit    $docsdir/spider.html        \
           $docsdir/spider_license.html $docsdir/spider_avail.html  \
           $docsdir/spi-register.html   $docsdir/spi-download.html  \
           $docsdir/*.pdf               $docsdir/*.css              \
           $docsdir/spider78.html                                   $wwwdocsdir/

echo " Copied non-headerized SPIDER doc files " >> $LOGFILE

# Copy SPIDER source files ------------------------------------------
echo " Copying SPIDER src files xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" 

senditr    --exclude="*.o"   --exclude="*.mod"      --exclude="jnk*"   \
           --exclude="*.a"   --exclude="ifort-jnk"  --exclude="gfort"  \
           $srcdir           $wwwhostdir/src

echo " Copied SPIDER source files " >> $LOGFILE


# Copy SPIDER proc files ------------------------------------------
echo  " Copying procs files xxxxxxxxxxxxxxxxxxxxxxxxxxx" 
$sendit  $procsdir/*     $wwwhostdir/proc  
echo  " Copied  pubsub  files "  >> $LOGFILE


# Copy spire files --------------------------------------------------
echo  " Copying spire files (except *.tar) xxxxxxxxxxxxxxxxxxxxx" 
$sendit -r --delete --exclude="*.tar.gz" --exclude="*.pyc"  \
                    --exclude="spire/bin/python2.5" --exclude="spire/bin/python2.7" --exclude="lib" \
           $spiredir   $wwwhostdir/spire  

$sendit  -r --delete --exclude="*.pyc"  \
           $spiredir/lib/python2.5/site-packages/Spider  $wwwhostdir/spire/lib/python2.5/site-packages/  
$sendit  -r --delete --exclude="*.pyc"  \
           $spiredir/lib/python2.5/site-packages/Spire   $wwwhostdir/spire/lib/python2.5/site-packages/  
$sendit  -r --delete --exclude="*.pyc"  \
           $spiredir/lib/python2.7/site-packages/Spider  $wwwhostdir/spire/lib/python2.7/site-packages/  

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

echo " Copying techs/supclass & techs/verify tar files xxxxxxxxxxxxxxxxxxxx" 

$sendit  $docsdir/techs/supclass/*.tar        $wwwdocsdir/techs/supclass/ 
$sendit  $docsdir/techs/verify/verify*.tar    $wwwdocsdir/techs/verify/ 

#echo " Copying SPIDER distribution file xxxxxxxxxxxxxxxxxxxxxxxxx" 
#$sendit -v $spider_send/spiderweb.*tar.gz        $wwwhostdir/download

echo " FINISHED"    >> $LOGFILE

echo "  "
echo "  "
echo "  "
echo " You may need to check: spider/download/spider! "
echo " Use: ssh spider-stage  cd /export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider"
exit 1



#############  Unused below here ##

# cp //usr16/software/send/spiderweb.tar.gz        $wwwhostdir/download/spiderweb.15.10.tar.gz
# scp -p spi-download.html $wwwhost':/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider/docs
# " scp -p reference_based.tar.gz nnewton:/usr3/WWW/wwwinternal/spider_doc/spider/download"
#  Legacy: sterecon stuff is in /~leith/www/sterecon
# scp -p spiderweb.21.00.tar.gz*   http-rev:/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider/download/spiderweb.21.00.tar.gz

# Copy SPIDER python tools files ------------------------------------------
echo " Copying  $toolsdir   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" 
ssh  $wwwhost mkdir -p $wwwhostdir/tools
$sendit    $toolsdir/readme   $toolsdir/install.sh $wwwhostdir/tools
$sendit    $toolsdir/*.html   $toolsdir/*.py $wwwhostdir/tools
$sendit    $toolsdir/bin/*py                 $wwwhostdir/tools/bin
senditr    $toolsdir/docs                    $wwwhostdir/tools
echo " Copied python tools files" >> $LOGFILE

# Copy SPIDER pubsub files ------------------------------------------
echo  " Copying pubsub & procs files xxxxxxxxxxxxxxxxxxxxxxxxxxx" 
set pubsubdir    = $spider_root/pubsub/ 
$sendit  $pubsubdir/*    $wwwhostdir/pubsub  
$sendit  $procsdir/*     $wwwhostdir/proc  
echo  " Copied  proc & pubsub  files "  >> $LOGFILE
