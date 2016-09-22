#! /bin/csh

# PURPOSE:  Create new python tools distribution (without SPIRE)

# ORIGINAL SOURCE: /usr8/spider/spire/spire_linux-1.5.5

# CURRENT SOURCE:  /usr8/spider/tools-src/tools-create-dist.csh

# Set some variables for input locations
set spiroot = /usr8/spider    
set srcdir  = $spiroot/tools-src

# Set some variables for output locations
set destdir  = $spiroot/tools-dist

# Set rsync  = verbose, compressed, update, 
#              preserve executability, preserve time, follow Symlinks

set sendit   = 'rsync -zuEtL --out-format="%n%L"  '

set excludes = "--exclude="RCS" --exclude="Attic" --exclude="\*.pyc" --exclude="howto-create-dist" "  

# Make necessary dir
mkdir -p $destdir 

# ------------------ Copy source files --------------------------------

echo  'Copying tools src files. '

$sendit -r $excludes --exclude="*.pyc" --exclude="orig" --exclude="tst" $srcdir/*  $destdir

# ------------------ Remove tools shell script files --------------------------------

#cd $destdir/bin/

# montage-spi is to overcome montage name collision

\rm $destdir/bin/montage-spi
 
\rm $destdir/bin/binarytree
\rm $destdir/bin/classavg
\rm $destdir/bin/ctfcircle
\rm $destdir/bin/ctfdemo
\rm $destdir/bin/ctfgroup
\rm $destdir/bin/ctfmatch
\rm $destdir/bin/emancoords2spiderdoc
\rm $destdir/bin/emanrctcoords2spiderdoc
\rm $destdir/bin/mkapps
\rm $destdir/bin/mkfilenums
\rm $destdir/bin/montagefromdoc
\rm $destdir/bin/montage
\rm $destdir/bin/pyplot
\rm $destdir/bin/qview
\rm $destdir/bin/scatter
\rm $destdir/bin/spiconvert
\rm $destdir/bin/verifybyview
\rm $destdir/bin/viewstack
\rm $destdir/bin/xmippsel2spiderdoc
\rm $destdir/bin/xplor

\rm $destdir/bin-python/python 


# ------------------ Create tar archive --------------------------------
cd $destdir
\rm -f tools.tar tools.tar.gz
tar cvf tools.tar *
tar --delete -vf  tools.tar use-spider-utils-create-tools-dist-to-update
gzip tools.tar  

echo ' 'Use:  /usr8/spider/utils/tosend.csh  to update distribution release

exit


# ------------------------------ File notes  ---------------------

# Master copy
#/usr8/spider/tools-src    79Mb    Has: Attic,RCS     Lacks: *.pyc

# Distribution copy created using: /usr8/spider/utils/create-tools-dist.csh
#/usr8/spider/tools-dist  105Mb    Has: tar=51Mb      Lacks: tst, *.pyc, RCS 

# Local execution copy
#/usr8/spider/tools        91Mb    Has tar.gz=21Mb

