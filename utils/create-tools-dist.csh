#! /bin/csh

# PURPOSE:  Create new python tools distribution (without SPIRE)

# ORIGINAL SOURCE: /usr8/spider/spire/spire_linux-1.5.5

# CURRENT SOURCE:  /usr16/software/spider/utils/create-tools-dist.csh

# Set some variables for input locations
set spiroot = /usr16/software/spider    
set srcdir  = $spiroot/tools

# Set some variables for output locations
set destdir  = $spiroot/tools-dist

# Set rsync  = verbose, compressed, update, 
#              preserve executability, preserve time, follow Symlinks

set sendit   = 'rsync -zuEtL --out-format="%n%L"  '

set excludes = "--exclude="RCS" --exclude="Attic" --exclude="\*.pyc" --exclude="howto-create-dist" "  

# Make necessary dir
mkdir -p $destdir 

##\rm *.pyc */*.pyc */*/*.pyc */*/*/*.pyc */*/*/*/*.pyc */*/*/*/*/*.pyc  */*/*/*/*/*/*.pyc */*/*/*/*/*/*/*.pyc */*/*/*/*/*/*/*/*.pyc

# ------------------ Copy source files --------------------------------

echo  'Copying tools src files. '

$sendit -r $excludes --exclude="*.pyc" --exclude="orig" --exclude="tst" $srcdir/*  $destdir

# ------------------ Create tar archive --------------------------------
cd $destdir
\rm -f tools.tar tools.tar.gz
tar cvf tools.tar *
gzip tools.tar  

echo ' 'Use:  /usr16/software/spider/utils/tosend.csh  to update distribution release

exit


# ------------------------------ File notes  ---------------------

# Master copy
#/usr16/software/spider/tools-src    79Mb    Has: Attic,RCS     Lacks: *.pyc

# Distribution copy created using: /usr16/software/spider/utils/create-tools-dist.csh
#/usr16/software/spider/tools-dist  105Mb    Has: tar=40Mb      Lacks: *.pyc, RCS 

# Local execution copy
#/usr16/software/spider/tools        91Mb    Has tar.gz=21Mb

