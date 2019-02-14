#! /bin/csh

# PURPOSE:  Create new Spire and python tools distribution from
#           contents of  !!spire-src!!  directory

# SOURCE:   /usr16/software/spider/utils/create-spire-dist.csh

# Set variable for spire root directory location
set spiroot  = /usr16/software/spider 
   
# Set some variables for input locations
set srcdir   = $spiroot/spire-src

# Set some variables for output locations 
set destdir  = $spiroot/spire-dist

# Set rsync  = compressed, update, preserve executability, 
#              preserve time, follow Symlinks, quiet

# Set rsync options
set sendit   = 'rsync -zuEtLq --out-format="%n%L"  '

# Weed out local directories not in release
set excludes  = "--exclude="RCS" --exclude="Attic" --exclude="utils" --exclude="tst" "  

# Make necessary output dir
mkdir -p $destdir 

echo ' ' ; echo ' Destination: '$destdir

##\rm *.pyc */*.pyc */*/*.pyc */*/*/*.pyc */*/*/*/*.pyc */*/*/*/*/*.pyc  */*/*/*/*/*/*.pyc */*/*/*/*/*/*/*.pyc */*/*/*/*/*/*/*/*.pyc

# ------------------ Copy source directory files ------------------------

echo  ' Copying:    ' $srcdir '  files. '

$sendit -r $excludes --exclude="*.pyc" --exclude="*jnk*" $srcdir/*  $destdir

\rm $destdir/make-changes* $destdir/see-utils*


# ------------------ Create tar archive --------------------------------

chdir  $destdir

echo  ' Creating:   ' $destdir 'tar files.'

\rm  -f spire.tar spire.tar.gz
tar cf spire.tar ./*
gzip spire.tar  

###gzip -k spire.tar  #Keep uncompressed tar file

# ------------------ Remove lib & libsys  directories (also in tar archive)  ------------------------

echo ' Removing:    lib & libsys  directories (left in tar archive)'

\rm -rf $destdir/lib-sys       $destdir/lib

echo ' Replacing:    site-packages'

mkdir -p $destdir/lib/python2.5/site-packages   $destdir/lib/python2.7/site-packages

$sendit -r $excludes --exclude="*.pyc"  $srcdir/lib/python2.5/site-packages/Spider  $destdir/lib/python2.5/site-packages
$sendit -r $excludes --exclude="*.pyc"  $srcdir/lib/python2.5/site-packages/Spire   $destdir/lib/python2.5/site-packages
$sendit -r $excludes --exclude="*.pyc"  $srcdir/lib/python2.7/site-packages/Spider  $destdir/lib/python2.7/site-packages

# ------------------ Finished  ------------------------
echo ' Use:         /usr16/software/spider/utils/tosend.csh  to update distribution release'

exit




# ------------------------------ File notes  ---------------------

# Master copy
#/usr16/software/spider/spire-src    Mb    Has: Attic,RCS,tst,utils     Lacks: *.pyc

# Distribution copy created using: /usr16/software/spider/utils/create-spire-dist.csh
#/usr16/software/spider/spire-dist   Mb    Has: tar= Mb      Lacks: Attic,*.pyc

# Local execution copy
#/usr16/software/spider/spire        Mb    Has tar.gz= Mb 

