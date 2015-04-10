#! /bin/csh

# PURPOSE:  Create new python tools distribution (without SPIRE)

# ORIGINAL SOURCE: /usr8/spider/spire/spire_linux-1.5.5

# CURRENT SOURCE:  /usr8/spider/tools-src/-create-dist.csh

# Set some variables for input locations
set spiroot = /usr8/spider    
set srcdir  = $spiroot/tools-src

# Set some variables for output locations
set destdir  = $spiroot/tools-dist

# Set rsync  = verbose, compressed, update, 
#             preserve executability, preserve time, follow Symlinks

set sendit   = 'rsync -zuEtL --out-format="%n%L"  '

set excludes = "--exclude="RCS" --exclude="Attic" --exclude="\*.pyc" --exclude="howto-create-dist" "  

# Make necessary dir
mkdir -p $destdir 

# ------------------ Copy source files --------------------------------

echo  'Copying tools src files. '

$sendit -r $excludes --exclude="*.pyc" --exclude="orig" --exclude="tst" $srcdir/*  $destdir

# ------------------ Remove tools shell script files --------------------------------

#cd $destdir/bin/

\rm $destdir/bin/binarytree           $destdir/bin/ctfmatch   $destdir/bin/montage-spi 
\rm $destdir/bin/verifybyview         $destdir/bin/classavg   $destdir/bin/scatter
\rm $destdir/bin/emancoords2spiderdoc $destdir/bin/pyplot     $destdir/bin/viewstack   
\rm $destdir/bin/xmippsel2spiderdoc   $destdir/bin/ctfdemo    $destdir/bin/montage        
\rm $destdir/bin/montagefromdoc       $destdir/bin/xplor      $destdir/bin/ctfgroup    
\rm $destdir/bin/spiconvert           $destdir/bin/qview    
\rm $destdir/bin/ctfcircle            $destdir/bin/mkfilenums          

\rm $destdir/bin-python/python 


# ------------------ Create tar archive --------------------------------
cd $destdir
tar cvf tools.tar *
gzip tools.tar  

echo ' 'Use:  /usr8/spider/utils/tosend.csh  to update distribution release

exit

