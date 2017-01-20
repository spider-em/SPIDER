#!/bin/sh
#
# SOURCE: /usr8/spider/utils/send-fftw.sh
#
# PURPOSE: Copy fftw3 files to distribution without duplicating a lot of stuff

# du -ak | grep '\.a' | grep -v '\.am'

# Source of fftw3 files
fftw_src="/usr8/spider/fftw"

# Distribution destination
fftw_dest="/usr8/send/spider/fftw"

# FFTW3 targets included in copy
#fftw_sources="fftw3-intel fftw3-32 fftw3-opt64 fftw3-intel64 fftw3-osx-32 fftw3-osx-64  fftw3-osx-32-pgi fftw3-osx-64-pgi " apr15
fftw_sources="fftw3-intel fftw3-32 fftw3-opt64 fftw3-intel64" 

# Name of FFTW3 distribution directory in distribution destination
build_dir="${fftw_dest}/FFTW3_dist"
   
# Copy the target FFTW3 object library files
echo Copying: config, howto, AUTHORS COPY, .libs lib files   
for i in ${fftw_sources}
  do
  cd ${fftw_src}/${i}
  mkdir -p $fftw_dest/${i} 
  cp -rpuf --parents config* how* AUTHORS COPY* .libs lib $fftw_dest/${i}   
done
echo " "

# Copy my README file
echo Copying: ${fftw_src}/README 
cp -pu ${fftw_src}/README ${fftw_dest}

echo Dir. of:  $fftw_dest
ls  $fftw_dest

exit # ------------------------------------

# Copy and clean the FFTW3 distribution directory (only needs to be done once)
cd ${fftw_src}
echo Copying: ${fftw_src}/fftw3-opt64 to FFTW3_distribution dir  
cp -rpu ${fftw_src}/fftw3-opt64 ${build_dir}  
cd ${build_dir}
make clean

echo " "

exit 


# may want to remove .libs directories??
