#! /bin/sh Makefile
# On GNU/Linux you must: alias make "/bin/sh Makefile $1" for this to work!!
#
# PURPOSE:      Chooses specific Makefile      
#
# NORMAL USAGE: make [version-selector]
#
# SAMPLE USAGE: make 
#               make [mp][f][d] e.g. 
#               on Linux
#                  make mp (to create spider_linux_mp)
#                  make all (to create all current executables)
#                  make new (to update all current executables version)

TARGET="$1"
echo TARGET: $TARGET

DEST="DEST=_tmp"

DEBUG=

if expr "$TARGET" : '.*d' > /dev/null
then # Want debugger suitable output, set flags & edit TARGET
  DEBUG="QFFLAGS=-g QLFLAGS=-g SUFFIX=d"
  # echo DEBUG: $DEBUG
  TARGET=`echo $TARGET | tr -d 'd'`
fi

UNAME=`uname`
#echo UNAME: $UNAME
UNAMEm=`uname -m`
#echo UNAMEm: $UNAMEm

case "$UNAME" in
#  Find current operating system
   Linux)  
      PLATFORM=_linux 
      #VERSIONS=" _linux _linux_mp_opt64  _linux_mp_intel _linux_mpi_opt64 _linux_mp_intel64 " 
      VERSIONS=" _linux _linux_mp_opt64  _linux_mp_intel  _linux_mp_intel64 " 
     
      if expr "$TARGET" : '.*mp' > /dev/null
      then #  Want mp with  version
         WHICH="mp"
      else
         if expr "$TARGET" : '.*f' > /dev/null
         then #  Want  version
            WHICH=""
         fi
      fi;;
esac


case "$TARGET" in
   new) # Update all VERSIONS to Next highest version number

        for SELECT in $VERSIONS
        do
           FILE=../bin/spider$SELECT"_tmp"
           if [ -e $FILE ] 
           then
              echo "Creating: `./Nextversion ../bin/spider$SELECT`"
              cp $FILE ../bin/`./Nextversion ../bin/spider$SELECT`
           else
              echo No: $FILE
           fi
        done
        exit ;;

   all) # Make all spiders to _tmp

        for SELECT in $VERSIONS
        do
           echo Running: make -f Makefile$SELECT $DEST
           make -f Makefile$SELECT $DEST
           echo " "
        done
        exit ;;
esac

#echo WHICH: $WHICH
 
# echo WHICH makefile selection to screen
echo Running: make -f Makefile$WHICH $DEBUG $DEST

# activate proper Makefile
make -f Makefile$WHICH $DEBUG $DEST

exit






