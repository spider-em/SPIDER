#!/usr/local/bin/perl
#
# perl script to convert x** or x* to [v**] 
#
#to run:
#  perl  subx.perl < file.spi > junk ; 

$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator

while (<>) 
    {
    chop;                               # strip record separator
    if ( m/[xX](\d+)/ )
      {
      #print 'found:' . $_  ;
      s/[xX](\d\d)/[v$1]/g;
      #print 'subdd:' . $_  ;
      }
    print $_;
    }


    #s/(x[\d|\dd]
    # s/(\.)(\d)(\d)(\d)+/$1$2$3$4/g;    
    # s/ambientIntensity (.*)/ambientColor $1 $1 $1 /g;     
    # s/(.*)\](.*)/$1}$2/g;

    #print 'found $_' ;
