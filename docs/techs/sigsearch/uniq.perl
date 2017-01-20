#!/usr/bin/perl
#!/usr/local/bin/perl
#
# perl script to create search doc file for unique locations and
# save location with highest CC
#
#to run:
#  uniq.perl < ryrsort.dat > ryrsortu

$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator

# change nsam, nrow and nslice to reflect the size of the target volume
$nsam   = 384;
$nrow   = 384;
$nslice = 110;

# peaks surrounding the largest peak inside a sphere of radius = $r are 
# supressed. mostly, it is the rdius of the motif.
$r      = 6;             # Search radius

# of unique peaks sought (~ # of expected occurence of the motif inside the 
# target volume)
$nwant  = 1000 ;         # Number of points wanted

# ----------------------------- code ---------------------------------

$r3   = 8 * $r * $r * $r;
$len  = $nsam * $nrow * $nslice;
$#buf = $len;

#foreach $item (@buf)
#   {$item = 0;}

$n = 0;
line : while (<>) 
   {
   chop;       # strip record separator
 
   if (!/;.*/) 
       { 
       #print $_ ;

       @Fld = split(' ', $_, 9999);  #   split fields on blanks

       $x  = $Fld[3];
       $y  = $Fld[4];
       $z  = $Fld[5];
       $cc = $Fld[6];
       #print  $len ;

       $iloc = ($z-1) * $nsam * $nrow + ($y-1) * $nsam + ($x - 1);
       if (($iloc < $len) && 
           ( !defined $buf[$iloc]))
           ####($buf[$iloc] == 0))
          { # No point at this location yet
          #print $iloc, $buf[$iloc];

          # See if there is any existing point within radius: r
          for ($iz=-$r; $iz<=$r; $iz++)
             {
             for ($iy=-$r; $iy<=$r; $iy++)
                {
                for ($ix=-$r; $ix<=$r; $ix++)
                   {
                   $rnow = $ix * $ix + $iy * $iy + $iz * $iz;
                   #print 'rnow' . $rnow;
                   if ( $rnow > $r3 ) 
                      {goto fini;}

                   # Within radius now 
                   $ixt =  $x + $ix;  $iyt =  $y + $iy;  $izt =  $z + $iz;

                   $ilocr = ($izt - 1) * $nsam * $nrow + 
                            ($iyt - 1) * $nsam +
                            ($ixt - 1);

              #print  $buf[$ilocr];

                   if ( ($ilocr >= 0) &&
                        ($ilocr < $len) &&
                        ($buf[$ilocr] > 1) )
                      { # Do not copy, another point with > CC within radius  
                      goto fini;
                      }
                   }
                }
             }  # end of: for
          # no point already within radius
          print $_ ;

          if (++$n >= $nwant)
             {exit;}

          fini : $buf[$iloc] = 5.0;   # set this point 
          }   #end of: if 
       } 
   }
