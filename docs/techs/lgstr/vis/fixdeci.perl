#!/usr/local/bin/perl
#
# Author: ArDean Leith  1 June 2000
#
# file: /net/java/marko/decimate/fixdeci.perl
# tests: /net/java/marko/decimate/
#
# Active Source location:/usr/local/spider/bin/fixdeci.perl 
#
# Purpose: perl script to convert a SPIDER "PP LL" output file which has
# been run through "Decimate" back into original scaling units.
#
# Usage:  fixdeci.perl input1 input2 output 
# Test:   fixdeci.perl tube.geo tube_gm1_e.iv jnkout.iv 
#

if (@ARGV[0] eq "-h"  || @ARGV[0] eq "?" || @ARGV[0] eq "")
   {    #help requested
   printf(" Purpose: Correctly scales Inventor file created from Decimate after 'PP LL'  \n");
   printf("  \n");
   printf(" Usage: fixdeci.perl deci-inputfile deci-outputfile  outputfile  \n");
   printf("        First input file  is Inventor file from Explorer after \'PP LL\'.\n");
   printf("        Second input file is Inventor file from Decimate.\n");
   printf("        Output file is scaled Inventor file from Decimate.\n"); 
   exit;
   }

#set default input file name
$ppfilename = "tube.geo";
if (@ARGV[0])
   {$ppfilename = @ARGV[0];}

open(INFILE, "< $ppfilename") 
     || die "Can't read: $ppfilename\n";

#set second default input file name
$decifilename = "tube_gm1_e.iv";
if (@ARGV[1])
   {$decifilename = @ARGV[1];}

open(INFILED, "< $decifilename") 
     || die "Can't read: $decifilename\n";

print 'Opened: ' . $ppfilename . " & " . $decifilename ."\n";

$keepgo = 1;

$minx   = 10e30;
$miny   = $minx;
$maxx   = -$minx;
$maxy   = $maxx;

line : while (($_ = <INFILE>)  &&  $keepgo == 1) 
    {
    chop;       # strip record separator

    if (/Transform/) 
        { $Transform =  $_ ; }

    elsif (/scaleFactor/) 
        { $scaleFactor =  $_ ; }

    elsif (/vertex/) 
        {
#        print 'vertex : ' . $_ . "\n";
        }

    elsif ((/[\d\.e-]+ [\d\.e-]+ [\d\.e-]+,/) ||
           (/[\d\.e-]+ [\d\.e-]+ [\d\.e-]+ \]/))

        {
#        print 'got : ' . $_ . "\n";
        $lastmatch = $_;

#       @fields = split(/[\s,]+/, $_); could not get to work
        @fields = /([\d.-e]+) ([\d.-e]+) ([\d.-e]+)/;

#        print '  0= ' . @fields[0] . '  1= ' . @fields[1] . '  2= ' . @fields[2] . "\n";

        $minx = @fields[0] < $minx ? @fields[0] : $minx;
        $miny = @fields[1] < $miny ? @fields[1] : $miny;
        $minz = @fields[2] < $minz ? @fields[2] : $minz;

        $maxx = @fields[0] > $maxx ? @fields[0] : $maxx;
        $maxy = @fields[1] > $maxy ? @fields[1] : $maxy;
        $maxz = @fields[2] > $maxz ? @fields[2] : $maxz;
        if (/[\d\.]+ [\d\.]+ [\d\.]+ \]/)
           { $keepgo = 0;  }
        }
    next line;
    } 

#   print 'lastmatch: ' . $lastmatch . "\n";
#   print '  0= ' . @fields[0] . '  1= ' . @fields[1] . '  2= ' . @fields[2] . "\n";

$avgx  = ($maxx + $minx) / 2.0;
$avgy  = ($maxy + $miny) / 2.0;
$avgz  = ($maxz + $minz) / 2.0;

$fx    = ($maxx - $minx) / 2.0;
$fy    = ($maxy - $miny) / 2.0;
$fz    = ($maxz - $minz) / 2.0;

$scale = $fx;
$scale = $fy > $scale ? $fy : $scale;
$scale = $fz > $scale ? $fz : $scale;

print '  minx:  ' . $minx  . '  miny: ' . $miny . '  minz: ' . $minz . "\n";
print '  maxx:  ' . $maxx  . '  maxy: ' . $maxy . '  maxz: ' . $maxz . "\n";
print '  fx=    ' . $fx    . '  fy=   ' . $fy   . '  fz=   ' . $fz   . "\n";
print '  scale= ' . $scale . "\n";

$inpts = 0;

# rewind the decimate input file
seek(INFILED,0,0);

#set  default output file name
$outfilename = "jnkout.iv";
if (@ARGV[2])
   {$outfilename = @ARGV[2];}

open(OUTFILE, "> $outfilename") 
     || die "Can't open: $outfilename\n";

print 'Opened: ' . $outfilename ."\n";

@factors[0] = 1.0;
@factors[1] = 1.0;
@factors[2] = 1.0;

line : while (($_ = <INFILED>)) 
    {
    chop;       # strip record separator

    if (/scaleFactor/) 
        { 
        @factors = /([\d.-]+) ([\d.-]+) ([\d.-]+e*[\d\.-]+)/;
# print 'factors= ' . @factors[0] . '  , ' . @factors[1] . '  , ' . @factors[2] . "\n";
        print  OUTFILE $_ . "\n";
        }

    elsif (/V1.0 ascii/) 
        {     # version 2.* needed for vertexOrdering change
        print OUTFILE "#Inventor V2.0 ascii   " .   "\n";
        }
    elsif (/ambientColor/ || /specularColor/) 
        {     # deletes some color info
        print OUTFILE " " .   "\n";
        }
    elsif (/Coordinate3/) 
        {     # must copy transform saved from decimate input file
        print OUTFILE 
           "   ShapeHints { vertexOrdering  COUNTERCLOCKWISE" .   "\n";
        print OUTFILE "       creaseAngle 2.0 }" .   "\n";

        print OUTFILE "      " . $Transform .   "\n";
        print OUTFILE "          " . $scaleFactor . "\n";

        print  OUTFILE $_ . "\n";
        $inpts = 1;
        }

    elsif ( ($inpts > 0) &&
           ((/[\d\.-]+e*[\d\.-]* [\d\.-]+e*[\d\.-]* [\d\.-]+e*[\d\.-]*,/) ||
            (/[\d\.-]+e*[\d\.-]* [\d\.-]+e*[\d\.-]* [\d\.-]+e*[\d\.-]* \]/)))
        {  # this is a coordinate triplet

# printf  "%s \n", $_;

        # save any prefix or suffix on the triplet
        $prefix = $_ ;
        $prefix =~ /(\s*[a-zA-Z]*\s*\[*)/;
# print " prefix: " . $1 . ":". "\n";

        # find the x, y & z values
        @fields = /([\d\.-]+e*[\d\.-]*) ([\d\.-]+e*[\d\.-]*) ([\d\.-]+e*[\d\.-]*)/;
# print '  0= ' . @fields[0] . '  1= ' . @fields[1] . '  2= ' . @fields[2] . "\n";

        # scale the x, y & z values
        $newx = (@fields[0] * $scale + $avgx) / @factors[0];
        $newy = (@fields[1] * $scale + $avgy) / @factors[1];
        $newz = (@fields[2] * $scale + $avgz) / @factors[2];

        # print  the scaled x, y & z values
# printf  "%s %9.3f %9.3f %9.3f , %s\n", $1,$newx,$newy,$newz;
        printf  OUTFILE "%s %9.3f %9.3f %9.3f", $1,$newx,$newy,$newz;

        if (/\]/)
           { # end of coordinate listing
           $inpts = 0;  
           print  OUTFILE " ]\n";
           }
        else
           { print OUTFILE ", \n"; }

        }
      else
        {  # just copy this line to output
        print  OUTFILE $_ . "\n";
        }
   }

     
exit;


