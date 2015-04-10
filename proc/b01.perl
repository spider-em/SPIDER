#!/usr/bin/perl
# b01.perl    (Should be executable) ArDean Leith July 2001
# Usage: b01.perl | spider dat

$\ = "\n";                           # Set Perl output record separator

$first = 1;                          # Flag for first pipe return

#spi("MD");           
#spi("TERM OFF");                    # Divert terminal output to RESULTS file

$pipename = "TMP_SPIDER_PIPE.pipe";  # Pipe name
print STDERR "Opening pipe: $pipename \n"; 
$FIFO = &openregpipe($pipename);     # Open pipe from SPIDER
###print STDERR "Opened pipe: $pipename \n"; 

spi("[testvar]=12");                 # Set contents of SPIDER variable: [testvar]
$regval = getreg("[testvar]");       # Retrieve SPIDER variable: [testvar]
print STDERR "  [testvar] = $regval \n";

spi("[testvar1]=98765.43");          # Set contents of SPIDER variable: [testvar1]
$regval = getreg("[testvar1]");      # Retrieve SPIDER variable: [testvar1]
print STDERR "  [testvar1] = $regval \n";
	
spi("IQ FI [exists]");               # SPIDER pipe test
spi("b01.perl");                     # Test existance of this file
	
$regval = getreg("[exists]");        # [exists] = 1 if b01.perl exists
print STDERR " [exists] = $regval \n";
		
for ($ifile=1; $ifile < 3; $ifile++) # File loop example
   {
   $iret=&b02($ifile);               # Invoke a Perl SPIDER procedure (b02)
   print STDERR "  perl subroutine b02 returned: $iret \n";
   }
spi("EN");                           # End SPIDER session
exit;
	
	

sub b02 # A Perl SPIDER procedure example ----------------------------------
   {                                # INPUT:  File number
   local($ival,$str);
   $str = sprintf("%3.3d",$_[0]);   # SPIDER likes 001 type names
   spi("IQ FI [yes-exists]");       # Recovers information in [yesexists]
   spi("jnk$str");                  # Appends "jnk" to file number to make filename
   $ival = getreg("[yes-exists]");  # Ask SPIDER for value of variable: [yes-exists]
   return $ival;                    # Returns: value of variable: [yes-exists]
   }
 	
	
# My Common Perl support routines for SPIDER Usage
	
sub openregpipe             # Opens FIFO for SPIDER register input  ------------
  {                         # INPUT:  pipe name    (Argument #1)
  use IO::Handle;
  local($pipename,$iret);
  STDOUT->autoflush(1);
  $pipename = $_[0];        # Get pipe name from argument #1
  if (! -p $pipename)
      {$iret = system("mkfifo $pipename"); }
  open(FIFO, "+<".$pipename ) || die  $!;

  spi("MD");
  spi("PIPE");              # Opens output pipe in SPIDER
  spi($pipename);           # Name of pipe

  return FIFO;              # Returns pipe id
  }


	
sub getreg                 # Gets SPIDER register variable value from pipe -------
  {                        # INPUT: register number (argument # 1)
  local($reg,$regno,$regval);
  $\ = "\n";               # Set output record separator

  $reg = $_[0];            # Get register from argument #1 to this subroutine

  spi("PI REG");           # Tell SPIDER to put register variable value on pipe
  spi("$reg");             # Register variable wanted

  #{print STDERR "Reading from pipe now for reg: $reg \n" ;}
  #($regno,$regval) = unpack("Lfc",<FIFO>); # Read register value from pipe
  #{print STDERR "Read string :$t0:$t1:$t2:$t3:$t4:$t5:$t6:$t7:$t8:$t9:$t10:\n" ;}

  # FOR SGI USE FOLLOWING  LINE
  #($regval) = unpack("f",<FIFO>);       # Read register value from pipe

  # FOR INTEL BASED LINUX WITH SPIDER COMPILED WITH: -byteswapio USE FOLLOWING 5 LINES
  ($t0,$t1,$t2) = unpack("Nff",<FIFO>);   # Read register value from pipe
  if ($first)              # I do not know why this is different or what t0 & t1 are
     { $regval = $t1; $first = 0; }
  else
     { $regval = $t2; }

  #($t0,$t1,$t2,$t3,$t4,$t5) = unpack("Nffffc",<FIFO>); # Read register value from pipe
  #{print STDERR "Read string :$t0:$t1:$t2:$t3:$t4:$t5:\n" ;}
  #{print STDERR "Returning:  $regval \n" ;}

  return $regval;
  }
	
	
sub spi
  {   # Pipes argument to SPIDER after variable substitution ---------------
  #print STDERR "GOT:" . $_[0];
  local($string,$ret);
  s/"/\\"/g;                      # Substitutes: \" for: " 
  $string = $_[0];
  $ret    = eval qq/"$string"/;
  print $ret;                     # This sends string down the pipe
  }

 
