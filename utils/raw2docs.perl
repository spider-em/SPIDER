#!/usr/bin/perl
#
# SOURCE: /usr16/software/spider/utils/raw2docs.perl
#
# CHANGES:
#         Added exa buttons and links            ArDean Leith  Jan. 2013
#         Moved raw man dir                      ArDean Leith  Mar. 2015
#         usr16                                  ArDean Leith  Apr. 2017
#
# PURPOSE: Puts buttons, header, and trailer on html documents. If -wadsworth
#          copies files to external Web site
#          Most HTML files (except if links) under the directory docs have 
#          their replacements done. The $rel_path is suitably modified  enroute..
#
# Locations and definitions previously in: htdefs.ph;
$spider_root   = qq($ENV{'SPIDER_DIR'});
print " SPIDER source dir: $spider_root \n";

$raw_docs_dir  = qq($spider_root/rawdocs);   # Unheaderized HTML docs directory
$docs_dir      = qq($spider_root/docs);      # Final headerized HTML docs directory
$raw_mans_dir  = qq(jnk_raw_mans);           # Unheaderized man pages directory created by: old2raw.perl
$html_man_dir  = qq($spider_root/docs/man);  # Final HTML manual pages directory
$exa_path      = "exa";                      # Examples dir

$wads_docs_dir = 'spider-stage:/export/apache/vhosts/spider.wadsworth.org/htdocs/spider_doc/spider/docs/';    

# For Wadsworth copy command compress, update, preserve Exec., preserve time,recursive,verbose
$sendit = qq(rsync -zuEtrv);

# --------------------------- Main routine ---------------------------

#print "sendit: $sendit \n";

$fnum      = 0;
$wadsworth = 0; 
if (@ARGV[0] =~ /wadsworth/)
   {  # Special destination for Wadsworth site
   $wadsworth = 1;
   print (" For Wadsworth WWW site \n");

   # Retarget $docs_out_dir & $docs_dir for external files
   $docs_dir     = qq(jnk_raw_www_docs);    # Created in utils directory
   $docs_out_dir = $wads_docs_dir;
   system("mkdir -p $docs_dir $docs_dir/man $docs_dir/exa/images") ;
   }


$rel_path = ".";
while (<$raw_docs_dir/*.html>)
   {
   # Remove: $raw_docs_dir from $fn 
   s/$raw_docs_dir\///e;
   $fn = $_;

   #print " Opening fn: $raw_docs_dir/$fn \n";
   
   open(fro, "< $raw_docs_dir/$fn");

   # Record input file's atime & mtime
   ($atime, $mtime) = (stat("$raw_docs_dir/$fn"))[8,9];

   #print " Opening fn: $docs_dir/$fn \n";
   open(too, "> $docs_dir/$fn")|| 
      die "Can't open: $docs_dir/$fn \n";

   # Search for 'nobuttons' flag in source html comment
   system("grep '<!-- NOBUTTONS -->' $raw_docs_dir/$fn > /dev/null");
   #system("grep 'garbage' $raw_docs_dir/$fn");
   
   $wantbuttons = 0;
   $wantbuttons = 0;
   if ($? == -1) 
      { print " Failed to execute: $! \n"; } 
   elsif ($? & 127) 
      { printf " Grep died with signal: %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without'; } 
   else 
      { 
      #printf " Child exited with value: %d \n", $? >> 8; 
      $wantbuttons = $? >> 8;
      }
   $wanttrailer = $wantbuttons;
      
   #print " Grep of: $raw_docs_dir/$fn Returned: $?  Wantbuttons: $wantbuttons \n";
   
   #if ($wantbuttons)
      #{ print " Want to add buttons and trailer to: $fn \n";}
      
   if ($wantbuttons == 0) 
      { print " Not adding buttons or trailer to: $fn\n";}

   if ($fn=~/^[A-Z]/) 
      {  # Filename starts with uppercase letter
      while ($linn=<fro>)  {print too $linn;}
      }
   else 
      {&installit;}
   close(fro);
   close(too);
   #print " Created: $docs_dir/$fn \n";

   # Alter the output files atime & mtime so same as input
   utime($atime, $mtime, "$docs_dir/$fn");
   $fnum++;
   }   

print " Created: $fnum  HTML SPIDER doc files in: $docs_dir \n";
print " ------------------------  \n";


# Now go one level deeper into the HTML manual directory
$rel_path = "..";
$fnum     = 0;

# Start with $raw_mans_dir
# Could not use earlier strategy - get "Arguments too long"  error msg.

opendir(raw_dir,"$raw_mans_dir") || 
    die "Can't open the unheaderized manual directory: $raw_mans_dir \n";
@files = readdir(raw_dir);		# Read in file names
closedir(raw_dir);

manf:
foreach $fn (@files)
   {			# Drop ma1 here if you want to
   # Only process .html files			
   unless ($fn =~/.*\.html/)
      {next manf;}

   #print " Opening fn: $raw_mans_dir/$fn \n";
   open(fro, "< $raw_mans_dir/$fn");

   # Record input file's atime & mtime
   ($atime, $mtime) = (stat("$raw_mans_dir/$fn"))[8,9];

   open(too, "> $docs_dir/man/$fn") || 
       die "Can't open: $docs_dir/man/$fn \n";;
   &installit;
   close(fro);
   close(too);

   # Alter the output files atime & mtime so same as input
   utime($atime, $mtime, "$docs_dir/man/$fn");
   $fnum++;
   }
    
print " Created: $fnum  HTML SPIDER man files in: $docs_dir/man \n";
print " ---------------------- \n";


# Process the exa example directory
$rel_path="..";
$fnum=0;

opendir(raw_dir,"$raw_docs_dir/$exa_path") || 
    die "Can't open the exa directory: $raw_docs_dir/$exa_path \n";
#print "Opened: $raw_docs_dir/$exa_path \n";

@files = readdir(raw_dir);		# Read in exa file names
@files = sort(@files);

closedir(raw_dir);

exaf:
foreach $fn (@files)
   {
   # Only process .html files			
   unless ($fn =~/.*\.html/)    
      {next exaf;}

   #print " Processing: $raw_docs_dir/$exa_path/$fn \n";
   
   open(fro, "< $raw_docs_dir/$exa_path/$fn");  # Open exa file

   # Record input file's atime & mtime
   ($atime, $mtime) = (stat("$raw_docs_dir/$exa_path/$fn"))[8,9];

   open(too, "> $docs_dir/exa/$fn") || 
       die "Can't open: $docs_dir/exa/$fn \n";;
       
   &installit ("exa") ;
   
   close(fro);
   close(too);

   # Alter the output files atime & mtime so same as input
   utime($atime, $mtime, "$docs_dir/exa/$fn");
   $fnum++;
   }   # End of: foreach ...
   
 
print " Processed: $fnum  EXA html files in: $docs_dir/exa \n";
print " ---------------------- \n";

# -------------------------- installit subroutine --------------------

sub installit
   {
   #if ($_[0] ne "exa")   removed mar 2015 al
   #   {print too "<html> \n";}
   
   while(<fro>)
      {
      if ((/<\/HEAD/ || /<\/head/)) 
         {   # Insert scripts at end of HEAD
         #print too qq(<!-- Begin Scripts -->\n);
         $line = $_  ;

         if ($wantbuttons) 
            {# Insert buttons scripts
            $dum=&include_stuff("scripts.ph");
            }
	    
         #print too qq(<!-- End Scripts -->\n);
         print too $line;
         }

      elsif (/<BODY/ || /<body/) 
         {   # Beginning of body
         print too $_;

         if ($wantbuttons) 
            {# # Insert Buttons at beginning of body 
            print too qq(<!-- Begin Buttons -->\n);
            # Insert local buttons
            $dum=&include_stuff("activebuttons.ph");
          
            print too qq(<br><hr>\n);
            print too qq(<!-- End Buttons -->\n);
	    }   
         print too qq(\n);
         }
      elsif ((/<\/BODY>/ || /<\/body>/) && ($_[0] ne "exa"))
         {   # End of body 
         if ($wanttrailer)
	    { # Insert trailer before end of body
            print too qq(<!-- Begin Trailer -->\n);
            $trailer = &print_trailer;
            print too $trailer;

            print too qq(<!-- End Trailer -->\n);
	    }
         print too $_;
         }
      else 
         {			#Otherwise.. MOST LINES
         $opline = $_;
	 
         if (($_[0] eq "exa") && $_ =~ /- Operation:/)
	    {     # Add link to operation;
	    #print qq( Opline: $opline );
	    
            #if ($fn =~/^pk/) {print qq( Opline: $opline );}
	    
	    # Extract OPERATION (w/wo quotes)
            #m/- Operation: ('?\w+ ?+\w+ ?+\w+'?)/;
            m/- Operation: ('?\w+ ?\w* ?\w* *'*)/;
	    
	    $op = $1; $op =~ s/'//g;   # Remove quotes if any
	    
            #if ($fn =~/^pk/) {print qq( op: $op \n); }
	    	    
	    $opfil = $op;
	    #print qq(Op: $opfil \n);
	    
	    $opfil =~ s/ //g;         # Remove blanks
	    #print qq( Opfil: $opfil \n);
	    
            $opfil =~ tr/[A-Z]/[a-z]/;
	    #print qq( Opfil lc: $opfil \n);

	    $nlet = length($opfil);
	    
            while($nlet >= 2)
	       {           # Try shorter name
	       #print  " Checking for: $html_man_dir/$opfil.html \n"; 
		    
	       if ( -e "$html_man_dir/$opfil.html") 
                  { 
	          #print qq( VOne: $1 \n);
		  
	          $opline =~ s#$op#<a href="../man/$opfil.html">$op</a>#;
	          #print qq( New opline: $opline \n);
                  last;
                  }

	       $nlet  = $nlet -1;
	       $opfil = substr($opfil,0,$nlet);
               }
            } #end of: if (($_[0] eq "exa"
	    
         print too $opline;
         } #end of: else
	 
      } #end of: while(<fro>)
      
      if ($_[0] ne "exa") 
         {print too "</html>";}
         
   } #end of: installit


#---------------------------------- include_stuff -------------------------
sub include_stuff
   {
   open(FILE,$_[0]);
   ##print " input file: $_[0] \n";

   while (<FILE>)
      {
      s/\$rel_path/$rel_path/g;
      print too $_;
      }	
   }

#---------------------------------- print_trailer ------------------------
sub print_trailer
   {
   return qq(<hr> <small>&copy; <a href="$rel_path/copyright.html"> \
   Copyright Notice</a> /           &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  \
   Enquiries: <a href="mailto:spider\@wadsworth.org">spider\@wadsworth.org</a> </small></a>\n);
   }

#---------------------------------------------------------------------------

