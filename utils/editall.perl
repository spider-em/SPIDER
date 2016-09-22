#!/usr/bin/perl
#
# SOURCE: /usr8/spider/utils/editall.perl
#
# PURPOSE:  Edit all files in a directory 
#
# CHANGES:         Author:    ArDean Leith June 2013
# USAGE:   /usr8/spider/utils/editall.perl 

$edit_dir = qq(/usr8/spider/rawdocs/exa/images);  # Files to be edited 
$new_dir  = qq($edit_dir/jnk);                    # Dir for edited files

	
# Actual program begins here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

opendir(EDIT_DIR,$edit_dir) 
   || die "Can't open $exa_dir\n";

@files = readdir(EDIT_DIR);	# Read in file names
@files = sort(@files);

closedir(EDIT_DIR);

$fnum = 0;			# Number of html files processed

foreach $file (@files)
  {
  $oldfile = "$edit_dir/$file";
  $newfile = "$new_dir/$file";
 
  #print " File: $file \n";
  # Check for .spi files, do not edit them
  if ($file !~ /$\.spi/ )    
     { 
     #print  " Not spi: $file\n";
     next;
     }

  unless (-T "$oldfile") 
     {print " No such File: $oldfile \n"; next;}	# Not a valid file
  
  #print " File: $file \n"; 

  open(FROM, " $oldfile ") || die "Can't open $oldfile \n";
  #print " Opened: $oldfile \n";

  open(TOO, "> $newfile") || die "Can't write to $newfile \n";
  #print " Opened: $newfile \n";

  @lines = <FROM>;   # Load whole file contents

  while(@lines)
     {
     $line = shift(@lines);

     if ($line =~ /BATCH/)
        { 
        $line =~ s/^-/ ; -/;
        #print  " Batch: $line";
	}
	
     if ($line =~ /^;/)
        { 
        $line =~ s/^;/ ;/;
        #print  " Semi: $line";
	}
	
     # Replace this filename with better filename
     $line =~ s/ip_001/sav_pp/g;
    	 
     print TOO "$line";
     }

  close(FROM);
  close(TOO);
  }
