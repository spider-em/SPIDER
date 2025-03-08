#!/usr/bin/env perl
#
# SOURCE: /home/dean/Software/spider-2020-2025/utils/file-redate.perl
#
# PURPOSE:  Edite all files in a directory to fix buttons & remove trailer
#
# CHANGES:     Author:    ArDean Leith Feb 2025
# USAGE:       /home/dean/Software/file-edit.perl 


#------- active

 #mkdir old_exa new_exa
 # cp -a spider-github-2025/docs/exa/*html old_exa
 # ./file-edit.perl
 # diff -aN -C 0 old_exa/acms.html  new_exa/acms.html
 # touch -r=spider-github-2020/README.md/new_exa/ac* 
 
 
 #$old_dir = qq(old_docs);      # Bad  input dir
 #$new_dir = qq(new_docs);      # Good output dir
 
 #$old_dir = qq(old_mans);    # Bad  input dir
 #$new_dir = qq(new_mans);    # Good output dir

 $old_dir = qq(old_exa);    # Bad  input dir
 $new_dir = qq(new_exa);    # Good output dir
 
 $tips    = qq(/tips/);
 $webadd  = qq(../../web/docs/web.html);
 $newweb  = qq(  <td><a href="https://spider-em.github.io/Web" id="web"> </a></td>); 
 $begin   = qq(!-- Begin Trailer -->)  ;     
 $copy    = qq(copyright.html)         ;  
 $notice  = qq(Copyright)              ;      
 $enq     = qq(Enquiries:)             ;
 $end     = qq(<!-- End Trailer -->)   ;

     	
# Actual program begins here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   my $fnum = 0;
     
   opendir(OLDDIR,$old_dir) || 
      die "Can't open directory: $old_dir \n";
   #print " Opened old dir: $old_dir \n";
   
   @all_files = readdir(OLDDIR);    # Read in file names

   allf:
   foreach $file (@all_files)
      {
      $oldfile = "$old_dir/$file" ; 
      # Check for non src files, do not edit them
      #print " File: $oldfile \n";

      if ($oldfile !~ /$\.html/ )    
         { next ; }    #print "Not doc: $file\n"; 
 
      open(FROM, " $oldfile ") || die "Can't open $oldfile \n";
      #print " Opened: $oldfile \n";

      $newfile = "$new_dir/$file";
      open(TOO, "> $newfile") || die "Can't write to $newfile \n";
      print " Opened: $newfile \n";
 
      @lines = <FROM>;   # Load whole file contents

      allines:
      while(@lines)
        {
        $line = shift(@lines);

        if ($line =~ /$tips/)    {next;} #"{ print  " Skip:  $line"; next;}
        if ($line =~ /$webadd/)
           { 
           #print " Sub:    $newweb \n" ;
           print  TOO "$newweb \n" ;
	   next;
	   }
	
        if ($line =~ /$begin/)   {next;} # { print  " Skip:  $line"; next;}
        if ($line =~ /$copy/)    {next;} # { print  " Skip:  $line"; next;}
        if ($line =~ /$notice/)  {next;} # { print  " Skip:  $line"; next;}
        if ($line =~ /$enq/)     {next;} # { print  " Skip:  $line"; next;}
        if ($line =~ /$end/)     {next;} # { print  " Skip:  $line"; next;}
	
        print TOO "$line";
        }
     # exit;
     $fnum++;
     #if ($fnum > 1) {exit;}
     }
     
     
# <td><a href="../../web/docs/web.html"               id="web">       </a></td>
# <!-- Begin Trailer -->
# <hr> <small>&copy; <a href="./copyright.html"> 
# Copyright Notice</a> /           &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  
# Enquiries: <a href="mailto:spiderem.software@gmail.com">spiderem.software@gmail.com</a> </small></a>
# <!-- End Trailer -->
