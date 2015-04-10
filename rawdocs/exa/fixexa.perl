#!/usr/local/bin/perl
#
# fixexa.perl
# fixexa
#
#USAGE: ls -1 *html >! junk
#        perl  fixexa.perl < junk > junk1

$\ = "\n";              # set output record separator

line: while (<>) 
    {
    $exafile = $_;

    $oper = uc $exafile;
    $oper =~s/\.HTML//g;
   #print ' Oper:'. $oper.""; 
   
   # Open  source file 
   open(EXAFILE,"<$exafile") ||  die "Can't open $exafile \n";
   #print ' opened:'. $exafile.""; 

   # Open temp output file
   open(EXAOUT,"+>jnkdir/$exafile") ||  die "Can't open jnkdir/$exaout \n";
   
   while (<EXAFILE>)
       {
       $line = $_;
       
       $line =~s/\<\/span\>//g;     # Remove <span>
       
       $line =~s/\t/  /g;           # Collapse tabs
       $line =~s/;/; &nbsp;  /g;           # Add space after ;
       
       # New classname for  for operation line      
       $line =~s/\<div class="com"\>/<div class="oplin">/g;
 
       # make class for operation in use      
       $line =~s/OPERATION.*\<td.*class="res/OPERATION  <\/td> <td class="opres/g;
              
       $line =~s/\<td\>\<span class="op"\>/<td class="op">/g;
       $line =~s/\<td\>\<span class="pr"\>/<td class="pr">/g;
       $line =~s/\<td\>\<span class="res"\>/<td class="res">/g;
       $line =~s/\<td\>\<span class="com"\>/<td class="com">/g;
        
       $line =~s/\<th\>\<span class="lab_in"\>/<th class="lab_in">/g; 
       $line =~s/\<th\>\<span class="lab_out"\>/<th class="lab_out">/g; 

       $line =~s/\<td\>\<span class="img_in"\>/<td class="img_in">/g; 
       $line =~s/\<td\>\<span class="img_out"\>/<td class="img_out">/g; 

       $line =~s/\<td\>\<span class="nam_in"\>/<td class="nam_in">/g; 
       $line =~s/\<td\>\<span class="nam_out"\>/<td class="nam_out">/g; 
      
       $line =~s/\<td\>\<span class="op"\>/<td class="op">/g; 
       $line =~s/\<\/span\>//g; 
       
       
       $line =~s|\.\./|./images/|g; 
       
       printf EXAOUT "%s", $line ;
       if ($line =~ m/\<body\>/)
           {
           printf EXAOUT '<br />' ;
           printf EXAOUT '<br />' ;
           printf EXAOUT "<h2>Usage Example - Operation: $oper</h2>" ;
           printf EXAOUT '<hr />' ;
	   }
	   
      }
   close(EXAOUT)         || die "Can't close jnkdir/$exafile \n";
   close(EXAFILE)        || die "Can't close $exafile \n";
   }

