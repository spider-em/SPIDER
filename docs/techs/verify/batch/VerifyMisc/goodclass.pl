#!/usr/sbin/perl

if ($ARGV[0] eq "!" )
{
        $back = "\n";
        print "FILE NAME  :   goodclass.pl".$back;
        print "CREATE DATA:   07/2004".$back;
        print "PURPOSE    :   Used in particle picking procedure, creating 'goodclass' doc file from 'firstgoodclass'".$back;
        print "FORMAT     :   goodclass.pl  <firstgoodparticle> <goodclass> ".$back;
        exit;
}
open(A,"$ARGV[0]");
open(B,">$ARGV[1]");

$count = 0;

while($linea =<A>)
{
        chomp($linea);
        $begina = substr($linea,1,1);
        if ($begina eq ";")
        {
        }
        else
        {
         $count = $count + 1;
         $begina = substr($linea,1,4);
         $begina = $begina + 0;
         printf B (" %04d 2         %4d            1\n", $count,$begina);
        }
}
close(B);
print "END \n ";

