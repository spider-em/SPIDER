; rescales RNA histograms for viewing comparison using gnuplot
vm
\rm hist_rna_calc.dat hist_rna.dat


;X20 = 297   ; HISTOGRAM NO. 1  MAX : output/hist_rna.dat 
X21 = 0     ; HISTOGRAM NO. 1  MIN

;X30 = 1129   ; HISTOGRAM NO. 2  MAX    :histo-RNAcore170mu083-g001.dat
X31 =  0     ; HISTOGRAM NO. 2  MIN

ud n, x32,x33
output/hist_rna.dat 

x20 = -99      ; HISTOGRAM NO. 1  MAX 

do lb1 x40 = 1,x32 
   ud ic, x40,x41,x42
   output/hist_rna.dat 
   if (x42 .gt. x20) then 
      x20 = x42
   endif
lb1


do lb3 x40 = 1,x32 
   ud ic, x40,x41,x42
   output/hist_rna.dat 
   x42 = x42/x20
   
   sd x40,x41,x42
   hist_rna.dat
   
lb3

sd e
hist_rna.dat


ud n, x52,x53
output/histo-RNAcore120mu090-g003.dat


x30 = -99      ; HISTOGRAM NO. 1  MAX 

do lb2 x40 = 1,x52 
   ud ic, x40,x41,x42
   output/histo-RNAcore120mu090-g003.dat
   if (x42 .gt. x30) then 
      x30 = x42
   endif
lb2


do lb4 x40 = 1,x52 
   ud ic, x40,x41,x42
   output/histo-RNAcore120mu090-g003.dat
   x42 = x42/x30
   
   sd x40,x41,x42
   hist_rna_calc.dat
   
lb4

sd e
 hist_rna_calc.dat




en


