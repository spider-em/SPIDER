; rescales RNA/protein histograms for viewing comparison using gnuplot
vm
\rm hist_rna_calc hist_pro_calc


;X20 = 297   ; HISTOGRAM NO. 1  MAX : output/hist_rna 
X21 = 0     ; HISTOGRAM NO. 1  MIN

;X30 = 1129   ; HISTOGRAM NO. 2  MAX    :histo-RNAcore170mu083-g001
X31 =  0     ; HISTOGRAM NO. 2  MIN

ud n, x32,x33
output/hist_rna 

x20 = -99      ; HISTOGRAM NO. 1  MAX 

do lb1 x40 = 1,x32 
   ud ic, x40,x41,x42
   output/hist_rna 
   if (x42 .gt. x20) then 
      x20 = x42
   endif
lb1


do lb3 x40 = 1,x32 
   ud ic, x40,x41,x42
   output/hist_rna 
   x42 = x42/x20
   
   sd x40,x41,x42
   hist_rna_calc
   
lb3

sd e
hist_rna_calc


ud n, x52,x53
output/hist_pro


x30 = -99      ; HISTOGRAM NO. 1  MAX 

do lb2 x40 = 1,x52 
   ud ic, x40,x41,x42
   output/hist_pro
   if (x42 .gt. x30) then 
      x30 = x42
   endif
lb2


do lb4 x40 = 1,x52 
   ud ic, x40,x41,x42
   output/hist_pro
   x42 = x42/x30
   
   sd x40,x41,x42
   hist_pro_calc
   
lb4

sd e
 hist_pro_calc




en


