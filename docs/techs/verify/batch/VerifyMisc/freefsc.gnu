pxsize = 2.82
idim   = 130

set grid xtics
set ytics 0,.1,1
set xzeroaxis
set xlabel 'Resolution (A^-1)'
set ylabel 'FSC'

set title 'Cross-validated alignment, variable #particles'

plot [0:0.1][-.05:1.0675] \
'combires.dat'                    using ($3/pxsize):5        title 'Reconstruction 1' with lines, \
'../Reconstruction2/combires.dat' using ($3/pxsize):5        title 'Reconstruction 2' with lines, \
'docmaskfreefsc.dat'              using ($1/(pxsize*idim)):3 title 'mask profile'     with lines

pause -1 'Press RETURN to continue.'

