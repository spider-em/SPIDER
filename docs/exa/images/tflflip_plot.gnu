set title "CTF Plots"
set xzeroaxis
set xlabel "Resolution, Angstroms^-1"
set ylabel "Transfer"
plot [0:0.3730][-1.1:1.4] \
"tflflip_doc.dat" using 7:3 title "flipped CTF" with lines, \
"tflflip_doc.dat" using 7:4 title "straight CTF" with lines, \
"tflflip_doc.dat" using 7:5 title "trapped CTF" with lines
