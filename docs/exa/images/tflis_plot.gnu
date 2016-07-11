set ylabel "Amplitude"
set xlabel "Radius (1/A)"
set title "Amplitude vs Radius"
plot  \
"tmp-tflis-nohead" using  4:6 title "-Straight CTF" with lines, \
"tmp-tflis-nohead" using  4:7 title "Flipped CTF"  with lines, \
"tmp-tflis-nohead" using  4:8 title "Trapped CTF"  with lines
