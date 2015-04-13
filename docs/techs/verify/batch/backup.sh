ln -s $BACKUP Bak

date
echo "Copying batch files"
cp -a *.tar* Bak/
tar cf Bak/allbatch.tar */*spi 

echo "Copying top-level parameters"
cp -a *.dat Bak/

echo "Copying micrographs (slow)"
mkdir Bak/Micrographs
cp -av Micrographs/mic* Bak/Micrographs/

echo "Copying defocus parameters"
tar cf  Bak/powerspectra.tar Power_Spectra/*.dat

echo "Copying particle coordinates"
tar cf Bak/coords.tar Particles/coords/sndc*

echo "Copying particles-per-micrograph"
tar cf  Bak/orderpicked.tar Particles/*dat

echo "Copying alignment files"
tar cvf Bak/align.tar Alignment/align_01*

echo "Copying selected particles"
tar cf  Bak/select.tar Particles/good/sel*

echo "Copying initial-reconstruction data"
mkdir Bak/Reconstruction
cp -a Reconstruction/vol001.dat Bak/Reconstruction/
cp -a Reconstruction/combires.dat Bak/Reconstruction/
cp -a Reconstruction/resolution.dat Bak/Reconstruction/

echo; echo "Done"; date
