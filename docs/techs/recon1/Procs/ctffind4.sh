#!/bin/csh -x
#
#   ctffind2
#   input micrograph in MRC format
#   output power spectrumfile mic_number
#   CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]
#   Box, ResMin[A], ResMax[A], dFMin[A], dFMax[A], FStep
#
#               The output image file to check the result of the fitting
#               shows the filtered average power spectrum of the input
#               image in one half, and the fitted CTF (squared) in the
#               other half. The two halfs should agree very well for a
#               sucessfull fit.
#
#               CS: Spherical aberration coefficient of the objective in mm
#               HT: Electron beam voltage in kV
#               AmpCnst: Amount of amplitude contrast (fraction). For ice
#                        images 0.07, for negative stain about 0.15.
#               XMAG: Magnification of original image
#               DStep: Pixel size on scanner in microns
#               Box: Tile size. The program devides the image into square
#                    tiles and calculates the average power spectrum. Tiles
#                    with a significatly higher or lower variance are 
#                    excluded; these are parts of the image which are unlikely
#                    to contain useful information (beam edge, film number 
#                    etc). IMPORTANT: Box must have a value of power of 2.
#               ResMin: Low resolution end of data to be fitted.
#               ResMaX: High resolution end of data to be fitted.
#               dFMin: Starting defocus value for grid search in Angstrom. 
#                      Positive values represent an underfocus. The program
#                      performs a systematic grid search of defocus values 
#                      and astigmatism before fitting a CTF to machine 
#                      precision.
#               dFMax: End defocus value for grid search in Angstrom.
#               FStep: Step width for grid search in Angstrom.
#
#
# EDIT FOLLOWING LINE FOR LOCATION OF: ctffind4.exe ON YOUR SYSTEM
# IF IT IS NOT ALREADY IN YOUR PATH

time /home/tapu/local/ctffind-4.0.13/bin/ctffind << eof
$1
$2
$3
$4
$5
$6
$7
$8
$9
$10
$11
$12
$13
$14
eof
#
