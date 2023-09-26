
-------------------------------------------------------

**SPIDER**  (**S**ystem for **P**rocessing **I**mage **D**ata from **E**lectron microscopy and **R**elated fields) 
is an image processing system for electron microscopy, especially usefull for single-particle reconstruction. 
**SPIDER** has been in use since 1978 and contains 
over 210,000 lines of Fortran code and 7400 files. The lastest release from July 2020 (which will probably be the 
final release) contains support for reading/writing MRC files.  

This latest release is available [here](http://github.com/spider-em/SPIDER/releases). 

* Uses:
               
   1.  Averaging of single particle objects.
  
   2.  Multivariate statistical analysis.
   
   3.  3D reconstruction from projections.
   
   4.  Electron tomography.
   
   5.  Object segmentation from volumes.
       
* Features:
               
   1.   Interactive command line interface
      
   3.   Hierarchical modular design for scripting.  

   4.   Graphical user interface, [*Web*](http://github.com/spider-em/web) 
        for  visualizing and interacting with images.  

   6.   Uses a File format which is interchangeable with other electron microscopy imaging systems.  

   7.   Extensive documentation of  operations and techniques.  
                  
   8.   Includes source code and executables for use on Linux. 
 
* History:

   1.   Originated in by 
           [Joachim Frank](http://joachimfranklab.org) who
           shared the 2017 Nobel Prize in Chemistry for his pioneering cryo-em research done at the 
           [Wadsworth Center](http://www.wadsworth.org) in Albany, NY.
       
   2.  Contributors (in chronological order): 
          J. Frank,    B. Shimkin,   H. Dowse,       L. Miranda,   W. Goldfarb,  S. Basu,       R. Banerjee,
          C. Mannella, J. P. Bretaudiere, 
          M. V. Heel,  A. Verschoor, M. Radermacher, A. Leith,     J. M. Carazo, P. Penczek,    S. Sibal, 
          L. Odesanya, Y. H. Li,     M. Ladjadj,     Y. W. Chen,   K. R. Lata,   J. Zhu,        W. P. Liu,   B. Rath, 
          C. Yang,     B. Baxter,    R. Hegerl,      A. Frangakis, L. Joyeux,    Z. Huang,      R. J. Renka, 
          T. Shaikh,   J. Sengupta,  J. LeBarron,    N. Boisset,   H. Gao,       G. Kishchenko, J. M. Kennedy, A. Gunggong. 
                
          In addition SPIDER is dependent on the Fourier package [FFTW](http://www.fftw.org) by M. Frigo and SG Johnson.  

   3.   SPIDER has often been used in conjunction with other software e.g. the [iterative helical real space reconstruction methods](http://doi.org/10.1016/j.jsb.2006.05.015) developed by Edward H. Egelman. 
             
   4.   Maintained since 2007 by: ArDean Leith & Tanvir (Tapu) Shaikh. During these 16 years they have added numerous 
        operations, procedures, parallelizations, and over 500 pages of documentation.  There have been no 
        significant contributors other than Leith, Shaikh, Kischenko, and Kennedy since 2007.

  5.    In 2023 we are happy to see that citations of the use of SPIDER continue to occur in the methodology sections of significant new research publications. 

-----------------------------------------------------

[Copyright under GPL].   

