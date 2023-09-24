
-------------------------------------------------------

**SPIDER**  (**S**ystem for **P**rocessing **I**mage **D**ata from **E**lectron microscopy and **R**elated fields) 
is an image processing system for electron microscopy, especially usefull for single-particle reconstruction. 
[*SPIDER*](./spider.html) has been in use since 1978 and contains 
over 210,000 lines of Fortran code and 7400 files. The lastest release from July 2020 (which will probably be the 
final release) contains  support for reading/writing MRC files.  

This latest release is available [here](http://github.com/spider-em/SPIDER/releases). 

* Uses:
               
   1.  [Averaging of single particle objects](http://spider.wadsworth.org/spider_doc/spider/docs/philosophy.html) 
  
   2.  [Multivariate statistical analysis](http://spider.wadsworth.org/spider_doc/spider/docs/techs/MSA/index.html)
  
   3.  [3D reconstruction](http://spider.wadsworth.org/spider_doc/spider/docs/strategies.html)

   4.  [Electron tomography](http://spider.wadsworth.org/spider_doc/spider/docs/techs/lgstr/tomo/tomo.html)
                     
   5.  [Object segmentation from volumes](http://spider.wadsworth.org/spider_doc/spider/docs/techs/segment/segment.html)

        
* Features:
               
   1.   Interactive [command line interface](http://spider.wadsworth.org/spider_doc/spider/docs/user_doc.html).  

   2.   Hierarchical modular design for scripting.  

   3.   Graphical user interface, [*Web*](http://spider.wadsworth.org/spider_doc/web/docs/web.html), for 
        visualizing and interacting with images.  

   4.   Uses a [File format](http://spider.wadsworth.org/spider_doc/spider_doc/spider/docs/image_doc.html) 
        interchangeable with other electron microscopy imaging systems.  

   5.   Extensive [documentation](http://spider.wadsworth.org/spider_doc/spider/docs/documents.html) of 
                     [operations](http://spider.wadsworth.org/spider_doc/spider/docs/operations_doc.html) and
                     [techniques](http://spider.wadsworth.org/spider_doc/spider/docs/techniques.html).  
                  
   6.   Includes [source code](http://spider.wadsworth.org/spider_doc/spider/src) and executables for use on Linux. 
 
* History:

   1.   Originated in [1978](http://spider.wadsworth.org/spider_doc/spider/docs/spider78.html) by 
           [Joachim Frank](http://franklab.cpmc.columbia.edu/franklab) who
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

   3.   SPIDER has often been used in conjunction with other software e.g. the often  cited methods developed by Edward H. Egelman. See:  [The iterative helical real space reconstruction method: Surmounting the problems posed by real polymers](http://doi.org/10.1016/j.jsb.2006.05.015).
             
   4.   Maintained since 2007 by: ArDean Leith & Tanvir Shaikh. During these 16 years they added numerous 
        operations, procedures, parallelizations, and over 500 pages of documentation.  There have been no 
        significant contributors other than Leith, Shaikh, Kischenko, and Kennedy since 2007.

  5.    In 2023 we are happy to see that citations of the use of SPIDER continue to occur in the methodology sections of significant new research publications. 

-----------------------------------------------------

[Copyright under GPL].   

