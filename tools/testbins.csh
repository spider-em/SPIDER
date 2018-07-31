# SOURCE:  tools/testbins.csh
# PURPOSE: Test python tools

 unsetenv PYTHONPATH
 cd $SPIDER_DIR/tools/tst

 #  # Use:  bin-python/use-python2.7 -------------------------------------
 
 # montagefromdoc   Display montage using doc file numbers  WORKS!!
 $SPIDER_DIR/tools/bin/pytest 

 # montagefromdoc   Display montage using doc file numbers  WORKS!!
 $SPIDER_DIR/tools/bin/montagefromdoc $SPIDER_DIR/tools/tst/sel_micrograph.tst  $SPIDER_DIR/tools-try/tst/pw_avg* &

 # qview            View SPIDER image                       WORKS!!
 $SPIDER_DIR/tools/bin/qview           $SPIDER_DIR/tools/tst/savavg000.tst                           
 $SPIDER_DIR/tools/bin/qview           $spid/save_face_16.dat 
 $SPIDER_DIR/tools/bin/qview           $SPIDER_DIR/tools/tst/bpass.jpg             
 $SPIDER_DIR/tools/bin-python/python   $SPIDER_DIR/tools/tst/qview.py  prj001/classavg002.tst 

 # spiconvert        Convert to/from SPIDER                  WORKS!!
 $SPIDER_DIR/tools/bin/spiconvert      $SPIDER_DIR/tools/tst/savavg000.tst   $SPIDER_DIR/jnk.tif
 display $SPIDER_DIR/jnk.tif &         (use Imagemagick to display)

 # montage           Display montage of files               WORKS!!
 $SPIDER_DIR/tools/bin/montage         $SPIDER_DIR/tools/tst/savavg***.tst &
 $SPIDER_DIR/tools/bin/montage         $SPIDER_DIR/tools/tst/savavg000.tst &

 # verifybyview      Show montage and select images        WORKS!!
 $SPIDER_DIR/tools/bin/verifybyview    $SPIDER_DIR/tools/tst/prj001

 # mkfilenums        Make list of file numbers              WORKS!!   
 $SPIDER_DIR/tools/bin/mkfilenums      $SPIDER_DIR/tools/tst/jnk_outdocfile.tst $SPIDER_DIR/savavg***.tst
 cat jnk_outdocfile.tst

 # classavg                                                WORKS!!
 $SPIDER_DIR/tools/bin/classavg        $SPIDER_DIR/savavg***.tst 
 # Select: savavg***.tst

 # pyplot            Plots doc file                          WORKS!!
 $SPIDER_DIR/tools/bin/pyplot          $SPIDER_DIR/tools/tst/fscdoc_m_02.tst 

 # xplor             Shows dir contents                      WORKS!!
 $SPIDER_DIR/tools/bin/xplor    
 # Select: a directory to view

 # viewstack         Sequentially step thru stack           WORKS!!
 $SPIDER_DIR/tools/bin/viewstack       $SPIDER_DIR/tools/tst/savribstk.tst & 

 # ctfcircle         Draws circle on image                     WORKS!!
 $SPIDER_DIR/tools/bin/ctfcircle       $SPIDER_DIR/tools/tst/savavg000.tst

 # Use:  bin-python/use-python2.5 kludge to overcome lack of PMW- Blt extension -------------------------
 #       which is needed by these tools   (python2.5 fails with error:: name "::blt::graph) 

 # ctfdemo           Study effect of ctf parameters         
 $SPIDER_DIR/tools/bin/ctfdemo 
 
 # ctfmatch          Find ctf equation                                          
 $SPIDER_DIR/tools/bin/ctfmatch        $SPIDER_DIR/tools/tst/ctf*   &                         

 # ctfgroup                                         Has default input???                       
 $SPIDER_DIR/tools/bin/ctfgroup                                     

 #emancoords2spiderdoc                                        
 $SPIDER_DIR/tools/bin/emancoords2spiderdoc    

 # scatter.py         Plots doc file?              No output
 $SPIDER_DIR/tools/bin/scatter fscdoc_m_03.tst              

 # Starts --------------------- Need test data, etc ----------
  
 $SPIDER_DIR/tools/bin/binarytree             Syntax: tree.py node_img001.ext {selectfile.ext max_depth margin_width canvas_width}

 $SPIDER_DIR/tools/bin/xmippsel2spiderdoc 

 import sys
 print '\n'.join(sys.path)
