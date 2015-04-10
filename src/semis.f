C++*********************************************************************
C
C SEMIS.F
C                
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2013  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* General Public License for more details.                           *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C                                                                      *
C  PURPOSE: CLASSIFICATION OF OBJECTS ACCORDING TO THEIR               *
C           FACTORIAL COORDINATES                                      *         
C                                                                      *         
C   ROUTINES CALLED - NOYAU - PARST - RETIR - SEN3A                    *         
C                                   - CLAST                            *         
C                                   - STABK - SHELK                    *         
C                           - DEUCL                                    *         
C                           - CHAVA                                    *         
C                           - ARBRE                                    *         
C                           - COUPE                                    *         
C                                                                      *         
C **********************************************************************

      SUBROUTINE SEMIS(Q, MAXMEM, NGUS, NGRI, MODE,                           
     &               ICARD, NFAC, KFAC, KV, NBASE, NITER, NCLAS, NKLA)          

      DIMENSION  Q(MAXMEM), KV(16)
                                              
      COMMON  / ENSOR / LEC,IMP                                                 

      WRITE (IMP,2000)                                                          
 2000 FORMAT (23X,' ** STEP: CLASSY **',///,
     &         1X,72('-'))                     

      WRITE(IMP, 2100) (KV(K), K=1, KFAC)                                      
 2100 FORMAT (/,' SPECIFICATIONS FOR: CLASSY',/,                       
     &          '   FACTORS USED :',/' ',4X, 16I4)
                                         
      WRITE(IMP, 2150) NBASE, NITER, NCLAS, NKLA                               
 2150 FORMAT(/,' ',2X,'NBASE=',I6,5X,'NITER=',I6,5X,                          
     &                'NCLAS=',I6,5X,'NKLA =',I6)
                                              
      IF (KFAC <= 0) THEN
         CALL ERRT(102,'KFAC IS <= 0',KFAC)                              
         RETURN
      ENDIF

C     MEMORY SEGMENTATION                                                  

      N2DIM = MAX (NKLA, 2*NCLAS)                                              
      KDIM  = MAX (ICARD, NCLAS**NBASE)                                        
      KDIM  = MAX (KDIM,  2*NKLA - 1)                                          
      MDIM  = NKLA*(NKLA-1) / 2                                                 
      LDIM  = NKLA                                                              
      L2DIM = 2*NKLA - 1                                                        

      NKLAS = 1 + ICARD                                                         
C     NKLAS = 1                                                                 
      ND    = NKLAS + ICARD                                                     
      NU    = ND   + MDIM                                                         
      NJV   = NU   + KDIM                                                         
      NJW   = NJV  + KDIM                                                        
      NIDK  = NJW  + KDIM                                                        
      NCI   = NIDK + NKLA                                                       
      NGT   = NCI  + NFAC                                                        
      NNUM  = NGT  + N2DIM*KFAC                                                  
      NLA   = NNUM + LDIM                                                       
      NLB   = NLA  + LDIM                                                        
      NIV   = NLB  + LDIM                                                        
      NIW   = NIV  + LDIM                                                        
      NV    = NIW  + LDIM                                                        
      NW    = NV   + LDIM                                                         
      NNT   = NW   + LDIM                                                         
      NVAL  = NNT  + L2DIM                                                       
      NPK   = NVAL + L2DIM                                                      
      NNO   = NPK  + L2DIM                                                       
      NFIN  = NNO  + L2DIM
                                                       
      IF (NFIN > MAXMEM) THEN
         WRITE (IMP,2500)  MAXMEM, NFIN                                           
 2500    FORMAT (/,' UNLABELED COMMON BLOCK MEMORY AVAILABLE:',I10,5X,
     &              'YOU NEED: ',I11,/)
         CALL ERRT(6,'SEMIS',NDUM)                              
         RETURN
      ENDIF

C     CALL PRINCIPAL ROUTINE FOR CALCULATIONS                             
      CALL NOYAU (N2DIM,KDIM,MDIM,LDIM,L2DIM,                                   
     &             NFAC,KFAC,ICARD,KV,NBASE,NITER,NCLAS,NKLA, Q(1),             
     &  Q(NKLAS),Q(ND),Q(NU),Q(NJV),Q(NJW),Q(NIDK),Q(NCI),Q(NGT),               
     &  Q(NNUM),Q(NLA),Q(NLB),Q(NIV),Q(NIW),Q(NV),Q(NW),                        
     &  Q(NNT),Q(NVAL),Q(NPK),Q(NNO) ,NGUS,NGRI, MODE)
                          
      WRITE (IMP,2600)                                                          
 2600 FORMAT (//,' ',72('-'),//,
     &           ' ',5X,' ** END OF STEP:  CLASSY **',//,
     &           ' ',72('-'))
                     
       END                                                                    
