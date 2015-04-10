
C **********************************************************************
C                                                                      *
C NOYAU.F        ORIGINAL                                    31-JAN-86 *
C                DENDRO CALL ADDED TO PLOT DENDROGRAM        NOV 86 al *
C                DENDROGRAM DOC FILE BUGGY REMOVED           DEC 03 al *
C                CUTOFF CHANGED                              FEB 04 al *
C                FORMATTING CHANGED                          DEC 07 al *
C                COSMETIC OUTPUT CHANGES                     DEC 08 al * 
C                DOC FILE *                                  MAY 09 al *
C                DENDRO REWRITE                              MAY 09 al *
C                KV ADDED                                    NOV 11 al *
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
C*---------------------------------------------------------------------*         
C                                                                      *         
C PURPOSE:  CLASSIFICATION ACCORDING TO FACTORIAL COORDINATES          *
C                                                                      *
C INPUT  ...  DATA AND WEIGHTS IN IMC FILE ON: LUNI (FORMATTED)        *         
C      N2DIM     = SUP (NKLA, 2*NCLAS)                                 *         
C      KDIM      = SUP (NUMIM, NCLAS**NBASE, 2*NKLA - 1)               *         
C      MDIM_8    = SUP (NKLA*(NKLA-1)/2)                               *         
C      LDIM      = SUP (NKLA)                                          *         
C      L2DIM     = SUP (2*NKLA - 1)                                    *         
C      NFAC      = NUMBER OF COORDINATES IN RECORD FOLLOWED BY WEIGHT  *         
C      KFAC      = NUMBER OF COORDINATES USED FOR CALCULATION          *         
C      NUMIM     = NUMBER OF OBJECTS TO BE CLASSIFIED                  *         
C      NBASE NITER ... PARAMETERS DEFINED IN SUBROUTINE PARST          * 
C      NCLAS NKLA  ... PARAMETERS DEFINED IN SUBROUTINE PARST          *
C                                                                      *
C WORKING ARRAYS:  GT(*,*) D(*) KLAS(*) CI(*) U(*) JV(*) JW(*)         * 
C                  U(), JV(), JW() EQUIVALENCED TO D()                 *         
C                                                                      *         
C  CALL TREE:   SCLASSI - NOYAU - PARST  - RETIR                       *
C                  |        |            - CLAST                       *
C                  |        |            - STABK - SHELK               *
C                  |        |                                          *
C                  |      NOYAU - DEUCL                                * 
C                  |      NOYAU - CHAVA                                *
C                  |      NOYAU - DENDRO - DENLST                      * 
C                  |      NOYAU - COUPE                                *
C                  |                                                   *
C               SCLASSI - RGRI                                         *    
C                                                                      *         
C*---------------------------------------------------------------------*         

      SUBROUTINE NOYAU(N2DIM,KDIM,MDIM,LDIM,L2DIM,NFAC,KFAC,NUMIM,          
     &                 KV,NBASE, NITER, NCLAS, NKLA, IDI,                
     &                 KLAS, D, U, JV, JW, IDK, CI, GT,               
     &                 NUM,LA,LB,IV,IW,V,W,NT,VAL,PK,NO,LUNI,LUNK,
     &                 LUNDOC)

         
        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 
        !COMMON /UNITS/LUN,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

        INTEGER *8 :: MCARD_8
        INTEGER    :: N2DIM,KDIM,LDIM,L2DIM,NFAC,KFAC,NUMIM
        INTEGER    :: KV(KFAC)
        INTEGER    :: NBASE,NITER,NCLAS,NKLA
        INTEGER    :: IDI(NUMIM),KLAS(NUMIM)
        REAL       :: GT(N2DIM,KFAC),D(MDIM),CI(NFAC),U(KDIM)
        INTEGER    :: JV(KDIM), JW(KDIM), IDK(NKLA)
        INTEGER    :: NUM(LDIM),LA(LDIM),LB(LDIM),NT(L2DIM),               
     &                NO(L2DIM),IV(LDIM),IW(LDIM)                    
        REAL       :: VAL(L2DIM),PK(L2DIM),W(LDIM),V(LDIM)                    
        INTEGER    :: LUNI,LUNK, LUNDOC

        IRTFLG = 0

C       AGGREGATION AROUND MOBILE CENTERS AND STABLE CLUSTERING             
        CALL PARST(N2DIM,KDIM,NFAC,NUMIM,KFAC,NBASE,NITER,NCLAS,NKLA,            
     &                       KLAS,JW,U,CI,JV,GT, KFIN,LUNI,LUNK,IDI,KV)              
                                                                               
C       HIERARCHICAL CLASSIFICATION OF THE CLUSTER GRAVITY                  
C       CENTERS ACCORDING TO THE WARD'S VARIANCE CRITERION                  

C       IMPORTANT CHANGE: LIMIT NUMBER OF CLUSTERS USED IN HAC. 8/25/86

        WRITE(NDAT,*)' CLASS ASSIGNMENT FOR EACH IMAGE:'
        WRITE(NDAT,'(10I5)') (KLAS(K),K=1,NUMIM)

        IRANGE1 = MINVAL(JV(1:NKLA))  ! JV IS OCCUPANCY LEVEL
        IRANGE2 = MAXVAL(JV(1:NKLA))

        IF (NDAT .NE. NOUT) WRITE(NOUT,91)IRANGE1,IRANGE2,NKLA
        WRITE(NDAT,91) IRANGE1,IRANGE2, NKLA
91      FORMAT(/'  CLASS OCCUPANCY: ',I7,'....',I7,'   CLASSES: ',I6,/)

        ILEVL = IRANGE1
        CALL RDPRI1S(ILEVL,NOT_USED,
     &      'OCCUPANCY LEVEL FOR CLASS CUTOFF (<CR> = NO CUTOFF)',
     &      IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (ILEVL .NE. IRANGE1) THEN
           DO I=1,NKLA
              IF (JV(I) < ILEVL) THEN
                 NKLA = I - 1
                 EXIT
              ENDIF
            ENDDO
         ENDIF

        IF (NDAT .NE. NOUT) WRITE(NOUT,92)NKLA
        WRITE(NDAT,92) NKLA
92      FORMAT(/,'  USING: ',I6,' CLASSES.')

C       READS FROM _IMC TO CONSTRUCT MONO-INDEXED TABLE OF DISTANCES 
C       BETWEEN NKLA CLASSES IN THE EUCLIDEAN SPACE CHARACTERIZED 
C       BY THE KFAC FIRST FACTORIAL COORDINATES. READS FROM: LUNI                             *         
C       OUTPUTS: D, GT & PK

        IF (NKLA > KFIN) NKLA = KFIN                              
        MCARD = NKLA * (NKLA-1) / 2 
                             
        CALL DEUCL(N2DIM,NKLA,MCARD,NUMIM,KFAC,NFAC,KLAS,
     &             D,GT,PK,CI,LUNI)

C       CHAVA OVERWRITES AND THEN READS FROM: LUNK
        JFIN    = 2 * NKLA - 1                                
        MCARD_8 = MCARD
        CALL CHAVA(NKLA,MCARD_8,JFIN,D,PK,VAL,LA,LB,NT,NO,LUNK,5)                  
                                          
C       SAVES CLASSIFICATION TO FILE BY OVERWRITING: LUNK                                                                                                                     
        CALL REW(LUNK,0)                                                         
        WRITE(LUNK) NUMIM, NFAC, NKLA, KFAC                 
        WRITE(LUNK) (KLAS(I), I=1,NUMIM), 
     &              (IDI(I), I=1,NUMIM),                   
     &              (PK(L),L=1,NKLA), 
     &              ((GT(L,J),L=1,NKLA),J=1,KFAC)              
                                                                               
C     GENERATION OF IDENTIFIERS, ADDED 1/2001 al SEEMS TO BE
C     NEEDED IN DENDRO BUT REMOVED BY pp SOMETIME ??
      DO I = 1,NKLA
         IDK(I) = I
      ENDDO
       
#ifdef NEVER   
      write(0,*) ' l2dim,nkla,jfin: ',l2dim,nkla,jfin
#endif
          
C     DRAW CLASSIFICATION TREE OF THE NKLA CENTERS                             
      CALL DENDRO(NKLA, JFIN, VAL, LA, LB, PK, IDK, KLAS,NUMIM,IDI,
     &              .TRUE., NO,NUM,NT,IV,IW,V,W)             

C     SUCCESSIVE TRUNCATIONS OF THE CLASS. 
      DO KPRO =2,NKLA-1                                                   
         CALL COUPE(NUMIM,NKLA,KPRO, PK,LA,LB,IV,KLAS,NT, IW)                     

C        TRUNCATED TREE OUTPUT APPENDED TO LUNK 
         WRITE(LUNK) (IV(J), J=1,NKLA)                                           
      ENDDO

      END






