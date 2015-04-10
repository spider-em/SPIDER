
C **********************************************************************
C                                                                      *
C SCLASSI.F                                                            *
C             ORIGINAL CODE WRITTEN                02.09.81            *
C             MODIFIED                             31-JAN-86           *
C             DENDRO CALL ADDED                    NOV 86 ARDEAN LEITH *
C             LONG FILE NAMES                      JAN 89 ARDEAN LEITH *
C             INCLUDED FILES FOR SCLASSY, SEMIS    MAR 02 ARDEAN LEITH *
C             NEW IMC FORMAT                       OCT 02 ARDEAN LEITH *
C             EXCESSIVE PARTITION TRAP             DEC 05 ARDEAN LEITH *
C             NFAC VS KFAC BUG                     DEC 07 ARDEAN LEITH *
C             IPALIGN & REFACTORING                DEC 08 ARDEAN LEITH *
C             KV FACTOR SELECTION                  NOV 11 ARDEAN LEITH *
C                                                                      *
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
C AUTHOR:  JEAN-PIERRE BRETAUDIERE (DECEASED)                          *
C      THE UNIVERSITY OF TEXAS HEALTH SCIENCE CENTER AT HOUSTON        *
C      MEDICAL SCHOOL - DEPT. OF PATHOLOGY AND LABORATORY MEDICINE     *
C      P.O. BOX 20708, HOUSTON, TX 77225.                              *
C                                                                      *
C **********************************************************************
C                                                                      *         
C  SCLASSI(LUNI,LUNK,LUNDOC)                                           *
C                                                                      *         
C  PURPOSE:  PERFORMS FIRST STEP AN AUTOMATIC CLUSTERING OF            *         
C            OBJECTS BY AGGREGATION AROUND MOBILE CENTERS AND IN A 2ND *
C            STEP, A HIERARCHIC ASCENDENT CLASSIFICATION OF THE GRAVITY*         
C            CENTERS OF THE CLUSTERS DETERMINED IN THE FIRST STEP.     *
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE SCLASSI(LUNI,LUNK,LUNDOC)

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER(LEN=MAXNAM) :: CLUSFILE,IMCFILE,FILPRE

#ifndef SP_32
        INTEGER *8            :: IBIG
#else
        INTEGER *4            :: IBIG
#endif

        CHARACTER(LEN=1)      :: NULL = CHAR(0)

        REAL, ALLOCATABLE     :: Q(:)
                                                                                                          
C       MARCH 02 al NKLA IS REDEFINED LATER! SO REPEAT FAILED
        NKLA   = 100

        CALL FILERD(FILPRE,NLET,NULL,
     &          'CORAN/PCA FILE PREFIX (e.g.. CORAN_01_)~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IMCFILE = FILPRE(1:NLET) // '_IMC'//NULL

C       GET CLUSTER OUTPUT FILE NAME
        CALL FILERD(CLUSFILE,NLET,DATEXC,'CLUSTER OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       OPEN & READ HEADER OF _IMC FILE FORM='FORMATTED'
        CALL OPAUXFILE(.FALSE.,IMCFILE,DATEXC,LUNI,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        READ(LUNI,*) NUMIM, NFAC,IDUM,IDUM,IDUM,IDUM                                       

        WRITE(NOUT,90) NFAC,NUMIM
90      FORMAT(/,'  FACTORS AVAILABLE:',I5,'  OBJECTS (IMAGES):',I6)

        KFAC   = NFAC
        MINFAC = 1
        MAXFAC = NFAC

        WRITE(NOUT, *) ' WARNING. INPUT ALTERED 2011.  LIST',
     &                  ' ALL FACTORS WANTED, E.G. 1-7, NOT: 7'
                                      
        CALL RDPRAI(INUMBR,NIMAX,KFAC, MINFAC,MAXFAC,
     &              'FACTOR NUMBERS TO BE USED',NULL,IER)
        IF (IRTFLG .NE. 0) RETURN
        IF (KFAC <= 0) THEN
           CALL ERRT(102,'NO FACTORS REQUESTED',KFAC)                              
           RETURN
        ENDIF

C       SORT THE INUMBER FACTOR LIST IN ASCENDING ORDER
        IF (KFAC > 1) CALL SORTI(INUMBER,KFAC)

        NITER = 5
        NCLAS = 5
        CALL RDPRIS(NITER,NCLAS,NOT_USED,
     &       'NUMBER OF ITERATIONS & CENTERS',IRTFLG)
        IF (NITER .LE. 0) THEN
           CALL ERRT(102,'ILLEGAL: NUMBER OF ITERATIONS',NITER)                              
           RETURN
        ELSEIF (NCLAS .LE. 0) THEN
           CALL ERRT(102,'ILLEGAL: NUMBER OF CENTERS',NCLAS)                              
           RETURN
        ENDIF

        NBASE = 4
        CALL RDPRI1S(NBASE,NOT_USED,'NUMBER OF PARTITIONS',IRTFLG)
        IF (NBASE .LE. 0) THEN
           CALL ERRT(102,'ILLEGAL: NUMBER OF PARTITIONS',NBASE)                              
           RETURN
        ENDIF
                                                        
C       OPEN NEW CLUSTER FILE (ALREADY HAS EXTENSION ON IT)
        CALL OPAUXFILE(.FALSE.,CLUSFILE,NULL,-LUNK,0,
     &                       'N', ' ',.TRUE.,IRTFLG)

        WRITE(NDAT, 2100) (INUMBR(K), K=1,KFAC)                                      
 2100   FORMAT (/,'  FACTORS USED:',(10I6))
                                         
        WRITE(NDAT, 2150) NBASE, NITER, NCLAS, NKLA                               
 2150   FORMAT('  PARTITIONS:',I6,'  ITERATIONS:      ',I6,/,                          
     &         '  CENTERS:   ',I6,'  INITIAL CLASSES: ',I6,/)
                                              
C       MEMORY SEGMENTATION                                                  
        N2DIM = MAX(NKLA, 2*NCLAS)                                              
        KDIM  = MAX(NUMIM, NCLAS**NBASE)                                        
        KDIM  = MAX(KDIM,  2*NKLA - 1)                                          
        MDIM  = NKLA*(NKLA-1) / 2                                                 
        LDIM  = NKLA                                                              
        L2DIM = 2*NKLA - 1                                                        

        IBIG  = NCLAS
        IBIG  = (IBIG**NBASE) * 3     
        IBIG  = IBIG + (LDIM * 7) + (L2DIM * 4) + 2 * NUMIM +
     &          MDIM  + (3 *KDIM)+ NKLA + NFAC + N2DIM * KFAC    
        IBIG4 = HUGE(IBIG4)
        IF (IBIG .GE. IBIG4) THEN
           WRITE(NOUT,*)' *** MUST REDUCE NUMBER OF PARTITIONS'
           CALL ERRT(102,'EXCESSIVE MEMORY NEEDED>',IBIG4)
           GOTO 9999
        ENDIF

        N1    = 1
        NK    = IPALIGN64(N1   + NUMIM)                                                    
        ND    = IPALIGN64(NK   + NUMIM)                                                    
        NU    = IPALIGN64(ND   + MDIM)                                                        
        NJV   = IPALIGN64(NU   + KDIM)                                                        
        NJW   = IPALIGN64(NJV  + KDIM)                                                       
        NIDK  = IPALIGN64(NJW  + KDIM)                                                       
        NCI   = IPALIGN64(NIDK + NKLA)                                                      
        NGT   = IPALIGN64(NCI  + NFAC)                                                        
        NNUM  = IPALIGN64(NGT  + N2DIM * KFAC)                                                 
        NLA   = IPALIGN64(NNUM + LDIM)                                                      
        NLB   = IPALIGN64(NLA  + LDIM)                                                       
        NIV   = IPALIGN64(NLB  + LDIM)                                                        
        NIW   = IPALIGN64(NIV  + LDIM)                                                        
        NV    = IPALIGN64(NIW  + LDIM)                                                        
        NW    = IPALIGN64(NV   + LDIM)                                                         
        NNT   = IPALIGN64(NW   + LDIM)                                                        
        NVAL  = IPALIGN64(NNT  + L2DIM)                                                       
        NPK   = IPALIGN64(NVAL + L2DIM)                                                      
        NNO   = IPALIGN64(NPK  + L2DIM)                                                       
        NFIN  = IPALIGN64(NNO  + L2DIM)
                                                       
        ALLOCATE (Q(NFIN),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'SCLASSI; Q',NFIN)
           GOTO 9999
        ENDIF

C       CLASSIFICATION OF OBJECTS ACCORDING TO FACTORIAL COORDINATES          
        CALL NOYAU(N2DIM,KDIM,MDIM,LDIM,L2DIM,                                   
     &   NFAC,KFAC,NUMIM,INUMBR,NBASE,NITER,NCLAS,NKLA,
     &   Q(N1), Q(NK),  Q(ND),  Q(NU), Q(NJV),Q(NJW),Q(NIDK),
     &   Q(NCI),Q(NGT), Q(NNUM),Q(NLA),Q(NLB),Q(NIV),Q(NIW),Q(NV),Q(NW),                        
     &   Q(NNT),Q(NVAL),Q(NPK), Q(NNO), LUNI,LUNK,LUNDOC)
                                               
C       PRINT OUT LIST OF CLUSTER MEMBERS, LIST OF CENTERS OF GRAVITY, 
C       AND RE-CLASSIFICATION LOOKUP TABLE

        REWIND(LUNK)                                                                        
        READ(LUNK) NUMIM, NFAC, NKLA, KFAC  
                
        NK    = 1                                                                 
        NI    = IPALIGN64(NK  + NUMIM )                                                       
        NPK   = IPALIGN64(NI  + NUMIM)                                                        
        NGT   = IPALIGN64(NPK + NKLA )                                                     
        NIV   = IPALIGN64(NGT + NKLA * KFAC)                                               
                                                     
C       READS: LUNI AND LUNK (WRITTEN FROM NOYAU)                
        CALL RGRI(NUMIM, KFAC, NKLA,                                              
     &            Q(NK),Q(NI),Q(NPK),Q(NGT),Q(NIV), 
     &            LUNK,LUNI,NFAC,INUMBR)

9999    CLOSE(LUNK)
        CLOSE(LUNI)

        END












