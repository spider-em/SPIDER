
C **********************************************************************
C
C  GETCOOT
C                         REWRITTEN            SEP 2003 ARDEAN LEITH
C                         AVERAGE WEIGHTP      MAY 2004 ARDEAN LEITH
C                         PIXEL FILE FDUM      JUN 2009 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C
C  GETCOOT(NFAC,NUMIM,NPIX,INUMBR,USE_PCA,
C          EVECTS,EVALS,WEIGHTI,WEIGHTP,SUMW, CO, BLU, BLW, 
C          LUNS, LUNI, LUNP)
C
C  PURPOSE: GET COORDINATES FOR TRANSPOSED DATA
C
C  PARAMETERS:           
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED              (INPUT)
C	NPIX      NUMBER OF PIXELS PER IMAGE                    (INPUT)
C       NUMIM	  NUMBER OF IMAGES                              (INPUT)
C       INUMBR()  IMAGE NUMBER LIST                             (INPUT)
C	USE_PCA   CORAN VS PCA FLAG                             (INPUT)
C       EVECTS    EIGENVECTORS (COLUMN)  OF X'X AND             (INPUT)
C                 X(I,*)= BLU() J=1,JTOT W/ I=1,ITOT
C       EVALS()    EIGENVALUES                                  (INPUT)
C	WEIGHTI() SUM OF PIXEL VALUES FOR THIS IMAGE            (INPUT)
C	WEIGHTP() SUM OF PIXEL VALUES AT THIS PIXEL             (INPUT)
C	SUMW      SUM OF ALL THE PIXEL VALUES IN ALL IMAGES     (INPUT)
C       CO()      WORKING ARRAY  
C       BLU()     WORKING ARRAY FOR INPUTS
C       BLW()     WORKING ARRAY FOR OUTPUTS
C	LUNS	  SEQUENTIAL IMAGE I/O UNIT ( FOR INPUT)        (INPUT)
C	LUNI      IMAGE COORDINATE I/O UNIT (FOR OUTPUT)        (INPUT)
C	LUNP      PIXEL COORDINATE I/O UNIT (FOR OUTPUT)        (INPUT)
C
C  NOTE: FOR CORAN WEIGHTP IS RETURNED UNCHANGED
C        FOR CORAN WEIGHTI IS RETURNED DIVIDED BY SUMW
C        FOR CORAN _PIX    IS WEIGHTP  DIVIDED BY SUMW
C        FOR CORAN _IMC    IS WEIGHTI  DIVIDED BY SUMW
C
C        FOR PCA   WEIGHTP IS RETURNED UNCHANGED
C        FOR PCA   WEIGHTI IS RETURNED UNCHANGED
C        FOR PCA   _PIX    IS WEIGHTP  UNCHANGED
C        FOR PCA   _IMC    IS WEIGHTI  UNCHANGED
C
C **********************************************************************

        SUBROUTINE GETCOOT(NFAC,NUMIM,NPIX,INUMBR,USE_PCA,
     &                  EVECTS,EVALS,WEIGHTI,WEIGHTP,SUMW, CO, BLU, BLW, 
     &                  LUNS, LUNI, LUNP)
            
        REAL    :: EVECTS(NUMIM, NUMIM), EVALS(NUMIM), CO(NUMIM)
        REAL    :: BLU(NUMIM), BLW(NFAC)
        INTEGER :: INUMBR(NUMIM)
        REAL    :: WEIGHTI(NUMIM)
        REAL    :: WEIGHTP(NPIX)
        LOGICAL :: USE_PCA

C       SKIP HEADERS ON RELEVANT I/O FILES          
        REWIND(LUNI)
        READ(LUNI,*)  IDUM, IDUM, IDUM, IDUM, IDUM, IDUM

        REWIND(LUNP)
        READ(LUNP,*) IDUM, IDUM, IDUM, IDUM, IDUM
      
        IF (USE_PCA) THEN
           CODUM = 0.0 
           DOR   = 0.0    
        ELSE
C          INITIALIZE ARRAYS  
           CO      =  0.0
           WEIGHTI =  WEIGHTI / SUMW 
        ENDIF      
        ACT  = 1
        FDUM = 0.0

C       WRITE PIX  DATA
        DO I = 1,NPIX

C         READ THE WHOLE PIXEL  IN BLU ARRAY FROM _SEQ FILE.
          READ(LUNS,REC=I,IOSTAT=IERR) BLU

          FPIX = I
          PIA  = WEIGHTP(I) 

C         COORDINATES TO ORIGIN FOR THE ROWS (SINCE TRANSPOSED)
          DO K=1,NFAC                                                   
            BLW(K) = 0.0 

            IF (USE_PCA) THEN
                DO J=1, NUMIM
                   BLW(K) = BLW(K) + (BLU(J) * EVECTS(J,K))
                ENDDO
             ELSE
                DO J=1, NUMIM
                   BLW(K) = BLW(K) + (BLU(J) * EVECTS(J,K)) / PIA
                ENDDO
             ENDIF
          ENDDO


          IF (USE_PCA) THEN
C            FOR PCA ANALYSIS
     
C            WRITE  DATA TO PIX FILE 
             FPIX = I

C            GET AVERAGE FOR _PIX FILE (MAY 04)
             PIAT = PIA / NPIX
             WRITE(LUNP,90) (BLW(K),K=1,NFAC), PIAT, DOR, FPIX, FDUM 
90           FORMAT(10(1PG12.5,' '))

          ELSE
C            FOR CORAN FIND DISTANCES TO ORIGIN FOR THE ROWS     
             DOR = 0.0   
             DO J=1,NUMIM
                DOR = DOR + (BLU(J)/PIA - WEIGHTI(J))**2 / WEIGHTI(J)
             ENDDO

             PI = PIA / SUMW   
             DO J=1,NUMIM
                CO(J) = CO(J) + (BLU(J)/(SUMW * WEIGHTI(J)) - PI)**2/PI
             ENDDO

C            WRITE  DATA TO PIX FILE 
C            GET AVERAGE FOR _PIX FILE (MAY 04)
             PIT = PI * SUMW / NPIX
             WRITE(LUNP,90) (BLW(K),K=1,NFAC), PIT, DOR, FPIX, FDUM
          ENDIF
        ENDDO

C       WRITE IMC DATA
C       COORDINATES  FOR THE COLS (SINCE TRANSPOSED)
        DO J=1,NUMIM
          DO K=1,NFAC 
C            IF  IMAGES HAVE SOME CROSS-CORRELATION (I.E; (PART OF) ONE 
C	     IMAGE IS THE SAME AS (PART OF) ANOTHER), ONE OF THE FACTORS
C	     WILL BE 0 (SO A NEGATIVE EIGENVALUE MAY EXIST).
             EVALS(K) = MAX(EVALS(K), 1.0E-9)

             BLW(K) = EVECTS(J,K) * SQRT(EVALS(K))            
          ENDDO

C         WRITE IMC DATA TO IMC FILE 
          FIM = INUMBR(J)
          IF (USE_PCA) THEN
             WRITE(LUNI,90) (BLW(K),K=1,NFAC),WEIGHTI(J),CODUM,FIM,ACT
          ELSE
             WRITE(LUNI,90) (BLW(K),K=1,NFAC),WEIGHTI(J),CO(J),FIM,ACT
          ENDIF
        ENDDO

        RETURN
        END
