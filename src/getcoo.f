
C **********************************************************************
C  GETCOO 
C                         REWRITTEN         SEP 2003 ARDEAN LEITH
C                         AVERAGE WEIGHTP   MAY 2004 ARDEAN LEITH
C
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
C  GETCOOT(NFAC,NPIX,NUMIM,INUMBR,USE_PCA,
C          EVECTS,EVALS,WEIGHTI,WEIGHTP,SUMW,CO, BLU, BLW, 
C          LUNS, LUNI, LUNP)
C
C  PURPOSE: GET COORDINATES FOR NON-TRANSPOSED DATA
C
C  PARAMETERS:           
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED              (INPUT)
C	NPIX      NUMBER OF PIXELS PER IMAGE                    (INPUT)
C       NUMIM	  NUMBER OF IMAGES                              (INPUT)
C       INUMBR()  IMAGE NUMBER LIST                             (INPUT)
C	USE_PCA   CORAN VS PCA FLAG                             (INPUT)
C       EVECTS()  EIGENVECTORS (COLUMN)  OF X'X AND             (INPUT)
C                 X(I,*)= BLU() J=1,JTOT W/ I=1,ITOT
C       EVALS()   EIGENVALUE ARRAY                           (INPUT/OUTPUT)
C	WEIGHTI() SUM OF PIXEL VALUES FOR THIS IMAGE            (INPUT)
C	WEIGHTP() SUM OF PIXEL VALUES AT THIS PIXEL             (INPUT)
C       CO()      WORKING ARRAY  
C	SUMW      SUM OF ALL THE PIXEL VALUES IN ALL IMAGES     (INPUT)
C       BLU()     WORKING ARRAY FOR INPUTS
C       BLW()     WORKING ARRAY FOR OUTPUTS
C	LUNS	  SEQUENTIAL IMAGE I/O UNIT ( FOR INPUT)        (INPUT)
C	LUNI      IMAGE COORDINATE I/O UNIT (FOR OUTPUT)        (INPUT)
C	LUNP      PIXEL COORDINATE I/O UNIT (FOR OUTPUT)        (INPUT)
C
C	S(,) HAS THE EIGENVECTORS (COLUMN)  OF X'X AND D() HAS THE
C       EIGENVALUES. X(I,*)= U() J=1,JTOT W/ I=1,ITOT
C
C  NOTE: FOR CORAN WEIGHTP IS RETURNED UNCHANGED
C        FOR CORAN WEIGHTI IS RETURNED DIVIDED BY SUMW
C        FOR CORAN _PIX    IS WEIGHTP  UNCHANGED
C        FOR CORAN _IMC    IS WEIGHTI  DIVIDED BY SUMW
C
C        FOR PCA   WEIGHTP IS RETURNED UNCHANGED
C        FOR PCA   WEIGHTI IS RETURNED UNCHANGED
C        FOR PCA   _PIX    IS WEIGHTP  UNCHANGED
C        FOR PCA   _IMC    IS WEIGHTI  UNCHANGED
C
C **********************************************************************
      
        SUBROUTINE GETCOO(NFAC,NPIX,NUMIM,INUMBR,USE_PCA,
     &             EVECTS,EVALS,WEIGHTI,WEIGHTP,SUMW,CO, BLU, BLW, 
     &             LUNS, LUNI, LUNP)

        REAL    :: EVECTS(NPIX,NPIX), EVALS(NPIX)
        REAL    :: CO(NPIX), BLU(NPIX), BLW(NFAC)
        INTEGER :: INUMBR(NUMIM)
        REAL    :: WEIGHTP(NPIX)
        REAL    :: WEIGHTI(NUMIM)
        LOGICAL :: USE_PCA

C       POSITION _SEQ, _IMC, & _PIX FILES
        CALL REW(LUNS,  1)

        REWIND(LUNI)
        READ(LUNI,*)  IDUM, IDUM, IDUM, IDUM, IDUM, IDUM

        REWIND(LUNP)
        READ(LUNP,*) IDUM, IDUM, IDUM, IDUM, IDUM
        
C       INITIALIZE ARRAYS  
        CO     =  0.0

        IF (USE_PCA) THEN
           DOR = 0.0
           PI  = 0.0       
        ELSE
C          WHOLE ARRAY
           WEIGHTP =  WEIGHTP / SUMW 
        ENDIF      
        ACT = 1.0

C       WRITE _IMC  DATA

        DO I  =  1, NUMIM
C         READ THE WHOLE IMAGE  IN BLU ARRAY FROM _SEQ FILE.
          READ(LUNS,IOSTAT=IERR) BLU,FIM

C         COORDINATES  TO ORIGIN FOR THE ROWS
C         [W] = [X -<X>].[EIGENVECTORS]   MATIX NOTATION. 
C         <XI> AVERAGE VALUE OF PIXELS = SUM (XIJ), J=1,NUMIM
 
          PIA = WEIGHTI(I)

          DO  K=1,NFAC                                                   
             BLW(K) = 0.0 

             IF (USE_PCA) THEN
                DO J=1,NPIX
                   BLW(K) = BLW(K) + (BLU(J) - WEIGHTP(J) / NUMIM) * 
     &                      EVECTS(J, K)
                ENDDO
             ELSE
                DO J=1,NPIX
                   BLW(K) = BLW(K) + (BLU(J) * EVECTS(J, K)) / PIA
                ENDDO
             ENDIF
          ENDDO

C         DISTANCES TO ORIGIN FOR THE ROWS                    
          DOR  =  0.0 
          FIM =  INUMBR(I)
          IF (USE_PCA) THEN
             DO J=1, NPIX
                DOR =  DOR + (BLU(J) - WEIGHTP(J))**2
             ENDDO

             WRITE(LUNI,90) (BLW(K),K=1,NFAC), PIA, DOR, FIM, ACT
          ELSE
             DO J=1,NPIX
                DOR = DOR + (BLU(J)/PIA - WEIGHTP(J))**2 / WEIGHTP(J)
             ENDDO

             PI = PIA / SUMW   
             DO J=1,NPIX
                CO(J) = CO(J) + (BLU(J)/(SUMW * WEIGHTP(J)) - PI)**2/PI
             ENDDO

             WRITE(LUNI,90) (BLW(K), K = 1, NFAC), PI, DOR, FIM, ACT
          ENDIF
        ENDDO

C       WRITE _PIX  DATA

C       COORDINATES  FOR THE COLS

        FDUM = 0.0   ! FOR COMPATIBILITY WITH IMC
        DO J=1,NPIX

          DO K=1,NFAC
C            IF IMAGES HAVE SOME CROSS-CORRELATION (I.E; (PART OF) ONE 
C	     IMAGE IS THE SAME AS (PART OF) ANOTHER), ONE FACTOR
C	     WILL BE 0 (SO A NEGATIVE EIGENVALUE MAY EXIST).

             EVALS(K) = MAX(EVALS(K), 1.0E - 9)
             BLW(K)   = EVECTS(J, K) * SQRT(EVALS(K))            
          ENDDO

C         WRITE  DATA TO PIX*** FILE 
          FPIX = J
          PIAT = WEIGHTP(J) / NPIX
          WRITE(LUNP,90) (BLW(K),K=1,NFAC),PIAT, CO(J), FPIX, FDUM 
90        FORMAT(10(1PG12.5,' '))
        ENDDO

        RETURN
        END
