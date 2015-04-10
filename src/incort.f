
C***********************************************************************
C INCORT
C                         REWRITTEN      SEP 2003 ARDEAN LEITH
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
C  INCORT(NPIX,NUMIM,NFAC,LUNT,USE_PCA,
C         MATS,EVALS,U,WEIGHTI,WEIGHTP,TRACE,SUMW,IRTFLG)
C
C  PARAMETERS:
C       NPIX	  NUMBER OF PIXELS                              (INPUT)
C	NUMIM     NUMBER OF IMAGES                              (INPUT)
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED              (INPUT)
C	LUNT	  TRANSPOSED IMAGE I/O UNIT                     (INPUT)
C	USE_PCA   CORAN VS PCA FLAG                             (INPUT)
C	MATS()    MATRIX & EIGENVECTOR ARRAY                    (OUTPUT) 
C	EVALS()   EIGENVALUE ARRAY                              (OUTPUT)
C	BLU()     INPUT BUFFER                                  (WORK)
C	WEIGHTI() SUM OF PIXEL VALUES FOR THIS IMAGE            (INPUT)
C	WEIGHTP() SUM OF PIXEL VALUES AT THIS PIXEL             (INPUT)
C	TRACE	  SUM OF THE ELEMENTS ARRAY DIAGONAL            (OUTPUT)
C	SUMW      SUM OF ALL THE PIXEL VALUES IN ALL IMAGES     (INPUT)
C	IRTFLG    ERROR FLAG                                    (OUTPUT)
C
C  NOTE:  IMAGE ROW & COLUMNS HAVE BEEN TRANSPOSED ALREADY ON UNIT LUNT
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
        
        SUBROUTINE INCORT(NPIX,NUMIM,NFAC,LUNT,USE_PCA,MATS,
     &                    EVALS,BLU,WEIGHTI,WEIGHTP,TRACE,SUMW,IRTFLG)
                                                  
        INCLUDE 'CMBLOCK.INC'

        REAL                              :: MATS(NUMIM,NUMIM)  
        REAL                              :: EVALS(NUMIM), BLU(NUMIM)
        REAL                              :: WEIGHTI(NUMIM) 
        REAL                              :: WEIGHTP(NPIX)
        REAL, ALLOCATABLE, DIMENSION(:,:) :: EVECTS

C       AUTOMATIC ARRAYS
        INTEGER,PARAMETER                 :: LWORK=100
        REAL                              :: WORK(LWORK*NPIX)
        INTEGER                           :: IWORK(5*NUMIM)
        INTEGER                           :: IFAIL(NUMIM)

        LOGICAL                           :: USE_PCA

C       SET MATS ARRAY = 0.0
        MATS     = 0.0

C       READ ALL THE ROWS. EACH ROW CONTAINS VALUES FROM ONE PIXEL
        DO I = 1, NPIX

C           READ THE WHOLE PIXEL  IN BLU ARRAY.
            READ(LUNT,REC=I,IOSTAT=IERR) BLU

C           WE ARE ASSUMING ARRAY IS SYMMETRICAL.
C	    COMPUTE  MATS(NUMIM,NUMIM) = MATS'. MATS WHERE MATS' IS MATS TRANSPOSE.

            IF (USE_PCA) THEN

C              SUBSTRACT AV. PIXEL VALUE.
               BLU = BLU - (WEIGHTP(I) / NPIX)

               DO J=1,NUMIM
                  DO JJ =1,J
                     MATS(J, JJ) = MATS(J,JJ) + (BLU(J) * BLU(JJ))   
                     MATS(JJ,J)  = MATS(J,JJ)   
                  ENDDO
               ENDDO
            ELSE
               PIA = WEIGHTP(I)

               DO J=1,NUMIM
                  DO JJ =1,J
                     MATS(J,JJ) = MATS(J,JJ) + (BLU(J) * BLU(JJ)) / PIA 
                  ENDDO
               ENDDO
            ENDIF
        ENDDO


C       ALL THE IMAGES HAVE BEEN READ.

        IF (.NOT. USE_PCA) THEN
C          WEIGHTI(K) IS THE SUM OF ALL PIXELS IN ONE IMAGE
           DO J =1,NUMIM
              DO JJ =1,J
                 AAA        =  SQRT(WEIGHTI(J) * WEIGHTI(JJ))
                 MATS(J,JJ) =  MATS(J, JJ) / AAA  -  AAA / SUMW
              ENDDO
           ENDDO
        ENDIF

C       "TRACE" OF A MATRIX IS THE SUM OF THE ELEMENTS ON DIAGONAL
        TRACE = 0.0
        DO J=1, NUMIM
           TRACE = TRACE + MATS(J, J)
        ENDDO

C       COMPUTE EIGENVALUES AND EIGENVECTORS.

#ifdef  SP_LAPACK
C       USE LAPACK EIGENVECTOR ROUTINE

        ALLOCATE (EVECTS(NUMIM,NUMIM),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = NUMIM**2
            CALL ERRT(46,'EVECTS',MWANT)
            RETURN
        ENDIF

        LWORKT = LWORK * NPIX 
        WRITE(NOUT,*) ' Using LAPACK eigenvalues routine.'
        CALL ssyevx('V', 'A', 'L', NUMIM, MATS, NUMIM, 
     &              DUM,DUM, IDUM,IDUM,
     &              0.0, NGOT, EVALS, EVECTS, NUMIM, 
     &              WORK, LWORKT, IWORK, IFAIL, IRTFLG)

         WRITE(NOUT,*)' RETURNED:',NGOT, 'EIGENVALUES.'

         WRITE(NOUT,*)
     &        ' USED LWORK: ',LWORKT,'  OPTIMAL LWORK: ', WORK(1)

C        REVERSE ORDER OF EIGENVALUES & EIGENVECTORS 
         IEND = NGOT / 2

         MATS = EVECTS
         IF (ALLOCATED(EVECTS)) DEALLOCATE(EVECTS)

         DO I=1,IEND
            IT         = NGOT-I+1
            TMP        = EVALS(I)
            EVALS(I)   = EVALS(IT)
            EVALS(IT)  = TMP

            BLU        = MATS(:,I)
            MATS(:,I)  = MATS(:,IT)
            MATS(:,IT) = BLU(:)
         ENDDO 

#else
C       USE SPIDER EIGENVECTOR ROUTINE
        CALL VPROP(NUMIM,  NUMIM,  MATS, EVALS,  BLU,  IRTFLG)
#endif       
        IF (IRTFLG .NE. 0)  THEN
           CALL ERRT(102,'DIAGONALIZATION FAILURE',IRTFLG)
           RETURN
        ENDIF

C       EVALS() HOLDS THE EIGENVALUES 
C       WHILE MATS() HOLDS THE CORRESPONDING EIGENVECTORS.

c       WRITE(NOUT,*) ' EIGENVALUES:'
c       WRITE(NOUT,90) (EVALS(K),K=1,NFAC)
c90     FORMAT(5(1PG12.5,' '))
c       WRITE(NOUT,*) ' EIGENVECTORS:'
c       DO I=1,NFAC
c            WRITE(NOUT,92) (MATS(K,I),K=1,NUMIM)
c       ENDDO
c92     FORMAT(5(1PG12.5,' '),/,5(1PG12.5,' '),/,
c     &        5(1PG12.5,' '),/,5(1PG12.5,' '),/,1(1PG12.5,' ')/)

     
        IF (.NOT. USE_PCA) THEN
           DO L=1,NFAC
              DO K=1,NUMIM
                MATS(K,L) = MATS(K,L) * SQRT(SUMW / WEIGHTI(K))
              ENDDO
           ENDDO
        ENDIF
  
        RETURN
        END

#ifdef NEVER
original old SPIDER
 1.63677E-06  1.78826E-06  1.87376E-06  1.93852E-06  2.14307E-06  
 2.21764E-06  2.31187E-06  2.35837E-06  2.63528E-06  4.35528E-06
#endif
