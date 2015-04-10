
C***********************************************************************
C INCOR3
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
C  INCOR3(NUMIM,NPIX,NFAC,USE_PCA,
C         MATS,EVALS,WEIGHTI,WEIGHTP,TRACE,SUMW,IRTFLG)
C
C  PARAMETERS:
C       NUMIM	  NUMBER OF IMAGES                              (INPUT)
C	NPIX      NUMBER OF PIXELS PER IMAGE                    (INPUT)
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED              (INPUT)
C	USE_PCA   CORAN VS PCA FLAG                             (INPUT)
C	MATS	  EIGENVECTOR ARRAY                             (OUTPUT)                                    (OUTPUT)
C	EVALS	  EIGENVALUE ARRAY                              (OUTPUT)
C	U	  INPUT BUFFER                                  (WORK)
C	WEIGHTI	  SUM OF PIXEL VALUES FOR THIS IMAGE             (INPUT)
C	WEIGHTP	  SUM OF PIXEL VALUES AT THIS PIXEL              (INPUT)
C	TRACE	  SUM OF THE ELEMENTS ARRAY DIAGONAL            (OUTPUT)
C	SUMW      SUM OF ALL THE PIXEL VALUES IN ALL IMAGES     (INPUT)
C	IRTFLG    ERROR FLAG                                    (OUTPUT)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
        
        SUBROUTINE INCOR3(NUMIM,NPIX,NFAC,LUNS,USE_PCA,MATS,
     &                    EVALS,BLU,WEIGHTI,WEIGHTP,TRACE,SUMW,IRTFLG)
                                                  
        INCLUDE 'CMBLOCK.INC'

        REAL                              :: MATS(NPIX,NPIX)  
        REAL                              :: EVALS(NPIX)
        REAL                              :: WEIGHTP(NPIX)
        REAL                              :: WEIGHTI(NUMIM)
        REAL                              :: BLU(NPIX)
        REAL, ALLOCATABLE, DIMENSION(:,:) :: EVECTS

C       AUTOMATIC ARRAYS
        INTEGER,PARAMETER                 :: LWORK=100
        REAL                              :: WORK(LWORK*NPIX)
        INTEGER                           :: IWORK(5*NPIX)
        INTEGER                           :: IFAIL(NPIX)

        LOGICAL                           :: USE_PCA

C       SET MATS ARRAY = 0.0
        MATS     = 0.0

C       POSITION _SEQ FILE
        REWIND(LUNS)
        READ(LUNS,IOSTAT=IERR) NUMIMT,NPIXT

C       READ ALL THE ROWS. EACH ROW CONTAINS VALUES FROM ONE IMAGE

C       WE ARE ASSUMING ARRAY IS SYMMETRICAL.
C       COMPUTE  MATS(NUMIM,NUMIM) = MATS' . MATS WHERE MATS' IS MATS TRANSPOSE.

        DO I = 1, NUMIM
C          READ THE WHOLE IMAGE IN BLU ARRAY.
           READ(LUNS,IOSTAT=IERR) BLU,FIM

           IF (USE_PCA) THEN
              DO J=1,NPIX
                 DO JJ =1,J
                    MATS(J, JJ) = MATS(J,JJ) + (BLU(J) * BLU(JJ))   
                    MATS(JJ,J)  = MATS(J,JJ)   
                 ENDDO
              ENDDO
           ELSE
              PIA = WEIGHTI(I)

              DO J=1,NPIX
                 DO JJ =1,J
                    MATS(J,JJ) = MATS(J,JJ) + (BLU(J) * BLU(JJ)) /PIA 
                 ENDDO
              ENDDO
           ENDIF
        ENDDO

C       ALL THE IMAGES HAVE BEEN READ INTO MATS.

C       WEIGHTP(K) IS THE SUM OF THIS PIXEL IN ALL IMAGES.
        IF (USE_PCA) THEN
C          PCA
           DO J =  1, NPIX
             DO JJ  =  1, J
               MATS(J, JJ) =  (MATS(J, JJ) - 
     &                      (WEIGHTP(J) * WEIGHTP(JJ)) / NUMIM) 
               MATS(JJ, J) = MATS(J, JJ)                                 
             ENDDO
           ENDDO

        ELSE
C          CORAN
           DO J =1,NPIX
              DO JJ =1,J
                 AAA        =  SQRT(WEIGHTP(J) * WEIGHTP(JJ))
                 MATS(J,JJ) =  MATS(J, JJ) / AAA  -  AAA / SUMW
              ENDDO
           ENDDO
        ENDIF

C       THE "TRACE" OF A MATRIX IS THE SUM OF THE ELEMENTS ON DIAGONAL
        TRACE = 0.0
        DO J = 1, NPIX
           TRACE = TRACE + MATS(J, J)
        ENDDO

C       COMPUTE EIGENVALUES AND EIGENVECTORS.

#ifdef  SP_LAPACK

        ALLOCATE (EVECTS(NPIX,NPIX),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            MWANT = NPIX**2
            CALL ERRT(46,'EVECTS',MWANT)
            RETURN
        ENDIF

        LWORKT = LWORK * NPIX 
        WRITE(NOUT,*) ' Using LAPACK eigenvalues routine.'
        CALL ssyevx('V', 'A', 'L', NPIX, MATS, NPIX, 
     &              DUM,DUM, IDUM,IDUM,
     &              0.0, NGOT, EVALS, EVECTS, NPIX, 
     &              WORK, LWORKT, IWORK, IFAIL, IRTFLG)

         WRITE(NOUT,*)' RETURNED:',NGOT, 'EIGENVALUES.'

         WRITE(NOUT,*)
     &         ' USED LWORK: ',LWORKT,'  OPTIMAL LWORK: ', WORK(1)

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
         CALL VPROP(NPIX,  NPIX,  MATS, EVALS,  BLU,  IRTFLG)
#endif       
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(101,'DIAGONALIZATION FAILURE',IRTFLG)
           RETURN
        ENDIF

C       EVALS() HOLDS THE EIGENVALUES 
C       WHILE MATS() HOLDS THE CORRESPONDING EIGENVECTORS.

c       WRITE(NOUT,*) ' EIGENVALUES:'
c       WRITE(NOUT,90) (EVALS(K),K=1,NFAC)
c90     FORMAT(5(1PG12.5,' '))
c       WRITE(NOUT,*) ' EIGENVECTORS:'
c       DO I=1,NFAC
c            WRITE(NOUT,92) (MATS(K,I),K=1,NPIX)
c       ENDDO
c92     FORMAT(5(1PG12.5,' '),/,5(1PG12.5,' '),/,
c     &        5(1PG12.5,' '),/,5(1PG12.5,' '),/,1(1PG12.5,' ')/)

     
        IF (.NOT. USE_PCA) THEN
           DO L=1,NFAC
              DO K=1,NPIX
                MATS(K,L) = MATS(K,L) * SQRT(SUMW / WEIGHTP(K))
              ENDDO
           ENDDO
        ENDIF
  
        RETURN
        END

