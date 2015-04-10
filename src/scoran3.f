
C **********************************************************************
C
C   SCORAN3               USED ALLOCATE       JAN 2001 ARDEAN LEITH
C                         ENLARGED MADAR      MAR 2002 ARDEAN LEITH
C                         TSIZE BUG           AUG 2002 ARDEAN LEITH
C                         REWRITTEN           SEP 2003 ARDEAN LEITH
C                         LUNE HEADER SPACE   MAR 2009 ARDEAN LEITH
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
C *********************************************************************
C
C   SCORAN3(NUMIM, NFAC, NPIX, INUMBR,LUNS, LUNI, LUNP, LUNE, 
C           WEIGHTI,WEIGHTP, TRANSPOSE)
C
C       NUMIM	  NUMBER OF IMAGES                               (INPUT)
C	NFAC 	  NUMBER OF EIGENVECTORS REQUESTED               (INPUT)
C	NPIX      NUMBER OF ACTIVE PIXELS PER IMAGE              (INPUT)
C       INUMBR	  IMAGE NUMBER LIST                              (INPUT)
C	USE_PCA   CORAN VS PCA FLAG                              (INPUT)
C	LUNS	  SEQUENTIAL IMAGE FILE I/O UNIT (INPUT FILE)    (INPUT)
C	LUNI      IMAGE COORDINATE FILE I/O UNIT (OUTPUT FILE)   (INPUT)
C	LUNP      PIXEL COORDINATE FILE I/O UNIT (OUTPUT FILE)   (INPUT)
C	LUNE	  EIGENVALUE  FILE I/O UNIT (OUTPUT FILE)        (INPUT)
C	WEIGHTI	  SUM OF PIXEL VALUES FOR THIS IMAGE             (INPUT)
C	WEIGHTP	  SUM OF PIXEL VALUES AT THIS PIXEL              (INPUT)
C	TRANSPOSE TRANSPOSE FLAG                                 (INPUT)
C
C       ALL FILES ARE FORMATTED EXCEPT FOR LUNS!
C
C  CALL TREE:
C     JPMSK1 ---> SCORAN3 --->  INCOR3
C                             | GETCOO
C                             | GETCOOT
C
C **********************************************************************

        SUBROUTINE SCORAN3(NUMIM, NFAC, NPIX, INUMBR, USE_PCA,
     &                     LUNS, LUNI, LUNP, LUNE, 
     &                     WEIGHTI,WEIGHTP,
     &                     TRANSPOSE, SUMW)
 
        INCLUDE 'CMBLOCK.INC'

        INTEGER, DIMENSION(NUMIM)          :: INUMBR
        REAL, DIMENSION(NUMIM)             :: WEIGHTI
        REAL, DIMENSION(NPIX)              :: WEIGHTP

        REAL, ALLOCATABLE, DIMENSION(:,:)  :: MATS
        REAL, ALLOCATABLE, DIMENSION(:)    :: EVALS
        REAL, ALLOCATABLE, DIMENSION(:)    :: BLU
        REAL, ALLOCATABLE, DIMENSION(:)    :: BLCO

        REAL, DIMENSION(NFAC)              :: BLW
        LOGICAL                            :: TRANSPOSE, USE_PCA

        IF (TRANSPOSE) THEN
           WRITE(NOUT,*) ' IN-CORE SOLUTION, TRANSPOSED DATA --- '

           MWANT = NUMIM**2 + 3*NUMIM
           ALLOCATE(MATS(NUMIM,NUMIM),EVALS(NUMIM),BLU(NUMIM),
     &              BLCO(NUMIM),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'MATS',MWANT)
              RETURN
           ENDIF

           CALL INCORT(NPIX, NUMIM, NFAC, LUNS, USE_PCA, MATS,  
     &                EVALS, BLU, WEIGHTI,WEIGHTP, TRACE, SUMW,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL GETCOOT(NFAC, NUMIM, NPIX, INUMBR, USE_PCA, MATS, 
     &		     EVALS,  WEIGHTI, WEIGHTP, SUMW, BLCO, BLU, BLW,
     &		     LUNS, LUNI, LUNP)

         ELSE
           WRITE(NOUT,*) ' IN-CORE SOLUTION --- '

           MWANT = NPIX**2 + 3*NPIX
           ALLOCATE (MATS(NPIX,NPIX), EVALS(NPIX), BLU(NPIX),
     &               BLCO(NPIX), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
               CALL ERRT(46,'MATS',MWANT)
               RETURN
           ENDIF

           CALL INCOR3(NUMIM, NPIX, NFAC, LUNS, USE_PCA, MATS,  
     &                 EVALS, BLU, WEIGHTI,WEIGHTP, TRACE, SUMW,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           CALL GETCOO(NFAC, NPIX, NUMIM, INUMBR, USE_PCA, MATS,
     &		     EVALS,  WEIGHTI, WEIGHTP, SUMW, BLCO, BLU, BLW,
     &		     LUNS, LUNI, LUNP)
        ENDIF

C       FIND SIZE OF EIGENVECTS ARRAY
	N = NPIX
        IF (TRANSPOSE) N = NUMIM

        KIND_PCA = 0
        IF (USE_PCA) KIND_PCA = 1

C       WRITE EIGEN FILE HEADER
        WRITE(LUNE,90) NFAC, SUMW, TRACE, KIND_PCA, N
90      FORMAT(I10,' ',1PG12.5,' ',1PG12.5,' ',I10,' ',I10)

C       SAVE EIGENVALUES TO _EIG FILE (ONE FACTOR PER LINE)
	CUL = 0.0
	DO I = 1, NFAC
           PER = 100.0 * EVALS(I) / TRACE
           CUL = CUL + PER

           WRITE(LUNE,91) EVALS(I), PER, CUL
 91        FORMAT(1PE12.5,'  ',E12.5,'  ', E12.5)
        ENDDO

C       SAVE EIGENVECTORS ARRAY TO _EIG FILE 
	DO I = 1, N
           WRITE(LUNE,92) (MATS(I,J),J=1,N)
92         FORMAT(10(1PG12.5,' '))
        ENDDO

9999    IF (ALLOCATED(MATS))   DEALLOCATE(MATS)
        IF (ALLOCATED(EVALS))  DEALLOCATE(EVALS)
        IF (ALLOCATED(BLU))    DEALLOCATE(BLU)
        IF (ALLOCATED(BLCO))   DEALLOCATE(BLCO)

        RETURN
        END 

      
 
