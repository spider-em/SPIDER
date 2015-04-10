
C ++********************************************************************
C                                                                      *
C SSNRB.F                                                               *
C              OPFILEC                             FEB  03 ARDEAN LEITH
C              OPFILES REWRITE                     DEC  10 ARDEAN LEITH
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C SSNRB()
C
C PARAMETERS:     
C
C OPERATION: 'RF SN'
C
C PURPOSE:  COMPUTE THE SPECTRAL SIGNAL-TO-NOISE RATIO (SSNR) OF A  
C           SERIES OF IMAGES.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE SSNRB()

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         CHARACTER(LEN=MAXNAM)         :: FILPAT1,FILMASK
         DOUBLE PRECISION, ALLOCATABLE :: FPOWER(:,:)
         DOUBLE COMPLEX,   ALLOCATABLE :: FAMPL(:,:)
         REAL, ALLOCATABLE             :: B(:,:), RMSK(:,:)
         CHARACTER(LEN=MAXNAM)         :: DOCNAM

         REAL                          :: DLIST(5)
         INTEGER, ALLOCATABLE          :: LR(:) 
         DOUBLE PRECISION, ALLOCATABLE :: SIGNAL(:), SIGDIF(:)
         CHARACTER (LEN=1)             :: NULL = CHAR(0)
         CHARACTER (LEN=100)           :: COMMENT
         LOGICAL                       :: NEWFILE
         LOGICAL                       :: FOUROK = .FALSE.
         LOGICAL                       :: ASKNAM = .TRUE.

         INTEGER, PARAMETER            :: LUN1   = 20
         INTEGER, PARAMETER            :: LUNM   = 21
         INTEGER, PARAMETER            :: LUNDOC = 80
         INTEGER, PARAMETER            :: LUNXM1 = 81

	 REAL                          :: SNR
	 REAL                          :: SIGNALSUM
	 REAL                          :: SIGDIFSUM

C        OPEN FIRST INPUT FILE, NOT FOURIER
         MAXIM1  = 0
         NILMAX  = NIMAX         ! INUMBR FROM CMLIMIT
         CALL OPFILES(0,LUN1,LUNDOC,LUNXM1,  
     &             ASKNAM,FILPAT1,NLET1, 'O',
     &             IFORM1,NX1,NY1,NZ1,MAXIM1,
     &             NULL,
     &             FOUROK, INUMBR,NILMAX, 
     &             NDUM,NGOT1,IMG1, IRTFLG) 
         IF (IRTFLG .NE. 0) RETURN

C        OPEN MASK FILE (NOT FOURIER)
         MAXIMM = 0.0
         CALL OPFILEC(0,ASKNAM,FILMASK,LUNM,'O',IFORMM,
     &                NXM,NYM,NZM,
     &                MAXIMM,'MASK',FOUROK,IRTFLG)
         IF (IRTFLG .NE. 0)  RETURN

         CALL SIZCHK(UNUSED,NX1,NY1,NZ1,0,
     &                      NXM,NYM,NZM,0,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 999

C        TOTAL NUMBER OF IMAGES
         WRITE(NOUT,2001) NGOT1
 2001    FORMAT('  NUMBER OF IMAGES: ',I5)

         WI = 1.0
         CALL RDPRM1S(WI,NOT_USED,'RING WIDTH',IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 999

         Y1     = FLOAT(MAX(NX1,NY1))
         INC    = INT(Y1/WI) / 2 + 1
         NXLD   = NX1 + 2 - MOD(NX1,2)

         ALLOCATE (RMSK  (NXLD,  NY1),
     &             FPOWER(NXLD/2,NY1), 
     &             FAMPL (NXLD/2,NY1), 
     &             B     (NXLD,  NY1),  
     &             LR    (INC),  
     &             SIGNAL(INC),  
     &             SIGDIF(INC), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = NX1*NY1  + 2 * NXLD/2*NY1 +  
     &             NXLD*NY1 + 3 * INC
           CALL ERRT(46,'RF SN, RMSK...',MWANT)
           GOTO 999
        ENDIF

C       READ MASK INTO: RMSK
        CALL READV(LUNM,RMSK,NXLD,NY1, NXM,NYM,1)

C       INITIALIZE THE SUMS
        FPOWER = 0.0D0
        FAMPL     = CMPLX(0.0,0.0)

        NINDX1 = 1

        DO                      ! LOOP OVER ALL IMAGES 

C          READ  CURRENT IMAGE INTO: B
           CALL READV(LUN1,B,NXLD,NY1, NX1,NY1,1)

C          MASK MULTIPLY
c$omp      parallel do private(i,j)
           DO J=1,NY1
              DO I=1,NX1
                 B(I,J) = B(I,J) * RMSK(I,J)
              ENDDO
           ENDDO

C          FORWARD FFT
           INV = +1
           CALL FMRS_2(B,NX1,NY1,INV)

           IF (INV .NE. 1) GOTO 999

c$omp      parallel do private(i,j)
           DO J=1,NY1

              DO I=0,NXLD/2-1

              FAMPL (I+1,J) = FAMPL(I+1,J) + 
     &                        CMPLX(B(2*I+1,J),B(2*I+2,J))
              FPOWER(I+1,J) = FPOWER(I+1,J) + 
     &                        B(2*I+1,J) * DBLE(B(2*I+1,J)) + 
     &                        B(2*I+2,J) * DBLE(B(2*I+2,J))
              ENDDO
           ENDDO

           !write(6,*)'FPOWER(15,15) =',FPOWER(15,15)
           !write(6,*)'FAMPL(15,15)  =',REAL(FAMPL(15,15)),AIMAG(FAMPL(15,15))

C          OPEN NEXT INPUT FILE, UPDATES NINDX1 
           CALL NEXTFILE(NINDX1, INUMBR, 
     &                   FOUROK, LUNXM1,
     &                   NGOT1,  MAXIM1,  
     &                   LUN1,0, FILPAT1, 'O',
     &                   IMG1,   IRTFLG)
           IF (IRTFLG .EQ. -1) EXIT      ! END OF INPUT STACK
           IF (IRTFLG .NE. 0) GOTO 999   ! ERROR  

           CALL LUNGETSIZE(LUN1,NX,NY,NZ,IRTFLG)
           CALL SIZCHK(UNUSED,NX1,NY1,  1,0,
     &                        NX, NY, NZ,0,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999   ! BAD SIZE

        ENDDO

        SIGNAL = 0.0D0  ! ARRAY INIT
        SIGDIF = 0.0D0  ! ARRAY INIT
        LR     = 0

        DO  J=1,NY1
           JJ = J-1
           IF (JJ > NY1/2)  JJ = JJ - NY1

           DO  I=1,NXLD/2
              IF (I == 1 .AND. JJ < 0) CYCLE

              PII = 0.5*SQRT((FLOAT(JJ)/ FLOAT(NY1/2))**2 +
     &                       (FLOAT(I-1)/FLOAT(NX1/2))**2)

              IF (PII <= 0.5)  THEN
                 L         = MIN(MAX(NINT(PII*Y1/WI)+1,1),INC)

                 LR(L)     = LR(L) + 1

                 SIGNAL(L) = SIGNAL(L) + ABS(FAMPL(I,J))**2

                 SIGDIF(L) = SIGDIF(L) + FPOWER(I,J) - 
     &                       ABS(FAMPL(I,J))**2 / NGOT1
              ENDIF
           ENDDO            
        ENDDO            

C       SAVE RESULTS
        CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOC,NDOC,ASKNAM,
     &           'OUTPUT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,
     &            NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999  

C                  123456789 123456789 123456789 123456789 
        COMMENT = 'RADIUS NORM.RAD        SSNR   '     //
     &            '      NVALS       VARIANCE    '
        CALL LUNDOCPUTCOM(NDOC,COMMENT(1:60),IRTFLG)

        DO L=1,INC
           IF (LR(L) .NE. 0)  THEN
              DLIST(1) = L

              DLIST(2) = FLOAT(L-1) / FLOAT(INC-1) * 0.5

C              DLIST(3) = AMAX1(0.0, SNGL((FLOAT(NGOT1-1) / NGOT1**2) *
C     &                   SIGNAL(L) / SIGDIF(L)) - 1.0)
              DLIST(3) = AMAX1(0.0, SNGL((FLOAT(NGOT1-1) / NGOT1) *
     &                   SIGNAL(L) / SIGDIF(L)) - 1.0)

              DLIST(4) = LR(L)

C             N1       = LR(L)
C             N2       = (NGOT1-1) * LR(L)
C             DLIST(5) = SQRT(2.*((N2+N1-2) * (1+2*DLIST(3)) +
C     &                   N1 * DLIST(3)**2) / FLOAT(N1*(N2-4)))

              DLIST(5) = SQRT((2 + 4*DLIST(3))/NGOT1 +
     &        (2 + 4*DLIST(3) + 
     &         2*DLIST(3)**2)/(NGOT1*(LR(L)-1)))


              CALL LUNDOCWRTDAT(NDOC,L,DLIST(2),4,IRTFLG)
           ENDIF
        ENDDO

	SIGNALSUM = 0
	SIGDIFSUM = 0

        SNR = 0
        DO L=1,INC
	   SIGNALSUM = SIGNALSUM + SIGNAL(L)
	   SIGDIFSUM = SIGDIFSUM + SIGDIF(L)
        ENDDO

        SNR = (SIGNALSUM/SIGDIFSUM) / FLOAT(NGOT1)
        WRITE(NOUT,*) '  Integral SNR= ',SNR

999     CLOSE(LUNDOC)
        CLOSE(LUNM)
        CLOSE(LUN1)

        IF (ALLOCATED(FPOWER)) DEALLOCATE (FPOWER) 
        IF (ALLOCATED(FAMPL))  DEALLOCATE (FAMPL)
        IF (ALLOCATED(RMSK))   DEALLOCATE (RMSK) 
        IF (ALLOCATED(B))      DEALLOCATE (B)

        END 
