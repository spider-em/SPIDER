C++*********************************************************************
C
C FALB
C              B-SPLINE INTERPOLATION INTRODUCED              09/21/89
C              RESTRICTION OF THE INTERPOLATION FIELD         10/13/89
C              SUBTRACTION OF ONE IMAGE                       03/27/91 
C              QUADRATIC INTERPOLATION USED AS AN OPTION      06/24/91
C              SCRATCH FILE ON THE DISK                       08/01/91
C              PROMPTS                             JAN 02 ARDEAN LEITH
C              OPFILEC                             FEB 03 ARDEAN LEITH
C              FALB_P INSERTED                     MAR 12 ARDEAN LEITH
C
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
C
C         FALB
C         FALB_P(BUF,ILIST,NX,NY,LSAM,LROW,NIMA,
C         ANG(RKK,MODE)
C         ENER(CIRC,LCIRC,NRING,NUMR,MODE)
C         ALPRBS(NUMR,NRING,LCIRC,MODE)
C         ALRQ
C         UPDTC(CIRC1,CIRC2,LCIRC,NRING,NUMR,TOT,MAXRIN,IS)
C         OUTRNG
C         CROSRNG
C         FOURING(CIRC,LCIRC,NUMR,NRING,E,MODE)
C         LOG2(N)
C         PRB1D(B,NPOINT,POS)
C         FFTR_D(X,NV)   
C         FFTC_D(BR,BI,LN,KS)
C
C         USES NON-FFTW FOURIER
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE FALB

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'   

        INTEGER, ALLOCATABLE  :: NUMR(:,:)

        CHARACTER(LEN=MAXNAM) :: FINPAT,FINPIC

        INTEGER               :: MAXRIN
        CHARACTER(LEN=1)      :: MODE,ASK
        CHARACTER(LEN=1)      :: NULL = CHAR(0) 
        LOGICAL               :: FOUROK = .FALSE.

        INTEGER, PARAMETER    :: INPIC = 77
        INTEGER,PARAMETER     :: LUNXM = 0  ! SELFILE NOT ALLOWED
        INTEGER               :: IMGNUM

C       OPEN INPUT IMAGE(S)
        MAXIM  = 0
        NILMAX = NIMAX
        IMGNUM = 0    ! Passing uninitialized variable may result in 'INVALID IMAGE NUMBER' error
        CALL OPFILES(0,INPIC,LUNDOC,LUNXM, 
     &             .TRUE.,FINPAT,NLET, 'O',
     &             IFORM,NX,NY,NZ,MAXIM,
     &             'INPUT FILE TEMPLATE (E.G. PIC****)~',
     &             FOUROK,INUMBR,NILMAX, 
     &             NDUM,NIMA,IMGNUM, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)

C       NIMA - TOTAL NUMBER OF IMAGES
        IF (NIMA > 0)  THEN
           WRITE(NOUT,2001) NIMA
2001       FORMAT('  Number of images: ',I5)
        ELSE
           CALL ERRT(100,'NO IMAGES',NDUM)
           RETURN
        ENDIF
	 
        MR    = 5
        NRAD  = MIN(NX/2-1, NY/2-1)
        NR    = NRAD
        ISKIP = -1    ! TRAP FOR OLD STYLE INPUT

        CALL RDPRI3S(MR,NR,ISKIP,NOT_USED,
     &       'FIRST, LAST RING, & RING SKIP,',IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        IF (ISKIP < 0) THEN
            CALL RDPRI1S(ISKIP,NOT_USED,'RING SKIP',IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 9999
        ENDIF
        ISKIP = MAX(1,ISKIP)

	IF (MR <= 0) THEN
	   CALL ERRT(101,'FIRST RING MUST BE > 0',NE)
	   GOTO 9999

	ELSEIF (NR < MR)  THEN 
	   CALL ERRT(102,'LAST RING MUST BE > ',MR)
	   GOTO 9999

	ELSEIF (NR >= NRAD)  THEN 
	   CALL ERRT(102,'LAST RING MUST BE < ',NRAD)
	   GOTO 9999
        ENDIF

        CALL  RDPRMC(ASK,NA,.TRUE.,
     &               'ANALYZE FULL OR HALF RING? (F/H)',
     &                NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (ASK == 'F')  THEN
           MODE = 'F'
        ELSEIF (ASK == 'H')  THEN
           MODE = 'H'
        ELSE
	   CALL ERRT(101,'INVALID RESPONSE',NRAD)
	   GOTO 9999
        ENDIF

C       CALL  RDPRMI(JACUP,NDUMP,NOT_USED,
C     &         'Precision of peak location (0..100)')
C       JACUP=MAX0(0,MIN0(100,JACUP))
        JACUP=0

C       FIND TOTAL NUMBER OF RINGS
        NRING = 0
        DO I=MR,NR,ISKIP
           NRING = NRING + 1
        ENDDO

        ALLOCATE (NUMR(3,NRING), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'FALB; NUMR',3*NRING)
           GOTO 9999
        ENDIF

C       FILL RINGS POINTER
        NRING = 0
        DO I=MR,NR,ISKIP
           NRING         = NRING + 1
           NUMR(1,NRING) = I
        ENDDO
 
C       CALCULATION OF ACTUAL DIMENSION OF AN IMAGE TO BE INTERPOLATED
C       2*(NO. OF RINGS)+(0'TH ELEMENT)+2*(MARGIN OF 1)

        NRA  = MIN(((NX-1)/2)*2+1, ((NY-1)/2)*2+1, 2*NR+3)
        LSAM = NX
        LROW = NY
        NX   = NRA
        NY   = NRA

C       FILL NUMR
        CALL ALPRBS(NUMR,NRING,LCIRC,MODE)

        MAXRIN = NUMR(3,NRING)  ! CORRECT AS NOT FFTW al
 
        CALL  FALB_P(INUMBR,NX,NY,LSAM,LROW,NIMA,
     &               NRING,LCIRC,MAXRIN,JACUP,NUMR,MODE, FINPAT,NLET)

9999    IF (ALLOCATED(NUMR)) DEALLOCATE(NUMR)

        END


C++*********************************************************************
C
C FALB_P.F                ROT FIXED & RANDOMIZED JULY 2000 ARDEAN LEITH
C
C **********************************************************************
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE  FALB_P(ILIST,NX,NY,LSAM,LROW,NIMA,
     &          NRING,LCIRC,MAXRIN,JACUP,NUMR,MODE, FINPAT,NLET)

C       BUFIN,ROT,CIROLD,CIRNEW,CIRTMP,TEMP,EC ARE AUTOMATIC ARRAYS
	
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'   

        REAL, ALLOCATABLE  :: X(:,:),CIRC(:,:) 

        INTEGER            :: MAXRIN,MAXRI,NUMR(3,NRING)
        REAL               :: BUFIN(LSAM),DIST(NIMA)
        REAL               :: CIROLD(LCIRC),CIRNEW(LCIRC),CIRTMP(LCIRC)
        DOUBLE PRECISION   :: TEMP(MAXRIN,2)
        DOUBLE PRECISION   :: ENER,TOTMIN,SOLD,SNEW,EAV,EC(NIMA)
        INTEGER            :: ILIST(NIMA)
        REAL               :: DLIST(5),ROT(NIMA)
        LOGICAL            :: NEWFILE

        INTEGER            :: NLET
        LOGICAL*1          :: CH_ANG
        CHARACTER(LEN=1)   :: MODE

        CHARACTER(LEN=MAXNAM) :: FINPIC,FINPAT,OUTDOC
        CHARACTER(LEN=MAXNAM) :: MSG

        COMMON  /MXR/ MAXRI
C       MXR in: ang.f,gali.f,hali_p.f

        INTEGER,PARAMETER    :: INPIC  = 77
        INTEGER,PARAMETER    :: LUNDOC = 80   

        MAXRI = MAXRIN

        LQ  = LROW/2+1
        LR1 = (NY-1)/2
        LR2 = LQ+LR1
        LR1 = LQ-LR1
        LQ  = LSAM/2+1
        LS1 = (NX-1)/2
        LS2 = LQ+LS1
        LS1 = LQ-LS1

        ALLOCATE (X(NX,NY), CIRC(LCIRC,NIMA), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           MWANT = NX*NY + LCIRC*NIMA 
           CALL ERRT(46,'FALB_P; X',MWANT)
           RETURN
        ENDIF

        DO K1=1,NIMA

           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K1),INTFLAG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'FILE DOES NOT EXIST',NE)
              GOTO 9999
           ENDIF

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,
     &           NXT,NYT,NSL,MAXIM,'DUMMY',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999

           DO K2=LR1,LR2
              CALL  REDLIN(INPIC,BUFIN,LSAM,K2)
              DO K3=LS1,LS2
                 X(K3-LS1+1,K2-LR1+1) = BUFIN(K3)
              ENDDO
           ENDDO
       
           CLOSE(INPIC)

           CALL ALRQ(X,NX,NY,NUMR,CIRC(1,K1),LCIRC,NRING,MODE,K1)
           CALL FOURING(CIRC(1,K1),LCIRC,NUMR,NRING,EC(K1),MODE)

        ENDDO

C       BUILD FIRST AVERAGE

C       TWO ESTIMATION OF INITIAL AVERAGE ARE USED
C       ONLY ONE !!!  11/06/91

C          DIST IS USED HERE FOR THE RANDOM CHOOSING OF IMAGES

           DO IMI=1,NIMA
              DIST(IMI) = 0.0
           ENDDO

           CALL RANDOM_NUMBER(CIID)
           IMI = MIN(NIMA,MAX0(1, INT(CIID*NIMA+0.5)))

           DO I=1,LCIRC
              CIROLD(I) = CIRC(I,IMI)
           ENDDO

           ROT(IMI)  = 1.0
           DIST(IMI) = 1.0

C          write(nout,*) 'rot(',imi,'):',rot(imi),ciid,nima,dist(imi)
           DO KTN=2,NIMA

804           CALL RANDOM_NUMBER(CIID) 
              M = MIN(NIMA, MAX(1, INT(CIID*(NIMA-KTN+1)+0.5)))

              IMI = 0
              DO I=1,NIMA
                 IF (DIST(I) .NE. 1.0)  THEN
                    IMI = IMI + 1
                    IF (IMI == M)  GOTO 810
                 ENDIF
809           CONTINUE
              ENDDO
              GOTO  804

810           IMI       = I
              DIST(IMI) = 1.0

              CALL CROSRNG(CIROLD,CIRC(1,IMI),LCIRC,NRING,TEMP,
     &             TEMP(1,2),MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)

              ROT(IMI) = TOT

C             write(nout,*) 'rot(',imi,'):',rot(imi),ciid,m,dist(imi),i
              CALL UPDTC(CIROLD,CIRC(1,IMI),LCIRC,NRING,NUMR,
     &                   TOT,MAXRIN,KTN)
           ENDDO

           SOLD = ENER(CIROLD,LCIRC,NRING,NUMR,MODE)
           WRITE(NOUT,2037)  SOLD
2037      FORMAT('  Random approximation: 1,  Squared sum:    ',1PE12.5)

           DO I=1,LCIRC
              CIRNEW(I) = CIROLD(I)
           ENDDO
                 
C       PRINT  *,SOLD*FLOAT(NIMA)/(NIMA-1)
C       PRINT  2001,(ANG(ROT(J),MODE),J=1,NIMA)
2001    FORMAT(8(1X,F8.3))

C       ITERATIONS TO GET BETTER APPROXIMATION

        ITER = 0
901     CONTINUE

        ITER   = ITER+1
        CH_ANG = .FALSE.

        DO IMI=1,NIMA
           CALL OUTRNG(CIROLD,CIRC(1,IMI),CIRTMP,LCIRC,NRING,
     &                 NUMR,ROT(IMI),MAXRIN,NIMA)

           CALL CROSRNG(CIRTMP,CIRC(1,IMI),LCIRC,NRING,TEMP,
     &               TEMP(1,2),MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)

           IF (ROT(IMI) .NE. TOT)  THEN
              CH_ANG   = .TRUE.
              ROT(IMI) = TOT
           ENDIF

          CALL UPDTC(CIRNEW,CIRC(1,IMI),LCIRC,NRING,NUMR,TOT,MAXRIN,IMI)

        ENDDO

        SNEW = ENER(CIRNEW,LCIRC,NRING,NUMR,MODE)

        WRITE(NOUT,2030)  ITER,SNEW
2030    FORMAT('  Iteration:         ',I4,',  New squared sum:',1PE12.5)

        IF (SNEW .GE. SOLD .AND. CH_ANG)  THEN
           DO  I=1,LCIRC 
              CIROLD(I) = CIRNEW(I)
           ENDDO
           SOLD = SNEW
           GOTO  901
        ENDIF

        !IF (VERBOSE) WRITE(NOUT,*1) ' ANGLES:'
        !IF (VERBOSE) WRITE(NOUT,2001) (ANG(ROT(J),MODE),J=1,NIMA)

C       OPEN ALIGNMENT OUTPUT DOC FILE 
        CALL OPENDOC(OUTDOC,.TRUE.,NLET,LUNDOC,LUNDOCO,.TRUE.,
     &                   'CLASSIFICATION & ALIGNMENT DOC FILE',
     &              .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

C            123456789 123456789 123456789 123456789 123456789
        MSG='         FILE#         ANGLE      DIST/MIR         GROUP'
        CALL LUNDOCPUTCOM(LUNDOCO,MSG,IRTFLG)

        DLIST(4) = 1.0   ! FOR CONSISTENCY WITH 'AP C'

        DO   IMI=1,NIMA
C          CALCULATE DISTANCES
           CALL  OUTRNG(CIROLD,CIRC(1,IMI),CIRTMP,LCIRC,NRING,
     &           NUMR,ROT(IMI),MAXRIN,NIMA)

           CALL CROSRNG (CIRTMP,CIRC(1,IMI),LCIRC,NRING,TEMP,
     &          TEMP(1,2), MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)

           EAV      = ENER(CIRTMP,LCIRC,NRING,NUMR,MODE)

           DLIST(1) = ILIST(IMI)
           DLIST(2) = ANG(ROT(IMI),MODE)  ! FUNCTION CALL
           DLIST(3) = EAV + EC(IMI) - 2.0 * TOTMIN

           CALL LUNDOCWRTDAT(LUNDOC,IMI,DLIST,4,IRTFLG)

           !!CALL SAVD(NDOC,DLIST,5,IRTFLG)
        ENDDO

        CLOSE(LUNDOC)

9999    IF (ALLOCATED(X))    DEALLOCATE(X)
        IF (ALLOCATED(CIRC)) DEALLOCATE(CIRC)

        END
