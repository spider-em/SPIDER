C **********************************************************************
C
C HALI.F
C              RESTRICTION OF THE INTERPOLATION FIELD 10/13/89
C              QUADRATIC INTERPOLATION                07/09/93
C              OPFILEC                                02/24/03 al
C              CHKMIRROR                              06/18/08 al
C              HALI_P INSERTED                        03/21/12 al
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
C HALI
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE HALI(ASKMIRROR)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        LOGICAL                :: ASKMIRROR

        INTEGER, ALLOCATABLE   :: NUMR(:,:)

        CHARACTER(LEN=MAXNAM)  :: FINPAT 
        CHARACTER(LEN=1)       :: MODE,ASK
        LOGICAL                :: NORM
        LOGICAL                :: CHKMIRROR

        LOGICAL                :: FOUROK = .FALSE.

        CHARACTER(LEN=1)       :: NULL = CHAR(0)

        INTEGER,PARAMETER      :: INPIC  = 21
        INTEGER,PARAMETER      :: LUNDOC = 80   
        INTEGER,PARAMETER      :: LUNXM  = 0  ! SELFILE NOT ALLOWED

        NILMAX = NIMAX

C       OPEN INPUT IMAGE(S)
        MAXIM = 0
        CALL OPFILES(0,INPIC,LUNDOC,LUNXM, 
     &             .TRUE.,FINPAT,NLET, 'O',
     &             IFORM,NX,NY,NZ,MAXIM,
     &             'INPUT FILE TEMPLATE (E.G. PIC****)~',
     &             FOUROK,INUMBR,NILMAX, 
     &             NDUM,NIMA,IMGNUM, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN

        !CLOSE(LUNIN) bug?? july 2014
        CLOSE(INPIC)

C       NIMA - TOTAL NUMBER OF IMAGES
        IF (NIMA > 0)  THEN
           WRITE(NOUT,2001) NIMA
2001       FORMAT('  Number of images: ',I5)
        ELSE
           CALL ERRT(101,'NO IMAGES',NDUM)
           GOTO 9999
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

        CALL RDPRMC(ASK,NA,.TRUE.,'NORMALIZE UNDER MASK? (N/Y)',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        NORM = (ASK == 'Y')   

        JACUP = 0   ! UNKNOWN PURPOSE al

        CALL RDPRI2S(NKMAX,MAXIT,NOT_USED,
     &          'NUMBER OF GROUPS, MAX. NUMBER OF ITERATIONS',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

C       FIND TOTAL NUMBER OF RINGS
        NRING = 0
        DO I=MR,NR,ISKIP
           NRING = NRING + 1
        ENDDO

        ALLOCATE(NUMR(3,NRING), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'HALI; NUMR',3*NRING)
           GOTO 9999
        ENDIF

C       FILL RINGS POINTER
        NRING = 0
        DO I=MR,NR,ISKIP
           NRING         = NRING + 1
           NUMR(1,NRING) = I
        ENDDO

        CHKMIRROR = .FALSE.
        IF (ASKMIRROR) THEN
C          MUST ASK IF WANT TO CHECK MIRROR
           CALL RDPRMC(ASK,NA,.TRUE.,
     &          'CHECK MIRRORED POSITIONS? (N/Y)',NULL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           CHKMIRROR = (ASK(1:1)  == 'Y'  .OR. 
     &                  ASK(1:1) .NE. '0') 
        ENDIF
  
C       CALCULATION OF ACTUAL DIMENSION OF AN IMAGE TO BE INTERPOLATED
C       2*(No. OF RINGS)+(0'TH ELEMENT)+2*(MARGIN OF 1)

        NRA  = MIN(((NX-1)/2)*2+1, ((NY-1)/2)*2+1, 2*NR+3)
        LSAM = NX
        LROW = NY
        NX   = NRA
        NY   = NRA
        CALL ALPRBS(NUMR,NRING,LCIRC,MODE)

        MAXRIN = NUMR(3,NRING)

        CALL HALI_P(INUMBR,NX,NY,LSAM,LROW,NIMA,NRING,LCIRC,
     &        MAXRIN,JACUP,NUMR,NKMAX,MAXIT,MODE,CHKMIRROR,NORM,NOUT,
     &        FINPAT,NLET)
     
9999    IF (ALLOCATED(NUMR)) DEALLOCATE(NUMR)

        END

C ++********************************************************************
C                                                                      
C HALI_P.F                                                                   
C              OPFILEC                                02/24/03 al
C              FINPAT PARAMETER                       06/18/08 al
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
C  HALI_P                                                                    
C                                                                      
C  PURPOSE:                                                            
C                                                                      
C  PARAMETERS:                                                         
C 
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE HALI_P(ILIST,NX,NY,LSAM,LROW,NIMA,NRING,
     &      LCIRC,MAXRIN,JACUP,NUMR,NKMAX,MAXIT,MODE,
     &      MIRROR,NORM,NOUT, FINPAT,NLET)

C       BUFIN,ROT,CIRNEW,TEMP,DIST,EC,CIROLD,CIRSEED,ES,E,IP,IQ
C       ARE AUTOMATIC ARRAYS

        INCLUDE 'CMLIMIT.INC'

        REAL, ALLOCATABLE    :: X(:,:),CIRC(:,:)
        INTEGER              :: MAXRIN,MAXRI,NUMR(3,NRING)
        REAL                 :: BUFIN(LSAM),CIROLD(LCIRC),CIRNEW(LCIRC)
	INTEGER              :: KLIST(1)
        DOUBLE PRECISION     :: TEMP(MAXRIN,2),ENER,TOTMIN,TOTMIM,SOLD
        DOUBLE PRECISION     :: SNEW,EAV,EC(NIMA)
        DIMENSION            :: DIST(NIMA),ROT(NIMA)
        DIMENSION            :: ILIST(NIMA),DLIST(5)

        CHARACTER(LEN=*)     :: FINPAT
        CHARACTER(LEN=MAXNAM):: FINPIC,FINP,OUTDOC,OUTPAT
        CHARACTER(LEN=MAXNAM):: MSG
        INTEGER              :: NLET

        COMMON  /MXR/  MAXRI  ! DANGER USED IN ANG( FUNCTION al

        LOGICAL              :: CH_ANG,NORM
        CHARACTER*1          :: MODE
        LOGICAL              :: MIRROR,NEWFILE

        INTEGER,PARAMETER    :: INPIC  = 22
        INTEGER,PARAMETER    :: LUNDOC = 80   

C  --------------------------------------------
C       USED ONLY IN  HKMC
        INTEGER*2            :: IP(NIMA),IQ(NKMAX)
        DIMENSION            :: CIRSEED(LCIRC,NKMAX)
        DOUBLE PRECISION     :: E(NKMAX),ES(NKMAX)
C  --------------------------------------------

        MAXRI = MAXRIN

        LQ  = LROW/2+1
        LR1 = (NY-1)/2
        LR2 = LQ+LR1
        LR1 = LQ-LR1
        LQ  = LSAM/2+1
        LS1 = (NX-1)/2
        LS2 = LQ+LS1
        LS1 = LQ-LS1

        ALLOCATE(X(NX,NY), 
     &           CIRC(LCIRC,NIMA), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = NX*NY + LCIRC*NIMA
           CALL ERRT(46,'AP C; X & CIRC',MWANT)
           RETURN
        ENDIF

        DO K1=1,NIMA
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K1),IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(101,'FILE DOES NOT EXIST',NE)
              GOTO 9999
           ENDIF

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPIC,'O',IFORM,
     &                  NXT,NYT,NSL,
     &                  MAXIM,'DUMMY',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 9999

           DO K2=LR1,LR2
              CALL REDLIN(INPIC,BUFIN,LSAM,K2)
              DO K3=LS1,LS2
                 X(K3-LS1+1,K2-LR1+1) = BUFIN(K3)
              ENDDO
           ENDDO
           CLOSE(INPIC)

C          NORMALIZE IF REQUESTED                                       
           IF (NORM) CALL NORMAS(X,-NX/2,NX/2,-NY/2,NY/2,
     &                          NUMR,NUMR(1,NRING))

           CALL ALRQ(X,NX,NY,NUMR,CIRC(1,K1),LCIRC,NRING,MODE,K1)

           CALL FOURING(CIRC(1,K1),LCIRC,NUMR,NRING,EC(K1),MODE)
        ENDDO

C       BUILD FIRST AVERAGE

C       DIST  IS USED HERE FOR THE RANDOM CHOOSING OF IMAGES
        DIST = 0.0

        CALL RANDOM_NUMBER(CIID)
        IMI = MIN(NIMA,MAX(1,INT(CIID*NIMA+0.5)))

        CIROLD    = CIRC(:,IMI)
        ROT(IMI)  = 1.0
        DIST(IMI) = 1.0
        SOLD      = 0.0D0
        EAV       = EC(IMI)

        DO KTN=2,NIMA

804        CALL RANDOM_NUMBER(CIID) 
           M = MIN(NIMA, MAX(1,INT(CIID*(NIMA-KTN+1)+0.5)))

           IMI = 0
           DO I=1,NIMA
              IF (DIST(I) .NE. 1.0)  THEN
                 IMI = IMI + 1
                 IF (IMI .EQ. M)  GOTO  810
              ENDIF
           ENDDO
           GOTO  804

810        IMI       = I
           DIST(IMI) = 1.0

           CALL CROSRNG(CIROLD,CIRC(1,IMI),LCIRC,NRING,TEMP,TEMP(1,2),
     &                  MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)
	   IF (MIRROR)  THEN
              CALL CROSRMG(CIROLD,CIRC(1,IMI),LCIRC,NRING,
     &                     TEMP,TEMP(1,2),
     &                     MAXRIN,JACUP,NUMR,TOTMIM,TMT,MODE)
              IF (TMT > TOT)  THEN
                 ROT(IMI) = -TMT
                 SOLD     = SOLD+EAV+EC(IMI)-2.0*TOTMIM

                 CALL UPDTM(CIROLD,CIRC(1,IMI),LCIRC,NRING,NUMR,TOT,
     &                      MAXRIN,KTN)
                 GOTO 151
              ENDIF
           ENDIF

           ROT(IMI) = TOT
           SOLD     = SOLD + EAV + EC(IMI) - 2.0 * TOTMIN

           CALL UPDTC(CIROLD,CIRC(1,IMI),LCIRC,NRING,NUMR,TOT,
     &                 MAXRIN,KTN)

 151       EAV = ENER(CIROLD,LCIRC,NRING,NUMR,MODE)
            
        ENDDO


        CIRNEW = CIROLD
        ROT = 0.0

C       WRITE(NOUT,*)    SOLD*FLOAT(NIMA)/(NIMA-1)
C       WRITE(NOUT,2001) (ANG(ROT(J),MODE),J=1,NIMA)

C       ITERATIONS TO GET BETTER APPROXIMATION

        ITER = 0

901     CONTINUE
        ITER   = ITER+1
        CH_ANG = .FALSE.
        SNEW   = 0.0D0
C
        DO IMI=1,NIMA
           CALL CROSRNG(CIROLD,CIRC(1,IMI),LCIRC,NRING,TEMP,TEMP(1,2),
     &            MAXRIN,JACUP,NUMR,TOTMIN,TOT,MODE)

	   IF (MIRROR)  THEN
              CALL CROSRMG(CIROLD,CIRC(1,IMI),LCIRC,NRING,TEMP,
     &                      TEMP(1,2),MAXRIN,JACUP,NUMR,TOTMIM,TMT,MODE)
              IF (TMT .GT. TOT)  THEN
                 IF (ROT(IMI) .NE. -TMT) THEN
                    CH_ANG   = .TRUE.
                    ROT(IMI) = -TMT
                 ENDIF
                 CALL UPDTM(CIRNEW,CIRC(1,IMI),LCIRC,NRING,NUMR,
     &                  TMT,MAXRIN,IMI)
                 TOTMIN = TOTMIM
                 GOTO  152
              ENDIF
           ENDIF

           IF (ROT(IMI) .NE. TOT) THEN
              CH_ANG   = .TRUE.
              ROT(IMI) = TOT
           ENDIF
           CALL UPDTC(CIRNEW,CIRC(1,IMI),LCIRC,NRING,NUMR,TOT,
     &                MAXRIN,IMI)
152        SNEW      = SNEW+EAV+EC(IMI)-2.0*TOTMIN
           DIST(IMI) = EAV+EC(IMI)-2.0*TOTMIN
        ENDDO

        WRITE(NOUT,2020) ITER,SNEW
2020    FORMAT('  Iteration #',I3,'  Sum of distances=',1PD13.6)

        IF (SNEW.LE.SOLD .AND. CH_ANG)  THEN
           CIROLD = CIRNEW
           EAV    = ENER(CIROLD,LCIRC,NRING,NUMR,MODE)
           SOLD   = SNEW
           GOTO 901
        ENDIF

cc        IF (VERBOSE)
cc     &    WRITE(NOUT,2001) (SIGN(ANG(ABS(ROT(J)),MODE),ROT(J)),J=1,NIMA)
cc2001    FORMAT(8(1X,F8.3))

        CALL SEEDS(CIRSEED,CIRC,DIST,NKMAX,LCIRC,IP,NIMA,NOUT)

        CALL HKMC(CIRSEED,CIRC,CIRNEW,NKMAX,LCIRC,IP,IQ,ES,EC,E,
     &             DIST,ROT,NRING,TEMP,MAXRIN,JACUP,NUMR,MAXIT,
     &             NIMA,MODE,SNEW,NOUT,MIRROR)

        NMAX = 0
        CALL FILSEQP(OUTPAT,NLET,KLIST,NMAX,NIXX,
     &          'GROUP SELECTION FILE TEMPLATE',IRTFLG)

C                  123456789 123456789 123456789 123456789 123456789
        MSG = '         FILE#'

	DO II=1,NKMAX
C           MAKE GROUP SELECTION FILENAME 
	    CALL FILGET(OUTPAT,OUTDOC,NLET,II,IRTFLG)

C           OPEN GROUP SELECTION OUTPUT DOC FILE 
            CALL OPENDOC(OUTDOC,.TRUE.,NLET,LUNDOC,LUNDOCO,.FALSE.,
     &                   ' ',
     &                   .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
            IF (IRTFLG .NE. 0)  GOTO 9999

            CALL LUNDOCPUTCOM(LUNDOCO,MSG,IRTFLG)

	    KEY = 0
	    DO IIII=1,NIMA
		IF (IP(IIII) == II) THEN
	    	   KEY      = KEY + 1
		   DLIST(1) = ILIST(IIII)
                   CALL LUNDOCWRTDAT(LUNDOC,KEY,DLIST,1,IRTFLG)

		   !CALL SAVDN1(LUNDOCO,FINP,DLIST,NLS,III-1,IAP)
		ENDIF
	   ENDDO
	   CLOSE(LUNDOC)
	ENDDO

C       OPEN ALIGNMENT OUTPUT DOC FILE 
        CALL OPENDOC(OUTDOC,.TRUE.,NLET,LUNDOC,LUNDOCO,.TRUE.,
     &                   'ROTATIONAL ALIGNMENT DOC FILE',
     &              .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

C            123456789 123456789 123456789 123456789 123456789
        MSG='         FILE#         ANGLE      DIST/MIR         GROUP'
        CALL LUNDOCPUTCOM(LUNDOCO,MSG,IRTFLG)

        I = 0
        DO IMI=1,NIMA
           I = I + 1
712        IF (ILIST(I) == -1)  THEN
              I = I + 1
              GOTO 712
           ENDIF
   
           DLIST(1) = ILIST(I)
           DLIST(2) = ANG(ABS(ROT(IMI)),MODE)  ! FUNCTION CALL
           DLIST(3) = SIGN(DIST(IMI),ROT(IMI))
           DLIST(4) = IP(IMI)

           CALL LUNDOCWRTDAT(LUNDOC,IMI,DLIST,4,IRTFLG)
        ENDDO
        CLOSE(LUNDOC)

9999    IF (ALLOCATED(X))     DEALLOCATE(X)
        IF (ALLOCATED(CIRC))  DEALLOCATE(CIRC)

        END


C       ---------------------- UPDTM ------------------------------

        SUBROUTINE UPDTM(CIRC1,CIRC2,LCIRC,NRING,NUMR,TOT,MAXRIN,IS)

        DIMENSION  :: CIRC1(LCIRC),CIRC2(LCIRC)
        INTEGER    :: NUMR(3,NRING),MAXRIN
        COMPLEX    :: C

        PI2 = 8.0D0*DATAN(1.0D0)

c$omp   parallel do private(i,j,nsirt,arg,c)
        DO I=1,NRING
           NSIRT = NUMR(3,I)

           CIRC1(NUMR(2,I)) =
     &        (CIRC1(NUMR(2,I))*(IS-1)+CIRC2(NUMR(2,I)))/REAL(IS)

           CIRC1(NUMR(2,I)+1) =
     &        (CIRC1(NUMR(2,I)+1)*(IS-1)+CIRC2(NUMR(2,I)+1)*
     &        COS(PI2*(TOT-1.0)/2.0
     &        *REAL(NSIRT)/REAL(MAXRIN)))/REAL(IS)

           DO J=3,NSIRT,2
              ARG = PI2*(TOT-1.0)*REAL(J/2)/REAL(MAXRIN)

              C   = CMPLX(CIRC2(NUMR(2,I)+J-1),-CIRC2(NUMR(2,I)+J))*
     &              CMPLX(COS(ARG),SIN(ARG))

              CIRC1(NUMR(2,I)+J-1) =
     &              (CIRC1(NUMR(2,I)+J-1)*(IS-1)+REAL(C))/REAL(IS)

              CIRC1(NUMR(2,I)+J) =
     &              (CIRC1(NUMR(2,I)+J)*(IS-1)+AIMAG(C))/REAL(IS)
           ENDDO
        ENDDO

        END
