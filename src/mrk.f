
C ++********************************************************************
C                                                                      
C   MRK
C          PROMPTS & DOC FILE HEADERS IMPROVED   FEB 2014 ARDEAN LEITH
C                                                                      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C   MRK(IDUM)
C
C   WARNING ! PARAMETERS LS AND LV ARE HARD-WIRED IN ALL SUBROUTINES.
C
C   PURPOSE:  MAIN SUBROUTINE FOR MARKER PROGRAM.
C     ATTEMPTS TO ALIGN A DATA SET USING FIDUCIARY MARKERS.
C
C     THE PROCEDURE IS TO FIRST ROUGHLY TRANSLATIONALLY
C     ALIGN THE IMAGES. USING THIS SET, RECONSTRUCT THE
C     COORDS OF THE POINTS IN THE 3D OBJECT. THEN PROJECT
C     THE POINTS BACK DOWN ONTO A PLANE AND COMPARE THEM
C     WITH THE POINTS FROM PREVIOUS TRIAL. CORRECT THE
C     POINTS AND ONCE AGAIN RECONSTRUCT THE 3D, ETC.
C
C     THE CORRECTIONS ACCUMULATE IN A FEW REGISTERS.
C     ALL COMPARISONS EXCEPT FOR THE ROUGH ALIGNMENT ARE
C     BETWEEN THE ORIGINAL PROJECTIONS AND THEIR PSUEDO-
C     PROJECTIONS.
C
C **********************************************************

      SUBROUTINE MRK(IDUM)

      INCLUDE     'CMBLOCK.INC'
      INCLUDE     'CMLIMIT.INC'

      INTEGER, PARAMETER   :: LV=300
      INTEGER, PARAMETER   :: LS=256
      INTEGER, PARAMETER   :: MAXLOOP = 19999
      REAL, PARAMETER      :: MAXPER  = 1.0E-7
      INTEGER, PARAMETER   :: NLIST   = 2 

      CHARACTER(LEN=MAXNAM):: DOCNAM

      DIMENSION         DLIST(NLIST)
      LOGICAL*1         FULL
      CHARACTER*1       ASKFULL,NULL

      LOGICAL           PARAMQ(4)
      LOGICAL           PTACTIVE(LS,LV)
      INTEGER           NUMPTS(LV)
      REAL              CIR(2)

      DIMENSION         XYPTS(2,LS,LV),PRJ(2,LS,LV)
      DIMENSION         ROT(3,3),ANGLE(3,LV),TPSI(LV)
      DIMENSION         P3D(3,LS),SCALE(LV)
      DIMENSION         NUMBER(LV),SCALEI(LV)
      DIMENSION         CVPT(3)
      DIMENSION         SHFT(2),SHIFT(2,LV),TSHIFT(2,LV)
      DIMENSION         ERVW(LV),ERPT(LS)

      LOGICAL           FIRST
      INTEGER           NTIMES
      DOUBLE PRECISION  TCH,ERTOT

      NDOUT = 70

      NTIMES = 0
      NULL   = CHAR(0)

C     GET COORDS OF MARKERS IN EACH PROJECTION

      CALL MRGETINFO(XYPTS,IREF,NUMBER,ANGLE,SCALE,TSHIFT,FIRST,NTTVL,
     &               PARAMQ,PTACTIVE,NUMPTS,NTVW,NTPT,CIR)

      DO  IVIEW=1,NTVW

        TPSI(IVIEW)    = ANGLE(1,IVIEW)
        SHIFT(1,IVIEW) = TSHIFT(1,IVIEW)
        SHIFT(2,IVIEW) = TSHIFT(2,IVIEW)

        CALL MRALIGN(XYPTS(1,1,IVIEW),
     &	      ANGLE(1,IVIEW),SHIFT(1,IVIEW),SCALE(IVIEW),
     &	      PTACTIVE(1,IVIEW),NTPT)
      ENDDO

C     GET 3D VOLUME
      CALL MR2TO3D(P3D, XYPTS, ANGLE,PTACTIVE,NUMPTS,NTVW,NTPT)
      CALL MRPROJ(P3D, PRJ, ANGLE,NTVW,NTPT)
      CALL MRCALERR(XYPTS,PRJ,ERTOT,ERVW,ERPT,PTACTIVE,NUMPTS,NTVW,NTPT)

      WRITE(NOUT,*) ' INITIAL INFO:'

      FULL = .FALSE.
      CALL MRPUTINFO(XYPTS,PRJ,ANGLE,P3D,
     &     NUMBER,TSHIFT,SCALE,ERVW,ERPT,FULL,PTACTIVE,NTVW,NTPT,CIR)

      TCH = ERTOT
      WRITE(NOUT,701) NTIMES,ERTOT
 701  FORMAT('  PASS #',I5,' TOTCHANGE: ',10('-'),' ERTOT:',G11.4,/)

C     CREATE PSUEDO-PROJECTIONS AND COMPARE TO LAST RUN

 200  NTIMES = NTIMES + 1

C     DO ONLY IF _IN PLANE_ ROTATION OR _SCALE_ CORRECTIONS REQUESTED.

      IF (PARAMQ(1) .OR. PARAMQ(3))  THEN
	 SHFT(1) = 0.0
	 SHFT(2) = 0.0

         DO  JVIEW=0,NTVW
C          THIS STRANGE TRICK IS TO MAKE SURE THAT FIRST PASS IS OVER 
C          REFERENCE IMAGE, IN SUBSEQUENT ITERATIONS REFERENCE IMAGE 
C          IS SKIPPED.
           IVIEW = JVIEW
           IF (JVIEW .EQ. IREF) CYCLE

           IF (JVIEW .EQ. 0) IVIEW = IREF

C          ASSUMPTION IS THAT TILT OF PROJ AND TILT OF ORIG ARE EQUAL
C          THIS IS THE ONLY REASON WE CAN REALLY COMPARE THE TWO

C          FIND PSI (IN PLANE) ANGLE
           IF (PARAMQ(1)) THEN
              CALL MRANG2(PRJ(1,1,IVIEW),XYPTS(1,1,IVIEW),IVIEW,PSI,
     &                    PTACTIVE,NTPT)
	   ELSE
	      PSI = 0.0
	   ENDIF

C          FIND SCALE
           IF (PARAMQ(3)) THEN
             CALL MRSCALE(PRJ(1,1,IVIEW),XYPTS(1,1,IVIEW),IVIEW,SCALEI,
     &                    PTACTIVE,NUMPTS)
             SCAL = SCALEI(IVIEW)
             IF (IVIEW .EQ. IREF) SCAD = 1 / SCALEI(IVIEW)
             SCALEI(IVIEW) = SCALEI(IVIEW) * SCAD
           ELSE
	     SCAL          = 1.0
             SCALEI(IVIEW) = 1.0
           ENDIF

C          UPDATE CORRECTION INFO IN REGISTERS, NO SHIFTS AT THIS MOMENT
           C               = COS(PSI)
           S               = SIN(PSI)
           QT              = TSHIFT(1,IVIEW)
           TSHIFT(1,IVIEW) = SCALEI(IVIEW) * ( QT*C+TSHIFT(2,IVIEW)*S)
           TSHIFT(2,IVIEW) = SCALEI(IVIEW) * (-QT*S+TSHIFT(2,IVIEW)*C)
           TPSI(IVIEW)     = PSI + TPSI(IVIEW)
           SCALE(IVIEW)    = SCALE(IVIEW) * SCALEI(IVIEW)
           ANGLE(1,IVIEW)  = PSI

C          APPLY CORRECTIONS TO POINTS WITH SHFT SET TO ZERO ...
           CALL MRALIGN(XYPTS(1,1,IVIEW),
     &	            PSI,SHFT,SCALEI(IVIEW),PTACTIVE(1,IVIEW),NTPT)

        ENDDO   ! ------------------------------------------------

C       NOW GOTO 3D
	CALL MR2TO3D(P3D, XYPTS, ANGLE,PTACTIVE,NUMPTS,NTVW,NTPT)
      ENDIF


      IF (PARAMQ(4)) THEN
C       FIND TRANSLATION (SHIFT)

C       FIND CENTER OF GRAVITY OF POINTS IN 3D
C       AND APPLY ANY SHIFTS FOUND TO THESE POINTS
	CALL MRCG3D(P3D,NTPT)

C       GOTO 2D
	CALL MRPROJ(P3D, PRJ, ANGLE,NTVW,NTPT)

C       FIND THE SHIFT BETWEEN POINTS AND SHIFT XYPTS
	CALL MRSHIFT(PRJ,XYPTS,SHIFT,PTACTIVE,NTVW,NTPT)

C       UPDATE PARAMETERS
	DO IVIEW=1,NTVW
           TSHIFT(1,IVIEW) = TSHIFT(1,IVIEW) + SHIFT(1,IVIEW)
           TSHIFT(2,IVIEW) = TSHIFT(2,IVIEW) + SHIFT(2,IVIEW)
	ENDDO

C       GOTO 3D AGAIN
        CALL MR2TO3D(P3D,XYPTS,ANGLE,PTACTIVE,NUMPTS,NTVW,NTPT)
      ENDIF


      IF (PARAMQ(2)) THEN
C        FIND THETA (TILT) ANGLE

C        GOTO 2D
	 CALL MRPROJ(P3D, PRJ, ANGLE,NTVW,NTPT)
	 DO JVIEW=0,NTVW

C          THIS STRANGE TRICK IS TO MAKE SURE THAT FIRST PASS IS OVER REFERENCE
C          IMAGE, IN SUBSEQUENT ITERATIONS REFERENCE IMAGE IS SKIPPED.
C          WON'T MODIFY XYPTS, ONLY THETA

           IVIEW = JVIEW
           IF (JVIEW .NE. IREF) THEN
             IF (JVIEW .EQ. 0) IVIEW=IREF
             THETA = ANGLE(2,IVIEW)

C            SET SCALE=1.0 AND PSI=0.0
             CALL MRTHETA(PRJ(1,1,IVIEW),XYPTS(1,1,IVIEW),
     &                    IVIEW,P3D,THETA,PTACTIVE,NUMPTS,NTPT)
             IF (IVIEW .EQ. IREF) THAD = THETA
             THETA = THETA - THAD

C            CORRECTION OF THETA IS FOUND, SO CALCULATE NEW THETA
             ANGLE(2,IVIEW) = THETA
	  ENDIF
	ENDDO

C       GOTO 3D AGAIN
	CALL MR2TO3D(P3D, XYPTS, ANGLE,PTACTIVE,NUMPTS,NTVW,NTPT)
      ENDIF

C     GOTO 2D AND GET THE ERROR ...
      CALL MRPROJ(P3D,PRJ,ANGLE,NTVW,NTPT)
      CALL MRCALERR(XYPTS,PRJ,ERTOT,ERVW,ERPT,PTACTIVE,NUMPTS,NTVW,NTPT)

      TTCH = SNGL((TCH-ERTOT)/TCH)
      TCH  = ERTOT
      IF (MOD(NTIMES,150) .EQ. 0)  THEN 
           WRITE(NOUT,700)NTIMES,TTCH,ERTOT
 700       FORMAT('  PASS #',I5,' TOTCHANGE= ',G10.3,' ERTOT= ',G11.4)

	   FULL = .FALSE.
     	   CALL MRPUTINFO(XYPTS,PRJ,ANGLE,P3D,
     &       NUMBER,TSHIFT,SCALE,ERVW,ERPT,FULL,PTACTIVE,NTVW,NTPT,CIR)
	ENDIF

C       DO THE NEXT ITERATION ?
	IF (ABS(TTCH) .GT. MAXPER .AND. NTIMES.LT.MAXLOOP) GOTO 200

C       LOOPS ABOVE --------------------------------------------------


	IF (NTIMES.GE.MAXLOOP)
     &	   WRITE(NOUT,*)'  EXCEEDED LOOP LIMIT, ENDING CYCLE'

C       JUST DO SOME CLEANING UP BEFORE PUTTING OUT VALUES
	TILR = ANGLE(2,IREF)

	DO  IVIEW=1,NTVW
           ANGLE(1,IVIEW) = TPSI(IVIEW)
           ANGLE(2,IVIEW) = ANGLE(2,IVIEW)-TILR
	ENDDO

	CALL  RDPRMC(ASKFULL,NUMC,.TRUE.,'FULL OUTPUT (Y/N)',NULL,IRT)
	FULL = (ASKFULL .NE. 'N')

	CALL MRPUTINFO(XYPTS,PRJ,ANGLE,P3D,
     &      NUMBER,TSHIFT,SCALE,ERVW,ERPT,FULL,PTACTIVE,NTVW,NTPT,CIR)

987     CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOUT,NDOUTT,.TRUE.,
     &               'ERROR PER VIEW OUTPUT DOC',
     &               .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO  51
        IF (NDOUTT .LE. 0) THEN
           CALL ERRT(102,'INCORE DOC FILE NOT ALLOWED HERE.',NE)
           GOTO 987
        ENDIF 
        CALL LUNDOCPUTCOM(NDOUTT,'   ERROR_PER_VIEW', IRTFLG)
       
	DO  IVIEW=1,NTVW
	   DLIST(1) = ERVW(IVIEW)
           CALL LUNDOCWRTDAT(NDOUTT,NUMBER(IVIEW),DLIST,1,IRTFLG)
	ENDDO

        CALL LUNDOCPUTCOM(NDOUTT,'        TOTAL_ERROR', IRTFLG)
C       DLIST(1) = -1
	DLIST(2) = ERTOT
        CALL LUNDOCWRTDAT(NDOUTT,-1,DLIST(2),1,IRTFLG)

51	CLOSE(NDOUT)

986     CALL OPENDOC(DOCNAM,.TRUE.,NLET,NDOUT,NDOUTT,.TRUE.,
     &               'ERROR PER MARKER OUTPUT DOC',
     &               .FALSE.,.FALSE.,.TRUE.,NEWFILE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO  52
        IF (NDOUT .LE. 0) THEN
           CALL ERRT(102,'INCORE DOC FILE NOT ALLOWED HERE.',NE)
           GOTO 986
        ENDIF 

        CALL LUNDOCPUTCOM(NDOUTT,'       TOTAL_ERROR', IRTFLG)
	DO  IPT=1,NTPT
	   DLIST(1) = ERPT(IPT)
           CALL LUNDOCWRTDAT(NDOUTT,IPT,DLIST,1,IRTFLG)
	ENDDO

        CALL LUNDOCPUTCOM(NDOUTT,'          ERROR_PER_PT', IRTFLG)
	DLIST(1) = ERTOT
        CALL LUNDOCWRTDAT(NDOUTT,-1,DLIST,1,IRTFLG)
52      CLOSE(NDOUT)

	CALL MRDOC(SCALE,TSHIFT,ANGLE,NUMBER,P3D,NTVW,NTPT)

	END

C------------- UNUSED BELOW ---------------------------

