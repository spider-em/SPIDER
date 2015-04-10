C++*********************************************************************
C
C    ORMD.F 
C             OPFILEC                          FEB 2003  ARDEAN LEITH
C             REWRITTEN                        JUN 2008  ARDEAN LEITH
C             APRINGS_INIT_PLANS PARAMS        JUN 2011  ARDEAN LEITH
C             ANGLE BUG                        JAN 2013  ARDEAN LEITH                      
C
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
C
C  ORMD()
C
C PURPOSE: DETERMINES ROTATIONAL ORIENTATION BETWEEN TWO IMAGES USING 
C          RESAMPLING INTO POLAR COORDINATES. SAME AS: 
C          'AP REF' OR 'AP SH' WITHOUT THE TRANSLATIONAL SEARCH AND
C          IT ONLY ACTS ON A SINGLE PAIR OF IMAGES. 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ORMD(ASKPEAKS)

        INCLUDE 'SETNUMRINGS.INC'
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=MAXNAM)    :: FILNAM,REFNAM,OUTANG
        INTEGER, ALLOCATABLE      :: NUMR(:,:)
        REAL, ALLOCATABLE         :: X(:,:)
        CHARACTER (LEN=1)         :: ASK,MODE
        LOGICAL                   :: NEWFILE 
	CHARACTER(LEN=80)         :: COMMEN
        LOGICAL                   :: ASKPEAKS
        CHARACTER (LEN=1)         :: NULL = CHAR(0)

        INTEGER, PARAMETER        :: NPLANS = 14
        INTEGER *8                :: FFTW_PLANS(NPLANS)

        INTEGER, PARAMETER        :: LUNREF   = 20
        INTEGER, PARAMETER        :: LUNEXP   = 21
        INTEGER, PARAMETER        :: LUNDOC   = 70 

C       ASK FOR INPUT FILES
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNEXP,'O',ITYPE,NX,NY,
     &               NZ,MAXIM,'EXPERIMENTAL IMAGE~',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,REFNAM,LUNREF,'O',ITYPE,NX1,NY1,
     &               NZ1,MAXIM,'REFERENCE IMAGE~',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL SIZCHK(UNUSED, NX,NY,NZ,0,  
     &                      NX1,NY1,NZ1,0, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ISKIP = 0
        CALL RDPRI3S(MR,NR,ISKIP,NOT_USED,
     &              'FIRST, LAST RING, & RING STEP',IRTFLG)
        IF (ISKIP .LE. 0) THEN
           ISKIP = 1
           CALL RDPRI1S(ISKIP,NOT_USED,'RING STEP',IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9999
        ISKIP = MAX(1,ISKIP)

        NRAD = MIN(NX/2-1, NY/2-1) -1

	IF (MR .LE. 0) THEN
	   CALL ERRT(101,'FIRST RING MUST BE > 0',NE)
	   GOTO 9999

	ELSEIF (NR .GE. NRAD)  THEN 
	   CALL ERRT(102,'LAST RING MUST BE < ',NRAD)
	   GOTO 9999
        ENDIF

        MODE = 'F'
        CALL RDPRMC(ASK,NA,.TRUE.,'FULL OR HALF CIRCLE (F/H)',
     &              NULL,IRT)
        IF (IRT .NE. 0) GOTO 9999
        IF (ASK == 'H') MODE = 'H'

        NPEAK = 1
        IF (ASKPEAKS) THEN 
           CALL RDPRI1S(NPEAK,NOT_USED,'NUMBER OF PEAKS',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           NPEAK = MAX(1,NPEAK)
        ENDIF
        
        IF (NPEAK > 1) THEN
C          OPEN PEAK OUTPUT DOC FILE (FOR APPENDING)
           NOUTANG = LUNDOC
           CALL OPENDOC(OUTANG,.TRUE.,NLET,LUNDOC,NOUTANG,.TRUE.,
     &           'OUTPUT ANGLE DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &            NEWFILE,IRTFLG)

           COMMEN = '       ANGLE,     PEAK-HEIGHT '
           CALL LUNDOCPUTCOM(LUNDOC,COMMEN,IRTFLG)
        ENDIF

C       CREATES NUMR ARRAY HOLDING THE SPECS FOR RADIAL RINGS
        CALL SETNUMRINGS(MR,NR,ISKIP,MODE, 1,
     &                   NUMR,NRING,LCIRC,
     &                   IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999
        MAXRAYS = NUMR(3,NRING) -2  ! ACTUAL LENGTH OF LONGEST RING

        ALLOCATE (X(NX,NY), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'ORMD, X',NX*NY)
           GOTO 9999
        ENDIF

C       INITIALIZE FFTW3 PLANS FOR USE WITHIN OMP || SECTIONS
        CALL APRINGS_INIT_PLANS(NUMR,NRING,
     &                          FFTW_PLANS,NPLANS,NX,NY,IRTFLG)

        CALL ORMD_P(NX,NY, 
     &      NRING,LCIRC,MAXRAYS,NUMR,X,NPEAK,MODE,FFTW_PLANS,
     &      LUNEXP,LUNREF,LUNDOC)

9999    CALL FFTW3_KILLPLANS(FFTW_PLANS,NPLANS,IRTFLG)

        CLOSE(LUNREF)
        CLOSE(LUNEXP)  
        CLOSE(LUNDOC)
  
        IF (ALLOCATED(NUMR)) DEALLOCATE (NUMR)
        IF (ALLOCATED(X))    DEALLOCATE (X)

        END
