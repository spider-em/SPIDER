C ++********************************************************************
C                                                                      *
C DECIMATE.F                                                           *
C                  OPFILEC                       ARDEAN LEITH   FEB 03 *
C                  CHECKED RANGES, REWRITE       ARDEAN LEITH   FEB 09 *
C                  SETPRS                        ARDEAN LEITH   NOV 10 *
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
C  DECIMATE
C
C  PURPOSE: DECIMATE 2-D OR 3-D REAL IMAGE  BY SKIPPING PIXELS  OR BY
C           AVERAGING OF PIXELS               
C                                                                 
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE DECIMATE

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)      :: FILNAM
	CHARACTER(LEN=1)           :: NULL = CHAR(0)
	REAL, ALLOCATABLE          :: Q(:),W(:)

	INTEGER, PARAMETER         :: LUN1 = 50
	INTEGER, PARAMETER         :: LUN2 = 51

	MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	CALL FILERD(FILNAM,NLETI,NULL,'OUTPUT',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

	MX = 0
	MY = 0
        MZ = 0
        IF (IFORM .EQ. 1) THEN
           CALL RDPRIS(MX,MY,NOT_USED,
     &                'DECIMATION FACTORS FOR X & Y',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999
	   MZ = 1
        ELSE
           CALL RDPRI3S(MX,MY,MZ, NOT_USED,
     &                 'DECIMATION FACTORS FOR X, Y, & Z',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 999

	   IF (MZ .LE. 0 ) THEN
              CALL RDPRI1S(MZ,NOT_USED,'DECIMATION FACTOR FOR Z',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 999
          ENDIF
	ENDIF
        IF (IRTFLG .NE. 0) GOTO 999

	IF (MY .EQ.0) MY = MX
	IF (MZ .EQ.0) MZ = MX

	IF (MX <= 0) THEN
           CALL  ERRT(102,'ILLEGAL X DECIMATION FACTOR',MX)
           GOTO 999
	ELSEIF (MY <= 0) THEN
           CALL  ERRT(102,'ILLEGAL Y DECIMATION FACTOR',MY)
           GOTO 999
	ELSEIF (MZ <= 0) THEN
           CALL  ERRT(102,'ILLEGAL Z DECIMATION FACTOR',MZ)
           GOTO 999
	ENDIF

	NSAMP   = NSAM   / MX
	NROWP   = NROW   / MY
	NSLICEP = NSLICE / MZ

	MAXIM = 0
        CALL OPFILEC(LUN1,.FALSE.,FILNAM,LUN2,'U',IFORM,
     &               NSAMP,NROWP,NSLICEP,
     &               MAXIM,NULL,.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

	ALLOCATE(Q(NSAM),STAT=IRTFLG)
	IF (IRTFLG .NE. 0) THEN
           CALL  ERRT(46,'DC, Q',NSAM)
           GOTO 999
        ENDIF

	IF (FCHAR(4:4) .EQ. 'S')  THEN
C          DECIMATE BY AVERAGING PIXELS
	   ALLOCATE(W(NSAMP), STAT=IRTFLG)
	   IF (IRTFLG .NE. 0) THEN
               CALL  ERRT(46,'DC, W',NSAMP)
               GOTO 998
           ENDIF

	   CALL DECIMS(LUN1,LUN2,Q,W,
     &		 NSAM,NROW,NSLICE, NSAMP,NROWP, MX,MY,MZ)
	   IF (ALLOCATED(W)) DEALLOCATE(W)

	ELSE
C          SIMPLE DECIMATION
	   CALL DECIM(LUN1,LUN2,Q,
     &		NSAM,NROW,NSLICE,NSAMP,NROWP, MX,MY,MZ)
	ENDIF

C       FILE HEADER FOR PIXSIZ HAS CHANGED
        SCALEX = FLOAT(NSAM)   / FLOAT(NSAMP)
        SCALEY = FLOAT(NROW)   / FLOAT(NROWP)
        SCALEZ = FLOAT(NSLICE) / FLOAT(NSLICEP) 
        !write(6,*) ' set scaless:',scalex,scaley,scalez

        SCALET  = SCALEX
        IF (SCALEY .NE. SCALEX) SCALET = 0.0   ! X NOT SAME AS Y
        IF (NSLICE > 1 .AND. SCALEZ .NE. SCALEX) SCALET = 0.0

C       UPDATE THE INCORE HEADER VALUE & FILE HEADER FOR PIXSIZ
        CALL SETPRMS(LUN2, SCALET,IRTFLG)

        !write(6,*) ' set scale:',scalet

998     IF (ALLOCATED(Q)) DEALLOCATE(Q)
999     CLOSE(LUN1)
	CLOSE(LUN2)

	END

C       --------------------------- DECIMS ----------------------------

	SUBROUTINE DECIMS(LUN1,LUN2,BI,BO,
     &		    NSAM,NROW,NSLICE, NSAMP,NROWP, MX,MY,MZ)

	REAL :: BI(NSAM) ,BO(NSAMP)
        
        FACM = 1.0 / REAL(MX * MY * MZ)

        DO K=1,NSLICE-MZ+1,MZ
           DO J=1,NROW-MY+1,MY 

              BO = 0.0             ! ZEROS WHOLE ARRAY

              DO KT=1,MZ
                 DO JT=1,MY
                    IRECIN = J + JT - 1 + NROW * (K+KT-2)
                    CALL REDLIN(LUN1,BI,NSAM,IRECIN)

                    DO I=1,NSAM-MX+1,MX
                       DO IT=0,MX-1
                          ILOC     = (I+MX-1)/MX
                          BO(ILOC) = BO(ILOC) +  BI(I+IT)
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO

              BO = BO * FACM      ! ARRAY OPERATION

              IRECOUT = (J+MY-1)/MY + NROWP * ((K+MZ-1)/MZ-1)
              CALL WRTLIN(LUN2,BO,NSAMP,IRECOUT)

           ENDDO
        ENDDO

	END	

C       --------------------------- DECIM ----------------------------

	SUBROUTINE  DECIM(LUN1,LUN2,BI,
     &		NSAM,NROW,NSLICE,NSAMP,NROWP,MX,MY,MZ)

	REAL ::  BI(NSAM)

        DO K=1,NSLICE-MZ+1,MZ
           DO J=1,NROW-MY+1,MY
              CALL  REDLIN(LUN1,BI,NSAM,J+NROW*(K-1))

              DO I=1,NSAM-MX+1,MX
                 BI((I+MX-1)/MX) = BI(I)
              ENDDO

              IRECOUT = (J+MY-1)/MY + NROWP * ((K+MZ-1)/MZ-1)
              CALL WRTLIN(LUN2,BI,NSAMP,IRECOUT)
           ENDDO
        ENDDO

        END	
