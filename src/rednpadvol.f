

C++*********************************************************************
C
C  REDNPADVOL.F   NEW                           JUL 2008 ARDEAN LEITH
C                 REDNPADVOL_1P                 AUG 2011 ARDEAN LEITH
C                 REDNPADVOL_SEL                NOV 2011 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C REDNPADVOL(LUNVOL,PADVAL,NXIN,NYIN,NZIN,
C                         NXOUT,NYOUT,NZOUT, BUFVOL, IRTFLG)
C
C PURPOSE: LOAD AND PAD A VOLUME 
C
C PARAMETERS:
C       LUNVOL                       I/O UNITS                 (INPUT)
C       PADVAL                       PADDING VALUE             (INPUT)
C       NXIN,NYIN,NZIN               INPUT VOLUME DIM.         (INPUT)
C       NXOUT,NYOUT,NZOUT            VOLUME SIZE               (INPUT)
C       BUFVOL                       VOLUME                    (OUTPUT)
C       IRTFLG                       ERROR FLAG                (OUTPUT)
C
C **********************************************************************

        SUBROUTINE REDNPADVOL(LUNVOL,PADVAL, 
     &                       NXIN,NYIN,NZIN,
     &                       NXOUT,NYOUT,NZOUT,
     &                       BUFVOL, IRTFLG)

        IMPLICIT NONE

	INTEGER, INTENT(IN)   :: LUNVOL 
	REAL,    INTENT(IN)   :: PADVAL
	INTEGER, INTENT(IN)   :: NXIN,NYIN,NZIN 
	INTEGER, INTENT(IN)   :: NXOUT,NYOUT,NZOUT
 	REAL,    INTENT(OUT)  :: BUFVOL(NXOUT,NYOUT,NZOUT)
	INTEGER, INTENT(OUT)  :: IRTFLG

	INTEGER               :: IREC,IZ,IY
        LOGICAL, PARAMETER    :: MPIBCAST = .FALSE.

        IREC = 1

        DO IZ = 1,NZIN
           DO IY = 1,NYIN
              CALL REDLIN_SEL(LUNVOL,NXIN,IREC,MPIBCAST,
     &                        BUFVOL(1,IY,IZ),IRTFLG)

C             FILL REMAINING PADDING COLS 
              BUFVOL(NXIN+1:NXOUT, IY, IZ) = PADVAL

              IREC = IREC + 1
           ENDDO

C          FILL REMAINING PADDING ROWS 
           BUFVOL(1:NXOUT, NYIN+1:NYOUT, IZ) = PADVAL
        ENDDO

C       FILL REMAINING PADDING SLICES 
        BUFVOL(1:NXOUT, 1:NYOUT, NZIN+1:NZOUT) = PADVAL

        END

C **************************** REDNPADVOL_SEL ***********************

        SUBROUTINE REDNPADVOL_SEL(LUNVOL,PADVAL, 
     &                       NXIN,NYIN,NZIN,
     &                       NXOUT,NYOUT,NZOUT, MPIBCAST,
     &                       BUFVOL, IRTFLG)

        IMPLICIT NONE

	INTEGER, INTENT(IN)   :: LUNVOL 
	REAL,    INTENT(IN)   :: PADVAL
	INTEGER, INTENT(IN)   :: NXIN,NYIN,NZIN 
	INTEGER, INTENT(IN)   :: NXOUT,NYOUT,NZOUT
        LOGICAL, INTENT(IN)   :: MPIBCAST
 	REAL,    INTENT(OUT)  :: BUFVOL(NXOUT,NYOUT,NZOUT)
	INTEGER, INTENT(OUT)  :: IRTFLG

	INTEGER               :: IREC,IZ,IY

        IREC = 1

        DO IZ = 1,NZIN
           DO IY = 1,NYIN
              CALL REDLIN_SEL(LUNVOL,NXIN,IREC,MPIBCAST,
     &                        BUFVOL(1,IY,IZ),IRTFLG)

C             FILL REMAINING PADDING COLS 
              BUFVOL(NXIN+1:NXOUT, IY, IZ) = PADVAL

              IREC = IREC + 1
           ENDDO

C          FILL REMAINING PADDING ROWS 
           BUFVOL(1:NXOUT, NYIN+1:NYOUT, IZ) = PADVAL
        ENDDO

C       FILL REMAINING PADDING SLICES 
        BUFVOL(1:NXOUT, 1:NYOUT, NZIN+1:NZOUT) = PADVAL

        END


C **********************************************************************
C
C REDMASKNPADVOL(LUNVOL,PADVAL,RADI,NXIN,NYIN,NZIN,
C               NXOUT,NYOUT,NZOUT, BUFVOL,MPIBCAST, IRTFLG)
C
C PURPOSE: LOAD, MASK, AND PAD A VOLUME 
C
C PARAMETERS:
C       LUNVOL                       I/O UNITS                 (INPUT)
C       PADVAL                       PADDING VALUE             (INPUT)
C       RADI                         MASK RADIUS               (INPUT)
C       NXIN,NYIN,NZIN               INPUT VOLUME DIM.         (INPUT)
C       NXOUT,NYOUT,NZOUT            VOLUME SIZE               (INPUT)
C       MPIBCAST                     MPI BROADCAST BUFVOL      (INPUT)
C       BUFVOL                       VOLUME                    (OUTPUT)
C       IRTFLG                       ERROR FLAG                (OUTPUT)
C
C **********************************************************************

        SUBROUTINE REDMASKNPADVOL(LUNVOL,PADVAL,RADI,
     &                       NXIN, NYIN, NZIN,
     &                       NXOUT,NYOUT,NZOUT,
     &                       MPIBCAST, BUFVOL,IRTFLG)

        IMPLICIT NONE

	INTEGER, INTENT(IN)  :: LUNVOL 
	REAL,    INTENT(IN)  :: PADVAL
	REAL,    INTENT(IN)  :: RADI
	INTEGER, INTENT(IN)  :: NXIN, NYIN, NZIN 
	INTEGER, INTENT(IN)  :: NXOUT,NYOUT,NZOUT
        LOGICAL, INTENT(IN)  :: MPIBCAST 
 	REAL,    INTENT(OUT) :: BUFVOL(NXOUT,NYOUT,NZOUT)
	INTEGER, INTENT(OUT) :: IRTFLG

	INTEGER              :: IREC,IZ,IY,IX
        REAL                 :: RADISQ,RSQ,XCEN,YCEN,ZCEN

        RADISQ = RADI **2

        XCEN   = NXIN/2 + 1
        YCEN   = NYIN/2 + 1
        ZCEN   = NZIN/2 + 1

        IREC = 1

        DO IZ = 1,NZIN
           DO IY = 1,NYIN

              CALL REDLIN_SEL(LUNVOL,NXIN,IREC,MPIBCAST,
     &                        BUFVOL(1,IY,IZ),IRTFLG)

C             MASK OUTSIDE CIRCLE
              IF (NZIN == 1) THEN    ! FASTER FOR IMAGES
                 DO IX = 1,NXIN
                    RSQ = (FLOAT(IX) - XCEN) **2 + 
     &                    (FLOAT(IY) - YCEN) **2  
                    IF (RSQ > RADISQ) BUFVOL(IX,IY,1) = PADVAL 
                 ENDDO
              ELSE
                 DO IX = 1,NXIN
                    RSQ = (FLOAT(IX) - XCEN) **2 + 
     &                    (FLOAT(IY) - YCEN) **2 + 
     &                    (FLOAT(IZ) - ZCEN) **2
                    IF (RSQ > RADISQ) BUFVOL(IX,IY,IZ) = PADVAL 
                 ENDDO
              ENDIF

C             FILL REMAINING PADDING COLS 
              BUFVOL(NXIN+1:NXOUT, IY, IZ) = PADVAL

              IREC = IREC + 1
           ENDDO

C          FILL REMAINING PADDING ROWS 
           BUFVOL(1:NXOUT, NYIN+1:NYOUT, IZ) = PADVAL
        ENDDO

C       FILL REMAINING PADDING SLICES 
        BUFVOL(1:NXOUT, 1:NYOUT, NZIN+1:NZOUT) = PADVAL

        END


