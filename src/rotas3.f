C++*********************************************************************
C
C ROTAS3.F                        REWRITTEN APRIL 2002 ArDean Leith
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
C IMAGE_PROCESSING_ROUTINE
C
C  ROTAS3(LUN1,LUN2,NSAM,NROW,NSLICE,MODE)
C
C  PARAMETERS:     MODE                                        (SENT)
C                  3D   - QUADRATIC INTERPOLATION
C                  3A   - QUADRATIC WITH SPECIFIED ORIGIN
C                  3DL  - LINEAR INTERPOLATION
C                  3AL  - LINEAR WITH SPECIFIED ORIGIN
C
C  PURPOSE:        3D ROTATION USING EULER ANGLES OF VOLUME. 
C                  CHOOSES LINEAR OR TRI-QUADRATIC INTERPOLATION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE ROTAS3(LUN1,LUN2,NSAM,NROW,NSLICE,MODE) 
    
        INCLUDE 'CMBLOCK.INC'
  
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: AIMG
        CHARACTER(LEN=3)                    ::  MODE

C       3-D  ROTATION

        ALLOCATE (AIMG(NSAM,NROW,NSLICE), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'RT 3, AIMG',IER)
           RETURN
        ENDIF

        CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,AIMG,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       GET EULER ANGLES
        PHI   = 0.0
        THETA = 0.0
        PSI   = HUGE(PSI)

        CALL  RDPRM3S(PHI,THETA,PSI,NOT_USED,'PHI, THETA, & PSI',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (PSI .EQ. HUGE(PSI)) THEN
           PSI = 0.0
           CALL  RDPRM1S(PSI,NOT_USED,'PSI',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF


        IF (MODE(2:2) .EQ. 'A')  THEN
C          USER SPECIFIES CENTER OF ROTATION
           NZ = HUGE(NZ)
           CALL  RDPRI3S(NX,NY,NZ,NOT_USED,
     &                 'X, Y & Z FOR CENTER OF ROTATION',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (NZ .EQ. HUGE(NZ)) THEN
              CALL  RDPRI1S(NZ,NOT_USED,'Z FOR CENTER OF ROTATION',
     &                      IRTFLG)
              IF (IRTFLG .NE. 0) RETURN
           ENDIF

           KLX = -NX+1
           KNX = NSAM-NX
           KLY = -NY+1
           KNY = NROW-NY
           KLZ = -NZ+1
           KNZ = NSLICE-NZ
 
        ELSE
           IF (MOD(NSAM,2) .EQ. 0)  THEN
              KNX = NSAM/2-1
           ELSE
              KNX = NSAM/2
           ENDIF
           KLX = -NSAM/2
           IF (MOD(NROW,2) .EQ. 0)  THEN
              KNY = NROW/2-1
           ELSE
              KNY = NROW/2
           ENDIF
           KLY = -NROW/2
           IF (MOD(NSLICE,2) .EQ. 0)  THEN
              KNZ = NSLICE/2-1
           ELSE
              KNZ = NSLICE/2
           ENDIF
           KLZ = -NSLICE/2
        ENDIF

        IF (THETA.EQ.0.0 .AND .PHI.EQ.0.0 .AND. PSI.EQ.0.0)  THEN
C          NO ROTATION NEEDED
           CALL WRTVOL(LUN2,NSAM,NROW,1,NSLICE,AIMG,IRTFLG)

        ELSEIF (MODE(3:3) .EQ. 'S') THEN
C          LINEAR INTERLPOLATION
           CALL ROTS3(LUN2,AIMG,KLX,KNX,KLY,KNY,KLZ,KNZ,PSI,THETA,PHI)

        ELSE
C          TRI-QUADRATIC INTERLPOLATION (DEFAULTED IN MAY 2002)
           CALL ROTS3Q(LUN2,AIMG,KLX,KNX,KLY,KNY,KLZ,KNZ,PSI,THETA,PHI)
        ENDIF

        IF (ALLOCATED(AIMG)) DEALLOCATE(AIMG)
        
        RETURN
        END


