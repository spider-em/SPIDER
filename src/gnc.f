C ++********************************************************************
C                                                                      *
C  GNC                                                                 *
C                                                                      *
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
C                                                                      *
C  GNC                                                                 *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C
C IMAGE_PROCESSING_ROUTINE
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE GNC(LUN1,LUN2,NSAM,NROW)

        INCLUDE 'CMBLOCK.INC'
      
        REAL, ALLOCATABLE, DIMENSION(:,:) :: AIMG
        REAL, ALLOCATABLE, DIMENSION(:,:) :: BIMG

	
        CALL RDPRMI  (L,IDUM, NOT_USED,'LAMBDA')
        QL=L

        CALL RDPRM  (H0, NOT_USED,'H0')

        CALL RDPRM  (EPS, NOT_USED,'EPS')
 
        ALLOCATE (AIMG(NSAM,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'GNC, AIMG',IER)
           RETURN
        ENDIF

        ALLOCATE (BIMG(NSAM,NROW), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'GNC, BIMG',IER)
           RETURN
        ENDIF

        NSLICE = 1
        CALL READV(LUN1,AIMG,NSAM,NROW,NSAM,NROW,NSLICE)

        CALL  GNC2S(AIMG,BIMG,NSAM,NROW,QL,H0,EPS,NDAT)

        CALL WRITEV(LUN2,BIMG,NSAM,NROW,NSAM,NROW,NSLICE)

        DEALLOCATE (AIMG)
        DEALLOCATE (BIMG)
        END 
