C ++********************************************************************
C                                                                      *
C   MEDIAN           K WAS WRONG MAR 01 ARDEAN LEITH                                                        *
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
C  MEDIAN(LUN1,LUN2,NSAM,NROW,NSLICE)                                                                    *
C                                                                      *
C  PURPOSE: MEDIAN FILTRATION USING BOX OR CROSS PATTERN KERNAL                                                           *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C IMAGE_PROCESSING_ROUTINE                                             *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE MEDIAN(LUN1,LUN2,NSAM,NROW,NSLICE)

        INCLUDE 'CMBLOCK.INC'
 
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: AIMG
        CHARACTER(LEN=1)                     :: MODE,NULL

        NULL = CHAR(0)

        IF (NSLICE .LE. 1) NSLICE=1

        LENGTH = 3
10      CALL RDPRI1S(LENGTH,NOT_USED,'LENGTH OF FILTER',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (LENGTH .LE. 1) THEN
           CALL ERRT(102,'LENGTH MUST BE GREATER THAN',2) 
           GOTO 10
        ELSEIF (MOD(LENGTH,2) .EQ. 0) THEN
           LENGTH = LENGTH + 1
           WRITE(NOUT,90) LENGTH 
90         FORMAT(' EFFECTIVE LENGTH OF FILTER:',I5)
        ENDIF

        CALL RDPRMC(MODE,NA,.TRUE.,'BOX OR CROSS (B/C)?',NULL,IRTFLG)
        IF (IRTFLG .EQ. -1) GOTO 10
        IF (IRTFLG .NE. 0) RETURN

        IF (MODE .EQ. 'B') THEN
           K = LENGTH * LENGTH
           IF (NSLICE .GT. 1) K = K * LENGTH
        ELSE
           K = 2 * LENGTH - 1
           IF (NSLICE .GT. 1) K = K + LENGTH - 1
        ENDIF
        
        ALLOCATE (AIMG(NSAM,NROW,NSLICE), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'CE MED, AIMG',IER)
           RETURN
        ENDIF

        CALL READV(LUN1,AIMG,NSAM,NROW,NSAM,NROW,NSLICE)

        IF (NSLICE .EQ. 1)  THEN
           CALL MD2(AIMG,NSAM,NROW,LENGTH,K,MODE,LUN2)
        ELSE
           CALL MD3(AIMG,NSAM,NROW,NSLICE,LENGTH,K,MODE,LUN2)
        ENDIF

        DEALLOCATE(AIMG)
       
        END

