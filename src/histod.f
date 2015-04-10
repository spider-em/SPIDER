C++*********************************************************************
C
C HISTOD.F
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
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C--*********************************************************************

        SUBROUTINE HISTOD

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
 
        CHARACTER(LEN=MAXNAM)   ::  FILNAM
        REAL, ALLOCATABLE, DIMENSION(:,:,:) :: AIMG
	        
        CHARACTER *1  NULL

        DATA  LUN1/11/

       
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',ITYPE,
     &          NSAM,NROW,NSLICE,MAXIM,
     &          'IMAGE TO CORRECT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO  130
  
111     IF(IFORM.LT.1) GOTO 145
        NSR=NSAM*NROW*NSLICE
        CALL  RDPRMI(LENH,ITRMAX,NOT_USED,'HISTOGRAM LENGTH')
	
        ALLOCATE (AIMG(NSAM,NROW,NSLICE), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'CE OD, AIMG',IER)
           CLOSE(LUN1)
           RETURN
        ENDIF

       
        CALL READV(LUN1,AIMG,NSAM,NROW,NSAM,NROW,NSLICE)
        CLOSE(LUN1)

        CALL  HISTODC(AIMG,NSR,LENH)

        DEALLOCATE(AIMG)
        RETURN

130     CALL ERRT(4,'CE OD  ',NE)
        GOTO 5
145     CALL ERRT(2,'CE OD  ',NE)
5       CLOSE(LUN1)
        END
