
C++*********************************************************************
C
C  REDREC8.F  -- CREATED OCT 98 ArDean Leith
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
C  REDREC8(LUN,LBUF,NB,NREC,IRTFLG)
C
C  PURPOSE:  BASIC RAW INPUT ROUTINE
C
C        0         2         3         4         5         6         7     
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE REDREC8(LUN,LBUF,LENREC,NREC,IRTFLG)

      INCLUDE 'CMBLOCK.INC'


      LOGICAL * 1   LBUF(LENREC)

      READ(LUN,REC=NREC,IOSTAT=IERR) LBUF

      IF (IERR .NE. 0) THEN
         WRITE(NOUT,*) '*** ERROR:',IERR, ' READING INPUT RECORD:',NREC
         CALL ERRT(100,'REDREC8',NE)
         IRTFLG = 1
         RETURN
      ENDIF

      IRTFLG = 0
      RETURN
      END

