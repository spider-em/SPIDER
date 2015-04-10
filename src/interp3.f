C ++********************************************************************
C                                                                      *
C INTERP3     ADDE ALLOCATE                        MAY 00 ARDEAN LEITH *
C             USED ALLOCATE ALWAYS                 OCT 11 ARDEAN LEITH *
C                                                                      *
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
C     INTERP3(LUN1,LUN2,NSAM,NROW,NSLICE,NSAM1,NROW1,NSLICE1,IDUM)
C
C     PURPOSE:  3 D INTERPOLATION
C                                                                     *
C **********************************************************************

      SUBROUTINE INTERP3(LUN1,LUN2,NSAM,NROW,NSLICE,
     &                             NSAM1,NROW1,NSLICE1,IDUM)
 
      INCLUDE 'CMBLOCK.INC' 
      INCLUDE 'CMLIMIT.INC'
 
      COMMON /IOBUF/ BUF(NBUFSIZ)

      REAL, ALLOCATABLE, DIMENSION(:) :: QPT
 
C     ALLOCATE MEMORY FOR INPUT VOLUME
      MEMWANT = NSAM * NROW * NSLICE

      ALLOCATE (QPT(MEMWANT), STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CALL ERRT(46,'INTERP3, QPT',NDUM)
         RETURN
      ENDIF
      WRITE(NOUT,*) 'ALLOCATED MEMORY FOR VOLUME:',MEMWANT

C     LOAD VOLUME
      CALL READV(LUN1,QPT,NSAM,NROW,NSAM,NROW,NSLICE)

C     CALL INTERPOLATION ROUTINE
      CALL IRP3(QPT,BUF,NSAM,NROW,NSLICE,
     &          NSAM1,NROW1,NSLICE1,LUN2)
 
      DEALLOCATE(QPT)

      RETURN
      END
