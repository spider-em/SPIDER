
C ++********************************************************************
C                                                                      *
C                                                                      *
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
C      CNSTUFF.F         LAST UPDATE               8/27/96 ARDEAN LEITH
C                                                  3/27/92 ARDEAN LEITH
C                                                  2/04/98 ARDEAN LEITH
C                                      REWRITTEN   3/10/99 ARDEAN LEITH
C
C **********************************************************************
C
C    PURPOSE:      PLACE CONTOUR FROM CNTRCE IN OUTPUT FILE
C
C    PARAMETERS:   LUN     OUTPUT FILE UNIT                       (SENT)
C                  X,Y     DATA ARRAY                             (SENT)
C                  NPTS    NUMBER OF POINTS IN DATA               (SENT)
C                  MULTIZ  LOGICAL FLAG FOR STERECON INPUT FILE   (SENT)
C                  MAXPTS  THIN LIMIT                             (SENT)
C                  IPEN    COLOR FOR CONTOUR (UNUSED)             (SENT)
C                  IRTFLG  ERROR FLAG                             (RET.)
C
C    CALLED BY:    CNTRCE
C                                                          ...CNCALC
C    CALL TREE:    PLOT1...CNINT3...CNTUR...CNSCAN...CNTRCE...
C                                                          ...CNSTUFF
C                                                             SSPUSH
C
C--********************************************************************

      SUBROUTINE CNSTUFF(LUN,X,Y,NPTS,MULTIZ,MAXPTS,IPENT,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      COMMON /POLY/  MINPTS,ISLICE

      PARAMETER      (NSIZE = 2000)
      DIMENSION      DATA(3,NSIZE),WORK(3,NSIZE)
      EQUIVALENCE    (BUF(1),DATA(1,1)), (BUF(6400),WORK(1,1))

      COMMON /IOBUF/ BUF(NBUFSIZ)

      DIMENSION      X(NSIZE),Y(NSIZE)
      LOGICAL        THINED

      IRTFLG = 0

      DO I2 = 1, NPTS
         DATA(1,I2) = X(I2)
         DATA(2,I2) = Y(I2)
      ENDDO

      IF (NPTS .GT. MAXPTS) THEN 
C       TOO MANY POINTS ON THIS CONTOUR, TRY TO THIN IT
        N0      = NPTS
        RETANG  = 178.0
        MAXTRYS = 2
        CALL RATHIN(DATA,NPTS,0.0333,RETANG,MAXPTS,MAXTRYS,WORK,IRTFLGR)
      ENDIF
      THINED = (NPTS .LT. N0) 

      IF (NPTS .LT. MINPTS) THEN
C        ABANDON THIS TOO-SHORT CONTOUR
         WRITE(NOUT,1112) NPTS
 1112    FORMAT(' SHORT CONTOUR WITH: ',I4,' POINTS ABANDONED.')

      ELSEIF (FCHAR(4:4) .EQ. 'S') THEN
C        STORE THE CONTOUR IN STERECON INPUT FILE
         CALL SSPUSH(LUN,DATA,NPTS,ISLICE,IRTFLG)
      ELSE
C        STORE THE CONTOUR IN POSTSCRIPT FILE
         CALL POARAYF(LUN,DATA,NPTS,.TRUE.,.FALSE.)
      ENDIF
        
      IF (THINED) THEN
         WRITE(NOUT,101) ISLICE,N0,NPTS
  101    FORMAT(' CONTOUR:',I5,' WITH: ',I4,' --> ',I4,
     &          ' POINTS STORED.')
      ELSE
         WRITE(NOUT,103) ISLICE,NPTS
  103    FORMAT(' CONTOUR:',I5, ' WITH: ',I4,' POINTS STORED.')
      ENDIF

      RETURN
      END
