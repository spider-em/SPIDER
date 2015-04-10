
C++*********************************************************************
C
C DUST.F
C                 SETPRMB PARAMETERS       ARDEAN LEITH         5/19/09
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
C    HIST(LUN,NX,NY,NZ,HSIG,HMODE,IRTFLG)
C
C    PURPOSE:    REMOVE DATA THAT ARE 
C                OUT OF A SPECIFIED STATISTICAL RANGE.
C
C    PARAMETERS:  LUN        IO UNIT NUMBER OF IMAGE FILE
C                 NX,NY  DIMENSIONS OF IMAGE
C                 NZ     DIMENSIONS OF IMAGE
C                 HSIG       HISTOGRAM STANDARD DEVIATION
C                 HMODE      HISTOGRAM MODE
C                 IRTFLG     UNUSED
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE DUST(LUN,NX,NY,NZ,HSIG,HMODE,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      REAL, ALLOCATABLE :: REDBUF(:)

      FMULT = 2.0
      CALL RDPRM1S(FMULT,NOT_USED,'STANDARD DEVIATION FACTOR',IRTFLG)
      IF (IRTFLG .EQ. -1) RETURN

      ISIDE = 3
      CALL RDPRI1S(ISIDE,NOT_USED,
     &             'BOTTOM=(1),  TOP=(2), OR BOTH SIDES=(3)',IRTFLG)
      IF (IRTFLG .EQ. -1) RETURN

      IF (FCHAR(4:4) == 'V') THEN
         CALL RDPRM1S(VALUE,NOT_USED,'VALUE TO BE SUBSTITUTED',IRTFLG)
         IF (IRTFLG .EQ. -1) RETURN
      ENDIF

      BOT  = FMIN
      TOP  = FMAX
      YSUB = HSIG * FMULT

      IF (ISIDE == 1 .OR. ISIDE == 3) BOT = HMODE - YSUB
      IF (ISIDE == 2 .OR. ISIDE == 3) TOP = HMODE + YSUB

C     FOR 'DU' CASE
      TT  = TOP
      TB  = BOT

C     FOR 'DU V' CASE
      IF (FCHAR(4:4) == 'V') THEN
         TT = VALUE
         TB = VALUE
      ENDIF

      IF (ISIDE == 1 .OR. ISIDE == 3) THEN
         WRITE(NOUT,90) '  REPLACING PIXELS < ',BOT,'  WITH: ',TB
      ENDIF

      IF (ISIDE == 2 .OR. ISIDE == 3) THEN
         WRITE(NOUT,90) '  REPLACING PIXELS > ',TOP,'  WITH: ',TT
 90      FORMAT(A,ES10.3, A,ES10.3)
      ENDIF

      ALLOCATE (REDBUF(NX), STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN 
         CALL ERRT(46,'DUST , REDBUF',NX)
         RETURN
      ENDIF

      DO  I=1,NY*NZ
         CALL REDLIN(LUN,REDBUF,NX,I)

         DO  K=1,NX
            IF (REDBUF(K) > TOP) REDBUF(K) = TT
            IF (REDBUF(K) < BOT) REDBUF(K) = TB
         ENDDO
         CALL WRTLIN(LUN,REDBUF,NX,I)
      ENDDO

C     SET STATISTICS TO UNDETERMINED
      CALL SETPRMB(LUN, 0.0,0.0, 0.0,0.0)

      DEALLOCATE (REDBUF)

      END

