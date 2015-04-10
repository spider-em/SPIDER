C ++********************************************************************
C
C           INTENSITY ADDED                      FEB 2014 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C MCIRCL(LUN,NX,NY,RP,IDIM)
C
C PURPOSE: MAKES A CIRCLE WITHIN A FILE
C
C PARAMETERS:
C       LUN :   LOGICAL UNIT NUMBER
C       NX,NY : FILE DIMENSIONS
C       RP :    VALUE TO BE CORRECTED
C       IDIM :  1 = CIRCLE LINE
C               2 = FILLED CIRCLE 
C               IDIM < 0 : USES COOS OF 3 POINTS TO DETERMINE CIRCLE
C               IDIM > 0 : USES CENTER COOS AND RADIUS TO DETERMINE CIRCLE
C
C  NOTES: INEFFICIENT CODING al
C
C **********************************************************************

        SUBROUTINE MCIRCL(LUN,NX,NY,RP,IDIM)

        INCLUDE 'CMBLOCK.INC'

        INTEGER :: LUN,NX,NY,IDIM
        REAL    :: RP
        REAL    :: BUF(NX)

        IF (IDIM < 0) THEN
C          ENTER 3 POINTS TO DETERMINE CIRCLE

           IDIM = IABS(IDIM)
           CALL RDPRMI(IX1,IY1,NOT_USED,'COOS OF FIRST  POINT')
           CALL RDPRMI(IX2,IY2,NOT_USED,'COOS OF SECOND POINT')
           CALL RDPRMI(IX3,IY3,NOT_USED,'COOS OF THIRD  POINT')

           IF (IY1 == IY2 .AND. IY2 == IY3) GOTO 9000
           
           X1  = FLOAT(IX2-IX1)
           Y1  = FLOAT(IY2-IY1)
           X2  = FLOAT(IX3-IX1)
           Y2  = FLOAT(IY3-IY1)

           XM1 = FLOAT(IX1)+X1/2.
           YM1 = FLOAT(IY1)+Y1/2.
           XM2 = FLOAT(IX1)+X2/2.
           YM2 = FLOAT(IY1)+Y2/2.

           WRITE(NOUT,9999) XM1,YM1,XM2,YM2
9999       FORMAT(5X,'(',F5.1,',',F5.1,')',5X,'(',F5.1,',',F5.1,')')

           IF (IY1 == IY2) THEN
111           AM22 = -X2/Y2
              X    = XM1
              Y    = AM22*(X-XM2)+YM2

           ELSEIF (IY1 == IY3) THEN
112           X    = XM2
              AM11 = -X1/Y1
              Y    = AM11*(X-XM1)+YM1
           ELSE

              AM11 = -X1/Y1
              AM22 = -X2/Y2
              WRITE(NOUT,9998) AM11,AM22
9998          FORMAT(' AM11= ',F6.2,5X,' AM22= ',F6.2)

              IF (AM11 == AM22) GOTO 9000

              X = (YM2-YM1+AM11*XM1-AM22*XM2) / (AM11-AM22)
              Y = AM11 * (X-XM1)+YM1
           ENDIF

           IX = IFIX(X+.5)
           IY = IFIX(Y+.5)

           WRITE(NOUT,9997) X,Y,IX,IY
9997       FORMAT(1X,' X= ',F5.1,2X,' Y= ',F5.1,2X,' (',I2,',',I2,')')

           R = SQRT( (X-FLOAT(IX1))**2 + (Y-FLOAT(IY1))**2 )

           WRITE(NOUT,9996) R
9996       FORMAT(' RADIUS: ',F12.2)

        ELSE
           CALL RDPRM3S(X,Y,RP,NOT_USED,
     &        'CENTER COORDINATES (NX,NY) & OPTIONAL INTENSITY',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IX = NINT(X)
           IY = NINT(Y)

           CALL RDPRM1S(R,NOT_USED,'RADIUS',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (R <= 0.0) GOTO 9000
        ENDIF


9       IR     = IFIX(R+0.5)
        I0     = IY-IR
        I1     = IY+IR
        IYSTRT = MAX(1,I0)
        IYEND  = MIN(NY,I1)
        IF (IYSTRT > NY .OR. IYEND <= 0) GOTO 9000

        I0     = IX-IR
        I1     = IX+IR
        IXSTRT = MAX(1,I0)
        IXEND  = MIN(NX,I1)
        IF (IXSTRT > NX .OR. IXEND <= 0) GOTO 9000

        DO  I=IYSTRT,IYEND
           CALL REDLIN(LUN,BUF,NX,I)

           DO J=IXSTRT,IXEND
              T  = FLOAT(J-IX)**2 + FLOAT(I-IY)**2
              IT = IFIX(SQRT(T))
              IF (IDIM == 2 .AND. IT <= IR) BUF(J) = RP
              IF (IDIM == 1 .AND. IT == IR) BUF(J) = RP
           ENDDO
           CALL WRTLIN(LUN,BUF,NX,I)
        ENDDO

        RETURN

9000    CALL ERRT(101,'INCONSISTENT INPUT PARAMETERS',NF)
        
        END
