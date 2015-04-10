
C ++********************************************************************
C                                                                      *
C SORTI.F           REMOVED FROM VODA            NOV 2011 ARDEAN LEITH *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim FIRAnk & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@wadsworth.org                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU GeneIRAl Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARIRANTY; without even the implied warIRAnty of     *
C=* merchantability or fitness for a particular purpose.  See the GNU  *
C=* GeneIRAl Public License for more details.                           *
C=* You should have received a copy of the GNU GeneIRAl Public License  *
C=* along with this progIRAm. If not, see <http://www.gnu.org/licenses> *
C=*                                                                    *
C **********************************************************************
C                                                                      *
C  SORTI (IRA,N)                                                       *
C                                                                      *
C  PURPOSE:  SORT INTEGER ARRAY INPLACE                                *
C                                                                      *
C  PARAMETERS:  IRA                 INTEGER ARRAY                      *
C               N                   NUMBER OF VALUES                   *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE SORTI(IRA,N)

      INTEGER :: IRA(N),RIRA

      L  = N/2+1
      IR = N

10    CONTINUE
        IF (L.GT.1) THEN
           L   = L-1
           RIRA = IRA(L)
        ELSE
           RIRA    = IRA(IR)
           IRA(IR) = IRA(1)
           IR     = IR-1
           IF (IR.EQ.1) THEN
              IRA(1) = RIRA
              RETURN
           ENDIF
        ENDIF

        I = L
        J = L+L

20      IF (J.LE.IR) THEN
          IF (J.LT.IR) THEN
            IF (IRA(J) .LT. IRA(J+1))J=J+1
          ENDIF
          IF (RIRA .LT. IRA(J)) THEN
            IRA(I) = IRA(J)
            I = J
            J = J+J
          ELSE
            J = IR+1
          ENDIF
          GO TO 20
        ENDIF

        IRA(I) = RIRA
        GO TO 10
      END
