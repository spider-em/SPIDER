
C++*********************************************************************
C
C TRIQUAD.F                            NEW APRIL 4 2002 ArDean Leith
C                                 F(1) RST BUG JUL 2004 ArDean Leith
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
C  TRIQUAD(R,S,T,F)
C 
C  PURPOSE:        TRI-QUADRATIC INTERPOLATION
C
C  NOTE:        PROBABLY CREATED BY GENERALIZATION FROM QUADRI.F
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      REAL FUNCTION TRIQUAD(R,S,T,F)

      DOUBLE PRECISION               :: R,S,T
      REAL, DIMENSION(27),INTENT(IN) :: F

      REAL, PARAMETER :: C2 = 1.0 / 2.0
      REAL, PARAMETER :: C4 = 1.0 / 4.0
      REAL, PARAMETER :: C8 = 1.0 / 8.0

      RS   = R * S
      ST   = S * T
      RT   = R * T
      RST  = R * ST

      RSQ  = 1-R**2
      SSQ  = 1-S**2
      TSQ  = 1-T**2

      RM1  = (1-R)
      SM1  = (1-S)
      TM1  = (1-T)

      RP1  = (1+R)
      SP1  = (1+S)
      TP1  = (1+T)

      TRIQUAD = 
     &  (-C8) * RST * RM1  * SM1  * TM1 * F( 1) + 
     &  ( C4) * ST  * RSQ  * SM1  * TM1 * F( 2) + 
     &  ( C8) * RST * RP1  * SM1  * TM1 * F( 3) + 
     &  ( C4) * RT  * RM1  * SSQ  * TM1 * F( 4) + 
     &  (-C2) * T   * RSQ  * SSQ  * TM1 * F( 5) + 
     &  (-C4) * RT  * RP1  * SSQ  * TM1 * F( 6) + 
     &  ( C8) * RST * RM1  * SP1  * TM1 * F( 7) + 
     &  (-C4) * ST  * RSQ  * SP1  * TM1 * F( 8) + 
     &  (-C8) * RST * RP1  * SP1  * TM1 * F( 9) + 

     &  ( C4) * RS  * RM1  * SM1  * TSQ * F(10) + 
     &  (-C2) * S   * RSQ  * SM1  * TSQ * F(11) + 
     &  (-C4) * RS  * RP1  * SM1  * TSQ * F(12) + 
     &  (-C2) * R   * RM1  * SSQ  * TSQ * F(13) + 
     &                RSQ  * SSQ  * TSQ * F(14) + 
     &  ( C2) * R   * RP1  * SSQ  * TSQ * F(15) + 
     &  (-C4) * RS  * RM1  * SP1  * TSQ * F(16) + 
     &  ( C2) * S   * RSQ  * SP1  * TSQ * F(17) + 
     &  ( C4) * RS  * RP1  * SP1  * TSQ * F(18) +
 
     &  ( C8) * RST * RM1  * SM1  * TP1 * F(19) + 
     &  (-C4) * ST  * RSQ  * SM1  * TP1 * F(20) + 
     &  (-C8) * RST * RP1  * SM1  * TP1 * F(21) + 
     &  (-C4) * RT  * RM1  * SSQ  * TP1 * F(22) + 
     &  ( C2) * T   * RSQ  * SSQ  * TP1 * F(23) + 
     &  ( C4) * RT  * RP1  * SSQ  * TP1 * F(24) + 
     &  (-C8) * RST * RM1  * SP1  * TP1 * F(25) + 
     &  ( C4) * ST  * RSQ  * SP1  * TP1 * F(26) + 
     &  ( C8) * RST * RP1  * SP1  * TP1 * F(27)  

         END
