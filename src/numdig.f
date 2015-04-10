
C++*********************************************************************
C
C    NUMDIG.F   -- CREATED OCT 88 ArDean Leith
C                  ADDED SUPPORT FOR NEGATIVES AUG 2002 ARDEAN LEITH
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
C    NUMDIG(IVALUE,MINVAL)
C
C    PURPOSE:  FIND NUMBER OF DIGITS IN IVALUE
C
C    PARAMETERS:  
C         IVALUE    INPUT NUMBER
C         MINVAL    MIN. VALUE RETURNED FOR IVALUE
C
C **********************************************************************

          FUNCTION NUMDIG(IVALUE,MINVAL)

          IVALUET = ABS(IVALUE)

          LENI  = 1
          LENT  = 9
 3        IF (IVALUET .GT. LENT) THEN
             LENT = 9 * (10**LENI) + LENT
             LENI = LENI + 1
             GOTO 3
          ENDIF

C         WHAT IF < 0
          IF (IVALUE .LT. 0) LENI = LENI + 1

C         MINVAL IS MINIMUM LENGTH
          NUMDIG = MAX(MINVAL,LENI)         

          RETURN
          END

       
