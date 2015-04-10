 
 
#ifdef SP_SUN4 
 
C   THIS ROUTINE NOT AVAILABLE AT SUN SITES
 
       SUBROUTINE vax32u
 
       COMMON /UNITS/LUNC,NIN,NOUT
 
       WRITE(NOUT,*) 'DUMMY CALL: vax32u'
       RETURN
       END
#else
 
C++*********************************************************************
C
C  VAX32U.FOR  --  JUL 93 al
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
C    VAX32U(FVAL)
C
C    PURPOSE:  QUICKLY CONVERT A VAX FLOATING POINT NUMBER TO A UNIX 
C              FLOATING POINT NUMBER.  THIS SUBROUTINE IS MEANT TO BE
C              RUN ON A UNIX SYSTEM (NOT ON A VAX!)
C
C    PARAMETERS:   IVAL       FLOATING POINT NUMBER!!!
C
C    NOTES:        IVALIN     --  --  --  --
C                  LVALIN      1   2   3   4
C                             IN1 IN2 IN3 IN4
C
C
C                  IVALOUT    --  --  --  --
C                  LVALOUT     1   2   3   4
C                  I2         -----1
C
C **********************************************************************

      SUBROUTINE VAX32U(IVAL)

      LOGICAL * 1  IO1,IO2,IO3,IO4,IN1,IN2,IN3,IN4,LVALIN(4),LVALOUT(4)

      INTEGER * 2  I2

      INTEGER * 4  IVALIN,IVAL,IVALOUT
      EQUIVALENCE                 (LVALIN,IVALIN),(LVALOUT,IVALOUT)
      
      EQUIVALENCE  (IN1,LVALIN(1)), (IN2,LVALIN(2)), (IN3,LVALIN(3)),
     &             (IN4,LVALIN(4))
      EQUIVALENCE  (IO1,LVALOUT(1)),(IO2,LVALOUT(2)),(IO3,LVALOUT(3)),
     &             (IO4,LVALOUT(4))

#ifdef __osf__
      EQUIVALENCE  (I2,LVALOUT(3))
#else
      EQUIVALENCE  (I2,IVALOUT)
#endif
    
      INTEGER * 2  P32640
      P32640 = 32640

      IVALIN = IVAL

C     CHANGE THE EXPONENT TO EXCESS 127 NOTATION FOR UNIX

C     FLIP THE BYTES IN EACH OF THE WORDS
#ifndef __osf__
      IO1 = IN2
      IO2 = IN1
      IO3 = IN4
      IO4 = IN3
#else
      IO1 = IN3
      IO2 = IN4
      IO3 = IN1
      IO4 = IN2
#endif

C     32640 = 0111 1111 1000 0000
C     256 = 0000 0001 0000 0000

      IF (IAND(I2,P32640) .NE. 0)  I2 = I2 - 256

      IVAL =  IVALOUT

      RETURN
      END
   
#endif
