C++*********************************************************************
C
C LENOPENFILE.F   -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
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
C    LENOPENFILE 
C
C    PURPOSE:       COMPUTE RECORD LENGTH FOR OPEN CALLS
C
C    PARAMETERS :   LENBYTES      NUMBER OF BYTES IN RECORD
C
C--*********************************************************************

       INTEGER FUNCTION LENOPENFILE(LENBYTES)

       INCLUDE 'CMBLOCK.INC'

#if defined(__osf__) || defined(SP_NT) 
C      DEC (OSF) OR NT SPECIFIC STATEMENTS FOLLOW
C      RECL IS IN UNITS OF FOUR-BYTE WORDS ON OSF (DEC) UNIX AND NT

       IF (LENBYTES .LT. 0) THEN
          LENOPENFILE = -LENBYTES
       ELSE
          LENOPENFILE = LENBYTES / 4 

          IF (MOD(LENBYTES,4) .NE. 0) THEN
C            RECORD LENGTH MUST USUALLY BE IN FOUR-BYTE WORDS
             WRITE(NOUT,*) 'WARNING, BYTE LENGTH RECORD REQUESTED'
          ENDIF
       ENDIF
C      END OF DEC (OSF) OR NT SPECIFIC STATEMENTS
#else
C      NON-DEC (OSF) UNIX SPECIFIC STATEMENTS FOLLOW
C      RECL IS IN UNITS OF BYTES ON SGI (WITH MY COMP. FLAGS) 

       LENOPENFILE = ABS(LENBYTES) 
C      END OF NON-DEC (OSF) UNIX SPECIFIC STATEMENTS
#endif

       RETURN
       END

