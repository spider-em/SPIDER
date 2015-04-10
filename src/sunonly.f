 

#ifndef SP_SUN4
 
C      THESE ROUTINES ONLY NEEDED ON SUN
 
       SUBROUTINE SUNONLY
 
       COMMON /UNITS/LUNC,NIN,NOUT
 
       WRITE(NOUT,*) ' THIS MESSAGE SHOULD NEVER BE SEEN!!'
       RETURN
       END
 

#else

C ++********************************************************************
C
C SUNONLY.F
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
C  THIS FILE CONTAINS INTRINSIC FUNCTIONS WHICH ARE NOT RECOGNIZED 
C  BY THE THE SUN4 FORTRAN COMPILER
C
C **********************************************************************

        SUBROUTINE DATE(STRING)

 

        CHARACTER *(*)  STRING
        CHARACTER *24   CDATE,ctime
        INTEGER         time

        ISTIME      = time()
        CDATE       = ctime(ISTIME)

        STRING(1:2) = CDATE(9:10)
        STRING(3:3) = '-'
        STRING(4:6) = CDATE(5:7)
        STRING(3:3) = '-'
        STRING(8:9) = CDATE(23:24)
        ILEN = LEN(STRING)
        IF (ILEN .GT. 9) STRING(10:ILEN) = ' '
 
        RETURN
        END

        INTEGER FUNCTION JIABS(IT)
        INTEGER IT
        JIABS = ABS(IT)
        RETURN
        END

        INTEGER FUNCTION JMOD(IT1,IT2)
        JMOD = MOD(IT1,IT2)
        RETURN
        END

        INTEGER FUNCTION JINT(FVAL)
        INTEGER FVAL
        JINT = INT(FVAL)
        RETURN
        END

        INTEGER FUNCTION JNINT(FVAL)
        REAL FVAL
        JNINT = NINT(FVAL)
        RETURN
        END

        INTEGER FUNCTION IISHFT(IT1,IT2)
        INTEGER * 2 IT1,IT2
        IISHFT = ISHFT(IT1,IT2)
        RETURN
        END

#endif
