
C++*********************************************************************
C
C ERRF2.F  -- CREATED JULY 15 1987
C
C **********************************************************************
C *  AUTHOR:  ARDEAN LEITH
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
C    ERRF2(FNUM1,FNUM2,NVAL,FLOW1,FHI1,FLOW2,FHI2)
C
C    PARAMETERS:    NVAL         NUMBER OF VALUES TO BE CHECKED
C                   FNUM1        FIRST VALUE
C                   FNUM2        SECOND VALUE
C                   FLOW1        LOWEST VALUE FOR  FNUM1
C                   FHI1         HIGHEST VALUE FOR FNUM1
C                   FLOW2        LOWEST VALUE FOR  FNUM2
C                   FHI2         HIGHEST VALUE FOR FNUM2
C
C    CALLED BY:     CSAXIS (ONLY)
C
C        0         2         3         4         5         6         7
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       LOGICAL FUNCTION  ERRF2(FNUM1,FNUM2,NVAL,FLOW1,FHI1,FLOW2,FHI2)

 

       COMMON /UNITS/ LUNDOC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

        ERRF2 = .TRUE.
        IF (FNUM1 .LT. FLOW1 .OR. FNUM1 .GT. FHI1) THEN
            WRITE(NOUT,92)FLOW1,FHI1
   92       FORMAT(' ERROR, FIRST INPUT RANGE: (',1PG11.3,'....',
     &               1PG11.3,')',/)
            RETURN
        ENDIF

        IF (NVAL .GE. 2 .AND. (FNUM2 .LT. FLOW2 .OR. FNUM2 .GT. FHI2)) 
     &      THEN
            WRITE(NOUT,97) FLOW2,FHI2
   97       FORMAT(' ERROR, SECOND INPUT RANGE: (',1PG11.3,'....',
     &               1PG11.3,')',/)
            RETURN
        ENDIF

        ERRF2 = .FALSE.

        RETURN
        END
