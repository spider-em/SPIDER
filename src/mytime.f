C ++********************************************************************
C                                                                      *
C MYTIME                                  NEW APRIL 98 FOR F90 al        *
C                        RENAMED FROM TIME DEC 00  ARDEAN LEITH
C                                                          *
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
C  MYTIME(TIMEVAR)
C
C  PARAMETERS:   TIMEVAR    CHAR. VARIABLE FOR TIME        (RETURNED)
C                           FORMAT IS HH-MM-SS
C
C --*********************************************************************

       SUBROUTINE  MYTIME(TIMEVAR)

C      USUAL RETURNED LENGTH OF TIMEVAR IF CTIM IS 8

       CHARACTER  *(*)       :: TIMEVAR

C      CHARACTER(LEN=8)      :: DATET (fails on altix)
       CHARACTER(LEN=64)     :: DATET
       CHARACTER(LEN=10)     :: TIMET
       CHARACTER(LEN=5)      :: ZONE
       INTEGER, DIMENSION(8) :: VALUES

C      CALL DATE_AND_TIME(DATET,TIMET,ZONE,VALUES) (fails on altix)
       CALL DATE_AND_TIME(DATET,TIMET)

       TIMEVAR(1:8) = TIMET(1:2) // ':' // TIMET(3:4) // ':' // 
     &                TIMET(5:6) 

       IF (LEN(TIMEVAR) .GE. 9) TIMEVAR(9:9) = CHAR(0)

       RETURN
       END

