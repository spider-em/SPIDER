C ++********************************************************************
C                                                                      *
C DATE                                  NEW APRIL 98 FOR F90 al        *
C                                                                      *
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
C  DATE(DATEVAR)
C
C  PARAMETERS:   DATEVAR    CHAR. VARIABLE FOR DATE        (RETURNED)
C                           FORMAT IS DD-MON-YY
C
C --*********************************************************************

       SUBROUTINE DATE_DUMMY(DATEVAR)

       CHARACTER*(*)         :: DATEVAR

       CHARACTER(LEN=8)      :: DATET
       CHARACTER(LEN=10)     :: TIME
       CHARACTER(LEN=5)      :: ZONE
       INTEGER, DIMENSION(8) :: VALUES
 
       CALL DATE_AND_TIME(DATET,TIME,ZONE,VALUES)

       DATEVAR(1:3) = DATET(7:8) // '-'

       IF     (DATET(1:2) .EQ. '01') THEN
           DATEVAR(4:7) = 'JAN-'
       ELSEIF (DATET(1:2) .EQ. '02') THEN
           DATEVAR(4:7) = 'FEB-'
       ELSEIF (DATET(1:2) .EQ. '03') THEN
           DATEVAR(4:7) = 'MAR-'
       ELSEIF (DATET(1:2) .EQ. '04') THEN
           DATEVAR(4:7) = 'APR-'
       ELSEIF (DATET(1:2) .EQ. '05') THEN
           DATEVAR(4:7) = 'MAY-'
       ELSEIF (DATET(1:2) .EQ. '06') THEN
           DATEVAR(4:7) = 'JUN-'
       ELSEIF (DATET(1:2) .EQ. '07') THEN
           DATEVAR(4:7) = 'JUL-'
       ELSEIF (DATET(1:2) .EQ. '08') THEN
           DATEVAR(4:7) = 'AUG-'
       ELSEIF (DATET(1:2) .EQ. '09') THEN
           DATEVAR(4:7) = 'SEP-'
       ELSEIF (DATET(1:2) .EQ. '10') THEN
           DATEVAR(4:7) = 'OCT-'
       ELSEIF (DATET(1:2) .EQ. '11') THEN
           DATEVAR(4:7) = 'NOV-'
       ELSEIF (DATET(1:2) .EQ. '12') THEN
           DATEVAR(4:7) = 'DEC-'
       ENDIF
      
       DATEVAR(8:9) = DATET(3:4) 

       IF (LEN(DATEVAR) .GE. 10) DATEVAR(10:10) = CHAR(0)

       RETURN
       END

