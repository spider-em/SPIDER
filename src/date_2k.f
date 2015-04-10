C++*********************************************************************
C
C DATE_2K.F   -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
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
C    DATE_2K(DATEVAR)
C
C    PURPOSE:       RETURNS DATE IN: 21-JAN-1999   FORMAT
C
C    PARAMETERS :   DATEVAR  CHAR*12 CONTAINING DATE        (RETURNED)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE DATE_2K(DATEVAR)

C      DATEVAR SHOULD BE 12 CHAR. LONG!
       CHARACTER *(*)    DATEVAR

       CHARACTER *8      DATET
       CHARACTER *10     TIME
       CHARACTER *5      ZONE
       INTEGER           IVALUES(8)

C      Y2K SAFE CALL

C       CALL DATE_AND_TIME(DATET,TIME,ZONE,IVALUES) (fails on altrix)
       CALL DATE_AND_TIME(DATET)

C      PUT IN DAY
       DATEVAR(1:3) = DATET(7:8) // '-'

C      PUT IN MONTH
       IF     (DATET(5:6) .EQ. '01') THEN
           DATEVAR(4:7) = 'JAN-'
       ELSEIF (DATET(5:6) .EQ. '02') THEN
           DATEVAR(4:7) = 'FEB-'
       ELSEIF (DATET(5:6) .EQ. '03') THEN
           DATEVAR(4:7) = 'MAR-'
       ELSEIF (DATET(5:6) .EQ. '04') THEN
           DATEVAR(4:7) = 'APR-'
       ELSEIF (DATET(5:6) .EQ. '05') THEN
           DATEVAR(4:7) = 'MAY-'
       ELSEIF (DATET(5:6) .EQ. '06') THEN
           DATEVAR(4:7) = 'JUN-'
       ELSEIF (DATET(5:6) .EQ. '07') THEN
           DATEVAR(4:7) = 'JUL-'
       ELSEIF (DATET(5:6) .EQ. '08') THEN
           DATEVAR(4:7) = 'AUG-'
       ELSEIF (DATET(5:6) .EQ. '09') THEN
           DATEVAR(4:7) = 'SEP-'
       ELSEIF (DATET(5:6) .EQ. '10') THEN
           DATEVAR(4:7) = 'OCT-'
       ELSEIF (DATET(5:6) .EQ. '11') THEN
           DATEVAR(4:7) = 'NOV-'
       ELSEIF (DATET(5:6) .EQ. '12') THEN
           DATEVAR(4:7) = 'DEC-'
       ENDIF
    
C      PUT IN YEAR  
       DATEVAR(8:11) = DATET(1:4) 

       DATEVAR(12:12) = CHAR(0)

       RETURN
       END

