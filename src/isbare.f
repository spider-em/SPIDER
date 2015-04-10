
C++*********************************************************************
C
C    ISBARE.F                                                NEW JAN 98
C                                                  AUTHOR: ARDEAN LEITH
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
C      PURPOSE: FIND IF A FILENAME DENOTES A BARE STACK
C
C      PARAMETER: FILNAM   INPUT FILENAME
C
C **********************************************************************

       LOGICAL FUNCTION ISBARE(FILNAM)

 

       CHARACTER *(*) FILNAM

C      FIND LOCATION OF @ IN FILNAM
       ILOCAT = INDEX(FILNAM,'@')

C      FIND LOCATION OF LAST NON-NULL NON-BLANK CHARACTER IN FILE 
       LENB   = lnblnk(FILNAM)
       LENN   = INDEX(FILNAM,CHAR(0)) 
       IF (LENN .GT. 1 .AND. LENN .LE. LENB) LENB = LENN - 1

C      IF IT IS A BARE STACK ILOCAT IS LENB
       ISBARE =  ILOCAT .GT. 1 .AND. ILOCAT .EQ. LENB

       RETURN
       END
