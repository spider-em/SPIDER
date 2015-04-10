
C++*********************************************************************
C
C  CHKINPQ.FOR           MODIFIED FOR CHAR. FROM CHKINP.FOR AUG 89 al
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
C     CHKINPQ(SMENU,WHICH,IOUT)
C
C     PARAMETERS:    SMENU         MENU TO BE SEARCHED
C                    WHICH         STRING BEING SEARCHED FOR
C                    IOUT          LOCATION OF STRING IN MENU (RETURNED)
C
C--*******************************************************************

       SUBROUTINE CHKINPQ(SMENU,WHICH,IOUT)

 

       CHARACTER *(*) SMENU,WHICH
       CHARACTER *80 TEMP

       IOUT   = 1
       LENMEN = LEN(SMENU)
       LENDOL = INDEX(SMENU,'$')
       IF (LENDOL .GT. 0 .AND. LENDOL .LE. LENMEN) LENMEN = LENDOL - 1
       DO  LENUSE = LENMEN,1,-1
          IF (SMENU(LENUSE:LENUSE) .NE. ' ' ) GOTO 5
       ENDDO

C      MENU IS ALL BLANK
       RETURN

5      IF (WHICH .NE. ' ') THEN
C         SEARCH STRING IS NOT BLANK, FIND USED LENGTH OF SEARCH STRING 
          LENS   = LEN(WHICH)
          IBLANK = INDEX(WHICH,' ')
          IF (IBLANK .EQ. 0) IBLANK = LENS + 1
          ICOMMA = INDEX(WHICH,',')
          IF (ICOMMA .EQ. 0) ICOMMA = LENS + 1
          LENGTH = MIN(IBLANK-1,ICOMMA-1,LENS)

C         FIND LOCATION OF SEARCH STRING IN MENU
          IGO  = 1
          IOUT = 1

8         IOUT = IOUT + 1

C         FIND RANGE OF NEXT MENU CHOICE
          IEND = INDEX(SMENU(IGO:),',')
          IF (IEND .EQ. 0) THEN
C            LAST MENU CHOICE IN THE LIST
             IEND = LENMEN
          ELSE
             IEND = IGO + IEND - 2
          ENDIF

C         DO UPPER CASE COMPARISON TO MENU...
          TEMP(1:LENGTH) = WHICH(1:LENGTH)
          CALL SSUPCAS(TEMP)
          IF (SMENU(IGO:IEND) .EQ. WHICH(1:LENGTH)) THEN
C            HAVE FOUND WHICH IN SMENU
             RETURN
          ELSE
C            THIS MENU CHOICE IS NOT WHICH, TRY NEXT MENU CHOICE, IF ANY
             IGO = IEND + 2
             IF (IGO .LE. LENUSE) GOTO 8
          ENDIF
       ENDIF

C      WHICH NOT FOUND IN SMENU
       IOUT = 1
       RETURN
       END

