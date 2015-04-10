
C++*********************************************************************
C
C  CHARINSIDE.F                 NEW            MAY 1997 ARDEAN LEITH
C                               ADDED FROMBACK JAN 2001 ARDEAN LEITH
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
C CHARINSIDE(STRING,CLEFT,CRIGHT,INCLUSIVE,FROMBACK,IGO,IEND,NCHAR)
C
C PARAMETERS:    STRING     INPUT STRING                 (SENT)
C                CLEFT      START DELIMITER              (SENT)
C                CRIGHT     END   DELIMITER              (SENT)
C                INCLUSIVE  FLAG TO EXCLUDE DELIMITERS   (SENT)
C                FROMBACK   FLAG TO FIND FROM BACK       (SENT)
C                IGO        STARTING LOCATION            (RETURNED)
C                IEND       ENDING LOCATION              (RETURNED)
C                NCHAR      NUMBER OF CHAR.              (RETURNED)
C 
C PURPOSE:     FINDS FIRST OCCURANCE OF SUBSTRING DELIMINATED BY 
C              CLEFT AND CRIGHT INSIDE OF STRING.  IF LOGICAL FLAG
C              FROMBACK IS TRUE THEN HUNTS FOR LAST OCCURANCE.
C              IF INCLUSIVE THEN THE POSTIONS RETURNED IN IGO..IEND
C              DO NOT INCLUDE THE CLEFT & CRIGHT DELIMITERS.  NCHAR
C              IS NUMBER OF CHARACTERS IN IGO...IEND.  IF THERE
C              ARE MULTIPLE OCCURANCES OF THE DELIMITERS WITHIN
C              THE SUBSTRING THEN THE INNERMOST SET IS CHOSEN
C              E.G. IF FROMBACK IS FALSE AND STRING IS CLEFT, CLEFT,
C              CRIGHT THEN THE RETURNED STRING IS CLEFT,CRIGHT ONLY
C
C--*********************************************************************

       SUBROUTINE CHARINSIDE(STRING,CLEFT,CRIGHT,INCLUSIVE,FROMBACK,
     &                       IGO,IEND,NCHAR)

       INCLUDE 'CMBLOCK.INC'

       CHARACTER *(*) STRING
       CHARACTER *(*) CLEFT
       CHARACTER *(*) CRIGHT

       LOGICAL        INCLUSIVE,FROMBACK

C      SUBSTRING NOT PRESENT (IS DIFFERENT FROM EMPTY SUBSTRING)
       NCHAR = -1

       LENT  = LEN(STRING)

       IF (.NOT. FROMBACK) THEN
C         SEARCH FROM FRONT END OF STRING

C         FIND BEGINNING OF SUBSTRING (CLEFT) INSIDE STRING
          IGO = INDEX(STRING(1:LENT),CLEFT)

C         RETURN IF  START OF SUBSTRING (CLEFT) NOT PRESENT
          IF (IGO .LE. 0 ) RETURN

C         RETURN IF START OF SUBSTRING (CLEFT) IS AT END OF STRING
          IF (IGO .GE. LENT) RETURN   

C         FIND END OF SUBSTRING (CRIGHT) FOLLOWING IGO
          IEND  = INDEX(STRING(IGO+1:),CRIGHT) + IGO

C         RETURN IF END OF SUBSTRING  (CRIGHT) WAS NOT PRESENT 
          IF (IEND .LE. IGO) RETURN


C         FIND LAST BEGINNING OF SUBSTRING (CLEFT) INSIDE NEW STRING
          IGO2 = INDEX(STRING(IGO+1:IEND-1),CLEFT,BACK=.TRUE.)
          IGO  = IGO + IGO2

       ELSE
C         SEARCH FROM BACK END OF STRING

C         FIND END OF SUBSTRING (CRIGHT) INSIDE STRING
          IEND = INDEX(STRING(1:LENT),CRIGHT,BACK=FROMBACK)

C         RETURN IF  END OF SUBSTRING (CRIGHT) NOT PRESENT
          IF (IEND .LE. 1) RETURN

C         FIND BEGINNING OF SUBSTRING (CLEFT) PRECEDING IEND
          IGO = INDEX(STRING(:IEND-1),CLEFT,BACK=FROMBACK)

C         RETURN IF BEGINNING OF SUBSTRING (CLEFT) NOT PRESENT 
          IF (IGO .LE. 0) RETURN

C         FIND FIRST END OF SUBSTRING (CRIGHT) INSIDE NEW STRING
          IEND2 = INDEX(STRING(IGO:IEND-1),CRIGHT)
          IF (IEND2 .GT. 0) IEND  = IGO + IEND2 - 1

       ENDIF

       IF (INCLUSIVE) THEN
C         DO NOT WANT TO INCLUDE DELIMITER LOCATIONS
          IGO  = IGO  + 1
          IEND = IEND - 1
       ENDIF

C      FIND NUMBER OF CHARACTERS IN RETURNED SUBSTRING
       NCHAR = IEND - IGO + 1

       RETURN
       END


