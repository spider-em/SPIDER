

C++*************************************************************************
C
C  SUBSYMPAR.F    CREATED  FROM SPIDER.F         Sep 2000  ARDEAN LEITH 
C                 MULTIPLE VARIABLE SUBSTITUTION JAN 2001  ARDEAN LEITH
C                 ALPHABETICAL VARIABLES         JUN 2002  ARDEAN LEITH
C                 NESTED VARIABLES               SEP 2002  ARDEAN LEITH
C                 [] DEFAULT FOR VARIABLES       OCT 2005  ARDEAN LEITH
C                 CVAR                           OCT 2006  ARDEAN LEITH
C                 DONOTRECURSE
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
C    SUBSYMPAR(INPUT,OUTPUT,NCHAR,ILEVELT,IRTFLG)
C
C    PURPOSE:       RUN-TIME VARIABLE SUBSTITUTION FOR ALL [ID] IN
C                   INPUT STRING AT THIS STACK LEVEL
C
C    PARAMETERS:    INPUT     INPUT LINE CONTAINING [STRING..] (SENT)
C                   OUTPUT    SUBSTITUTED OUTPUT LINE         (RET.)
C                   NCHAR     LAST NON_BLANK CHAR BEFORE ;    (RET.)
C                   ILEVELT   NESTING LEVEL                   (SENT)
C                   IRTFLG    RETURN FLAG (0 IS NORMAL)       (RETURNED)
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE SUBSYMPAR(INPUT,OUTPUT,NCHAR,ILEVELT,IRTFLG)

      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*), INTENT(IN)  :: INPUT
      CHARACTER(LEN=*), INTENT(OUT) :: OUTPUT
      LOGICAL, PARAMETER            :: DONOTRECURSE = .FALSE.

C     FOR ISTOP 
      INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
      COMMON /QSTR_STUFF1/ ISTOP,ITI,ITIN,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      ILEVEL = ILEVELT
      IF (ILEVEL .LE. 0) ILEVEL = ISTOP

      CALL SYMPAR_SUB(INPUT,OUTPUT,NCHAR,ILEVEL,DONOTRECURSE,IRTFLG)

      END
   
