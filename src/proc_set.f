C++*********************************************************************
C
C PROC_SET.F   -- NEW JAN 2001                   AUTHOR: ARDEAN LEITH
C                 PLINE LENGTHENED MAR 2001 ARDEAN LEITH
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
C  PROC_SET(NUMPRC,NCHARS,NLINES,PLINEGO,PDATA,IRTFLG) 
C  PROC_GETPLINE(IPLINE,IPROCNUMT,PLINE,NCHAR,IRTFLG)
C
C  PURPOSE:  HANDLES ALL INTERACTIONS WITH SPIDER PROCEDURE FILE
C            CONTAINS SUBROUTINES STARTING  WITH PREFIX: PROC_
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      MODULE PROC_INFO

        INCLUDE 'CMLIMIT.INC'

C       ARRAY OF POINTERS TO PROCEDURE STORAGES
        TYPE CHARA_POINTER
           CHARACTER, DIMENSION(:), POINTER :: CPTR 
        END TYPE CHARA_POINTER

        TYPE INTA_POINTER
           INTEGER, DIMENSION(:), POINTER :: IPTR 
        END TYPE INTA_POINTER

        TYPE(CHARA_POINTER), DIMENSION(MAXPRCNAM) :: PROCSTOREPTRS
        TYPE(INTA_POINTER) , DIMENSION(MAXPRCNAM) :: PROCLINESPTRS


      END MODULE PROC_INFO

C      ******************** PROC_SET *****************************

       SUBROUTINE PROC_SET(NUMPRC,NCHARS,NLINES,PLINEGO,PDATA,IRTFLG)

       USE PROC_INFO

       INCLUDE 'CMBLOCK.INC'

       INTEGER, DIMENSION(NLINES)   :: PLINEGO
       CHARACTER, DIMENSION(NCHARS) :: PDATA
 
       CHARACTER, DIMENSION(:), POINTER :: PROCSTORE 
       INTEGER,   DIMENSION(:), POINTER :: PROCLINES

C      ALLOCATE NEW PROC STORAGE ARRAYS FOR THIS LENGTH PROCEDURE
       ALLOCATE(PROCSTORE(NCHARS),PROCLINES(NLINES+1),STAT=IRTFLG)
       IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'SPIDER,PROCSTORE...',NE)
           RETURN
       ENDIF

C      KEEP POINTERS
       PROCSTOREPTRS(NUMPRC)%CPTR => PROCSTORE
       PROCLINESPTRS(NUMPRC)%IPTR => PROCLINES

C      FILL ARRAY HOLDING START LOCATION OF EACH LINE
       DO I = 1, NLINES
          PROCLINES(I) = PLINEGO(I)
       ENDDO
       PROCLINES(NLINES+1) = -NCHARS

C      FILL ARRAY HOLDING PROC. LINES
       DO I = 1, NCHARS
          PROCSTORE(I) = PDATA(I) 
       ENDDO

       IRTFLG = 0
       RETURN
       END



C      ******************** PROC_GETPLINE *****************************

       SUBROUTINE PROC_GETPLINE(IPLINE,IPROCNUMT,PLINE,NCHAR,IRTFLG)

       USE PROC_INFO

       INCLUDE 'CMBLOCK.INC'

       COMMON /PROC_STUFF/ NUMPROCNOW

       CHARACTER(LEN=*)                 :: PLINE
       CHARACTER, DIMENSION(:), POINTER :: PROCSTORE 
       INTEGER,   DIMENSION(:), POINTER :: PROCLINES

C      RETRIEVES CURRENT PROCEDURE LINE

       IRTFLG    = 1
       NCHAR     = 0

       IF (IPROCNUMT .GT. 0) THEN
          IPROCNUM  = IPROCNUMT
       ELSE
          IPROCNUM  = NUMPROCNOW
       ENDIF

       PROCSTORE => PROCSTOREPTRS(IPROCNUM)%CPTR
       PROCLINES => PROCLINESPTRS(IPROCNUM)%IPTR

       IGO = PROCLINES(IPLINE)
C      CHECK TO SEE IF VALID LINE
       IF (IGO .LE. 0) RETURN

       IEND  = ABS(PROCLINES(IPLINE+1))
       NCHAR = IEND - IGO

       IRTFLG = 0
       IF (NCHAR .LE. 0) THEN
          PLINE(1:) = ' '
          RETURN
       ENDIF

       NCHARMAX = LEN(PLINE)
       NCHAR    = MIN(NCHAR,NCHARMAX)
       DO I = 1,NCHAR
          PLINE(I:I) = PROCSTORE(IGO+I-1)
       ENDDO
       IF (NCHAR .LT. NCHARMAX) PLINE(NCHAR+1:) = ' '

       RETURN
       END


