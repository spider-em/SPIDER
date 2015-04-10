
C++*********************************************************************
C
C GETSTACK.F  -- CREATED 10 SEP 02 ARDEAN LEITH
C                  USED GETOLDSTACK, GETNEWSTACK  APRIL 99 ARDEAN LEITH
C                  INDEXED STACK                  JAN 02   ARDEAN LEITH
C                  GETNEWSTACK PARAM.             FEB 03   ARDEAN LEITH
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
C  GETSTACK(LUN1,LUN2,IMGNUM,MAXIM,SAYIT,LOADIT,BUF,NORMIT,IRTFLG)
C
C  PURPOSE:  LOCATE NEXT IMAGE/VOL. IN INPUT & OUTPUT STACK & OPTIONALLY
C            LOAD IT INTO BUF.
C 
C  PARAMETERS:  IMGNUM    LAST IMAGE                        (SENT/RET.)
C                 USE -3 ON INPUT FOR NON STACK LOOP END         (SENT)
C               MAXIM     STACK FLAG                             (SENT)
C               SAYIT     REPORT FILE OPENING INFO               (SENT)
C               LOADIT    LOAD IMAGE/VOL. IN BUF                 (SENT)
C               BUF       IMAGE/VOLUME AREA                      (RET.)
C               NORMIT    FIND IMAGE/VOL. MIN,MAX....            (SENT)
C               IRTFLG    ERROR FLAG                             (RET.)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE GETSTACK(LUN1,LUN2,IMGNUM,MAXIM,
     &                     SAYIT,LOADIT,BUF,NORMIT,IRTFLG)

       INCLUDE    'CMBLOCK.INC' 

       LOGICAL  ::  WANTNEXT,MUSTGET,SAYIT,LOADIT,NORMIT

       WANTNEXT = .TRUE.
       MUSTGET  = .FALSE.

C      GET CURRENT NSAM,.......
       CALL LUNGETSIZE(LUN1,NSAM,NROW,NSLICE,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       IF (MAXIM .LT. 0) THEN
C         NOT A WHOLE STACK OPERATION
          IMGNUM = MAXIM 

       ELSE
C         NEGATIVE IMGNUM SET TO ZERO
          IMGNUM = MAX(IMGNUM+1,1)

C         GET INPUT IMAGE FROM INPUT STACK            
          CALL GETOLDSTACK(LUN1,NSAM,IMGNUM,
     &                     WANTNEXT,MUSTGET,SAYIT,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

C         CREATE OUTPUT IMAGE IN OUTPUT STACK
          CALL GETNEWSTACK(LUN1,LUN2,NSAM,IMGNUM,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

       ENDIF

C      SET CURRENT FMIN.......
       CALL LUNGETSTAT(LUN1,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)

       IF (IMAMI.NE.1 .AND. NORMIT) THEN
          CALL NORM3(LUN1,NSAM,NROW,NSLICE,FMAX,FMIN,AV)
          IF (IRTFLG .NE. 0) RETURN
       ENDIF

       IF (LOADIT) THEN
C         LOAD IMAGE/VOL. IN BUF
          CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,BUF,IRTFLG)
       ENDIF

       RETURN
       END
