
C++*********************************************************************
C
C GETOLDSTACK.F    NEW                           7 JAN 99 ARDEAN LEITH
C                  LUNSETCOMMON BUG FIXED       17 OCT 02 ARDEAN LEITH
C                  INDEXED STACK                   JAN 02 ARDEAN LEITH
C                  MAXNAM                          JUL 14 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C    GETOLDSTACK(LUN,IMGNUM,WANTNEXT,MUSTGET,SAYIT,IRTFLG)
C
C    PURPOSE:       TO OPEN A SPECIFIED IMAGE WITHIN STACK FOR RANDOM 
C                   ACCESS READING/WRITING.
C
C    PARAMETERS:
C        LUN        LOGICAL UNIT NUMBER FOR FILNAM.               (SENT)
C        IMGNUM     IMAGE NUMBER WANTED AND FOUND            (SENT/RET.) 
C        WANTNEXT   LOGICAL FLAG TO GET NEXT IMAGE >= IMGNUM      (SENT) 
C        MUSTGET    LOGICAL FLAG TO STOP IF IMAGE DOES'NT EXIST   (SENT) 
C        SAYIT      SAY FILE OPENIN INFO.                         (RET.)
C        IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                   IRTFLG =  3    ERROR, END OF FILE BEFORE IMGNUM
C                   IRTFLG =  2    ERROR, MUSTGET BUT IMGNUM UNUSED
C                   IRTFLG =  1    ERROR, REDHED FAILED
C                   IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                   IRTFLG = -1    NOT A STACK
C                   IRTFLG = -2    IMAGE NOT IN USE, MUSTGET IS FALSE
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE GETOLDSTACK(LUN,NSAM,IMGNUM,WANTNEXT,MUSTGET,
     &                         SAYIT,IRTFLG)

        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=1)      :: DSP
        CHARACTER(LEN=MAXNAM) :: FILNAM

        LOGICAL               :: WANTNEXT,MUSTGET,SAYIT,ISBARET

        CALL LUNGETISBARE(LUN,ISBARET,IRTFLG)
        IF (.NOT. ISBARET) THEN
C          INPUT IMAGE IS SIMPLE OR A SPECIFIC STACKED IMAGE
C          IF SAYIT, WRITE OUT FILE OPENING INFO 
           IF (SAYIT) CALL LUNSAYINFO(LUN,IRTFLG)
           IRTFLG = -1
           RETURN
        ENDIF

C       V---- ONLY FOR INLINE OR REGULAR IMAGE STACK NOW ------------V

C       GET SPECIFIED IMAGE HEADER FROM STACK FILE LOCATION
C       DO NOT CALL ERRT IF RUNS OFF END OF FILE
C       MUST LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)

8888    CALL LUNREDHED(LUN,NSAM,0,.FALSE.,IRTFLG)
        CALL LUNREDHED(LUN,NSAM,IMGNUM,.FALSE.,IRTFLG)
        IF (IRTFLG .GT. 0) THEN
C          PROBABLY RAN OFF END OF STACK FILE, ERRT NOT CALLED
           IF (MUSTGET) 
     &         CALL ERRT(102,'THIS IMAGE NOT USED IN STACK',IMGNUM) 
           IRTFLG = 3
           RETURN
        ELSEIF (IRTFLG .EQ. 0) THEN
C          NEED IMUSED FROM THIS STACKED IMAGE
           CALL LUNGETINUSE(LUN,IMUSED,IRTFLG)
        ELSE
           IMUSED = 0
        ENDIF
      
        IF (IMUSED .EQ. 0) THEN
C          THIS IMAGE NOT AN EXISTING IMAGE WITHIN STACK!
           IF (WANTNEXT) THEN
C             INCREMENT IMGNUM AND TRY AGAIN
              IMGNUM = IMGNUM + 1
              GOTO 8888

           ELSEIF (MUSTGET) THEN
              CALL ERRT(102,'THIS IMAGE NOT USED IN STACK',IMGNUM)
              IRTFLG = 2
              RETURN
           ELSE
              IRTFLG = -2
              RETURN
           ENDIF
        ENDIF

C       GET FILENAM FROM CURRENT HEADER OBJECT
        CALL LUNGETFILE(LUN,FILNAM,NLET,DSP,IRTFLG)

C       APPEND CURRENT STACKED IMAGE NUMBER TO STACK FILE NAME
C       (INTTOCHAR ALSO RETURNS NEW VALUE FOR NLET) 
        LENAT = INDEX(FILNAM,'@')
        CALL INTTOCHAR(IMGNUM,FILNAM(LENAT+1:),NLET,0)
        NLET = NLET + LENAT

C       SET FILENAME IN HEADER OBJECT
        CALL LUNSETFILE(LUN,FILNAM(1:NLET),'O',IRTFLG)

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUN,IMGNUM,NSAM,IRTFLG)

C       WRITE OUT FILE OPENING INFO 
        CALL LUNSAYINFO(LUN,IRTFLG)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON(LUN,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
 
        RETURN
	END



