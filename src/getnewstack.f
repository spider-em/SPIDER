
C++*********************************************************************
C
C GETNEWSTACK.F  NEW                            7 JAN 99   ARDEAN LEITH
C                SET INUSE                        OCT 99   ARDEAN LEITH
C                INDEXED STACK                    JAN 03   ARDEAN LEITH
C                HEADER COPY                      FEB 03   ARDEAN LEITH
C                COPYSTAT PARAM                   OCT 2010 ARDEAN LEITH
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
C    GETNEWSTACK(LUNIN,LUNOUT,COPYSTAT,IMGNUM,IRTFLG)
C
C    PURPOSE:       TO OPEN A NEW IMAGE WITHIN STACK AND TIE TO LUNOUT
C
C    PARAMETERS:
C        LUNIN      UNIT TO COPY HEADER VALUES FROM               (SENT)
C        LUNOUT     LOGICAL UNIT NUMBER FOR FILE.                 (SENT)
C        COPYSTAT   COPY FMIN.. FROM INPUT FILE                   (SENT)
C        IMGNUM     IMAGE NUMBER WANTED                           (SENT) 
C        IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                   IRTFLG =  1    ERROR, REDHED FAILED
C                   IRTFLG =  0    NORMAL RETURN, IMAGE IS STACKED
C                   IRTFLG = -1    OK, BUT FILE REQUESTED IS NOT STACK
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE GETNEWSTACK(LUNIN,LUNOUT,COPYSTAT,NSAM,IMGNUM,IRTFLG)

        INTEGER            :: LUNIN,LUNOUT,NSAM,IMGNUM,IRTFLG
        LOGICAL            :: COPYSTAT

        CHARACTER (LEN=71) :: FILNAM
        CHARACTER (LEN=1)  :: DSP

        LOGICAL            :: ISBARET

        CALL LUNGETISBARE(LUNOUT,ISBARET,IRTFLG)
        IF (.NOT. ISBARET) THEN
C          INPUT IMAGE IS SIMPLE OR A SPECIFIC STACKED IMAGE
           IRTFLG = -1
           RETURN
        ENDIF

C       ONLY FOR INLINE OR REGULAR IMAGE STACK NOW -------------V
C       MAKE A NEW IMAGE WITHIN EXISTING STACK 

C       GET OVERALL HEADER FROM THE STACK FILE
        CALL LUNREDHED(LUNOUT,NSAM,0,.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           IRTFLG = 1
           RETURN
        ENDIF

C       RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
        CALL LUNGETMAXIM(LUNOUT,MAXIM,IRTFLG)

        IF (IMGNUM .GT. MAXIM) THEN
C          MUST UPDATE OVERALL HEADER WITH MAXIMUM IMAGE NUMBER
           CALL LUNSETMAXIM(LUNOUT,IMGNUM,IRTFLG)
           CALL LUNSETMAXALL(LUNOUT,IMGNUM,IRTFLG)
        ENDIF

C       NEED SOME VALUES TO PUT IN STACKED IMAGE HEADER OBJECT
        CALL LUNGETSIZE(LUNOUT,NSAM,NROW,NSLICE,IRTFLG)
        CALL LUNGETTYPE(LUNOUT,ITYPE,IRTFLGT)
        CALL LUNCOPYSTK(LUNOUT,ISTACK,IRTFLGT)

        IF (ISTACK .LT. 0) THEN
C          MAKING A NEW INDEXED STACKED FILE, UPDATE INDX LOCATION
           CALL LUNWRTINDX(LUNOUT,IMGNUM,NSAM,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN
        ENDIF

        IF (IMGNUM .GT. MAXIM .OR. ISTACK .GT. 2) THEN
C          SAVE OVERALL HEADER NOW TO PRESERVE MAXIM & LASTINDX
           CALL LUNWRTHED(LUNOUT,NSAM,0,IRTFLGT)
        ENDIF

C       GET FILENAM FROM CURRENT HEADER OBJECT
        CALL LUNGETFILE(LUNOUT,FILNAM,NLET,DSP,IRTFLG)

C       APPEND CURRENT STACKED IMAGE NUMBER TO STACK FILE NAME
C       (INTTOCHAR ALSO RETURNS NEW VALUE FOR NLET) 
        LENAT = INDEX(FILNAM,'@')
        CALL INTTOCHAR(IMGNUM,FILNAM(LENAT+1:),NLET,0)
        NLET = NLET + LENAT

C       INITIALIZE HEADER OBJECT FOR NEW STACKED IMAGE FROM LUNIN
C         (I THINK THIS COULD BE MOSTLY SKIPPED??)
        CALL LUNSETHDR(LUNIN,LUNOUT,NSAM,NROW,NSLICE,
     &                 ITYPE,ISTACK,IRTFLG)

        IF (COPYSTAT) THEN
C          COPY STATISTICS FROM INPUT FILE
           CALL LUNGETSTAT(LUNIN, IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)
           CALL LUNSETSTAT(LUNOUT,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)
        ENDIF

C       PUT IMAGE INUSE IN HEADER OBJECT
        CALL LUNSETINUSE( LUNOUT,IMGNUM,IRTFLG)
        CALL LUNSETISTACK(LUNOUT,0,IRTFLG)

C       PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
        CALL LUNSETCOMMON(LUNOUT,IRTFLG)

C       PUSH HEADER OBJECT INFO INTO NEW STACKED FILE
        CALL LUNWRTHED(LUNOUT,NSAM,IMGNUM,IRTFLG)

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUNOUT,IMGNUM,NSAM,IRTFLGT)

C       WRITE OUT FILE OPENING INFO TO SCREEN
        CALL LUNSAYINFO(LUNOUT,IRTFLG)

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
 
        RETURN
	END


