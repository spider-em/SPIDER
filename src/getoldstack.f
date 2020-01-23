
C++*********************************************************************
C
C GETOLDSTACK.F  NEW                            7 JAN 99 ArDean Leith
C                LUNSETCOMMON BUG FIXED        17 OCT 02 ArDean Leith
C                INDEXED STACK                    JAN 02 ArDean Leith
C                MAXNAM                           JUL 14 ArDean Leith
C                MRC SUPPORT                      SEP 19 ArDean Leith
C                LOCAT BUG IN MRC              JAN 2020   ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  GETOLDSTACK(LUN,IMGNUM,WANTNEXT,MUSTGET,SAYIT,IRTFLG)
c  GETOLDSTACK_MRC(LUN,NX,IMGNUM,WANTNEXT,IRTFLG)
C
C  PURPOSE:       OPEN A SPECIFIED IMAGE WITHIN STACK FOR RANDOM 
C                 ACCESS READING/WRITING.
C
C  PARAMETERS:
C      LUN        LOGICAL UNIT NUMBER FOR FILNAM.               (SENT)
C      IMGNUM     IMAGE NUMBER WANTED AND FOUND            (SENT/RET.) 
C      WANTNEXT   LOGICAL FLAG TO GET NEXT IMAGE >= IMGNUM      (SENT) 
C      MUSTGET    LOGICAL FLAG TO STOP IF IMAGE DOESN'T EXIST   (SENT) 
C      SAYIT      SAY FILE OPENIN INFO.                         (RET.)
C      IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                 IRTFLG =  3    ERROR, END OF FILE BEFORE IMGNUM
C                 IRTFLG =  2    ERROR, MUSTGET BUT IMGNUM UNUSED
C                 IRTFLG =  1    ERROR, REDHED FAILED
C                 IRTFLG =  0    NORMAL RETURN, IMAGE IS STACK
C                 IRTFLG = -1    NOT A STACK
C                 IRTFLG = -2    IMAGE NOT IN USE, MUSTGET IS FALSE
C
C  NOTE:      THE THREE SUBROUTINES: GETSTACK, GETOLDSTACK, AND
C             GETNEWSTACK WERE WRITTEN TO HANDLE WHOLE STACK (@STACK)
C             OPERATIONS. THEY ARE ONLY USED IN THE OPERATIONS:
C             'AF' (UNDOCUMENTED), 'CP TO STK','DE','CE AD', 
C             'CE' (MANY?), 'AR SCA',
C             'AD D2', 'AD DF', 'AD DR', 'AD F', 'AD R', 'DIV 2', 
C             'MUL 2', 'MUL O', 'MU O', & 'SU B2'.
C             I DECIDED THAT WHOLE STACK OPERATION SUPPORT WAS LESS
C             IMPORTANT THAN SUPPORTING OPERATIONS ON SELECTED 
C             STACKED IMAGES SO I STOPPED CONVERTING OPERATIONS TO
C             HANDLING WHOLE STACKS WITHOUT SELECTED STACK SUPPORT.
C             THE NEW PARADIGM USES OPTILES AND NEXTFILE WHICH
C             HANDLES BOTH SELECTED AND WHOLE STACKS. 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE GETOLDSTACK(LUN,NX,IMGNUM,WANTNEXT,MUSTGET,
     &                         SAYIT,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'    ! NEED: MAXNAM

	INTEGER               :: LUN,NX,IMGNUM,IRTFLG
        LOGICAL               :: WANTNEXT,MUSTGET,SAYIT

        CHARACTER(LEN=1)      :: DSP
        CHARACTER(LEN=MAXNAM) :: FILNAM

	INTEGER               :: IMUSED, NZ,MZ,NSTACK,NLET,LOCAT
        INTEGER               :: LOCSLASH,IGO,LENAT

        LOGICAL               :: ISBARET,IS_MRC
    
    
        CALL LUNGETISBARE(LUN,ISBARET,IRTFLG)     ! MRC OK
        IF (.NOT. ISBARET) THEN
C          INPUT IMAGE IS SIMPLE OR A SPECIFIC STACKED IMAGE
C          IF SAYIT, WRITE OUT FILE OPENING INFO 
           IF (SAYIT) CALL LUNSAYINFO(LUN,IRTFLG)
           IRTFLG = -1       ! NOT A STACK
           RETURN
        ENDIF

C       ONLY FOR INLINE OR REGULAR IMAGE STACK NOW 

C       DETERMINE IF MRC OR SPIDER 
        CALL LUNGETIS_MRC(LUN,IS_MRC,IRTFLG)

        IF (IS_MRC) THEN
C          MRC STACK
           CALL GETOLDSTACK_MRC(LUN,IMGNUM,WANTNEXT,IRTFLG)
           RETURN
        ENDIF


C       GET SPECIFIED IMAGE HEADER FROM STACK FILE LOCATION
C       DO NOT CALL ERRT IF RUNS OFF END OF FILE
C       MUST LOAD OVERALL HEADER FIRST FOR LUNREDHED (MAY BE MT NOW!)

8888    CALL LUNREDHED(LUN,NX,0,.FALSE.,IRTFLG)
        CALL LUNREDHED(LUN,NX,IMGNUM,.FALSE.,IRTFLG)

        IF (IRTFLG > 0) THEN
C          PROBABLY RAN OFF END OF STACK FILE, ERRT NOT CALLED
           IF (MUSTGET) 
     &       CALL ERRT(102,'THIS IMAGE NOT USED IN STACK',IMGNUM) 
             IRTFLG = 3
           RETURN
        ELSEIF (IRTFLG == 0) THEN
C          NEED IMUSED FROM THIS STACKED IMAGE
           CALL LUNGETINUSE(LUN,IMUSED,IRTFLG)
        ELSE
           IMUSED = 0
        ENDIF
      
        !write(6,*) ' In getoldstack, imgnum 3:',imgnum,irtflg,imused

        IF (IMUSED == 0) THEN
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

        !write(6,*) ' In getoldstack, filnam 1:',filnam(:nlet)
C       APPEND CURRENT STACKED IMAGE NUMBER TO STACK FILE NAME
        
C       (INTTOCHAR ALSO RETURNS NEW VALUE FOR NLET) 
        LENAT = INDEX(FILNAM,'@')
        CALL INTTOCHAR(IMGNUM,FILNAM(LENAT+1:),NLET,0)
        NLET = NLET + LENAT

C       SET FILENAME IN HEADER OBJECT
        CALL LUNSETFILE(LUN,FILNAM(1:NLET),'O',IRTFLG)  

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUN,IMGNUM,NX,IRTFLG)

C       WRITE OUT FILE OPENING INFO 
        CALL LUNSAYINFO(LUN,IRTFLG)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON(LUN,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
 
	END


C       ---------------------- GETOLDSTACK_MRC ------------------------  

	SUBROUTINE GETOLDSTACK_MRC(LUN,IMGNUM,WANTNEXT,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'    ! NEED: MAXNAM

	INTEGER               :: LUN,IMGNUM,IRTFLG
        LOGICAL               :: WANTNEXT

        CHARACTER(LEN=1)      :: DSP
        CHARACTER(LEN=MAXNAM) :: FILNAM,FIL_AT
        CHARACTER(LEN=10)     :: CNUM

	INTEGER               :: IMUSED, NX,NY,NZ,MZ,NSTACK,NLET,IDUM
        INTEGER               :: LOCAT,LOCSLASH,IGO,LENAT
        
        INTEGER               :: lnblnkn   ! FUNCTION

C       ONLY FOR REGULAR MRC IMAGE STACK 

C       DO NOT CALL ERRT IF RUNS OFF END OF FILE

C       SEE IF IMAGNUM IS BEYOND LENGTH OF STACK
        CALL LUNGETNSTACK_MRC(LUN, NX,NY,NZ,MZ,NSTACK,IRTFLG)
        IF (IMGNUM > NSTACK) THEN
            IRTFLG = -2     !  END OF FILE BEFORE IMGNUM
            RETURN
        ENDIF

C       GET FILENAME FROM CURRENT HEADER OBJECT
        CALL LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)

C       CONVERT IMGNUM TO CHARACTER STRING
        CALL INTTOCHAR(IMGNUM,CNUM,NLET,1)

C       APPEND CURRENT STACKED IMAGE NUMBER TO STACK FILE NAME
        LOCAT  = INDEX(FILNAM,'@')
        FIL_AT = CNUM(1:NLET) // '@' // FILNAM(LOCAT+1:)
        FILNAM = FIL_AT
        LOCAT  = INDEX(FILNAM,'@')
        !write(6,*)' In getnewstack_mrc - imgnum,locat:',imgnum,locat

        IF (LOCAT == 0) THEN
C          HANDLES: 'stk'
           CALL ERRT(101,'NO @ IN FILENAME',IDUM)
           IRTFLG = 1
           RETURN

        ELSEIF (LOCAT == 1)  THEN
C          HANDLES: '@stk'
           FILNAM = FIL_AT

        ELSEIF (FILNAM(LOCAT-1:LOCAT-1) == '/')  THEN
C          HANDLES: '/@stk'
           FILNAM = FILNAM(:LOCAT-1) // FIL_AT

        ELSE 
C          HANDLES: 'nnn@stk' AND '/nnn@stk'
           LOCSLASH = INDEX(FIL_AT(1:LOCAT),'/',BACK=.TRUE.)
           IF (LOCSLASH == 0) THEN
              FILNAM = FIL_AT
           ELSE
              FILNAM = FILNAM(:LOCSLASH-1) // FIL_AT
           ENDIF
        ENDIF
        NLET = LNBLNKN(FILNAM)

        !write(6,*) 'In getoldstack - filnam: ',filnam(:nlet)

C       SET FILENAME IN HEADER OBJECT
        CALL LUNSETFILE_MRC(LUN,FILNAM(1:NLET),'O',IRTFLG) 

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETPOS_MRC(LUN,IMGNUM,IRTFLG)

C       WRITE OUT FILE OPENING INFO 
        CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON_MRC(LUN,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
 
	END



