
C++*********************************************************************
C
C GETNEWSTACK.F  NEW                         7 JAN   99   ArDean Leith
C                SET INUSE                     OCT   99   ArDean Leith
C                INDEXED STACK                 JAN 2003   ArDean Leith
C                HEADER COPY                   FEB 2003   ArDean Leith
C                COPYSTAT PARAM                OCT 2010   ArDean Leith
C                MRC SUPPORT                   SEP 2019   ArDean Leith
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
C  GETNEWSTACK(LUN_CP,LUNOUT,COPYSTAT,NX,IMGNUM,IRTFLG)
C
C  PURPOSE:       OPEN A NEW IMAGE WITHIN STACK AND TIE TO LUNOUT
C
C  PARAMETERS:
C      LUN_CP     UNIT TO COPY HEADER VALUES FROM               (SENT)
C      LUNOUT     LOGICAL UNIT NUMBER FOR FILE.                 (SENT)
C      COPYSTAT   COPY FMIN.. FROM INPUT FILE                   (SENT)
C      NX         IMAGE X DIMENSION                             (SENT) 
C      IMGNUM     IMAGE NUMBER WANTED                           (SENT) 
C      IRTFLG     ERROR RETURN FLAG.                            (RET.)
C                 IRTFLG =  1    ERROR, REDHED FAILED
C                 IRTFLG =  0    NORMAL RETURN, IMAGE IS STACKED
C                 IRTFLG = -1    OK, BUT FILE REQUESTED IS NOT STACK
C
C  NOTE:      THE THREE SUBROUTINES: GETSTACK, GETOLDSTACK, AND
C             GETNEWSTACK WERE WRITTEN TO HANDLE WHOLE BARE 
D             STACK (@STACK) OPERATIONS. THEY ARE ONLY USED IN 
D              THE OPERATIONS:
C             'AF' (NO DOCS), 'CP TO STK','DE','CE AD','AR SCA', 
C             'CE' (MANY?),
C             'AD D2', 'AD DF', 'AD DR', 'AD F', 'AD R', 'DIV 2', 
C             'MUL 2', 'MUL O', 'MU O', & 'SU B2'.
C             I DECIDED THAT BARE WHOLE STACK OPERATION SUPPORT 
C             WAS LESS IMPORTANT THAN SUPPORTING OPERATIONS ON 
C             SELECTED STACKED IMAGES SO I STOPPED CONVERTING 
C             OPERATIONS TO HANDLING WHOLE STACKS WITHOUT SELECTED 
C             STACK SUPPORT. THE NEW PARADIGM USES OPFILES AND 
C             NEXTFILE WHICH HANDLES BOTH SELECTED AND WHOLE STACKS 
C             BUT I HAVE NOT IMPLEMENTED IT YET IN MOST OPERATIONS AND
C             AM RUNNING OUT OF TIME TO EVER FINISH IMPLEMENTING 
C             ANYTHING MORE 
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

	SUBROUTINE GETNEWSTACK(LUN_CP,LUNOUT,COPYSTAT,NX,IMGNUM,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'     ! FOR: MAXNAM
       
        INTEGER                :: LUN_CP,LUNOUT,NX,IMGNUM,IRTFLG
        LOGICAL                :: COPYSTAT

        CHARACTER (LEN=MAXNAM) :: FILNAM
        CHARACTER (LEN=1)      :: DSP

        INTEGER                :: NY,NZ,ITYPE,ISTACK,MAXIM,LENAT,NLET
        INTEGER                :: IMAMI
        REAL                   :: FMIN,FMAX,AV,SIG
        INTEGER                :: IRTFLGT
        LOGICAL                :: ISBARET,IS_MRC

C       SE IF WHOLE STACK OPERATION
        CALL LUNGETISBARE(LUNOUT,ISBARET,IRTFLG)    ! MRC OK
        IF (.NOT. ISBARET) THEN
C          INPUT IMAGE IS SIMPLE OR A SPECIFIC STACKED IMAGE
           IRTFLG = -1
           RETURN
        ENDIF

C       DETERMINE IF MRC OR SPIDER STACK FILE
        CALL LUNGETIS_MRC(LUNOUT,IS_MRC,IRTFLG)

        IF (IS_MRC) THEN
C          MRC STACK
	   CALL GETNEWSTACK_MRC(LUNOUT,COPYSTAT,IMGNUM,IRTFLG)
           RETURN
        ENDIF

C       MAKE A NEW IMAGE WITHIN EXISTING INLINE OR REGULAR STACK 

C       GET OVERALL HEADER FROM THE STACK FILE
        CALL LUNREDHED(LUNOUT,NX,0,.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           IRTFLG = 1
           RETURN
        ENDIF

C       RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM OVERALL HEADER
        CALL LUNGETMAXIM(LUNOUT,MAXIM,IRTFLG)

        IF (IMGNUM > MAXIM) THEN
C          MUST UPDATE OVERALL HEADER WITH MAXIMUM IMAGE NUMBER
           CALL LUNSETMAXIM(LUNOUT, IMGNUM,IRTFLG)
           CALL LUNSETMAXALL(LUNOUT,IMGNUM,IRTFLG)
        ENDIF

C       NEED SOME VALUES TO PUT IN STACKED IMAGE HEADER OBJECT
        CALL LUNGETSIZE(LUNOUT,NX,NY,NZ,IRTFLG)
        CALL LUNGETTYPE(LUNOUT,ITYPE,IRTFLGT)
        CALL LUNCOPYSTK(LUNOUT,ISTACK,IRTFLGT)

        IF (ISTACK < 0) THEN
C          MAKING A NEW INDEXED STACKED FILE, UPDATE INDX LOCATION
           CALL LUNWRTINDX(LUNOUT,IMGNUM,NX,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN
        ENDIF

        IF (IMGNUM > MAXIM .OR. ISTACK > 2) THEN
C          SAVE OVERALL HEADER NOW TO PRESERVE MAXIM & LASTINDX
           CALL LUNWRTHED(LUNOUT,NX,0,IRTFLGT)
        ENDIF

C       GET FILENAM FROM CURRENT HEADER OBJECT
        CALL LUNGETFILE(LUNOUT,FILNAM,NLET,DSP,IRTFLG)

C       APPEND CURRENT STACKED IMAGE NUMBER TO STACK FILE NAME
C       (INTTOCHAR ALSO RETURNS NEW VALUE FOR NLET) 
        LENAT = INDEX(FILNAM,'@')
        CALL INTTOCHAR(IMGNUM,FILNAM(LENAT+1:),NLET,0)
        NLET = NLET + LENAT

C       INITIALIZE HEADER OBJECT FOR NEW STACKED IMAGE FROM LUN_CP
C         (I THINK THIS COULD BE MOSTLY SKIPPED??)
        CALL LUNSETHDR(LUN_CP,LUNOUT,NX,NY,NZ,
     &                 ITYPE,ISTACK,IRTFLG)

        IF (COPYSTAT) THEN
C          COPY STATISTICS FROM INPUT FILE
           CALL LUNGETSTAT(LUN_CP,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)
           CALL LUNSETSTAT(LUNOUT,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)
        ENDIF

C       PUT IMAGE INUSE IN HEADER OBJECT
        CALL LUNSETINUSE( LUNOUT,IMGNUM,IRTFLG)
        CALL LUNSETISTACK(LUNOUT,0,IRTFLG)

C       PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
        CALL LUNSETCOMMON(LUNOUT,IRTFLG)

C       PUSH HEADER OBJECT INFO INTO NEW STACKED FILE
        CALL LUNWRTHED(LUNOUT,NX,IMGNUM,IRTFLG)

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUNOUT,IMGNUM,NX,IRTFLGT)

C       WRITE OUT FILE OPENING INFO TO SCREEN
        CALL LUNSAYINFO(LUNOUT,IRTFLG)

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
 	END


C       ---------------------- GETNEWSTACK_MRC ------------------------  

	SUBROUTINE GETNEWSTACK_MRC(LUN,COPYSTAT,IMGNUM,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMLIMIT.INC'    ! NEED: MAXNAM
       
        INTEGER                :: LUN,IMGNUM,IRTFLG
        LOGICAL                :: COPYSTAT     !   UNUSED

        CHARACTER (LEN=MAXNAM) :: FILNAM,FIL_AT
        CHARACTER (LEN=1)      :: DSP
        CHARACTER (LEN=10)     :: CNUM

        INTEGER                :: NX,NY,NZ,ITYPE,ISTACK,NLET,IDUM
        INTEGER                :: MZ,NSTACK,LOCAT,LOCSLASH,IGO,LENU

        INTEGER                :: lnblnkn ! function

C       MAKE A NEW IMAGE WITHIN EXISTING MRC STACK 

C       RETRIEVE CURRENT MAXIMUM IMAGE NUMBER FROM MRC HEADER
        CALL LUNGETNSTACK_MRC(LUN, NX,NY,NZ,MZ,NSTACK,IRTFLG)

        !write(6,*) 'In getnewstack_mrc- nz,mz,nstack:',nz,mz,nstack

        IF (IMGNUM > NSTACK) THEN
C          MUST UPDATE OVERALL HEADER WITH MAXIMUM IMAGE NUMBER
           NSTACK = IMGNUM
           CALL LUNSETNSTACK_MRC(LUN, NZ,NSTACK,IRTFLG)
        ENDIF

C       GET FILENAM FROM HEADER OBJECT FOR EARLIER FILE IN SERIES
        CALL LUNGETFILE_MRC(LUN,FILNAM,NLET,DSP,IRTFLG)

        !write(6,*) 'In getnewstack_mrc - filnam1:',filnam(:nlet)

C       SUBSTITUTE @ IN NAME WITH  CURRENT STACKED IMAGE NUMBER AND @
C       TO STACK FILE NAME

C       CONVERT IMGNUM TO CHARACTER STRING
        CALL INTTOCHAR(IMGNUM,CNUM,NLET,1)

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

        !call subchar(cnum(1:nlet+1),filnam,locslash,lenu,irtflg)
        !write(6,*)' In getnewstack_mrc - nlet,cnum:',nlet,cnum

C       SET FILENAME IN HEADER OBJECT
        CALL LUNSETFILE_MRC(LUN,FILNAM(1:NLET),'N',IRTFLG)  
        
C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETPOS_MRC(LUN,IMGNUM,IRTFLG)

C       PUSH HEADER OBJECT INFO INTO OVERALL STACK HEADER IN FILE
        CALL LUNWRTHED_MRC(LUN,IRTFLG)

C       WRITE OUT FILE OPENING INFO TO SCREEN
        CALL LUNSAYINFO_MRC(LUN,.TRUE.,IRTFLG)

C       PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
        CALL LUNSETCOMMON_MRC(LUN,IRTFLG)

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
 	END





