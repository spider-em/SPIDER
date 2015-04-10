
C++*********************************************************************
C
C OPENINSTK.F   -- NEW NOVEMBER 1996   AUTHOR: ArDean Leith
C                  CHANGED IRTFLG DEFAULT        JAN 99 -- ArDean Leith
C                  REWRITTEN                     JAN 99 -- ArDean Leith 
C                  USED LUNHDR                   FEB 99 -- ArDean Leith 
C                  INLINE VOLUME STACKS          AUG 02 -- ArDean Leith                         
C                  INDEXED STACKS                JAN 03 -- ArDean Leith
C                  HEADER COPY                   FEB 03 -- ArDean Leith
C                  OPENINLY *8                   OCT 10 -- ArDean Leith
C                  ERROR MSG                     FEB 12 -- ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C    OPENINSTK(LUNT,FILNAM,LUN,NX,NY,NZ,NSTACK,ITYPE, DISP,IRTFLG)
C
C    PURPOSE:       TO OPEN A NEW OR EXISTING INLINE STACK. CAN OPEN
C                   OVERALL (BARE) STACK OR SPECIEFIED IMAGE WITHIN
C                   AN INLINE STACK
C
C    PARAMETERS:
C        LUNT       UNIT TO COPY HEADER VALUES FROM              (SENT)
C        FILNAM     CHARACTER ARRAY, CONTAINING FILE NAME        (SENT)
C        LUN        LOGICAL UNIT NUMBER FOR FILNAM I/O.          (SENT)
C        NX,NY      DIMENSIONS OF FILE                      (SENT/RET.)
C        NZ         NUMBER OF PLANES                        (SENT/RET.)
C        NSTACK     STACK INDICATOR                         (SENT/RET.)
C                   ON INPUT:
C                      >0 : REGULAR STACK FILE (IF NEW)
C                      <0 : INDEXED STACK FILE (IF NEW)
C                   ON OUTPUT:                               
C                      -2 : NOT STACK = ERROR
C                      -1 : STACKED IMAGE
C                       0 : REGULAR BARE STACK, CONTAINS NO IMAGE(S)
C                      >0 : INDEXED BARE STACK, VALUE IS MAX. IMAGE
C        ITYPE      IFORM                                   (SENT/RET.)                    
C        DISP       FILE DISPOSITION, SEE OPFIL FOR VALUES       (SENT)
C        IRTFLG     ERROR RETURN FLAG.                      (SENT/RET.)
C                   IRTFLG = 0    NORMAL RETURN
C                   IRTFLG = 1    ERROR RETURN
C
C    CALL TREE:  SEE OPFIL
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
C  2,195,511,375 10

	SUBROUTINE OPENINSTK(LUNT,FILNAM,LUN,NX,NY,NZ,
     &                       NSTACK,ITYPE,DISP,IRTFLG)

C       USE INLINE BUFFER COMMON AREA
        INCLUDE 'INLN_INFO.INC'

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'LABLOCK.INC'

        INTEGER                      :: LUNT
        CHARACTER (LEN=*)            :: FILNAM,DISP
        INTEGER                      :: LUN,NX,NY,NZ,NSTACK,ITYPE,IRTFLG

	LOGICAL                      :: STACKOPN,ISDIGI,CALLERRT

        INTEGER, PARAMETER           :: I_8 = SELECTED_INT_KIND(12)
	INTEGER(KIND=I_8)            :: NWORDS8   ! FROM INLN_INFO
        INTEGER(KIND=I_8), PARAMETER :: ZERO_8 = 0


C       SHOULD NOT STOP IF DISP == 'Z' AND REDHED FAILS
        CALLERRT =  (DISP(1:1) .NE. 'Z') 

C       SET ERROR RETURN
        IRTFLG   = 1

        NSTACKIN = NSTACK

C       MAKE SURE USER WANTS TO USE A INLINE STACK      
        ILOCAT = INDEX(FILNAM,'@')

        IEND = ILOCAT + 1
        IF (ISDIGI(FILNAM(ILOCAT + 2:ILOCAT+2))) IEND = IEND + 1
        IF (ILOCAT .LE. 2 .OR. FILNAM(1:1) .NE. '_') THEN
           WRITE(NOUT,*)'*** BAD SYNTAX FOR INLINE STACK:',FILNAM
           CALL ERRT(100,'OPENINSTK',NE)
           RETURN
        ENDIF

C       RETRIVE INLINE BUFFER NUMBER FROM FILE NAME
        CALL INLNBUF(FILNAM(1:IEND),NLET,INLNED,IRTFLGT)
        IF (IRTFLGT .NE. 0)  RETURN

        IF (ISDIGI(FILNAM(ILOCAT + 1:ILOCAT + 1))) THEN
C          FIND IMAGE NUMBER WITHIN STACK FILE 
           CALL GETFILENUM(FILNAM(ILOCAT:),IMGNUM,NDIGITS,
     &                     .TRUE.,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN
           IF (IMGNUM .LE. 0) THEN
              CALL ERRT(101,'STACKS START WITH IMAGE: 1',NE)
              RETURN
           ENDIF

C          FOR SPECIFIC IMAGE RETURN NSTACK = -1
           NSTACK = -1
        ELSE
C          SET DEFAULT IMGNUM FOR BARESTACK 
           IMGNUM   = 0
C          (BARE STACKS RETURN NSTACK >= 0)
           NSTACK   = 0
        ENDIF

C       SEE IF INLINE STACK EXISTS NOW
        STACKOPN = (NSAMBUF(INLNED) > 0)


C -------------------------------- NEW --------------------------------

10	IF (DISP(1:1)  ==  'U' .OR. DISP(1:1)  ==  'N') THEN
C         WANT TO MAKE A NEW STACK OR NEW IMAGE WITHIN EXISTING STACK
   
          IF (.NOT. STACKOPN) THEN
C            INLINE STACK DOES NOT EXIST YET, CREATE NEW INLINE STACK 

             IF (IMGNUM  ==  0 .AND. STACKOPN) THEN
C	        ACCESS THE WHOLE STACK, NOT A PARTICULAR IMAGE BUT SINCE
C               DISP IS NEW AND STACK IS ALREADY OPEN THIS IS REALLY
C               ASKING TO CREATE A NEW STACK IN AN EXISTING STACK. 

C               CLOSE OLD STACK (THIS DESTROYS DATA IN IT!!)
                CALL OPENINLN(LUN,INLNED,.TRUE.,0,ZERO_8,
     &                       .FALSE.,IRTFLGT)
             ENDIF

             CALL RDPRI1S(NIMAGE,NOT_USED,
     &              'NUMBER OF IMAGES/VOLUMES ALLOWED IN STACK',IRTFLGT)
             IF (IRTFLGT .NE. 0) RETURN

             IF (NSTACKIN .LT. 0) THEN
                 CALL RDPRI1S(ISTACK,NOT_USED,
     &           'HIGHEST IMAGE/VOLUME NUMBER ALLOWED IN STACK',IRTFLGT)
                 IF (IRTFLGT .NE. 0) RETURN
                 IF (ISTACK .LT. 1) THEN
                     CALL ERRT(101,'MIN. NO. FOR INDEXED STACK IS 1',NE)
                     RETURN                        
                  ENDIF
                  ISTACK = - ISTACK
             ELSE
                  ISTACK = 2
             ENDIF

C            CREATE NEW OVERALL HEADER 
             CALL LUNSETHDR(0,LUN,NX,NY,NZ,
     &                      ITYPE,ISTACK,IRTFLGT)
             IF (IRTFLGT .NE. 0) RETURN

C            SAVE ISTACK & MAXIM IN STATIC PART OF HEADER OBJECT
             CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT) ! sets (259) & (260)
             CALL LUNSETMAXALL(LUN,0,IRTFLGT)

C            GET RECORD INFO
             CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,LABBYT,
     &                      LENBYT,IRTFLGT)
             IF (IRTFLGT .NE. 0) RETURN

C            SET UP INLINE BUFFER AND TIE IT TO LUN
C            NWORDS = (LABREC + INDXREC +NIMAGE  * NRECS) * (LENBYT / 4) 
             NWORDS8 = (LABREC + INDXREC +NIMAGE  * NRECS) 
             NWORDS8 = NWORDS8  * (LENBYT / 4) 
 
             CALL OPENINLN(LUN,INLNED,.TRUE.,NX,NWORDS8,
     &                     .TRUE.,IRTFLGT)
             IF (IRTFLGT .NE. 0)  RETURN

C            WRITE OVERALL HEADER RECORD(S) INTO INLINE STACK	    
             CALL LUNWRTHED(LUN,NX,0,IRTFLGT)
             IF (IRTFLGT .NE. 0) RETURN

             IRECBUF(INLNED)   = IREC  
             LABRECBUF(INLNED) = LABREC  

             IF (ISTACK .LT. 0) THEN
C               CLEAR STACK INDEX IN NEW FILE
                CALL LUNCLRINDX(LUN,NX,IRTFLGT)
             ENDIF

             IF (IMGNUM  ==  0) THEN
C               DO NOT HAVE AN IMAGE TO PLACE IN THE STACK YET
C	        ACCESS THE WHOLE STACK NOT A PARTICULAR IMAGE

                CALL LUNSETISBARE(LUN,.TRUE.,IRTFLGT)
                IF (IRTFLGT .NE. 0) RETURN
                GOTO 7777
             ENDIF

          ELSE 
C            PUTTING IMAGE INTO EXISTING INLINE FILE 

C            USE EXISTING INLINE BUFFER, TIE IT TO LUN & GET NX
             CALL OPENINLN(LUN,INLNED,.FALSE.,NX,ZERO_8,
     &                    .TRUE.,IRTFLGT)
             IF (IRTFLGT .NE. 0)  RETURN
          ENDIF
       
C         RECOVER OVERALL HEADER FROM INLINE STACK FILE
          CALL LUNREDHED(LUN,NX,0,.TRUE.,IRTFLGT)
          IF (IRTFLGT .NE. 0) RETURN

          CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG)   !gets (24) & (26)

C         COPY ISTACK & MAXIM INTO STATIC PART OF HEADER OBJECT 
          CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT)     !sets (259) & (260)
          CALL LUNSETMAXALL(LUN,MAXIM,IRTFLGT)

          CALL LUNGETTYPE(LUN,ITYPEF,IRTFLGT)

          CALL LUNGETSIZE(LUN,NXF,NYF,NZF,IRTFLGT)

          CALL LUNGETLAB(LUN,LABRECF,INDXREC,NRECF,NDUM,NDUM,IRTFLGT)

          IF (ISTACK  ==  0) THEN
C             INLINE BUFFER DOES NOT CONTAIN A STACK
              CALL ERRT(101,'INLINE BUFFER DOES NOT CONTAIN A STACK',NE)
              RETURN

          ELSEIF (NXF   .NE. NX .OR. NYF .NE. NY .OR.
     &            NZF .NE. NZ) THEN
C             EXISTING FILE HAS DIFFERING DIMENSIONS
              CALL ERRT(102,
     &                  'IMAGE DIMENSIONS NOT SAME AS IN INLINE STACK',
     &                  INLNED)
              RETURN

          ELSEIF (ITYPE .NE. ITYPEF) THEN
C            EXISTING STACK FILE FORMAT NOT SAME AS IMAGE FORMAT
             WRITE(NOUT,96) ITYPE,ITYPEF
96           FORMAT('IMAGE FORMAT: ',I5,
     &              ' NOT COMPATIBLE WITH EXISTING STACK FORMAT: ',I5)
             CALL ERRT(100,'OPENINSTK',NE)
             RETURN
          ENDIF

          IF (IMGNUM  ==  0) THEN
C	     ACCESS THE WHOLE STACK NOT A PARTICULAR IMAGE
             CALL LUNSETISBARE(LUN,.TRUE.,IRTFLGT)
             IF (IRTFLGT .NE. 0) RETURN
             GOTO 7777
          ENDIF

          IF (IMGNUM .GT. MAXIM) THEN
C            UPDATE OVERALL HEADER WITH MAX. IMAGE NUMBER IN USE
             CALL LUNSETMAXIM(LUN,IMGNUM,IRTFLGT)  ! sets(26)
             IF (IRTFLGT .NE. 0) RETURN
          ENDIF

          IF (ISTACK .LT. 0) THEN
C            MAKING A NEW INDEXED STACKED FILE, SET STORAGE LOCATION
             CALL LUNWRTINDX(LUN,IMGNUM,NX,IRTFLGT)
             IF (IRTFLGT .NE. 0) RETURN
          ENDIF

          IF (IMGNUM .GT. MAXIM .OR. ISTACK .LT. 0) THEN
C            SAVE OVERALL HEADER TO PRESERVE IMGNUM & LASTINDX
             CALL LUNWRTHED(LUN,NX,0,IRTFLGT)
          ENDIF

C         SET MAXIM VALUE IN HEADER OBJECT STATIC AREA
          MAXIM = MAX(IMGNUM,MAXIM)
          CALL LUNSETMAXALL(LUN,MAXIM,IRTFLG)

C         CREATE HEADER FOR NEW STACKED IMAGE NOW
C         KEEPS STATIC HEADER SETTINGS
          ISTACK = 0
          CALL LUNSETHDR(LUNT,LUN,NX,NY,NZ,ITYPE,ISTACK,IRTFLGT)
          IF (IRTFLGT .NE. 0) RETURN

          CALL LUNSETMAXIM(LUN,0,IRTFLGT)

C         SET IMUSED FLAG FOR STACKED IMAGE HEADER
          CALL LUNSETINUSE(LUN,IMGNUM,IRTFLGT)
          IF (IRTFLGT .NE. 0) RETURN

C         PLACE NEW IMAGE HEADER INTO PROPER STACK LOCATION
          CALL LUNWRTHED(LUN,NX,IMGNUM,IRTFLGT)
          IF (IRTFLGT .NE. 0) RETURN

C -------------------------------- OLD --------------------------------
           
	ELSEIF (DISP(1:1)  ==  'O' .OR. DISP(1:1)  ==  'K' .OR.
     &          DISP(1:1)  ==  'Z' .OR. 
     &          DISP(1:1)  ==  'E' .OR. DISP(1:1)  ==  'M') THEN
C          WANT AN EXISTING IMAGE FROM EXISTING STACK OR AN
C          EXISTING BARE STACK HEADER

           IF (.NOT. STACKOPN) THEN
C             INLINE STACK DOES NOT EXIST YET, ERROR
              WRITE(NOUT,*) '*** INLINE STACK DOES NOT EXIST'
C	      FOR DISP=Z, DO NOT STOP BATCH JOBS BY CALLING ERRT
              IF (CALLERRT) CALL ERRT(100,'OPENINSTK',NE)
              RETURN
           ENDIF

C          USE EXISTING INLINE BUFFER, TIE IT TO LUN & GET NX
           CALL OPENINLN(LUN,INLNED,.FALSE.,NX,ZERO_8,.TRUE.,IRTFLGT)
           IF (IRTFLGT .NE. 0)  RETURN

C          GET OVERALL HEADER 
           CALL LUNREDHED(LUN,NX,0,CALLERRT,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN
          
C          RECOVER MAXIM FROM OVERALL HEADER & PLACE IN STATIC AREA 
           CALL LUNCOPYMAXIM(LUN,MAXIM,IRTFLGT)

C          RECOVER ISTACK FROM OVERALL HEADER & PLACE IN STATIC AREA 
           CALL LUNCOPYSTK(LUN,ISTACK,IRTFLGT)

           IF (IMGNUM .NE. 0) THEN
C             GET SPECIFIED IMAGE HEADER 
              CALL LUNREDHED(LUN,NX,IMGNUM,CALLERRT,IRTFLGT)
              IF (IRTFLGT .NE. 0) THEN
                 IF (CALLERRT) THEN
                    CALL ERRT(102,'STACK LACKS IMAGE NUMBER',IMGNUM)
                 ENDIF
                 RETURN
              ENDIF
           ENDIF

C          RECOVER IMAGE PARAMETERS FROM STACKED IMAGE HEADER
           CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLGT)

           IF (IMGNUM  ==  0) THEN
C             FOR OVERALL STACK RETURN NSTACK = MAX IMAGE (BUF(26))
              NSTACK = MAXIM
           ELSE

C             SEE IF THIS IMAGE IS USED IN THE STACK
              CALL LUNGETINUSE(LUN,IMGNUMT,IRTFLGT)

              IF (IMGNUM .NE. IMGNUMT) THEN
C                NO EXISTING IMAGE WITHIN STACK, (THIS IMAGE UNUSED)
                 IF (CALLERRT) THEN
                    CALL ERRT(102,'THIS IMAGE NOT USED IN STACK',IMGNUM)
                 ENDIF
                 RETURN
              ENDIF 
           ENDIF
           CALL LUNGETTYPE(LUN,ITYPE,IRTFLG)

        ELSE
          WRITE(NOUT,*)'*** PGM. ERROR: UNKNOWN DISP IN OPENINSTK:',DISP
          CALL ERRT(100,'OPENINSTK',NE)
          RETURN
        ENDIF

C ------------------------------- BOTH --------------------------------

C       SET OFFSETS FOR REDLIN/WRTLIN ON THIS LUN
        CALL LUNSETIMGOFF(LUN,IMGNUM,NX,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       WRITE OUT FILE OPENING INFO
7777    CALL LUNSAYINFO(LUN,IRTFLGT)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON(LUN,IRTFLGT)

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0
        RETURN

	END

