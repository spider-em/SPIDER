
C++*********************************************************************
C
C OPENSTK.F -- CREATED                         DEC 96 -- ArDean Leith
C              USED LUNHDR                     FEB 99 -- ArDean Leith
C              INDEXED STACKS                  JAN 03 -- ArDean Leith
C              HEADER COPY                     FEB 03 -- ArDean Leith
C              OPENFIL PARAMETERS              APR 04 -- ArDean Leith
C              BAD IRTFLG RETURN               AUG 04 -- ArDean Leith
C              ERROR MSG                       DEC 10 -- ArDean Leith
C              MPI ERROR MSG                   MAR 11 -- ArDean Leith
C              MSG                             FEB 12 -- ArDean Leith
C              'ST' MISSING STACKED IMAGE      AUG 14 -- ArDean Leith
C              'FS' ON MISSING STACKED IMAGE   SEP 16 -- ArDean Leith
C              NSTK_FLG -> NSTK_FLG            APR 26 -- ArDean Leith
C              DEBUG OUTPUT ADDED              APR 26 -- ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2026  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email:                                                             *
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
C  OPENSTK(LUNT,FILNAM,LUN,NX,NY,NZ,NSTK_FLG,ITYPE,DISP,IRTFLG)
C
C  PURPOSE:       TO OPEN A NEW OR EXISTING STACK FILE.  
C                 NOT FOR INLINE STACKS
C
C  PARAMETERS:
C        LUNT       UNIT TO COPY HEADER VALUES FROM              (SENT)
C        FILNAM     CHARACTER ARRAY, CONTAINING FILE NAME        (SENT)
C        LUN        LOGICAL UNIT NUMBER FOR FILNAM               (SENT)
C        NX,NY      DIMENSIONS OF FILE                       (SENT/RET)
C        NZ         NUMBER OF PLANES                         (SENT/RET)
C    
C        NSTK_FLG   STACK INDICATOR                          (SENT/RET)
C                   ON INPUT:
C                      >0 : REGULAR STACK FILE (IF NEW)
C                      <0 : INDEXED STACK FILE (IF NEW)
C                   ON OUTPUT:                               
C                      -2 : NOT STACK = ERROR
C                      -1 : STACKED IMAGE
C                       0 : REGULAR BARE STACK, CONTAINS NO IMAGE(S)
C                      >0 : INDEXED BARE STACK, VALUE IS MAX. IMAGE
C                       5 : NOT SPIDER FILE?
C
C        ITYPE      IFORM                                    (SENT/RET)
C        DISP       FILE DISPOSITION, SEE OPFIL FOR VALUES       (SENT)
C        IRTFLG     ERROR RETURN FLAG.                            (RET)
C                   IRTFLG = 0    NORMAL RETURN
C
C  CALL TREE:  SEE OPFILC
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE OPENSTK(LUNT,FILNAM,LUN,NX,NY,NZ,
     &                     NSTK_FLG,ITYPE,DISP,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER                  :: LUNT
        CHARACTER (LEN=*)        :: FILNAM
        INTEGER                  :: LUN,NX,NY,NZ,NSTK_FLG,ITYPE,IRTFLG
        CHARACTER (LEN=*)        :: DISP

        CHARACTER (LEN=MAXNAM)   :: FILNOAT,FILNPE
        CHARACTER (LEN=2*MAXNAM) :: MSG
	LOGICAL                  :: EX,ISDIGI,CALLERRTRED,INDXD
        INTEGER                  :: IMUSED,IMGNUM,NSTK_FLG_IN
        INTEGER                  :: NLET,NLETE,NLEN
        INTEGER                  :: ILOCAT,NDIGITS,IRTFLGT
        INTEGER                  :: ITYPEIN
        INTEGER                  :: NDUM1,NDUM2,NDUM3,NDUM4,IVAL,MAXIM
        INTEGER                  :: ISTACK,NXF,NYF,NZF,IER,NE,NOT_USED
        INTEGER                  :: NRECS
        INTEGER                  :: ICOMM,MYPID,MPIERR

        INTEGER                  :: lnblnkn


#ifdef USE_MPI
        include 'mpif.h'
#endif

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       SHOULD NOT STOP IF DISP == 'Z' AND REDHED FAILS
        CALLERRTRED = (DISP(1:1) .NE. 'Z')

C       SET ERROR RETURN
        IRTFLG      = 1
        NSTK_FLG_IN = NSTK_FLG

        ILOCAT = INDEX(FILNAM,'@')      

        IF (ISDIGI(FILNAM(ILOCAT + 1:ILOCAT + 1))) THEN
C          FIND IMAGE NUMBER WITHIN STACK FILE 
           CALL GETFILENUM(FILNAM(ILOCAT:),IMGNUM,NDIGITS,.TRUE.,IER)
           IF (IER .NE. 0) RETURN

           IF (IMGNUM <= 0) THEN
              NLEN = lnblnkn(FILNAM)
              WRITE(NOUT,*) 'STACK NAME:',FILNAM(1:NLEN)
              CALL ERRT(102,'STACKS START WITH IMAGE: 1 NOT',IMGNUM)
              RETURN
           ENDIF
        ELSE
C          SET IMGNUM FOR BARESTACK
           IMGNUM = 0
        ENDIF

C       GET FILENAME WITHOUT @ AND DATEXC
        FILNOAT = FILNAM(1:ILOCAT-1) // CHAR(0) 

C       CREATE STACK FILE NAME WITHOUT '@' BUT WITH EXTENSION
        CALL FILNAMANDEXT(FILNOAT,DATEXC,FILNPE,NLET,.TRUE.,IRTFLGT)
	IF (IRTFLGT .NE.0) RETURN

C       SEE IF STACK FILE ALREADY EXISTS NOW

#ifdef USE_MPI
        INQUIRE(FILE=FILNPE,IOSTAT=IER,EXIST=EX)
        CALL MPI_BCAST(IER,1,MPI_INTEGER,0,ICOMM,MPIERR)
        IF (MPIERR .NE. 0) THEN
           WRITE(0,*) ' OPENSTK: FAILED TO BCAST IER'
           STOP
        ENDIF

        CALL MPI_BCAST(EX,1,MPI_LOGICAL,0,ICOMM,MPIERR)
        IF (MPIERR .NE. 0) THEN
           WRITE(0,*) ' OPENSTK: FAILED TO BCAST EX'
           STOP
        ENDIF
#else
        INQUIRE(FILE=FILNPE,IOSTAT=IER,EXIST=EX)
#endif
        IF (IER .NE. 0) THEN
           WRITE(NOUT,*) '*** FILE INQUIRY ERROR: ',FILNPE(1:NLET)
           CALL ERRT(100,'OPENSTK',NE)
           RETURN
        ENDIF


C       !!!!!!!!!!!!!!!!!!!!!!!!!!!
cc#define SP_DBUGIO
#if defined (SP_DBUGIO)  
        write(3,*)'  '
        write(3,*)' In openstk  11 ; filnpe:      ', trim(filnpe)
        write(3,*)' In openstk  11 ; ier:         ', ier
        write(3,*)' In openstk  11 ; ex:                 ', ex
        write(3,*)' In openstk  11 ; imgnum :     ', imgnum
        write(3,*)' In openstk  11 ; nstk_flg_in: ', nstk_flg_in
        write(3,*)' In openstk  11 ; disp:              |',disp,'|'
#endif

 
	IF (DISP(1:1) == 'U' .OR. DISP(1:1) == 'N') THEN
C          WANT TO MAKE A NEW STACK OR NEW IMAGE WITHIN EXISTING STACK

C ------------------------ NEW FILE OR IMAGE IN STACK ----------------

           IF (.NOT. EX .OR. IMGNUM  ==  0) THEN
C             STACK FILE DOES NOT EXIST YET, OR MUST BE REPLACED

              IF (NSTK_FLG_IN < 0) THEN
C                FLAG FOR INDEXED STACK

                 CALL RDPRI1S(NSTK_FLG,NOT_USED,
     &           'HIGHEST IMAGE/VOLUME NUMBER ALLOWED IN STACK',IRTFLGT)
                 IF (IRTFLGT .NE. 0) RETURN

                 IF (NSTK_FLG < 1) THEN
                     CALL ERRT(101,'HIGHEST NUMBER MUST BE > 0',NE)
                     RETURN                        
                 ENDIF
                 NSTK_FLG = -NSTK_FLG
              ELSE
C                REGULAR NEW STACK
                 NSTK_FLG = NSTK_FLG_IN
              ENDIF

#if defined (SP_DBUGIO)
              write(3,*)' '
              write(3,*)' In openstk  12 ; imgnum:      ', imgnum
              write(3,*)' In openstk  12 ; nstk_flg_in: ', nstk_flg_in
              write(3,*)' In openstk  12 ; nstk_flg:    ', nstk_flg
              write(3,*)' In openstk  12 ; disp:            |', disp,'|'
              write(3,*)' In openstk  12 ; Calling openfil now'
#endif

C             CREATE NEW STACK FILE, OPENFIL WILL RETURN NSTK_FLG =0
	      CALL OPENFIL(0,FILNOAT,LUN, NX,NY,NZ,NSTK_FLG,
     &                     ITYPE,DISP,.FALSE.,IRTFLG)

#if defined (SP_DBUGIO)
              write(3,*)' '
              write(3,*)' In openstk 12b ; nstk_flg_in: ', nstk_flg_in
              write(3,*)' In openstk 12b ; nstk_flg:    ', nstk_flg
              write(3,*)' In openstk 12b ; lun:         ', lun
              write(3,*)' In openstk 12b ; imgnum:      ', imgnum
              write(3,*)' In openstk 12b ; irtflg:      ', irtflg
#endif

              IF (IRTFLG .NE. 0) RETURN

              IF (NSTK_FLG_IN < 0) THEN
C                CLEAR STACK INDEX IN NEW FILE
                 CALL LUNCLRINDX(LUN,NX,IRTFLGT)
              ENDIF

              IF (IMGNUM <= 0) THEN
C                ONLY WANT TO OPEN NEW BARE STACK
                 IRTFLG  = 0
                 RETURN
              ENDIF

#if defined (SP_DBUGIO)
              write(3,*)' '
              write(3,*)' In openstk 12c ; nstk_flg_in: ', nstk_flg_in
              write(3,*)' In openstk 12c ; lun:         ', lun
              write(3,*)' In openstk 12c ; imgnum:      ', imgnum
              write(3,*)' In openstk 12c ; irtflg:      ', irtflg
#endif


           ELSE
C             OPEN EXISTING STACK FILE TO APPEND A NEW STACKED IMAGE

#if defined (SP_DBUGIO)
              write(3,*)' In openstk  13; Calling openfil'
              write(3,*)' '
#endif
              ITYPEIN = ITYPE
	      CALL OPENFIL(0,FILNOAT,LUN, NXF,NYF,NZF,NSTK_FLG,
     &                     ITYPE,'O',.FALSE.,IRTFLGT)
              IF (IRTFLGT .NE. 0)  RETURN

C             OPENFIL WILL RETURN NUMBER OF IMAGES IN STACK, OR -1
C             IF THIS IS A SPECIFIC IMAGE WITHIN THE STACK, -2 IS
C             FOR NON-STACK IMAGE.

#if defined (SP_DBUGIO)
              write(3,*)' '
              write(3,*)' In openstk  14; nstk_flg: ',nstk_flg
#endif

              IF (NSTK_FLG <= -2) THEN
C                EXISTING FILE IS NOT A STACK
                 CALL ERRT(101,'EXISTING FILE IS NOT A STACK',NE)
                 RETURN

              ELSEIF (NXF .NE. NX .OR. NYF .NE. NY .OR.
     &                NZF .NE. NZ) THEN
C                EXISTING FILE HAS DIFFERING DIMENSIONS
                 MSG = 'IMAGE DIMENSIONS NOT SAME AS IN STACK '//
     &                  FILNOAT(1:NLET) 

                 CALL ERRT(101,MSG,NE)
                 RETURN

              ELSEIF (ITYPEIN .NE. ITYPE) THEN
C                EXISTING STACK FILE FORMAT NOT SAME AS IMAGE FORMAT
                 CALL ERRT(101,
     &              'IMAGES IN STACK MUST HAVE SAME FILE FORMAT',NE)
                 RETURN
              ENDIF
           ENDIF


#if defined (SP_DBUGIO)
           write(3,*)' '
           write(3,*)' In openstk 15; Calling lungetstk 4444444444444'
           write(3,*)'  '
#endif

C          RECOVER MAXIM & ISTACK FROM OVERALL HEADER 
           CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN

#if defined (SP_DBUGIO)
           write(3,*)' In openstk 16; imgnum:   ', imgnum
           write(3,*)' In openstk 16; istack:   ', istack
           write(3,*)' In openstk 16; maxim:    ', maxim
           write(3,*)'  '
#endif

           IF (IMGNUM > MAXIM) THEN
C             UPDATE OVERALL HEADER WITH MAXIMUM IMAGE NUMBER IN USE NOW
              CALL LUNSETMAXIM(LUN,IMGNUM,IRTFLGT)
              CALL LUNSETMAXALL(LUN,IMGNUM,IRTFLGT)
           ENDIF

           IF (ISTACK < 0) THEN
C             NEW INDEXED STACKED FILE, UPDATE INDX LOCATION
              CALL LUNWRTINDX(LUN,IMGNUM,NX,IRTFLGT)
              IF (IRTFLGT .NE. 0) RETURN
           ENDIF

           IF (IMGNUM > MAXIM .OR. ISTACK < 0) THEN
C             SAVE OVERALL HEADER NOW TO PRESERVE MAXIM & LASTINDX
              CALL LUNWRTHED(LUN,NX,0,IRTFLGT)
           ENDIF



C          CREATE HEADER FOR NEW STACKED IMAGE,
C          KEEPS STATIC ISBARE SETTING, MAXIM, AND STKALL SETTING
           ISTACK = 0
           CALL LUNSETHDR(LUNT,LUN,NX,NY,NZ,
     &                    ITYPE,ISTACK,IRTFLGT)

C          SET IMGNUM FOR THIS CURRENT IMAGE
           CALL LUNSETINUSE(LUN,IMGNUM,IRTFLGT)

#if defined (SP_DBUGIO)
           write(3,*)' In openstk 17; imgnum:   ', imgnum
           write(3,*)' In openstk 17; istack:   ', istack
           write(3,*)' In openstk 17; maxim:    ', maxim
           write(3,*)' In openstk 17; Writing new lunwrthed'
           write(3,*)'  '
#endif


C          PLACE NEW STACKED IMAGE HEADER INTO PROPER STACK LOCATION
           CALL LUNWRTHED(LUN,NX,IMGNUM,IRTFLGT)

C          SET PROPER OFFSET INTO LUNSTK FOR IMGNUM
           CALL LUNSETIMGOFF(LUN,IMGNUM,NX,IRTFLGT)

C          RETURNS NSTK_FLG = -1 TO SIGNIFY THIS IS STACKED IMAGE
           NSTK_FLG  = -1

#if defined (SP_DBUGIO)
           write(3,*)' In openstk 17; lunt:      ', lunt
           write(3,*)' In openstk 17; lun:       ', lun
           write(3,*)' In openstk 17; imgnum:    ', imgnum
           write(3,*)' In openstk 17; maxim:     ', maxim
           write(3,*)' In openstk 17; istack:    ', istack
           write(3,*)' In openstk 17; nstk_flg:  ', nstk_flg
           write(3,*)'  '
#endif

C -------------------------------- OLD --------------------------------
           
	ELSEIF (DISP(1:1) == 'O' .OR. DISP(1:1) ==  'K' .OR.
     &          DISP(1:1) == 'Z' .OR. 
     &          DISP(1:1) == 'E' .OR. DISP(1:1) ==  'M') THEN
C          WANT AN EXISTING IMAGE FROM EXISTING STACK OR AN
C          EXISTING BARE STACK HEADER

           NLETE = lnblnkn(FILNOAT)
           IF (.NOT. EX) THEN
C             STACK FILE DOES NOT EXIST YET, ERROR
              WRITE(NOUT,*)'*** STACK FILE NOT FOUND: ',FILNOAT(:NLETE)
C	      FOR DISP=Z, DO NOT STOP THE BATCH JOB BY CALLING ERRT
              IF (DISP .NE. 'Z') CALL ERRT(100,'OPENSTK',NE)
              RETURN
           ENDIF

C          OPEN EXISTING OVERALL STACK FILE, RETURNS MAXIM IN NSTK_FLG

#if defined (SP_DBUGIO)
           write(3,*)'    '
           write(3,*)' In openstk  100; nstk_flg: ', nstk_flg
           write(3,*)' In openstk  100; filnoat:  ', filnoat(:nlete)
           write(3,*)' In openstk  100; Calling openfil'
           write(3,*)'    '
#endif

	   CALL OPENFIL(0,FILNOAT,LUN, NX,NY,NZ,NSTK_FLG,
     &                  ITYPE,'O',.FALSE.,IRTFLGT)

#if defined (SP_DBUGIO)
           write(3,*)'    '
           write(3,*)' In openstk  101; After openfil'
           write(3,*)' In openstk  101; nstk_flg:  ',nstk_flg
           write(3,*)' In openstk  101; imgnum:    ',imgnum
           write(3,*)' In openstk  101; irtflgt:   ',irtflgt
#endif
           IF (IRTFLGT .NE. 0)  RETURN

           IF (NSTK_FLG <= -2) THEN
C             EXISTING FILE IS NOT A STACK FILE
              CALL ERRT(101,'EXISTING FILE IS NOT A STACK',NE)
              RETURN

           ELSEIF (IMGNUM <= 0) THEN
C             JUST WANT BARE STACK, RETURN NSTK_FLG = MAX IMAGE IN STACK 
              IRTFLG = 0
              RETURN

           ELSEIF (IMGNUM > NSTK_FLG) THEN
C             STOP IF REQUESTED IMAGE NOT IN STACK 
              IF (DISP .NE. 'Z') THEN
                 CALL ERRT(102,'THIS IMAGE NOT USED IN STACK',IMGNUM)
              ENDIF
              RETURN
           ENDIF


C          SET OFFSET INTO LUNSTK FOR THIS STACKED IMAGE

           CALL LUNSETIMGOFF(LUN,IMGNUM,NX,IRTFLGT)
           !write(3,*)' In openstk; lun,imgnum,nx: ',lun,imgnum,nx 
           !write(3,*)' In openstk; irtflgt: ',  irtflgt
           IF (IRTFLGT .NE. 0) RETURN

C          GET SPECIFIED IMAGE HEADER FROM STACK FILE LOCATION
           CALL LUNREDHED(LUN,NX,IMGNUM,CALLERRTRED,IRTFLGT)

           !write(3,*)' In openstk; lun,imgnum:     ',lun,imgnum
           !write(3,*)' In openstk; callerrtred+++: ',callerrtred

           IF (IRTFLGT .NE. 0) RETURN

C          RECOVER IMAGE PARAMETERS FROM SPECIFIC IMAGE HEADER

C          GET IMUSED FOR THIS CURRENT IMAGE
           IMUSED = IMGNUM
           CALL LUNGETINUSE(LUN,IMUSED,IRTFLGT)

           !write(3,*)' In openstk; lun,imgnum: ', lun,imgnum
           !write(3,*)' In openstk; imused:    ', imused

           IF (IMUSED .NE. IMGNUM) THEN
C             NO EXISTING IMAGE WITHIN STACK??

              IF (IMUSED == 0) THEN
C                SOME VERY OLD STACKS DID NOT HAVE IMGNUM IN THEM
                 CALL LUNGET25(LUN,IVAL,IRTFLGT)

                 !write(3,*)' After lunget25; lun,ival,irtflgt:', 
                 ! lun,ival,irtflgt

                 IF (IVAL .NE. 1) THEN
C                   STACK MAY LACK IMAGE OR HAS BAD EM2EM HEADER

C                   GET RECORD INFO (CAN BE FROM OVERALL HEADER)
                    CALL LUNGETLAB(LUN,NDUM1,NDUM2,NRECS,
     &                                 NDUM3,NDUM4,IRTFLGT)

                    !write(3,*)' In openstk, imgnum,nrecs,disp:'
                    ! ,imgnum,nrecs,disp(1:1)

                    IF (IMGNUM > 0 .AND. NRECS <= 0 .AND.
     &                   DISP == 'Z') THEN
C                      MAY BE BAD IREC IN STACKS CREATED WITH EM2EM!
c                      WRITE(NOUT,'(/,A,A)')
c    &                 ' *** IF IMAGE FROM EM2EM',
c    &                 ' USE OPERATION: <ST EM2> TO FIX HEADER.'

                       CALL LUNSETINUSE(LUN,IMGNUM,IRTFLGT)
                       IRTFLG = 5 
                       RETURN

                    ELSEIF (DISP == 'Z') THEN
C                      DO NOT STOP ON ERROR
                       WRITE(NOUT,'(A,I0)') ' *** MISSING IMAGE',IMGNUM
                       IRTFLG = 5 
                       RETURN

                    ELSE
                       CALL ERRT(102,'MISSING IMAGE',IMGNUM)
                       IRTFLG = 1
                       RETURN 
                   ENDIF
                   GOTO 999

                 ENDIF

                 IMUSED = IMGNUM
                 CALL LUNSETINUSE(LUN,IMUSED,IRTFLGT)
              ENDIF
           ENDIF

C          RETURN NSTK_FLG = -1 (FOR STACKED IMAGE)
           NSTK_FLG = -1

        ELSE
           CALL ERRT(101,'PGM. ERROR: UNKNOWN DISP IN OPENSTK',NE)
           RETURN
        ENDIF

C ------------------------------- BOTH --------------------------------

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0

999     CONTINUE
     
C       WRITE OUT FILE OPENING INFO
        CALL LUNSAYINFO(LUN,IRTFLGT)

C       SET COMMON BLOCK VARIABLES
        CALL LUNSETCOMMON(LUN,IRTFLGT)

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0

#if defined (SP_DBUGIO)
           write(3,*)'    '
           write(3,*)' In openstk  999; RETURNING -----------'
           write(3,*)' In openstk  999; nstk_flg:  ',nstk_flg
           write(3,*)' In openstk  999; imgnum:    ',imgnum
           write(3,*)' In openstk  999; irtflg:    ',irtflg
#endif

	END


