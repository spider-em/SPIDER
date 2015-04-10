
C++*********************************************************************
C
C REG_SET.F                                          AUTHOR: A. LEITH
C                    REMOVED FILNAMSUB CALL IN READPQ APR 01 A. Leith
C                    ADDED REGLPIPE                   JUL 01 A. Leith
C                    READPQ MOD                       NOV 05 A. Leith
C                    REG_SET_BANKED                   NOV 05 A. Leith
C                    REWRITE                          NOV 05 A. Leith
C                    CHANGED TO < TAGEND              JAN 06 A. Leith
C                    DECREASED MAXRSTRQ               JAN 06 A. Leith
C                    ERRT USAGE                       FEB 06 A. Leith
C                    REG_GET_NAME BUG                 MAR 06 A. Leith
C                    RECREATE OLD BUGGY BEHAVIOUR     APR 06 A. Leith
C                    REMOVED REDUNDANT ERROR MSG      JUN 09 A. Leith
C                    GLO VAR SET BUG                  DEC 09 A. Leith
C                    REG_FIND_IREG                    JAN 10 A. Leith
C                    REG_GET_SEL FOR GLO VAR SEND     MAY 14 A. Leith
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
C    CONTAINS SUBROUTINES FOR CREATING REGISTER BANKS,
C    QUERYING REGISTERS AND SETTING REGISTERS
C
C    REG_INIT(IBANK,IRTFLG)
C    PURPOSE: INITIALIZES A BANK OF REGISTERS   
C
C    REG_FIND(IBANKT,NAME,VALUE,IREG,IRTFLG)
C    PURPOSE: FINDS IF VARIABLE EXISTS, RETURNS VALUE & IREG  
C
C   REG_FIND_IREG(TYPE,STRING,ISGLOBAL,IREG,IRTFLG)
C    PURPOSE: FINDS IF VARIABLE EXISTS, RETURNS OR CREATES VARIABLE  
C
C    REG_NEW(IBANKT,NAME,VALUE,IREG,IRTFLG)
C    PURPOSE: CREATES A NEW REGISTER VAR & VALUE   
C
C    REG_SET_VAR(IBANK,STRING,CREATE,VALUE,IREG,IRTFLG)
C    PURPOSE: SETS FIRST VARIABLE FROM STRING
C
C    REG_SET_BYNUM(IREG,VALUE,IRTFLG)
C    PURPOSE: SETS REGISTER IREG=VALUE
C
C    REG_SET(IXREG,CXREG,VALUE,IRTFLG)
C    PURPOSE: SETS X REGISTER VALUE FOR: IXREG OR CXREG
C
C    REG_GET_BYNUM(IREG,VALUE,IRTFLG)
C    PURPOSE: GETS REGISTER VALUE FOR: IREG
C
C    REG_GET_VAR(IBANK,STRING,CREATE,VALUE,IREG,IEND,IRTFLG)
C    PURPOSE: GETS FIRST VARIABLE FROM STRING
C
C    REG_GET(IBANK,IXREG,CXREG,VALUE,IREGRET,IRTFLG)
C    PURPOSE: GETS A CURRENT REGISTER VALUE FROM IXREG OR CXREG INPUT   
C
C    REG_GET_SEL(IBANK,STRING,CREATE,WANTGLO,NREG,IRTFLG)
C    PURPOSE: PARSES REGISTER LIST INTO NSELREG
C
C    REG_GET_SELS(ISELS,NREG,IRTFLG)
C    PURPOSE: RETURNS REGISTER NUMBERS (NOT CONTENTS) FROM  NSELREG.
C
C    REG_LIST_COPY()
C    PURPOSE: COPIES LISTIN REGISTER VALUES TO LISTOUT REGISTERS
C
C    REG_SET_NSEL(IGO,VAL0,VAL1,VAL2,VAL3,VAL4,IRTFLG)
C    PURPOSE: SETS A REGISTER SPECIFIED IN NSELREG(NVAL) TO VALUE   
C
C    REG_SET_NSELA(NREG,VALUES,FILLALL,IRTFLG)
C    PURPOSE: SETS REGISTERS SPECIFIED IN NSEL TO VALUES   
C
C    REG_GET_NAME(IPOS,NAME,NLET,IRTFLG)
C    PURPOSE: REVERSE LOOKUP OF REGISTER(S) SPECIFIED BY REGVALUES  
C
C    REG_GET_NUMS(IREGS)
C    PURPOSE: GETS TOTAL CURRENT NUMBER OF REGISTERS 
C
C    REG_OPENPIPE(CXNUM,IRTFLG)
C    PURPOSE: OPENS PIPE FOR REGISTERS  
C
C    REG_PIPE(NAME,IRTFLG)
C    PURPOSE: SENDS REGISTER VALUE DOWN LUNREGPIPE   
C
C    REG_REINIT()
C    PURPOSE: RESIZES REGISTER SPACE   
C 
C     REGVALUES   CONTAINS CONTENTS OF REGISTERS X0.....
C     NSELREG     CONTAINS THE NUMBERS OF REGISTERS SELECTED IN THE 
C                 OPERATION LINE. PK X12,X20 WOULD RETURN 13,21 IN 
C                 NSELREG(1) AND NSELREG(2)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      MODULE REG_STUFF

         SAVE

         INTEGER, PARAMETER :: NUMREGLOOP = 26

         LOGICAL            :: REGPIPE  = .FALSE.
         INTEGER, PARAMETER :: LUNREGPIPE  = 302

C        DANGER NPARG=MAXNSEL IS ALSO SET IN spider.f (NPARG) AND readrq!! 
         INTEGER, PARAMETER          :: MAXNSEL = 24  ! REGISTER LIST 
         INTEGER, DIMENSION(MAXNSEL) :: NSELREG 
         INTEGER                     :: NSEL_USED = 0

C        DANGER MAXPRC IS ALSO SET IN spider.f!!
         INTEGER, PARAMETER :: MAXPRC = 20      ! NO. OF LEVELS

         INTEGER, DIMENSION(MAXPRC) :: IGORSTRQ,   IENDRSTRQ
         INTEGER, DIMENSION(MAXPRC) :: IGOREGNUM,  IENDREGNUM
         INTEGER                    :: IMAXREGNUM1,IMAXRSTRQ1

C        THIS SHOULD BE RE-DONE WITH ALLOCATABLE CHAR. ARRAY??
C        GLOBAL (BANK 1) AND OTHER NAMESPACES ARE CONCATENATED!
         INTEGER, PARAMETER      :: MAXRSTRQG = 1600
         INTEGER, PARAMETER      :: MAXRSTRQ  = 16000
         CHARACTER(LEN=MAXRSTRQ) :: RSTRQ = ' '

C        GLOBAL (BANK 1) AND OTHER REGISTERS ARE ALSO CONCATENATED!
         INTEGER, PARAMETER      :: NUMREGISG_ORIG = 1000

         INTEGER, PARAMETER      :: NUMREGIS_ORIG  = 16000
         INTEGER                 :: NUMREGIS       = 0

         DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: REGVALUES 

      END MODULE REG_STUFF


C++*********************************************************************
C
C REG_INIT                    NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_INIT(IBANK,IRTFLG)
C
C    PURPOSE:     INITIALIZES  BANK OF REGISTERS   
C
C    PARAMETERS:  IBANK     BANK OF REGISTERS                   (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_INIT(IBANK,IRTFLG)

      USE REG_STUFF

      IF (IBANK > MAXPRC) THEN
         IT     = MAXPRC
         CALL ERRT(102,'MAXPRC OVERFLOW',IT)
         IRTFLG = 1
         RETURN
      ENDIF

      IF (NUMREGIS <= 0) THEN
C        CREATE THE REGISTER STORAGE ARRAY (ONLY OCCURS ONCE)
         ALLOCATE (REGVALUES(NUMREGIS_ORIG), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'UNABLE TO ALLOCATE REGISTERS:',NUMREGIS_ORIG)
            RETURN
         ENDIF
         NUMREGIS  = NUMREGIS_ORIG

C        SET FLAG FOR MT BANKS (ARRAY OPS!)
         IGOREGNUM   = 0
         IENDREGNUM  = 0

      ENDIF

      IF (IBANK < 1) THEN
         CALL ERRT(102,'ILLEGAL VARIABLE BANK:',IBANK)
         IRTFLG = 1
         RETURN

       ELSEIF (IBANK == 1) THEN
C        INITIALIZE BANK ONE (WHICH IS ALSO THE GLOBAL BANK)
         IGORSTRQ(IBANK)   = 1
         IENDRSTRQ(IBANK)  = 0
         IGOREGNUM(IBANK)  = 1
         IENDREGNUM(IBANK) = 0

         IMAXREGNUM1       = NUMREGISG_ORIG 
         IMAXRSTRQ1        = MAXRSTRQG

      ELSEIF (IBANK == 2) THEN
C        SECOND  BANK
         IGORSTRQ(IBANK)   = IMAXRSTRQ1       + 1
         IENDRSTRQ(IBANK)  = IMAXRSTRQ1 
         IGOREGNUM(IBANK)  = IMAXREGNUM1      + 1
         IENDREGNUM(IBANK) = IMAXREGNUM1

      ELSE
C        THIRD, ..... BANK
         IGORSTRQ(IBANK)   = IENDRSTRQ(IBANK-1)  + 1
         IENDRSTRQ(IBANK)  = IGORSTRQ(IBANK)     - 1 
         IGOREGNUM(IBANK)  = IENDREGNUM(IBANK-1) + 1
         IENDREGNUM(IBANK) = IGOREGNUM(IBANK)    - 1
      ENDIF

C     INITIAL < IN RSTRQ
      RSTRQ(IGORSTRQ(IBANK):IGORSTRQ(IBANK)) = '<'  
      IENDRSTRQ(IBANK) = IGORSTRQ(IBANK) 

C     SET LOOP REG NONE IN BANK
      CALL REG_SET_VAR(IBANK,'[_0]',.TRUE.,0.0,IDUM,IRTFLG)

C     SET ERROR FLAG TO NONE IN BANK
      CALL REG_SET_VAR(IBANK,'[_9]',.TRUE.,0.0,IDUM,IRTFLG)

      IRTFLG = 0
 
      RETURN
      END

C++*********************************************************************
C
C REG_REINIT                    NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_REINIT()
C
C    PURPOSE:     RESIZES   REGISTER SPACE   
C
C    PARAMETERS:  IRTFLG    ERROR FLAG                          (RET.)
C
C    YES, I KNOW THAT IT SHOULD BE WRITTEN USING POINTERS BUT
C    I DOUBT ANYONE WILL EVER USE THIS!! al
C
C--*******************************************************************

      SUBROUTINE REG_REINIT(IRTFLG)

      USE REG_STUFF

      INCLUDE 'CMBLOCK.INC'
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: REGVALUEST

      CALL REG_GET_NUMS(NREG,NCHAR)
      WRITE(NOUT,90)NREG,NCHAR  
90    FORMAT(' CURRENT REGISTERS: ',I7,' NAME CHARACTERS: ',I8)
  
      CALL RDPRI1S(NREGN,NOT_USED,'NUMBER OF REGISTERS WANTED',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (NREGN > NREG) THEN
         ALLOCATE (REGVALUEST(NREG), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'UNABLE TO INCREASE REGISTERS:',NREG)
            RETURN
         ENDIF
              
         REGVALUEST(1:NUMREGIS) = REGVALUES
         DEALLOCATE(REGVALUES)

         ALLOCATE (REGVALUES(NREGN), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'UNABLE TO INCREASE REGISTERS:',NREGN)
            RETURN
         ENDIF

         REGVALUES(1:NUMREGIS) = REGVALUEST
         DEALLOCATE(REGVALUEST)
         
         NUMREGIS = NREGN

      ENDIF

      IRTFLG = 0
 
      RETURN
      END


C++*********************************************************************
C
C REG_BANK_OK                    NEW            AUG 2005 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_BANK_OK(IBANKIN,IBANKOUT,IRTFLG)
C
C    PURPOSE:     FINDS IF REGISTER BANK EXISTS, RETURNS BANK  
C
C    PARAMETERS:  IBANKIN   STACK LEVEL (0 IS CURRENT ISTOP)    (SENT)
C                 IBANKOUT  STACK LEVEL                         (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_BANK_OK(IBANKIN,IBANKOUT,IRTFLG)

      USE REG_STUFF

      INCLUDE 'CMBLOCK.INC'
      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      IF (IBANKIN == -9999) THEN
         IBANKOUT = 1      ! UNUSED CAPABILITY FOR ALLBANK
      ELSEIF (IBANKIN < 0) THEN
         IBANKOUT = ISTOP + IBANKIN
      ELSEIF (IBANKIN == 0) THEN
         IBANKOUT = ISTOP
      ELSE
         IBANKOUT = IBANKIN
      ENDIF

      IF (IBANKOUT < 0 .OR. IBANKOUT > MAXPRC) THEN
         CALL ERRT(102,'ILLEGAL REGISTER VARIABLE BANK',IBANKOUT)
         IRTFLG = 1
         RETURN

      ELSEIF (IBANKOUT > ISTOP) THEN
         CALL ERRT(102,'BANK NOT IN CURRENT USE',IBANKOUT)
         IRTFLG = 1
         RETURN
      ENDIF

      IRTFLG   = 0
      END

     
C++*********************************************************************
C
C REG_FIND                       NEW            AUG 2005 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_FIND(IBANKT,NAME,VALUE,IREG,IRTFLG)
C
C    PURPOSE:     FINDS IF VARIABLE EXISTS, RETURNS VALUE & IREG  
C
C    PARAMETERS:  IBANK     STACK LEVEL () IS CURRENT TOP       (SENT)
C                 NAME      REGISTER NAME, WITH [...]           (SENT)
C                 VALUE     VALUE                               (RET.)
C                 IREG      REGISTER NUMBER                     (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_FIND(IBANKT,NAME,VALUE,IREG,IRTFLG)

      USE REG_STUFF

      INCLUDE 'CMBLOCK.INC'

C     COMMON NEEDED FOR: ISTOP
      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      CHARACTER(LEN=*)   :: NAME
      CHARACTER(LEN=160) :: NAMET

C     SEE IF REGISTER VARIABLE EXISTS YET

C     <> ARE USED AS VARIABLE ID DELIMITERS IN RSTRQ
      LENT             = LEN(NAME)
      NAMET            = '<' // NAME(2:LENT-1) // '>' 

C     FIND AND CHECK BANK NUMBER
      CALL REG_BANK_OK(IBANKT,IBANK,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IGOQ   = IGORSTRQ(IBANK)
      IENDQ  = IENDRSTRQ(IBANK)
      IRTFLG = 0
      IREG   = 0

      IF (IENDQ > IGOQ) THEN

C        SEARCH FOR THIS VARIABLE AT THIS LEVEL
         CALL GETREGVAR(RSTRQ(IGOQ:IENDQ),NAMET(1:LENT),'<',
     &                   IGO,IEND,IRTFLG)

         IF (IRTFLG == 0) THEN
C           VARIABLE EXISTS, FIND REGISTER NUMBER
            READ(RSTRQ(IGOQ+IGO-1:IGOQ+IEND-1),*) IREG
            VALUE = REGVALUES(IREG)
         ENDIF
      ENDIF

!     if (ireg > 0) then
!        write(6,90) namet(1:lent),ibank,ireg,value
!90       format('  regfind:',t15,a,t27,i8,i6,f9.3)
!      else
!        write(6,90),namet(1:lent),ibank,ireg
!      endif

      END



C++*********************************************************************
C
C REG_NEW                         NEW            AUG 2005 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_NEW(IBANK,NAME,VALUE,IREG,IRTFLG)
C
C    PURPOSE:     CREATES A NEW REGISTER VAR & VALUE   
C
C    PARAMETERS:  IBANK     STACK LEVEL                         (SENT)
C                 NAME      REGISTER NAME, WITH [...]           (SENT)
C                 VALUE     VALUE                               (SENT)
C                 IREG      REGISTER NUMBER                     (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C     NOTES: REGISTER VARIABLE SHOULD NOT ALREADY EXIST (NOT CHECKED)
C
C--*******************************************************************

      SUBROUTINE REG_NEW(IBANKT,NAME,VALUE,IREG,IRTFLG)

      USE REG_STUFF

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)  :: NAME
      CHARACTER(LEN=80) :: CREG

C     IF REGISTER VARIABLE ALREADY EXISTS (ERROR)

C     CHECK BANK NUMBER
      CALL REG_BANK_OK(IBANKT,IBANK,IRTFLG)

C     VARIABLE DOES NOT EXIST, CREATE IT
C     WRITE(NDAT,*)' CREATING REGISTER VARIABLE: ',NAME

C     INCREMENT REGISTER NUMBER
      IREG = IENDREGNUM(IBANK) + 1
      CALL INTTOCHAR(IREG,CREG,LENR,1)

      IF (IBANK == 1 .AND. IREG > IMAXREGNUM1) THEN
C        OVER-RUN OF GLOBAL REGISTER VALUE ARRAY
         IT     = IMAXREGNUM1
         CALL ERRT(102,'TOO MANY GLOBAL REGISTERS REQUESTED, LIMIT',IT)
         IRTFLG = 1
         RETURN

      ELSEIF (IREG > NUMREGIS) THEN
C        OVER-RUN OF REGISTER VALUE ARRAY
         IT = NUMREGIS

         WRITE(nout,*) '  VARIABLE NAME: ',NAME
         write(nout,*) ' BANK:      ',IBANK
         write(nout,*) ' IGORSTRQ:  ',IGORSTRQ
         write(nout,*) ' IGOREGNUM: ',IGOREGNUM
         write(nout,*) ' IENDREGNUM: ',IENDREGNUM
         do i=2401,3900,60
            write(nout,*) I,' RSTRQ(I:I+60): ',RSTRQ(I:i+59)
         enddo
       write(nout,*)'RSTRQ: ',RSTRQ(IGORSTRQ(IBANK):IGORSTRQ(IBANK)+60)

         WRITE(NOUT,*) 'igo,iend,nchar: ',igo,iend,nchar
         CALL ERRT(102,'TOO MANY REGISTERS REQUESTED, LIMIT',IT)
         IRTFLG = 1
         RETURN
      ENDIF

C     PLACE VARIABLE STRING IN RSTRQ ARRAY
      LENT    = LEN(NAME)
      LENADD  = LENT + LENR 

      IENDQAT = IENDRSTRQ(IBANK) + LENADD 

      IF (LENT > 80) THEN
         WRITE(NOUT,*) '  VARIABLE NAME: ',NAME
         CALL ERRT(102, 'OVERLENGTH VARIABLE NAME',LENT)
         IRTFLG = 1
         RETURN
      ELSEIF (LENR > 80) THEN
         WRITE(NOUT,*) '  VARIABLE SELECTOR: ',CREG
         CALL ERRT(102, 'OVERLENGTH REGISTER SELECTOR',LENR)
         IRTFLG = 1
         RETURN
      ELSEIF (LENADD > 92) THEN
         WRITE(NOUT,*) '  VARIABLE NAME: ',NAME
         WRITE(NOUT,*) '  VARIABLE SELECTOR: ',CREG
         CALL ERRT(102, 'OVERLENGTH TOTAL REG. SELECTOR',LENADD)
         IRTFLG = 1
         RETURN
      ENDIF

      IF (IBANK == 1 .AND. IREG > IMAXRSTRQ1) THEN
C        OVER-RUN OF GLOBAL REGISTER VALUE ARRAY
         CALL ERRT(102,
     &         'RSTRQ OVERFLOW, GLOBAL VARIABLE NAMES ARRAY IS FULL',
     &         IMAXRSTRQ1)
         IRTFLG = 1
         RETURN

      ELSEIF (IENDQAT >= MAXRSTRQ) THEN
C        OVER-RUN OF RSTRQ NAMESPACE ARRAY

       write(nout,*) ' BANK:      ',IBANK
       write(nout,*) ' IGORSTRQ:  ',IGORSTRQ
       write(nout,*) ' IGOREGNUM: ',IGOREGNUM
       write(nout,*) ' IENDREGNUM: ',IENDREGNUM
       do i=2401,3900,60
          write(nout,*) I,' RSTRQ(I:I+60): ',RSTRQ(I:i+59)
       enddo
       write(nout,*) 'RSTRQ: ',RSTRQ(IGORSTRQ(IBANK):IGORSTRQ(IBANK)+60)

       WRITE(NOUT,*) 'igo,iend,nchar: ',igo,iend,nchar
       stop
   
         CALL ERRT(102,
     &         'RSTRQ OVERFLOW, VARIABLE NAMES ARRAY IS FULL',IENDQAT)
         IRTFLG = 1
         RETURN
      ENDIF

C     FIND START FOR VARIABLE DEFINITION SEQ.
      IGOQ = IENDRSTRQ(IBANK) 

      IGOQ = IENDRSTRQ(IBANK) + 1
      RSTRQ(IGOQ:IGOQ+LENADD-1) =  NAME(2:LENT-1) // '>' // 
     &                                CREG(1:LENR)   // '<' 
       
      IENDRSTRQ(IBANK)  = IENDQAT 

      REGVALUES(IREG)   = VALUE
      IENDREGNUM(IBANK) = IREG

      IRTFLG = 0

      RETURN
      END






C++*********************************************************************
C
C REG_SET_VAR.F   
C                 NATIVE []                       NOV 2005 ARDEAN LEITH
C **********************************************************************
C
C    REG_SET_VAR(IBANK,STRING,CREATE,VALUE,IREG,IRTFLG)
C
C    PURPOSE:          SET REGISTER VAR FROM STRING, RETURNS THE
C                      REGISTER NUMBER AND IT'S VALUE
C
C    PARAMETERS:       IBANK     STACK LEVEL                     (SENT)
C                      STRING    OPERATION STRING                  SENT
C                      CREATE    CREATE FLAG                       SENT
C                      VALUE     VALUE                             SENT
C                      IREG      NO. OF REGISTER VARIABLE          RET.
C                      IRTFLG    ERROR FLAG                        RET.
C
C--*******************************************************************

      SUBROUTINE REG_SET_VAR(IBANK,STRING,CREATE,VALUE,IREG,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR NALPH
      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)   :: STRING
      LOGICAL :: CREATE
      LOGICAL :: ISDIGI

      NCHAR  = LEN(STRING)
      IFIRST = 1

c     GET REG. VAR. (CHAR. STRING DELIMITED BY [])  FROM STRING
      CALL GETNEXTVAR(STRING,IFIRST,IGO,IEND,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     VARIABLE STRING FOUND, FIND REG. NUMBER FOR LIST
      CALL REG_FIND(IBANK,STRING(IGO:IEND),VALDUM,IREG,IRTFLG)
      IF (IREG <= 0) THEN
C        REGISTER DOES NOT EXIST
         IF (CREATE) THEN
C           CREATE THE REGISTER
            CALL REG_NEW(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)
         ELSE
C           DO NOT WANT TO CREATE REG. IF DOES NOT ALREADY EXIST
            CALL ERRT(8,STRING(IGO:IEND),NE)

#ifdef NEVER
            IF ((STRING(IGO:IGO) == '_') .AND.
     &          ISDIGI(STRING(IGO+1:IGO+1))) THEN
                VALPREV = 0.0
                IF (IBANK > 1) THEN
C                  TRY TO COPY REGISTER FROM PREVIOUS BANK
                   CALL REG_FIND(IBANK-1,STRING(IGO:IEND),VALPREV,
     &                        IREG,IRTFLG)
                ENDIF
                CALL REG_NEW(IBANK,STRING(IGO:IEND),VALPREV,IREG,IRTFLG)

                WRITE(NOUT,*)' *** UNDEFINED REGISTER VARIABLE: X',
     &                   STRING(IGO+1:IEND-1)
                WRITE(NOUT,*)' *** PLEASE FIX THIS AS IT MAY NOT BE',
     &                       ' ACCEPTED IN FUTURE SPIDER RELEASES'
            ENDIF
#endif
         ENDIF
      ELSE
         REGVALUES(IREG) = VALUE 
      ENDIF

      RETURN
      END

C++*********************************************************************
C
C REG_SET_BYNUM.F   
C                 NATIVE []                       NOV 2005 ARDEAN LEITH
C **********************************************************************
C
C    REG_SET_BYNUM(IREG,VALUE,IRTFLG)
C
C    PURPOSE:          SET REGISTER VAR IREG = VALUE
C
C    PARAMETERS:       IREG      NO. OF REGISTER VARIABLE          RET.
C                      VALUE     VALUE                             SENT
C                      IRTFLG    ERROR FLAG                        RET.
C
C--*******************************************************************

      SUBROUTINE REG_SET_BYNUM(IREG,VALUE,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR NALPH
      INCLUDE 'CMBLOCK.INC'
      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      IF (IREG <= 0) THEN
C        REGISTER CAN NOT EXIST
         CALL ERRT(102,'BAD REGISTER NUMBER',IREG)
         IRTFLG = 1

      ELSEIF (IREG > IENDREGNUM(ISTOP)) THEN
C        REGISTER DOES NOT EXIST (BEYOND END)
         CALL ERRT(102,'BAD REGISTER NUMBER',IREG)
         IRTFLG = 1

      ELSE
         REGVALUES(IREG) = VALUE 
         IRTFLG = 0
      ENDIF

      RETURN
      END


C++*********************************************************************
C
C REG_GET_BYNUM.F   
C                 NATIVE []                       NOV 2005 ARDEAN LEITH
C **********************************************************************
C
C    REG_GET_BYNUM(IREG,VALUE,IRTFLG)
C
C    PURPOSE:          SET REGISTER VAR IREG = VALUE
C
C    PARAMETERS:       IREG      NO. OF REGISTER VARIABLE          SENT
C                      VALUE     VALUE                             RET.
C                      IRTFLG    ERROR FLAG                        RET.
C
C--*******************************************************************

      SUBROUTINE REG_GET_BYNUM(IREG,VALUE,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR NALPH
      INCLUDE 'CMBLOCK.INC'
      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                       IPNUMSTACK,IPARNUM

      IF (IREG <= 0) THEN
C        REGISTER CAN NOT EXIST
         CALL ERRT(102,'BAD REGISTER NUMBER',IREG)
         IRTFLG = 1

      ELSEIF (IREG > IENDREGNUM(ISTOP)) THEN
C        REGISTER DOES NOT EXIST (BEYOND END)
         CALL ERRT(102,'REGISTER NUMBER NOT DEFINED',IREG)
         IRTFLG = 1

      ELSE
         VALUE  = REGVALUES(IREG)
         IRTFLG = 0
      ENDIF

      RETURN
      END

C++*********************************************************************
C
C REG_GET_VAR.F   
C                 NATIVE []                       NOV 2005 ARDEAN LEITH
C **********************************************************************
C
C    REG_GET_VAR(IBANK,STRING,CREATE,VALUE,IREG,IEND,IRTFLG)
C
C    PURPOSE:          EXTRACTS REGISTER VAR FROM STRING, RETURNS THE
C                      REGISTER NUMBER AND IT'S VALUE
C
C    PARAMETERS:       IBANK     STACK LEVEL                       SENT 
C                      STRING    OPERATION STRING                  SENT
C                      CREATE    CREATE FLAG                       SENT
C                      VALUE     VALUE                             RET.
C                      IREG      NUMBER FOR REGISTER VARIABLE      RET.
C                      IEND      LAST VARIABLE POSITION IN STRING  RET.
C                      IRTFLG    ERROR FLAG                        RET.
C
C--*******************************************************************

      SUBROUTINE REG_GET_VAR(IBANK,STRING,CREATE,VALUE,IREG,IEND,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR: NOUT
      INCLUDE 'CMBLOCK.INC'

C     COMMON NEEDED FOR: ISTOP
      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      CHARACTER(LEN=*) :: STRING
      LOGICAL          :: CREATE,ISDIGI

      NCHAR  = LEN(STRING)
      IFIRST = 1

c     GET REG. VAR. (CHAR. STRING DELIMITED BY [])  FROM STRING
      CALL GETNEXTVAR(STRING,IFIRST,IGO,IEND,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     VARIABLE STRING FOUND, FIND REG. NUMBER FOR LIST
      CALL REG_FIND(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)

!      if (ireg > 0) then
!         write(6,90) string(igo:iend),istop,ireg,value
!90       format('  get_var:',t15,a,t27,i8,i6,f9.3)
!      else
!         write(6,90) string(igo:iend),ibank,ireg
!      endif

      IF (IREG > 0) RETURN     ! REGISTER VARIABLE FOUND OK

C     SEE IF FOUND IN GLOBAL BANK  
      CALL REG_FIND(1,STRING(IGO:IEND),VALUE,IREG,IRTFLG)

!      ibankz = 1
!      if (ireg > 0) then
!         ibankz = 1
!         write(6,91) string(igo:iend),ibankz,ireg,value
!91       format('  get_var1:',t15,a,t27,i8,i6,f9.3)
!      else
!         write(6,91) string(igo:iend),ibankz,ireg
!      endif


      IF (IREG <= 0) THEN
C        NOT FOUND IN IBANK OR IN GLOBAL BANK

         IF (CREATE) THEN
C           CREATE THE VARIABLE
            CALL REG_NEW(IBANK,STRING(IGO:IEND),0.0,IREG,IRTFLG)

         ELSE
C           DO NOT CREATE, WANT EXISTING VARIABLE ONLY
            WRITE(NOUT,*)'  '

            IF ((STRING(IGO+1:IGO+1) == '_') .AND.
     &          ISDIGI(STRING(IGO+2:IGO+2))) THEN
                VALUE = 0.0
C               TRY TO COPY REGISTER FROM PREVIOUS BANK
                CALL REG_FIND(-1,STRING(IGO:IEND),VALUE,
     &                        IREG,IRTFLG)

                CALL REG_NEW(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)
                WRITE(NOUT,*)' ***  UNDEFINED REGISTER VARIABLE: X',
     &                       STRING(IGO+2:IEND-1)
                WRITE(NOUT,*)'      OK NOW, BUT USE OF UNDEFINED ',
     &                       'REGISTER VARIABLES '

                WRITE(NOUT,*)'      MAY NOT BE ACCEPTED IN FUTURE ',
     &                       'SPIDER RELEASES'
                IRTFLG = 0
            ELSE 
                CALL ERRT(8,STRING(IGO:IEND),NE)
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END


C++*********************************************************************
C
C REG_FIND_IREG.F   
C                 FROM REG_GET_VAR                JAN 2010 ARDEAN LEITH
C **********************************************************************
C
C    REG_FIND_IREG(TYPE,STRING,ISGLOBAL,IREG,IRTFLG)
C
C    PURPOSE:  EXTRACTS REGISTER VAR, NAME FROM STRING, RETURNS THE
C              GLOBAL REGISTER NUMBER, IF NO GLOBAL REG. EXISTS AND
C             , CREATES SPECIFIED GLOBAL REG. OTHERWISE 
C              RETURNS LOCAL REG. IF EXISTS OR CREATES IT. 
C
C    PARAMETERS:       TYPET      GLOBAL STRING                    SENT 
C                      ISGLOBAL   GLOBAL FLAG                      RET. 
C                      STRING     OPERATION STRING                 SENT
C                      IREG       NUMBER FOR REGISTER VARIABLE     RET.
C                      IRTFLG     ERROR FLAG                       RET.
C
C--*******************************************************************

      SUBROUTINE REG_FIND_IREG(TYPET,STRING,ISGLOBAL,IREG,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR: ISTOP
      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

C     COMMON NEEDED FOR: FCHAR
      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=3) :: TYPE,TYPET
      INTEGER          :: IREG,IBANK,IRTFLG
      LOGICAL          :: ISGLOBAL

      NCHAR  = LEN(STRING)
      IFIRST = 1

c     GET REG. VAR. (CHAR. STRING DELIMITED BY [])  FROM STRING
      CALL GETNEXTVAR(STRING,IFIRST,IGO,IEND,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      TYPE = TYPET
      IF (FCHAR(1:2) == 'UD') THEN
C        ALWAYS WANT A LOCAL VARIABLE CREATION
         TYPE = 'LOC'
      ENDIF

      ISGLOBAL = .FALSE.

      IF (TYPE == 'GLO') THEN
C        USE GLOBAL VARIABLE, CREATE IF IT DOES NOT EXIST

         IBANK = 1

C        SEE IF FOUND IN GLOBAL BANK (IREG WILL BE > 0)  
         CALL REG_FIND(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)

         IF (IREG <= 0) THEN
            CALL REG_NEW(IBANK,STRING(IGO:IEND),0.0,IREG,IRTFLG)
         ENDIF
         ISGLOBAL = .TRUE.

      ELSEIF (TYPE == 'LOC') THEN
C        SEE IF NAMED REGISTER IS IN LOCAL BANK
         IBANK = 0
         CALL REG_FIND(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)

         IF (IREG <= 0) THEN
C           NOT FOUND IN LOCAL BANK, CREATE LOCAL VARIABLE
            IBANK = 0
            CALL REG_NEW(IBANK,STRING(IGO:IEND),0.0,IREG,IRTFLG)
        ENDIF

      ELSEIF (TYPE == 'RED') THEN
C        SEE IF NAMED REGISTER IS IN LOCAL BANK
         IBANK = 0
         CALL REG_FIND(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)

         !write(6,*) ' type:',type, ' ibank:',ibank,'string:',string,
         !' value:',value,' ireg:',ireg,' irtflg:',irtflg

         IF (IREG <= 0) THEN
C           NOT FOUND IN LOCAL BANK, USE GLOBAL VARIABLE IF EXISTS

C           SEE IF FOUND IN GLOBAL BANK (IREG WILL BE > 0) 
            IBANK = 1
            CALL REG_FIND(IBANK,STRING(IGO:IEND),VALUE,IREG,IRTFLG)
         
            IF (IREG <= 0) THEN
C              NOT FOUND IN GLOBAL BANK, CREATE LOCAL VARIABLE INSTEAD
               IBANK = 0
               CALL REG_NEW(IBANK,STRING(IGO:IEND),0.0,IREG,IRTFLG)            ISGLOBAL = .TRUE.
               ISGLOBAL = .FALSE.
            ELSE
               ISGLOBAL = .TRUE.
            ENDIF
         ENDIF
      ENDIF

!      if (isglobal) then
!         write(6,90) string(igo:iend),ibank,ireg,type
!90       format('  find_ireg:',t15,a,t27,i8,i6,'          :',a,':')
!      else
!         write(6,90) string(igo:iend),istop,ireg,type
!      endif

      RETURN
      END

C++*********************************************************************
C
C REG_GET_SEL.F   ADAPTED FROM READP.FOR FOR CHAR. AUG 1989 ARDEAN LEITH
C                 NATIVE []                        NOV 2005 ARDEAN LEITH
C                 REG_FIND_IREG                    JAN 2009 ARDEAN LEITH
C **********************************************************************
C
C    REG_GET_SEL(IBANK,STRING,CREATE,WANTGLO,NREG,IRTFLG)
C
C    PURPOSE:          PARSES REGISTER LIST FROM INPUT LINE.  PLACES
C                      REGISTER NUMBERS (UP TO MAXNSEL NUMBERS) IN 
C                      NSELREG.
C
C    PARAMETERS:       IBANK     STACK LEVEL                       SENT 
C                      STRING    OPERATION STRING                  SENT
C                      CREATE    FLAG TO CREATE VAR.               SENT
C                      WANTGLO   WANT EXISTIN GLOBAL REG           SENT
C                      NREG      NO. OF REGISTER VARIABLES         RET.
C                      IRTFLG    ERROR FLAG                        RET.
C
C--*******************************************************************

      SUBROUTINE REG_GET_SEL(IBANK,STRING,CREATE,WANTGLO,NREG,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR: NOUT
      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)   :: STRING
      CHARACTER(LEN=10)  :: REGNAME
      LOGICAL            :: CREATE,WANTGLO,ISGLOBAL

      NCHAR = LEN(STRING)
      NREG   = 0
      IFIRST = 1
      IRTFLG = 1

      DO
c        GET NEXT VAR. (CHAR. STRING DELIMITED BY [] ) FROM STRING
         CALL GETNEXTVAR(STRING,IFIRST,IGO,IEND,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         IF (IGO <= 0 .OR. IEND <= IGO) THEN
C           ALL TOKENS FROM STRING HAVE BEEN EVALUATED
            IRTFLG = 0
            EXIT
         ENDIF

         IF (IGO  > 0 .AND. IEND <= IGO) THEN
C           ERROR
            write(nout,*) ' IFIRST:     ',IFIRST,IGO,IEND
            write(nout,*) ' BANK:       ',IBANK
            write(nout,*) ' NCHAR:      ',NCHAR
            write(nout,*) ' STRING:     ',STRING
            write(nout,*) ' IGORSTRQ:   ',IGORSTRQ
            write(nout,*) ' IENDRSTRQ:  ',IENDRSTRQ
            write(nout,*) ' IGOREGNUM:  ',IGOREGNUM
            write(nout,*) ' IENDREGNUM: ',IENDREGNUM
            write(nout,*) ' NSEL_USED:  ',NSEL_USED
            write(nout,*) ' NSELREG:    ',NSELREG
            do i=2401,IENDRSTRQ(IBANK),60
               write(nout,*) I,' RSTRQ(I:I+60): ',RSTRQ(I:i+59)
            enddo
            CALL GETNEXTVARbug(STRING,IFIRST,IGO,IEND,IRTFLG)
            CALL ERRT(102,'TOO MAY REGISTERS ON THIS LINE',NREG)
            NSEL_USED = NREG -1
            IRTFLG = 1
         stop
            RETURN
         ENDIF

C        VARIABLE FOUND, FIND REG. NUMBER FOR LIST

         NREG = NREG + 1
         IF (NREG > MAXNSEL) THEN
C           ERROR
            write(nout,*) ' IFIRST:     ',IFIRST,IGO,IEND
            write(nout,*) ' BANK:       ',IBANK
            write(nout,*) ' NCHAR:      ',NCHAR
            write(nout,*) ' STRING:     ',STRING
            write(nout,*) ' IGORSTRQ:   ',IGORSTRQ
            write(nout,*) ' IENDRSTRQ:  ',IENDRSTRQ
            write(nout,*) ' IGOREGNUM:  ',IGOREGNUM
            write(nout,*) ' IENDREGNUM: ',IENDREGNUM
            write(nout,*) ' NSEL_USED:  ',NSEL_USED
            write(nout,*) ' NSELREG:    ',NSELREG
            do i=2401,IENDRSTRQ(IBANK),60
               write(nout,*) I,' RSTRQ(I:I+60): ',RSTRQ(I:i+59)
            enddo
            CALL GETNEXTVARbug(STRING,IFIRST,IGO,IEND,IRTFLG)
            CALL ERRT(102,'TOO MANY REGISTERS ON THIS LINE',NREG)
            NSEL_USED = NREG -1
            IRTFLG = 1
            STOP
            RETURN
         ENDIF

         IF (STRING(IGO:IGO) .NE. '[') THEN
C           JUST A NUMBER, NOT A REGISTER VARIABLE
            WRITE(NOUT,*) 'GETNEXVAR RETURNED:',STRING(IGO:IEND)
            WRITE(NOUT,*) 'SHOULD HAVE BEEN VARIABLE',IGO,IEND
            CALL ERRT(100,'REG_GET_SEL',NE)
            STOP
         ENDIF

         IF (IBANK == 1) THEN

C           ACCESS EXISTING GLOBAL REGISTER OR CREATE NEW GLOBAL ONE 
            CALL REG_FIND_IREG('GLO',STRING(IGO:IEND),
     &                         ISGLOBAL,IREG,IRTFLG)
         ELSEIF (WANTGLO) THEN

C           ACCESS EXISTING LOCAL/GLOBAL REGISTER OR CREATE NEW LOCAL ONE 
            CALL REG_FIND_IREG('RED',STRING(IGO:IEND),
     &                         ISGLOBAL,IREG,IRTFLG)
         ELSE

C           ACCESS EXISTING LOCAL REGISTER OR CREATE NEW LOCAL ONE 
            CALL REG_FIND_IREG('LOC',STRING(IGO:IEND),
     &                         ISGLOBAL,IREG,IRTFLG)
         ENDIF

         NSELREG(NREG) = IREG
         IFIRST        = IEND + 1
         IF (IFIRST > NCHAR) EXIT
          
      ENDDO

      NSEL_USED = NREG
      IRTFLG    = 0

      END

C++*********************************************************************
C
C REG_GET_SELS.F   NATIVE []                      NOV 2005 ARDEAN LEITH
C **********************************************************************
C
C    REG_GET_SELS(IREGSELS,NLEN,NREG,IRTFLG)
C
C    PURPOSE: RETRIEVES REGISTER NUMBERS FROM NSELREG.
C
C    PARAMETERS: IREGSELS     SELECTED REGISTER LIST.             RET.
C                NLEN         LENGTH OF SELECTED REGISTER LIST.   SENT
C                NREG         NO. OF REGISTER VARIABLES           RET.
C                IRTFLG       ERROR FLAG (UNUSED)                 RET.
C
C--*******************************************************************

      SUBROUTINE REG_GET_SELS(IREGSELS,NLEN,NREG,IRTFLG)

      USE REG_STUFF

C     COMMON NEEDED FOR NALPH
      INCLUDE 'CMBLOCK.INC'

      INTEGER, DIMENSION(*) :: IREGSELS

      IRTFLG = 0
      NREG   = NSEL_USED

      IF (NREG > NLEN) THEN
         CALL ERRT(102,'TOO MANY REGISTERS SPECIFIED',NREG)
         NREG   = NLEN
         IRTFLG = 1
      ENDIF

      DO I = 1,NREG
         IREGSELS(I) = NSELREG(I)
      ENDDO 

      END


C++*********************************************************************
C
C REG_LIST_COPY.F   NATIVE []                       NOV 2005 ARDEAN LEITH
C **********************************************************************
C
C    REG_LIST_COPY(NL,LISTIN,LISTOUT)
C
C    PURPOSE:  COPIES LISTIN REGISTER VALUES TO LISTOUT REGISTERS
C
C    PARAMETERS: NL        NUMBER OF VALUES IN LIST             SENT
C                LISTIN    LIST OF INPUT REG. NUMBERS           SENT
C                LISTOUT   LIST OF OUTPUT REGISTER NUMBERS      RET.
C
C--*******************************************************************

      SUBROUTINE REG_LIST_COPY(NL,LISTIN,LISTOUT)

      USE REG_STUFF

      INTEGER, DIMENSION(NL) :: LISTIN,LISTOUT

      DO I = 1,NL
         REGVALUES(LISTOUT(I)) =  REGVALUES(LISTIN(I)) 
      ENDDO

      RETURN
      END




C++*********************************************************************
C
C REG_SET_NSEL                     NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_SET_NSEL(IGO,VAL0,VAL1,VAL2,VAL3,VAL4,IRTFLG)
C
C    PURPOSE:     SETS A REGISTER SPECIFIED IN NSELREG(NVAL) TO VALUE   
C
C    PARAMETERS:  IGO       STARTING REGISTER NUMBER           (SENT)
C                 VAL...    VALUES                             (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_SET_NSEL(IGO,NVAL,VAL0,VAL1,VAL2,VAL3,VAL4,IRTFLG)

      USE REG_STUFF

C     NSELREG IS STILL CARRYING ADJUSTED REG NUMBER (+1)

      NVALS = MIN(NVAL+IGO-1, NSEL_USED)
      
      IF (NVALS >= IGO+0) THEN
         REGVALUES(NSELREG(IGO+0)) = VAL0
         IF (NVALS >= IGO+1) THEN
            REGVALUES(NSELREG(IGO+1)) = VAL1
            IF (NVALS >= IGO+2) THEN
               REGVALUES(NSELREG(IGO+2)) = VAL2
               IF (NVALS >= IGO+3) THEN
                  REGVALUES(NSELREG(IGO+3)) = VAL3
                  IF (NVALS >= IGO+4) THEN
                     REGVALUES(NSELREG(IGO+4)) = VAL4
                  ENDIF 
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END


C++*********************************************************************
C
C REG_SET_NSELA                    NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_SET_NSELA(NREG,VALUES,FILLALL,IRTFLG)
C
C    PURPOSE:     SETS REGISTERS SPECIFIED IN NSEL TO VALUES   
C
C    PARAMETERS:  NREGT     NUMBER OF VARIABLES TO SET          (SENT)
C                 VALUES    VALUES ARRAY                        (SENT)
C                 FILLALL   FLAG TO ZERO REST OF SEL. VAR.      (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_SET_NSELA(NREGT,VALUES,FILLALL,IRTFLG)

      USE REG_STUFF

      REAL,DIMENSION(*) :: VALUES
      LOGICAL           :: FILLALL

C     ONLY SET A MAX OF: NSEL_USED REGISTERS 
      NREG = MIN(NREGT,NSEL_USED)

      IF (NREG > 0) THEN
C        NSELREG CONTAINS: NSEL_USED REGISTER NUMBERS
         DO IREG=1,NREG
            REGVALUES(NSELREG(IREG)) = VALUES(IREG)
         ENDDO
 
         IF (NREG < NSEL_USED) THEN
            DO IREG=NREG,NSEL_USED
               REGVALUES(NSELREG(IREG)) = 0.0
            ENDDO
         ENDIF
      ENDIF

      
      RETURN
      END

C++*********************************************************************
C
C REG_GET_NSELA                    NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_GET_NSELA(NREG,VALUES,FILLALL,IRTFLG)
C
C    PURPOSE:     GETS VALUES FROM REGISTER(S) LISTED IN NSEL    
C
C    PARAMETERS:  NREGT     NUMBER OF VARIABLES TO GET          (SENT)
C                 VALUES    VALUES ARRAY                        (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_GET_NSELA(NREGT,VALUES,IRTFLG)

      USE REG_STUFF

      REAL,DIMENSION(*) :: VALUES

C     ONLY GET A MAX OF: NSEL_USED REGISTERS 
      NREG = MIN(NREGT,NSEL_USED)

      DO IREG = 1,NREG
         CALL REG_GET_BYNUM(NSELREG(IREG),VALUES(IREG),IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDDO
 
     
      RETURN
      END




C++*********************************************************************
C
C REG_GET                         NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_GET(IBANK,IXREG,CXREG,VALUE,IREGRET,IRTFLG)
C
C    PURPOSE:     GETS A CURRENT REGISTER VALUE FROM X OR CXREG INPUT   
C
C    PARAMETERS:  IBANK     IBANK NUMBER                        (SENT)
C                 IXREG     REGISTER NUMBER X#                  (SENT)
C                 CXREG     INDEX REGISTER                      (SENT.)
C                 VALUE     VALUE FOR REGISTER                  (RET.)
C                 IREGRET   REGISTER #                          (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_GET(IBANK,IXREG,CXREG,VALUE,IREGRET,IRTFLG)

      USE REG_STUFF

      CHARACTER(LEN=1)  :: CXREG
      CHARACTER(LEN=80) :: REGNAME

      REGNAME(1:2) = '[_'
      IF (IXREG >= 0) THEN
         CALL INTTOCHAR(IXREG,REGNAME(3:),NLET,1)
         NLET = NLET + 3
      ELSE
         REGNAME(3:3) = CXREG
         NLET         = 4
      ENDIF
      REGNAME(NLET:NLET) = ']'

      CALL REG_GET_VAR(IBANK,REGNAME,.FALSE.,VALUE,IREGRET,IEND,IRTFLG)
 
      RETURN
      END

C++*********************************************************************
C
C REG_SET                         NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_SET(IXREG,VALUE,CXREG,IRTFLG)
C
C    PURPOSE:     SETS A CURRENT REGISTER VALUE   
C
C    PARAMETERS:  IXREG      REGISTER NUMBER X#                  (SENT)
C                 VALUE     VALUE FOR REGISTER                  (SENT)
C                 CXREG      INDEX REGISTER                      (SENT.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_SET(IXREG,VALUE,CXREG,IRTFLG)

      USE REG_STUFF

      CHARACTER(LEN=1)  :: CXREG
      CHARACTER(LEN=80) :: REGNAME

      REGNAME(1:2) = '[_'
      IF (IXREG >= 0) THEN
         CALL INTTOCHAR(IXREG,REGNAME(3:),NLET,1)
         NLET = NLET + 3
      ELSE
         REGNAME(3:3) = CXREG
         NLET         = 4
      ENDIF
      REGNAME(NLET:NLET) = ']'

      CALL REG_SET_VAR(0,REGNAME(:NLET),.TRUE.,VALUE,IREGRET,IRTFLG)

      RETURN
      END


C++*********************************************************************
C
C REG_GET_USED               NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_GET_USED(NREG)
C
C    PURPOSE:     GETS NUMBER OF REGISTER USED IN NSELREG   
C
C    PARAMETERS:  NREG      REGISTER NUMBERS IN USE             (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_GET_USED(NREG)

      USE REG_STUFF

      NREG = NSEL_USED

      RETURN
      END

      SUBROUTINE REG_SET_USED(NREG)

      USE REG_STUFF

      NSEL_USED = NREG

      RETURN
      END


C     ------------------------- REG_DOC_PARSE -----------------------------
C
C    REG_DOC_PARSE(CCHAR,COMOUT,IKEY,NLIST,IRTFLG)
C
C    PURPOSE:    SUBROUTINE TO PARSE UD & SD  TYPE LINE WHERE
C                VALUES AFTER THE FIRST ARE ALL REFERENCES TO REGISTERS
C                AND ARE RETURNED IN ILIST AS RAW REGISTER NUMBERS 
C
C    PARAMETERS:   CCHAR   INPUT LINE                             (SENT)
C                  COMOUT  COMMENT INDICATOR                      (RET.)
C                  IKEY    NUMBER OF FIRST VALUE IN CCHAR         (RET.)
C                  ILIST   ARRAY REGISTER LIST                    (RET.)
C                  NMAX    MAX LENGTH OF ARRAY REGISTER LIST      (SENT)
C                  NLIST   NUMBER OF ELEMENTS IN ARRAY            (RET.)
C                  IRTFLG  ERROR FLAG                             (RET.)
C
C--*********************************************************************

      SUBROUTINE REG_DOC_PARSE(CCHAR,COMOUT,IKEY,NLIST,IRTFLG)

      USE REG_STUFF

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*) :: CCHAR
      CHARACTER(LEN=1) :: CTEMP
      LOGICAL          :: COMOUT,ISCHAR       

C       PARSE REGISTER LINE, CHECK FOR ',' OR '[' ---------------------
       
        IRTFLG = 1
        COMOUT = .FALSE.

C       CHECK IF FIRST ENTRY IS A REGISTER, INTEGER, OR COMMENT /
        ILEN = LNBLNKN(CCHAR)

C       FIND FIRST NON-BLANK, NON-COMMA CHAR IN CCHAR
        K = VERIFY(CCHAR(1:ILEN),', ')
        IF (K <= 0) THEN
           WRITE(NOUT,90) CCHAR(1:ILEN)
          IF (NDAT .NE. 6) WRITE(6,90) CCHAR(1:ILEN)
90         FORMAT('  *** UNABLE TO PARSE REGISTER VARIABLE LINE: ',A)
           CALL ERRT(100,'REG_DOC_PARSE',NE)
           RETURN
        ENDIF

        CTEMP  = CCHAR(K:K)

        IF (CTEMP == '/') THEN
C          JUST WANT TO PUT A COMMENT IN THE DOC FILE.
           COMOUT = .TRUE.
           IRTFLG = 0
           RETURN

        ELSEIF (CTEMP == '[') THEN
C         FIRST ENTRY IS A REGISTER VARIABLE. PUT VAR. CONTENTS IN IKEY
          
C         FIND THE REGISTER NUMBER
          CALL REG_GET_NSEL(1,FKEY,FDUM,FDUM,FDUM,FDUM,IRTFLG)
          IF (IRTFLG .NE. 0) THEN
             WRITE(NOUT,90) CCHAR(1:ILEN)
             IF (NDAT .NE. 6) WRITE(NOUT,90) CCHAR(1:ILEN)
             CALL ERRT(100,'REG_DOC_PARSE',NE)
             RETURN
          ENDIF
          IKEY = FKEY

C         REGISTER VARIABLE NUMBERS WERE RETRIEVED IN RDPR.F
C         BUT FIRST REGISTER IS INTERPRETED AS THE KEY!!
C         THROW AWAY FIRST REGISTER VARIABLE

          IF (NSEL_USED > 1) THEN
             DO I = 2,NSEL_USED
                NSELREG(I-1) = NSELREG(I)
             ENDDO
             NSEL_USED = NSEL_USED - 1
          ENDIF

        ELSEIF (ISCHAR(CTEMP)) THEN
C          FIRST ENTRY IS A OLD DO LOOP INDEX, PUT ITS VALUE INTO IKEY.

           CALL REG_GET(0,-1,CTEMP,FKEY,IREGRET,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           IKEY = FKEY

        ELSE
C          FIRST ENTRY IS A ,NUMBER. SO PUT THE NUMBER IN IKEY.
           NCHAR = INDEX(CCHAR(K:),',') - 1
           READ(CCHAR(K:K+NCHAR-1),*,IOSTAT=IER) IKEY

           IF (IER .NE. 0) THEN
              WRITE(NOUT,90) CCHAR(1:ILEN)
              IF (NDAT .NE. 6) WRITE(6,90) CCHAR(1:ILEN)
              CALL ERRT(100,'REG_DOC_PARSE',NE)
              RETURN
           ENDIF
        ENDIF

        IF (IKEY == 0) THEN
           CALL ERRT(101,'*** INVALID KEY NUMBER: 0',NE)
           RETURN
        ENDIF

        NLIST = NSEL_USED
        IRTFLG = 0

        RETURN
        END



C++*********************************************************************
C
C REG_GET_NSEL                     NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_GET_NSEL(IGO,VAL0,VAL1,VAL2,VAL3,VAL4,IRTFLG)
C
C    PURPOSE:     GETS VALUES OF CONTENTS OF REGISTER(S) SPECIFIED 
C                 IN NSELREG 
C
C    PARAMETERS:  IGO       STARTING REGISTER NUMBER           (SENT)
C                 VAL...    VALUES                              (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_GET_NSEL(IGO,VAL0,VAL1,VAL2,VAL3,VAL4,IRTFLG)

      USE REG_STUFF

      IF (NSEL_USED >= IGO+0) THEN
C        NSEL IS STILL CARRYING ADJUSTED REG NUMBER (+1)
         CALL REG_GET_BYNUM(NSELREG(IGO+0),VAL0,IRTFLG)
         IF (NSEL_USED >= IGO+1) THEN
            CALL REG_GET_BYNUM(NSELREG(IGO+1),VAL1,IRTFLG)
            IF (NSEL_USED >= IGO+2) THEN
               CALL REG_GET_BYNUM(NSELREG(IGO+2),VAL2,IRTFLG)
               IF (NSEL_USED >= IGO+3) THEN
                  CALL REG_GET_BYNUM(NSELREG(IGO+3),VAL3,IRTFLG)
                  IF (NSEL_USED >= IGO+4) THEN
                     CALL REG_GET_BYNUM(NSELREG(IGO+4),VAL4,IRTFLG)
                  ENDIF 
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END


C++*********************************************************************
C
C REG_GET_NAME                 NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_GET_NAME(IREG,NAME,NLET,IRTFLG)
C
C    PURPOSE:     REVERSE LOOKUP OF REGISTER(S) SPECIFIED IN IREG 
C
C    PARAMETERS:  IREG      REGISTER POSITION IN REGVALUES     (SENT.)
C                 NAME      REGISTER NAME                      (RET.)
C                 NLET      CHARS IN REGISTER NAME             (RET.)
C                 IRTFLG    ERROR FLAG                         (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_GET_NAME(IREG,NAME,NLET,IRTFLG)

      USE REG_STUFF

      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      CHARACTER(LEN=*)  :: NAME
      CHARACTER(LEN=80) :: SEARCH
      LOGICAL           :: ISDIGI

      CALL INTTOCHAR(IREG,SEARCH(2:),NLET,1)
      SEARCH(1:1)           = '>' 
      SEARCH(NLET+2:NLET+2) = '<' 
      ILOC = INDEX(RSTRQ(IGORSTRQ(ISTOP):IENDRSTRQ(ISTOP)), 
     &                   SEARCH(1:NLET+2))

      IEND = IGORSTRQ(ISTOP) + ILOC - 2
      IGO  = INDEX(RSTRQ(1:IEND),'<',.TRUE.) + 1

      IF (RSTRQ(IGO:IGO) == '_') THEN
         IF (ISDIGI(RSTRQ(IGO+1:IGO+1))) THEN
C           OLD FASHIONED X REGISTER NAME
            NAME = 'X' // RSTRQ(IGO+1:IEND) // CHAR(0)
            NLET = IEND - (IGO + 1) + 2
         ELSE
C           INDEX LETTER
            NAME = RSTRQ(IGO+1:IGO+1) // CHAR(0)
            NLET = 1
         ENDIF
      ELSE
C        MODERN REGISTER VARIABLE NAME
         NAME = '[' // RSTRQ(IGO:IEND) // ']' // CHAR(0)
         NLET = IEND - (IGO + 1) + 4
      ENDIF

      IRTFLG = 0
      RETURN
      END

C++*********************************************************************
C
C REG_GET_NUMS                    NEW            AUG 2000 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_GET_NUMS(IREGS)
C
C    PURPOSE:     GETS TOTAL CURRENT NUMBER OF REGISTERS 
C
C    PARAMETERS:  IREGS     NUMBER OF REGISTER                  (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_GET_NUMS(IREGS,NCHAR)

      USE REG_STUFF

      IREGS = NUMREGIS
      NCHAR = MAXRSTRQ
      END


C++*********************************************************************
C
C REG_OPENPIPE                   NEW            JULY 2001 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_OPENPIPE(CXNUM,IRTFLG)
C
C    PURPOSE:     OPENS PIPE FOR REGISTERS   
C
C    PARAMETERS:  PIPENAME  PIPE NAME                           (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE REG_OPENPIPE(PIPENAME,IRTFLG)

      USE REG_STUFF

      CHARACTER (LEN=80) ::    PIPENAME
      CHARACTER (LEN=80+24) :: MSG

      OPEN(UNIT=LUNREGPIPE, FILE=PIPENAME,
     &    FORM='UNFORMATTED',
     &    ACCESS='SEQUENTIAL',
     &    STATUS='OLD',
     &    ACTION='WRITE',
     &    IOSTAT=IRTFLG)

      IF (IRTFLG .NE. 0) THEN 
         MSG = 'FAILED TO OPEN PIPE: ' // PIPENAME
         CALL ERRT(101,MSG,IRTFLG)
         RETURN
      ENDIF

      REGPIPE = .TRUE.
      IRTFLG  = 0

      RETURN
      END

C++*********************************************************************
C
C REG_PIPE                         NEW           JULY 2001 ARDEAN LEITH
C
C **********************************************************************
C
C    REG_PIPE(NAME,IRTFLG)
C
C    PURPOSE:    SENDS REGISTER VALUE DOWN LUNREGPIPE   
C
C    PARAMETERS:  NAME  REGISTER NAME                           (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*********************************************************************

      SUBROUTINE REG_PIPE(NAME,IRTFLG)

      USE REG_STUFF

      CHARACTER(LEN=*) :: NAME
      CHARACTER(LEN=8) :: CREG

      REAL           RVAL
      INTEGER * 1    I1VAL(4), I1TMP1,I1TMP2,I1TMP3,I1TMP4
      EQUIVALENCE    (RVAL,I1VAL(1))

      INOT = VERIFY(NAME,'0123456789')
      IF (INOT == 0) THEN
C        GOT AN OLD FASHIONED NUMBER
         CREG = '[_' // NAME // ']' // CHAR(0)
         CALL REG_GET_VAR(0,CREG,.FALSE.,
     &                    VALUE,IREG,IENDVAR,IERR)
      ELSE
         CALL REG_GET_VAR(0,NAME,.FALSE.,
     &                    VALUE,IREG,IENDVAR,IERR)
      ENDIF
      IF (IRTFLG .NE. 0) RETURN

      RVAL = VALUE

#ifdef __linux__
      I1TMP1   = I1VAL(1)     !THIS COULD BE SIMPLIFIED
      I1TMP2   = I1VAL(2)
      I1TMP3   = I1VAL(3)
      I1TMP4   = I1VAL(4)

      I1VAL(1) = I1TMP4
      I1VAL(2) = I1TMP3
      I1VAL(3) = I1TMP2
      I1VAL(4) = I1TMP1

C     write(0,*) ' VALUE: ',VALUE,'  --> ', RVAL
      VALUE    = RVAL
#endif

      IF (REGPIPE) THEN
C        WRITE REG NUMBER & VALUE TO NAMED PIPE
C        LINE_FEED IS NECESSARY (EVEN ON LINUX), DO NOT ASK ME WHY!

         WRITE(LUNREGPIPE,IOSTAT=IRTFLG) RVAL,CHAR(10)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'PIPING REGISTER',IREG)
            RETURN
         ENDIF

C        THIS LINE_FEED IS NECESSARY, DO NOT ASK ME WHY!
C        WRITE(LUNREGPIPE,IOSTAT=IRTFLG) CHAR(10)

#ifdef __linux__
         CALL FLUSHFILE(LUNREGPIPE)
#endif
      ELSE
         CALL ERRT(102,'NO PIPE OPEN ON LUNREGPIPE',LUNREGPIPE)
         IRTFLG = 1
      ENDIF
 
      RETURN
      END

C--*************************** GETNEXTVAR *****************************

      SUBROUTINE GETNEXTVAR(STRING,IFIRST,IGO,IEND,IRTFLG)

C     COMMON NEEDED FOR WRITE
      INCLUDE 'CMBLOCK.INC'

C     VAR DELIMITERS ARE []

      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=1) :: CTEMP

C     SET DEFAULT RETURN VALUES
      IGO       = 0
      IEND      = 0
      IRTFLG    = 0

C     FIND LAST CHAR POSITION IN STRING
      ILAST = LEN(STRING)

      DO I = IFIRST,ILAST
         CTEMP = STRING(I:I)

         IF (IGO == 0 .AND. CTEMP == '[') THEN
C           ARE STARTING A VARIABLE, SET IGO 
            IGO  = I

         ELSEIF (CTEMP == ';' .AND. IGO > 0) THEN
C           START OF COMMENT BUT NO ENDING ] FOR VARIABLE
            WRITE(NOUT,*)'  *** NO ENDING ] FOR LAST VARIABLE IN: ',STRING
            CALL ERRT(101,'BAD VARIABLE SYNTAX',NE)
            IRTFLG = I
            EXIT

        ELSEIF (CTEMP ==  ']' .AND. IGO > 0) THEN
C           ARE ENDING A VAR.
            IEND = I
            EXIT

         ELSEIF (CTEMP == ';') THEN
C           START OF COMMENT, THIS IS LINE END
            EXIT

         ELSEIF (IGO > 0 .AND. CTEMP < CHAR(31)) THEN
C           ILLEGAL CHAR INSIDE A REG. VAR. NAME
            WRITE(NOUT,*)'  *** ILLEGAL CHAR.: ',CTEMP,' AT POSITION: ',
     &           I,' IN REGISTER VAR. STRING: ',STRING
            CALL ERRT(101,'BAD CHAR. IN REGISTER VAR. NAME',NE)
            IRTFLG = I

         ELSEIF (CTEMP < CHAR(31)) THEN
C           ILLEGAL CHAR OUTSIDE A REG. VAR. NAME
            WRITE(NOUT,*)'  *** ILLEGAL CHAR.: ',CTEMP,' AT POSITION: ',
     &           I,' IN STRING: ',STRING
         ENDIF
      ENDDO

C     VARIABLE FOUND OR RAN OFF END OF STRING

      RETURN
      END





      SUBROUTINE GETNEXTVARbug(STRING,IFIRST,IGO,IEND,IRTFLG)

C     COMMON NEEDED FOR WRITE
      INCLUDE 'CMBLOCK.INC'

C     VAR DELIMITERS ARE []

      CHARACTER(LEN=*) :: STRING
      CHARACTER(LEN=1) :: CTEMP

C     SET DEFAULT RETURN VALUES
      IGO       = 0
      IEND      = 0
      IRTFLG    = 0

       write(NOUT,*)'  STRING: ',STRING

C     FIND LAST CHAR POSITION IN STRING
      ILAST = LEN(STRING)
       write(NOUT,*)'  STRING: ',STRING
       write(NOUT,*)'  IFIRST,ILAST: ',IFIRST,ILAST

      DO I = IFIRST,ILAST
         CTEMP = STRING(I:I)
       write(NOUT,*)'  CTEMP: ',CTEMP,' AT: ',I,' IGO: ',IGO

         IF (IGO == 0 .AND. CTEMP == '[') THEN
C           ARE STARTING A VAR, SET IGO 
            IGO  = I

         ELSEIF (CTEMP == ';') THEN
C           START OF COMMENT, THIS IS LINE END
            EXIT

        ELSEIF (CTEMP ==  ']' .AND. IGO > 0) THEN
C           ARE ENDING A VAR.
            IEND = I
            EXIT

         ELSEIF (IGO > 0 .AND. CTEMP < CHAR(31)) THEN
C           ILLEGAL CHAR INSIDE A REG. VAR. NAME
            WRITE(NOUT,*)'  *** ILLEGAL CHAR.: ',CTEMP,' AT POSITION: ',
     &           I,' IN REGISTER VAR. STRING: ',STRING
            CALL ERRT(101,'BAD CHAR. IN REGISTER VAR. NAME',NE)
            IRTFLG = I

         ELSEIF (CTEMP < CHAR(31)) THEN
C           ILLEGAL CHAR OUTSIDE A REG. VAR. NAME
            WRITE(NOUT,*)'  *** ILLEGAL CHAR.: ',CTEMP,' AT POSITION: ',
     &           I,' IN STRING: ',STRING
         ENDIF
      ENDDO

C     VARIABLE FOUND OR RAN OFF END OF STRING

      RETURN
      END



C++*********************************************************************
C
C  GETREGVAR.F                NEW JUNE 2002 ARDEAN LEITH
C
C **********************************************************************
C * SPIDER - MODULAR IMAGE PROCESSING SYSTEM.    AUTHOR: J.FRANK       *
C * COPYRIGHT (C)1985, 2002. HEALTH RESEARCH INCORPORATED (HRI),       *
C * ONE UNIVERSITY PLACE, RENSSELAER, NY 12144-3455.                   *
C * THE CONTENTS OF THIS DOCUMENT ARE PROPRIETARY TO HRI AND ARE NOT   *
C * TO BE DISCLOSED TO OTHERS OR USED FOR PURPOSES OTHER THAN INTENDED *
C * WITHOUT WRITTEN APPROVAL OF HRI.                                   *
C **********************************************************************
C
C  GETREGVAR(QSTRQ,QFIND,QEND,IGO,IEND,IRTFLG)
C
C  PARAMETERS:  
C               IRTFLG                                            (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE GETREGVAR(QSTRQ,QFIND,QEND,IGO,IEND,IRTFLG)

      CHARACTER(LEN=*), INTENT(IN)  ::      QSTRQ,QFIND,QEND
      INTEGER, INTENT(OUT)  ::              IGO,IEND,IRTFLG

      IRTFLG = 1

C     FIND LENGTH OF SEARCH STRING
      LENFIND  = LEN(QFIND)
       
C     FIND STARTING LOCATION OF SEARCH STRING IN QSTRQ
      IGO = INDEX(QSTRQ,QFIND)

C     RETURN IF NO SEARCH STRING IN QSTRQ
      IF (IGO <= 0) THEN
C          THIS WAS DUE TO XLRF90 MAC OPTIMIZER BUG
           IRTFLG = 1
           RETURN
      ENDIF

C     FIND START OF ASSOCIATED VALUE 
      IGO    = IGO + LENFIND

C     FIND END OF ASSOCIATED VALUE
      IEND   = IGO + INDEX(QSTRQ(IGO:),QEND) - 2 

C     RETURN IF NO ASSOCIATED VALUE OR QEND
      IF (IEND < IGO) RETURN

      IRTFLG = 0

      RETURN
      END




