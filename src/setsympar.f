

C++*********************************************************************
C
C  SETSYMPAR.F -- CREATED 2/20/01 FROM SPIDER.F     ARDEAN LEITH 
C                 ?<ANYTHING? BUG              FEB 2002 ARDEAN LEITH
C                 REWRITTEN                    JUN 2002 ARDEAN LEITH
C                 [] DEFAULT FOR VARIABLES     OCT 2005 ARDEAN LEITH
C                 GLOBAL DUPLICATES OK         FEB 2006 ARDEAN LEITH
C                 CVAR                         OCT 2006 ARDEAN LEITH
C                 NCHAR IN SYMPAR_SUB          NOV 2009 ARDEAN LEITH
C                 DONOTRECURSE IN SYMPAR_SUB   MAY 2010 ARDEAN LEITH
C                 ISSYMPAR IRTFLG BUG          MAY 2010 ARDEAN LEITH
C                 UNDEFINED VARS FIXED         DEC 2010 ARDEAN LEITH
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
C    SETSYMPAR(SYMPARID,SYMPARVAL,LOCAL,IRTFLG)
C
C    PURPOSE:       RECORD VARIABLE IN CVAR ARRAY
C
C    PARAMETERS:    SYMPARID   ID (WITH <>)                    (SENT)
C                   SYMPARVAL  VALUE                           (SENT)
C                   LOCAL      FLAG FOR LOCAL                  (SENT)
C                   IRTFLG     RETURN FLAG (0 IS NORMAL)   (RETURNED)
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************



      MODULE SYMPAR_STUFF

         SAVE

C        DANGER, MAXPRC IS ALSO SET IN spider.f!!
         INTEGER, PARAMETER :: MAXPRC = 20      ! NO. OF PROC. LEVELS
         INTEGER, PARAMETER :: MAXVARNAME_LEN  = 80

C        BANK ZERO IS FOR GLOBAL BANK
         INTEGER, DIMENSION(0:MAXPRC) :: IGOCSTRQ,   IENDCSTRQ
         INTEGER, DIMENSION(0:MAXPRC) :: IGOCVARNUM, IENDCVARNUM
         INTEGER, DIMENSION(0:MAXPRC) :: LENCSTRQ,   NCVAR

C        THIS SHOULD BE RE-DONE WITH ALLOCATABLE CHAR. ARRAY??
         INTEGER, PARAMETER         :: MAXCSTRQG     = 6000
         INTEGER, PARAMETER         :: MAXCSTRQ      = 16000

         INTEGER, PARAMETER         :: NUMCVARG_ORIG = 2000
         INTEGER, PARAMETER         :: NUMCVAR_ORIG  = 5300

         INTEGER                    :: NUMCVAR       = 0
         INTEGER                    :: NUMCVARG      = 0

         CHARACTER(LEN=MAXCSTRQ+MAXCSTRQG) :: CSTRQ  = ' '

         INTEGER, PARAMETER         :: MAXCVAR_LEN = 120
         CHARACTER(LEN=MAXCVAR_LEN), ALLOCATABLE, 
     &                   DIMENSION(:) :: CVARVALUES 

      END MODULE SYMPAR_STUFF

C------------------------- END MODULE ---------------------------------

C CONTENTS:
C       ISSYMPAR(NAME,IBANK,ICVAR,IRTFLG)
C       SETSYMPAR(SYMPARID,SYMPARVAL,LOCAL,IRTFLG)
C       SYMPAR_NEW(IBANK,NAME,CVALUE,ICVAR,IRTFLG)
C       SYMPAR_OLD(ICVAR,CVALUE,IRTFLG)
C       SYMPAR_INIT(IBANK,IRTFLG)
C       SYMPAR_FIND(IBANK,NAME,ICVAR,IRTFLG)
C       SYMPAR_SUB(INPUT,OUTPUT,NCHAR,ILEVEL,IRTFLG)
C       SYMPAR_REINIT(IRTFLG)
C       SYMPAR_GET_NUMS(ICVARS,NCHAR)
C       ASSOCARRAY(QSTRQ,QFIND,IGO,IEND,IRTFLG)


C------------------------- SETSYMPAR ---------------------------------


      SUBROUTINE SETSYMPAR(SYMPARID,SYMPARVAL,LOCAL,IRTFLG)

      COMMON /QSTR_STUFF1/ ISTOP,NSTDOUTP,NSTDINP,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      CHARACTER(LEN=*)  :: SYMPARID,SYMPARVAL
      LOGICAL           :: LOCAL
      CHARACTER(LEN=1)  :: NULL

      IBANK = 0
      IF (LOCAL) IBANK = ISTOP

      NULL = CHAR(0)
      IF (SYMPARID .EQ. NULL) THEN
C        SET INITIAL VARIABLE INFO (BEST MIGRATED INTO SPIDER CALL)
	 CALL SYMPAR_INIT(IBANK,IRTFLG)
         IRTFLG = 0
         RETURN
      ENDIF

C     SEARCH FOR AN EXISTING VARIABLE OF SAME NAME
      CALL SYMPAR_FIND(IBANK,SYMPARID,ICVAR,IRTFLG)

      IF (ICVAR .LE. 0) THEN
C        CREATE NEW VARIABLE 
         CALL SYMPAR_NEW(IBANK,SYMPARID,SYMPARVAL,ICVAR,IRTFLG)

      ELSE
C        REPLACE VARIABLE VALUE 
         CALL SYMPAR_OLD(ICVAR,SYMPARVAL,IRTFLG)
      ENDIF

      RETURN
      END

C++*********************************************************************
C
C SYMPAR_NEW                         NEW            OCT 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    SYMPAR_NEW(IBANK,NAME,CVALUE,ICVAR,IRTFLG)
C 
C    PURPOSE:     CREATES A NEW CVAR VARIABLE & ASSOCIATED VALUE   
C
C    PARAMETERS:  IBANK     STACK LEVEL                         (SENT)
C                 NAME      VARIABLE NAME, WITH [...]           (SENT)
C                 CVALUE    VARIABLE VALUE                      (SENT)
C                 ICVAR     VARIABLE NUMBER                     (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C     NOTES: VARIABLE SHOULD NOT ALREADY EXIST (NOT CHECKED)
C
C--*******************************************************************

      SUBROUTINE SYMPAR_NEW(IBANK,NAME,CVALUE,ICVAR,IRTFLG)

      USE SYMPAR_STUFF

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)  :: NAME,CVALUE
      CHARACTER(LEN=20) :: CICVAR

C     NOTE: IF VARIABLE ALREADY EXISTS THIS GIVES ERROR
      IRTFLG = 1

C     VARIABLE DOES NOT EXIST, CREATE IT
C     WRITE(NDAT,*)' CREATING VARIABLE: ',NAME

C     INCREMENT VARIABLE NUMBER
      ICVAR = IENDCVARNUM(IBANK) + 1

      IF (IBANK .EQ. 0 .AND. ICVAR .GT. NCVAR(IBANK)) THEN
C        OVER-RUN OF GLOBAL VARIABLE VALUE ARRAY
         IT = NCVAR(IBANK)
         CALL ERRT(102,'TOO MANY GLOBAL VARIABLES REQUESTED, LIMIT',IT)
         RETURN

      ELSEIF (ICVAR .GT. NCVAR(IBANK)) THEN
C        OVER-RUN OF CVAR VARIABLE VALUE ARRAY

         write(6,*) ' NEW VARIABLE NAME:  ',NAME
         write(6,*) ' BANK:               ',IBANK
         write(6,*) ' IGOCSTRQ(IBANK):    ',IGOCSTRQ(IBANK)
         write(6,*) ' IGOCVARNUM(IBANK):  ',IGOCVARNUM(IBANK)
         write(6,*) ' IENDCVARNUM(IBANK): ',IENDCVARNUM(IBANK)
         write(6,*) 'CSTRQ(IBANK)(1:60): ',
     &          CSTRQ(IGOCSTRQ(IBANK):IGOCSTRQ(IBANK)+60)

         IT = NCVAR(IBANK)
         CALL ERRT(102,'TOO MANY LOCAL VARIABLES REQUESTED, LIMIT',IT)
         RETURN
      ENDIF

C     PLACE VARIABLE NAME STRING IN CSTRQ ARRAY
      CALL INTTOCHAR(ICVAR,CICVAR,LENNUM,1)  ! LENGTH OF VARIABLE NUM.
      LENVAR  = LEN(CVALUE)         ! LENGTH OF VARIABLE VALUE
      LENNAM  = LEN(NAME)           ! LENGTH OF VARIABLE NAME

      IF (LENNAM .GT. 80) THEN
         WRITE(NOUT,*) '  VARIABLE NAME: ',NAME
         CALL ERRT(102, 'OVERLENGTH VARIABLE NAME, CHARS.',LENNAM)
         RETURN
      ELSEIF (LENVAR .GE. MAXCVAR_LEN) THEN
         WRITE(NOUT,*) '  VARIABLE NAME: ',NAME
         CALL ERRT(102,'OVERLENGTH VARIABLE, CHARS.',LENVAR)
         RETURN
      ENDIF

      LENADD  = LENNAM + LENNUM          ! LENGTH OF CSTRQ ENTRY
      IGOQ    = IENDCSTRQ(IBANK) + 1
      IENDQAT = IENDCSTRQ(IBANK) + LENADD 

      IF (IENDQAT .GE. LENCSTRQ(IBANK)) THEN
C        OVER-RUN OF CSTRQG ARRAY
         CALL ERRT(102,'CSTRQG OVERFLOW, TOO MANY VARIABLES',IENDQAT)
         RETURN
      ENDIF

#ifdef DEBUGD
      write(6,*) ' BANK:               ',IBANK
      write(6,*) ' IENDCSTRQ(IBANK):  ',IENDCSTRQ(IBANK)
      write(6,*) ' LENVAR:  ',LENVAR
      write(6,*) ' LENADD:  ',LENADD
      write(6,*) ' LENNAM:  ',LENNAM
      write(6,*) ' LENNUM:  ',LENNUM
      write(6,*) ' NEW VARIABLE NAME:  ',NAME(2:LENNAM-1)
      write(6,*) ' NEW VARIABLE VALUE:  ',CVALUE
      write(6,*) ' IGOQ:  ',IGOQ,'   IENDQAT:  ',IENDQAT
      write(6,*) ' CICVAR(1:LENNUM):  ',CICVAR(1:LENNUM)
      write(6,*) ' ICVAR:  ',ICVAR
#endif

      CSTRQ(IGOQ:IGOQ+LENADD-1) = NAME(2:LENNAM-1) // '>' // 
     &                            CICVAR(1:LENNUM) // '<' 
      CVARVALUES(ICVAR) = CVALUE(1:LENVAR) 

#ifdef DEBUGD
      write(6,*) ' CSTRQ(1:50): ',CSTRQ(1:50)
      write(6,*) ' CSTRQ(6000:6050): ',CSTRQ(6000:6050)
      write(6,*) ' CVARVALUES(ICVAR): ',CVARVALUES(ICVAR)(1:LENVAR)
#endif

      IENDCSTRQ(IBANK)   = IENDQAT 
      IENDCVARNUM(IBANK) = ICVAR
      IRTFLG             = 0

      RETURN
      END

C++*********************************************************************
C
C SYMPAR_OLD                         NEW            OCT 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    SYMPAR_OLD(ICVAR,CVALUE,IRTFLG)
C 
C    PURPOSE:     REPLACES EXISTING CVAR VARIABLE'S ASSOCIATED VALUE   
C
C    PARAMETERS:  ICVAR     VARIABLE NUMBER                     (SENT)
C                 CVALUE    VARIABLE VALUE                      (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C     NOTES: VARIABLE SHOULD NOT ALREADY EXIST (NOT CHECKED)
C
C--*******************************************************************

      SUBROUTINE SYMPAR_OLD(ICVAR,CVALUE,IRTFLG)

      USE SYMPAR_STUFF
         
      CHARACTER(LEN=*) :: CVALUE 

C     IF VARIABLE NOT ALREADY EXISTS (ERROR)
      CVARVALUES(ICVAR) = CVALUE

      IRTFLG = 0

      END

C++*********************************************************************
C
C SYMPAR_INIT                         NEW            OCT 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    SYMPAR_INIT(IBANK,IRTFLG)
C 
C    PURPOSE:     INITIALIZES CVAR VARIABLE'S    
C
C    PARAMETERS:  IBANK     PROCEDURE BANK NUMBER               (SENT)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C
C--*******************************************************************

      SUBROUTINE SYMPAR_INIT(IBANK,IRTFLG)

      USE SYMPAR_STUFF

      IRTFLG = 1
      IF (IBANK .LT. 0) THEN
         CALL ERRT(102,'ILLEGAL VARIABLE BANK:',IBANK)
         RETURN
      ELSEIF (IBANK .GT. MAXPRC) THEN
         IT     = MAXPRC
         CALL ERRT(102,'VARIABLE BANK EXCEEDS MAXPRC',IT)
         RETURN
      ENDIF

      IF (NUMCVAR .LE. 0) THEN
C        CREATE THE CVAR STORAGE ARRAY (ONLY OCCURS ONCE)
         MWANT = NUMCVAR_ORIG + NUMCVARG_ORIG
         ALLOCATE (CVARVALUES(MWANT), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'UNABLE TO ALLOCATE CVARS:',MWANT)
            RETURN
         ENDIF
         NUMCVAR  = NUMCVAR_ORIG
         NUMCVARG = NUMCVARG_ORIG

         NCVAR(0)           = NUMCVARG
         NCVAR(1:MAXPRC)    = MWANT

         LENCSTRQ(0)        = MAXCSTRQ
         LENCSTRQ(1:MAXPRC) = MAXCSTRQ + MAXCSTRQG
      ENDIF

      IF (IBANK .EQ. 0) THEN
C        INITIALIZE GLOBAL BANK ZERO, (SHOULD ONLY OCCUR ONCE)

         IGOCSTRQ(IBANK)    = 1
         IGOCVARNUM(IBANK)  = 1
         IENDCVARNUM(IBANK) = 0

      ELSEIF (IBANK .EQ. 1) THEN
         IGOCSTRQ(IBANK)    = MAXCSTRQG         + 1
         IGOCVARNUM(IBANK)  = NUMCVARG          + 1
         IENDCVARNUM(IBANK) = IGOCVARNUM(IBANK) - 1
      ELSE
         IGOCSTRQ(IBANK)    = IENDCSTRQ(IBANK-1)   + 1
         IGOCVARNUM(IBANK)  = IENDCVARNUM(IBANK-1) + 1
         IENDCVARNUM(IBANK) = IGOCVARNUM(IBANK)    - 1
      ENDIF

C     PUT INITIAL < IN CSTRQ
      CSTRQ(IGOCSTRQ(IBANK):IGOCSTRQ(IBANK)) = '<'  
      IENDCSTRQ(IBANK) = IGOCSTRQ(IBANK) 

#ifdef DEBUGD
      write(6,*) 'initial IBANK:',IBANK 
      write(6,*) 'initial IENDCSTRQ(IBANK):',IENDCSTRQ(IBANK) 
#endif
      IRTFLG = 0
 
      RETURN
      END

C++*********************************************************************
C
C ISSYMPAR                          NEW            JUN 2009 ARDEAN LEITH
C
C **********************************************************************
C
C    ISSYMPAR(NAME,IBANK,ICVAR,IRTFLG)
C 
C    PURPOSE:     SEES IF A CVAR VARIABLE ALREADY EXISTS   
C
C    PARAMETERS:  IBANK     STACK LEVEL                         (SENT)
C                 NAME      VARIABLE NAME, WITH [...]           (SENT)
C                 ICVAR     VARIABLE NUMBER                     (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C
C--*******************************************************************

      SUBROUTINE ISSYMPAR(NAME,IBANK,ICVAR,IRTFLG)

      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER, INTENT(IN)           :: IBANK
      INTEGER, INTENT(OUT)          :: ICVAR,IRTFLG
 
C     FOR ISTOP 
      INTEGER, DIMENSION(MAXPRC) :: IPSTACK,IPNUMSTACK,IPARNUM
      COMMON /QSTR_STUFF1/ ISTOP,ITI,ITIN,IWHERE,IPSTACK,
     &                     IPNUMSTACK,IPARNUM

      ILEVEL = IBANK
      IF (ILEVEL .LT. 0) ILEVEL = ISTOP

      CALL SYMPAR_FIND(ILEVEL,NAME,ICVAR,IRTFLG)
      IF (IRTFLG .EQ. 1 .OR. ICVAR .LE. 0) THEN
C        NOT AT THIS BANK, TRY GLOBAL
         IF (ILEVEL .NE. 0) THEN
            ILEVEL = 0
            CALL SYMPAR_FIND(ILEVEL,NAME,ICVAR,IRTFLG)
         ENDIF
      ENDIF
 
      END

C++*********************************************************************
C
C SYMPAR_FIND                       NEW            OCT 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    SYMPAR_FIND(IBANKT,NAME,ICVAR,IRTFLG)
C
C    PURPOSE:     FINDS IF VARIABLE EXISTS, RETURNS VARIABLE NUMBER  
C
C    PARAMETERS:  IBANK     STACK LEVEL () IS CURRENT TOP       (SENT)
C                 NAME      VARIABLE NAME, WITH [...]           (SENT)
C                 ICVAR     VARIABLE NUMBER                     (RET.)
C                 IRTFLG    ERROR FLAG                          (RET.)
C
C--*******************************************************************

      SUBROUTINE SYMPAR_FIND(IBANK,NAME,ICVAR,IRTFLG)

      USE SYMPAR_STUFF

      CHARACTER(LEN=*)                :: NAME
      CHARACTER(LEN=MAXVARNAME_LEN+2) :: NAMET

C     SEE IF VARIABLE EXISTS YET

C     <> ARE USED AS VARIABLE ID DELIMITERS IN CSTRQ
      LENT  = LEN(NAME)
      NAMET = '<' // NAME(2:LENT-1) // '>' 

      IGOQ   = IGOCSTRQ(IBANK)
      IENDQ  = IENDCSTRQ(IBANK)
      IRTFLG = 0
      ICVAR  = 0

      IF (IENDQ .GT. IGOQ) THEN
C        SEARCH FOR THIS VARIABLE AT THIS LEVEL
         CALL ASSOCARRAY(CSTRQ(IGOQ:IENDQ),NAMET(1:LENT),
     &                   IGO,IEND,IRTFLG)

#ifdef DEBUGD
         write(6,*) ' HUNTING, BANK:',IBANK,':', NAMET(1:LENT)
         write(6,*) ' IN: ',CSTRQ(IGOQ:IENDQ)
#endif

         IF (IRTFLG .EQ. 0) THEN
C           VARIABLE EXISTS, FIND VARIABLE NUMBER
#ifdef DEBUGD
            write(6,*) ' FOUND VAR NUMBER: ', CSTRQ(IGO:IEND)
            write(6,*) ' AT CSTRQ: ',IGO,'..',IEND
#endif
            READ(CSTRQ(IGOQ+IGO-1:IGOQ+IEND-1),*) ICVAR
            LENVAR = lnblnk(CVARVALUES(ICVAR))

#ifdef DEBUGD
            !write(6,*) ' INTEGERIZED VAR. NUM: ', ICVAR
#endif

#ifdef DEBUGD
        ELSE
            write(6,*) ' SYMPAR_FIND, IRTFLG: ', IRTFLG
#endif
        ENDIF


      ENDIF

      RETURN
      END



C++*************************************************************************
C
C  SYMPAR_SUB.F   CREATED  FROM SPIDER.F         Sep 2000  ARDEAN LEITH 
C                 MULTIPLE VARIABLE SUBSTITUTION JAN 2001  ARDEAN LEITH
C                 ALPHABETICAL VARIABLES         JUN 2002  ARDEAN LEITH
C                 NESTED VARIABLES               SEP 2002  ARDEAN LEITH
C                 [] DEFAULT FOR VARIABLES       OCT 2005  ARDEAN LEITH
C                 CVAR                           OCT 2006  ARDEAN LEITH
C                 DONOTRECURSE                   JUN 2010  ARDEAN LEITH
C              
C **********************************************************************
C
C    SYMPAR_SUB(INPUT,OUTPUT,NCHAR,ILEVELT,DONOTRECURSE,IRTFLG)
C
C    PURPOSE:       RUN-TIME VARIABLE SUBSTITUTION FOR ALL [ID] IN
C                   INPUT STRING AT THIS STACK LEVEL
C
C    PARAMETERS:    INPUT     INPUT LINE CONTAINING [STRING..] (SENT)
C                   OUTPUT    SUBSTITUTED OUTPUT LINE         (RET.)
C                   NCHAR     LAST NON_BLANK CHAR BEFORE ;    (RET.)
C                   ILEVEL    NESTING LEVEL                   (SENT)
C                   IRTFLG    RETURN FLAG (0 IS NORMAL)       (RETURNED)
C   
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE SYMPAR_SUB(INPUT,OUTPUT,NCHAR,ILEVEL,
     &                      DONOTRECURSE,IRTFLG)

      USE SYMPAR_STUFF

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=*)           :: INPUT,OUTPUT
      CHARACTER(LEN=MAXCVAR_LEN) :: CVALUE
      LOGICAL                    :: LOCAL
      INTEGER, DIMENSION(160)    :: DONOT
      LOGICAL                    :: DONOTRECURSE

      INOT   = 0     
      LENT   = LEN(INPUT)
      OUTPUT = INPUT
      NCHAR  = lnblnk(INPUT)

      IDONE  = 0
      IP2N   = INDEX(OUTPUT,']',BACK=.TRUE.)

      DO
C        REPLACE ALL SYMBOL [] SETS IN POSITIONS 1...LENT

C        FIND START AND END OF FIRST VARIABLE STRING IN OUTPUT 
         CALL CHARINSIDE(OUTPUT(1:LENT),'[',']',.FALSE.,.FALSE.,
     &                  IP1,IP2,NCT)

C        GET OUT OF LOOP IF NO SYMBOL SET [] FOUND
         IF (NCT .LE. 0) EXIT

C        FIND SYMBOL STRING FOR SUBSTITUTION.

         IRTFLG = 1

C        <> ARE USED AS SYMBOL ID DELIMITERS IN CSTRQ
         OUTPUT(IP1:IP1) = '<'
         OUTPUT(IP2:IP2) = '>'

         IGOQ            = IGOCSTRQ(ILEVEL)
         IENDQ           = IENDCSTRQ(ILEVEL)

         IF (IENDQ .GT. IGOQ) THEN
C           SEARCH FOR THIS LOCAL VARIABLE
            IGOQ   = IGOCSTRQ(ILEVEL)
            IENDQ  = IENDCSTRQ(ILEVEL)
            CALL SYMPAR_FIND(ILEVEL,OUTPUT(IP1:IP2),ICVAR,IRTFLG)

            IF (ICVAR .GT. 0) THEN
C              VARIABLE WAS FOUND
#ifdef DEBUGD
               write(6,*) ' ILEVEL: ',ILEVEL
               write(6,*) ' OUTPUT(IP1:IP2) :',OUTPUT(IP1:IP2)
               write(6,*) ' IENDCSTRQ(s):',(IENDCSTRQ(I),I=1,ILEVEL)
               write(6,*) ' ICVAR: ',ICVAR
#endif
               CVALUE   = CVARVALUES(ICVAR)
               LENVALUE = lnblnk(CVALUE)

#ifdef DEBUGD
               write(6,*) ' CVALUE: ',CVALUE(:LENVALUE)
               write(6,*) ' CSTRQ(1:50): ',CSTRQ(1:50)
#endif

C              COPY CORRESPONDING CSTRQ STRING TO OUTPUT & UPDATE NCHAR
               CALL SUBCHAR(CVALUE(:LENVALUE),OUTPUT,
     &                      IP1,IP2,NCHAR,IRTFLG)
            ENDIF
         ENDIF

         IF (IRTFLG .NE. 0 .AND. IENDCSTRQ(0) .GT. 0) THEN
C           SEARCH FOR A GLOBAL VARIABLE
            IGOQ   = IGOCSTRQ(ILEVEL)
            IENDQ  = IENDCSTRQ(ILEVEL)
            CALL SYMPAR_FIND(0,OUTPUT(IP1:IP2),ICVAR,IRTFLG)

            IF (ICVAR .GT. 0) THEN
C              FOUND A GLOBAL VARIABLE
#ifdef DEBUGD
               write(6,*) ' ILEVEL: ',ILEVEL
               write(6,*) ' OUTPUT(IP1:IP2) :',OUTPUT(IP1:IP2)
               write(6,*) ' IENDCSTRQ(s):',(IENDCSTRQ(I),I=1,ILEVEL)
               write(6,*) ' ICVAR: ',ICVAR
#endif

               CVALUE   = CVARVALUES(ICVAR)
               LENVALUE = lnblnk(CVALUE)

#ifdef DEBUGD
               write(6,*) ' CVALUE: ',CVALUE(:LENVALUE)
               write(6,*) ' CSTRQ(1:50): ',CSTRQ(1:50)
#endif
C              COPY CORRESPONDING CSTRQ STRING TO OUTPUT & UPDATE NCHAR
               CALL SUBCHAR(CVALUE(:LENVALUE),OUTPUT,
     &                      IP1,IP2,NCHAR,IRTFLG)
            ENDIF
        ENDIF

C        END OF SUBSTITUTION
         IP2N = IP2N + (IENDQ - IGOQ) - (IP2 - IP1)

         IF (ICVAR .EQ. 0) THEN
C           NO SUBSTITUTION, PROBABLY A REGISTER VARIABLE
            INOT          = INOT + 1
            DONOT(INOT)   = IP1
            INOT          = INOT + 1
            DONOT(INOT)   = IP2
         ENDIF

C        NEXT SEARCH IS OVER WHOLE STRING
         LENT  = NCHAR
         IDONE = IDONE + 1
         IF (IDONE .GT. 100000) THEN
            CALL ERRT(102,'RECURSIVE VARIABLE ???? LOOPS',IDONE)
            RETURN
         ENDIF
         IF (DONOTRECURSE) EXIT       ! NON-RECURSIVE SUBSTITUTION
      ENDDO

C     END SYMBOL SUBSTITUTION 

C     FIX <> BEFORE RETURNING
      IF (INOT .GT. 0) THEN
         DO IV=1,INOT,2
            I = DONOT(IV)
            IF (OUTPUT(I:I) .EQ. '<') THEN
               OUTPUT(I:I) = '['
            ELSE
               CALL ERRT(102,'LOCATION IS NOT <',I)
            ENDIF
            I = DONOT(IV+1)
            IF (OUTPUT(I:I) .EQ. '>') THEN
               OUTPUT(I:I) = ']'
            ELSE
               CALL ERRT(102,'LOCATION IS NOT >',I)
            ENDIF
         ENDDO
      ENDIF

C     SET NORMAL RETURN FLAG
      IRTFLG = 0

      END


C++*********************************************************************
C
C SYMPAR_REINIT                    NEW            OCT 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    SYMPAR_REINIT()
C
C    PURPOSE:     RESIZES   CVAR SPACE   
C
C    PARAMETERS:  IRTFLG    ERROR FLAG                          (RET.)
C
C    YES, I KNOW THAT IT SHOULD BE WRITTEN USING POINTERS BUT
C    I DOUBT ANYONE WILL EVER USE THIS!! AL
C
C--*********************************************************************

      SUBROUTINE SYMPAR_REINIT(IRTFLG)

      USE SYMPAR_STUFF

      INCLUDE 'CMBLOCK.INC'
      CHARACTER(LEN=MAXCVAR_LEN), ALLOCATABLE, DIMENSION(:)::CVARVALUEST

      NCVAR     = NUMCVAR
      NCHARCVAR = MAXCSTRQ
      NCHAREACH = MAXCVAR_LEN

      WRITE(NOUT,90)NCVAR,NCHARCVAR,NCHAREACH  
90    FORMAT(' CURRENT NUMBER OF VARIABLES: ',I6,
     &       ' CHARACTERS FOR VARIABLE NAMES: ',I8,/,
     &       ' CHARACTERS / VARIABLE NAME: ',I4)
  
      CALL RDPRI1S(NCVARN,NOT_USED,
     &            'NUMBER OF VARIABLES WANTED',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (NCVARN .GT. NCVAR(MAXPRC)) THEN
         ALLOCATE (CVARVALUEST(NCVARN), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'UNABLE TO INCREASE VARIABLES:',NCVAR)
            RETURN
         ENDIF
              
         CVARVALUEST(1:NUMCVAR) = CVARVALUES
         DEALLOCATE(CVARVALUES)

         ALLOCATE (CVARVALUES(NCVARN), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'UNABLE TO INCREASE VARIABLES:',NCVARN)
            RETURN
         ENDIF

         CVARVALUES(1:NUMCVAR) = CVARVALUEST
         DEALLOCATE(CVARVALUEST)
         
         NUMCVAR = NCVARN
      ENDIF

      IRTFLG = 0
 
      RETURN
      END

C++*********************************************************************
C
C SYMPAR_GET_NUMS                    NEW            AUG 2006 ARDEAN LEITH
C
C **********************************************************************
C
C    SYMPAR_GET_NUMS(ICVARS,NCHAR)
C
C    PURPOSE:     RETRIEVE VAR INFO   
C
C    PARAMETERS:  ICVARS                              (RET.)
C                 NCHAR                               (RET.)
C
C
C--*******************************************************************

       SUBROUTINE SYMPAR_GET_NUMS(ICVARS,NCHAR)

       USE SYMPAR_STUFF

       ICVARS =  NUMCVARG  + NUMCVAR
       NCHAR  =  MAXCSTRQG + MAXCSTRQ

       END

C++*********************************************************************
C
C  ASSOCARRAY.F                NEW JUNE 2002               ARDEAN LEITH
C
C **********************************************************************
C
C  ASSOCARRAY(QSTRQ,QFIND,IGO,IEND,IRTFLG)
C
C  PURPOSE:     RETURN POSITION OF ASSOCIATED VARIABLE FOR QFIND IN QSTRQ
C               QFIND WILL BE SURROUNDED IN <>
C               QFIND WILL BE FOLLWED BY ASSOCIATED VARIABLE AND < 
C
C  PARAMETERS:  
C               QSTRQ                                             (SENT)
C               QFIND                                             (SENT)
C               IGO                                               (RET.)
C               IEND                                              (RET.)
C               IRTFLG                                            (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE ASSOCARRAY(QSTRQ,QFIND,IGO,IEND,IRTFLG)

      CHARACTER(LEN=*), INTENT(IN)  :: QSTRQ,QFIND
      INTEGER, INTENT(OUT)          :: IGO,IEND,IRTFLG
      CHARACTER(LEN=1)              :: TAGEND

      IRTFLG = 1

C     FIND LENGTH OF SEARCH STRING
      LENFIND  = LEN(QFIND)
       
C     FIND STARTING LOCATION OF SEARCH STRING IN QSTRQ
      IGO = INDEX(QSTRQ,QFIND)

C     RETURN IF NO SEARCH STRING IN QSTRQ  xlf90 compiler bug hack
      IF (IGO .LE. 0) THEN
          IRTFLG = 1
          RETURN
      ENDIF

C     FIND START OF ASSOCIATED VALUE 
      IGO = IGO + LENFIND

C     FIND END OF ASSOCIATED VALUE
      TAGEND = '<' 
      IEND   = IGO + INDEX(QSTRQ(IGO:),TAGEND) - 2 

C     RETURN IF NO ASSOCIATED VALUE xlf90 compiler bug hack
      IF (IEND .LT. IGO) THEN
          IRTFLG = 1
          RETURN
      ENDIF

      IRTFLG = 0

      RETURN
      END


