C++*********************************************************************
C
C  INQUIREIF1.F    NEW ROUTINE                      SEP 97 al
C                  F90 CHANGES                      MAR 98 al
C                  REMOVED ifdef sgi                NOV 01 al
C                  INDEXED STACK                    JAN 03 ARDEAN LEITH
C                  LUNRED                           FEB 03 ARDEAN LEITH
C                  NLET = 0 BUG                     FEB 05 ARDEAN LEITH
C                  OPENINLN KIND                    OCT 10 ARDEAN LEITH
C                  _4@2 NON  BUG                    JAN 11 ARDEAN LEITH
C                  TYPET ifort                      MAY 12 ARDEAN LEITH
C                  ==, nsam                         SEP 16 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C   INQUIREIF1(LUN,FILNAM,TYPET,EX,ISOPEN,LUNOP,INLNED,IMGNUM,IRTFLG)    
C
C   PURPOSE:     DETERMINES IF A FILE EXISTS
C
C   PARAMETERS:  LUN                                          (SENT)
C                FILNAM                                       (SENT)
C                TYPET   SEARCH TYPE                          (SENT)
C                EX                                           (RET.)
C                ISOPEN                                       (RET.)
C                LUNOP                                        (RET.)
C                INLNED                                       (RET.)
C                IMGNUM                                       (RET.)
C                IRTFLG                                       (RET.)
C
C--*********************************************************************

        SUBROUTINE INQUIREIF1(LUN,FILNAM,TYPET,EX,ISOPEN,LUNOP,
     &                        INLNED,IMGNUM,IRTFLG)

C       USE INLINE BUFFER COMMON AREA
        INCLUDE 'INLN_INFO.INC'

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC' 

        INTEGER                      :: LUN
        CHARACTER (LEN=*)            :: FILNAM
        CHARACTER (LEN=*)            :: TYPET
        LOGICAL                      :: EX,ISOPEN
        INTEGER                      :: LUNOP,INLNED,IMGNUM,IRTFLG

        CHARACTER (LEN=MAXNAM)       :: FILDUM
        CHARACTER (LEN=1)            :: NULL,FIRSTC
        LOGICAL                      :: STACKOPN

        INTEGER, PARAMETER           :: I_8 = SELECTED_INT_KIND(12)
        INTEGER(KIND=I_8), PARAMETER :: ZERO_8 = 0

        NULL    = CHAR(0)

        NLET    = lnblnk(FILNAM)
        FIRSTC  = FILNAM(1:1)
        ILOCAT  = INDEX(FILNAM,'@')
        NAMEND  = NLET
        IF (ILOCAT > 1) NAMEND = ILOCAT - 1
        EX      = .FALSE.
        ISOPEN  = .FALSE.
        INLNED  = 0
        IMGNUM  = 0
        IRTFLG  = 0

C       CHECK FOR ANONMOLOUS INPUT
        IF (NLET <= 0 .OR. NAMEND <= 0) RETURN

        IF (FIRSTC .NE. '_' .AND. ILOCAT <= 0) THEN
C          NO LEADING '_' AND NO '@' MEANS THAT IT IS A REGULAR 
C          FILE_BASED NON-STACK IMAGE OR OTHER FILE (SUCH AS A
C          DOCUMENT FILE)

#if defined (__INTEL_COMPILER)
           IF (TYPET == 'DIR') THEN
C             SEE IF THIS DIR EXISTS, (RETURNS EX)
              INQUIRE(DIRECTORY=FILNAM(1:NAMEND),EXIST=EX,ERR=999)
           ELSE
C             SEE IF THIS FILE EXISTS, (RETURNS EX, ISOPEN, LUNOP)
              INQUIRE(FILE=FILNAM(1:NAMEND),EXIST=EX,OPENED=ISOPEN,
     &             NUMBER=LUNOP,ERR=999)
           ENDIF
#else
C          SEE IF THIS FILE EXISTS, (RETURNS EX, ISOPEN, LUNOP)
           INQUIRE(FILE=FILNAM(1:NAMEND),EXIST=EX,OPENED=ISOPEN,
     &             NUMBER=LUNOP,ERR=999)
#endif

        ELSEIF (FIRSTC .NE. '_') THEN
C          CHECK TO SEE IF IMAGE EXISTS IN THIS EXISTING STACK
C          NO LEADING '_' MEANS THIS IS FILE_BASED STACK OR BARE STACK 

           IFOUND = -4   ! SET IFOUND TO DECREASE OPENING OUTPUT INFO
           MAXIM  =  1   ! A BARE STACK FILE IS OK

           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'Z',IFORM,NX,NY,NZ,
     &              MAXIM,' ',.TRUE.,IRTFLG)
           IFOUND = 0    ! RESET NEEDED, THIS IS IN COMMON /UNITS/

           EX     = (IRTFLG == 0)
           !write(6,*) ' in inquirif1,  irtflg,ex:',irtflg,ex

        ELSEIF (FIRSTC == '_') THEN
C          INLINE IMAGE FILE OR OVERALL INLINE STACK ACCESS WANTED

C          RETRIVE INLINE BUFFER NUMBER FROM FILE NAME
           CALL INLNBUF(FILNAM,NLET,INLNED,IRTFLG)
           IF (IRTFLG .NE. 0)  RETURN

C          SEE IF INLINE STACK EXISTS NOW
           STACKOPN = (NSAMBUF(INLNED) > 0)
           IF (.NOT. STACKOPN) RETURN

           IF (ILOCAT == 0) THEN
C              SIMPLE INLINE IMAGE OR OVERALL INLINE STACK
               EX = .TRUE.
               RETURN
           ENDIF

C          FIND IMAGE NUMBER WITHIN STACK 
C          READ(FILNAM(ILOCAT+1:),*,IOSTAT=IER) IMGNUM -changed for osf-liy
           CALL FILCAD(FILNAM(ILOCAT:),FILDUM,IMGNUM,IER)
           IF (IER .NE. 0) THEN
              CALL ERRT(101,'UNABLE TO INTERPRET IMAGE NUMBER',NE)
              RETURN
           ENDIF
                 
           IFOUND = -4   ! SET IFOUND TO DECREASE OPENING OUTPUT INFO
           MAXIM  =  1   ! A BARE STACK FILE IS OK

           CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'Z',IFORM,NX,NY,NZ,
     &              MAXIM,' ',.TRUE.,IRTFLG)
           IFOUND = 0    ! RESET NEEDED, THIS IS IN COMMON /UNITS/

           EX     = (IRTFLG == 0)

        ENDIF

        CLOSE(LUN)
        RETURN


999     WRITE(NOUT,*)'*** ERROR INQUIRING ABOUT FILE: ',FILNAM(1:NLET)
        CALL ERRT(100,' ',NE)
        EX = .FALSE.
        RETURN

        END





