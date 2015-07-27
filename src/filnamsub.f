
C++*********************************************************************
C
C    FILNAMSUB.F                             NEW  JUL 1997 ArDean Leith
C                  CHARINSIDE PARAMETERS CHANGED  JAN 2001 ArDean Leith
C                  SGI LEAK ON INTERNAL FMT       AUG 2002 ArDean Leith
C                  BETTER ERROR MSG               SEP 2002 ArDean Leith
C                  ADDED {%f3.4%x11} SUPPORT      JAN 2004 ArDean Leith
C                  ***[] SUPPORT FROM SUBSYM      DEC 2005 ArDean Leith
C                  ***/xxx***{v1][v2] SUPPORT     JAN 2007 ArDean Leith
C                  TOLERATES abc{***} TEMPLATE    MAY 2007 ArDean Leith
C                  %f% val < 0 bug                MAY 2012 ArDean Leith
C                  %I% SUPPORTED                  MAY 2012 ArDean Leith
C                  %I% BUG (not %i%               APR 2013 ArDean Leith
C                  aa**bb**[v1][v2] BUG           OCT 2013 ArDean Leith
C                  ff***[r1]  [ss] BUG            OCT 2013 ArDean Leith
C                  1-[r1]  [ss] BUG               NOV 2013 ArDean Leith
C                  igo blank format BUG           JUL 2015 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2015  Health Research Inc.,                         *
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
C    FILNAMSUB(FILNAM,NLET,IBANK,IRTFLG)
C
C    PARAMETERS: 
C           FILNAM     CHARACTER STRING FOR SUBSTITUION  (SENT/RETURNED)
C           NLET       NO. OF CHARACTERS IN STRING       (SENT/RETURNED)
C           IBANK      SET OF REG. VALUES TO BE USED     (SENT)
C           IRTFLG     ERROR FLAG, ZERO IS NORMAL        (RETURNED)
C
C    PURPOSE:    SUBSTITUTE REGISTERS AND DO LOOP INDICES INTO
C                CERTAIN POSITIONS OF A STRING (NOT NECESSARILY A
C                FILENAME, IT CAN BE ANY STRING). ALSO SUBSTITUTES FOR
C                ENVIRONMENTAL VARIABLES ${...}.
C
C    DESCRIPTION: SUBSTITUTES FOR REGISTER AND LOOP INDEX. 
C                 EXAMPLES IF ([v1]=1, [v11]=11, [v1111]=1111, INDEX:A=44)
C
C                 avg***[v1]           ---> avg001
C                 dir*/avg***[v1][v11] ---> dir1/avg011
C
C                 avg000{---[v1]}      ---> avg001
C                 avg000{----[iv1]}    ---> av0001 (no error reported)
C                 avg000{---[v1111]}   ---> av1111 (no error reported)
C                 avg000{---a}         ---> avg044 
C
C                 avg{***[v1]}         ---> avg001
C                 avg{****[v1]}        ---> avg0001
C                 avg{**a}             ---> avg44 
C                 avg{****}            ---> avg***  (no error reported)
C                 avg{***[v1111]}      ---> ERROR DAMAGES INVARIANT PART 
C                                           OF FILE NAME
C                 EXAMPLES IF ([v2]=-1
C                 avg{***[v2]}         ---> avg0-1
C                 {%i3%[v2]}           --->  -1
C                 {%f3%[v2]}           ---> -1.0000000
C                 {%f4.2%[v2]}         ---> -1.0
C                 {%g14.1%[v2]}        ---> -1.
C
C  NOTE:         CALLS ERRT IF ERROR OCCURS
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE FILNAMSUB(FILNAM,NLET,IBANK,IRTFLG)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      CHARACTER(LEN=*)      :: FILNAM
      INTEGER               :: NLET,IBANK,IRTFLG

      CHARACTER(LEN=1)      :: NULL = CHAR(0)
      CHARACTER(LEN=1)      :: JCHAR,CTEMP
      CHARACTER(LEN=4)      :: NAME
      CHARACTER(LEN=12)     :: SUBPAT
      LOGICAL               :: ISCHAR,CHKMULTIPLE
C     CHARACTER(LEN=80)     :: SUBENV,CVAL,FMT,CSUB JULY 2014 al
      CHARACTER(LEN=MAXNAM) :: SUBENV,CVAL,FMT,CSUB

      DATA SUBPAT/'{***********'/

      IRTFLG      = 1
      CHKMULTIPLE = .TRUE.

      IGOM   = 1    ! THIS ENSURES MULTIPLE TERMINAL [] ONLY SEARCHED 1X

C     SUBSTITUTION LOOP STARTS HERE
10    CONTINUE

C     SEE IF "$(...)" IS PRESENT FOR ENVIRONMENTAL VARIABLES ---------
      CALL CHARINSIDE(FILNAM(1:NLET),'${','}',.FALSE.,.FALSE.,
     &                IGOP,IENDP,NCHARP)

      IF (NCHARP > 0) THEN   
C        ENVIRONMENTAL VARIABLE SUBSTITUTION NEEDED 
         CALL MYGETENV(FILNAM(IGOP+2:IENDP-1),SUBENV,NCHARS,NULL,IER)
         CALL SUBCHAR(SUBENV(1:NCHARS),FILNAM,IGOP,IENDP,NLET,IRTFLG)

C        SEE IF ANY MORE ENV. VARIABLE SUBSTITUTIONS ARE NEEDED
         GOTO 10
      ENDIF

C     HANDLE A FILENAME TERMINATING WITH REG. OR LOOP ---------------
C     I.E.: ABC***[reg]  OR   ABC000[reg]  OR   ABC***i (at very end)

      NVAR = 0
      DO I= 1,NLET
         IF (FILNAM(I:I) == '[') NVAR = NVAR + 1
      ENDDO

      IF (NVAR == 1 .AND. 
     &    NLET > 3  .AND. 
     &    FILNAM(NLET:NLET) == ']' .AND.
     &    FILNAM(1:1) .NE. '?') THEN    ! ? NEEDED FOR 'FR' USE

C        PROBABLE TERMINAL REGISTER (STRING VARIABLES ALREADY REPLACED)

         IGOVAR = SCAN(FILNAM, '[', .TRUE.) 
         IENDA  = IGOVAR - 1                 ! POSSIBLE END OF ***'s

         IF (IENDA > 0) THEN
C           DEFINATE TERMINAL REGISTER: ....[REG]

C           FIND LAST NON-ASTERICK 
            INOTA  = VERIFY(FILNAM(1:IENDA), '*', BACK=.TRUE.)
            IASTGO = INDEX(FILNAM(1:IENDA),'*')

            !write (6,*)'inota,ienda,iastgo:',inota,ienda,iastgo

            IF (INOTA == IENDA .AND. IASTGO == 0) THEN
C              NO *(s) ANYWHERE, USE TERMINAL DIGITS OR OVERWRITE

               INOTDIG =VERIFY(FILNAM(1:IENDA),'0123456789',BACK=.TRUE.)
               !write(6,*)' inota,ienda:inotdig' ,inota,ienda,inotdig ! 4,4,4 
                  
               IF (INOTDIG .NE. IENDA) THEN
C                 HAS TERMINAL DIGITS, REPLACE THEM WITH *(s) 
                  IGOA = INOTDIG + 1
                  NAST = IGOVAR - IGOA
                  !write (6,*)'set0 =',igoa
                  !write(6,*)' csub:',CSUB(1:NCHARSUB),':',FILNAM(1:20)

                  CSUB = SUBPAT(1:NAST+1)//
     &                          FILNAM(IGOVAR:NLET)// '}'//NULL 
                  NCHARSUB = NAST+1+NLET-IGOVAR+1+1

                  CALL SUBCHAR(CSUB(1:NCHARSUB),FILNAM,
     &                         IGOA,NLET,NLET,IRTFLG)
              ENDIF
            ENDIF
        ENDIF

      ELSEIF (NLET > 1 .AND. ISCHAR(FILNAM(NLET:NLET)) ) THEN
C        POSSIBLE TERMINAL LOOP INDEX: ..**I
         IGOVAR = NLET 
         IENDA  = IGOVAR - 1 
       
         IF (IENDA > 0) THEN
C           ASSUME TERMINAL LOOP INDEX: ....[REG]
            INOTA = VERIFY(FILNAM(1:IENDA), '*', BACK=.TRUE.)

            IF (INOTA .NE. IENDA) THEN
C              HAS  *(s), SUBSTITUTE FOR *(s)
               IGOA = INOTA + 1
               NAST = IGOVAR - IGOA 
               CALL SSUPCAS(FILNAM(NLET:NLET))
               CSUB = SUBPAT(1:NAST+1) // '[_' // FILNAM(NLET:NLET) //
     &                       ']}' // NULL 

               NCHARSUB = NAST+1+2+1+2
               CALL SUBCHAR(CSUB(1:NCHARSUB),FILNAM,
     &                      IGOA,NLET,NLET,IRTFLG)
            ENDIF
         ENDIF
      ENDIF
      !write(6,*)' filnam:',FILNAM(1:20),IGOA,NLET


C     HANDLE FILNAM WITH MULTIPLE *** SUBSTITUTIONS  ----Jan07-----
C     E.G.  ABC***/DEF***[reg1][reg2]  
C     THIS WOULD BE MUCH SIMPLER IF FORTRAN HAD REGULAR EXPRESSIONS!

      !write(6,*) ' after filnam:', filnam(1:nlet),' igom:',igom, '  nlet:',nlet

      DO WHILE (CHKMULTIPLE .AND. IGOM <= (NLET - 2))
         CHKMULTIPLE = .FALSE.

         !write(6,*) ' Remaining filnam: ', filnam(igom:nlet)

C        FIND FIRST '*' IN REMAINING PART OF FILENAME
         IGOAST = INDEX(FILNAM(IGOM:NLET),'*')
         IF (IGOAST <= 0) EXIT

C        FOR THIS SUBSTITUTION WE CAN NOT HAVE A '{' JUST BEFORE FIRST '*'
         IF (IGOAST > 1 .AND.
     &       FILNAM(IGOM+IGOAST-2:IGOM+IGOAST-2) == '{')EXIT

C        FIND NUMBER OF '*' IN THIS ***
         NAST = VERIFY(FILNAM(IGOM+IGOAST-1:NLET),'*') - 1

C        NEED SOMETHING AFTER LAST '*' IN THIS ***
         IF (NAST <= 0) EXIT

         IGOVAR = INDEX(FILNAM(IGOM:NLET), '[') + IGOM -1 
         !write(6,*)' FILNAM(IGOM:NLET): ', FILNAM(IGOM:NLET)

C        NEED [.] REG. VARIABLE SOMEWHERE AFTER *** SEQUENCE
         IF (IGOVAR <= 0) EXIT

C        FIND NUMBER OF CHAR. IN THIS [] REG. VARIABLE
         NVAR = INDEX(FILNAM(IGOVAR:NLET),']')

C        NEED VALID [.] REG. VARIABLE AFTER *** SEQUENCE
         IF (NVAR <= 2) EXIT

         !write(6,*)' igom,igoast,nast,igovar,nvar: ', igom,igoast,nast,igovar,nvar
          
C        OK GO AHEAD AND SUBSTITUTE
         LOC1 = IGOM + IGOAST - 1
         LOC2 = LOC1 + NAST   - 1
         LOC3 = IGOVAR 
         LOC4 = LOC3 + NVAR   - 1
         NSUB = 1 + (LOC2- LOC1 + 1) + (LOC4 - LOC3 + 1) + 1 

          
         CSUB = '{'//FILNAM(LOC1:LOC2) // 
     &               FILNAM(LOC3:LOC4) // '}' // NULL
         !write(6,*)' nsub,loc1,loc2,loc3,loc4:',nsub,loc1,loc2,loc3,loc4
         !write(6,*)' csub: ', csub(1:nsub)

C        BLANK THE MATCHED VARIABLE [..] OUT                     
         FILNAM(LOC3:) = FILNAM(LOC4+1:) // NULL

C        SUBSTITUTE *** WITH {***[reg]}
         CALL SUBCHAR(CSUB(1:NSUB),FILNAM, LOC1,LOC2,NLET,IRTFLG)
         IF (IRTFLG .NE. 0) EXIT

C        CONTINUE LOOP TO FIND FURTHER *** STRINGS TO BE SUBSTITUTED
         IGOM        = LOC1 + NSUB
         CHKMULTIPLE = .TRUE.
       
         !write(6,*)' subed: ', filnam(1:nlet)
         !write(6,*)' igom,loc1,nsub: ', igom,loc1,nsub
 
      ENDDO
         

C     SEE IF "{...}" IS PRESENT  -------------------------------------    
      CALL CHARINSIDE(FILNAM(1:NLET),'{','}',.FALSE.,.FALSE.,
     &                IGO,IEND,NCHAR)

      IF (NCHAR > 0) THEN
C        HAVE {...}  NEED TO SUBSTITUTE INTO FILNAM

C        HUNT FOR ASTERICKS,- or %'s INSIDE {}
         NMINUS = 0
         NAST   = 0
         NPER   = 0

         J      = IGO + 1
         DO WHILE (J  <  IEND .AND. FILNAM(J:J) == '*')
            NAST = NAST + 1
            J    = J + 1
         ENDDO

         IF (NAST <= 0) THEN
C           HUNT FOR MINUS SIGNS INSIDE {}
            J    = IGO + 1
            DO WHILE (J  <  IEND .AND. FILNAM(J:J) == '-')
               NMINUS = NMINUS + 1
               J      = J + 1
            ENDDO
         ENDIF

         IF (NAST <= 0 .AND. NMINUS <= 0) THEN
C           SEE IF "(%...%)" IS PRESENT 
            CALL CHARINSIDE(FILNAM(J:NLET),'%','%',.FALSE.,.FALSE.,
     &                IGOPER,IENDPER,NPER)
            J = J + IENDPER
         ENDIF

         IF (NAST <= 0 .AND. NMINUS <= 0 .AND. NPER <= 0) THEN
C           NO ASTERICKS, MINUS, OR,PERCENT SIGNS IS AN ERROR
            WRITE(NOUT,*)'*** NO *,- OR, % IN SUBSTITUTION STRING: ',
     &                   FILNAM
            CALL ERRT(100,'FILNAMSUB',NE)
            RETURN

         ELSEIF ((NAST + 2) == NCHAR .OR. 
     &           (NMINUS > 0 .AND. (NMINUS +2) == NCHAR)) THEN
C           SUBSTITUTION STRING HAS {***..} OR {---} ONLY, REMOVE {}

            IF (IEND < NLET) THEN
               FILNAM = FILNAM(1:IGO-1) // 
     &                  FILNAM(IGO+1:IEND-1) // 
     &                  FILNAM(IEND+1:NLET)  // ' '
            ELSE
               FILNAM = FILNAM(1:IGO-1) // 
     &                  FILNAM(IGO+1:IEND-1) // ' '
            ENDIF

            NLET   = NLET - 2
            IRTFLG = 0
            RETURN
         ENDIF

         JCHAR   = FILNAM(J:J)

         IF (JCHAR == '[') THEN
C           REGISTER SUBSTITUTION WANTED

            CALL CHARINSIDE(FILNAM(J:),'[',']',.FALSE.,.FALSE.,
     &                       IGOBRAK,IENDBRAK,NCHARBRAK)

            CALL REG_GET_VAR(IBANK,FILNAM(J+IGOBRAK-1:J+NCHARBRAK-1),
     &                       .FALSE.,RVALT,IREG,IENDVAR,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            IF (NPER <= 0) THEN
C              MAKE REG. CONTENTS INTO INTEGER FOR FILE NUMBER
               IF (RVALT >= 0 ) THEN
                  IVALT = RVALT + 0.5
               ELSE
                  IVALT = RVALT - 0.5
               ENDIF
            ENDIF
            !write(6,*) ' fil:',filnam(J+IGOBRAK-1:J+NCHARBRAK-1) 
            !write(6,*) ' ivalt:',ibank,ireg,iendvar,irtflg,rvalt,ivalt

         ELSEIF (ISCHAR(JCHAR)) THEN
C           OLD STYLE LOOP COUNTER SUBSTITUTION WANTED
            CALL SSUPCAS(JCHAR)
            NAME = '[_' // JCHAR // ']'

            CALL REG_GET_VAR(0,NAME,.FALSE.,RVALT,IREG,IENDVAR,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
            IF (RVALT >= 0 ) THEN
               IVALT = RVALT + 0.5
            ELSE
               IVALT = RVALT - 0.5
            ENDIF

         ELSE
            CALL ERRT(101,'UNDECIPHERABLE SUBSTITUTION REQUEST',NE)
            RETURN
         ENDIF

C        FIND DIGITS IN IVALT
         NDIGITS = NUMDIG(IVALT,0)
         IF (NAST >  0 .AND. NDIGITS > NAST) THEN
            
            WRITE(NOUT,90) IVALT, FILNAM(1:NLET)
90          FORMAT(' *** SUBSTITUTING: ',I7, ' INTO: ',A,/,
     &             '     DAMAGES INVARIANT PART OF STRING')
            CALL ERRT(100,'FILNAMSUB',NE)
            RETURN
         ENDIF

C        FIND LOCATION FOR SUBSTITUTION WITHIN FILENAME
         NSUB      = NAST
         IF (NMINUS > 0) THEN
C           USING {---x?} SUBSTITUTION
            NSUB = MAX(NMINUS,NDIGITS)
C           EAT AWAY SUPPLIED FILENAME WHEN USING - TYPE SUBSTITUTION
            IGO = IGO - NSUB

C           EAT AWAY ANY PRECEEDING ASTERICKS IN FILENAME ALSO
            DO WHILE (IGO >= 2 .AND. FILNAM(IGO-1:IGO-1) == '*')
               IGO  = IGO - 1
               NSUB = NSUB + 1
            ENDDO
 
            IF (IGO <= 0) THEN
               WRITE(NOUT,*) '*** SUBSTITUTION NUMBER: ',IVALT,
     &                       ' WILL NOT FIT IN STRING'
               CALL ERRT(100,'FILNAMSUB',NE)
               RETURN
            ENDIF

         ELSEIF (NPER > 0) THEN
C           USING {%...%x?} SUBSTITUTION
            CALL CHARINSIDE(FILNAM(IGO:IGO+NCHAR-1),'%','%',
     &                      .TRUE.,.FALSE.,IGOPER,IENDPER,NPER)
            FMT = '(' // FILNAM(IGO+IGOPER-1:IGO+IENDPER-1) // 
     &                 ')' //CHAR(0)

            !write(6,*) 'Fmt:',fmt(1:15),':'

            IF (FMT(2:2) == 'i' .OR. FMT(2:2) == 'I') THEN
C              WRITE REG. CONTENT AS INTEGER
               IF (RVALT >= 0 ) THEN
                  IVALT = RVALT + 0.5
               ELSE
                  IVALT = RVALT - 0.5
               ENDIF
               WRITE(CVAL,FMT) IVALT

            ELSE
C              WRITE REG. CONTENT AS REAL
               WRITE(CVAL,FMT) RVALT
            ENDIF

            CVAL = adjustl(CVAL)
            NSUB = lnblnkn(CVAL) 
         ENDIF

         NLET  = IGO + NSUB + (NLET - IEND) - 1
         NLEN  = LEN(FILNAM)

         IF (NLET > NLEN) THEN
            CALL ERRT(102,'SUBSTITUTION IS TOO LONG',NLET)
            RETURN
         ENDIF

C        PRESERVE END OF FILENAME AFTER SUBSTITUTION AREA    
         FILNAM(IGO+NSUB:) = FILNAM(IEND+1:)

C        WRITE SUBSTITUTION CONTENTS INTO THE  STRING
         IF (NPER > 0) THEN
            FILNAM(IGO:IGO+NSUB-1) = CVAL(1:NSUB)


         ELSE

            CALL INTTOCHAR(IVALT,FILNAM(IGO:IGO+NSUB-1),NNN,NSUB)

            IF (NNN < 0) THEN
               CALL ERRT(102,'SUBSTITUTING NUMBER INTO STRING',
     &                  IVALT)
               RETURN
            ENDIF
         ENDIF

        !write(6,*) ' final: ',nper,':',filnam(1:25)

C        SEE IF ANY MORE SUBSTITUTIONS ARE NEEDED
         GOTO 10
      ENDIF

30    IF (INDEX(FILNAM(1:NLET), '$') > 0) THEN
C       CHECK FOR PRESENCE OF $DATEXT IN FILNAM -----------------
        IDATEXT = INDEX(FILNAM(1:NLET),'$DATEXT')
        IF (IDATEXT <= 0) IDATEXT = INDEX(FILNAM(1:NLET),'$datext')

        DO WHILE (IDATEXT > 0) 
C          SUBSTITUTE CURRENT DATEXC FOR $DATEXT
           FILNAM = FILNAM(1:IDATEXT-1) // DATEXC(1:3) // 
     &              FILNAM(IDATEXT+7:NLET)
           NLET      = NLET - 4 
           IDATEXT = INDEX(FILNAM(1:NLET),'$DATEXT')
           IF (IDATEXT <= 0)IDATEXT = INDEX(FILNAM(1:NLET),'$datext')
        END DO

C       CHECK FOR PRESENCE OF $PRJEXT IN FILNAM -----------------
C       DO NOT CHANGE CASE OF PRJEXT HERE!!!!

        IPRJEXT = INDEX(FILNAM(1:NLET),'$PRJEXT')
        IF (IPRJEXT <= 0) IPRJEXT = INDEX(FILNAM(1:NLET),'$prjext')

        DO WHILE (IPRJEXT > 0) 
C          SUBSTITUTE CURRENT PRJEXC FOR $PRJEXT
           FILNAM = FILNAM(1:IPRJEXT-1) // PRJEXC(1:3) // 
     &              FILNAM(IPRJEXT+7:NLET)
           NLET      = NLET - 4 
           IPRJEXT = INDEX(FILNAM(1:NLET),'$PRJEXT')
           IF (IPRJEXT <= 0) IPRJEXT = INDEX(FILNAM(1:NLET),'$prjext')
        END DO
 
C        SEE IF ANY MORE SUBSTITUTIONS ARE NEEDED
         GOTO 30
      ENDIF
      !    write(6,*) ' Return: ',filnam(1:25)

      IRTFLG = 0

      END
