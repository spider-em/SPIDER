
C++*********************************************************************
C
C OPFILEC.F                AUTHOR: ArDean Leith                        
C                          REMOVED OPENALL CALL    JAN 99 ArDean Leith
C                          USED LUNHDR             FEB 99 ArDean Leith
C                          CAN KEEP EXTENSION      NOV 02 ArDean Leith
C                          INDEXED STACKS          JAN 03 ArDean Leith
C                          HEADER COPY             FEB 03 ArDean Leith
C                          OPFIL --> OPFILEC       FEB 03 ArDean Leith  
C                          REMOVED IRTFLG INPUT    APR 04 ArDean Leith
C                          DISP 'B'                APR 04 ArDean Leith
C                          FOURIER ERRT            DEC 10 ArDean Leith
C                          KEEPEXT LOGIC           JUN 14 ArDean Leith
C                          INDEX(FILNAM,'.' BUG    AUG 14 ArDean Leith
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
C  OPFILEC(LUNT,ASKNAM,FILNAM,LUN,DISPT,NX,NY,ITYPE,NZ,
C          MAXIM,PROMPT,FOUROK,IRTFLG)
C
C  PURPOSE:         SOLICITS FILE NAME AND OPENS FILE
C  
C  PARAMETERS:
C 
C        LUNT       UNIT TO COPY HEADER VALUES FROM               (SENT)
C        ASKNAM     LOGICAL FLAG TO QUERY NAME                    (SENT) 
C        FILNAM     FILENAME (WITOUT EXTENSION)               (SENT/RET)
C        LUN        UNIT TO OPEN FILE ON                          (SENT)
C        DISPT      CHARACTER CONTAINING ONE OF THE               (SENT) 
C                   FOLLOWING DISPOSITION SPECIFICATIONS:
C                   'O'   -  FILE IS ASSUMED TO EXIST.  DIMENSIONS,
C                            ITYPE AND HEADER INFO (IN COMMON) ARE 
C                            RETURNED TO THE CALLING PROGRAM. 
C                   'B'   -  SAME AS OLD BUT NO LIMIT ON BUFFER LENGTH
C                            FOR OPENCHK. 
C                   'Z/E' -  THE FILE IS ASSUMED TO EXIST.
C                            IF FILE DOES NOT EXIST, THEN BATCH DOES
C                            NOT STOP. (ONLY DIFFERENCE FROM 'O'). 
C                   'N'  -   WANT NEW FILE. NX, NY, NZ, AND
C                            ITYPE MUST BE SENT.
C                   'U'  -   IT IS NOT KNOWN IF THE FILE EXISTS.  
C                            NX, NY, NZ, AND ITYPE MUST 
C                            BE SENT.  IF THE FILE ALREADY EXISTS, IT 
C                            WILL BE REPLACED.
C                   'K & M'-    NO LONGER USED.
C        ITYPE      IFORM FOR FILE                         (SENT OR RET) 
C        NX         IMAGE SIZE                             (SENT OR RET)
C        NY         IMAGE SIZE                             (SENT OR RET)
C        NZ         IMAGE SIZE                             (SENT OR RET)
C        MAXIM      STACK INDICATOR                           (SENT/RET)
C                   ON INPUT (IF NEW):
C                        0 IS FOR SPECIFIC IMAGE           
C                       +1 STACK                               
C                       -1 INDEXED STACK                               
C                   ON INPUT (IF EXISTING):
C                        0 IS FOR SPECIFIC IMAGE ONLY         
C                       <0 or >0 ALLOWS WHOLE STACK OPERATION                              
C                   ON OUTPUT:
C                       -2 NON-STACK IMAGE                
C                       -1 STACKED IMAGE                  
C                     >= 0 IS CURRENT MAX. IMAGE NO. FOR STACK             
C        PROMPT     PROMPT FOR FILNAME                            (SENT)
C                       AT END:  ~  SKIPS FILE ON PROMPT
C                                ~7 CAN OPEN A STACK WITHOUT @
C                                ~9 KEEPS INCOMING EXTENSION
C                                ~6 KEEPS OLD DATE/TIME
C        FOUROK     CAN USE EXISTING FOURIER FILES?               (SENT)
C        IRTFLG     ERROR FLAG (0 IS NORMAL)                      (RET)
C                        -1 GOTO PREVIOUS QUESTION
C                         0 NORMAL RETURN
C                         2 CAN'T USE AN EXISTING FOURIER FILE
C                         3 NO @ ON A BARE STACK FILE
C                         4 OPERATION DOESN'T WORK ON WHOLE STACKS
C
C    CODING:   BASED ON PARAMETERS NX,NY, & NZ A
C              NEW FILE IS OPENED WITH IREC RECORDS, EACH NX*4 
C              BYTES LONG.  IREC ALLOWS SPACE FOR THE 2-D OR 3-D 
C              IMAGE  PLUS HEADER.  A STACK FILE CONTAINS AN OVERALL
C              HEADER PLUS MAXIM * IREC RECORDS. EACH IMAGE IN THE
C              STACK HAS ITS OWN HEADER RECORD(S) WHOSE FORMAT IS THE
C              SAME AS THE OVERALL HEADER RECORDS
C
C    COMMON VARIABLES:
C              IFORM (TYPE)  FILE TYPE SPECIFIER.                ( RET.)
C               +1    R     2-D IMAGE FILE
C               +3    R3    3-D IMAGE (VOLUME) FILE
C               -9    FS    3-D SIMPLE FORMAT FOURIER (MR'S FORMAT)
C               -11   FO    2-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -12   FE    2-D FOURIER TRANSFORM, MIXED RADIX EVEN
C               -21   FE    3-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -22   FE    3-D FOURIER TRANSFORM, MIXED RADIX EVEN
C
C     CALL TREE:
C                                                        
C                                          (regular)
C                           (not stack)     (file)
C     OPFILE  ---------?------------------> OPENFIL -------->
C                      |                      |
C                      |                      `----------> INLNBUF
C                      |                        (inline)   OPENINLN
C                      |                        (file  )
C                      |
C                      |
C                      |    (inline stack)
C                      |-------------------> OPNINSTK -->      
C                      |
C                      |    (file stack)
C                      ` ------------------>  OPENSTK --> OPENFIL -->
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE OPFILEC(LUNT,ASKNAM,FILNAM,LUN,DISPT,ITYPE,
     &                  NX,NY,NZ,MAXIM,PROMPT,FOUROK,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=*) :: FILNAM,PROMPT,DISPT
        CHARACTER (LEN=1) :: DSP,DISP
        LOGICAL           :: ASKNAM,FOUROK,OPSTKNOAT,KEEPEXT
        CHARACTER (LEN=1) :: NULL = CHAR(0)

C       HACK TO OPEN NEW FILES GREATER THAN NBUFSIZ (DANGEROUS)
        DISP = DISPT
        IF (DISP(1:1) == 'B') DISP = 'U'

        LENP = lnblnkn(PROMPT)
C       CAN PASS ~7  AT END OF PROMTP TO OPEN* TO OPEN A STACK WITHOUT @
        OPSTKNOAT = (PROMPT(LENP-1:LENP) == '~7')

C       CAN PASS ~9  AT END OF PROMTP TO OPEN* TO KEEP INCOMING EXTENSION
        KEEPEXT = (PROMPT(LENP-1:LENP) == '~9')

        KEEPEXT = (KEEPEXT      .AND. 
     &             .NOT. ASKNAM .AND.
     &             (INDEX(FILNAM,'.',BACK = .TRUE.)) > 
     &             (INDEX(FILNAM,'/',BACK = .TRUE.)))

C       SET DEFAULT ERROR RETURN IRTFLG
        IRTFLG = 1

        IF (LUN <= 0 .OR. LUN > 100) THEN
           CALL ERRT(102,'PGM. ERROR: LUN MUST BE 1...100',LUN)
           RETURN
        ENDIF

        IF (ASKNAM) THEN
C          SOLICIT FILE NAME, CAN KEEP EXTENSION IF PROMPT ENDS WITH: ~9 
           CALL FILERD(FILNAM,NLETI,NULL,PROMPT,IRTFLG)
           IF (IRTFLG == -1) RETURN
        ENDIF

        DSP = 'O'
        IF (DISP(1:1) .NE. 'O' .AND. 
     &      DISP(1:1) .NE. 'Z' .AND.
     &      DISP(1:1) .NE. 'E') THEN
C          WILL OPEN A NEW FILE, CHECK THAT NECESSARY INFO IS HERE

C          HACK TO OPEN FILES GREATER THAN NBUFSIZ (DANGEROUS)
           NBUFSIZT = NBUFSIZ
           IF (DISPT(1:1) == 'B') NBUFSIZT = HUGE(NBUFSIZT)

           !write(6,*)'   Calling opench with nbufsizt: ',NBUFSIZT

           CALL OPENCHK(NX,NY,NZ,ITYPE,NBUFSIZT,IRTFLGT)
           IF (IRTFLGT .NE. 0) RETURN    
           DSP = 'N'
        ENDIF

C       CREATE A NEW HEADER OBJECT FOR THIS LUN
        CALL LUNNEWHDR(LUN,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       PUT FILENAME AND DSP IN OFF-FILE AREA OF THE HEADER OBJECT
        CALL LUNSETFILE(LUN,FILNAM,DSP,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       PUT ISBARE IN STATIC AREA OF HEADER OBJECT
        CALL LUNBAREFILE(LUN,FILNAM,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       MAKE SURE THE STACK OFFSET IS ALWAYS ZEROED, LUNARB SET ...
        CALL LUNSETLUNS(LUN,0,0,LUN,0,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

        IRTFLG = 0
        IF (INDEX(FILNAM,'@') == 0) THEN

C          NOT AN IMAGE STACK (MAY BE AN INLINE FILE).  NOTE:
C          THIS IS THE PATH FOR 'REGULAR' SPIDER IMAGE FILES.

           NSTACK = 0
	   CALL OPENFIL(LUNT,FILNAM,LUN,NX,NY,NZ,NSTACK,
     &                  ITYPE,DISP(1:1),KEEPEXT,IRTFLG)

           !write(6,*) ' opfilec, openfil returned irtflg:',irtflg,filnam

           IF (IRTFLG .NE. 0) RETURN

C          RETURNS NSTACK -2 FOR NON-STACK or MAXIM

           IF (NSTACK >= 0 .AND. .NOT. OPSTKNOAT ) THEN 
C             BARE STACK REFERENCE ALLOWED WITHOUT '@' FROM 'ST' & 'LI'
C             BARE STACK REFERENCE NOT ALLOWED WITHOUT '@' NORMALLY
              CALL ERRT(101,'STACK INDICATOR (@) MISSING',NE)
              IRTFLG = 3
              RETURN
           ENDIF

        ELSEIF (FILNAM(1:1) == '_') THEN
C           INLINE IMAGE STACK ACCESS WANTED
            NSTACK = MAXIM
            CALL OPENINSTK(LUNT,FILNAM,LUN,NX,NY,NZ,
     &                     NSTACK,ITYPE,DISP(1:1),IRTFLG)

        ELSE
C          WANT TO ACCESS A FILE BASED IMAGE STACK
           NSTACK = MAXIM
           CALL OPENSTK(LUNT,FILNAM,LUN,NX,NY,NZ,
     &                  NSTACK,ITYPE,DISP(1:1),IRTFLG)

        ENDIF

C----------------------------------------

C       RETURN IF THERE IS ANY ERROR OPENING FILE IN OPEN*
        IF (IRTFLG .NE. 0) RETURN

        IF ((DISP(1:1) == 'O'  .OR. 
     &       DISP(1:1) == 'Z'  .OR. 
     &       DISP(1:1) == 'E') .AND.
     &          ITYPE < 0 .AND. .NOT. FOUROK) THEN
C          CAN NOT USE EXISTING FOURIER FILE
           CALL ERRT(101,'OPERATION DOES NOT ACCEPT FOURIER FILES',NE)
           IRTFLG = 2
           RETURN

        ELSEIF (.NOT. OPSTKNOAT .AND. 
     &          MAXIM == 0 .AND. NSTACK >= 0) THEN
C          THIS OPERATION DOES NOT ACCEPT WHOLE STACKS 
           CALL ERRT(101,'OPERATION DOES NOT WORK ON WHOLE STACKS',NE)
           IRTFLG = 4
           RETURN
        ENDIF

        IFORM = ITYPE

C       RETURN NSTACK IN MAXIM 
        MAXIM = NSTACK

        END
