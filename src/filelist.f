
C++*********************************************************************
C
C    FILELIST.F         LONG FILENAMES         JUL 1999 ARDEAN LEITH
C                       OPENDOC PARAMETERS     DEC 2000 ARDEAN LEITH
C                       DOC FILE SLICING       APR 2001 ARDEAN LEITH
C                       INCORE OPENDOC         JUL 2003 ARDEAN LEITH
C                       KEYED ILIST            SEP 2003 ARDEAN LEITH
C                       XMIPP SELFILE          DEC 2010 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C
C FILELIST(GETTEMPLATE,LUNDOC,FILPAT,NLETP,ILIST,NMAX,NUM,PROMPT,IRTFLG)
C
C     PURPOSE:      INPUTS FILE NAME TEMPLATE AND NUMBERS FOR FILE
C                   NAME LOOP.  USUALLY USED WITH FILGET.FOR
C                   I.E.  CALL FILGET(FILPAT,FILNAM,NLET,INUM,IRTFLG)
C                         
C     PARAMETERS:   GETTEMPLATE   FLAG TO INPUT TEMPLATE      (SENT)
C                   LUNDOC        DOC FILE I/O UNIT           (SENT)
C                   FILNAM        FILE NAME PATTERN           (RETURNED)
C                   NLETP         LENGTH OF FILNAM            (RETURNED)
C                   ILIST         ARRAY FOR NUMBERS           (RETURNED)
C                   NMAX          MAX. LENGTH OF ILIST        (SENT)
C                                  IF ZERO ONLY GETS FILPAT NOT ILIST
C                                  IF < ZERO GETS KEYED ILIST
C                   NUM           NUMBER OF VALUES IN ILIST   (RETURNED)
C                                   IF < ZERO, SELFILE ON LUNDOC
C                   PROMPT        PROMPT                      (SENT)
C                   IRTFLG        ERROR FLAG; 0 IS NORMAL     (RETURNED)
C
C--*********************************************************************

        SUBROUTINE FILELIST(GETTEMPLATE,LUNDOC,FILPAT,NLETP,
     &                      ILIST,NMAX,NUM,PROMPT,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'

        LOGICAL           :: GETTEMPLATE
        INTEGER           :: LUNDOC
        CHARACTER(LEN=*)  :: FILPAT
        INTEGER           :: NLETP

C       ILIST IS DIMENSIONED AS (*) HERE SO NMAX=0 IS ACCEPTED
C**	INTEGER*4     ILIST(NMAX)      ! ACTUAL SIZE
        INTEGER           :: ILIST(*)
        INTEGER           :: NMAX,NUM
        CHARACTER(LEN=*)  :: PROMPT
        INTEGER           :: IRTFLG

        CHARACTER(LEN=1)  :: NULL
        LOGICAL           :: GOTAST
        INTEGER           :: NSEL

        NULL = CHAR(0)

        IF (GETTEMPLATE) THEN
C          GET FILE NAME TEMPLATE 
           CALL FILELISTA(FILPAT,NLETP,PROMPT,LUNDOC,NSEL,IRTFLG)
           IF (IRTFLG .NE. 0 .AND. NSEL > 0) THEN
C             XMIPP SELFILE OPENED ON LUNDOC
              NUM = -NSEL     ! # OF SELECTED IMAGES IN LIST
              RETURN
           ENDIF

C          SEE IF INPUT CONTAINS A TEMPLATE
           GOTAST = (INDEX(FILPAT,'*') > 0) 
        ELSE
           GOTAST = .TRUE.
        ENDIF

        IF (GOTAST) THEN
           IF (NMAX .LT. 0) THEN
C             FILL THE NUMBERS ARRAY ALSO
              CALL FILELISTC(LUNDOC,ILIST,-NMAX,NUM,NULL,IRTFLG)

           ELSEIF (NMAX .GT. 0) THEN
C             FILL THE NUMBERS ARRAY ALSO
              CALL FILELISTB(LUNDOC,ILIST,NMAX,NUM,NULL,IRTFLG)
           ENDIF
        ENDIF

        RETURN
        END

C       ********************* FILELISTA *******************************
      
        SUBROUTINE FILELISTA(FILPAT,NLETP,PROMPT,LUNXM,NSEL,IRTFLG)

        IMPLICIT NONE
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*)  :: FILPAT   ! NAME PATTERN           (RET)
        INTEGER           :: NLETP    ! CHAR IN PATTERN        (RET)
        CHARACTER(LEN=*)  :: PROMPT   ! OPTIONAL PROMPT        (SENT)
        INTEGER           :: LUNXM    ! LUN FOR SELFILE        (SENT)
        INTEGER           :: NSEL     ! # FILES IN SELFILE     (RET)
        INTEGER           :: IRTFLG   ! ERROR FLAG             (RET)

        CHARACTER(LEN=81)     :: PROMPT2
        CHARACTER(LEN=1)      :: NULL
        INTEGER               :: lnblnkn    
        INTEGER               :: NLET,LOCAST,LOCAT
        CHARACTER(LEN=MAXNAM) :: FILNAM  

        NULL    = CHAR(0)

        PROMPT2 = PROMPT
        IF (PROMPT == NULL) THEN
           PROMPT2 = 'TEMPLATE FOR FILENAMES (E.G. PIC@****)'
        ENDIF

C       DO NOT CHANGE CASE OF THE RDPRMC INPUT
        IRTFLG = -999

C       READ IN FILE NAME TEMPLATE
        NLET = lnblnkn(PROMPT2)
        CALL RDPRMC(FILPAT,NLETP,.TRUE.,PROMPT2(:NLET),NULL,IRTFLG)
        IF (IRTFLG == -1) RETURN

        LOCAST = INDEX(FILPAT(1:NLETP),'*')    
        LOCAT  = INDEX(FILPAT(1:NLETP),'@')    

        IF (LOCAST == 0 .AND. LOCAT == 0) THEN
C          CHECK FOR XMIPP SELFILE
           CALL OPENXMSELFILE(FILPAT(:NLETP),LUNXM,
     &                        FILNAM,NLET,NSEL,IRTFLG)

        ELSEIF (NLETP == 3 .AND. 
     &          FILPAT(NLETP:NLETP) .NE. '*' .AND.
     &          FILPAT(1:1) .NE. '_') THEN
C          MAKE NEW STYLE TEMPLATE
           WRITE(NOUT,*) ' *** POSSIBLE: ARCHAIC 3 CHAR. TEMPLATE? ',
     &                    ' NO LONGER SUPPORTED.'
           WRITE(NOUT,*) ' IF THIS IS 3 CHAR. TEMPLATE, APPEND ***'
        ENDIF

        FILPAT(NLETP+1:NLETP+1) = NULL

        RETURN
        END

C       ********************* FILELISTB *******************************

        SUBROUTINE FILELISTB(LUNDOCT,ILIST,NMAX,NUM,PROMPT,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*) ::      PROMPT
        CHARACTER(LEN=MAXNAM) :: FILLIST,PROMPT2
        CHARACTER(LEN=1) ::      NULL
        LOGICAL      ::          ISCHAR,ISFILENAME,NEWFILE
        REAL, DIMENSION(2) ::    PLIST

C       ILIST IS DIMENSIONED AS (1) HERE SO NMAX=0 IS ACCEPTED
C**	INTEGER*4     ILIST(NMAX)      ! ACTUAL SIZE
        INTEGER, DIMENSION(*) :: ILIST

        NULL = CHAR(0)
        
        PROMPT2 = PROMPT
        IF (PROMPT(1:1) == NULL) THEN
           PROMPT2 = 'FILE NUMBERS OR SELECTION DOC. FILE NAME' 
        ENDIF
        LENPROM = LNBLNKN(PROMPT2)

C       FILL THE NUMBERS ARRAY 

C       GET SELECTION FILENAME OR FILE NUMBER LIST
C       IRTFLG OF -999 SAYS DO NOT UPPERCASE RDPRMC INPUT
        IRTFLG = -999
        CALL RDPRMC(FILLIST,NLET,.TRUE.,PROMPT2(:LENPROM),NULL,IRTFLG)

        IF (ISFILENAME(FILLIST,NLET)) THEN
C          FILLIST IS A SELECTION DOC FILE NAME

C          CHECK FOR SLICING (X?? X?? SEPARATED FROM NAME)
           LOCB = INDEX(FILLIST(1:NLET),' ')
           
           IF (LOCB .LE. 0) THEN
C             FILL THE NUMBERS ARRAY (ILIST) FROM SELECTION FILE 
              FILLIST(NLET+1:) = NULL
              CALL OPENDOC(FILLIST,.TRUE.,NLET,LUNDOCT,LUNDOC,.FALSE.,
     &                     ' ', .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL LUNDOCREDSEL(LUNDOC,ILIST,NMAX,NUM,MAXGOTY,IRTFLG)
             
           ELSE
C             PROBABLY WANT A SLICE FROM DOC FILE
C             TILDE SIGN AT BEGINNING OF PROMPT OVERRIDES INPUT IN RDPRAI
              PROMPT2 = '~' // FILLIST(LOCB:)

C             GET LIST OF NUMBERS CONTAINED IN SLICE (MAY BE X??) 
C             MAXIMUM VALUE PLACED IN ILIST IS 9999999 CURRENTLY
              NUM = 2
              CALL RDPRAI(ILIST,NMAX,NUM,0,9999999,PROMPT2(1:NLET+1),
     &                    NULL,IRTFLG)

C             OPEN DOC FILE
              NLET   = LOCB - 1
              CALL OPENDOC(FILLIST(1:NLET),.TRUE.,NLET,LUNDOCT,LUNDOC,
     &              .FALSE.,' ', .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

C             RETRIEVE REGISTER 1 VALUES FROM SPECIFIED KEYS 
              IGOY  = ILIST(1)
              IENDY = ILIST(2)
CCC           CALL LUNDOCREDSLI(LUNDOC,ILIST,NMAX,.TRUE.,1,IGO,IEND,NUM,IRTFLG)
              CALL LUNDOCREDSLC(LUNDOC,.TRUE.,ILIST,DUM,1,NMAX,
     &              .FALSE.,.FALSE.,1,1, IGOY,IENDY, NUM,MAXGOTY,IRTFLG)
           ENDIF
           CLOSE(LUNDOCT)             
        ELSE

C          TILDE SIGN AT BEGINNING OF PROMPT OVERRIDES INPUT IN RDPRAI
           PROMPT2 = '~' // FILLIST(1:MAXNAM-1)

C          SET NUM TO NMAX FOR NUMBER OF FILES ALLOWED
           NUM = NMAX

C          MAXIMUM VALUE PLACED IN ILIST IS 9999999 CURRENTLY
           CALL RDPRAI(ILIST,NMAX,NUM,0,9999999,PROMPT2(1:NLET+1),
     &                 NULL,IRTFLG)
        ENDIF

        RETURN
        END

C       ********************* FILELISTC *******************************

        SUBROUTINE FILELISTC(LUNDOCT,ILIST,NMAX,NUM,PROMPT,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=*) ::      PROMPT
        CHARACTER(LEN=MAXNAM) :: FILLIST,PROMPT2
        CHARACTER(LEN=1) ::      NULL
        LOGICAL      ::          ISCHAR,ISFILENAME,NEWFILE
        REAL, DIMENSION(2) ::    PLIST

C       ILIST IS DIMENSIONED AS (1) HERE SO NMAX=0 IS ACCEPTED
C**	INTEGER*4     ILIST(NMAX)      ! ACTUAL SIZE
        INTEGER, DIMENSION(*) :: ILIST
        INTEGER, ALLOCATABLE,DIMENSION(:) :: ILISTT

        NULL = CHAR(0)
        
        PROMPT2 = PROMPT
        IF (PROMPT(1:1) == NULL) THEN
           PROMPT2 = 'FILE NUMBERS OR SELECTION DOC. FILE NAME' 
        ENDIF
        LENPROM = LNBLNKN(PROMPT2)

C       FILL THE NUMBERS ARRAY 

C       GET SELECTION FILENAME OR FILE NUMBER LIST
C       IRTFLG OF -999 SAYS DO NOT UPPERCASE RDPRMC INPUT
        IRTFLG = -999
        CALL RDPRMC(FILLIST,NLET,.TRUE.,PROMPT2(:LENPROM),NULL,IRTFLG)

        IF (ISFILENAME(FILLIST,NLET)) THEN
C          FILLIST IS A SELECTION DOC FILE NAME

C          CHECK FOR SLICING (X?? X?? SEPARATED FROM NAME)
           LOCB = INDEX(FILLIST(1:NLET),' ')
           
           IF (LOCB .LE. 0) THEN
C             FILL THE NUMBERS ARRAY (ILIST) FROM SELECTION FILE 
              FILLIST(NLET+1:) = NULL
              CALL OPENDOC(FILLIST,.TRUE.,NLET,LUNDOCT,LUNDOC,.FALSE.,
     &                     ' ', .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

              CALL LUNDOCREDSLC(LUNDOC,.TRUE.,ILIST,DUM,1,NMAX,
     &             .TRUE.,.FALSE. ,1,1, 1,MAXY,NUM,MAXGOTY,IRTFLG)

           ELSE
C             PROBABLY WANT A SLICE FROM DOC FILE
C             TILDE SIGN AT BEGINNING OF PROMPT OVERRIDES INPUT IN RDPRAI
              PROMPT2 = '~' // FILLIST(LOCB:)

C             GET LIST OF NUMBERS CONTAINED IN SLICE (MAY BE X??) 
C             MAXIMUM VALUE PLACED IN ILIST IS 9999999 CURRENTLY
              NUM = 2
              CALL RDPRAI(ILIST,NMAX,NUM,0,9999999,PROMPT2(1:NLET+1),
     &                    NULL,IRTFLG)

C             OPEN DOC FILE
              NLET   = LOCB - 1
              CALL OPENDOC(FILLIST(1:NLET),.TRUE.,NLET,LUNDOCT,LUNDOC,
     &              .FALSE.,' ', .TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
              IF (IRTFLG .NE. 0) RETURN

C             RETRIEVE REGISTER 1 VALUES FROM SPECIFIED KEYS 
              IGOY  = ILIST(1)
              IENDY = ILIST(2)
CCC           CALL LUNDOCREDSLI(LUNDOC,ILIST,NMAX,.TRUE.,1, CCC                                          IGO,IEND,NUM,IRTFLG)
              CALL LUNDOCREDSLC(LUNDOC,.TRUE.,ILIST,DUM,1,NMAX,
     &              .TRUE.,.FALSE.,1,1, IGOY,IENDY, NUM,MAXGOTY,IRTFLG)
           ENDIF
           CLOSE(LUNDOCT)             
        ELSE

C          TILDE SIGN AT BEGINNING OF PROMPT OVERRIDES INPUT IN RDPRAI
           PROMPT2 = '~' // FILLIST(1:MAXNAM-1)

C          SET NUM TO NMAX FOR NUMBER OF FILES ALLOWED
           NUM = NMAX

           ALLOCATE (ILISTT(NUM),STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN
              CALL ERRT(46,'FILELISTC, ILISTT....',NUM)
              RETURN
           ENDIF

C          MAXIMUM VALUE PLACED IN ILIST IS 9999999 CURRENTLY
           CALL RDPRAI(ILISTT,NMAX,NUM,0,9999999,PROMPT2(1:NLET+1),
     &                 NULL,IRTFLG)

           ILIST(1:NMAX) = 0
           DO I = 1,NUM
              IT = ILISTT(I)
              IF (IT .LT. NMAX) THEN
                 ILIST(IT) = NMAX
              ENDIF
           ENDDO  
        ENDIF

        RETURN
        END


C     ********************* ISFILENAME  *******************************

      LOGICAL FUNCTION ISFILENAME(STRING,NLET)

C     CRITERION IS THAT A FILENAME MUST ALWAYS HAVE A ALPHABETIC
C     CHARACTER WHICH IS NOT "X" OR "x", AND WHICH IS NOT 
C     WITHIN {}'s OR []'s.

      CHARACTER(LEN=*) :: STRING  
      LOGICAL          :: INBRAK,ISCHAR,INSQBRAK
      CHARACTER(LEN=1) :: CTEMP

      ISFILENAME = .FALSE.
      INBRAK     = .FALSE.
      INSQBRAK   = .FALSE.

      DO I=1,NLET
         CTEMP = STRING(I:I)

         IF (ISCHAR(CTEMP)) THEN
C           CHAR. (A..Z) SAYS STRING MAY BE A FILE NAME ?
            IF (.NOT. INBRAK .AND. .NOT. INSQBRAK) THEN
               IF (CTEMP .NE. 'X' .AND. CTEMP .NE. 'x') THEN
                  ISFILENAME = .TRUE.
                  RETURN
               ENDIF
            ENDIF
         ELSE
            IF (CTEMP == '{') INBRAK   = .TRUE.
            IF (CTEMP == '[') INSQBRAK = .TRUE.

            IF (CTEMP == '}') INBRAK   = .FALSE.
            IF (CTEMP == ']') INSQBRAK = .FALSE.
         ENDIF
      ENDDO

      RETURN
      END



