
C++*********************************************************************
C                                                                      
C  STAR2DOC.F  XMIPP STAR FILE SUPPORT           APR 2013 ARDEAN LEITH 
C              LINE END BUG                      FEB 2014 ARDEAN LEITH 
C              REWRITE WITH VAR NAMES            FEB 2014 ARDEAN LEITH 
C              NREGSTAR = MIN(NVARSTARNAME BUG   MAY 2014 ARDEAN LEITH
C              MAXLENVAR = 40 --> 120            OCT 2015 ARDEAN LEITH
C              NVARSTAR BUG                      OCT 2015 ARDEAN LEITH
C              INTEGER OVERFOW HACK              NOV 2015 ARDEAN LEITH
C              _ NOT _rln                        SEP 2016 ARDEAN LEITH
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
C     STAR2DOC(LUNXM,LUNDOC)
C
C     PURPOSE:      CRUDE CONVERTOR FROM STAR FILE TO DOC FILE
C                         
C     PARAMETERS:   LUNXM,LUNDOC I/O UNITS                      (SENT)
C
C     VARIABLES:     
C        
C     NOTE:         NEEDS MPI MODS IF CALLED UNDER MPI
C
C--*********************************************************************

      SUBROUTINE STAR2DOC(LUNXM,LUNDOC)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER                  :: LUNXM   
      INTEGER                  :: LUNDOC  

      CHARACTER(LEN=MAXNAM)    :: FILXM,DOCNAM 

      INTEGER, PARAMETER       :: MAXLENVAR = 120
      INTEGER, PARAMETER       :: NVARSMAX  = 50

      INTEGER                  :: NVARHEAD,NVARSTAR,IHEAD,ISTAR
      INTEGER                  :: NVARSPI, NVARSTARNAME,NREGSTAR 
      INTEGER                  :: NVARSPINAME,IGOT,NSPI

      CHARACTER(LEN=MAXLENVAR) :: VARHEAD     (NVARSMAX)
      CHARACTER(LEN=MAXLENVAR) :: VARSTAR     (NVARSMAX)
      CHARACTER(LEN=MAXLENVAR) :: VARSTARNAME (NVARSMAX)
      INTEGER                  :: NCOLSTAR    (NVARSMAX)
      INTEGER                  :: NCOLSPI     (NVARSMAX)
      INTEGER                  :: NCOLSPINAME (NVARSMAX)
      INTEGER                  :: NCOLSTARNAME(NVARSMAX)
      CHARACTER(LEN=MAXLENVAR) :: VARTMP      (NVARSMAX)
      REAL                     :: DLIST       (NVARSMAX)
    ! DOUBLE PRECISION         :: DLIST       (NVARSMAX)
      CHARACTER(LEN=MAXLENVAR) :: CTOKEN     
      CHARACTER(LEN=240)       :: RECLIN 
      CHARACTER(LEN=240)       :: RECNAM 
      CHARACTER(LEN=10)        :: VALID = '0123456789'
      CHARACTER(LEN=100)       :: FORMOUT

      INTEGER                  :: ILABEL(NVARSMAX) 
      integer                  :: idlist,NBAD

      INTEGER                  :: IENDVAR,IGO,IEND,NCHAR2        
      INTEGER                  :: IFIRST,NCHAR,NGOT,IPUT,IWHICH       
      INTEGER                  :: NLET,ICOMM,MYPID,MPIERR
       
      CHARACTER(LEN=1)         :: NULL = CHAR(0)  
      CHARACTER(LEN=1)         :: CDUM  
      INTEGER                  :: LUNDOCT,ILOC,IT,NLIST,ILAST,NDUM
      INTEGER                  :: IHI,I,IKEY,IRTFLG,LENREC,IVAR,ITOK

      LOGICAL                  :: ISOPEN,EX,CALLERRT,ASKNAM
      LOGICAL                  :: ADDEXT,GETNAME
      LOGICAL                  :: ISOLD,APPEND,MESSAGE,NEWFILE

      INTEGER                  :: lnblnkn
      LOGICAL                  :: ischar,isdigi

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C     OPEN STAR FILE
      LENREC   = 0              ! FORMATTED, SEQ. ACCESS
      ASKNAM   = .TRUE.
      CALLERRT = .TRUE.
      CALL OPAUXFILE(ASKNAM,FILXM,NULL,LUNXM,LENREC,
     &               'O', 'STAR',CALLERRT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IKEY     = 0
      NVARHEAD = 0
      WRITE(NOUT,*) ' STAR HEADER VARIABLES -------------'

      DO 
C        READ SINGLE HEADER LINE FROM STAR FILE
         READ(LUNXM,'(A)',IOSTAT=IRTFLG) RECLIN
         IF (IRTFLG .NE. 0) RETURN

         NLET = lnblnkn(RECLIN)
         !write(6,*) '  reclin: ', reclin(:nlet)
         IF (NLET <= 0) CYCLE         ! BLANK LINE

C        PARSE STAR FILE HEADER LINE INTO TOKENS
         CALL GET_TOKENS(RECLIN,NLET,MAXLENVAR,2,.FALSE.,' ',
     &                   VARTMP, NGOT, IRTFLG)
         IF (IRTFLG .NE. 0 .OR. NGOT < 1) EXIT

         IF (VARTMP(1)(1:1) == '_') THEN
C           GOT A STAR REGISTER VARIABLE NAME
            NVARHEAD          = NVARHEAD + 1
            NCHAR             = lnblnkn(VARTMP(1))
            NCHAR2            = lnblnkn(VARTMP(2))
            VARHEAD(NVARHEAD) = VARTMP(1)(1:NCHAR)
            WRITE(NOUT,90) '  ',VARHEAD(NVARHEAD)(1:NCHAR),
     &                          VARTMP(2)(1:NCHAR2)
90          FORMAT(A,A,1X,A,I0)

         ELSE
            IF ( INDEX('dl_#',VARTMP(1)(1:1)) == 0) EXIT

         ENDIF              
      ENDDO
      WRITE(NOUT,*) ' '

      IF (NVARHEAD <= 0) THEN
         CALL ERRT(101, 
     &      'NO STAR VARIABLES FOUND IN FILE (I.E. NO _*** VARIABLES)',
     &      NDUM)
         RETURN
      ENDIF



C     GET LIST OF STAR VARIABLES WANTED --------------------------
      IRTFLG = -999   ! DO NOT UPPERCASE RECLIN
      CALL RDPRMC(RECLIN,NLET,.TRUE., 
     &        'STAR FILE REGISTER VARIABLES WANTED', CDUM,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      !write(6,*) ' reclin:',reclin(1:nlet)

C     PARSE STAR VARIABLES LINE INTO TOKEN ARRAY
      CALL GET_TOKENS(RECLIN,NLET,MAXLENVAR,NVARSMAX,.FALSE.,' ,',
     &                VARSTAR, NVARSTAR, IRTFLG)
      IF (IRTFLG .NE. 0 .OR. NGOT  < 1) GOTO 9999

      !write(6,*) ' ngot:',ngot,varstar(1:ngot)(1:10)

C     ASSOCIATE STAR VARIABLES WITH COLUMNS IN STAR FILE ---------
      NCOLSTAR = 0     ! ARRAY ZERO
      DO ISTAR = 1,NVARSTAR
         DO IHEAD = 1,NVARHEAD
            IF (VARSTAR(ISTAR) == VARHEAD(IHEAD)) THEN
               NCOLSTAR(ISTAR) = IHEAD
               NCHAR = lnblnkn(VARSTAR(ISTAR))
               WRITE(NOUT,91) '  ',VARSTAR(ISTAR)(1:NCHAR),
     &                        ' From Col: ',NCOLSTAR(ISTAR)
91             FORMAT(A,A,T40,A,I0)
            ENDIF
         ENDDO
         IF (NCOLSTAR(ISTAR) ==  0) THEN
            CALL ERRT(102, 'NO DOC REGISTER FOR STAR VARIABLE',ISTAR)
            GOTO 9999
         ENDIF
      ENDDO
      WRITE(NOUT,*) ' '

C     GET LIST OF SPIDER REGISTERS TO BE FILLED  -----------------
      NVARSPI = NVARSTAR
      IHI     = HUGE(NSPI)
      CALL RDPRAI(NCOLSPI,NVARSMAX,NVARSPI, 1,IHI,
     &        'CORRESPONDING SPIDER DOC REGISTER NUMBERS', CDUM,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      IF (NVARSPI .NE.  NVARSTAR) THEN
         CALL ERRT(102,
     &     'NO. OF DOC REGISTERS MUST = NUMBER OF STAR VARIABLES',
     &      NVARSTAR)
         GOTO 9999
      ENDIF
      WRITE(NOUT,*) ' '

C     GET LIST OF STAR FILENAME VARIABLES WANTED -----------------
      IRTFLG = -999   ! DO NOT UPPERCASE RECLIN
      CALL RDPRMC(RECLIN,NLET,.TRUE., 
     &  'STAR FILE NAME VARIABLES WANTED (e.g. _rlnImageName)',
     &  CDUM,IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      NVARSPINAME  = 0
      NVARSTARNAME = 0
      IF (NLET > 0 .AND. RECLIN(1:1) .NE. '*') THEN
C        PARSE STAR NAME VARIABLES LINE INTO ARRAY
         CALL GET_TOKENS(RECLIN,NLET,MAXLENVAR,NVARSMAX,
     &                .FALSE.,' ,', VARSTARNAME, NVARSTARNAME, IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
 
C        ASSOCIATE STAR FILENAME VARIABLES WITH COLS IN STAR FILE ----
         NVARSPINAME = NVARSTARNAME

         DO ISTAR = 1,NVARSTARNAME
            DO IHEAD = 1,NVARHEAD
               IF (VARSTARNAME(ISTAR) == VARHEAD(IHEAD)) THEN
                  NCOLSTARNAME(ISTAR) = IHEAD
                  NCHAR               = lnblnkn(VARSTARNAME(ISTAR))
                  WRITE(NOUT,91) '  ',VARSTARNAME(ISTAR)(1:NCHAR),
     &                        ' From Col: ',NCOLSTARNAME(ISTAR)
              ENDIF
            ENDDO
         ENDDO
         WRITE(NOUT,*) ' '
 
C        GET LIST OF SPIDER REGISTERS TO BE FILLED  -----------------
         IHI         = HUGE(NVARSPINAME)
         NVARSPINAME = NVARSMAX
         NCOLSPINAME = 0   ! ARRAY ZERO
         CALL RDPRAI(NCOLSPINAME,NVARSMAX,NVARSPINAME, 1,IHI,
     &       'CORRESPONDING SPIDER DOC REGISTER NUMBERS', CDUM,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

      ENDIF
      WRITE(NOUT,*) ' '

C     OPEN SPIDER OUTPUT DOC FILE
      ADDEXT  = .TRUE.
      GETNAME = .TRUE.
      ISOLD   = .FALSE.
      APPEND  = .TRUE.
      MESSAGE = .TRUE.
      IRTFLG  = -8         ! NO IC USE

      CALL OPENDOC(DOCNAM,ADDEXT,NLET,LUNDOC,LUNDOCT,GETNAME,
     &           'SPIDER DOC',ISOLD,APPEND,MESSAGE,
     &            NEWFILE,IRTFLG)
      IF (IRTFLG == -1) RETURN


      REWIND(LUNXM)           ! RETURN TO BEGINNING

      NREGSTAR = 0
      NLIST          = MAXVAL(NCOLSPI(1:NVARSPI))
      I              = MAXVAL(NCOLSPINAME(1:NVARSPINAME))
      NLIST          = MAX(NLIST,I)
      DLIST(1:NLIST) = 0.0
      IKEY           = 0

      ILABEL = 0  ! ARRAY ZERO
      DO I = 1,NVARSPI
         ILOC         = NCOLSPI(I)
         ILABEL(ILOC) = NCOLSTAR(I)
         !write(6,*) ' Col:',iloc,' Label#:',ilabel(iloc)
      ENDDO

      !write(6,*) 'Nvarspi:',    Nvarspi 
      !write(6,*) 'Nvarspiname:',Nvarspiname
      !write(6,*) 'Ncolstarname(1):', ncolstarname(1) 

      DO I = 1,NVARSPINAME             !  1:3
         ILOC         = NCOLSPINAME(I) !  1,2,3
         ILABEL(ILOC) = NCOLSTARNAME(1)  ! WRONG IF > 1 namevar!!!
         !write(6,*) 'Ncolstarname:',i, ncolstarname(i) 
      ENDDO

      !write(6,'("star col:",100(i5,9x))') ilabel(1:nlist)

      WRITE(RECLIN,'("STAR COL:",100(I5,9X))') ILABEL(1:NLIST)
      NCHAR = lnblnkn(RECLIN)
      CALL LUNDOCPUTCOM(LUNDOCT,RECLIN(:NCHAR),IRTFLG)

      NBAD    = 0

C     LOOP OVER ALL DATA LINES IN STAR FILE
      DO
 
C        READ SINGLE LINE FROM STAR FILE
         READ(LUNXM,'(A)',IOSTAT=IRTFLG) RECLIN
         IF (IRTFLG .NE. 0) RETURN

         !write(6,*) '  reclin: ', reclin(:nlet)
         RECLIN = ADJUSTL(RECLIN)
         NLET   = lnblnkn(RECLIN)
         IF (NLET <= 0) CYCLE         ! SKIP BLANK LINE

         IGOT = SCAN(RECLIN(1:1),'_dl#!')
         IF (IGOT > 0) CYCLE         ! SKIP HEADER LINE

C        PARSE STAR FILE REGISTER LINE INTO TOKENS
         CALL GET_TOKENS(RECLIN,NLET,MAXLENVAR,NVARSMAX,.FALSE.,' ',
     &                   VARTMP, NREGSTAR, IRTFLG)
         IF (IRTFLG .NE. 0) EXIT
         IF (NREGSTAR < 1)  CYCLE        ! SKIP COMMENT LINE

C        INTERPRET TOKENS AS FLOATS, PLACE IN DLIST
         DO ITOK = 1,NVARSPI
            ILOC   = NCOLSPI(ITOK)

            !write(6,*) ' ncolstar:',ncolstar(1:20)
            !write(6,*) ' itok,iloc:',itok,iloc,ncolstar(itok)
            CTOKEN = TRIM(VARTMP(NCOLSTAR(ITOK)))

            !write(6,*) ' Token:',ctoken
            READ(CTOKEN,'(f10.0)') DLIST(ILOC) 
            !write(6,*) ' Token:',ctoken,' dlist:',dlist(ITOK)
         ENDDO

C        EXTRACT ALL INTEGERS FIELDS FROM FILENAME TOKENS

         ! NOTE: NVARSTARNAME MAY BE ZERO
         IT      = 1
         FORMOUT = '(I7,1X,I2,1X,50(G14.7,1X))'

         !write(6,*) ' nvarstarname:',nvarstarname

         DO ITOK = 1,NVARSTARNAME
            IGO    = 0
            RECNAM = ADJUSTL(TRIM(VARTMP(NCOLSTARNAME(ITOK))))
            NCHAR  = lnblnkn(RECNAM)
            !write(6,*) ' recnam:',recnam(1:nchar)

C           PARSE STAR FILE TOKEN INTO SUB-TOKENS
            CALL GET_TOKENS_INTS(RECNAM,NCHAR,MAXLENVAR,NVARSMAX,
     &                     .TRUE.,VALID, VARTMP, NREGSTAR, IRTFLG)
            IF (IRTFLG .NE. 0) EXIT
            IF (NREGSTAR < 1)  EXIT        ! SKIP 
            !write(6,*) ' vartmp(1):', adjustl(trim(vartmp(1)))
            !write(6,*) ' nregstar:',nregstar,'  it:',it

C           INTERPRET TOKENS AS FLOATS, PLACE IN DLIST
            NREGSTAR = MIN(NREGSTAR,(NVARSPINAME - IT +1))
            !write(6,*)' nregstar:',nregstar

            IPUT = 0
            DO I = IT,IT+NREGSTAR-1
               IWHICH = I - IT + 1
               ILOC   = NCOLSPINAME(I)
               CTOKEN = ADJUSTL(VARTMP(IWHICH))
               NCHAR  = lnblnkn(CTOKEN)

               READ(CTOKEN(1:NCHAR),'(I)') IDLIST
               DLIST(ILOC) = IDLIST
 
               !write(6,*)' Token:',ctoken(1:nchar),i,iloc,idlist

               IPUT = IPUT + 1
            ENDDO
            IT = IT + IPUT

         ENDDO       ! END OF: DO I = 1,NVARSTARNAME

C        WRITE OUTPUT LINE IN DOC FILE
         IKEY = IKEY + 1
         !CALL LUNDOCWRTDAT(LUNDOCT,IKEY,DLIST,NLIST,IRTFLG)
         CALL LUNDOCWRTDATF(LUNDOCT,IKEY,DLIST,NLIST,FORMOUT,IRTFLG)
         !write(6,*)'  dlsit:',dlist(1:7)

      ENDDO          ! END OF: DO 

 
9999  CLOSE(LUNDOC)

C     PUT NUMBER OF KEYS IN OPERATION LINE REGISTER
      CALL REG_SET_NSEL(1,1,FLOAT(IKEY),0,0,0,0,IRTFLG)

      IRTFLG =  0

      END
 


