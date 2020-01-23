
C++*********************************************************************
C
C UTIL2SUP.F   NEW                              8/1/97    ArDean Leith  
C              REWRITTEN                        MAR 99    ArDean Leith
C              USED REDVOL                      DEC 2000  ArDean Leith
C              USED OPFILEC                     FEB 2003  ArDean Leith
C              ADDFAC                           MAR 2003  ArDean Leith
C              GETNEWSTAK PARAMS.               OCT 2010  ArDean Leith
C              REAL * FOURIER ALLOWED           JUL 2011  ArDean Leith
C              LUN23 != 23 CLOSURE              APR 2014  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  UTIL2SUP (PROMPT1,PROMPT2,PROMPT3, LUN1,LUN2,LUN3, SIGN)
C  UTIL2SUPL(PROMPT1,PROMPT2,PROMPT3, SIGN) 
C
C  PARAMETERS:      
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE UTIL2SUP(PROMPT1,PROMPT2,PROMPT3,
     &                      LUN1,LUN2,LUN3,
     &                      SIGN) 

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        CHARACTER (LEN=*)      :: PROMPT1,PROMPT2,PROMPT3
        INTEGER                :: LUN1,LUN2,LUN3
        REAL                   :: SIGN

        INTEGER                :: LUNA,LUNB,LUNOUT
        CHARACTER(LEN=MAXNAM)  :: FILNAM1,FILNAM2,FILNAM3

        CHARACTER (LEN=1)      :: NULL = CHAR(0)
        REAL, ALLOCATABLE      :: VOLBUF(:)
        LOGICAL                :: ASKNAME,MUSTGET,WANTNEXT
        LOGICAL                :: BARE1,BARE2,ISBARE

C       IN CASE LUN1,... ARE CONSTANTS

        LUNA   = LUN1
        LUNB   = LUN2
        LUNOUT = LUN3

        IPVOL  = 0

        CALL FILERD(FILNAM1,NLETI,NULL,PROMPT1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL FILERD(FILNAM2,NLETI,NULL,PROMPT2,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        BARE1 = ISBARE(FILNAM1)
        BARE2 = ISBARE(FILNAM2)

          
        IF (SIGN >= 1000) THEN
           CALL RDPRM2S(FACT1,FACT2,NOT_USED,
     &                  'FACTORS FOR FIRST & SECOND FILES',IRTFLG)
	   IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

        IF (.NOT.(BARE1 .OR. BARE2)) THEN

           MAXIM1 = 0
           MAXIM2 = 0
           MAXIM3 = 0

C          NOT A STACKS OPERATION, OPEN FIRST INPUT FILE ON LUNA
           CALL OPFILEC(0,.FALSE.,FILNAM1,LUNA,'O',IFORM1,
     &                  NX,NY,NZ,
     &                 MAXIM1,PROMPT1,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

C          ALLOCATE SPACE IN VOLBUF
           ALLOCATE(VOLBUF(NX*NY*NZ), STAT=IRTFLGT)
           IF (IRTFLGT .NE. 0) THEN
              MWANT = NX*NY*NZ
              CALL ERRT(46,'UTIL2SUP; VOLBUF',MWANT)
              GOTO 9999
           ENDIF

C          LOAD VOLUME FROM FIRST FILE INTO VOLBUF
           CALL REDVOL(LUNA,NX,NY,1,NZ,VOLBUF,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

           CLOSE(LUNA)
           ITER = 0

C          OPEN 2ND... INPUT FILE ON LUNB

10         ITER    = ITER + 1
           ASKNAME = (ITER > 1)
           CALL OPFILEC(0,ASKNAME,FILNAM2,LUNB,'O',IFORM2,
     &                NXT,NYT,NZT,MAXIM2,PROMPT2,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998
           
           CALL SIZCHK(NULL,NXT,NYT,NZT, 0,
     &                      NX ,NY ,NZ , 0, IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

           !IF (IFORM2 .NE. IFORM1) CALL ERRT(40,'UTIL2SUP',NE)

           ! KLUDGE TO ALLOW MANIPULATION OF FOURIER FILES
           IFORMT = IFORM1
           IF ( IFORM1 .NE. IFORM2 .AND.
     &        ((IFORM1 < 0 .AND. IFORM2 > 0) .OR.
     &         (IFORM1 > 0 .AND. IFORM2 < 0))) THEN
              IFORMT = MAX(IFORM1,IFORM2)
           ENDIF
       
           IF (SIGN < 1000) THEN
C             ADD, ETC SECOND FILE TO STORED VOLUME
              CALL ADD(VOLBUF,LUNB,IFORMT,NX,NY,NZ,SIGN)
           ELSE
C             CARRY OUT ADDITION, ETC
              CALL ADDFAC(VOLBUF,LUNB,IFORMT,NX,NY,NZ,SIGN,
     &                    FACT1,FACT2)
           ENDIF

C          CLOSE SECOND FILE (IN CASE OUTPUT IS SAME FILE)
           CLOSE(LUNB)

C          OPEN OUTPUT FILE ON LUNOUT
           ASKNAME = (ITER <= 1)
           CALL OPFILEC(LUNA,ASKNAME,FILNAM3,LUNOUT,'U',IFORM1,
     &                  NX,NY,NZ,MAXIM3,PROMPT3,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

C          PUT SUM,ETC. IN OUTPUT FILE ON LUNOUT
           CALL WRTVOL(LUNOUT,NX,NY, 1,NZ,VOLBUF,IRTFLG)

           CLOSE(LUNOUT)

C          CONTINUE UNTIL '*' IS INPUT
           IF (SIGN < 1000) GOTO 10 

        ELSE
C          STACKS OPERATION

           MAXIM1 = -1
           MAXIM2 = -1
           MAXIM3 = -1
           IMGNUM = 0

C          OPEN FIRST INPUT FILE ON LUNA
           CALL OPFILEC(0,.FALSE.,FILNAM1,LUNA,'O',IFORM1,
     &                 NX,NY,NZ,
     &                 MAXIM1,PROMPT1,.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

C          ALLOCATE SPACE IN VOLBUF
           ALLOCATE(VOLBUF(NX*NY*NZ), STAT=IRTFLGT)
           IF (IRTFLGT .NE. 0) THEN
              CALL ERRT(46,'UTIL2SUP; VOLBUF',NX*NY*NZ)
              GOTO 9999
           ENDIF

C          OPEN SECOND INPUT FILE ON LUNB (IF NECESSARY)
           IF (FILNAM2 == FILNAM1) THEN
              LUNB = LUNA
           ELSE
              CALL OPFILEC(0,.FALSE.,FILNAM2,LUNB,'O',IFORM2,
     &              NXT,NYT,NZT,MAXIM2,PROMPT2,.TRUE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9998

              CALL SIZCHK(NULL,NXT,NYT,NZT, 0,
     &                         NX ,NY ,NZ , 0, IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9998

              !IF (IFORM2 .NE. IFORM1)  CALL ERRT(40,'UTIL2SUP',NE)
           ENDIF

C          FIND TOTAL NUMBER OF COMMON IMAGES IN STACKS
           IF (MAXIM1 > 0) NIMAGE = MAXIM1
           IF (MAXIM2 > 0) NIMAGE = MAXIM2
           IF (MAXIM1 > 0 .AND. MAXIM2 > 0)
     &         NIMAGE = MIN(MAXIM1,MAXIM2)
 
C          FIND OUTPUT STACK NAME
           CALL FILERD(FILNAM3,NLETI,NULL,PROMPT2,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

C          OPEN OUTPUT STACK ON LUNOUT (IF NECESSARY)
           IF (FILNAM3 == FILNAM2) THEN
              LUNOUT = LUNB
           ELSEIF (FILNAM3 == FILNAM1) THEN
              LUNOUT = LUNA
           ELSE
C             OUTPUT IS DIFFERENT STACK FROM EITHER INPUT
              CALL OPFILEC(LUNA,.FALSE.,FILNAM3,LUNOUT,'U',IFORM1,
     &              NX,NY,NZ,MAXIM3,PROMPT3,.TRUE.,IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9998
           ENDIF

C          IF FIRST FILE IS NOT A STACK CAN SKIP TO NEXTIMAGE
           WANTNEXT = (MAXIM1 < 0)
 
C          IF FIRST FILE IS STACK MUSTGET SPECIFIED IMGNUM
           MUSTGET = .NOT. WANTNEXT 

 
20         IMGNUM = IMGNUM + 1
           IF (VERBOSE) WRITE(NOUT,*) ' '

C          GET IMGNUM FROM FIRST INPUT
           CALL GETOLDSTACK(LUNA,NX,IMGNUM,.TRUE.,.FALSE.,
     &                     .TRUE.,IRTFLG)
           IF (IRTFLG > 0 .AND. IMGNUM <= 1) GOTO 9998
           IF (IRTFLG > 0) GOTO 9998
        
C          LOAD VOLUME FROM FIRST INPUT FILE INTO VOLBUF 
           CALL REDVOL(LUNA,NX,NY,1,NZ,VOLBUF,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9998

C          WHAT HAPPENS IF STACK + STACK BUT ONE STACKED IMAGE
C          NOT PRESENT?? CURRENTLY PART OF STACK COULD BE 
C          INCONSISTENT WITH PREVIOUS PART OF STACK SINCE IT WILL
C          HALT ON ERROR!!!!!!!!!!!!!!al

           CALL GETOLDSTACK(LUNB,NX,IMGNUM,WANTNEXT,MUSTGET,
     &                     .TRUE.,IRTFLG)
           IF (IRTFLG > 0) GOTO 9998

           ! KLUDGE TO ALLOW MANIPULATION OF FOURIER FILES
           IFORMT = IFORM1
           IF ( IFORM1 .NE. IFORM2 .AND.
     &        ((IFORM1 < 0 .AND. IFORM2 > 0) .OR.
     &         (IFORM1 > 0 .AND. IFORM2 < 0))) THEN
              IFORMT = MAX(IFORM1,IFORM2)
           ENDIF

           IF (SIGN < 1000) THEN
C             CARRY OUT ADDITION, ETC
              CALL ADD(VOLBUF,LUNB,IFORMT,NX,NY,NZ,SIGN)
           ELSE
C             CARRY OUT ADDITION, ETC
              CALL ADDFAC(VOLBUF,LUNB,IFORMT,NX,NY,NZ,SIGN,
     &                    FACT1,FACT2)
           ENDIF

c          STACK OPERATION, POINT TO NEXT STACKED IMAGE
           CALL GETNEWSTACK(LUNA,LUNOUT,.FALSE.,NX,IMGNUM,IRTFLG)
           IF (IRTFLG .GT. 0) GOTO 9998

C          PUT SUM,ETC. IN OUTPUT FILE ON LUNOUT
           CALL WRTVOL(LUNOUT,NX,NY, 1,NZ,VOLBUF,IRTFLG)

C          CONTINUE UNTIL LAST STACKED IMAGE REACHED, DO NOT CLOSE 

           IF (IMGNUM < NIMAGE) GOTO 20
        ENDIF


C       DEALLOCATE VOLBUF
9998    IF (ALLOCATED(VOLBUF)) DEALLOCATE(VOLBUF)

9999    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUN3)


        RETURN
        END


C++*********************************************************************
C
C UTIL2SUPL.F    FROM UTIL2SUP                  OCT 2012  ArDean Leith  
C                IMAGE/STACK SERIES SUPPORT     OCT 2012  ArDean Leith
C                OVERWRITING LOGIC              OCT 2014  ArDean Leith
C
C **********************************************************************
C
C  UTIL2SUPL(PROMPT1,PROMPT2,PROMPT3, SIGN) 
C
C  PARAMETERS:      
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE UTIL2SUPL(PROMPT1,PROMPT2,PROMPT3,
     &                     SIGN) 

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC' 

      CHARACTER (LEN=*)      :: PROMPT1,PROMPT2,PROMPT3
      INTEGER                :: IT
      REAL                   :: SIGN

      CHARACTER(LEN=MAXNAM)  :: FILPAT1,FILPAT2,FILPAT3

      REAL,    ALLOCATABLE   :: VOLBUF(:,:,:),VOLBUF1(:,:,:)
      INTEGER, ALLOCATABLE   :: ILIST1(:),ILIST2(:),ILIST3(:)
      LOGICAL                :: ISBAREIN,ISBAREOUT
      INTEGER                :: LOCAT1,LOCAT2,LOCAT3
      INTEGER                :: LOCAST1,LOCAST2,LOCAST3
      INTEGER                :: NILMAX,NDUM
      INTEGER                :: MAXIM1,MAXIM2,MAXIM3,NOT_USED
      INTEGER                :: NLET1, NLET2, NLET3
      INTEGER                :: ITYPE1,ITYPE2,ITYPE3,ITYPE
      INTEGER                :: NINDX1,NINDX2,NINDX3
      INTEGER                :: NLIST1,NLIST2,NLIST3
      LOGICAL                :: ISTACK1,ISTACK2,ISTACK3
      INTEGER                :: IMGNUM1,IMGNUM2,IMGNUM3
      INTEGER                :: NX, NY, NZ
      INTEGER                :: NX2,NY2,NZ2
      INTEGER                :: I,ITER, LUNCP,ILOC,MWANT,IRTFLG
      INTEGER                :: LUN2T,LUN3T
      REAL                   :: FACT1,FACT2
      LOGICAL                :: BARE1,BARE2,BARE3
      LOGICAL                :: SAME1,SAME2,SAME3
      LOGICAL                :: OVERWRITE1,OVERWRITE2
      LOGICAL                :: ISOPEN 

      CHARACTER (LEN=1)      :: NULL = CHAR(0)

      LOGICAL                :: FOUROK
      LOGICAL, PARAMETER     :: ASKNAME  = .TRUE.
      LOGICAL, PARAMETER     :: ASKNAME3 = .FALSE.

      INTEGER, PARAMETER     :: LUN1    = 21
      INTEGER, PARAMETER     :: LUN2    = 22
      INTEGER                :: LUN3    = 23
      INTEGER, PARAMETER     :: LUNDOC  = 81
      INTEGER, PARAMETER     :: LUNDOC2 = 82
      INTEGER, PARAMETER     :: LUNDOC3 = 83
      INTEGER, PARAMETER     :: LUNXM1  = 84
      INTEGER, PARAMETER     :: LUNXM2  = 85
      INTEGER, PARAMETER     :: LUNXM3  = 86

      LUN3 = 23

      NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
      ALLOCATE(ILIST1(NILMAX),
     &         ILIST2(NILMAX),
     &         ILIST3(NILMAX),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'UTIL2SUPL; ILIST1....',3*NILMAX)
          RETURN
      ENDIF

      FOUROK = .TRUE.
      IF (SIGN >= 1000) THEN
         FOUROK = .FALSE.
      ENDIF
 
C     OPEN FIRST INPUT IMAGE(S)
      CALL OPFILES(0,LUN1,LUNDOC,LUNXM1, ASKNAME,
     &             FILPAT1,NLET1, 'O',
     &             ITYPE1,NX,NY,NZ,MAXIM1,
     &             PROMPT1,
     &             FOUROK, ILIST1,NILMAX, 
     &             NOT_USED,NLIST1,IMGNUM1, IRTFLG) 
      IF (IRTFLG .NE. 0) RETURN
      LOCAT1  = INDEX(FILPAT1,'@')
      LOCAST1 = INDEX(FILPAT1,'*')
      ISTACK1 = (MAXIM1 > 0)                     ! USING A STACK
      BARE1   = (LOCAT1 > 0 .AND. LOCAST1 == 0)  ! BARESTACK
      SAME1   = (LOCAST1 == 0 .AND. .NOT. BARE1) ! IMAGE IS CONSTANT
      !write(6,*)'1 same1:',LOCAST1,BARE1,SAME1

      IF (NLIST1 > 0) ILIST2 = ILIST1            ! COPY STACK FOR ILIST2

      !write(6,*)'1 maxim,nlist,num:',maxim1,nlist1,imgnum1,filpat1(:11)
      !write(6,*) ' opened lun1:',lun1,filpat1(1:20)

C     OPEN SECOND INPUT IMAGE(S)
      CALL OPFILES(0,LUN2,LUNDOC,LUNXM2, ASKNAME,
     &             FILPAT2,NLET2, 'O',
     &             ITYPE2,NX2,NY2,NZ2,MAXIM2,
     &             PROMPT2,
     &             FOUROK, ILIST2,NILMAX, 
     &             NOT_USED,NLIST2,IMGNUM2, IRTFLG) 
      IF (IRTFLG .NE. 0) GOTO 9999
      LOCAT2  = INDEX(FILPAT2,'@')
      LOCAST2 = INDEX(FILPAT2,'*')
      ISTACK2 = (MAXIM2 > 0)                     ! USING A STACK
      BARE2   = (LOCAT2 > 0 .AND. LOCAST2 == 0)  ! BARESTACK
      SAME2   = (LOCAST2 == 0 .AND. .NOT. BARE2) ! IMAGE IS CONSTANT

c     write(6,*)' 2 maxim,nlist,num:', maxim2,nlist2,imgnum2,filpat2(:11)
c     write(6,*)' 2 at2,ast2,stk2,bar2,sam2:',
c     &              locat2,locast2,istack2,bare2,same2
c     write(6,*) ' opened lun2:',lun2,filpat2(1:20)
           
      CALL SIZCHK(NULL,NX, NY, NZ, ITYPE1, 
     &                 NX2,NY2,NZ2,ITYPE2, IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      IF (SIGN >= 1000.0) THEN
         CALL RDPRM2S(FACT1,FACT2,NOT_USED,
     &                'FACTORS FOR FIRST & SECOND FILES',IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
      ENDIF

C     FIND OUTPUT IMAGE NAME
      CALL FILERD(FILPAT3,NLET3,NULL,PROMPT3,IRTFLG)
      IF ( IRTFLG .NE. 0 ) GOTO 9999

      LOCAT3     = INDEX(FILPAT3,'@')
      LOCAST3    = INDEX(FILPAT3,'*')
      ISTACK3    = (LOCAT3 > 0)                     ! USING A STACK
      BARE3      = (LOCAT3 > 0 .AND. LOCAST3 == 0)  ! BARESTACK
      SAME3      = (LOCAST3 == 0 .AND. .NOT. BARE3) ! IMAGE IS CONSTANT

      OVERWRITE1 = .FALSE.
      OVERWRITE2 = .FALSE.

      IF (LOCAT3 == 0 .AND. LOCAST3 == 0 .AND.
     &    NLET3 == NLET1  .AND.
     &    FILPAT3(1:NLET3) == FILPAT1(1:NLET1) ) THEN
C         SIMPLE OUTPUT FILE - OVERWRITING INPUT FILE
       
         WRITE(NOUT,*) ' WARNING - OVERWRITING INPUT FILE: ', 
     &                   FILPAT1(1:NLET1)

         OVERWRITE1 = .TRUE.  
         LUN3       = LUN1
         NLIST3     = NLIST1

       ELSEIF(LOCAT3 == 0 .AND. LOCAST3 == 0 .AND.
     &    NLET3 == NLET2  .AND.
     &    FILPAT3(1:NLET3) == FILPAT2(1:NLET2) ) THEN
C         SIMPLE OUTPUT FILE - OVERWRITING INPUT FILE 
       
         OVERWRITE2 = .TRUE.  
 
         WRITE(NOUT,*) ' WARNING - OVERWRITING INPUT FILE: ', 
     &                   FILPAT2(1:NLET2)

         OVERWRITE2 = .TRUE.  
         LUN3       = LUN2
         NLIST3     = NLIST2

      ELSE

C        OPEN OUTPUT IMAGE(S)
         ITYPE3  = ITYPE1       ! IMAGE TYPE
         MAXIM3  = -1           ! ALLOW BARE STACK
         IMGNUM3 = IMGNUM1      ! IMAGE # WANTED
         LUNCP   = LUN1

         CALL OPFILES(LUNCP,LUN3,LUNDOC,LUNXM3,ASKNAME3,
     &             FILPAT3,NLET3, 'U',
     &             ITYPE3,NX,NY,NZ,MAXIM3,
     &             FILPAT3,
     &             FOUROK, ILIST3,NILMAX, 
     &             NOT_USED,NLIST3,IMGNUM3, IRTFLG) 
C        IRTFLG = -2, FILE ALREADY OPEN
         IF (IRTFLG .NE. 0 .AND. IRTFLG .NE. -2) GOTO 9999

         !write(3,*)' filpat3:',filpat3(1:nlet3)
         !write(3,*)' nlist3,ilist3: ',nlist3,ilist3(1:nlist3)
         !write(3,*)' nilmax,imgnum3:',nilmax,imgnum3

         IF (NLET3 == NLET1  .AND.
     &       FILPAT3(1:NLET3) == FILPAT1(1:NLET1) .AND.
     &       NLIST3 == NLIST1) THEN
       
            OVERWRITE1 = .TRUE.  
            DO I=1,NLIST3
               IF (ILIST1(I) .NE. ILIST3(I)) THEN
                  OVERWRITE1 = .FALSE.
                  EXIT
               ENDIF
            ENDDO

C           write(6,*) ' nlet1,nlet3,nlist1,nlist3:',nlet1,nlet3,nlist1,nlist3
            IF (OVERWRITE1) THEN
C              DUPLICATE FILE NAME!
               WRITE(NOUT,*) ' WARNING - OVERWRITING INPUT FILE: ', 
     &                        FILPAT1(1:NLET1)
               LUN3   = LUN1
               NLIST3 = NLIST1
            ENDIF

         ELSEIF (NLET3 == NLET2  .AND.
     &           FILPAT3(1:NLET3) == FILPAT2(1:NLET2) .AND.
     &           NLIST3 == NLIST2) THEN

c           write(6,*) ' nlet2,nlet3,nlist2,nlist3:',nlet2,nlet3,nlist2,nlist3
      
            OVERWRITE2 = .TRUE.  
            DO I=1,NLIST2
               IF (ILIST2(I) .NE. ILIST3(I)) THEN
                  OVERWRITE2 = .FALSE.
                  EXIT
               ENDIF
            ENDDO

            IF (OVERWRITE2) THEN
C              DUPLICATE FILE NAME!
               WRITE(NOUT,*) ' WARNING - OVERWRITING INPUT FILE: ', 
     &                       FILPAT2(1:NLET2)
               LUN3   = LUN2
               NLIST3 = NLIST2
            ENDIF
         ENDIF
      ENDIF

      !write(6,*) ' opened lun3:',lun3,filpat3(1:20)
      !write(6,*)'3 maxim,nlist,num:',maxim3,nlist3,imgnum3,filpat3(:11)
      !write(6,*) ' now lun3:',lun3,nlet2,nlet3,overwrite2


C     ALLOCATE SPACE IN VOLBUF
      ALLOCATE(VOLBUF(NX,NY,NZ), STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          MWANT = NX*NY*NZ
          CALL ERRT(46,'UTIL2SUPL; VOLBUF',MWANT)
          GOTO 9999
      ENDIF

      IF (SAME1) THEN 
C        LOAD VOLUME FROM FIRST FILE INTO VOLBUF
         ALLOCATE(VOLBUF1(NX,NY,NZ), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'UTIL2SUPL; VOLBUF1',MWANT)
            GOTO 9999
         ENDIF

         CALL REDVOL(LUN1,NX,NY,1,NZ,VOLBUF1,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
         CLOSE(LUN1)
      ENDIF

      ! KLUDGE TO ALLOW MANIPULATION OF FOURIER FILES
      ITYPE = ITYPE1
      IF ( ITYPE1 .NE. ITYPE2 .AND.
     &   ((ITYPE1 < 0 .AND. ITYPE2 > 0) .OR.
     &    (ITYPE1 > 0 .AND. ITYPE2 < 0))) THEN
           ITYPE = MAX(ITYPE1,ITYPE2)
      ENDIF

      ITER       = 0
      NINDX1     = 1
      NINDX2     = 1
      NINDX3     = 1

        LUN2T = LUN2
        IF (SAME2 .OR. OVERWRITE2) THEN
C          DO NOT OPEN NEXT SECOND INPUT FILE 
           LUN2T = 0
        ENDIF
        LUN3T = LUN3
        IF (OVERWRITE1 .OR. OVERWRITE2) THEN
C          DO NOT OPEN NEXT THIRD INPUT FILE 
           LUN3T = 0
        ENDIF


      DO 

      !write(3,*)' In util2supl, num,indx:',imgnum1,nindx1,nindx2,nindx3

        ITER = ITER + 1

C       LOAD VOLUME FROM FIRST FILE INTO VOLBUF
        IF (SAME1) THEN
c          REUSE SAME FIRST INPUT BUFFER
           VOLBUF = VOLBUF1 
        ELSE 
           !write(6,*)' Loading buf1,lun1:',lun1
           CALL REDVOL(LUN1,NX,NY,1,NZ,VOLBUF,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (.NOT. ISTACK1 .AND..NOT. BARE1) CLOSE(LUN1)
        ENDIF

        IF (SIGN < 1000) THEN
C          COMBINE FIRST & SECOND FILE --> VOLBUF
           !write(6,*)' Adding bufs,lun2:',lun2

           CALL ADD(VOLBUF,LUN2, ITYPE, NX,NY,NZ, SIGN)
        ELSE
C          ADD FIRST & SECOND FILE --> VOLBUF WITH FACTORS
           CALL ADDFAC(VOLBUF,LUN2, ITYPE, NX,NY,NZ, SIGN,
     &                 FACT1,FACT2)
        ENDIF

C       CLOSE SECOND FILE (IN CASE OUTPUT IS SAME FILE)
        IF (.NOT. ISTACK2 .AND. 
     &      .NOT. BARE2   .AND. 
     &      .NOT. SAME2)  CLOSE(LUN2)

C       PUT SUM, ETC. IN OUTPUT FILE ON LUN3
        CALL WRTVOL(LUN3,NX,NY, 1,NZ,VOLBUF,IRTFLG)

        IF ((NINDX1 >= NLIST1 .AND. NINDX2 >= NLIST2)) then
           !write(3,*)' nindx1 >= nlist1:', nindx1,' >= ',nlist1
           !write(3,*)' nindx2 >= nlist2:', nindx2,' >= ',nlist2
           EXIT      ! END OF INPUT LIST
        ENDIF
      
        !write(3,*)' In util2sup, same1,overwrite1:',same1,overwrite1
        IF (.NOT. SAME1) THEN
C          OPEN NEXT FIRST INPUT FILE 
 
           CALL NEXTFILE(NINDX1, ILIST1, 
     &                   FOUROK,LUNXM1,
     &                   NLIST1,MAXIM1,   
     &                   LUN1,LUNCP,
     &                   FILPAT1,'O',
     &                   IMGNUM1, IRTFLG) 
c           write(3,'(A,5i6)')
c     &         ' In util2supl 1, num1,list1,indx1,irtflg:',
c     &                         imgnum1,nlist1,nindx1,irtflg
           IF (IRTFLG < 0)  EXIT      ! END OF INPUT FILES
        ENDIF

c       write(6,*)' same2,overwrite2:',same2,overwrite2

        IF (.NOT. SAME2) THEN
C          OPEN NEXT SECOND INPUT FILE 
           !write(6,*)' same2:',same2,nindx2,nlist2,filpat2(1:10)

           CALL NEXTFILE(NINDX2, ILIST2, 
     &                   FOUROK,LUNXM2,
     &                   NLIST2,MAXIM2,   
     &                   LUN2,LUNCP,
     &                   FILPAT2,'O',
     &                   IMGNUM2, IRTFLG) 

c           write(3,'(A,5i6)')
c     &         ' In util2supl 2  num2,nlist2,nindx2,irtflg:',
c     &                        imgnum2,nlist2,nindx2,irtflg
           IF (IRTFLG < 0)  EXIT       ! END OF INPUT FILES
        ENDIF


        !write(6,*)' same3,overwrite1 2:',same3,overwrite1,overwrite2
        IF (.NOT. SAME3 .AND. 
     &      .NOT. OVERWRITE1 .AND. 
     &      .NOT. OVERWRITE2) THEN

C          OPEN NEXT OUTPUT FILE 
           IF (BARE3) THEN
              NLIST3  = NLIST1
              IF (BARE2) NLIST3 = MAX(NLIST1,NLIST2)
           ENDIF

c          write(3,*)' In util2sup 3, bare3,nindx3,nlist3: ',
c     &                               bare3,nindx3,nlist3
c          write(3,*)' In util2sup, nlist3,ilist3(1..),maxim3,imgnum3:',
c     &                            nlist3,ilist3(1:nlist3),maxim3,imgnum3
c          write(3,*)' In util2sup,filpat3: ',filpat3(1:10)
c          !!if (nindx1 > ngot1 == nlist3) then overflow

           CALL NEXTFILE(NINDX3, ILIST3, 
     &                   FOUROK,LUNXM1,
     &                   NLIST3,MAXIM3,   
     &                   LUN3,LUNCP,
     &                   FILPAT3,'N',
     &                   IMGNUM3, IRTFLG) 

c           write(3,'(A,5i6)')
c     &         ' In util2supl 3  num3,nlist3,nindx3,irtflg:',
c     &                        imgnum3,nlist3,nindx3,irtflg

           IF (BARE3) IMGNUM3 = IMGNUM3 + 1   

         ENDIF

         IF (IRTFLG == -99) THEN
             CALL ERRT(102,'INSUFFICIENT OUTPUT FILE NAMES',NLIST3)
             EXIT         
         ELSEIF (IRTFLG < 0) THEN
             EXIT         ! END OF INPUT FILES
         ENDIF
         IF (IRTFLG .NE. 0) GOTO 9999    ! ERROR

      ENDDO

9999  IF (ALLOCATED(VOLBUF))    DEALLOCATE(VOLBUF)
      IF (ALLOCATED(ILIST1))    DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2))    DEALLOCATE(ILIST2)
      IF (ALLOCATED(ILIST3))    DEALLOCATE(ILIST3)

      IF (VERBOSE) WRITE(NOUT,*) ' '

      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
      CLOSE(23)        ! IN CASE REDEFINED ABOVE
      CLOSE(LUNDOC)
      CLOSE(LUNXM1)
      CLOSE(LUNXM2)
      CLOSE(LUNXM3)
     
      END

#ifdef NEVER
        !inquire(unit=lun1,opened=isopen, name=filopend)
        !write(6,*) 'isopen1a:',isopen,filopend(1:12)
#endif


