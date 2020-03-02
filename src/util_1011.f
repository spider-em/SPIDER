
C++*********************************************************************
C
C  UTIL_1011.F   DERIVED FROM ROTQSS              2/14/12 ARDEAN LEITH 
C                CENT_ROD BUGGY  - REMOVED        8/08/14 ARDEAN LEITH
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
C   UTIL_1011()
C
C   INPUT/OUTPUT:  INPUT FILE,  OPTIONAl SEL FILE, 
C                  OUTPUT FILE, OPTIONAL OUTPUT DOC FILE
C
C   PURPOSE:  CENTER IMAGES
C             DRIVER FOR OPERATIONS WHICH BOTH READ AND WRITE FILES
C             CAN USE A SELECTION DOC FILE OR XMIPP SELFILE
C             CAN OPEN AN OUTPUT DOC FILE WITH ALIGNMENT PARAMETERS
C
C--*********************************************************************

      SUBROUTINE UTIL_1011()

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      REAL,    ALLOCATABLE   :: BUF1(:,:),BUF2(:,:)
      INTEGER, ALLOCATABLE   :: ILIST1(:),ILIST2(:)

      REAL                   :: DLIST(5)
      CHARACTER (LEN=MAXNAM) :: FILPATIN,FILPATOUT,FILNAM
      CHARACTER (LEN=MAXNAM) :: DOCNAM,OUTDOC,MSG

      CHARACTER (LEN=1)      :: NULL = CHAR(0)

      REAL                   :: THETA
      REAL                   :: XC,YC,ZC,OFFX,OFFY,OFFZ
      INTEGER                :: NILMAX,MAXIMOUT,ITYPE,NX,NY,NZ
      INTEGER                :: NUMB,MAXX,MAXY,NLET,NDUM,NIMGOUT
      INTEGER                :: NIMG,IMGNUM,IMGNUMOUT,MWANT
      INTEGER                :: NINDX1,NINDX2,ICOUNT,ISLICE
      INTEGER                :: MAXIMIN
      INTEGER                :: LOCAT,LOCAST,NIMGT,NLETF
      INTEGER                :: NXLD,IKEY,IDUM
      LOGICAL                :: ISBAREIN,ISBAREOUT,NEWFILE,WANTDOCOUT
      INTEGER                :: IRTFLG,NE
      INTEGER                :: LUNDOCT,NLIST
      LOGICAL                :: FOURIER = .TRUE.
      LOGICAL                :: SHIFTIT = .TRUE.

      INTEGER,PARAMETER      :: LUNIN     = 21 
      INTEGER,PARAMETER      :: LUNOUT    = 22
      INTEGER,PARAMETER      :: LUNDOCSEL = 81
      INTEGER,PARAMETER      :: LUNXM1    = 82
      INTEGER,PARAMETER      :: LUNXM2    = 83
      INTEGER,PARAMETER      :: LUNDOC    = 84

      INTEGER                :: lnblnkn
 
      NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
      ALLOCATE(ILIST1(NILMAX),
     &         ILIST2(NILMAX),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'UTIL_1011; ILIST1....',2*NILMAX)
          RETURN
      ENDIF

C     OPEN INPUT IMAGE(S)
      CALL OPFILES(0,LUNIN,LUNDOCSEL,LUNXM1, 
     &             .TRUE.,FILPATIN,NLET, 'O',
     &             ITYPE,NX,NY,NZ,MAXIMIN,
     &             NULL,
     &             .FALSE.,ILIST1,NILMAX, 
     &             NDUM,NIMG,IMGNUM, IRTFLG) 
      IF (IRTFLG .NE. 0) RETURN

      LOCAT  = INDEX(FILPATIN,'@')
      LOCAST = INDEX(FILPATIN,'*')
      IF (NIMG > 0) ILIST2 = ILIST1

      !write(6,*)  ' Input image #: ',IMGNUM,NIMG,maximin,locast
  
C     SET SINGLE IMAGE OPERATION 
      IF (NIMG <= 1) IMGNUM  = 0
      
C     OPEN OUTPUT IMAGE(S)
      MAXIMOUT  = -1           ! ALLOW BARE STACK
      IMGNUMOUT = IMGNUM       ! IMAGE # WANTED
      NIMGT     = -NIMG        ! DO NOT ASK FOR ANOTHER LIST

      CALL OPFILES(LUNIN,LUNOUT,LUNDOCSEL,LUNXM2,
     &            .TRUE.,FILPATOUT,NLET, 'U',
     &            ITYPE,NX,NY,NZ,MAXIMOUT,
     &            'OUTPUT FILE NAME OR TEMPLATE (E.G. IMG@****)~',
     &            .FALSE., ILIST2,NIMGT, 
     &            NDUM,NIMGOUT, IMGNUMOUT, IRTFLG) 

      IF (IRTFLG .NE. 0) GOTO 9999
      !write(6,*)'output Image #:',IMGNUMOUT,NIMG,maximout,filpatout(:11)

      WANTDOCOUT = .FALSE.
      IF (NIMG > 1) THEN
C        OPEN OPTIONAL OUTPUT DOC FILE (FOR APPENDING)
         WANTDOCOUT = .TRUE.
         LUNDOCT    = LUNDOC
         CALL OPENDOC(OUTDOC,.TRUE.,NLET,LUNDOC,LUNDOCT,.TRUE.,
     &           'OUTPUT ALIGNMENT DOCUMENT',.FALSE.,.TRUE.,.TRUE.,
     &            NEWFILE,IRTFLG)
         IF (IRTFLG == -1) THEN
C           DO NOT WANT OUTPUT DOC FILE
            LUNDOCT    = 0
            WANTDOCOUT = .FALSE.
         ELSEIF (IRTFLG .NE. 0) THEN
            GOTO 9999
         ENDIF
      ENDIF


      SELECT CASE(FCHAR(4:4))
      CASE ('P')  !  ------- PH --------------------------------- 'PH' 
         NXLD = NX + 2 - MOD(NX,2)
         ALLOCATE (BUF1(NXLD,NY*NZ),
     &             BUF2(NX*6,1), 
     &             STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = NXLD*NY*NZ + NX*6 
            CALL ERRT(46,'UTIL_1011; BUF1..',MWANT)
            GOTO 9999
         ENDIF  

      CASE ('R','S') !  --------------------------------- 'ROD' & 'SYM' 
         NXLD = NX + 2 - MOD(NX,2)
         ALLOCATE (BUF1(NXLD,NY),
     &             BUF2(NXLD,NY), 
     &             STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = NX*NY + NXLD*NY 
            CALL ERRT(46,'UTIL_1011; BUF1..',MWANT)
            GOTO 9999
         ENDIF  

      CASE DEFAULT
          CALL ERRT(101,'UNKNOWN CASE,PGM ERROR',IDUM)
      END SELECT

 
      IF (WANTDOCOUT) THEN
C               123456789 123456789 123456789 123456789 123456789 123
         MSG = '          IMG#         ANGLE            SX           '//
     &         ' SY          SZ' 
         CALL LUNDOCPUTCOM(LUNDOCT,MSG,IRTFLG)
      ENDIF

      THETA  = 0.0
      OFFZ   = 0     
      NINDX1 = 1
      NINDX2 = 1
      IKEY   = 1
      DO 
         !write(6,*)  'Image #:',IMGNUM

         SELECT CASE(FCHAR(4:4))

         CASE ('P') !  ------- PH -------------------------------- 'PH' 
 
            IF (NZ > 1) THEN
               CALL CENT_3PH(LUNIN,NX,NY,NZ, XC,YC,ZC)
               OFFZ  = ZC - (NZ/2+1)
            ELSE
C              READ INPUT IMAGE, NEEDS: BUF1(NX,NY)
               CALL READV(LUNIN,BUF1, NX,NY, NX,NY,1)

C              FIND CENTER OFFSET
               CALL CENT_PH(BUF1, NX,NY, XC,YC)
            ENDIF

            OFFX  = XC - (NX/2+1)  ! RELATIVE TO SPIDER CENTER
            OFFY  = YC - (NY/2+1)

            !write(6,*) ' Center:',xc,yc,offx,offy,offz

C           CENTER IMAGE/VOLUME. NEEDS: BUF1(NXLD,NY,NZ), BUF2(NX*6)

            CALL SHIFTS(FOURIER, LUNIN,LUNOUT, BUF1,BUF2,
     &                  NX,NY,NZ, -OFFX,-OFFY,-OFFZ, IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

         CASE ('S') !  ------- SYM ------------------------------ 'SYM' 

C           READ INPUT IMAGE.  NEEDS: BUF1(NXLD,NY)
            CALL READV(LUNIN,BUF1, NXLD,NY, NX,NY,1)

C           CENTER THE IMAGE. NEEDS: BUF1(NXLD,NY), BUF2(NX,NY)
            CALL CENT_SYM(BUF1, BUF2, SHIFTIT, NXLD,NX,NY, 
     &                    OFFX,OFFY, IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

C           WRITE OUTPUT IMAGE/VOLUME
            CALL WRITEV(LUNOUT,BUF2, NX,NY,NX,NY, NZ)

         CASE ('R') !  ------------------------------------------ 'ROD' 

            CALL ERRT(101,'OPERATION REMOVED DUE TO UNRESOLVED BUG',NE)
            GOTO 9999

#ifdef NEVER
C           READ INPUT IMAGE.  NEEDS: BUF1(NXLD,NY)
            CALL READV(LUNIN,BUF1, NXLD,NY, NX,NY,1)

C           CENTER THE IMAGE NEEDS: BUF1(NXLD,NY), BUF2(NX,NY)
            CALL CENT_ROD(BUF1, BUF2, NXLD, NX,  NY,
     &                       THETA, OFFX,OFFY, IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999

C           WRITE OUTPUT IMAGE/VOLUME
            CALL WRITEV(LUNOUT,BUF2, NX,NY,NX,NY, NZ)
#endif

         END SELECT



         SELECT CASE(FCHAR(4:4))
         CASE ('P','S') !  --------------------------------- 'PH','SYM' 

            IF (NZ == 1) THEN    ! IMAGE
               NLIST = 4
               IF (IMGNUM > 0) THEN
                  IF (VERBOSE) WRITE(NOUT,90) IMGNUM,OFFX,OFFY
90                FORMAT('  IMAGE:',I6,
     &                   '  SHIFT:(',F9.3,',',F9.3,')')
               ELSE 
                  WRITE(NOUT,91) OFFX,OFFY 
91                FORMAT('  SHIFT:(',F9.3,',',F9.3,')')
               ENDIF
            ELSE   ! VOLUME
               NLIST = 5
               IF (IMGNUM > 0) THEN
                  IF (VERBOSE) WRITE(NOUT,93) IMGNUM,
     &                           OFFX,OFFY,OFFZ
93                FORMAT('  VOLUME:',I6,
     &                   '  SHIFT:(',F9.3,',',F9.3,',',F9.3')')
               ELSE 
                  WRITE(NOUT,94) OFFX,OFFY,OFFZ 
94                FORMAT('  SHIFT:(',F9.3,',',F9.3,',',F9.3,')')
               ENDIF
            ENDIF

C           FILL COMMAND LINE REGISTERS IF WANTED
            CALL REG_SET_NSEL(1,3,OFFX,OFFY,OFFZ,0.0,0.0,IRTFLG)

         CASE ('R') !  ----------------------------------------- 'ROD' 

            NLIST = 4
            IF (IMGNUM > 0) THEN
               IF (VERBOSE)  WRITE(NOUT,92)IMGNUM,THETA,OFFX,OFFY
92             FORMAT('  IMAGE:',I6,
     &                '  ANGLE:',F8.3,
     &                '  SHIFT:(',F9.3,',',F9.3,')')
            ELSE 
               !write(6,91) theta,offx,offy 
               WRITE(NOUT,95) THETA,OFFX,OFFY 
95             FORMAT('  ANGLE:',F8.3,
     &                '  SHIFT:(',F9.3,',',F9.3,')')

C              FILL COMMAND LINE REGISTERS IF WANTED
               CALL REG_SET_NSEL(1,3,THETA,OFFX,OFFY,0.0,0.0,IRTFLG)
            ENDIF
         END SELECT

         IF (WANTDOCOUT) THEN
            DLIST(1) = IMGNUM
            DLIST(2) = THETA
            DLIST(3) = OFFX
            DLIST(4) = OFFY
            DLIST(5) = OFFZ
            CALL LUNDOCWRTDAT(LUNDOC,IKEY,DLIST,NLIST,IRTFLG)
         ENDIF

C        IMAGE STATISTICS HAVE CHANGED
         CALL SETPRMB(LUNOUT, 0,0, 0,0)

         IF (NINDX1 >= NIMG) EXIT      ! END OF INPUT LIST
             
C        OPEN NEXT SET OF I/O FILES 
         !write(6,*) NINDX1, NINDX2,NIMG,NIMGOUT,MAXIMIN,MAXIMOUT
         !write(6,*) LUNIN,LUNIN,LUNOUT, FILPATIN,FILPATOUT,IMGNUM,IMGNUMOUT
         CALL NEXTFILES(NINDX1, NINDX2, ILIST1,ILIST2, 
     &                  .FALSE.,LUNDOC,LUNXM2,
     &                  NIMG,NIMGOUT,   
     &                  MAXIMIN,MAXIMOUT,   
     &                  LUNIN,LUNIN,LUNOUT, FILPATIN,FILPATOUT,
     &                  IMGNUM,IMGNUMOUT, IRTFLG) 

c          write(6,'(A,5i6)')
c     &       ' Nextfiles imgnum,imgnumout,nindx1,nindx2,irtflg:',
c     &                   imgnum,imgnumout,nindx1,nindx2,irtflg

         IF (IRTFLG == -99) THEN
             CALL ERRT(102,'INSUFFICIENT OUTPUT FILE NAMES',NINDX2)
             EXIT         
         ELSEIF (IRTFLG < 0) THEN
             EXIT         ! END OF INPUT FILES
         ENDIF
         IF (IRTFLG .NE. 0) GOTO 9999    ! ERROR
      ENDDO

9999  IF (ALLOCATED(BUF1))      DEALLOCATE(BUF1)
      IF (ALLOCATED(BUF2))      DEALLOCATE(BUF2)
      IF (ALLOCATED(ILIST1))    DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2))    DEALLOCATE(ILIST2)

      IF (VERBOSE) WRITE(NOUT,*) ' '

      CLOSE(LUNIN)
      CLOSE(LUNOUT)
      CLOSE(LUNDOCSEL)
      CLOSE(LUNXM1)
      CLOSE(LUNXM2)
      CLOSE(LUNDOC)
      
      END

