
C++*********************************************************************
C
C  UTIL_1110.F  FROM ROTQSS                      MAR 2012 ARDEAN LEITH 
C               ZSH < 0 BUG                      APR 2014 ARDEAN LEITH
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
C   UTIL_1110()
C
C   PURPOSE:  SHIFT IMAGES.
C             CAN USE A SELECTION DOC FILE
C             CAN TAKE TRANSFORMATION PARAMETERS FROM A DOC FILE
C
C--*********************************************************************

      SUBROUTINE UTIL_1110()

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

C     DOC FILE POINTER
      INCLUDE 'F90ALLOC.INC'
      REAL, POINTER          :: PBUF(:,:)

      REAL,    ALLOCATABLE   :: BUF1(:),BUF2(:)
      INTEGER, ALLOCATABLE   :: ILIST1(:),ILIST2(:)

      CHARACTER (LEN=MAXNAM) :: FILPATIN,FILPATOUT,FILNAM
      CHARACTER (LEN=MAXNAM) :: DOCNAM,DOCNAME,MSG

      CHARACTER (LEN=1)      :: NULL = CHAR(0)
      CHARACTER (LEN=83)     :: PROMPT

      INTEGER                :: NILMAX,ITYPE,NX,NY,NZ
      INTEGER                :: NUMB,MAXX,MAXY,NLET,NDUM,NIMGOUT
      INTEGER                :: NIMG,IMGNUM,IMGNUMOUT,MWANT
      INTEGER                :: NINDX1,NINDX2,NC,NLETD
      INTEGER                :: MAXIMIN,MAXIMOUT
      INTEGER                :: LOCAT,LOCAST,NIMGT,NLETF
      INTEGER                :: NXLD,IDUM,NOT_USED,IREGX,IREGY,IREGZ
      INTEGER                :: ICOUNT, IRTFLG
      REAL                   :: XSH,YSH,ZSH
      LOGICAL                :: FOURIER

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
          CALL ERRT(46,'UTIL_1110; ILIST1....',2*NILMAX)
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
      IF (ITYPE .NE. 1 .AND. ITYPE .NE. 3) THEN
	 CALL ERRT(102,
     &     'OPERATION NOT IMPLEMENTED FOR DATA FORMAT',ITYPE)
         GOTO 9999
      ENDIF   ! END OF: IF (ITYPE....
   
C     SINGLE IMAGE OPERATION 
      IF (NIMG <= 1) IMGNUM = 0
            
C     OPEN OUTPUT IMAGE(S)
      PROMPT    = 'OUTPUT FILE NAME OR TEMPLATE (E.G. IMG@****)~'
      IF (LOCAST == 0)  PROMPT = 'OUTPUT'
      IF (LOCAST == 0 .AND. LOCAT > 0)  PROMPT = 'OUTPUT STACK'
 
      MAXIMOUT  = -1           ! ALLOW BARE STACK
      IMGNUMOUT = IMGNUM       ! IMAGE # WANTED
      NIMGT     = -NIMG        ! USE FIRST LIST FOR OUTPUT

      CALL OPFILES(LUNIN,LUNOUT,LUNDOCSEL,LUNXM2,
     &            .TRUE.,FILPATOUT,NLET, 'U',
     &            ITYPE,NX,NY,NZ,MAXIMOUT,
     &            PROMPT,
     &            .FALSE., ILIST2,NIMGT, 
     &            NDUM,NIMGOUT, IMGNUMOUT, IRTFLG) 

      IF (IRTFLG .NE. 0) GOTO 9999
      !write(6,*)'output Image #:',IMGNUMOUT,NIMG,maximout,filpatout(:11)
      !write(6,*) 'fchar:',fchar(1:10)

      FOURIER = (FCHAR(4:4) == 'F')

      NLETF = lnblnkn(FCHAR)
      SELECT CASE(FCHAR(1:NLETF))


      CASE ('SH','SH F') !--------------------------------------- 'SH' 
C        SH       SHIFT, SH F SHIFT USING FOURIER INTERP

         XSH  = 0
         YSH  = 0
         ZSH = -999999

         IF (NIMG <= 1) THEN
C           SINGLE IMAGE OPERATION, READ IN THE SHIFTS...FROM  USER

            IF (IFORM == 3) THEN
C              SHIFT VOLUME  ---------------------------------------- 3D

               CALL RDPRM3S(XSH,YSH,ZSH, NOT_USED,
     &                 'SHIFT COMPONENTS IN X, Y, & Z',IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9999
               IF (ZSH == -999999) THEN
C                 LEGACY INPUT
                  CALL RDPRM1S(ZSH,NOT_USED,
     &                         'SHIFT COMPONENT IN Z',IRTFLG)
                  IF (IRTFLG .NE. 0) GOTO 9999
               ENDIF

               NC   = 52
C                        123456789 123456789 123456789 123456789 12345 
               PROMPT = "REG. #'S FOR X, Y & Z SHIFTS " //
     &                  "(OR * FOR: 2,3,4)"
            ELSE
C              IMAGE  ----------------------------------------------- 2D
               CALL RDPRM2S(XSH,YSH,NOT_USED,
     &                 'SHIFT COMPONENTS IN X & Y',IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9999
               ZSH = 0

               NC  = 48
C                        123456789 123456789 123456789 123456789 12345 
               PROMPT = "REG. #'S FOR X & Y SHIFTS " //
     &                  "(OR * FOR: 2,3)"
            ENDIF   ! END OF: IF (IFORM == 3) THEN

         ELSE
C           HAVE MULTIPLE IMAGES, READ SHIFTS... FROM A DOC. FILE
         
            CALL FILERD(DOCNAM,NLETD,NULL,
     &		'SHIFT DOCUMENT',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           CALL FILNAMANDEXT(DOCNAM,DATEXC,DOCNAME,NLETD,
     &                       .TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (IFORM == 3) THEN
C              SHIFT VOLUME  ---------------------------------------- 3D

               NC   = 52
C                        123456789 123456789 123456789 123456789 12345 
               PROMPT = "REG. #'S FOR X, Y & Z SHIFTS " //
     &                  "(OR * FOR: 2,3,4)"

               IREGX = 2
               IREGY = 3
               IREGZ = 4
               CALL RDPRI3S(IREGX,IREGY,IREGZ,
     &                      NOT_USED,PROMPT(1:NC),IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9999

            ELSE
C              IMAGE  ----------------------------------------------- 2D
               NC   = 47
C                        123456789 123456789 123456789 123456789 12345 
               PROMPT = "REG. #'S FOR X & Y SHIFTS " //
     &                  "(OR * FOR: 2,3)"
               IREGX = 2
               IREGY = 3
               CALL RDPRI2S(IREGX,IREGY,NOT_USED,
     &                      PROMPT(1:NC),IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9999

            ENDIF   ! END OF: IF (IFORM == 3) THEN

C           MAXX IS 1 + NUM OF REGISTERS SINCE PBUF CONTAINS COUNT ALSO
            MAXX  = MAX(IREGX,IREGY) + 1
            MAXY  = 0
            CALL GETDOCDAT('',.FALSE.,DOCNAME,LUNDOC,
     &                      .TRUE.,MAXX, MAXY,PBUF,IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
         ENDIF ! END OF: IF (NIMG <= 1)

         NXLD  = NX+2-MOD(NX,2)
         MWANT = NX
         IF (FOURIER) MWANT = NXLD * NY * NZ

C        BUF2 ALWAYS EVEN IF FOUIER IN CASE INTEGER SHIFT
         ALLOCATE(BUF1(MWANT), BUF2(NX*6),STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            MWANT = NXLD*NY*NZ + NX*6
            CALL ERRT(46,'UTIL8, Q',MWANT)
            GOTO 9999
         ENDIF

      CASE DEFAULT
          CALL ERRT(101,'UNKNOWN CASE,PGM ERROR',IDUM)
      END SELECT


      NINDX1 = 1
      NINDX2 = 1
      NLETF  = lnblnkn(FCHAR)
      DO 
         !write(6,*)  'Image #:',IMGNUM

         SELECT CASE(FCHAR(1:NLETF))

         CASE ('SH','SH F') !------------------------------------ 'SH' 
C           SH= SHIFT, SH F= SHIFT USING FOURIER INTERP

            IF (NIMG > 1) THEN
               ICOUNT = PBUF(1,IMGNUM)
c                                 0      9    21   1              0
               !write(6,*)' ICOUNT: ',icount,maxx,maxy,imgnum,pbuf(:,imgnum)
               IF (ICOUNT .LT. (MAXX-1)) THEN
                  !write(6,*)' ICOUNT: ',icount,maxx,maxy,imgnum,pbuf(:,imgnum)
                  !write(nout,*)' ICOUNT: ',icount,maxx,maxy,imgnum,pbuf(:,imgnum)
                  CALL ERRT(102,
     &              'LACK ROTATE/SHIFT PARAMETERS FOR IMAGE',IMGNUM)
                  GOTO 9999
               ENDIF

               XSH  = PBUF(IREGX+1, IMGNUM)
               YSH  = PBUF(IREGY+1, IMGNUM)
            ENDIF

            CALL SHIFTS(FOURIER, LUNIN,LUNOUT, BUF1,BUF2,
     &                  NX,NY,NZ, XSH,YSH,ZSH, IRTFLG)

            IF (VERBOSE .AND. IMGNUM > 0) THEN
               IF (NZ == 1) THEN
                  WRITE(NOUT,90)IMGNUM,XSH,YSH 
90                FORMAT('  IMAGE:',I6,
     &                   '  SHIFT:(',F9.3,',',F9.3,')')
               ELSE
                  WRITE(NOUT,91)IMGNUM,XSH,YSH,ZSH 
91                FORMAT('  IMAGE:',I6,
     &                   '  SHIFT:(',F9.3,',',F9.3,',',F9.3,')')
               ENDIF
            ENDIF

C           IMAGE STATISTICS HAVE CHANGED
            CALL SETPRMB(LUNOUT, 0,0, 0,0)
         END SELECT

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

          !write(6,'(A,5i6)')
      !&       ' Nextfiles imgnum,imgnumout,nindx1,nindx2,irtflg:',
      !&                   imgnum,imgnumout,nindx1,nindx2,irtflg

         IF (IRTFLG == -99) THEN
             CALL ERRT(102,'INSUFFICIENT OUTPUT FILE NAMES',NINDX2)
             EXIT         
         ELSEIF (IRTFLG < 0) THEN
             EXIT         ! END OF INPUT FILES
         ENDIF
         IF (IRTFLG .NE. 0) GOTO 9999    ! ERROR
      ENDDO

9999  IF (ALLOCATED(BUF1))    DEALLOCATE(BUF1)
      IF (ALLOCATED(BUF2))    DEALLOCATE(BUF2)
      IF (ALLOCATED(ILIST1))  DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2))  DEALLOCATE(ILIST2)

      IF (VERBOSE) WRITE(NOUT,*) ' '

      CLOSE(LUNIN)
      CLOSE(LUNOUT)
      CLOSE(LUNXM1)
      CLOSE(LUNXM2)
      
      END









