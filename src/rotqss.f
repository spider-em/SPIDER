
C++*********************************************************************
C
C  ROTQSS.F   ADDED STACK OPERATION                    98 ARDEAN LEITH 
C             USED AUTO. ARRAYS                   5/02/00 ARDEAN LEITH
C             SELECTION FILE ADDED               12/15/06 ARDEAN LEITH
C             MULTIFILE SUPPORT ADDED            12/17/10 ARDEAN LEITH
C             MERGED WITH ROTQSS_DL               1/10/11 ARDEAN LEITH
C             ADDED RTKSQ                         5/20/11 ARDEAN LEITH
C             RENAMED RTKSQ --> RTFS              6/02/11 ARDEAN LEITH
C             RENAMED ROT2QS --> RTSQ            12/28/11 ARDEAN LEITH
C             NSAM --> NX, RTSQ CALL PARAM        1/04/12 ARDEAN LEITH
C             INSUFFICIENT OUTPUT FILE NAMES      1/18/12 ARDEAN LEITH
C             CHKMIRROR                           3/08/12 ARDEAN LEITH
C             NANG == 1 .AND. ILIST1(1) BUG       4/05/12 ARDEAN LEITH
C             MIRROR BUG                          5/16/12 ARDEAN LEITH
C             MIRROR = (FMIRROR <  0.0)           6/11/12 ARDEAN LEITH
C             NULLIFY(PBUF)                      12/07/12 ARDEAN LEITH
C             
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C   ROTQSS(LUNIN,LUNOUT,LUNDOC,LUNXM1,LUNXM2,OLDRTSQ,USERTSF, IRTFLG)
C
C   PURPOSE:  ROTATE, SCALE, AND SHIFT IMAGES.
C             CAN USE A SELECTION DOC FILE
C             CAN TAKE TRANSFORMATION PARAMETERS FROM A DOC FILE
C
C--*********************************************************************

      SUBROUTINE ROTQSS(LUNIN,LUNOUT,LUNDOC,LUNXM1,LUNXM2, 
     &                  OLDRTSQ,USERTSF,CHKMIRROR,IRTFLG)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER                :: LUNIN,LUNOUT,LUNDOC,LUNXM1,LUNXM2
      LOGICAL                :: OLDRTSQ, USERTSF,CHKMIRROR
      INTEGER                :: IRTFLG

C     DOC FILE POINTER
      INCLUDE 'F90ALLOC.INC'
      REAL, POINTER          :: PBUF(:,:)

      REAL,    ALLOCATABLE   :: XIMG(:,:),BUF2(:,:),BUFM(:,:)
      INTEGER, ALLOCATABLE   :: ILIST1(:),ILIST2(:)

      INTEGER                :: ILIST(5)
      CHARACTER (LEN=MAXNAM) :: FILPATIN,FILPATOUT,FILNAM
      CHARACTER (LEN=MAXNAM) :: DOCNAM
      CHARACTER (LEN=1)      :: ANS

      CHARACTER (LEN=1)      :: NULL = CHAR(0)
      CHARACTER (LEN=83)     :: PROMPT

      REAL                   :: THETA,SCLI,SHXI,SHYI,SCALE,FMIRROR
      INTEGER                :: NILMAX,MAXIMOUT,ITYPE,NX,NY,NZ
      INTEGER                :: NUMB,MAXX,MAXY,NLET,NDUM,NANGOUT
      INTEGER                :: NANG,IMGNUM,IMGNUMOUT,MWANT
      INTEGER                :: NINDX1,NINDX2,ICOUNT,ISLICE
      INTEGER                :: MAXIMIN,NOT_USED
      INTEGER                :: LOCAT,LOCAST,NANGT
      INTEGER                :: NXLD,NC,I,J,NXP1
      LOGICAL                :: ISBAREIN,ISBAREOUT,MIRROR,EVEN

      NULLIFY(PBUF)

      NILMAX  = NIMAXPLUS      ! FROM CMLIMIT
      ALLOCATE(ILIST1(NILMAX),
     &         ILIST2(NILMAX),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
          CALL ERRT(46,'ROTQSS; ILIST1....',2*NILMAX)
          RETURN
      ENDIF

C     OPEN INPUT IMAGE(S)
      CALL OPFILES(0,LUNIN,LUNDOC,LUNXM1,  .TRUE.,FILPATIN,NLET, 'O',
     &              ITYPE,NX,NY,NZ,MAXIMIN,
     &              NULL,
     &              .FALSE.,ILIST1,NILMAX, 
     &              NDUM,NANG,IMGNUM, IRTFLG) 
      IF (IRTFLG .NE. 0) RETURN
      LOCAT  = INDEX(FILPATIN,'@')
      LOCAST = INDEX(FILPATIN,'*')
      IF (NANG > 0) ILIST2 = ILIST1

      !write(6,*)' Input image #: ',IMGNUM,nang,maximin,locast
      !writewrite(6,*)' Input image #: ',IMGNUM,nang,ilist1(1),maximin,locast
  
      IF (OLDRTSQ .AND. LOCAST <= 0) THEN
C        OLD 'RT SQ' & NO '*' PROMPT ORDER. OPEN OUTPUT IMAGE(S) NOW
         MAXIMOUT = -1            ! ALLOW BARE STACK
         IMGNUMOUT = 1            ! IMAGE # WANTED
         NANGT     = NANG
         !write(6,*)  'old input Image #: ',IMGNUM,nang,maximin,locast
         CALL OPFILES(LUNIN,LUNOUT,LUNDOC,LUNXM2,
     &            .TRUE.,FILPATOUT,NLET, 'U',
     &            ITYPE,NX,NY,NZ,MAXIMOUT,
     &            'OUTPUT FILE NAME OR TEMPLATE (E.G. ROT@****)~',
     &            .FALSE., ILIST2,NANGT, 
     &            NDUM,NANGOUT, IMGNUMOUT, IRTFLG) 
         IF (IRTFLG .NE. 0) GOTO 9999
      ENDIF

      IF (NANG <= 0  .OR. 
     &    NANG == 1 .AND. ILIST1(1) == 0) THEN

C        SINGLE IMAGE OPERATION, READ IN THE ANGLES...
         THETA  = 0.0
         SCLI   = 1.0
         MIRROR = .FALSE.
         IF (CHKMIRROR) THEN
            FMIRROR = 0.0
            CALL RDPRM3S(THETA,SCLI,FMIRROR,NOT_USED,
     &          'ROTATION ANGLE, SCALE FACTOR, MIRROR (IF < 0)',
     &          IRTFLG)
            MIRROR = (FMIRROR <  0.0)
         ELSE
            CALL RDPRM2S(THETA,SCLI,NOT_USED,
     &                'ROTATION ANGLE, SCALE FACTOR',IRTFLG)
         ENDIF
	 IF (SCLI <= 0.0) SCLI = 1.0

         SHXI = 0.0
         SHYI = 0.0
	 CALL RDPRM2S(SHXI,SHYI,NOT_USED,
     &               'SHIFTS IN X AND Y',IRTFLG)
         IMGNUM = 0

      ELSE
C        HAVE MULTIPLE ANGLES, READ ANGLES... FROM A DOC. FILE
         IF (CHKMIRROR) THEN
             NUMB = 5
             NC   = 79
C                      123456789 123456789 123456789 123456789 12345 
             PROMPT = "ENTER REG. #'S FOR ANGLE, SCALE, X & Y " //
     &                "SHIFT, AND MIRROR (OR * FOR: 6,0,7,8,15)"
         ELSE
             NUMB = 4
             NC   = 64
C                      123456789 123456789 123456789 123456789 12345
             PROMPT = "ENTER REG. #'S FOR ANGLE, SCALE, X & Y " //
     &                "SHIFT (OR * FOR: 6,0,7,8)"
         ENDIF
         CALL RDPRAI(ILIST,5,NUMB,1,15,PROMPT(1:NC),'T',IRTFLG)
         IF (IRTFLG == -1) THEN
            ILIST(1) = 6  ! Usual  Reg. #s for angle, scale, & shift
            ILIST(2) = 0
            ILIST(3) = 7
            ILIST(4) = 8
            ILIST(5) = 15
            NUMB     = 5
            IRTFLG   = 0
         ENDIF

         IF (IRTFLG .NE. 0) GOTO 9999

C        MAXX IS 1 + NUM OF REGISTERS SINCE PBUF CONTAINS COUNT ALSO
         MAXX  = MAXVAL(ILIST) + 1
         MAXY  = 0
         CALL GETDOCDAT('ANGLE/SHIFT DOCUMENT',.TRUE.,DOCNAM,LUNDOC,
     &             .TRUE.,MAXX, MAXY,PBUF,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
      ENDIF
      
       
      IF (.NOT. OLDRTSQ .OR. LOCAST > 0 .OR. 
     &   (.NOT. OLDRTSQ .AND. USERTSF)) THEN
C        'RT SQ'   WITH SET OF STACKED IMAGES, OPEN OUTPUT IMAGE(S) NOW
C        'RTD SQ'  OPEN OUTPUT IMAGE(S)
C        'RT  SF'  OPEN OUTPUT IMAGE(S)
         MAXIMOUT = -1            ! ALLOW BARE STACK
         IMGNUMOUT = IMGNUM       ! IMAGE # WANTED
         NANGT     = NANG

         IF (OLDRTSQ .AND. LOCAST > 0) THEN
C           'RT SQ' WITH IMG@***, USE ILIST2, KLUDGE FOR BACKWARD COMPAT.
             NANGT = -NANGT
         ENDIF

         !write(6,*)'output Image #:',IMGNUMOUT,nang,maximout,filpatout(:11)

         CALL OPFILES(LUNIN,LUNOUT,LUNDOC,LUNXM2,
     &            .TRUE.,FILPATOUT,NLET, 'U',
     &            ITYPE,NX,NY,NZ,MAXIMOUT,
     &            'OUTPUT FILE NAME OR TEMPLATE (E.G. ROT@****)~',
     &            .FALSE., ILIST2,NANGT, 
     &            NDUM,NANGOUT, IMGNUMOUT, IRTFLG) 

         IF (IRTFLG .NE. 0) GOTO 9999
      ENDIF
 
      EVEN = (MOD(NX,2) == 0)
      NXLD = NX 
      IF (USERTSF) NXLD = NX + 2 - MOD(NX,2)

      ALLOCATE (XIMG(NXLD,NY),
     &          BUF2(NX,NY), 
     &          STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN 
         MWANT = NXLD*NY + NX*NY 
         CALL ERRT(46,'ROTQSS; XIMG..',MWANT)
         GOTO 9999
      ENDIF  

      IF (CHKMIRROR) THEN
         ALLOCATE (BUFM(NX,NY), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            CALL ERRT(46,'ROTQSS; BUFM..',NX*NY)
            GOTO 9999
         ENDIF  
      ENDIF

      NINDX1 = 1
      NINDX2 = 1
      DO 
         !write(6,*)  'Image #,nang:',IMGNUM,NANG

         IF (NANG > 1 .OR. (NANG == 1 .AND. ILIST1(1) > 0) ) THEN
C           MUST READ ANGLES FROM DOC FILE
            ICOUNT = PBUF(1,IMGNUM)
c                                 0      9    21   1              0
            !write(6,*)' ICOUNT: ',icount,maxx,maxy,imgnum,pbuf(:,imgnum)
            IF (ICOUNT .LT. (MAXX-1)) THEN
            write(6,*)' ICOUNT: ',icount,maxx,maxy,imgnum,pbuf(:,imgnum)
         write(nout,*)' ICOUNT: ',icount,maxx,maxy,imgnum,pbuf(:,imgnum)
               CALL ERRT(102,
     &           'LACK ROTATE/SHIFT PARAMETERS FOR IMAGE',IMGNUM)
               GOTO 9999
            ENDIF

            THETA = 0.0 
            IF (ILIST(1) > 0) THEN
               THETA = PBUF(ILIST(1)+1, IMGNUM)
            ENDIF

            SCLI  = 1.0
            IF (ILIST(2) > 0) THEN
               SCLI  = PBUF(ILIST(2)+1, IMGNUM)
	       IF (SCLI .EQ. 0.0) SCLI = 1.0
            ENDIF

            SHXI  = 0.0
            IF (ILIST(3) > 0) THEN
               SHXI  = PBUF(ILIST(3)+1, IMGNUM)
            ENDIF

            SHYI  = 0.0
            IF (ILIST(4) > 0) THEN
               SHYI  = PBUF(ILIST(4)+1, IMGNUM)
            ENDIF

            IF (CHKMIRROR) MIRROR = (PBUF(ILIST(5)+1, IMGNUM) < 0) 
         ENDIF

C        LOAD AND ROTATE VOLUME, SLICE BY SLICE
         DO ISLICE=1,NZ

C           READ INPUT IMAGE
            CALL READV(LUNIN,XIMG, NXLD,NY, NX,NY,ISLICE)

            IF (USERTSF) THEN
               !write(6,*)' Calling rtsf:',theta,scli,shxi,shyi
               CALL RTSF(XIMG, BUF2, 
     &                   NXLD, NX, NY,
     &                   THETA,SCLI,SHXI,SHYI,  IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 9999
            ELSE
               !write(6,*)' Calling rtsq:',theta,scli,shxi,shyi
    	       CALL RTSQ(XIMG, BUF2, 
     &                   NX,NY, NX,NY,
     &                   THETA,SCLI,SHXI,SHYI, IRTFLG)
            ENDIF

C           WRITE OUTPUT IMAGE
            IF (.NOT. CHKMIRROR .OR. .NOT. MIRROR) THEN
               CALL WRITEV(LUNOUT,BUF2, NX,NY,NX,NY, ISLICE)
            ELSE
               NXP1 = NX+1
               DO J = 1,NY
                  IF (EVEN) THEN
C                    EVEN LINE LENGTH
                     DO I = 1,NX
 	                BUFM(MOD(NXP1-I,NX)+1,J) = BUF2(I,J) ! 1-->1,60-->2
                     ENDDO
                  ELSE
C                    ODD LINE LENGTH
                     DO I = 1,NX
                        BUFM(I,J) = BUF2(NXP1-I,J)  ! 1-->60,2-->59
                     ENDDO
                  ENDIF
               ENDDO
              
               CALL WRITEV(LUNOUT,BUFM, NX,NY,NX,NY, ISLICE)
            ENDIF

         ENDDO
         CALL SETPRMB(LUNOUT, 0,0, 0,0)

        IF (VERBOSE .AND. IMGNUM > 0) THEN
            IF (.NOT. CHKMIRROR) THEN 
               WRITE(NOUT,90)IMGNUM,THETA,SCLI,SHXI,SHYI 
          
90             FORMAT('  IMAGE:',I6,
     &                '  ANGLE:',F8.3,
     &                '  SCALE:',F7.3,
     &                '  SHIFT:(',F9.3,',',F9.3,')')
            ELSE          

               WRITE(NOUT,91)IMGNUM,THETA,SCLI,SHXI,SHYI,MIRROR 
91             FORMAT('  IMAGE:',I6,
     &                '  ANGLE:',F8.3,
     &                '  SCALE:',F7.3,
     &                '  SHIFT:(',F9.3,',',F9.3,')',
     &                '  MIRROR:  ',L)
            ENDIF
 
        ELSEIF (VERBOSE) THEN
            IF (.NOT. CHKMIRROR) THEN 
               WRITE(NOUT,92) THETA,SCLI,SHXI,SHYI 
          
92             FORMAT('  ANGLE:',F8.3,
     &                '  SCALE:',F7.3,
     &                '  SHIFT:(',F9.3,',',F9.3,')')
            ELSE          
               WRITE(NOUT,93)THETA,SCLI,SHXI,SHYI,MIRROR 
93             FORMAT('  ANGLE:',F8.3,
     &                '  SCALE:',F7.3,
     &                '  SHIFT:(',F9.3,',',F9.3,')',
     &                '  MIRROR:  ',L)
            ENDIF
         ENDIF

         IF (NINDX1 >= NANG) EXIT      ! END OF INPUT LIST
             
C        OPEN NEXT SET OF I/O FILES 
         CALL NEXTFILES(NINDX1, NINDX2, ILIST1,ILIST2, 
     &                  .FALSE.,LUNDOC,LUNXM2,
     &                  NANG,NANGOUT,   
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

9999  IF (ALLOCATED(XIMG))      DEALLOCATE(XIMG)
      IF (ALLOCATED(BUF2))      DEALLOCATE(BUF2)
      IF (ALLOCATED(BUFM))      DEALLOCATE(BUFM)
      IF (ALLOCATED(ILIST1))    DEALLOCATE(ILIST1)
      IF (ALLOCATED(ILIST2))    DEALLOCATE(ILIST2)

C     DEALLOCATE DOC. FILE MEMORY
      IF (NANG > 1 .AND. ASSOCIATED(PBUF)) THEN
         DEALLOCATE(PBUF)
         NULLIFY(PBUF)
      ENDIF
      IF (VERBOSE) WRITE(NOUT,*) ' '

      CLOSE(LUNIN)
      CLOSE(LUNOUT)
      CLOSE(LUNXM1)
      CLOSE(LUNXM2)
      
      END

