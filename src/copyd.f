
C++*********************************************************************
C
C COPYD.F      CREATED                        23 DEC 87 ARDEAN LEITH
C              USED GETOLDSTACK, GETNEWSTACK  APRIL 99  ARDEAN LEITH
C              GETNEWSTACK PARAM.             FEB 03    ARDEAN LEITH
C              FLIPEND                        FEB 03    ARDEAN LEITH
C              MPI                            FEB 04    Chao Yang
C              OPFILES                        DEC 10    ARDEAN LEITH
C              INDEXED STACK BUG              JAN 11    ARDEAN LEITH
C              INDEXED STACK BUG              FEB 11    ARDEAN LEITH
C              NON SPIDER IMAGE COPY          MAY 14    ARDEAN LEITH
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
C  COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPOUT)
C
C  PURPOSE:  COPY A SPIDER IMAGE FILE TO ANOTHER FILE
C 
C  PARAMETERS: LUN1,LUN2      READ & WRITE UNIT                 (SENT)
C              LUNDOC         READ & WRITE UNIT                 (SENT)
C              LUNXM1,LUNXM2  READ & WRITE UNIT                 (SENT)
C              INDXD          CREATE INDXD STACK                (SENT)
C              FLIPOUT        CHANGE BYTE ORDER                 (SENT)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPOUT)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOERR/  IERR
        COMMON /IOBUF/  BUF(NBUFSIZ)

        INTEGER                  :: LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2
        LOGICAL                  :: INDXD,FLIPOUT

        INTEGER                  :: IFLIPOUT
        CHARACTER (LEN=MAXNAM)   :: PROMPT 
        CHARACTER (LEN=MAXNAM)   :: FILNAM1,FILNAM2
        CHARACTER (LEN=2*MAXNAM) :: COMMAN
        LOGICAL                  :: VERBOSE_SAVE,IFLIP
        INTEGER,ALLOCATABLE      :: ILIST1(:),ILIST2(:)
        CHARACTER (LEN=1)        :: NULL = CHAR(0)
        CHARACTER (LEN=1)        :: DISP

        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

        VERBOSE_SAVE = VERBOSE           ! SAVE CURRENT VERBOSITY

        NILMAX       = NIMAX             ! FROM CMLIMIT
        ALLOCATE(ILIST1(NIMAX),
     &           ILIST2(NIMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'COPYD; ILIST....',2*NIMAX)
           RETURN
        ENDIF

C       OPEN FIRST INPUT FILE, DISP = 'E' DOES NOT STOP ON ERROR
        MAXIM1 = 0
        PROMPT = 'INPUT FILE NAME OR TEMPLATE (E.G. STK@****)~~9'
        CALL OPFILES(0,LUN1,LUNDOC,LUNXM1,  
     &               .TRUE.,FILNAM1,NLET1, 'E',
     &               IFORM1,NX1,NY1,NZ1,NSTACK1,
     &               PROMPT,
     &              .TRUE., ILIST1,NILMAX, 
     &               NDUM,NGOT1,IMG1, IRTFLG) 

        IF (IRTFLG == 5) THEN

           CALL FILERD(FILNAM2,NLET2,NULL,'OUTPUT FILE NAME',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 

           !write(6,*) ' Non spider input: ', filnam1(1:nlet1)
           !write(6,*) ' Non spider output: ',filnam2(1:nlet2)

           LOCDOT1 = ( INDEX(FILNAM1(1:NLET1),'.',BACK=.TRUE.) )
           LOCDOT2 = ( INDEX(FILNAM2(1:NLET2),'.',BACK=.TRUE.) )

           IF ( LOCDOT2 <= 0 ) THEN
C             APPEND EXTENSION TO FILNAM2
              FILNAM2 = FILNAM2(1:NLET2) // '.' // FILNAM1(LOCDOT1+1:)
              NLET2   = lnblnkn(FILNAM2)
           ENDIF

C          USE SYSTEM COPY
           COMMAN = 'cp ' // FILNAM1(1:NLET1) //' '// FILNAM2(1:NLET2)

           WRITE(NOUT,'(2X,2A,2X,A)'),
     &           'cp ',FILNAM1(:NLET1),FILNAM2(:NLET2)
           IRET = system(COMMAN)
           GOTO 999
        ENDIF

        IF (IRTFLG .NE. 0) RETURN

        !write(6,'(a,8i5)')' In  nstack1,ngot1,img1:',nstack1,ngot1,img1

        NSTACK2 =  1 ! UNUSED
        DISP    = 'U'          
        IF (NSTACK1 .GT. 0 .AND. INDXD) THEN
C          INPUT IS A WHOLE STACK AND WANT INDEXED OUTPUT STACK
           DISP    = 'I'
           NSTACK2 = NSTACK1  ! MAX SIZE
        ENDIF

        IF (FLIPOUT) THEN
C          STANDARD COPY WITH FLIPPED ENDEDNESS
           CALL LUNGETFLIP(LUN1,IFLIPIN,IRTFLG)
           IF (IFLIPIN .NE. 1) IFLIPOUT = 1
        ENDIF
        IF (IFLIPOUT == 1) CALL LUNSETFLIP(LUN2,IFLIPOUT,IRTFLG)

C	OPEN FIRST OUTPUT FILE
        IMG2 = IMG1
        CALL OPFILES(LUN1,LUN2,LUNDOC,LUNXM2, 
     &             .TRUE.,FILNAM2,NLET2,DISP,
     &             IFORM1,NX1,NY1,NZ1,NSTACK2,
     &             NULL,
     &             .TRUE., ILIST2,NILMAX, 
     &             NDUM,NGOT2,IMG2, IRTFLG) 

        !write(6,'(A,4i6)')' Out nstack2,ngot2,img2:',nstack2,ngot2,img2

        IFLIPOUT = 0
        IF (FLIPOUT) THEN
C          STANDARD COPY WITH FLIPPED ENDEDNESS
           CALL LUNGETFLIP(LUN1,IFLIPIN,IRTFLG)
           IF (IFLIPIN .NE. 1) IFLIPOUT = 1
        ENDIF
        !write(6,'(A,4i6)')' IFLIPIN, IFLIPOUT:',IFLIPIN, IFLIPOUT 

       IF (IFLIPOUT == 1) THEN
C          TELL WRITLIN TO FLIP CONTENTS DURING I/O
           CALL LUNSETFLIP(LUN2,IFLIPOUT,IRTFLG)
C          REPLACE HEADER WITH BYTE-FLIPPED HEADER
           CALL LUNWRTHED (LUN2,NX1,0,IRTFLG)
        ENDIF
       !call lungetstat(lun2,imamit,fmint,fmaxt,avt,sigt,irtflg)
       !write(6,*)' stats3:',imamit,fmint,fmaxt,avt,sigt 

        NINDX1 = 1
        NINDX2 = 1
        DO                ! LOOP OVER ALL IMAGES/STACKS

C          DO NOT REPORT FILE INFO IF WHOLE STACK
	   IF (NSTACK1 > 0 .AND. NSTACK2 >= 0) VERBOSE = .FALSE. 

C          COPY THE DESIRED NUMBER OF DATA RECORDS
           DO IREC = 1,NY1 * NZ1
              CALL REDLIN(LUN1,BUF,NX1,IREC)
              CALL WRTLIN(LUN2,BUF,NX1,IREC)
           ENDDO

C          OPEN NEXT SET OF I/O FILES, UPDATES NINDX* 
           !write(6,*)  ' CALLING nextfiles'
           CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                    .FALSE., LUNXM1,LUNXM2,
     &                    NGOT1,NGOT2,    NSTACK1,NSTACK2,  
     &                    LUN1,LUN1,LUN2, FILNAM1,FILNAM2,
     &                    IMG1,IMG2, IRTFLG)
 
           !write(6,*) ' After nextfiles, irtflg',irtflg
           !if (irtflg .ne. 0) then
           !write(6,'(A,4i6)') 
!     &        ' Nextfiles img1,img2,irtflg:',img1,img2,irtflg
           !write(6,'(A,4i6)') 
!     &        ' Nextfiles nindx1,nindx2:',nindx1,nindx2
           !endif

           IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
       ENDDO

       IRTFLG = 0
   
999    CLOSE(LUN1)
       CLOSE(LUN2)

       VERBOSE = VERBOSE_SAVE          ! RESTORE VERBOSITY 
       IF (ALLOCATED(ILIST1)) DEALLOCATE(ILIST1)
       IF (ALLOCATED(ILIST2)) DEALLOCATE(ILIST2)

       END

