
C++*********************************************************************
C
C COPYD.F      CREATED                        23 DEC 87 ArDean Leith
C              USED GETOLDSTACK, GETNEWSTACK  APRIL 99  ArDean Leith
C              GETNEWSTACK PARAM.             FEB 03    ArDean Leith
C              FLIPEND                        FEB 03    ArDean Leith
C              MPI                            FEB 04    Chao Yang
C              OPFILES                        DEC 10    ArDean Leith
C              INDEXED STACK BUG              JAN 11    ArDean Leith
C              INDEXED STACK BUG              FEB 11    ArDean Leith
C              NON SPIDER IMAGE COPY          MAY 14    ArDean Leith
C              DEBUG OUTPUT ADDED             SEP 20    ArDean Leith
C              FLIP OUTPUT ADDED              MAR 26    ArDean Leith
C              DEBUG OUTPUT ADDED             APR 26    ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2026  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email:                                                             *
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
C  COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPEND)
C
C  PURPOSE:  COPY SPIDER IMAGE/VOL TO ANOTHER FILE INCLUDING MRC FIL
C 
C  PARAMETERS: LUN1,LUN2      READ & WRITE UNIT                 (SENT)
C              LUNDOC         READ & WRITE UNIT                 (SENT)
C              LUNXM1,LUNXM2  READ & WRITE UNIT                 (SENT)
C              INDXD          CREATE INDEXED STACK              (SENT)
C              FLIPEND        CHANGE BYTE ORDER                 (SENT)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE COPYD(LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2,INDXD,FLIPEND)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL            :: BUF
        INTEGER         :: IERR
        INTEGER         :: LUNARA,LUNSTK,LUNARB,LUNFLIP

        COMMON /IOERR/  IERR
        COMMON /IOBUF/  BUF(NBUFSIZ)

        common /lunara/ lunara(100),lunstk(100),lunarb(100),lunflip(100)

        INTEGER                  :: LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2
        LOGICAL                  :: INDXD,FLIPEND

        INTEGER                  :: IFLIPIN, IFLIPOUT
        CHARACTER (LEN=MAXNAM)   :: PROMPT 
        CHARACTER (LEN=MAXNAM)   :: FILNAM1,FILNAM2
        CHARACTER (LEN=2*MAXNAM) :: COMMAN
        LOGICAL                  :: VERBOSE_SAVE,IFLIP
        INTEGER,ALLOCATABLE      :: ILIST1(:),ILIST2(:)
        CHARACTER (LEN=1)        :: NULL = CHAR(0)
        CHARACTER (LEN=1)        :: DISP
        INTEGER                  :: NE,LOCDOT1,LOCDOT2,NLET1,NLET2,IRET
        INTEGER                  :: NINDX1,NINDX2,NLIST1,NLIST2,NSTACK1
        INTEGER                  :: NSTACK2,IMG1,IMG2,NDUM,NX1,NY1,NZ1
        INTEGER                  :: MAXIM1,MAXIM2,IRTFLG,NILMAX,MPIERR
        INTEGER                  :: IFORM1,IFORM2,ICOMM,MYPID

        integer                  :: lnblnkn,system

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
        PROMPT = 'INPUT FILE NAME OR TEMPLATE (E.G. STK@*)~~9'
        CALL OPFILES(0,LUN1,LUNDOC,LUNXM1,  
     &               .TRUE.,FILNAM1,NLET1, 'E',
     &               IFORM1,NX1,NY1,NZ1,NSTACK1,
     &               PROMPT,
     &              .TRUE., ILIST1,NILMAX, 
     &               NDUM,NLIST1,IMG1, IRTFLG) 

        IF (IRTFLG == 5) THEN
C          NOT A VALID SPIDER OR MRC INPUT DATA FILE ---------------

           CALL FILERD(FILNAM2,NLET2,NULL,'OUTPUT FILE NAME',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN 

           !write(3,*) ' In copyd: Non spider input: ', filnam1(1:nlet1)
           !write(3,*) ' In copyd: Non spider output: ',filnam2(1:nlet2)

           LOCDOT1 = ( INDEX(FILNAM1(1:NLET1),'.',BACK=.TRUE.) )
           LOCDOT2 = ( INDEX(FILNAM2(1:NLET2),'.',BACK=.TRUE.) )

           IF ( LOCDOT2 <= 0 ) THEN 
C             APPEND EXTENSION TO FILNAM2
              FILNAM2 = FILNAM2(1:NLET2) // '.' // FILNAM1(LOCDOT1+1:)
              NLET2   = lnblnkn(FILNAM2)
           ENDIF

C          USE SYSTEM COPY
           COMMAN = 'cp ' // FILNAM1(1:NLET1) //' '// FILNAM2(1:NLET2)

           WRITE(NOUT,'(2X,2A,2X,A)')
     &           'cp ',FILNAM1(:NLET1),FILNAM2(:NLET2)
           IRET = system(COMMAN)
           GOTO 999
	   
	ENDIF      ! END OF:  IF (IRTFLG == 5) -------------------
	
        IF (IRTFLG .NE. 0) GOTO 999   ! ABORT OPERATION ----------

        NSTACK2 =  1 ! UNUSED
        DISP    = 'U'          
        IF (NSTACK1 > 0 .AND. INDXD) THEN
C          INPUT IS A WHOLE STACK AND WANT INDEXED OUTPUT STACK
           DISP    = 'I'
           NSTACK2 = NSTACK1  ! MAX IMAGE IN STACK
        ENDIF

C       FIND CURRENT ENDIAN FLIP OF SPIDER INPUT     
        CALL LUNGETFLIP(LUN1,IFLIPIN,IRTFLG)


C       FIND ENDIANESS FLIP OF INPUT
     
        IF (FLIPEND) THEN
           IF ((NSTACK2 > -2) .AND.  (IFLIPOUT .NE. IFLIPIN)) THEN
C             WANT TO FLIP ENDIANESS BUT OUTPUT IS A STACKED IMAGE 
C             OR SERIES OF IMAGES. FLIP MUST MATCH CURRENT STACK FLIP
              CALL ERRT(101,
     &           'EXISTING OUTPUT STACK NEEDS SAME ENDIANESS',NE)
              RETURN
           ELSE

C             TELL WRTLIN TO FLIP FILE CONTENTS OUTPUT
C             FLIP ENDIANESS FOR THE WHOLE FILE

             IFLIPOUT = 1  
           ENDIF

        ELSE
C          NORMAL COPY. DO NOT FLIP ENDIANESS
           IFLIPOUT = 0
        ENDIF

        CALL LUNSETFLIP(LUN2,IFLIPOUT,IRTFLG)

c    !!!!!!!!!!!!!!!!!!!!!!!!! REMOVE
cc#undef SP_DBUGIO
cc#define SP_DBUGIO
#if defined (SP_DBUGIO)
        write(3,*)' ' 
        write(3,*)' In copyd 0; lun2:          ', lun2 
        write(3,*)' In copyd 0; lunflip(lun2): ', lunflip(lun2)
        write(3,*)' In copyd 0; iflipout:      ', iflipout 
        write(3,*)' ' 
#endif

C	OPEN FIRST OUTPUT FILE.  OUTPUT STACK FILE MAY EXIST!!
        IMG2   = IMG1
        NLIST2 = 0
        CALL OPFILES(LUN1,LUN2,LUNDOC,LUNXM2, 
     &             .TRUE.,FILNAM2,NLET2,DISP,
     &             IFORM1,NX1,NY1,NZ1,NSTACK2,
     &             NULL,
     &             .TRUE., ILIST2,NILMAX, 
     &             NDUM,NLIST2,IMG2, IRTFLG) 
	
#if defined (SP_DBUGIO)
        write(3,*)' ' 
        write(3,*)' In copyd 1; flipend:       ', flipend 

        write(3,*)' In copyd 1; lun1:          ', lun1
        write(3,*)' In copyd 1; lunflip(lun1): ', lunflip(lun1)
        write(3,*)' In copyd 1; iflipin:       ', iflipin
        write(3,*)' In copyd 1; nstack1:       ', nstack1
        write(3,*)' In copyd 1; nlist1:        ', nlist1
        write(3,*)' In copyd 1; img1:          ', img1
        write(3,*)' ' 
        write(3,*)' In copyd 1; lun2:          ', lun2
        write(3,*)' In copyd 1; lunflip(lun2): ', lunflip(lun2)
        write(3,*)' In copyd 1; iflipout:      ', iflipout 
        write(3,*)' In copyd 1; nstack2:       ', nstack2
        write(3,*)' In copyd 1; nlist2:        ', nlist2
        write(3,*)' In copyd 1; img2:          ', img2
        write(3,*)' In copyd 1; irtflg         ', irtflg
        write(3,*)' -------------- COPYING NOW -----------------------'
        write(3,*)' '
#endif
	
        NINDX1 = 1
        NINDX2 = 1

        DO       ! LOOP OVER ALL IMAGES/STACKS

C          DO NOT REPORT FILE INFO IF WHOLE STACK (VERBOSE IN COMMON)
	   !!!!!!!IF (NSTACK1 > 0 .AND. NSTACK2 >= 0) VERBOSE = .FALSE. 

C          COPY THE DESIRED NUMBER OF DATA RECORDS
           DO IREC = 1,NY1 * NZ1
              CALL REDLIN(LUN1,BUF,NX1,IREC)
              CALL WRTLIN(LUN2,BUF,NX1,IREC)
           ENDDO


C          OPEN NEXT SET OF I/O FILES, UPDATES NINDX1 & NINDX2 
#if defined (SP_DBUGIO)  
           write(3,*)'   ' 
           write(3,*)' In copyd 6; nindx1,nindx2: ',nindx1,nindx2 
           write(3,*)' In copyd 6; nlist1,nlist2: ',nlist1,nlist2 
           write(3,*)' In copyd 6; img1,img2:     ',img1,img2 
           write(3,*)' In copyd 6; CALLING NEXTFILES NOW NNNNNNNNNNNNNN'
#endif

           CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                    .FALSE., LUNXM1,LUNXM2,
     &                    NLIST1,NLIST2,    NSTACK1,NSTACK2,  
     &                    LUN1,LUN1,LUN2, FILNAM1,FILNAM2,
     &                    IMG1,IMG2, IRTFLG)

#if defined (SP_DBUGIO)
           write(3,*)'   ' 
           write(3,*)' In copyd 7; irtflg   ', irtflg
           write(3,*)' In copyd 7; img1:    ', img1
           write(3,*)' In copyd 7; img2:    ', img2
           write(3,*)' In copyd 7; nlist1:  ', nlist1
           write(3,*)' In copyd 7; nlist2:  ', nlist2
           write(3,*)' In copyd 7; nstack1: ', nstack1
           write(3,*)' In copyd 7; nstack2: ', nstack2
           write(3,*)' In copyd 7; LOOPING NOW -----------------------'
           write(3,*)'   '
#endif

           IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
       ENDDO

       IRTFLG = 0
   
999    CLOSE(LUN1)
       CLOSE(LUN2)

       VERBOSE = VERBOSE_SAVE          ! RESTORE VERBOSITY 
       IF (ALLOCATED(ILIST1)) DEALLOCATE(ILIST1)
       IF (ALLOCATED(ILIST2)) DEALLOCATE(ILIST2)

       END

