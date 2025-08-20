
C **********************************************************************
C
C COPYFROMMRC   MODIFIED FROM COPYMRC             FEB 02 ArDean Leith         
C               ISSWAB ADDED                      JUL 02 ArDean Leith
C               FLIP QUESTION                     MAR 03 ArDean Leith
C               BAD IRECMRC4 & FLIP               SEP 03 ArDean Leith
C               SCALING                           JAN 05 ArDean Leith
C               I*8                               SEP 08 ArDean Leith
C               NPIX8                             DEC 08 ArDean Leith
C               BOTLEFT OPTION                    MAY 12 ArDean Leith
C               STREAM IO                         FEB 13 ArDean Leith
C               VOL BUG                           JUN 13 ArDean Leith
C               VOL BUG FIXED                     JUL 13 ArDean Leith
C               MODE 6 STACK SUPPORT              SEP 14 ArDean Leith
C               IPOSMRC INTEGER *8                JAN 15 ArDean Leith
C               BOTLEFT DEFAULT                   JUL 15 ArDean Leith
C               2015 STACK SUPPORT                JUL 15 ArDean Leith
C               STACK END BUG                     OCT 15 ArDean Leith
C               'MRCV'                            DEC 15 ArDean Leith
C               REWRITE FOR NATIVE MRC SUPPORT    JAN 20 ArDean Leith
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020 Health Research Inc.,                          *
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
C COPYFROMMRC(LUN1,LUN2, LUNDOC,LUNXM1,LUNXM2)
C                                                                      
C PURPOSE: CONVERTS MRC IMAGES TO SPIDER IMAGES, CRUDELY WRITTEN!!!
C
C NOTES: DATA IN MRC FILE
C        MODE   TYPES OF PIXEL IN IMAGE
C               0 : INTEGER*1 BYTES (UNSIGNED) 
C               1 : INTEGER*2       (SIGNED) 
C               2 : 32 BIT REALS
C               6 : INTEGER*2       (UNSIGNED)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE COPYFROMMRC(LUN1,LUN2, LUNDOC,LUNXM1,LUNXM2)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL        :: BUF
        COMMON /IOBUF/ BUF(NBUFSIZ)

        INTEGER                  :: LUN1,LUN2,LUNDOC,LUNXM1,LUNXM2

        LOGICAL                  :: VERBOSE_SAVE,ASKNAM
        LOGICAL                  :: IS_MRC,IS_ACTUALLY_A_VOL,IS_BARE
        LOGICAL                  :: FOUROK = .FALSE.

        CHARACTER (LEN=MAXNAM)   :: PROMPT 
        CHARACTER (LEN=MAXNAM)   :: FILNAM1,FILNAM2
        CHARACTER (LEN=1)        :: NULL = CHAR(0)
        CHARACTER (LEN=1)        :: DISP
        CHARACTER (LEN=4)        :: CAXIS
        CHARACTER (LEN=12)       :: CSTR
        INTEGER                  :: NC,NCC,NE,NOT_USED
        INTEGER                  :: NILMAX,IRTFLG,NSTACK1,NLET1,NINDX1
        INTEGER                  :: IFORM1,NX1,NY1,NZ1,NDUM,NGOT1,NGOT2
        INTEGER                  :: IMG1,NSTACK2,NLET2,IDUM,NINDX2,IMG2
        INTEGER,ALLOCATABLE      :: ILIST1(:),ILIST2(:)

        VERBOSE_SAVE = VERBOSE       ! SAVE CURRENT VERBOSITY

        NILMAX  = NIMAX              ! FROM CMLIMIT
        ALLOCATE(ILIST1(NIMAX),
     &           ILIST2(NIMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'COPYFROMMRC; ILIST...',2*NIMAX)
           RETURN
        ENDIF

C       OPEN FIRST INPUT MRC FILE 

        PROMPT  = 'MRC INPUT FILE OR TEMPLATE (E.G. STK@****)~~9'
        CALL FILERD(FILNAM1,NLET1,NULL, PROMPT,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       APPEND .mrc IF NEEDED
        IF (INDEX(FILNAM1,'.mrc') <= 0 .AND.
     &      INDEX(FILNAM1,'.MRC') <= 0 .AND.
     &      INDEX(FILNAM1,'.mrcs') <= 0 .AND.
     &      INDEX(FILNAM1,'.MRCS') <= 0 .AND.
     &      INDEX(FILNAM1,'.map') <= 0 .AND.
     &      INDEX(FILNAM1,'.MAP') <= 0) THEN
            FILNAM1 = FILNAM1(1:NLET1) // '.mrc'
            NLET1   = NLET1 + 4
        ENDIF

        NSTACK1 =  1   ! PARAMETER UNUSED ON INPUT
        DISP    = 'E'  ! DISP = 'E' DOES NOT STOP ON ERROR
        ASKNAM  = .FALSE.  
C       PRINT *, "copyfrommrc.f : 117: Calling OPFILES"
        IMG1 = 0    ! Passing uninitialized variable may result in 'INVALID IMAGE NUMBER'
        CALL OPFILES(0,LUN1,LUNDOC,LUNXM1,
     &               ASKNAM,FILNAM1,NLET1, DISP,
     &               IFORM1,NX1,NY1,NZ1,NSTACK1,
     &               FILNAM1,
     &               FOUROK, ILIST1,NILMAX, 
     &               NDUM,NGOT1,IMG1, IRTFLG) 
        PRINT *, "copyfrommrc.f : 124: OPFILES: NSTACK1:", NSTACK1
        IF (IRTFLG .NE. 0) RETURN

C       NSTACK1 RET:  -2   IS NON-STACK IMAGE,  -1 IS STACKED IMG,                  
C                    >= 0 IS CURRENT MAX. IMG # FOR STACK             

        !write(3,'(a,5i5)')' In copyfrommrc, nstack1,ngot1,img1:',
        !&                                      nstack1,ngot1,img1
 
C       BE SURE INPUT IS MRC
C       PRINT *, "copyfrommrc.f : 134: Calling LUNGETIS_MRC"
        CALL LUNGETIS_MRC(LUN1,IS_MRC,IRTFLG)
        IF (.NOT. IS_MRC) THEN    
           CALL ERRT(101,'INPUT IS NOT MRC FILE',NE)
           GOTO 999
        ENDIF
C       PRINT *, "copyfrommrc.f : 140: LUNGETIS_MRC ok"

C	OPEN FIRST SPIDER OUTPUT FILE
        NSTACK2 =  1   ! PARAMETER NOT USED ON INPUT
        DISP    = 'U'          
        IMG2    = IMG1

        IF (NSTACK1 >= 0) THEN
           PROMPT = 'SPIDER OUTPUT FILE OR TEMPLATE (E.G. STK@****)~~9'
        ELSE
           PROMPT = 'SPIDER OUTPUT FILE~'
        ENDIF
        CALL FILERD(FILNAM2,NLET2,NULL, PROMPT,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        CALL LUNGETISBARE_MRC(LUN1,IS_BARE,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

C       TREAT THIS IMAGE AS A STACK NOT A VOLUME (WHERE RELEVENT)
        IS_ACTUALLY_A_VOL = .FALSE.

        IF (NSTACK1 >= 0) THEN
C           MAY HAVE ABERRANT MRC HEADER
            IF  (FCHAR(9:12) == 'MRCV') THEN
               IS_ACTUALLY_A_VOL = .TRUE.

            ELSEIF (IS_BARE) THEN
               IMG2   = 1
               CALL RDPRI1S(IMG2,NOT_USED,
     &                  'FIRST IMAGE NUMBER IN SPIDER STACK',IRTFLG)
               IF (IRTFLG .NE. 0) GOTO 999
               IF (IMG2 <= 0) IS_ACTUALLY_A_VOL = .TRUE. 
            ELSE
               IMG2 = 1
            ENDIF
         ENDIF

         IF (IS_ACTUALLY_A_VOL) THEN
C           TREAT THIS IMAGE AS A VOLUME NOT A STACK
            NZ1     = NSTACK1
            NSTACK1 = -2
            CALL LUNSETNSTACK_MRC(LUN1, NZ1,NSTACK1,IRTFLG)
        ENDIF

        ASKNAM = .FALSE.
        CALL OPFILES(LUN1,LUN2,LUNDOC,LUNXM2, 
     &             ASKNAM,FILNAM2,NLET2,DISP,
     &             IFORM1,NX1,NY1,NZ1,NSTACK2,
     &             FILNAM2,
     &             FOUROK, ILIST2,NILMAX, 
     &             NDUM,NGOT2,IMG2, IRTFLG)
 
C       NSTACK2 RET:  -2   IS NON-STACK IMAGE,  -1 IS STACKED IMG,                  
C                    >= 0 IS CURRENT MAX. IMG # FOR STACK             

        !write(6,'(A,4i6)')' Out nstack2,ngot2,img2:',nstack2,ngot2,img2

        CALL WHICH_HAND_MRC(LUN1,FILNAM1,CAXIS,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        !write(3,*) ' Mrc data axis: ',caxis

        IF (NZ1 <= 1) THEN       ! IMAGE(S) INPUT
           CALL RDPRMC(CSTR,NC,.TRUE.,'DATA ORIGIN CORNER (UL/LL)', 
     &              NULL,IRTFLG)
        ELSE                     ! VOLUME(S) INPUT
           CALL RDPRMC(CSTR,NC,.TRUE.,  
     &              'DATA ORIGIN CORNER (UL/LL) & HANDEDNESS (L/R)',
     &              NULL,IRTFLG)
        ENDIF
        IF (IRTFLG .NE. 0)  RETURN

        NCC = INDEX(CSTR(1:NC),',')
        IF (NCC > 0) NCC = NCC - 1
        IF (INDEX(CSTR(1:NCC),'Y') > 0) THEN   ! LEGACY FLIP INPUT
           CALL ERRT(101,
     &         'DATA BYTE FLIP NO LONGER AVAILABLE HERE',IDUM)
           GOTO 999
        ENDIF

        IF (INDEX(CSTR(1:NC),'LL') > 0) CAXIS(1:2) = 'LL'

        IF (NZ1 > 1) THEN    
C          VOLUME OUTPUT
           CAXIS(4:4) = 'L'
           IF (INDEX(CSTR(3:NC),'R')  > 0) CAXIS(4:4) = 'R'

           !WRITE(NOUT,*)' DATA ORIGIN CORNER & HANDEDNESS: (',CAXIS,')'
        ELSE
           !WRITE(NOUT,*)' DATA ORIGIN CORNER: (',CAXIS(1:2),')'
        ENDIF

        CALL LUNSETHAND_MRC(LUN1,CAXIS,IRTFLG)
        CALL LUNSETPOS_MRC (LUN1,NGOT1,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 999

        NINDX1 = 1
        NINDX2 = 1
        DO                ! LOOP OVER ALL IMAGES/STACKS

C          DO NOT REPORT FILE INFO IF WHOLE STACK (VERBOSE IN COMMON)
	   !!!IF (NSTACK1 > 0 .AND. NSTACK2 >= 0) VERBOSE = .FALSE. 

C          COPY THE DESIRED NUMBER OF DATA RECORDS (MRC OK)
           DO IREC = 1,NY1 * NZ1
              CALL REDLIN(LUN1,BUF,NX1,IREC)
              CALL WRTLIN(LUN2,BUF,NX1,IREC)
           ENDDO
           !write(3,*)' in copyfrommrc,ilist.:   ',ilist1(:3),ilist2(:3) 
           !write(3,*)' in copyfrommrc,filnam1:  ',filnam1(1:20)
           !write(3,*)' in copyfrommrc,filnam2:  ',filnam2(1:20)
           !write(3,*)' in copyfrommrc,ngot1..2: ',ngot1,ngot2 
           !write(3,*)' In copyfrommrc,nindx1..2:',nindx1,nindx2 
           !write(3,*)' in copyfrommrc,img1..2:  ',img1,img2 

C          OPEN NEXT SET OF I/O FILES, UPDATES NINDX1 & NINDX2 
           CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                    .FALSE., LUNXM1,LUNXM2,
     &                    NGOT1,NGOT2,    NSTACK1,NSTACK2,  
     &                    LUN1,LUN1,LUN2, FILNAM1,FILNAM2,
     &                    IMG1,IMG2, IRTFLG)

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

