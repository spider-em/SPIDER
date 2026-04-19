
C++*********************************************************************
C
C OPENFIL.F     ADAPTED FROM OPEN3.FOR           OCT 1988 ArDean Leith
C               LONG FILENAMES                   OCT 1988 ArDean Leith
C               MERGED WITH OPENF                AUG 1996 ArDean Leith
C               -1, -3 & -7 FOURIER REMOVED      AUG 1996 ArDean Leith
C               F90 CHANGES                      APR 1998 ArDean Leith
C               USED LUNHDR                      JAN 1999 ArDean Leith
C               INDEXED STACKS                   JAN 2003 ArDean Leith
C               HEADER COPY                      FEB 2003 ArDean Leith
C               ENDEDNESS                        FEB 2003 ArDean Leith
C               REMOVED IRTFLG INPUT             APR 2004 ArDean Leith 
C               MPI OPEN BUG & BCAST CHANGES     OCT 2008 ArDean Leith 
C               OPENINLN INTEGER*8               OCT 2010 ArDean Leith 
C               NZ LIMITED TO 100000 BUG         JAN 2015 ArDean Leith 
C               MRC SUPPORT                      MAY 2019 ArDean Leith
C               OPEN ON 2 LUNS RECODED           AUG 2019 ArDean Leith
C               FLIP FLAG IN LUNFLIP FIXED       SEP 2025 ArDean Leith
C               DEBUG OUTPUT ADDED               APR 2026 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith & Tapu Shaikh *
C=*                                                                    *
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
C    OPENFIL(LUNT,FILNAM,LUNSENT,NX,NY,NZ,NSTK_FLG,
C            ITYPE,DISP,KEEPEXT,IRTFLG)
C
C    PURPOSE:       OPEN NEW OR EXISTING DATA FILE FOR RANDOM 
C                   ACCESS READING/WRITING.  CALLED FOR REGULAR FILES
C                   AND OVERALL STACK HEADERS.  THIS ROUTINE IS NOT
C                   CALLED FOR CREATING STACKED IMAGES WITHIN A
C                   STACK (ONLY FOR THE OVERALL STACK).  NEITHER
C                   IS IT USED FOR INLINE STACKS! 
C
C    PARAMETERS:
C        LUNT       UNIT TO COPY HEADER VALUES FROM               (SENT)
C        FILNAM     CHARACTER ARRAY, CONTAINING FILE NAME         (SENT)
C        LUNSENT    LOGICAL UNIT NUMBER FOR FILNAM.               (SENT)
C        NX,NY      DIMENSIONS OF FILE                (SENT OR RETURNED)
C        NZ         NUMBER OF SLICES IN IMAGE         (SENT OR RETURNED)
C        NSTK_FLG   STACK INDICATOR                    (SENT / RETURNED)
C                   ON INPUT: 
C                      FOR NEW FILE  0 : NOT A STACK, 
C                                   >0 : REGULAR STACK
C                                   <0 : INDEXED STACK, MAX TOTAL IMAGES
C                   ON OUTPUT: 
C                                  -2  : NOT A STACK 
C                                  >=0 : MAX IMAGE IN STACK
C
C        ITYPE        
C        DISP       CHARACTER VARIABLE CONTAINING DISPOSITION     (SENT)
C                   SEE OPFIL.F FOR POSSIBLE VALUES            
C        KEEPEXT    TO KEEP EXTENSION                             (SENT)
C        IRTFLG     ERROR RETURN FLAG.                        (RETURNED)
C                     0 : NORMAL RETURN
C                     1 : ERROR RETURN
C                     5 : NOT A SPIDER INPUT IMAGE
C
C    CODING:   BASED ON IMAGE PARAMETERS NX,NY, & NZ,  A
C              NEW FILE IS OPENED WITH IREC RECORDS, EACH NX * 4 
C              BYTES LONG.  IREC ALLOWS SPACE FOR THE 3-D DISTRIBUTION,
C              PLUS NORMAL IMAGE HEADER.
C
C    VARIABLES: ITYPE (TYPE)  FILE TYPE SPECIFIER. 
C               +1    R     2-D IMAGE
C               +3    R3    3-D VOLUME FILE
C               -11   O2    2-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -12   E2    2-D FOURIER TRANSFORM, MIXED RADIX EVEN
C               -21   O3    3-D FOURIER TRANSFORM, MIXED RADIX ODD
C               -22   E3    3-D FOURIER TRANSFORM, MIXED RADIX EVEN
C
C    CALL TREE:  
C
C    CALLED BY: openstk.f (3X),  opfilec.f (1X)
C
C                  MRC file
C    OPFILEC ----? -----------> OPENFIL_MRC --> ...
C                |
C                |                      Regular
C                | Not a stack          file
C                ?------------> OPENFIL ------> ...
C                |                 |
C                |                 |    Inline
C                |                 |    file
C                |                 `----------> INLNBUF 
C                |                              OPENINLN
C                ? File stack
C                |-------------> OPENSTK -----> OPENFIL -->
C                |                 
C                ? Inline stack
C                |-------------> OPNINSTK ----> ...
C
C       NOTE: OPENFIL SHOULD NEVER BE CALLED WITH MRC FILES     
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
        SUBROUTINE OPENFIL(LUNT,FILNAM,LUNSENT,NX,NY,NZ,NSTK_FLG,
     &                     ITYPE,DISP,KEEPEXT,IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'LABLOCK.INC'

        INTEGER                      :: LUNARA,LUNSTK,LUNARB,LUNFLIP
        COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

        INTEGER                      :: LUNT,LUNSENT
        CHARACTER(LEN=*)             :: FILNAM
        INTEGER                      :: NX,NY,NZ,NSTK_FLG,ITYPE
        CHARACTER(LEN=1)             :: DISP
        LOGICAL                      :: KEEPEXT
        INTEGER                      :: IRTFLG

        CHARACTER(LEN=MAXNAM)        :: FILNM
        CHARACTER(LEN=2*MAXNAM)      :: MSG 
        LOGICAL                      :: EX,CALLERRT,OPENED,GOTEXT
        INTEGER, PARAMETER           :: I_8 = SELECTED_INT_KIND(12)
        INTEGER(KIND=I_8)            :: NWORDS_8
        INTEGER(KIND=I_8), PARAMETER :: ZERO_8 = 0
        INTEGER                      :: LUN,IFLIP,ISTACK,LOCDOT,LOCSLASH
        INTEGER                      :: INLNED,MAXIM

        INTEGER                      :: NLET,IER,NE,LENREC,LUNNOW
        INTEGER                      :: IFLIPIN,IVAL,NDUM1,NDUM2,I
        INTEGER                      :: LENBYT,LABBYT,NRECS,INDXREC
        INTEGER                      :: LUNOPENED,LABREC,IRTFLGT

        INTEGER                      :: lnblnk,lenopenfile

        LOGICAL                      :: IS_MRC

C       INCLUSION FOR OPTIONAL MPI INITIALIZATION. MAY ALTER MYPID  
        INTEGER                      :: MYPID = -1
#include "MPI_INIT.INC"

        LUN = LUNSENT            ! IN CASE LUNSENT IS CONSTANT

C       SET FLAG FOR ERRONEOUS RETURN
        IRTFLG = 1
       
C       CHECK IF USER WANTS TO USE THE IN-LINE BUFFER.
        INLNED = 0

        IF (FILNAM(1:1) == '_') THEN
C          USE AN IN-LINE BUFFER.
           CALL INLNBUF(FILNAM,NLET,INLNED,IRTFLGT)
           IF (IRTFLGT .NE. 0)  RETURN

           IF (DISP == 'U' .OR. DISP == 'N') THEN
              EX  = .FALSE.
           ELSE
              EX  = .TRUE.
           ENDIF
           FILNM  = FILNAM(1:NLET)
           INLNED = 1

        ELSE
C          USE REGULAR I/O FILE
           INLNED = 0
           LUNC   = LUN

C          NULLIFY THE INLINE POINTER (CLOSEINLN IS INSIDE OPENINLN.F)
           CALL CLOSEINLN(LUN,IRTFLGT)

           LOCDOT   = INDEX(FILNAM,'.',BACK=.TRUE.)
           LOCSLASH = INDEX(FILNAM,'/',BACK=.TRUE.)

           GOTEXT   = (( LOCDOT > 1 ) .AND. ( LOCDOT > LOCSLASH))

           IF (GOTEXT .AND. KEEPEXT) THEN
              FILNM = FILNAM
              NLET  = lnblnk(FILNM)

           ELSE
C             APPEND DATEXC EXTENSION TO FILNAM
              CALL FILNAMANDEXT(FILNAM,DATEXC,FILNM,NLET,.TRUE.,IRTFLGT)
              IF (IRTFLGT .NE.0) RETURN
           ENDIF

C          SEE IF FILE EXISTS NOW
           IF (MYPID <= 0) THEN
              INQUIRE(FILE=FILNM(1:NLET),IOSTAT=IER,EXIST=EX)
           ENDIF

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CC#define SP_DBUGIO
#if defined (SP_DBUGIO)
           write(3,*)'  '
           write(3,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(3,*)' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(3,*)' In openfil   0 ; filnam:   ', trim(filnam)
           write(3,*)' In openfil   0 ; locdot:   ', locdot
           write(3,*)' In openfil   0 ; locslash: ', locslash
           write(3,*)' In openfil   0 ; gotext:   ', gotext
           write(3,*)' In openfil   0 ; keepext:  ', keepext
           write(3,*)' In openfil   0 ; disp:     ', disp
           write(3,*)' In openfil   0 ; ex:       ', ex
           write(3,*)' In openfil   0 ; filnm:    ', filnm(:nlet) 
           write(3,*)' In openfil   0 ; nstk_flg: ', nstk_flg 
           write(3,*) '  '
#endif


#ifdef USE_MPI
           IF (ONLYONE_RED) THEN
              CALL BCAST_MPI('OPENFIL','EX', EX,1,'L',ICOMM)
              CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
           ENDIF
#endif
           IF (IER .NE. 0) THEN
              MSG =  'OPENFIL; INQUIRY ERROR: ' // FILNM(1:NLET)
              CALL ERRT(101,MSG,NE)
              RETURN
           ENDIF
        ENDIF

C       SET ACTUAL LUN FOR THIS FILE
C       LUNARB(LUN) = LUN
        CALL LUNSETLUNB(LUN,LUN,IRTFLG)

        CALLERRT = (DISP(1:1) .NE. 'Z' .AND. DISP(1:1) .NE. 'E')



10      IF (DISP == 'U' .OR. DISP == 'N') THEN
C         USER WANTS A NEW FILE TO WRITE INTO ---------------------- NEW

          IF (EX) THEN
C           FILE IS IS OPENED WITH 'UNKNOWN' BUT ALREADY EXISTS; IT WILL
C           BE REPLACED.  OLD FILES ARE DELETED FIRST 
            IF (MYPID <= 0) THEN
               OPEN(UNIT=LUN,FILE=FILNM(1:NLET),STATUS='OLD')
               CLOSE(UNIT=LUN,STATUS='DELETE',IOSTAT=IER)
            ENDIF
#ifdef USE_MPI
            IF (ONLYONE_RED)
     &         CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
#endif

            IF (IER .NE. 0) THEN
               MSG =  'OPENFIL; DELETING FILE: ' // FILNM(1:NLET)
               CALL ERRT(101,MSG,NE)
               RETURN
            ENDIF
          ENDIF

C         INITIALIZE HEADER OBJECT FOR NEW IMAGE 
          CALL LUNGETIS_MRC(LUNT,IS_MRC,IRTFLG)
          IF (IS_MRC) THEN 
C            DO NOT COPY FROM LUNT IF IT IS A MRC FILE
             CALL LUNSETHDR(0,   LUN,NX,NY,NZ,ITYPE,NSTK_FLG,IRTFLGT)
          ELSE
             CALL LUNSETHDR(LUNT,LUN,NX,NY,NZ,ITYPE,NSTK_FLG,IRTFLGT)
          ENDIF
#if defined (SP_DBUGIO)
          write(3,*) '  '
          write(3,*) ' In openfil   1 ; After  lunsethdr '
          write(3,*) ' In openfil   1 ; lunt,lun:   ', lunt,lun
#endif
 

C         SET FILE'S FMIN.... TO AVOID FLT. PT. ERROR ON DEC
C         DOES NOT ALTER STATS: FMIN... IN COMMON BLOCK
          CALL LUNSETSTAT(LUN,0, 0.0,0.0,0.0,0.0,IRTFLGT)

C         SET NATIVE-ENDDEDNESS OF NEW FILE
          IFLIP = 0  ;
          IF (FCHAR(1:7) == 'CP TO O') THEN  
C              'CP TO OP END'  EXCEPTION 
             IFLIP = 1
          ENDIF
          CALL LUNSETFLIP(LUN, IFLIP,IRTFLG)
 

          IF (INLNED == 0) THEN
C             REGULAR FILE, NOT INLINE FILE --------------------

             LENREC = LENOPENFILE(NX*4)
             IF (MYPID <= 0) THEN
                OPEN(LUN,FILE=FILNM(1:NLET),STATUS='NEW',
     &               FORM='UNFORMATTED',
     &               ACCESS='DIRECT',IOSTAT=IER,RECL=LENREC)

#if defined (SP_DBUGIO)
             write(3,*)''
             write(3,*)' In openfil   2 ; Opened filnm: ', filnm(1:nlet)
#endif
             ENDIF

#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
                CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
             ENDIF
#endif
             IF (IER .NE. 0) THEN        
                MSG = 'OPENFIL; OPENING FILE: ' // FILNM(1:NLET)
                CALL ERRT(101,MSG,NE)
                RETURN
             ENDIF

C            FOR USING NEW FILE ON OLDER SPIDER RELEASES
             CALL LUNSET25(LUN,-1,IRTFLGT)

             ISTACK = 0
             IF (NSTK_FLG > 0 ) ISTACK = NSTK_FLG 
             CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT)


#if defined (SP_DBUGIO)
             write(3,*)'  '
             write(3,*)' In openfil   3 ; nstk_flg: ', nstk_flg
             write(3,*)' In openfil   3 ; istack:   ', istack
#endif

          ELSE
C            SET UP INLINE BUFFER AND TIE IT TO LUN

C            SAVE ISTACK IN NON-VOLATILE PART OF HEADER OBJECT
             CALL LUNGETISTACK(LUN,ISTACK,IRTFLGT)
             CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT)

             CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,LABBYT,
     &                      LENBYT,IRTFLGT)
             IF (IRTFLGT .NE. 0) THEN
                write(6,*) 'nrecs,indxrec,lenbyt,errorflag:'
                write(6,*)  nrecs,indxrec,lenbyt,irtflgt
                RETURN
             ENDIF

             NWORDS_8 = (NRECS + INDXREC) * LENBYT / 4

             CALL OPENINLN(LUN,INLNED,.TRUE.,NX,NWORDS_8,
     &                     .TRUE.,IRTFLGT)
             IF (IRTFLGT .NE. 0)  RETURN
          ENDIF

#if defined (SP_DBUGIO)
          write(3,*)'  '
          write(3,*)' In openfil  04 ; Calling lunwrthed'
#endif


C         PUSH OVERALL HEADER OBJECT INFO INTO NEW FILE
          CALL LUNWRTHED(LUN,NX,0,IRTFLGT)

#if defined (SP_DBUGIO)
          write(3,*)' ' 
          write(3,*)' In openfil  10 ; Opened filnm: ', filnm(1:nlet)
          write(3,*)' In openfil  10 ; nstk_flg:     ', nstk_flg
          write(3,*)' In openfil  10 ; istack:       ', istack

          write(3,*)' In openfil  10 ; lun:          ', lun 
          write(3,*)' In openfil  10 ; lunflip(lun): ', lunflip(lun)
          write(3,*)' ' 
#endif

          GOTO 2000




C         --------------------------------------------------------- OLD
C         --------------------------------------------------------- OLD

        ELSEIF (DISP == 'O' .OR. DISP == 'K' .OR.
     &          DISP == 'Z' .OR. 
     &          DISP == 'E' .OR. DISP == 'M') THEN

C         FILE EXISTS, AND IS ACCESSED WITH 'OLD', OPEN THE FILE

          IF (.NOT. EX) THEN
C            ERROR -- FILE DOES NOT EXIST, BUT BEING OPENED WITH 'OLD
             IF (MYPID <= 0)
     &          WRITE(NOUT,*) '*** FILE NOT FOUND: ',FILNM(1:NLET)

C            FOR DISP=Z, DO NOT STOP THE BATCH JOB BY CALLING ERRT
             IF (CALLERRT) THEN
                MSG = 'OPENFIL; FILE NOT FOUND: ' // FILNM(1:NLET)
                CALL ERRT(101,MSG,NE)
             ENDIF
             RETURN
          ENDIF    ! END OF: .NOT. EX


          IF (INLNED .NE. 0) THEN
C            INLINE FILE ACCESS ---------------------------- OLD INLINE

C            USE EXISTING INLINE BUFFER, TIE IT TO LUN & RETRIEVE NX
             CALL OPENINLN(LUN,INLNED,.FALSE.,NX,ZERO_8,
     &                     CALLERRT,IRTFLGT)
             IF (IRTFLGT .NE. 0)  RETURN



          ELSE 

C            REGULAR SPIDER FILE ACCESS  ------------------ OLD REGULAR

C            SEE IF FILE ALREADY OPEN ON DIFFERENT LUN THAN CURRENT
             IF (MYPID <= 0) THEN
                INQUIRE(FILE=FILNM(1:NLET),IOSTAT=IRTFLG,OPENED=OPENED,
     &                  NUMBER=LUNOPENED)
             ENDIF   ! END OF: MYPID <= 0

#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
               CALL BCAST_MPI('OPENFIL','IRTFLG',   IRTFLG, 1,'I',ICOMM)
               CALL BCAST_MPI('OPENFIL','OPENED',   OPENED, 1,'L',ICOMM)
               CALL BCAST_MPI('OPENFIL','LUNOPENED',LUNOPENED,
     &                         1,'I',ICOMM)
             ENDIF 
#endif

             IF (OPENED) THEN
C               MUST REDIRECT OLD LUNS TO CURRENT LUN BEFORE OPENING 
                IF (MYPID <= 0)  CLOSE(LUNOPENED)
                DO I = 1,100
                   CALL LUNGETLUNB(I,LUNNOW,IRTFLG)
                   IF (LUNNOW == LUNOPENED) THEN
                      CALL LUNSETLUNB(I,LUN,IRTFLG)
                   ENDIF
                ENDDO                      
             ENDIF   ! END OF: IF(OPENED)
             

C            NOW WE CAN OPEN THIS FILE ON LUN
             IF (MYPID <= 0) THEN
                LENREC = LENOPENFILE(256*4)  ! USED FOR SPIDER FILE

                OPEN(LUN,FILE=FILNM(1:NLET),STATUS='OLD',
     &               ACCESS='DIRECT',
     &               FORM='UNFORMATTED',RECL=LENREC,IOSTAT=IER)
             ENDIF  
#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
                CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
             ENDIF
#endif

             IF (IER .NE. 0) THEN
C              UNKNOWN OPENING ERROR
               MSG = 'OPENFIL; OPENING FILE: ' // FILNM(1:NLET)
               CALL ERRT(101,MSG,NE)
               RETURN
             ENDIF  ! END OF: IER .NE. 0 

          ENDIF  ! END OF: OLD REGULAR SPIDER FILE ACCESS  ------------


#if defined (SP_DBUGIO)
          write(3,*)' In openfil 12 ; Calling lunredhed -- Main header'
          write(3,*)'  '
#endif

C         READ 1024 BYTES FROM OVERALL FILE HEADER (UNFLIPPED)
          CALL LUNSETFLIP(LUN,0,IRTFLG)
          CALL LUNREDHED(LUN,256,0,.TRUE.,IRTFLGT)
          IF (IRTFLGT .NE. 0 .AND. MYPID <= 0) THEN
             WRITE(NOUT,*) '*** ERROR READING HEADER OF: ',FILNM(:NLET)
             IF (IRTFLGT .NE. 0) RETURN
          ENDIF    ! END OF: IRTFLGT .NE. 0 ......

C         NEED ITYPE
          CALL LUNGETTYPE(LUN,ITYPE,IRTFLGT)
C         NEED NX,NY & NZ  VALUES 
          CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLGT)
          IF (IRTFLGT .NE. 0) RETURN

C         CHECK FOR VALID DATA IN  HEADER
          CALL CHECK_FILE_ACCESS(ITYPE, NX,NY,NZ, 
     &                           .FALSE.,VERBOSE, IRTFLG)
          IF (IRTFLG .EQ. 1) RETURN

#if defined (SP_DBUGIO)
          write(3,*)' In openfil 13 ; nx,ny:  ',nx,ny
          write(3,*)' In openfil 13 ; itype:  ',itype
          write(3,*)' In openfil 13 ; irtflg: ',irtflg
#endif
 
          IF (IRTFLG .NE. 0) THEN
C           EITHER A NON-NATIVE SPIDER IMAGE OR NOT A SPIDER/MRC FILE

C           FLIP BYTES IN CURRENT HEADER OBJECT 
            CALL LUNFLIPHDR(LUN,IRTFLGT)
            CALL LUNSETFLIP(LUN,1,IRTFLG)

C           GET ITYPE & SIZE FROM FLIPPED HEADER
            CALL LUNGETTYPE(LUN,ITYPE,IRTFLGT)
            CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLGT)
            IF (IRTFLGT .NE. 0) RETURN

C           CHECK FOR VALID DATA IN FLIPPED HEADER 
            CALL CHECK_FILE_ACCESS(ITYPE, NX,NY,NZ, 
     &                             .TRUE.,VERBOSE, IRTFLG)
#if defined (SP_DBUGIO)
            write(3,*)' '
            write(3,*)' In openfil 14 ; nx,ny:  ',nx,ny
            write(3,*)' In openfil 14 ; itype:  ',itype
            write(3,*)' In openfil 14 ; irtflg: ',irtflg
#endif

            IF (IRTFLG == 5 .OR. IRTFLG == 1) THEN
C              RETURN ERROR FLAG FOR INVALID OR NON-SPIDER IMAGE
               FILNAM = FILNM(1:NLET)
               RETURN
            ENDIF 

#if defined (SP_DBUGIO)
            write(3,*)' '
            write(3,*)' In openfil 15 ; nx,ny: ',nx,ny
            write(3,*)' In openfil 15 ; itype: ',itype
            write(3,*)' In openfil 15 ; lunflip(lun): ',lunflip(lun)
#endif
            IF (VERBOSE) 
     &         WRITE(NOUT,*)' NON-NATIVE BYTE ORDERED SPIDER FILE'

          ENDIF



          IF (INLNED == 0) THEN   ! NOT INLINE FILE
C            NOT INLINE FILE, REOPEN FILE WITH CORRECT RECORD LENGTH
             CLOSE(LUN)

             LENREC  = LENOPENFILE(NX*4)
             IFLIPIN = LUNFLIP(LUN)

             IF (MYPID <= 0) THEN

#if defined (SP_DBUGIO)
                write(3,*)' '
                write(3,*)' In openfil 16 ; iflipin: ',iflipin
                write(3,*)' In openfil 16 ; nx,:     ',nx
                write(3,*)' In openfil 16 ; lenrec:  ',lenrec
                write(3,*)' In openfil 16 ; mypid:   ',mypid
                write(3,*)' In openfil 16 ; Reopening file'
                write(3,*)' '
#endif

                OPEN(LUN,FILE=FILNM(1:NLET),STATUS='OLD',
     &               FORM='UNFORMATTED',
     &               ACCESS='DIRECT',IOSTAT=IER,RECL=LENREC)

                IF (IER .NE. 0) THEN        
                   MSG = 'OPENFIL; REOPENING FILE: ' // FILNM(1:NLET)
                   CALL ERRT(101,MSG,NE)
                   RETURN
                ENDIF
             ENDIF   

#ifdef USE_MPI
             IF (ONLYONE_RED) THEN
               CALL BCAST_MPI('OPENFIL','IER',IER,1,'I',ICOMM)
             ENDIF
#endif

C            READ HEADER FROM OVERALL SPIDER FILE HEADER (FLIPPED)
c             CALL LUNREDHED(LUN,256,LENRED,.TRUE.,IRTFLGT)
c             IF (IRTFLGT .NE. 0 .AND. MYPID <= 0) THEN
c              WRITE(NOUT,*) '*** ERROR READING HEADER OF:         ',FILNM(:NLET)
c              IF (IRTFLGT .NE. 0) RETURN
c             ENDIF    ! END OF: IRTFLGT .NE. 0 ......

C            FIND IF THIS FILE CONTAINS STACK    
             CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG) 


CC#define SP_DBUGIO_STK 
#if defined (SP_DBUGIO_STK)
             write(3,*)' In openfil 17 ; istack shud be > 0'
             write(3,*)' In openfil 17 ; istack:  ', istack
             write(3,*)' In openfil 17 ; maxim:   ', maxim 
             write(3,*)' '
#endif

C            SET FINAL LUNARA OFFSET VALUE FOR USE BY REDLIN/WRTLIN.
C            FOR INDEXED STACKS THIS INCLUDES INDEX OFFSET
             IF (ISTACK .NE. 0) THEN
                CALL LUNSETIMGOFF(LUN,0,NX,IRTFLGT)
             ELSE
                IVAL = 1
                CALL LUNSETIMGOFF(LUN,IVAL,NX,IRTFLGT)
             ENDIF

          ENDIF   ! END OF: INLNED == .FALSE.

        ELSE  
           MSG = 'OPENFIL; PGM. ERROR, UNKNOWN DISP: ' // DISP
           CALL ERRT(101,MSG,NE)
           RETURN
        ENDIF    ! END OF:  (DISP == 'O' .OR. ......
        
C       FIND NATIVE-ENDIANESS OF INPUT     
C       CALL LUNGETFLIP(LUN,IFLIPIN,IRTFLG)

#if defined (SP_DBUGIO_STK)
        write(3,*)' In openfil 20 ; lun:          ', lun
        write(3,*)' In openfil 20 ; iflipin:      ', iflipin
        write(3,*)' In openfil 20 ; lunflip(lun): ', lunflip(lun)
#endif
 

C       FIND IF THIS FILE CONTAINS STACK    
        CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG) 
         

2000    CONTINUE ! ---------------------------------------------
                 ! ---------------------------------------------



#if defined (SP_DBUGIO_STK)
        write(3,*)'   ' 
        write(3,*)' ------------ AFTER 2000 IN OPENFIL --------' 
        write(3,*)' In openfil 22: nstk_flg: ', nstk_flg 
        write(3,*)' In openfil 22: istack:   ', istack
        write(3,*)' In openfil 22: maxim:    ', maxim 
        write(3,*)'   ' 
#endif

C       GET MAXIM VALUE FROM HEADER OBJECT (26) 
C       SETS MAXIM IN NON-STATIC AREA (260)
        CALL LUNCOPYMAXIM(LUN,MAXIM,IRTFLGT)

C       SAVE ISTACK IN NON-VOLATILE PART OF HEADER OBJECT
        CALL LUNGETISTACK(LUN,ISTACK,IRTFLGT)
        CALL LUNSETSTKALL(LUN,ISTACK,IRTFLGT)

C       SAVE MAXIM IN NON-VOLATILE PART OF HEADER OBJECT
        MAXIM = MAX(ISTACK,NSTK_FLG)
        CALL LUNCOPYMAXIM(LUN,MAXIM,IRTFLGT)


C       SET FINAL LUNARA OFFSET VALUE FOR USE BY REDLIN/WRTLIN.
C       FOR INDEXED STACKS THIS INCLUDES INDEX OFFSET
        IF (ISTACK .NE. 0) THEN
           CALL LUNSETIMGOFF(LUN,0,NX,IRTFLGT)
        ELSE
           IVAL = 1
           CALL LUNSETIMGOFF(LUN,IVAL,NX,IRTFLGT)
        ENDIF

C       GET RECORD INFO (CAN BE FROM OVERALL HEADER)
        CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,NDUM1,NDUM2,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       PUT COMMON VALUES INTO COMMON AREA (NOT NEEDED IN FUTURE?)
C       SET STATS: FMIN... IN COMMON BLOCK AND FILE HEADER
        CALL LUNSETCOMMON(LUN,IRTFLGT)

C       WRITE OUT FILE OPENING INFO
        CALL LUNSAYINFO(LUN,IRTFLGT)
      
CC#define SP_DBUGIO_STK
#if defined (SP_DBUGIO_STK)
        write(3,*)'   ' 
        write(3,*)' In openfil 30 : nstk_flg:    ', nstk_flg 
        write(3,*)' In openfil 30 : istack:      ', istack
        write(3,*)' In openfil 30 : maxim:       ', maxim 
        write(3,*)' In openfil 30 ; labrec:      ', labrec
        write(3,*)' '
#endif

C       NEED TO RETURN NSTK_FLG & MAXIM
 
        IF (ISTACK == 0) THEN
C          NOT A STACK, RETURN -2
           NSTK_FLG = -2

#if defined (SP_DBUGIO_STK)
           write(3,*)' In openfil 31; istack:       ', istack
           write(3,*)' In openfil 31; nstk_flg:     ', nstk_flg 
           write(3,*)' '
#endif
        ELSE
C          RETURN CURRENT HIGHEST IMAGE NUMBER IN NSTK_FLG

           MAXIM    = MAX(MAXIM,ISTACK)
           NSTK_FLG = MAXIM

#if defined (SP_DBUGIO_STK)
           write(3,*)' In openfil 32 ; istack:      ', istack
           write(3,*)' In openfil 32 ; nstk_flg:    ', nstk_flg 
#endif
        ENDIF


#if defined (SP_DBUGIO)
        write(3,*)' '
        write(3,*)' In openfil 99 ; labrec:       ', labrec
        write(3,*)' In openfil 99 ; nstk_flg:     ', nstk_flg 
        write(3,*)' In openfil 99 ; istack:       ', istack
        write(3,*)' In openfil 99 ; maxim:        ', maxim 
        write(3,*)' In openfil 99 ; iflipin:      ', iflipin
        write(3,*)' In openfil 99 ; lunflip(lun): ', lunflip(lun)
        write(3,*)' Leaving openfil 99 --------------- '
        write(3,*)' '
#endif


C       SET FLAG FOR NORMAL RETURN        
        IRTFLG = 0

        END



C++*********************************************************************
C
C CHECK_FILE_ACCESS    EXTRACTED FROM OPENFIL MAR 2026 ArDean Leith
C
C **********************************************************************
C
C   CHECK_FILE_ACCESS(NX,NY,NZ,  ITYPE,CALLERR,VERBOSE,IRTFLG)
C
C    PURPOSE:       CHECK IF RECORD FROM FILE HAS IMPROPER ENNDIANESS
C                   FOR CURRENT EXECUTABLE OR IS NOT A SUPPORTED 
C                   FILE FORMAT 
C
C    PARAMETERS:
C        NX,NY,NZ   DIMENSIONS OF FILE                        (SENT)
C        ITYPE                                                (SENT)
C        CALLERR                                              (SENT)
C        VERBOSE                                              (SENT)
C        IRTFLG     ERROR RETURN FLAG.                        (RETURNED)
C                     0 : NORMAL RETURN
C                     1 : ERROR RETURN
C                     5 : ERROR RETURN: NOT A SPIDER/MRC  IMAGE
C
C    VARIABLES: ITYPE (TYPE)  FILE TYPE SPECIFIER. 
C               +1    R     2-D IMAGE
C               +3    R3    3-D VOLUME FILE 
C                  
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
       SUBROUTINE CHECK_FILE_ACCESS(ITYPE, NX,NY,NZ, 
     &                             CALLERR,VERBOSET,IRTFLG)

       IMPLICIT NONE
       INCLUDE 'CMBLOCK.INC'

       INTEGER  :: ITYPE,NX,NY,NZ
       LOGICAL  :: CALLERR , VERBOSET
       INTEGER  :: IRTFLG
  
       INTEGER  :: NE = 0

C       SET FLAG FOR NORMAL RETURN        
        IRTFLG = 0
 
         IF (ITYPE == 0 ) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(102,'INVALID HEADER, BAD FILE TYPE',ITYPE)
            ENDIF
            IRTFLG = 5   ! RETURN ERROR FLAG 

         ELSEIF (ITYPE == -1 .OR. ITYPE == -3 .OR. ITYPE == -7) THEN
C           READING OBSOLETE FORMAT FOURIER FILE
            WRITE(NOUT,96)
96          FORMAT(' *** CAN NOT READ OBSOLETE FOURIER FORMAT',/,
     &            '*** CONVERT FOURIER FILE TO REAL FORMAT USING ',
     &            'ORIGINAL VERSION OF SPIDER.'/)

            CALL ERRT(100, 'OPENFIL',NE)
            IRTFLG = 1   ! RETURN ERROR FLAG 
            RETURN

         ELSEIF ( NX <= 0 ) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(101,'INVALID HEADER, NX',NX)
            ENDIF
            IRTFLG = 5   ! RETURN ERROR FLAG 

         ELSEIF ( NY <= 0 ) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(101,'INVALID HEADER, NY',NY)
            ENDIF
            IRTFLG = 5   ! RETURN ERROR FLAG 

         ELSEIF ( NZ <= 0 ) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(101,'INVALID HEADER, NZ',NX)
            ENDIF
            IRTFLG = 5   ! RETURN ERROR FLAG 

         ELSEIF (NX > 1000000) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(102,'INVALID HEADER, NX (LIMIT=10000000)',NX)
            ENDIF
            IRTFLG = 5   ! RETURN ERROR FLAG 

         ELSEIF (NY > 1000000) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(102,'INVALID HEADER, NY (LIMIT=10000000)',NY)
            ENDIF
            IRTFLG = 5   !RETURN ERROR FLAG 

         ELSEIF (NZ > 1000000) THEN
C           NON-NATIVE ENDIAN OR NOT A SPIDER FILE
            IF (CALLERR) THEN
               CALL ERRT(102,'INVALID HEADER, NZ (LIMIT=10000000)',NZ)
            ENDIF
            IRTFLG = 5   ! RETURN ERROR FLAG 
         ENDIF
             
         IF (IRTFLG == 5) THEN
           IF (VERBOSET) THEN 
              WRITE(NOUT,*)' NON-NATIVE ENDIAN OR NOT A SPIDER FILE'
           ENDIF
         ENDIF

         END
