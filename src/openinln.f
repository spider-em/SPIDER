
C++*********************************************************************
C
C OPENINLN.F       AUTHOR                                  ARDEAN LEITH
C                  INCREASED NUMINLNT               AUG 02 ARDEAN LEITH
C                  SUPPORTS INTEGER 8               OCT 10 ARDEAN LEITH                            
C                  TOBEDEALLOC OVERFLOW BUG         JUL 12 ARDEAN LEITH                            
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
C    OPENINLN(LUN,NUMBUF,NEEDNEW,NSAM,NWORDS8,NEEDERRT,IRTFLG))
C
C    PURPOSE:       TO OPEN A INLINE BUFFER FOR READING/WRITING 
C                   VIA REDLIN AND WRTLIN.
C
C    PARAMETERS:    LUN      EQUIVALENT LUN                    (SENT)
C                   NUMBUF   BUFFER NUMBER                     (SENT)
C                   NEEDNEW  LOGICAL FLAG THAT WANT NEW BUFFER (SENT)
C                   NSAM     LENGTH OF ROW USED IN BUFFER      (SENT)
C                   NWORDS8  SIZE OF BUFFER                    (SENT)
C                   NEEDERRT NEED TO CALL ERRT AND END IF ERROR(SENT)
C                   IRTFLG   ERROR FLAG                        (RETURNED)
C
C    NOTE:  jan 2011 should we zero the array so as to better detect
C           non-existing images??
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      MODULE INLN_INFO

         SAVE

         INTEGER, PARAMETER    :: NUMINLNT = 99
         INTEGER, PARAMETER    :: NUMLUNST = 100

         TYPE REAL_POINTER
            REAL, POINTER      :: IPT(:) 
         END TYPE REAL_POINTER

         TYPE(REAL_POINTER)    :: LOCBUF(NUMINLNT)
         TYPE(REAL_POINTER)    :: LOCLUN(NUMLUNST)

         INTEGER               :: NSAMLUN(NUMLUNST)

         INTEGER, PARAMETER    :: I_8 = SELECTED_INT_KIND(12)
         INTEGER(KIND=I_8)     :: NWORDLUN(NUMLUNST)

      END MODULE INLN_INFO

C     ---------------- OPENINLN -----------------------------

      SUBROUTINE OPENINLN(LUN,NUMBUF,NEEDNEW,NSAM,NWORDS8,
     &                    NEEDERRT,IRTFLG)

      USE INLN_INFO
      INCLUDE 'INLN_INFO.INC'

      INCLUDE 'CMBLOCK.INC'

      LOGICAL           :: NEEDNEW,NEEDERRT,DUPED
      INTEGER           :: LUN,NUMBUF,NSAM,IRTFLG
      INTEGER(KIND=I_8) :: NWORDS8
      INTEGER           :: NDEALLOC
      REAL,  POINTER    :: RPOINTER(:)
      INTEGER(KIND=I_8) :: NWORDBUF(NUMINLN)
      INTEGER           :: TOBEDEALLOC(NUMINLN)

      SAVE NWORDBUF,TOBEDEALLOC,NDEALLOC

      DATA NDEALLOC/0/

      IRTFLG = 1
      IF (LUN .LE. 0 .OR. LUN .GT. 100) THEN
C        NO SUCH LUN
         CALL ERRT(102,'PGM ERROR, BAD LUN IN OPENINLN',LUN)
         RETURN

      ELSEIF (NUMBUF .LE. 0 .OR. NUMBUF .GT. NUMINLN) THEN
C        NO SUCH BUFFER
         ITEMP = NUMINLN
         WRITE(NOUT,90) ITEMP,NUMBUF
90       FORMAT(' *** INLINE BUFFERS AVAILABLE:  1...',I3,
     &          ' BUT YOU WANTED: ',I4)
         CALL ERRT(100,'OPENINLN',NE)
         RETURN
      ENDIF

      IF (NEEDNEW .AND. (NSAM == 0 .OR. NWORDS8 == 0))THEN
C        WANT TO DEALLOCATE AND CLOSE AN OPEN INLINE BUFFER,
C        BUT WE WILL WAIT TILL NEXT ALLOCATE TO ACTUALLY 
C        DEALLOCATE, SO AS TO MINIMIZE MEMORY LEAK

         IF (NDEALLOC > 0) THEN
            DO I = 1,NDEALLOC
               IF ( NUMBUF  == TOBEDEALLOC(I)) THEN
C                 ALREADY IN LIST, DO NOT ADD DUPLICATE
                  IRTFLG = 0
                  RETURN
               ENDIF
            ENDDO
         ENDIF

C        PUT BUFFER NUMBER ON DELETION LIST
         NDEALLOC               = NDEALLOC + 1
         TOBEDEALLOC(NDEALLOC)  = NUMBUF

C        ADDED DEC 06 FOR BUFFER IN AND OUT USE al
         NWORDBUF(NUMBUF)  = 0
         NSAMBUF(NUMBUF)   = 0
         IRECBUF(NUMBUF)   = 0
         LABRECBUF(NUMBUF) = 0
         IRTFLG            = 0
         RETURN

      ELSEIF (NDEALLOC > 0) THEN
C        DEALLOCATE ALL BUFFERS ON DELETION LIST

         DO I = 1,NDEALLOC
            NTEMP = TOBEDEALLOC(I)
            DUPED = .FALSE.
            IF (I > 1) THEN
               DO J=1,I-1
                  IF (TOBEDEALLOC(J) == NTEMP) DUPED = .TRUE.
               ENDDO
            ENDIF  
            IF (.NOT. DUPED) THEN     
               RPOINTER => LOCBUF(NTEMP)%IPT
               IF (ASSOCIATED(RPOINTER)) THEN
                  DEALLOCATE(RPOINTER,STAT=IRTFLGT)
               ENDIF

C              ZERO BUFFER SIZE, POINTER, NSAM, ETC...
               NULLIFY(LOCBUF(NTEMP)%IPT)
               NWORDBUF(NTEMP)  = 0
               NSAMBUF(NTEMP)   = 0
               IRECBUF(NTEMP)   = 0
               LABRECBUF(NTEMP) = 0
            ENDIF
         ENDDO
         NDEALLOC = 0
      ENDIF

1111  CONTINUE
      IF (NEEDNEW .AND. 
     &   (ASSOCIATED(LOCBUF(NUMBUF)%IPT)) .AND.
     &    (NWORDBUF(NUMBUF) .NE. NWORDS8)) THEN
C        FREE EXISTING SPACE ALLOCATION FOR THIS FILE
         RPOINTER => LOCBUF(NUMBUF)%IPT
         DEALLOCATE(RPOINTER,STAT=IRTFLGT)
         IF (IRTFLGT .NE. 0) RETURN

         NWORDBUF(NUMBUF) = 0
         NULLIFY(LOCBUF(NUMBUF)%IPT) 
      ENDIF

      IF (NEEDNEW .AND. .NOT. ASSOCIATED(LOCBUF(NUMBUF)%IPT)) THEN
C        MUST ALLOCATE SPACE FOR THIS BUFFER
         !print *,' allocating: ',nwords8

         ALLOCATE(RPOINTER(NWORDS8),STAT=IRTFLGT)
         IF (IRTFLGT .NE. 0) THEN
            WRITE(NOUT,92) NWORDS8,NUMBUF  
92          FORMAT(' *** CAN NOT ALLOCATE:',I12,
     &             '  WORDS FOR INLINE FILE: _',I1)  
            CALL ERRT(102,'ALLOCATION FAILED FOR INLINE FILE',NUMBUF)
            RETURN
         ENDIF

         LOCBUF(NUMBUF)%IPT => RPOINTER
         NWORDBUF(NUMBUF)   = NWORDS8
         NSAMBUF(NUMBUF)    = NSAM
         !print *,' nwordbuf(',numbuf,'):',NWORDBUF(NUMBUF)
      ENDIF

      RPOINTER        => LOCBUF(NUMBUF)%IPT
      LOCLUN(LUN)%IPT => RPOINTER
      ISINLINE(LUN)   =  .TRUE.

      IF (.NOT. NEEDNEW) THEN
C        RETRIEVE RECORD LENGTH FROM PREVIOUSLY SET VALUE
         NSAM = NSAMBUF(NUMBUF)
         IF (NSAM .LE. 0) THEN
C           NO PREVIOUS RUN-TIME FILE AVAILABLE
            WRITE(NOUT,93) NUMBUF
93          FORMAT(' *** NO EXISTING INLINE BUFFER: ',I3)
            IF (NEEDERRT) 
     &         CALL ERRT(102,'NO EXISTING INLINE BUFFER',NUMBUF)
            RETURN
         ENDIF
      ENDIF

C     SET LENGTH OF RECORD AND BUFFER FOR THIS LUN
      NSAMLUN(LUN)    = NSAM
      NWORDLUN(LUN)   = NWORDBUF(NUMBUF)
      LOCLUN(LUN)%IPT => LOCBUF(NUMBUF)%IPT
      ISINLINE(LUN)   =  .TRUE.
      !print *,' nwordlun(',lun,'):',NWORDLUN(LUN)

      IF (NWORDS8 > NWORDLUN(LUN)) THEN
C        NOT ENOUGH SPACE IN BUFFER
         WRITE(NOUT,94) NUMBUF,NWORDLUN(LUN),NWORDS8
94       FORMAT(' *** INLINE BUFFER: ',I4,
     &          '  LIMITED TO:  ',     I12,
     &          ' WORDS, WANTED: ',    I12)
         CALL ERRT(100,'INLINE BUFFER OVERFLOW',NE)
         RETURN
      ENDIF

      IRTFLG = 0

      RETURN
      END



C++********************************************************************
C   CLOSEINLN.F                                
C **********************************************************************
C
C     CLOSEINLN(LUN,IRTFLG)
C
C     PURPOSE:    CLOSE AN INLN BUFFER USE
C
C **********************************************************************

       SUBROUTINE CLOSEINLN(LUN,IRTFLG)

       USE INLN_INFO
       INCLUDE 'INLN_INFO.INC'

       NULLIFY(LOCLUN(LUN)%IPT)
       ISINLINE(LUN) = .FALSE.

       IRTFLG = 0
       RETURN
       END


C++********************************************************************
C   INLN_WRTLIN.F                                
C **********************************************************************
C
C     INLN_WRTLIN(LUN,BUF,NB,NREC)
C
C     PURPOSE:    READ A LINE OF FLOATING POINT NUMBERS FROM BUFFER
C
C     PARAMETERS:
C        LUN    LOGICAL UNIT NUMBER OF FILE BEING READ
C        BUF    BUFFER WHERE RECORD IS READ IN
C        NB     NUMBER OF VALUES IN RECORD TO BE READ
C        NREC   RECORD TO BE READ
C 
C        IERR   ERROR CODE 1 IS RETURNED IN CASE OF ERROR
C               IERR IS DEFINED IN COMMON /IOERR/IERR
C 
C       MANY COMMANDS IN SPIDER READ A FILE, MANIPULATE THE DATA, AND
C       WRITE THE RESULTS INTO A FILE, ONLY TO HAVE THE NEXT COMMAND
C       READ THE DATA BACK AGAIN IN CORE TO MANIPULATE IT AND STORE
C       IT BACK TO A FILE THEN START THE CYCLE AGAIN...
C       TO SPEED UP THE PROCESS, THE USER HAS THE ABILITY TO HAVE THAT 
C       DATA KEPT IN INCORE MEMORY UNTIL HE/SHE DECIDES THAT HE/SHE 
C       WANTS IT IN A FILE.
C
C       DATA FROM A FILE IS STORED EXACTLY THE SAME WAY AS IT IS IN A
C       DISK FILE (I.E: A HEADER FOLLOWED BY NROWS OF DATA OF LENGTH NSAMS
C       REPEATED NSLICE TIMES).
C 
C--*******************************************************************

      SUBROUTINE INLN_WRTLIN(LUN,BUF,NB,NREC)

C     USE INLINE BUFFER COMMON AREA
      USE INLN_INFO

      REAL   BUF(NB)

      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      REAL, POINTER :: IPTNOW(:)
      INTEGER *8    :: IGO

      IPTNOW => LOCLUN(LUN)%IPT

C     GET THE CORRECT OFFSET INTO THE BUFFER.
C     IGO = (LUNARA(LUN) + LUNSTK(LUN) + NREC - 1) * NSAMLUN(LUN)  
      IGO = (LUNARA(LUN) + LUNSTK(LUN) + NREC - 1) 
      IGO = IGO * NSAMLUN(LUN)  

      IF ((IGO + NB - 1 ) > NWORDLUN(LUN)) THEN
C        GOES OUT OF RESERVED BUFFER!
         WRITE(NOUT,"(A,I3,A,I6,A,I5,A,I4,A,I6,A,I6,A,I12)") 
     &      '*** INLN_WRTLIN ARRAY OVERFLOW, LUN: ',LUN,
     &      '   LENGTH:',        IGO + NB,
     &      '   NREC: ',         NREC,
     &      '   LUNARA(LUN) :',  LUNARA(LUN),
     &      '   LUNSTK(LUN):',   LUNSTK(LUN),  
     &      '   NSAMLUN(LUN) :', NSAMLUN(LUN),  
     &      '   NWORDLUN(LUN) :',NWORDLUN(LUN)  
         CALL ERRT(102,'INLINE BUFFER OVERFLOW ON LUN',LUN)
         IERR = 1
         RETURN

      ELSEIF (IGO < 0) THEN
C        BEFORE RESERVED BUFFER!
         WRITE(NOUT,"(A,I3,A,I6,A,I5,A,I5,A,I6,A,I6,A,I12)") 
     &      '*** INLN_WRTLIN ARRAY UNDERFLOW, LUN: ',LUN,
     &      '   IGO:',          IGO,
     &      '   NREC: ',        NREC,
     &      '   LUNARA(LUN) :', LUNARA(LUN),
     &      '   LUNSTK(LUN):',  LUNSTK(LUN),  
     &      '   NSAMLUN(LUN) :',NSAMLUN(LUN)  
         CALL ERRT(100,'WRTLIN',NE)
         IERR = 1
         RETURN
      ENDIF

C     WRITE THE BUFFER (THERE MAY BE A FASTER WAY TO DO THIS!!)
      DO I= 1,NB
         IPTNOW(IGO + I) = BUF(I) 
      ENDDO
      IERR = 0

      RETURN
      END




C++*********************************************************************
C   INLN_REDLIN.FOR                                 MAHIEDDINE LADJADJ
C **********************************************************************
C
C     INLN_REDLIN(LUN,BUF,NB,NREC)
C
C     PURPOSE:    READ A LINE OF FLOATING POINT NUMBERS FROM BUFFER
C
C     PARAMETERS:
C        LUN    LOGICAL UNIT NUMBER OF FILE BEING READ
C        BUF    BUFFER WHERE RECORD IS READ IN
C        NB     NUMBER OF VALUES IN RECORD TO BE READ
C        NREC   RECORD TO BE READ
C
C        IERR   ERROR CODE 1 IS RETURNED IN CASE OF ERROR
C               IERR IS DEFINED IN COMMON /IOERR/IERR
C 
C       MANY COMMANDS IN SPIDER READ A FILE, MANIPULATE THE DATA, AND
C       WRITE THE RESULTS INTO A FILE, ONLY TO HAVE THE NEXT COMMAND
C       READ THE DATA BACK AGAIN IN CORE TO MANIPULATE IT AND STORE
C       IT BACK TO A FILE THEN START THE CYCLE AGAIN...
C       TO SPEED UP THE PROCESS, THE USER HAS THE ABILITY TO HAVE THAT 
C       DATA KEPT IN INCORE MEMORY UNTIL HE DECIDES THAT HE WANTS IT IN
C       A FILE. FOR THAT, SPIDER OFFERS AN INCORE MEMORY IN A COMMON BLOCK
C       BUFFER:
C
C       THIS BUFFER CAN BE SUBDIVIDED INTO TWO PART TO HOLD TWO FILES, 
C       CALLED INFILES.
C       BY DEFAULT, THE INFILES ARE CALLED ___1 AND ___2. 
C
C       DATA FROM A FILE IS STORED EXACTLY THE SAME WAY AS IT IS IN 
C       DISK FILE (I.E: A HEADER FOLLOWED BY NROWS OF DATA OF LENGTH NSAMS
C       REPEATED NSLICE TIMES).
C 
C--*******************************************************************

      SUBROUTINE INLN_REDLIN(LUN,BUF,NB,NREC)

C     USE INLINE BUFFER COMMON AREA
      USE INLN_INFO

      REAL            BUF(NB)
      COMMON /IOERR/  IERR
      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      REAL,  POINTER :: IPTNOW(:)
      INTEGER * 8    :: IGO

      IPTNOW  => LOCLUN(LUN)%IPT

C     GET THE CORRECT OFFSET INTO THE BUFFER.
C     IGO = (LUNARA(LUN) + LUNSTK(LUN) + NREC - 1) * NSAMLUN(LUN)
      IGO = (LUNARA(LUN) + LUNSTK(LUN) + NREC - 1) 
      IGO = IGO * NSAMLUN(LUN)

      IF ((IGO + NB - 1) > NWORDLUN(LUN)) THEN
C        GOES OUT OF RESERVED BUFFER!
         WRITE(NOUT,"(A,I5,A,I12)") 
     &          '*** INLN_REDLIN READS BEYOND BUFFER, LUN: ',LUN,
     &          ' LENGTH:' ,IGO + NB
         IERR = 1
         RETURN
      ENDIF

C     READ THE BUFFER (THERE MAY BE A FASTER WAY TO DO THIS!!)
      DO I= 1,NB
         BUF(I) = IPTNOW(IGO + I)
      ENDDO
      IERR = 0

      RETURN
      END

C++*********************************************************************
C   INLN_WRTVOX.F           NEW                    SEPT 03 ARDEAN LEITH                     
C **********************************************************************
C
C  INLN_WRTVOX(LUN,NROW,BUF,NB,IX,IY,IZ,IRTFLG)
C
C   PURPOSE:    WRITE A FLOATING POINT NUMBERS INTO INCORE 'FILE'
C
C   PARAMETERS:
C        LUN      LOGICAL UNIT NUMBER OF FILE BEING WRITTEN
C        BUF      BUFFER WHERE VALUE IS KEPT
C        NB       NUMBER OF VALUES TO BE WRITTEN
C        IX,IY,IZ LOCATION
C 
C        IERR   ERROR CODE 1 IS RETURNED IN CASE OF ERROR
C 
C--*******************************************************************

      SUBROUTINE INLN_WRTVOX(LUN,NROW,BUF,NB,IX,IY,IZ,IRTFLG)

C     USE INLINE BUFFER COMMON AREA
      USE INLN_INFO

      REAL   BUF(NB)

      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      REAL, POINTER :: IPTNOW(:)
      INTEGER * 8   :: IGOM1

      IPTNOW => LOCLUN(LUN)%IPT

C     GET THE CORRECT OFFSET INTO THE BUFFER.
C     IGOM1 = (LUNARA(LUN) + LUNSTK(LUN) + (IZ - 1) * NROW + (IY-1)) * 
C    &       NSAMLUN(LUN) + IX - 1

      IGOM1 = (LUNARA(LUN) + LUNSTK(LUN) + (IZ - 1) * NROW + (IY-1)) 
      IGOM1 = IGOM1 * NSAMLUN(LUN) + IX - 1

      IF ((IGOM1 + NB) > NWORDLUN(LUN)) THEN
C        GOES OUT OF RESERVED BUFFER!
         WRITE(NOUT,"(A,3I6,A,I12,A,I12,A,I4)") 
     &      '*** INLN_WRTVOX OVERFLOW AT VOXEL: ', IX,IY,IZ,
     &      '   OVERFLOW LOCATION: ',              IGOM1+NB,  
     &      '   AVAILABLE LENGTH: ',               NWORDLUN(LUN)  
         CALL ERRT(102,'OVERFLOW IN INLN_WRTVOX, LUN',LUN)
         IERR = 1
         RETURN

      ELSEIF (IGOM1 < 0) THEN
C        BEFORE RESERVED BUFFER!
         WRITE(NOUT,"(A,3I6,A,I12,A,I12,A,I4)") 
     &      '*** INLN_WRTVOX UNDERFLOW AT VOXEL: ', IX,IY,IZ,
     &      '   OVERFLOW LOCATION: ',               IGOM1+NB,  
     &      '   AVAILABLE LENGTH: ',                NWORDLUN(LUN)  
         CALL ERRT(102,'UNDERFLOW IN INLN_WRTVOX, LUN',LUN)
         IERR = 1
         RETURN
      ENDIF

C     WRITE THE BUFFER (THERE MAY BE A FASTER WAY TO DO THIS!!)
      DO I= 1,NB
         IPTNOW(IGOM1 + I) = BUF(I) 
      ENDDO
      IERR = 0

      RETURN
      END




C++*********************************************************************
C   INLN_REDVOX.F           NEW                    SEPT 03 ARDEAN LEITH                     
C **********************************************************************
C
C     INLN_REDVOX(LUN,NROW,BUF,NB,IX,IY,IZ,IRTFLG)
C
C     PURPOSE:    READ A LINE OF FLOATING POINT NUMBERS FROM BUFFER
C
C     PARAMETERS:
C        LUN      LOGICAL UNIT NUMBER OF FILE BEING WRITTEN
C        BUF      BUFFER WHERE VALUE IS KEPT
C        NB       NUMBER OF VALUES TO BE WRITTEN
C        IX,IY,IZ LOCATION
C
C        IRTFLG   ERROR CODE 1 IS RETURNED IN CASE OF ERROR
C 
C--*******************************************************************

      SUBROUTINE INLN_REDVOX(LUN,NROW,BUF,NB,IX,IY,IZ,IRTFLG)

C     USE INLINE BUFFER COMMON AREA
      USE INLN_INFO

      REAL            BUF(NB)

      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      REAL, POINTER :: IPTNOW(:)
      INTEGER * 8   :: IGOM1

      IPTNOW  => LOCLUN(LUN)%IPT

C     GET THE CORRECT OFFSET INTO THE BUFFER.
C     IGOM1 = (LUNARA(LUN) + LUNSTK(LUN) + (IZ - 1) * NROW + (IY-1)) * 
C     &       NSAMLUN(LUN) + IX - 1

      IGOM1 = (LUNARA(LUN) + LUNSTK(LUN) + (IZ - 1) * NROW + (IY-1)) 
      IGOM1 = IGOM1 * NSAMLUN(LUN) + IX - 1 

      IF ((IGOM1 + NB) > NWORDLUN(LUN)) THEN
         WRITE(NOUT,"(A,3I6,A,I12,A,I12,A,I4)") 
     &      '*** INLN_REDVOX READS BEYOND BUFFER AT VOXEL:',IX,IY,IZ,
     &      '   LOCATION:',                             IGOM1+NB,  
     &      '   AVAILABLE LENGTH: ',                    NWORDLUN(LUN)  
         IRTFLG = 1
         RETURN
      ENDIF

C     READ THE BUFFER (THERE MAY BE A FASTER WAY TO DO THIS!!)
      DO I= 1,NB
         BUF(I) = IPTNOW(IGOM1 + I)
      ENDDO

C      CALL memcpy_(IPTNOW(IGOM1 + 1),BUF(1),NB*4) DOES NOT WORK!!

      IRTFLG = 0

      RETURN
      END



