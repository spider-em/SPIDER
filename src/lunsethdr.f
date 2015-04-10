
C++*********************************************************************
C
C LUNSETHDR.F   -- NEW JAN 1999                   AUTHOR: ARDEAN LEITH
C                  REPLACED ALLOCIT WITH ALLOCATE MAY 00  ARDEAN LEITH
C                  USED MYTIME                    DEC 00  ARDEAN LEITH
C                  KANGLE BUG FIXED               JAN 01  ARDEAN LEITH 
C                  GETFILENUM EXTRACTED           AUG 02  ARDEAN LEITH
C                  INDEXED STACKS                 JAN 03  ARDEAN LEITH 
C                  FORMATING IN LUNSAYINFO        MAY 04  ARDEAN LEITH
C                  HEADER BYTES MSG I7            FEB 05  ARDEAN LEITH
C                  HEADER INFO SIZE I6            JAN 06  ARDEAN LEITH
C                  LUNSETIMNUM REWRITE            JAN 07  ARDEAN LEITH
C                  GETLAB BAD VALUE TRAP          JAN 11  ARDEAN LEITH
C                  LABBYT >= 10000000) TRAP       MAR 11  ARDEAN LEITH
C                  MAXNAM                         JUL 14  ARDEAN LEITH
C                  EM2EM BUG RECOVERY SUPPORT     AUG 14  ARDEAN LEITH
C                    
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014 Health Research Inc.,                          *
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
C
C    PURPOSE:       HANDLES ALL INTERACTIONS WITH SPIDER IMAGE FILE
C                   HEADERS. CONTAINS NUMEROUS SUBROUTINES ALL STARTING
C                   WITH PREFIX: LUN
C
C    SUBROUTINES:
C     ------------------------- LUNSETHDR -----------------------------
C     ------------------------- LUNMTHDR ------------------------------
C     ------------------------- LUNREDHED -----------------------------
C     ------------------------- LUNWRTHED -----------------------------
C     ------------------------- LUNWRTCURHED --------------------------
C     ------------------------- LUNSETLAB -----------------------------
C     ------------------------- LUNGETLAB -----------------------------
C     ------------------------- LUNGETTYPE ----------------------------
C     ------------------------- LUNSETTYPE ----------------------------
C     ------------------------- LUNGETISTACK ----------------------------
C     ------------------------- LUNSETISTACK ----------------------------
C     ------------------------- LUNGETSTK -----------------------------
C     ------------------------- LUNSET25 ---------------------------
C     ------------------------- LUNSETINUSE ---------------------------
C     ------------------------- LUNGETINUSE ---------------------------
C     ------------------------- LUNGETMAXIM ---------------------------
C     ------------------------- LUNSETMAXIM ---------------------------
C     ------------------------- LUNSAVMAXIM ---------------------------
C     ------------------------- LUNCOPYMAXIM --------------------------
C     ------------------------- LUNGETSIZE ----------------------------
C     ------------------------- LUNSETSIZE ----------------------------
C     ------------------------- LUNGETSTAT ----------------------------
C     ------------------------- LUNSETSTAT ----------------------------
C     ------------------------- LUNGETANG -----------------------------
C     ------------------------- LUNSETANG -----------------------------
C     ------------------------- LUNGETFILE ----------------------------
C     ------------------------- LUNSETFILE ----------------------------
C     ------------------------- LUNSETIMNUM ---------------------------
C     ------------------------- LUNGETDATE ----------------------------
C     ------------------------- LUNGETTITLE ---------------------------
C     ------------------------- LUNSETTIME ----------------------------
C     ------------------------- LUNSETTITLE ---------------------------
C     ------------------------- LUNSAYINFO ----------------------------
C     ------------------------- LUNSETCOMMON --------------------------
C     ------------------------- LUNSETLUNS ----------------------------
C     ------------------------- LUNNEWHDR -----------------------------
C     ------------------------- LUNGETISBARE --------------------------
C     ------------------------- LUNSETISBARE --------------------------
C     ------------------------- LUNGETOBJ -----------------------------
C     ------------------------- LUNCLRINDX ----------------------------
C     ------------------------- LUNREDINDX ----------------------------
C     ------------------------- LUNWRTINDX ----------------------------
C     ------------------------- LUNGETINDXTOP -------------------------
C     ------------------------- LUNSETINDXTOP -------------------------
C     ------------------------- LUNSETVALS ----------------------------
C     ------------------------- LUNSETIMGOFF -------------------------
C      
C     STATIC LOCATION    257 -- IDSP
C                        258 -- ISBARE
C                        259 -- ISTACK (OVERALL)
C                        260 -- MAXIM
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************
 
      MODULE LUNHDR_INFO
         INTEGER, PARAMETER :: NUMLUNST = 100
         TYPE REAL_POINTER
            REAL, POINTER   :: IPT(:) 
         END TYPE REAL_POINTER

         TYPE(REAL_POINTER) :: LUNHDRBUF(NUMLUNST)
      END MODULE LUNHDR_INFO


C     ----------- LUNSETLUNS ---------------------------------------

      SUBROUTINE LUNSETLUNS(LUN,IVALA,IVALSTK,IVALLUN,IVALFLIP,IRTFLG)

C     SETS OFFSETS IN LUNSTK(LUN)....

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      LUNARA(LUN)  = IVALA
      LUNSTK(LUN)  = IVALSTK
      LUNARB(LUN)  = IVALLUN
      LUNFLIP(LUN) = IVALFLIP

      IRTFLG = 0

      END

C     ----------- LUNREDHED ---------------------------------------

      SUBROUTINE LUNREDHED(LUN,NX,IMGNUM,CALLERRT,IRTFLG)

C     READS IMAGE HEADER INTO HEADER OBJECT 
C     FOR STACKED IMAGES THIS MUST BE PRECEEDED BY A READ OF THE 
C     OVERALL HEADER TO ENSURE THAT LUNGETHEDOFF SUCCEEDS

#include "LUNHDR.INC"

      INCLUDE 'CMBLOCK.INC'

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      COMMON /IOERR/  IERR

      LOGICAL :: CALLERRT
      
C     SET PROPER  OFFSET IN LUNSTK
      LUNSTKSAV   = LUNSTK(LUN)
      LABRECSAV   = LUNARA(LUN)

C     FIND PROPER IMGNUM OFFSET 
      CALL LUNGETHEDOFF(LUN,NX,IMGNUM,
     &                  LUNARA(LUN),LUNSTK(LUN),IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     READ HEADER RECORDS FROM FILE INTO HEADER OBJECT
      IERR   = 0
      IRECT  = 1
      ILOC   = 1
      IRTFLG = 0

      DO WHILE (ILOC <= LENBUF)
C        LENT IS REMAINING LENGTH OF HEADER TO BE FILLED
         LENT = MIN(NX,LENBUF - ILOC + 1)

         CALL REDLIN(LUN,HEADER(ILOC),LENT,IRECT)
         IF (IERR .NE. 0) THEN
             IF (CALLERRT) THEN
                CALL ERRT(102,'I/O ERROR ON FILE HEADER',IERR)
             ENDIF
             IRTFLG = 1
             EXIT
          ENDIF
          ILOC  = ILOC + NX
          IRECT = IRECT + 1
      ENDDO

C     REPLACE OFFSETS IN LUNARA & LUNSTK    
      LUNARA(LUN) = LABRECSAV
      LUNSTK(LUN) = LUNSTKSAV 

      END

C     ------------------------- LUNGETHEDOFF -------------------------

      SUBROUTINE LUNGETHEDOFF(LUN,NX,IMGNUM,
     &           LUNARAOFF,LUNSTKOFF,IRTFLG)

C     SUPPORT ROUTINE TO RETURN RECORD OFFSET FOR HEADER REDLIN/WRTLIN

      INCLUDE 'CMBLOCK.INC'

      IF (IMGNUM == 0) THEN
C        OVERALL HEADER
         LUNARAOFF = 0
         LUNSTKOFF = 0
         IRTFLG    = 0
         RETURN
      ENDIF


C     GET RECORD INFO (CAN BE FROM OVERALL HEADER)
      CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,NDUM1,NDUM2,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (INDXREC > 0) THEN
C        INDEXED STACK IMAGE
         CALL LUNREDINDX(LUN,IMGNUM,INDX,NX,IRTFLG)
         IF (IRTFLG == 0 .AND. INDX <= 0) IRTFLG = -1
         IF (IRTFLG .NE. 0) RETURN

         LUNSTKOFF = (INDX-1) * NRECS + LABREC + INDXREC 
         LUNARAOFF = 0  

      ELSE

         CALL LUNGETSTKALL(LUN,ISTACK,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (ISTACK == 0 .AND. IMGNUM > 1) THEN
C           NORMAL IMAGE, IMGNUM CAN NOT BE > 1
            CALL ERRT(102,'NOT A STACK, NO IMAGE',IMGNUM)
            RETURN

         ELSEIF (ISTACK == 0) THEN
C           NORMAL IMAGE, IMGNUM IS 1
            LUNSTKOFF = 0
            LUNARAOFF = 0 

         ELSE
C           NORMAL STACKS HAVE ADDITIONAL OVERALL HEADER AT BEGINNING

            CALL LUNGETMAXIM(LUN,MAXIM,IRTFLG)
            IF (IMGNUM > MAXIM) THEN
               IRTFLG = 1
               RETURN
            ENDIF

            LUNSTKOFF = (IMGNUM-1) * NRECS + LABREC
            LUNARAOFF = 0 

         ENDIF
      ENDIF

      END


C     ----------- LUNNEWHDR -----------------------------------------

      SUBROUTINE LUNNEWHDR(LUN,IRTFLG)

C     CREATES STORAGE SPACE FOR A HEADER OBJECT

#include "LUNHDR.INC"
      INCLUDE 'CMBLOCK.INC'

      IPOINTER => LUNHDRBUF(LUN)%IPT 
      IF (.NOT. ASSOCIATED(IPOINTER)) THEN
C        ALLOCATE SPACE FOR THIS HEADER OBJECT
         ALLOCATE(IPOINTER(LENHDR),STAT=IRTFLG)

         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(46,'LUNNEWHDR, FILE HEADER',LENHDR)
            IRTFLG = 1
            RETURN
         ENDIF
         LUNHDRBUF(LUN)%IPT => IPOINTER
      ENDIF

      IRTFLG = 0

      END


C     ------------------------- LUNGETOBJ -------------------------

      SUBROUTINE LUNGETOBJ(LUN,IPOINTER,IRTFLG)

      USE LUNHDR_INFO

      REAL, DIMENSION(:), POINTER :: IPOINTER 

      PARAMETER        (NUMLUNS = 100)

C     POINT TO HEADER OBJECT
      IRTFLG   = 1

      IF (LUN <= 0 .OR. LUN > NUMLUNS) THEN
         CALL ERRT(102,'PGM ERROR, LUN OUT OF RANGE', LUN)
         RETURN
      ENDIF

      IPOINTER => LUNHDRBUF(LUN)%IPT
      IF (.NOT. ASSOCIATED(IPOINTER)) RETURN

      IRTFLG = 0

      END

C     ------------------------- LUNSETIMGOFF -------------------------

      SUBROUTINE LUNSETIMGOFF(LUN,IMGNUM,NX,IRTFLG)

C     SUPPORT ROUTINE TO SET RECORD OFFSET FOR REDLIN/WRTLIN
C     SHOULD ONLY BE CALLED FOR STACKS!

      INCLUDE 'CMBLOCK.INC'

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

C     GET RECORD INFO (CAN BE FROM OVERALL HEADER)
      CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,NDUM1,NDUM2,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (IMGNUM == 0) THEN
C        OVERALL HEADER
         LUNARAOFF = 0
         LUNSTKOFF = 0

      ELSEIF (INDXREC > 0) THEN
C        INDEXED STACK
         CALL LUNREDINDX(LUN,IMGNUM,INDX,NX,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         LUNSTKOFF = (INDX-1) * NRECS + INDXREC + LABREC 
         LUNARAOFF = LABREC
  
      ELSE
C        NORMAL IMAGE OR NORMAL STACK
         CALL LUNGETSTKALL(LUN,ISTACK,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (ISTACK == 0 .AND. IMGNUM > 1) THEN
C           NORMAL IMAGE, IMGNUM CAN NOT BE > 1
            CALL ERRT(101,'NOT A STACK, NO IMAGE',IMGNUM)
            RETURN

         ELSEIF (ISTACK == 0) THEN
C           NORMAL IMAGE, IMGNUM IS 1
            LUNSTKOFF = 0
            LUNARAOFF = LABREC 

         ELSE
C           NORMAL STACKS HAVE ADDITIONAL OVERALL HEADER AT BEGINNING

            LUNSTKOFF = (IMGNUM-1) * NRECS + LABREC 
            LUNARAOFF = LABREC 
         ENDIF
      ENDIF

      LUNARA(LUN) = LUNARAOFF
      LUNSTK(LUN) = LUNSTKOFF

      END

C     ------------------------- LUNGETTYPE ----------------------------

      SUBROUTINE LUNGETTYPE(LUN,ITYPE,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      ITYPE  = HEADER(5)

      END

C     ------------------------- LUNSETTYPE ----------------------------

      SUBROUTINE LUNSETTYPE(LUN,ITYPE,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET VALUES
      HEADER(5)  = ITYPE

      END

C     ------------------------- LUNGETISTACK ----------------------------

      SUBROUTINE LUNGETISTACK(LUN,ISTACK,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      ISTACK  = HEADER(24)

      END
C     ------------------------- LUNSETISTACK ----------------------------

      SUBROUTINE LUNSETISTACK(LUN,ISTACK,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET VALUES
      HEADER(24) = ISTACK

      END

C     ------------------------- LUNGETLAB ----------------------------

      SUBROUTINE LUNGETLAB(LUN,LABREC,INDXREC,NRECS,
     &                     LABBYT,LENBYT,IRTFLG)

#include "LUNHDR.INC"
      INCLUDE 'CMBLOCK.INC'

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      LABREC  = HEADER(13)
      NRECS   = HEADER(3)
      LABBYT  = HEADER(22)
      LENBYT  = HEADER(23)

      IF (LABREC <= 0 .OR. LABREC >= 100000 .OR.
     &    LABBYT <= 0 .OR. LABBYT >= 10000000) THEN
C        IMPOSSIBLE VALUE --> NON-EXISTING IMAGE??
         IRTFLG = 1
         RETURN
      ENDIF

      IF (LENBYT <= 0) THEN
C        CORRECT BAD LENBYT (UNAVAILABLE ON OLD VAX FILES)
         CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
         LENBYT = NX * 4
      ENDIF

C     CORRECT UNREASONABLE LABREC (BAD VALUE ONCE)        
      LABRECT = 1024 / LENBYT

      ITEMP = MOD(1024,LENBYT)
      IF (ITEMP .NE. 0) LABRECT = LABRECT + 1

      IF (LABRECT <= 0 .OR. LABREC .NE. LABRECT) THEN
C         UNREASONABLE LABREC NUMBER SO DEFAULT IT
          LABREC = LABRECT
      ENDIF

C     FIND  NUMBER OF INDX RECORDS IN INDEXED STACK HEADER
      ISTACK  = HEADER(259)
      IF (ISTACK >= 0) THEN 
C        NOT AN INDEXED STACK
         INDXREC  = 0
      ELSE
         ISTACK = - ISTACK
         FINDXREC = FLOAT(ISTACK) / FLOAT(LENBYT / 4)
         INDXREC  = ISTACK /  (LENBYT / 4)
         IF (FINDXREC > INDXREC) INDXREC  = INDXREC + 1 
      ENDIF

      IRTFLG = 0

      END

C     ------------------------- LUNSETLAB ----------------------------

      SUBROUTINE LUNSETLAB(LUN,LABREC,NRECS,LABBYT,LENBYT,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      HEADER(13) = LABREC  
      HEADER(3)  = NRECS   
      HEADER(22) = LABBYT  
      HEADER(23) = LENBYT  

      IRTFLG = 0

      END

C     ----------- LUNWRTCURHED -----------------------------------------

      SUBROUTINE LUNWRTCURHED(LUN,IRTFLG)

      INTEGER  :: LUN,IRTFLG

      INTEGER  :: NX,NY,NZ,IMGNUM

C     GET IMGNUM VALUE 
      CALL LUNGETINUSE(LUN,IMGNUM,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET NX VALUE 
      CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     REPLACE THE CURRENT HEADER BACK IN THE FILE
      CALL LUNWRTHED(LUN,NX,IMGNUM,IRTFLG)

      END

C     ----------- LUNWRTHED -----------------------------------------

      SUBROUTINE LUNWRTHED(LUN,NX,IMGNUM,IRTFLG)

C     WRITES HEADER OBJECT TO SPECIFIED IMAGE HEADER

#include "LUNHDR.INC"

      INTEGER  :: LUNARA,LUNSTK,LUNARB
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      INCLUDE 'CMBLOCK.INC'

      INTEGER  :: IERR
      COMMON /IOERR/ IERR

      INTEGER  :: LUN,NX,IMGNUM,IRTFLG

      INTEGER  :: LUNSTKSAV,LABRECSAV,IRECT,ILOC,LENT


C     SET PROPER IMGNUM OFFSET 

C     SAVE CURRENT FILE OFFSETS FOR LUNARA & LUNSTK
      LUNSTKSAV = LUNSTK(LUN)
      LABRECSAV = LUNARA(LUN)

      CALL LUNGETHEDOFF(LUN,NX,IMGNUM,
     &     LUNARA(LUN),LUNSTK(LUN),IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     WRITE HEADER RECORDS FROM LUNHDR INTO FILE 
      IERR   = 0
      IRECT  = 1
      ILOC   = 1
      IRTFLG = 0

      DO WHILE (ILOC <= LENBUF)
C        LENT IS REMAINING LENGTH OF HEADER TO BE WRITTEN
         LENT = MIN(NX,LENBUF - ILOC + 1)

         CALL WRTLIN(LUN,HEADER(ILOC),LENT,IRECT)
         IF (IERR .NE. 0) THEN
            CALL ERRT(102,'WRITING TO FILE HEADER, RECORD #',IRECT)
            IRTFLG = 1
            GOTO 999
         ENDIF

         ILOC  = ILOC + NX
         IRECT = IRECT + 1
      ENDDO

C     REPLACE OFFSETS IN LUNARA & LUNSTK    
999   LUNARA(LUN) = LABRECSAV
      LUNSTK(LUN) = LUNSTKSAV 

      END


C     ----------- LUNMTHDR -----------------------------------------
 
      SUBROUTINE LUNMTHDR(LUN,IRTFLG)

#include "LUNHDR.INC"

C     WANT TO DEALLOCATE AND CLOSE AN OPEN HEADER OBJECT

      IPOINTER => LUNHDRBUF(LUN)%IPT
      IF (ASSOCIATED(IPOINTER))  DEALLOCATE(IPOINTER)
      NULLIFY(LUNHDRBUF(LUN)%IPT)

C     CLEAR FILENAME ALSO
      LUNFILNAM(LUN) = CHAR(0)
      
      IRTFLG = 0

      END

C     ------------------------- LUNGETSTK ----------------------------

       SUBROUTINE LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      ISTACK = HEADER(24)
      MAXIM  = HEADER(26)

      IRTFLG = 0

      END

C     ------------------------- LUNGETINUSE ----------------------------

       SUBROUTINE LUNGETINUSE(LUN,INUSE,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
 
C     SET STACK RELATED LOCATIONS
      INUSE = HEADER(27) 

      END

C     ------------------------- LUNSETINUSE ----------------------------

       SUBROUTINE LUNSETINUSE(LUN,INUSE,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET STACK RELATED LOCATIONS
      HEADER(27) = INUSE 

      IRTFLG = 0

      END

C     ------------------------- LUNGETMAXIM ----------------------------

       SUBROUTINE LUNGETMAXIM(LUN,MAXIM,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET MAXIM VALUE FROM HEADER OBJECT STATIC AREA
      MAXIM   = HEADER(260)

      IRTFLG = 0

      END

C     ------------------------- LUNCOPYMAXIM -------------------------

       SUBROUTINE LUNCOPYMAXIM(LUN,MAXIM,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET MAXIM VALUE FROM HEADER OBJECT NON-STATIC AREA
      MAXIM       = HEADER(26)
      HEADER(260) = HEADER(26)

      IRTFLG = 0

      END

C     ------------------------- LUNSETMAXIM ----------------------------

       SUBROUTINE LUNSETMAXIM(LUN,MAXIM,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET MAXIM VALUE IN HEADER OBJECT 
      HEADER(26) = MAXIM

      IRTFLG = 0

      END


C     ------------------------- LUNSETMAXALL -------------------------

       SUBROUTINE LUNSETMAXALL(LUN,MAXIM,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET MAXIM VALUE IN HEADER OBJECT STATIC AREA
      HEADER(260)  = MAXIM

      IRTFLG = 0

      END


C     ------------------------- LUNSAVMAXIM ----------------------------

      SUBROUTINE LUNSAVMAXIM(LUN,NX,MAXIM,IRTFLG)

C     COMPLICATED SINCE I DO NOT WANT TO ALTER CURRENT FILE HEADER OBJECT

#include "LUNHDR.INC"

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
      COMMON /IOBUF/ BUF(NBUFSIZ)

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      COMMON /IOERR/  IERR

C     SET PROPER IMGNUM OFFSET IN LUNSTK
      IMGOFFSET   = 0
      LUNSTKSAV   = LUNSTK(LUN)
      LUNSTK(LUN) = IMGOFFSET

C     SET NO OFFSET FOR HEADER RECORDS IN LUNARA
      LABRECSAV   = LUNARA(LUN)
      LUNARA(LUN) = 0

C     READ HEADER RECORDS FROM FILE INTO BUFFER
      IERR  = 0
      IRECT = 1
      ILOC  = 1
      DO WHILE (ILOC <= LENBUF)
C        LENT IS REMAINING LENGTH OF HEADER TO BE FILLED
         LENT = MIN(NX,LENBUF - ILOC + 1)

         CALL REDLIN(LUN,BUF(ILOC),LENT,IRECT)
         IF (IERR .NE. 0) THEN
            CALL ERRT(102,'READING OVERALL HEADER, I/O ERROR',IERR)
            IRTFLG = 1
            GOTO 999
         ENDIF
         ILOC  = ILOC + NX
         IRECT = IRECT + 1
      ENDDO

C     SET MAXIM VALUE IN BUF
      BUF(26)  = MAXIM
      BUF(260) = MAXIM

C     WRITE HEADER RECORDS BACK INTO FILE 
      IERR  = 0
      IRECT = 1
      ILOC  = 1
      DO WHILE (ILOC <= LENBUF)
C        LENT IS REMAINING LENGTH OF HEADER TO BE WRITTEN
         LENT = MIN(NX,LENBUF - ILOC + 1)

         CALL WRTLIN(LUN,BUF(ILOC),LENT,IRECT)
         IF (IERR .NE. 0) THEN
            CALL ERRT(102,'WRITING OVERALL HEADER, I/O ERROR',IERR)
            IRTFLG = 1
            GOTO 999
         ENDIF

         ILOC  = ILOC + NX
         IRECT = IRECT + 1
      ENDDO

      IRTFLG = 0

C     REPLACE OFFSETS IN LUNARA & LUNSTK    
999   LUNARA(LUN) = LABRECSAV
      LUNSTK(LUN) = LUNSTKSAV 

      END


C     ------------------------- LUNGETSIZE ----------------------------

       SUBROUTINE LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      NX = HEADER(12)
      NY = HEADER(2)
      NZ = HEADER(1)
      IF (NZ < 0) NZ = -NZ

      IRTFLG = 0

      END

C     ------------------------- LUNSETSIZE ----------------------------

       SUBROUTINE LUNSETSIZE(LUN,NX,NY,NZ,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET VALUES
      HEADER(12)  = NX
      HEADER(2)   = NY
      HEADER(1)   = NZ

      IRTFLG = 0

      END



C     ------------------------- LUNGETSTAT ----------------------------

      SUBROUTINE LUNGETSTAT(LUN,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      IMAMI = HEADER(6) + 0.5
      FMIN  = HEADER(8)
      FMAX  = HEADER(7)
      AV    = HEADER(9)
      SIG   = HEADER(10)

      IRTFLG = 0

      END
C     ------------------------- LUNSETSTAT ----------------------------

      SUBROUTINE LUNSETSTAT(LUN,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      HEADER(6)  = IMAMI 
      HEADER(8)  = FMIN  
      HEADER(7)  = FMAX  
      HEADER(9)  = AV    
      HEADER(10) = SIG  

      IRTFLG = 0

      END

C     ------------------------- LUNGETANG ----------------------------

      SUBROUTINE LUNGETANG(LUN,IANGLE,PHI,THETA,PSI,XOFF,YOFF,ZOFF,
     &                     KANGLE,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RETURN VALUES
      IANGLE   = HEADER(14)
      PHI      = HEADER(15)
      THETA    = HEADER(16)
      PSI      = HEADER(17)
      XOFF     = HEADER(18)
      YOFF     = HEADER(19)
      ZOFF     = HEADER(20)
      KANGLE   = HEADER(30)

      IRTFLG = 0

      END

C     ------------------------- LUNGETVAL ----------------------------

       SUBROUTINE LUNGETVAL(LUN,ILOC,VAL,IRTFLG)

#include "LUNHDR.INC"

      INTEGER   :: LUN,ILOC,IRTFLG
      REAL      :: VAL

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET STACK RELATED LOCATIONS
      VAL = HEADER(25)

      IRTFLG = 0

      END


C     ------------------------- LUNGETVALS ----------------------------

      SUBROUTINE LUNGETVALS(LUN,IGO,NVAL,BUFOUT,IRTFLG)

#include "LUNHDR.INC"

      REAL, DIMENSION(NVAL)  :: BUFOUT

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IEND = IGO + NVAL - 1
      IF (IGO <= 0 .OR. IEND > 256) THEN
         CALL ERRT(102,'HEADER LOCATION MUST BE 0...256',IEND)
         IRTFLG = 1
         RETURN
      ENDIF

C     GET RETURN VALUES
      DO IVAL = IGO,IEND
         BUFOUT(IVAL-IGO+1) = HEADER(IVAL)
      ENDDO

      IRTFLG = 0

      END

C     ------------------------- LUNSETVALS ----------------------------

      SUBROUTINE LUNSETVALS(LUN,IGO,NVAL,BUFVALS,IRTFLG)

#include "LUNHDR.INC"

      REAL, DIMENSION(NVAL)  :: BUFVALS

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
  
      IEND = IGO+NVAL-1
      IF (IGO <= 0 .OR. IEND > 256) THEN
         CALL ERRT(102,'HEADER LOCATION MUST BE < 257',IEND)
         IRTFLG = 1
         RETURN
      ENDIF

C     SET VALUES IN HEADER OBJECT
      DO IVAL = IGO,IEND
         HEADER(IVAL) = BUFVALS(IVAL-IGO+1)   
      ENDDO

C     COPY HEADER OBJECT TO FILE
      CALL LUNWRTCURHED(LUN,IRTFLG)

      END

C     ------------------------- LUNSETANG ----------------------------

      SUBROUTINE LUNSETANG(LUN,IANGLE,PHI,THETA,PSI,XOFF,YOFF,ZOFF,
     &                     KANGLE,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET VALUES
      HEADER(14) = IANGLE
      HEADER(15) = PHI
      HEADER(16) = THETA
      HEADER(17) = PSI
      HEADER(18) = XOFF
      HEADER(19) = YOFF
      HEADER(20) = ZOFF
      HEADER(30) = KANGLE

C     COPY HEADER OBJECT TO FILE
      CALL LUNWRTCURHED(LUN,IRTFLG)

      END

C     ------------------------- LUNGETFILE ----------------------------

      SUBROUTINE LUNGETFILE(LUN,FILNAM,NLET,DSP,IRTFLG)

#include "LUNHDR.INC"

      CHARACTER (LEN=*) :: FILNAM
      CHARACTER (LEN=1) :: DSP

C     RETRIEVE CURRENT FILENAME
      FILNAM = LUNFILNAM(LUN)
      NLET   = LNBLNKN(FILNAM)

      CALL LUNGETDSP(LUN,DSP,IRTFLG)

      IRTFLG = 0

      END

C     ------------------------- LUNSETFILE ----------------------------

      SUBROUTINE LUNSETFILE(LUN,FILNAM,DSP,IRTFLG)


#include "LUNHDR.INC"

      CHARACTER *(*)   FILNAM
      CHARACTER *1     DSP

      IF (FILNAM(1:1) .NE. CHAR(0)) THEN
C        SET CURRENT FILENAME IN HEADER OBJECT 
         NLET           = LNBLNKN(FILNAM)
         LUNFILNAM(LUN) = FILNAM(1:NLET)
      ENDIF

      IF (DSP .NE. CHAR(0)) THEN
         CALL LUNSETDSP(LUN,DSP,IRTFLG)
      ENDIF

      IRTFLG = 0

      END

C     ------------------------- LUNSETIMNUM -------------------------

      SUBROUTINE LUNSETIMNUM(LUN,FILNAM,IMGNUM,DSP,IRTFLG)

      CHARACTER(LEN=*)       ::  FILNAM
      CHARACTER(LEN=1)       ::  DSP

C     APPENDS IMGNUM TO INPUT: FILNAM  AFTER @ OR 
C         SETS FILENAME IN: HEADER OBJECT
C         ALSO RETURNS: NEW FILENAME IN: FILNAM

C     APPEND IMAGE NUMBER TO BARE STACK FILE NAME
C     (INTTOCHAR ALSO RETURNS NEW VALUE FOR NLET)
 
      LENAT = INDEX(FILNAM,'@')
      IF (LENAT == 0) LENAT = INDEX(FILNAM,'*') -1
      IF (LENAT < 0) LENAT = 0

      CALL INTTOCHAR(IMGNUM,FILNAM(LENAT+1:),NLET,0)
      IF (NLET < 0) THEN
         IRTFLG = 1
         RETURN
      ENDIF
      NLET = NLET + LENAT

C     SET NEW FILENAME IN HEADER OBJECT
      CALL LUNSETFILE(LUN,FILNAM(1:NLET),'N',IRTFLG)

      END

C     ------------------------- LUNGETINDXTOP -------------------------

      SUBROUTINE LUNGETINDXTOP(LUN,INDXTOP,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT (MUST BE OVERALL HEADER)
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET INDXTOP FROM HEADER OBJECT
      INDXTOP = HEADER(28)
      IRTFLG = 0

      END

C     ------------------------- LUNSETINDXTOP -------------------------

      SUBROUTINE LUNSETINDXTOP(LUN,INDXTOP,IRTFLG)

C     SETS INDXTOP IN HEADER OBJECT

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT (MUST BE OVERALL HEADER!!)
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET INDXTOP IN HEADER OBJECT
      HEADER(28) = INDXTOP
      IRTFLG = 0

      END

C     ------------------------- LUNGETDSP -------------------------

      SUBROUTINE LUNGETDSP(LUN,DSP,IRTFLG)

#include "LUNHDR.INC"

      CHARACTER * 1    DSP

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET DSP FROM HEADER OBJECT
      IDSP = HEADER(257)

      DSP = 'O'
      IF (IDSP == 1) DSP = 'N'

      END

C     ------------------------- LUNSETDSP -------------------------

      SUBROUTINE LUNSETDSP(LUN,DSP,IRTFLG)

#include "LUNHDR.INC"

      CHARACTER * 1    DSP

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET DSP IN HEADER OBJECT
      HEADER(257) = 0
      IF ( DSP == 'N') HEADER(257) = 1

      END

C     ------------------------- LUNGETISBARE -------------------------

      SUBROUTINE LUNGETISBARE(LUN,ISBARE,IRTFLG)

#include "LUNHDR.INC"

      LOGICAL    ISBARE

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET ISBARE FROM HEADER OBJECT
      ISBARE = .FALSE.
      IF (HEADER(258) == 1.0) ISBARE = .TRUE.

      END


C     ------------------------- LUNSETISBARE -------------------------

      SUBROUTINE LUNSETISBARE(LUN,ISBARE,IRTFLG)

#include "LUNHDR.INC"

      LOGICAL          ISBARE

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET ISBARE IN HEADER OBJECT
      HEADER(258)             = 0.0
      IF (ISBARE) HEADER(258) = 1.0

      END


C     ------------------------- LUNBAREFILE -------------------------

      SUBROUTINE LUNBAREFILE(LUN,FILNAM,IRTFLG)

      CHARACTER *(*)   FILNAM
      LOGICAL          ISBARE

      LOCAT  = INDEX(FILNAM,'@')
      NLET   = LNBLNKN(FILNAM)
      ISBARE = (LOCAT > 0 .AND. LOCAT == NLET) 

C     SET ISBARE IN HEADER OBJECT
      CALL LUNSETISBARE(LUN,ISBARE,IRTFLG)

      END

C     ------------------------- LUNSETSTKALL -------------------------

      SUBROUTINE LUNSETSTKALL(LUN,ISTACK,IRTFLG)

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET ISTACK IN STATIC AREA OF HEADER OBJECT
      HEADER(259) = ISTACK

      RETURN
      END

C     ------------------------- LUNCOPYSTK -------------------------

      SUBROUTINE LUNCOPYSTK(LUN,ISTACK,IRTFLG)

#include "LUNHDR.INC"

      LOGICAL          ISBARE

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET ISTACK IN HEADER OBJECT STATIC AREA
      ISTACK      = HEADER(24) 
      HEADER(259) = ISTACK

      RETURN
      END

C     ------------------------- LUNGETSTKALL -------------------------

      SUBROUTINE LUNGETSTKALL(LUN,IVAL,IRTFLG)

#include "LUNHDR.INC"

      LOGICAL    ISBARE

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET IVAL FROM HEADER OBJECT
      IVAL = HEADER(259) 

      RETURN
      END

C     ------------------------- LUNGETTITLE ----------------------------

      SUBROUTINE LUNGETTITLE(LUN,FILETITLE,NLET,IRTFLG)

#include "LUNHDR.INC"

      CHARACTER(LEN=*)   :: FILETITLE

      CHARACTER(LEN=180) :: CLINE
      REAL,DIMENSION(45) :: ZBUF
      EQUIVALENCE    (CLINE,ZBUF)

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     MOVE ALPHA-NUMERIC PART OF HEADER TO ZBUF FOR CLINE ACCESS
      DO  I = 1,45 
         ZBUF(I) = HEADER(I+211) 
      ENDDO

      IF (CLINE(2:2) == '-') THEN
C        NOTE THAT ALPHA-NUMERICAL DATA (ABCD) WILL BE WRITTEN (DCBA)
         CALL REVERSEBYTES(CLINE,180,IRTFLG)
      ENDIF

C     RECOVER TITLE FROM CLINE
      FILETITLE = CLINE(21:180)
      NLET      = lnblnkn(FILETITLE)

      IRTFLG = 0

      END

C     ------------------------- LUNGETDATE ----------------------------

      SUBROUTINE LUNGETDATE(LUN,FILEDATE,FILETIME,IRTFLG)

#include "LUNHDR.INC"

      CHARACTER *(*) FILEDATE,FILETIME
      DIMENSION      ZBUF(45)
      CHARACTER *180 CLINE
      EQUIVALENCE    (CLINE,ZBUF)

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     MOVE ALPHA-NUMERIC PART OF HEADER TO ZBUF FOR CLINE ACCESS
      DO  I = 1,7 
         ZBUF(I) = HEADER(I+211) 
      ENDDO

      IF (CLINE(2:2) == '-') THEN
C        NOTE THAT ALPHA-NUMERICAL DATA (ABCD) WILL BE WRITTEN (DCBA)
         CALL REVERSEBYTES(CLINE,28,IRTFLG)
      ENDIF

C     READ DATE & TIME FROM HEADER OBJECT
      FILEDATE  = CLINE(1:11) // ' '
      IF (FILEDATE(10:10) == CHAR(0) .OR. FILEDATE(10:10) == ' ')
     &   THEN
C        2 DIGIT DATE, MAKE IT 4 DIGIT DATE
         FILEDATE(10:11) = FILEDATE(8:9)
         FILEDATE(8:9)   = '19'
      ENDIF
      FILETIME  = CLINE(13:20) 

      IRTFLG = 0

      END

C     ------------------------- LUNSETTIME ----------------------------

      SUBROUTINE LUNSETTIME(LUN,IRTFLG)

#include "LUNHDR.INC"

      DIMENSION      ZBUF(5)
      CHARACTER *20  CLINE
      EQUIVALENCE    (CLINE,ZBUF)

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     PUT CURRENT DATE AND TIME INTO THIS IMAGE HEADER
C     Y2K DATE TAKES 2 & 3/4 FLOATING POINT VARIABLES IN BUF (11 CHAR)
      CALL DATE_2K(CLINE)

C     PUT CURRENT TIME INTO THIS IMAGE HEADER
C     TIME TAKES 2 FLOATING POINT VARIABLES IN BUF (8 CHAR.)
      CALL MYTIME(CLINE(13:20))

C     COPY CLINE STUFF INTO HEADER OBJECT
      DO  I = 1,5 
          HEADER(I+211) = ZBUF(I)
      ENDDO

      END     
C     ------------------------- LUNSETTITLE ----------------------------

      SUBROUTINE LUNSETTITLE(LUN,FILETITLE,IRTFLG)

#include "LUNHDR.INC"
      CHARACTER *(*) FILETITLE
      DIMENSION      ZBUF(40)
      CHARACTER *160 CLINE
      EQUIVALENCE    (CLINE,ZBUF)

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (FILETITLE(1:1) .NE. CHAR(0)) THEN
C        TITLE TAKES 40 FLOATING POINT VARIABLES IN BUF (160 CHAR)
         CLINE(1:160) = FILETITLE(1:160)

C        COPY CLINE STUFF  INTO HEADER OBJECT
         DO  I = 1,40 
            HEADER(I+216) = ZBUF(I)
         ENDDO
      ENDIF

      IRTFLG = 0

      END     

C     ------------------------- LUNSAYINFO ----------------------------

      SUBROUTINE LUNSAYINFO(LUN,IRTFLG)

      INCLUDE 'CMLIMIT.INC'

      CHARACTER *1     DSP
      CHARACTER *104   CSTRING
      CHARACTER * 2    TYPE
      CHARACTER *12    CDAT
      CHARACTER *8     CTIM
      CHARACTER *160   CTIT

      CHARACTER(LEN=MAXNAM):: FILNAM

      LOGICAL       :: SILENT,VERBOSE,USE_SPIRE
      COMMON /IPRTT/   IDUM245,NTRACE,NALPH,VERBOSE,USE_SPIRE,SILENT
      COMMON /UNITS/   LUNC,NIN,NOUT,NECHO,IFOUND,NLOG,NDAT

#ifdef USE_MPI
      include 'mpif.h' 
      INTEGER MYPID, COMM, IERR
      COMM = MPI_COMM_WORLD 
      CALL MPI_COMM_RANK(COMM, MYPID, IERR) 
#else
      MYPID = -1
#endif 
      IRTFLG = 1

C     RETRIEVE ITYPE
      CALL LUNGETTYPE(LUN,ITYPE,IRTFLG)

C     RETRIEVE CURRENT HEADER ISTACK 
      CALL LUNGETISTACK(LUN,ISTACK,IRTFLG)

C     RETRIEVE SIZE
      CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)

C     RETRIEVE IMGNUM FROM HEADER OBJECT
      CALL LUNGETINUSE(LUN,IMGNUM,IRTFLG)

      IF (ITYPE == -2) THEN
            TYPE = 'P '
      ELSEIF (ITYPE == -9)  THEN
            TYPE = 'FS'
      ELSEIF (ITYPE == -11) THEN
            TYPE = 'O2'
      ELSEIF (ITYPE == -12) THEN
            TYPE = 'E2'
      ELSEIF (ITYPE == -21) THEN
            TYPE = 'O3'
      ELSEIF (ITYPE == -22) THEN
            TYPE = 'E3'
      ELSEIF (ITYPE == 0)   THEN
            TYPE = 'D '
      ELSEIF (ITYPE == 1 .AND. ISTACK < 0 .AND. IMGNUM <= 0)THEN
            TYPE = 'I2'
      ELSEIF (ITYPE == 3 .AND. ISTACK < 0 .AND. IMGNUM <= 0)THEN
            TYPE = 'I3'
      ELSEIF (ITYPE == 1 .AND. ISTACK > 0 .AND. IMGNUM <= 0)THEN
            TYPE = 'S2'
      ELSEIF (ITYPE == 3 .AND. ISTACK > 0 .AND. IMGNUM <= 0)THEN
            TYPE = 'S3'
      ELSEIF (ITYPE == 3)  THEN
            TYPE = 'R3'
      ELSE
            TYPE = 'R '
      ENDIF

C     RETRIEVE CURRENT FILENAME (INCLUDES IMGNUM IF STACKED IMAGE)
      CALL LUNGETFILE(LUN,FILNAM,NLET,DSP,IRTFLG)

C     RECOVER FILE DATE, TIME & TITLE FROM HEADER
      CALL LUNGETDATE(LUN,CDAT,CTIM,IRTFLG)
      CALL LUNGETTITLE(LUN,CTIT,LENTIT,IRTFLG)
      
C     GET RECORD INFO
      CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,LABBYT,LENBYT,IRTFLG)

      IF (USE_SPIRE .AND. DSP == 'N' .AND. FILNAM(1:1) .NE. '_') THEN
         CALL SPIREOUT(FILNAM(:NLET),IRTFLG)

         IF (ISTACK .NE. 0 .AND. IMGNUM == 0 .AND. NZ > 1)THEN
            CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG)
            WRITE(CSTRING,89)TYPE,NX,NY,NZ,MAXIM,CDAT,CTIM,
     &                       DSP,LABBYT
         ELSEIF (ISTACK .NE. 0 .AND. IMGNUM == 0)THEN
C           OVERALL STACKED IMAGE 
            CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG)
            WRITE(CSTRING,90)TYPE,NX,NY,MAXIM,CDAT,CTIM,DSP,LABBYT

         ELSE IF (IMGNUM > 0) THEN
C           STACKED IMAGE
            WRITE(NOUT,92)TYPE,NX,NY,IMGNUM,CDAT,CTIM,DSP
            WRITE(CSTRING,92)TYPE,NX,NY,IMGNUM,CDAT,CTIM,DSP

        ELSE IF (NZ > 1) THEN
C           SIMPLE VOLUME
            WRITE(CSTRING,93)TYPE,NX,NY,NZ,CDAT,CTIM,DSP,LABBYT

        ELSE IF (NZ > 1) THEN
C           SIMPLE VOLUME
            WRITE(CSTRING,93)TYPE,NX,NY,NZ,CDAT,CTIM,DSP,LABBYT

         ELSE
C           SIMPLE IMAGE
            WRITE(CSTRING,94)TYPE,NX,NY,CDAT,CTIM,DSP,LABBYT
         ENDIF

         CALL SPIREOUT(CSTRING,IRTFLG)
      ENDIF



      IF (VERBOSE .AND. IFOUND .NE. -4) THEN
C         PRINT STACK OPENING INFORMATION

         LENT = LENTIT + NLET
         IF (LENTIT <= 0 .AND. NLET > 0) THEN
C           FILENAME BUT NO TITLE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,*) ' ',FILNAM(:NLET)
            ENDIF


         ELSE IF (LENT > 0 .AND. LENT < 70) THEN
C           HAS FILENAME AND TITLE THAT FIT ON ONE LINE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,*) ' ',FILNAM(:NLET),'     /',CTIT(:LENTIT)
            ENDIF


         ELSEIF (LENT > 0) THEN
C           FILENAME AND TITLE DO NOT FIT ON SINGLE LINE
            IF (NLET > 0) WRITE(NOUT,*) ' ',FILNAM(:NLET)
            LENT = MIN(80,LENTIT)

            IF (MYPID <= 0) THEN
               WRITE(NOUT,*) ' ',CTIT(1:LENT)
               IF (LENTIT > 80)  WRITE(NOUT,*) ' ',CTIT(81:LENTIT)
            ENDIF

         ENDIF

         IF (ISTACK .NE. 0 .AND. IMGNUM == 0 .AND. NZ > 1)THEN
C           OVERALL STACKED VOLUME FILE
            CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG)

            IF (MYPID <= 0) THEN
               WRITE(NOUT,89)TYPE,NX,NY,NZ,  MAXIM,CDAT,CTIM,
     &                       DSP,LABBYT
89             FORMAT('  (',A,') ',3(I0,1X),' (.. ',I0,') CREATED ',
     &                A11, ' AT ',A, 2X,A,' HEADER BYTES: ',I0)
            ENDIF
            

         ELSEIF (ISTACK .NE. 0 .AND. IMGNUM == 0)THEN
C           OVERALL STACKED IMAGE 
            CALL LUNGETSTK(LUN,ISTACK,MAXIM,IRTFLG)

            IF (MYPID <= 0) THEN
               WRITE(NOUT,90)TYPE,NX,NY ,MAXIM,CDAT,CTIM,DSP,LABBYT
            ENDIF
90          FORMAT('  (',A,') ',2(I0,1X),' (.. ',I0,') CREATED ',A11,
     &             ' AT ',A, 2X,A,' HEADER BYTES: ',I0)


         ELSEIF (IMGNUM > 0 .AND. NZ > 1) THEN
C           STACKED VOLUME
            IF (MYPID <= 0) THEN
               WRITE(NOUT,91) TYPE,NX,NY,NZ, IMGNUM,CDAT,CTIM,DSP
            ENDIF
91          FORMAT('  (',A,') ',3(I0,1X),' (@',I0,')  CREATED ',A11,
     &             ' AT ',A, 2X,A)


         ELSE IF (IMGNUM > 0) THEN
C           STACKED IMAGE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,92)TYPE,NX,NY, IMGNUM,CDAT,CTIM,DSP
            ENDIF
92          FORMAT('  (',A,') ',2(I0,1X),' (@',I0,')  CREATED ',A11,
     &             ' AT ',A, 2X,A)


        ELSE IF (NZ > 1) THEN
C           SIMPLE VOLUME
            IF (MYPID <= 0) THEN
               WRITE(NOUT,93) TYPE,NX,NY,NZ, CDAT,CTIM,DSP,
     &                        LABBYT
93             FORMAT('  (',A,') ',3(I0,1X),' CREATED ',A11,' AT ',A,2X,
     &                A, ' HEADER BYTES: ',I0)
            ENDIF


         ELSE
C           SIMPLE IMAGE
            IF (MYPID <= 0) THEN
               WRITE(NOUT,94)TYPE,NX,NY,CDAT,CTIM,DSP,LABBYT
94             FORMAT('  (',A,') ',2(I0,1X),' CREATED ',A11,' AT ',A,2X,
     &                A,' HEADER BYTES: ',I0)
            ENDIF
         ENDIF
      ENDIF

      IRTFLG = 0

      END




C     ------------------------- LUNSETCOMMON ----------------------------

      SUBROUTINE LUNSETCOMMON(LUN,IRTFLG)

#include "LUNHDR.INC"

      INCLUDE 'CMBLOCK.INC'

      CHARACTER(LEN=12) :: CDUM
      CHARACTER(LEN=1)  :: NULL = CHAR(0)

C     LABLOCK SHOULD BE INLINED AND REMOVED IN FUTURE
      INCLUDE 'LABLOCK.INC'

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET RARE ANGLE HEADER VALUES
      KANGLE   = HEADER(30)
      PHI1     = HEADER(31)
      THETA1   = HEADER(32)
      PSI1     = HEADER(33)
      PHI2     = HEADER(34)
      THETA2   = HEADER(35)
      PSI2     = HEADER(36)

      DO I = 1,64
         HDR_VALS(I) = HEADER(36+I)
      ENDDO

C     RETRIEVE ISTACK
      CALL LUNGETISTACK(LUN,NSTACK,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RETRIEVE SIZE
      CALL LUNGETSIZE(LUN,NSAMC,NROWC,NSLICE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RETRIEVE IREC
      CALL LUNGETLAB(LUN,NDUM1,INDXREC,IREC,LABLEN,NDUM2,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      CALL LUNGETSTAT(LUN,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RECOVER FILE DATE, TIME & TITLE FROM HEADER
C     DATE NOT PASSED IN COMMON ANY MORE AS COMMON IS ONLY 10 CHAR.
      CALL LUNGETDATE(LUN,CDUM,CTIM,IRTFLG)
      CALL LUNGETTITLE(LUN,CTIT,LENTIT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     RECOVER ANGLES
      CALL LUNGETANG(LUN,IANGLE,PHI,THETA,PSI,XOFF,YOFF,ZOFF, 
     &               KANGLE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET REGISTER VALUES AS NEEDED
      CALL REG_SET(1,FLOAT(NSAMC),  NULL, IRTFLG)
      CALL REG_SET(2,FLOAT(NROWC),  NULL, IRTFLG)
      CALL REG_SET(3,FMAX,          NULL, IRTFLG)
      CALL REG_SET(4,FMIN,          NULL, IRTFLG)
      CALL REG_SET(5,AV,            NULL, IRTFLG)
      CALL REG_SET(6,SIG,           NULL, IRTFLG)
      CALL REG_SET(7,FLOAT(NSLICE), NULL, IRTFLG) 
      CALL REG_SET(8,FLOAT(NSTACK), NULL, IRTFLG) 

      END

C     -------------- LUNSETHDR --------------------------------------

      SUBROUTINE LUNSETHDR(LUNT,LUN,NX,NY,NZ,
     &                     ITYPE,ISTACK,IRTFLG)

C     INITIALIZES HEADER OBJECT. DOES NOT ZERO HEADER OBJECT LOCATIONS
C     BEYOND LENBUF!

C     SETS: NX,NY,NZ, ITYPE, ISTACK,  LENBYT
C           LABBYT, LABREC, NREC, DATE, TIME, IMGNUM, INUSE

#include "LUNHDR.INC"

      REAL, DIMENSION(LENBUF) :: BUF
      CHARACTER(LEN=1)        :: DSP
      REAL, POINTER           :: HEADERT(:) 

      IRTFLG = 1

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (LUNT > 0 .AND. LUNT <= 100) THEN
C        COPY ALL HEADER BUFFER SPACES UP TO LENBUF
C        POINT TO TRANSFER HEADER OBJECT
         CALL LUNGETOBJ(LUNT,HEADERT,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C        ZERO UNCOPIED HEADER BUFFER SPACES 
         DO I = 1,13
            HEADER(I) = 0.0
         ENDDO

         DO I = 14,21
            HEADER(I) = HEADERT(I)
         ENDDO
         DO I = 22,30
            HEADER(I) = 0.0
         END DO
         DO I = 31,LENBUF
            HEADER(I) = HEADERT(I)
         ENDDO

      ELSE
C        ZERO ALL HEADER BUFFER SPACES UP TO LENBUF
         DO I = 1,LENBUF
            HEADER(I) = 0.0
         ENDDO
      ENDIF

C     SET SIZE IN HEADER OBJECT
      CALL LUNSETSIZE(LUN,NX,NY,NZ,IRTFLG)

C     SET TYPE IN HEADER OBJECT
      CALL LUNSETTYPE(LUN,ITYPE,IRTFLG)

C     SET ISTACK IN HEADER OBJECT
      CALL LUNSETISTACK(LUN,ISTACK,IRTFLG)

C     SET RECORD ELEMENTS IN HEADER OBJECT
C     ADJUST NUMBER OF HEADER RECORDS TO HAVE >=256*4 BYTES IN HEADER
      LENBYT = NX * 4
      LABREC = 1024 / LENBYT
      IF (MOD(1024,LENBYT) .NE. 0) LABREC = LABREC + 1
      LABBYT = LABREC * LENBYT
   
C     SET TOTAL NUMBER OF RECORDS IN EACH IMAGE & HEADER
      NRECS  = NY * NZ + LABREC
      CALL LUNSETLAB(LUN,LABREC,NRECS,LABBYT,LENBYT,IRTFLG)

C     SET TIME, DATE & TITLE IN HEADER OBJECT
      CALL LUNSETTIME(LUN,IRTFLG)

      IRTFLG = 0

      END


C     ----------- LUNCLRINDX -----------------------------------------

      SUBROUTINE LUNCLRINDX(LUN,NX,IRTFLG)

C     CLEARS ALL INDEX RECORDS

#include "LUNHDR.INC"

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)

      INCLUDE 'CMBLOCK.INC'

      COMMON /IOERR/   IERR

C     AUTOMATIC ARRAY
      INTEGER :: IZEROBUF(NX)

C     SAVE CURRENT OFFSETS
      LUNSTKSAV = LUNSTK(LUN)
      LABRECSAV = LUNARA(LUN)

C     SET NO OFFSET FOR HEADER IN LUNARA
      LUNARA(LUN) = 0
C     SET NO OFFSET FOR STACKED IMAGE IN LUNSTK
      LUNSTK(LUN) = 0

C     GET THE NUMBER OF INDX RECORDS IN OVERALL HEADER 
      CALL LUNGETLAB(LUN,LABREC,INDXREC,NRECS,LABBYT,LENBYT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     NEED NX VALUE 
      CALL LUNGETSIZE(LUN,NX,NY,NZ,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     CLEAR THE INDX BUFFER
      IZEROBUF = 0.0

C     WRITE BLANK INDX RECORDS INTO OVERALL FILE HEADER
      IERR   = 0
      IRTFLG = 0
      ILOC   = LENHDR + 1

      DO IRECT = LABREC + 1, LABREC + INDXREC 

         CALL WRTLIN(LUN,IZEROBUF,NX,IRECT)
         IF (IERR .NE. 0) THEN
            CALL ERRT(102,'WHILE WRITING INDX HEADER, INDX REC.',IRECT)
            IRTFLG = 1
            EXIT
         ENDIF
      ENDDO

C     REPLACE OFFSETS IN LUNARA & LUNSTK    
      LUNARA(LUN) = LABRECSAV
      LUNSTK(LUN) = LUNSTKSAV 

      END

C     ----------- LUNWRTINDX --------------------------------------

      SUBROUTINE LUNWRTINDX(LUN,IMGNUM,NX,IRTFLG)

C     SETS INDEX FOR A SPECIFIED IMGNUM 

      INCLUDE 'CMBLOCK.INC'
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)
 
      INTEGER :: INDXBUF(NX)

C     SAVE CURRENT FILE OFFSETS
      LUNSTKSAV = LUNSTK(LUN)
      LABRECSAV = LUNARA(LUN)

C     SET NO OFFSET FOR HEADER IN LUNARA
      LUNARA(LUN) = 0
C     SET NO OFFSET FOR STACKED IMAGE IN LUNSTK
      LUNSTK(LUN) = 0

C     GET THE NUMBER OF INDX RECORDS IN OVERALL HEADER 
      CALL LUNGETLAB(LUN,LABREC,MAXNDXREC,NDUM2,NDUM3,NDUM4,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      INDXREC = IMGNUM / NX
      ILOC    = MOD(IMGNUM,NX)
      IF (ILOC .NE. 0) INDXREC = INDXREC + 1
      IF (ILOC == 0) ILOC = NX

      IF (INDXREC > MAXNDXREC) THEN
         CALL LUNGETSTKALL(LUN,ISTACK,IRTFLG)
         CALL ERRT(102,'IMAGE NUMBER EXCEEDS INDEX LIMIT',ISTACK)
         GOTO 999
      ENDIF

C     GET THE NUMBER OF CURRENT INDICES IN USE 
      CALL LUNGETINDXTOP(LUN,LASTINDX,IRTFLG)

      LASTINDX = LASTINDX + 1

C     SAVE THE NUMBER OF CURRENT INDICES IN USE 
      CALL LUNSETINDXTOP(LUN,LASTINDX,IRTFLG)

C     READ THE CORRECT INDEXS RECORD FROM THE FILE
      IRECT = LABREC + INDXREC
      CALL REDLIN(LUN,INDXBUF,NX,IRECT)

      INDXBUF(ILOC) = LASTINDX

C     WRITE THE CORRECT INDEXS RECORD BACK IN THE FILE
      CALL WRTLIN(LUN,INDXBUF,NX,IRECT)

C     REPLACE CURRENT FILE OFFSETS IN LUNARA & LUNSTK    
999   LUNARA(LUN) = LABRECSAV
      LUNSTK(LUN) = LUNSTKSAV 

C     BE SURE TO WRITE OVERALL HEADER TO FILE NOW TO SAVE INDXTOP!!

      END

C     ----------- LUNREDINDX -----------------------------------------

      SUBROUTINE LUNREDINDX(LUN,IMGNUM,INDX,NX,IRTFLG)

C     RETURNS INDEX FOR A SPECIFIED IMGNUM 

      INCLUDE 'CMBLOCK.INC'

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100)
 
C     AUTOMATIC ARRAY
      INTEGER  :: INDXBUF(NX)

C     SAVE CURRENT FILE OFFSETS FOR LUNARA & LUNSTK
      LUNSTKSAV = LUNSTK(LUN)
      LABRECSAV = LUNARA(LUN)

C     SET NO OFFSET FOR HEADER IN LUNARA
      LUNARA(LUN) = 0
C     SET NO OFFSET FOR STACKED IMAGE IN LUNSTK
      LUNSTK(LUN) = 0

C     GET THE NUMBER OF INDX RECORDS IN OVERALL HEADER 
      CALL LUNGETLAB(LUN,LABREC,MAXNDXREC,NDUM2,NDUM3,NDUM4,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      INDXREC = IMGNUM / NX
      ILOC    = MOD(IMGNUM,NX)
      IF (ILOC .NE. 0) INDXREC = INDXREC + 1
      IF (ILOC == 0) ILOC = NX

      IF (INDXREC > MAXNDXREC) THEN
         CALL LUNGETSTKALL(LUN,ISTACK,IRTFLG)
         WRITE(NOUT,*) 'IMAGE NUMBER:',IMGNUM,
     &                 'EXCEEDS INDEX LIMIT:',-ISTACK
         CALL ERRT(100,' ',NE)
         GOTO 999
      ENDIF

C     GET THE CORRECT INDEXS RECORD FROM THE FILE
      IRECT = LABREC + INDXREC
      CALL REDLIN(LUN,INDXBUF,NX,IRECT)

C     GET THE SPECIFIED INDEX VALUE
      INDX = INDXBUF(ILOC) 

C     REPLACE CURRENT FILE OFFSETS IN LUNARA & LUNSTK    
999   LUNARA(LUN) = LABRECSAV
      LUNSTK(LUN) = LUNSTKSAV 

      END

C     ------------------------- LUNSET25 ----------------------------

       SUBROUTINE LUNSET25(LUN,INUSE,IRTFLG)

C      ONLY USED FOR BACWARD FILE COMPATIBILITY WITH OLDER SPIDERS

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     SET IMUSED RELATED LOCATION
      HEADER(25) = -1 

      IRTFLG = 0

      END


C     ------------------------- LUNGET25 ----------------------------

       SUBROUTINE LUNGET25(LUN,IVAL,IRTFLG)

C      ONLY USED FOR BACWARD FILE COMPATIBILITY WITH OLDER SPIDERS

#include "LUNHDR.INC"

C     POINT TO HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     GET STACK RELATED LOCATIONS
      IVAL = HEADER(25)

      IRTFLG = 0

      END

C     ------------------------- FLIPBYTESI ----------------------------

      SUBROUTINE FLIPBYTESI(IBUF,NVAL,IRTFLG)

C     NEEDED FOR PGI COMPILER SINCE FLIPBYTES DOES NOT WORK THERE!
      INTEGER *4 IBUF(*)
      EQUIVALENCE (I,L(1)), (J,L(2)), (K,L(1))

      INTEGER * 2  I,J,L(2)
      INTEGER * 4  K

      DO N = 1,NVAL
          K       = IBUF(N)
          I       = ISHFTC(I,8,16)
          J       = ISHFTC(J,8,16)
          K       = ISHFTC(K,16,32)
          IBUF(N) = K
      ENDDO
      END

C     ------------------------- FLIPBYTES ----------------------------

      SUBROUTINE FLIPBYTES(IBUFIN,IBUFOUT,NVAL,IRTFLG)

      INTEGER *4 IBUFIN(*),IBUFOUT(*)
      EQUIVALENCE (I,L(1)), (J,L(2)), (K,L(1))

      INTEGER * 2  I,J,L(2)
      INTEGER * 4  K

      DO N = 1,NVAL
          K          = IBUFIN(N)
          I          = ISHFTC(I,8,16)
          J          = ISHFTC(J,8,16)
          K          = ISHFTC(K,16,32)
          IBUFOUT(N) = K
      ENDDO
      END

C     ------------------------- LUNSETFLIP ----------------------------

      SUBROUTINE LUNSETFLIP(LUN,IFLIP,IRTFLG)

C     SETS FLIP IN LUNFLIP(LUN)

      IMPLICIT NONE
      INTEGER  :: LUN,IFLIP,IRTFLG

      INTEGER  :: LUNARA,LUNSTK,LUNARB,LUNFLIP

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      LUNFLIP(LUN) = IFLIP

      IRTFLG       = 0

      END

C     ------------------------- LUNGETFLIP ----------------------------

      SUBROUTINE LUNGETFLIP(LUN,IFLIP,IRTFLG)

C     GETS FLIP FROM LUNFLIP(LUN)

      IMPLICIT NONE

      INTEGER  :: LUN,IFLIP,IRTFLG
      INTEGER  :: LUNARA,LUNSTK,LUNARB,LUNFLIP

      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      IFLIP  = LUNFLIP(LUN) 

      IRTFLG = 0

      END

C     ------------------------- LUNFLIPHDR ----------------------------

      SUBROUTINE LUNFLIPHDR(LUN,IRTFLG)

      !IMPLICIT NONE
#include "LUNHDR.INC"
      
      INTEGER       :: LUN,IRTFLG
      LOGICAL       :: FLIPEND
      REAL, POINTER :: HEADEROUT(:)

C     POINT TO INPUT HEADER OBJECT
      CALL LUNGETOBJ(LUN,HEADER,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     FLIP VALUES
      CALL FLIPBYTESI(HEADER,LENBUF,IRTFLG)

      IRTFLG = 0

      END

