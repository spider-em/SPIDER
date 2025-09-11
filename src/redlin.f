
C++*********************************************************************
C   REDLIN.F     USED IOSTAT                    12/8/1997 ArDean Leith
C                F90 CHANGES                    4/7/1998  ArDean Leith
C                ENDEDNESS                      FEB 2003  ArDean Leith
C                MPI CHANGES                    OCT 2008  ArDean Leith
C                REDLIN_SEL                     NOV 2011  ArDean Leith
C                MRC SUPPORT                    MAY 2019  ArDean Leith
C                MORE MRC SUPPORT               DEC 2019  ArDean Leith
C                REDLIN_SEL HANDLES MRC NOW     JUN 2020  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
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
C     REDLIN(LUNT,BUF,NB,IREC)
C
C     PURPOSE:  READ A LINE OF FLOATING POINT NUMBERS FROM FILE
C
C     PARAMETERS:
C        LUNT   LOGICAL UNIT NUMBER OF FILE BEING READ       (SENT)
C        BUF    BUFFER WHERE RECORD IS READ IN               (SENT/RET.)
C        NB     NUMBER OF VALUES IN RECORD TO BE READ        (SENT)
C        IREC   RECORD TO BE READ                            (SENT)
C
C     VARIABLES:
C        IERR   ERROR CODE > 0 IS RETURNED IN CASE OF ERROR. (RETURNED)
C               IERR IS PASSED BACK IN COMMON /IOERR/IERR (LEGACY USE)
C           
C        LUNARA   NO. OF LABEL RECORDS
C        LUNSTK   CURRENT STACKED FILE = (IMGNUM-1) * IRECS + LABREC
C        LUNARB   ALTERNATE LUN FOR DUAL OPENING KLUDGE
C        LUNFLIP  0/1 FLIPBYTES FLAG 
C      
C   MRC_MODE BITS   TYPE                      MY_NBYT     RANGE
C         0    8    INTEGER*1 (SIGNED BYTES)     1      -128 -->127 
C         1   16    INTEGER*2 (SIGNED)           2 
C         2   32    REALS                        4
C         3   16    INTEGERS  (COMPLEX) FFT    NOT SUPPORTED 
C         4   32    REALS     (COMPLEX) FFT    NOT SUPPORTED      
C         6   16    INTEGER*2 (UNSIGNED)        -2
C
C-**********************************************************************

      SUBROUTINE REDLIN(LUNT,BUF,NB,IREC)

      IMPLICIT NONE

      INCLUDE 'CMLIMIT.INC'

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      INTEGER         :: LUNT,NB,IREC
      REAL            :: BUF(NB)

      INTEGER         :: IERR
      INTEGER         :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /IOERR/  IERR
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER *8      :: LUNMRCPOS
      INTEGER         :: LUNMRCNBYT
      COMMON /LUNMRC/    LUNMRCPOS(100),LUNMRCNBYT(100)

      INTEGER         :: LUN,I,IRTFLG

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER         :: MYPID = -1
#include "MPI_INIT.INC"

      IF (ISINLINE(LUNT)) THEN
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_REDLIN(LUNT,BUF,NB,IREC)
         RETURN
      ENDIF

C     FOR LUN EQUIVALENT INPUT == OUTPUT
      LUN = LUNARB(LUNT)

      IF (LUNMRCPOS(LUN) == 0) THEN

C       DIRECT ACCESS SPIDER FILE
C       ADD LUNARA(LUN) (FOR LABEL REC) AND LUNSTK (FOR STACK OFFSET)
C       TO IREC TO GET THE CORRECT RECORD NUMBER
        I = IREC + LUNARA(LUN) + LUNSTK(LUNT)

C       USING IOSTAT; IERR IS SET TO ZERO ON EACH SUCCESSFUL READ. DEC 97 al
C       IF USING MPI .AND. ONLYONE_RED IS .FALSE. IT ALWAYS READS
C       OTHERWISE, USING MPI, READS ONLY ON PROCESSOR MYPID = 0 

        IF (MYPID <= 0) READ(LUN,REC=I,IOSTAT=IERR) BUF

      ELSE

C        DIRECT ACCESS MRC FILE
         CALL REDLIN_MRC(LUN,BUF,NB,IREC, MYPID,IERR)

      ENDIF

      IF (LUNFLIP(LUN) == 1) THEN
C        FLIP INPUT BYTES IN BUF ARRAY
         CALL FLIPBYTESII(BUF,NB,IRTFLG)
      ENDIF

#ifdef USE_MPI
      IF (ONLYONE_RED) THEN
C        ALWAYS BROADCASTS IF ONLYONE_RED IS .TRUE. WHEN USING MPI
         CALL BCAST_MPI('REDLIN','BUF', BUF,NB,'R', ICOMM,MPIERR)
         CALL BCAST_MPI('REDLIN','IERR',IERR,1,'I', ICOMM,MPIERR)
      ENDIF
#endif

      END

C----------------------------- REDLIN_MRC --------------------

      SUBROUTINE REDLIN_MRC(LUN,BUF,NB,IREC, MYPID,IERR)

C     PURPOSE:  READ STREAM ACCESS MRC IMAGE/VOLUME FILE INTO BUF

      IMPLICIT NONE

      INCLUDE 'CMLIMIT.INC'   ! NEED: NBUFSIZ

      INTEGER         :: LUN,NB,IREC
      REAL            :: BUF(NB)
      INTEGER         :: MYPID,IERR

      INTEGER *8      :: LUNMRCPOS
      INTEGER         :: LUNMRCNBYT
      COMMON /LUNMRC/    LUNMRCPOS(100),LUNMRCNBYT(100)
      INTEGER         :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/   LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      INTEGER *1      :: I1BUF(NBUFSIZ)
      INTEGER *2      :: I2BUF(NBUFSIZ)

      INTEGER         :: NBYT_PER_VAL,NBYT_PER_REC,MRCTYP,NX
      INTEGER *8      :: IPOSMRC

C     LUNARA SIGN DENOTES IMAGE ORIGIN
      NX = ABS(LUNARA(LUN))

C     DIFFERENT MRC MODES DIFFER IN DATA LENGTHS
      MRCTYP       = LUNMRCNBYT(LUN)   
      NBYT_PER_VAL = ABS(MRCTYP)
      NBYT_PER_REC = NBYT_PER_VAL * NX

C     FILE POSITION TO BEGIN READING DEPENDS ON IREC

      IF (LUNARA(LUN) > 0) THEN
C        ORIGIN IS LOWER LEFT
         IPOSMRC = LUNMRCPOS(LUN) - IREC * NBYT_PER_REC

      ELSEIF (LUNARA(LUN) < 0) THEN
C        ORIGIN IS UPPER LEFT
         IPOSMRC = LUNMRCPOS(LUN) + IREC * NBYT_PER_REC

      ELSE
         IERR = 1
         RETURN
      ENDIF

C     READ STREAM ACCESS MRC IMAGE/VOLUME FILE

      IF (NBYT_PER_VAL == 4) THEN

C        MRC FILE CONTAINS 32 BIT REAL VALUES
         IF (MYPID <= 0) READ(LUN, POS=IPOSMRC,IOSTAT=IERR) BUF(1:NB)

      ELSEIF (NBYT_PER_VAL == 2) THEN
C        MRC FILE CONTAINS 16 BIT, INTEGER*2 VALUES
        
         IF (MYPID <= 0) READ(LUN, POS=IPOSMRC,IOSTAT=IERR) I2BUF(1:NB)

         BUF(1:NB) = I2BUF(1:NB)        ! SPIDER NEEDS REALS

      ELSEIF (NBYT_PER_VAL == 1) THEN
C        MRC FILE CONTAINS 8 BIT INTEGER*1 VALUES
        
         IF (MYPID <= 0) READ(LUN, POS=IPOSMRC,IOSTAT=IERR) I1BUF(1:NB)

         BUF(1:NB) = I1BUF(1:NB)        ! SPIDER NEEDS REALS
 
      ENDIF
   
      END


c      if (irec==1) write(3,'(A,2X,3i9)')
c     &            '  In redlin, lunflip(lun):',lunflip(lun)
c      if (irec==1) write(3,'(A,2X,5i9)')
c     &            '  In redlin, irec,lunmrcpos(),iposmrc:',
c     &                          irec,lunmrcpos(lun),iposmrc
c      if (irec==2) write(6,*)' irec,mrctyp,nbyt_per_val,nbyt_per_rec:'
c     &                        ,irec,mrctyp,nbyt_per_val,nbyt_per_rec,nx
c      if (irec==1) write(6,*)' irec,lunmrcpos,ipos: ',
c     &                         irec,lunmrcpos(lun),iposmrc
c     if (irec==1) write(3,'(A,2X,i9,f8.1)')
c     &            '  In redlin, i1buf,buf:',i1buf(1),buf(1)



C     ----------------------------- REDLIN1P  ------------------------

C     SHOULD REPLACE THIS WITH: REDLIN_SEL!!!
C     CALLED BY: BP32F.F, BPCG.F, BPRP3.F, BPRP.F, READV.F, REDVOL.F


#ifdef USE_MPI
C     THE SAME AS REDLIN EXCEPT NO MPI_BCAST CALLED DEPRECATED!!! USE

      SUBROUTINE REDLIN1P(LUNT,BUF,NB,IREC)

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      INTEGER         :: LUNT,NB,IREC
      REAL            :: BUF(NB)

      INTEGER         :: IERR
      COMMON /IOERR/     IERR
      INTEGER         :: LUNARA,LUNSTK,LUNARB,LUNFLIP
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER  :: MYPID = -1
#include "MPI_INIT.INC"

      IF (ISINLINE(LUNT)) THEN
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_REDLIN(LUNT,BUF,NB,IREC)
         RETURN
      ENDIF
      LUN = LUNARB(LUNT)

C     ADD LUNARA(LUN)(FOR LABEL REC) AND LUNSTK( FOR STACK OFFSET)
C     TO IREC TO GET THE CORRECT RECORD NUMBER
      I = IREC + LUNARA(LUN) + LUNSTK(LUNT)

C     USING IOSTAT; IERR IS SET TO ZERO ON EACH SUCCESSFUL READ.  DEC 97 al

      IF (MYPID <= 0) READ(LUN,REC=I,IOSTAT=IERR) BUF
      IF (LUNFLIP(LUN) == 1) CALL FLIPBYTESII(BUF,NB,IRTFLG)

      END
#endif

C     ----------------------------- REDLIN_SEL  ------------------------

      SUBROUTINE REDLIN_SEL(LUNT,NB,IREC,MPIBCAST, BUF,IRTFLG)

      IMPLICIT NONE

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      INTEGER    :: LUNT,NB,IREC,IRTFLG
      REAL       :: BUF(NB)
      LOGICAL    :: MPIBCAST

      INTEGER    ::  IERR
      COMMON /IOERR/ IERR

C     INCLUSION FOR OPTIONAL MPI INITIALIZATION.  
      INTEGER    :: MYPID = -1
#include "MPI_INIT.INC"

C     USE REDLIN FOR BOTH SPIDER AND MRC FILES
      CALL REDLIN(LUNT,BUF,NB,IREC)

#ifdef USE_MPI
      IF (MPIBCAST .AND. .NOT. ONLYONE_RED) THEN
C        ALWAYS BROADCASTS IF ONLYONE_RED IS .TRUE. WHEN USING MPI
         CALL BCAST_MPI('REDLIN','BUF', BUF,NB,'R', ICOMM,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         CALL BCAST_MPI('REDLIN','IERR',IERR,1,'I', ICOMM,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF
#endif

      IRTFLG = IERR
      END






