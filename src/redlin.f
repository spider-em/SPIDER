
C++*********************************************************************
C   REDLIN.F     USED IOSTAT                    12/8/1997 ArDean Leith
C                F90 CHANGES                    4/7/1998  ArDean Leith
C                ENDEDNESS                      FEB 2003  ArDean Leith
C                MPI CHANGES                    OCT 2008  ArDean Leith
C                REDLIN_SEL                     NOV 2011  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C     REDLIN(LUNT,BUF,NB,NREC)
C
C     PURPOSE:  READ A LINE OF FLOATING POINT NUMBERS FROM FILE
C
C     PARAMETERS:
C        LUNT   LOGICAL UNIT NUMBER OF FILE BEING READ       (SENT)
C        BUF    BUFFER WHERE RECORD IS READ IN               (SENT/RET.)
C        NB     NUMBER OF VALUES IN RECORD TO BE READ        (SENT)
C        NREC   RECORD TO BE READ                            (SENT)
C
C     VARIABLES:
C        IERR   ERROR CODE > 0 IS RETURNED IN CASE OF ERROR. (RETURNED)
C               IERR IS PASSED IN COMMON /IOERR/IERR
C 
C--*******************************************************************

      SUBROUTINE REDLIN(LUNT,BUF,NB,NREC)

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      REAL            BUF(NB)
      COMMON /IOERR/  IERR
      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      LOGICAL            :: ISOPEN

#ifdef USE_MPI
      LOGICAL            :: ONLYONE_RED,ONLYONE_WRT
      COMMON /COMM_MPI/ONLYONE_RED,ONLYONE_WRT
      include 'mpif.h'

      ICOMM  = MPI_COMM_WORLD
      CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)

C     ONLYONE_RED IS ALWAYS .TRUE. EXCEPT FROM dsgr_p.f
      IF (.NOT. ONLYONE_RED) MYPID = -1
#else
      MYPID = -1
#endif

      IF (ISINLINE(LUNT)) THEN
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_REDLIN(LUNT,BUF,NB,NREC)
         RETURN
      ENDIF
      LUN = LUNARB(LUNT)

C     ADD LUNARA(LUN) (FOR LABEL REC) AND LUNSTK (FOR STACK OFFSET)
C     TO NREC TO GET THE CORRECT RECORD NUMBER
      I = NREC + LUNARA(LUN) + LUNSTK(LUNT)

C     USING IOSTAT; IERR IS SET TO ZERO ON EACH SUCCESSFUL READ. DEC 97 al
C     IF USING MPI .AND. ONLYONE_RED IS .FALSE. IT ALWAYS READS
C     OTHERWISE, USING MPI, READS ONLY ON PROCESSOR MYPID = 0 
      IF (MYPID .LE. 0)  READ(LUN,REC=I,IOSTAT=IERR) BUF

#ifdef USE_MPI
      IF (ONLYONE_RED) THEN
C        ALWAYS BROADCASTS IF ONLYONE_RED IS .TRUE. WHEN USING MPI
         CALL BCAST_MPI('REDLIN','BUF', BUF,NB,'R', ICOMM,MPIERR)
         CALL BCAST_MPI('REDLIN','IERR',IERR,1,'I', ICOMM,MPIERR)
      ENDIF
#endif

      IF (LUNFLIP(LUN) .EQ. 1) THEN
         CALL FLIPBYTESI(BUF,NB,IRTFLG)
      ENDIF

      END

C     ----------------------------- REDLIN1P  ------------------------

C     SHOULD REPLACE THIS WITH: REDLIN_SEL!!!

#ifdef USE_MPI
C     THE SAME AS REDLIN EXCEPT NO MPI_BCAST CALLED DEPRECATED!!! USE

      SUBROUTINE REDLIN1P(LUNT,BUF,NB,NREC)

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      REAL            BUF(NB)
      COMMON /IOERR/  IERR
      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      include 'mpif.h'
      ICOMM   = MPI_COMM_WORLD
      CALL MPI_COMM_RANK(ICOMM, MYPID, MPIERR)

      IF (ISINLINE(LUNT)) THEN
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_REDLIN(LUNT,BUF,NB,NREC)
         RETURN
      ENDIF
      LUN = LUNARB(LUNT)

C     ADD LUNARA(LUN)(FOR LABEL REC) AND LUNSTK( FOR STACK OFFSET)
C     TO NREC TO GET THE CORRECT RECORD NUMBER
      I = NREC + LUNARA(LUN) + LUNSTK(LUNT)

C     USING IOSTAT; IERR IS SET TO ZERO ON EACH SUCCESSFUL READ.  DEC 97 al

      IF (MYPID .LE. 0) READ(LUN,REC=I,IOSTAT=IERR) BUF
      IF (LUNFLIP(LUN) .EQ. 1) CALL FLIPBYTESI(BUF,NB,IRTFLG)

      END
#endif

C     ----------------------------- REDLIN_SEL  ------------------------

      SUBROUTINE REDLIN_SEL(LUNT,NB,NREC,MPIBCAST, BUF,IRTFLG)

      IMPLICIT NONE

C     USE INLINE BUFFER COMMON AREA
      INCLUDE 'INLN_INFO.INC'

      INTEGER    :: LUNT,NB,NREC
      REAL       :: BUF(NB)
      LOGICAL    :: MPIBCAST

      INTEGER    :: LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT,IERR
      INTEGER    :: LUNARA(100),LUNSTK(100),LUNARB(100),LUNFLIP(100)

      COMMON /IOERR/  IERR
      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      COMMON /LUNARA/ LUNARA,LUNSTK,LUNARB,LUNFLIP

      LOGICAL    :: ISOPEN
      INTEGER    :: MYPID,LUN,I,IRTFLG

#ifdef USE_MPI
      INTEGER    :: ICOMM
      LOGICAL    :: ONLYONE_RED,ONLYONE_WRT
      COMMON /COMM_MPI/ONLYONE_RED,ONLYONE_WRT
      include 'mpif.h'

      ICOMM  = MPI_COMM_WORLD
      CALL MPI_COMM_RANK(ICOMM, MYPID, IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     ONLYONE_RED IS ALWAYS .TRUE. EXCEPT FROM: dsgr_p.f
      IF (.NOT. ONLYONE_RED) MYPID = -1
#else
      MYPID = -1
#endif

      IF (ISINLINE(LUNT)) THEN
C        USE INLINED BUFFER FOR I/O (SEE OPENINLN.F)
         CALL INLN_REDLIN(LUNT,BUF,NB,NREC)
         IRTFLG = 0
         RETURN
      ENDIF
      LUN = LUNARB(LUNT)

C     ADD LUNARA(LUN) (FOR LABEL REC) AND LUNSTK (FOR STACK OFFSET)
C     TO NREC TO GET THE CORRECT RECORD NUMBER
      I = NREC + LUNARA(LUN) + LUNSTK(LUNT)

C     USING IOSTAT; IERR IS SET TO ZERO ON EACH SUCCESSFUL READ. DEC 97 al
C     IF USING MPI .AND. ONLYONE_RED IS .FALSE. IT ALWAYS READS
C     OTHERWISE, USING MPI, READS ONLY ON PROCESSOR MYPID = 0 
      IF (MYPID .LE. 0) READ(LUN,REC=I,IOSTAT=IERR) BUF
      IF (IERR .NE. 0) THEN
         IRTFLG = IERR
         RETURN
      ENDIF

#ifdef USE_MPI
      IF (ONLYONE_RED .OR. MPIBCAST) THEN
C        ALWAYS BROADCASTS IF ONLYONE_RED IS .TRUE. WHEN USING MPI
         CALL BCAST_MPI('REDLIN','BUF', BUF,NB,'R', ICOMM,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         CALL BCAST_MPI('REDLIN','IERR',IERR,1,'I', ICOMM,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF
#endif

      IRTFLG = 0
      IF (LUNFLIP(LUN) .EQ. 1) THEN
         CALL FLIPBYTESI(BUF,NB,IRTFLG)
      ENDIF

      END

