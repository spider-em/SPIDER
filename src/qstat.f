
C++*********************************************************************
C
C    QSTAT.F           LONG FILE NAMES JAN             89 ArDean Leith
C                      MODIFIED FOR STACKS             98 ArDean Leith
C                      ==,NSAM,IMPLICIT                16 ArDean Leith
C                      (IN HEADER) MSG                 19 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
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
C    QSTAT(LUN1,LUNM,LUNDOC,LUNXM1)
C
C    PARAMETERS:   LUN????      LOGICAL UNIT NUMBERS
C                  LUNM         LOGICAL UNIT NUMBER OF MASK
C                  NX,NY        X & Y DIMENSIONS OF IMAGE
C                  NZ           Z DIMENSION  OF IMAGE
C                  NSTACK       STACK/MAXIM INDICATOR
C
C--*******************************************************************

      SUBROUTINE QSTAT(LUN1,LUNM,LUNDOC,LUNXM1)

      IMPLICIT NONE

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER               :: LUN1,LUNM,LUNDOC,LUNXM1

      CHARACTER(LEN=MAXNAM) :: FILNAM1,FILNAMM

      LOGICAL               :: ISFIRST,FOUROK,NEWCALC
      INTEGER               :: ILIST1(NIMAX)
      INTEGER               :: NILMAX,NLET1,IFORM1,NX1,NY1,NZ1
      INTEGER               :: NDUM,NGOT1,IMG1,IRTFLG
      INTEGER               :: NXM,NYM,NZM,MAXIMM,ninndx1,npoint

      INTEGER               :: ICOMM, MYPID, MPIERR
      INTEGER               :: NSTACK1, IMAMI1, IFORMM, NINDX1, NLET 
      INTEGER               :: lnblnk 

      REAL                  :: UNUSED
      CHARACTER (LEN=1)     :: NULL = CHAR(0)

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C     OPEN INPUT FILE
      NILMAX = NIMAX
      FOUROK = .TRUE.
      CALL OPFILES(0,LUN1,LUNDOC,LUNXM1,  
     &             .TRUE.,FILNAM1,NLET1, 'O',
     &             IFORM1,NX1,NY1,NZ1,NSTACK1,
     &             NULL,
     &             FOUROK, ILIST1,NILMAX, 
     &             NDUM,NGOT1,IMG1, IRTFLG)
 
      !write(6,'(a,8i5)')' nstack1,ngot1,img1:',nstack1,ngot1,img1
      !write(3,'(a,8i5)')' In qstat; ngot,nstack,img:',nstack1,ngot1,img1
      !write(3,'(a,8i5)')' In qstat; ilist1:',ilist1(1:5)
      IF (IRTFLG .NE. 0) RETURN
      

      IMAMI1 = IMAMI   ! FROM CMBLOCK

      IF (FCHAR(4:4) == 'M') THEN
C        FIND STATISTICS UNDER A MASKED AREA

C        OPEN MASK INPUT FILE
         CALL OPFILEC(0,.TRUE.,FILNAMM,LUNM,'O',IFORMM,
     &             NXM,NYM,NZM,MAXIMM,'MASK',.FALSE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CALL SIZCHK(UNUSED,NX1,NY1,NZ1,IFORM1,
     &                      NXM,NYM,NZM,IFORMM,IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

      ENDIF

      ISFIRST = .TRUE.
       
      NINDX1  = 1
      DO                ! LOOP OVER ALL IMAGES/STACKS

        IF (FCHAR(4:4) == 'M') THEN
           CALL NORMM(LUN1,LUNM,NX1,NY1,NZ1,
     &                 FMAX,FMIN,AV,NPOINT)
 
           NEWCALC = .FALSE.
           IF (MYPID <= 0 .AND. ISFIRST) THEN
              NLET = lnblnk(FILNAMM)
              WRITE(NOUT,99) FILNAMM(1:NLET)
99            FORMAT('  STATISTICS RELATING TO MASK: ',A)

              WRITE(NOUT,98) NPOINT
98            FORMAT('  NUMBER OF POINTS: ',I8)
           ENDIF
           ISFIRST = .FALSE.
           NEWCALC = .TRUE.

        ELSEIF (IMAMI1 .NE. 1) THEN    ! IMAMI IS FROM CMBLOCK
           NEWCALC = .TRUE.
           CALL NORM3(LUN1,NX1,NY1,NZ1,FMAX,FMIN,AV)
        ENDIF

        !write(3,'(a,8i5)')' In qstat; ngot,nstack,img:',nstack1,ngot1,img1
        CALL REPORTSTAT(NEWCALC, FMIN,FMAX,AV,SIG)

C       OPEN NEXT INPUT FILE, UPDATE NINDX1 
        CALL NEXTFILE(NINDX1,   ILIST1, 
     &                FOUROK,  LUNXM1,
     &                NGOT1,   NSTACK1,  
     &                LUN1,    0,  
     &                FILNAM1, 'O',
     &                IMG1,    IRTFLG)
        IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
        IMAMI1 = IMAMI

      ENDDO

9999  CLOSE(LUN1)
      CLOSE(LUNM)
      CLOSE(LUNXM1)

      END

C     --------------- REPORTSTAT -----------------------------------

      SUBROUTINE REPORTSTAT(NEWCALC,FMINT,FMAXT,AVT,SIGT)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'

      LOGICAL           :: NEWCALC
      REAL              :: FMINT,FMAXT,AVT,SIGT
      INTEGER           :: IRTFLG

      CHARACTER(LEN=9)  :: MSG
      INTEGER           :: NC

      MSG = 'IN HEADER'
      NC  = 9
      IF (NEWCALC) THEN
         MSG = 'NEW'
         NC  = 3
      ENDIF

      WRITE(NOUT,90) FMINT,FMAXT,AVT,SIGT,MSG(1:NC)
90    FORMAT('  FMIN: ', 1PG10.3,
     &       '  FMAX: ', 1PG10.3,
     &       '  AV: ',   1PG12.5,
     &       '  SIG: ',  1PG12.5,
     &       '  (',A,')')

C     SET VALUES IN COMMAND LINE REGISTERS
      CALL REG_SET_NSEL(1,4,FMAXT,FMINT, AVT, SIGT, 0.0, IRTFLG)

C     SET VALUES IN LEGACY LOWER REGISTERS
      CALL REG_SET(3,FMAXT,CHAR(0), IRTFLG)
      CALL REG_SET(4,FMINT,CHAR(0), IRTFLG)
      CALL REG_SET(5,AVT,  CHAR(0), IRTFLG)
      CALL REG_SET(6,SIGT, CHAR(0), IRTFLG)

      END
