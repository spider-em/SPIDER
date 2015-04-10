
C++*********************************************************************
C
C    QSTATLOC.F      NEW                         SEP 2014 ArDean Leith
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
C  QSTATLOC(LUN1,LUNDOC,LUNXM)
C
C  PARAMETERS:   LUN*         I/O UNIT NUMBERS
C
C  VARIABLES:    NX,NY,NZ     DIMENSIONS  OF IMAGE
C                NSTK         STACK/MAXIM INDICATOR
C
C--*********************************************************************

      SUBROUTINE QSTATLOC(LUN1,LUNDOC,LUNXM)

      IMPLICIT NONE
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'

      INTEGER               :: LUN1,LUNDOC,LUNXM

      CHARACTER(LEN=MAXNAM) :: FILNAM

      LOGICAL               :: FOUROK
      INTEGER               :: ILIST(NIMAX)  ! NIMAX COMES FROM CMBLOCK
      INTEGER               :: NILMAX,NLET,ITYPE, NX,NY,NZ
      INTEGER               :: NDUM,NGOT,IMG1,IRTFLG
      INTEGER               :: IMAMIT,NINDX,NLOCMIN,NLOCMAX
      REAL                  :: FVALMAX,FVALMIN,FVAL, PERMIN,PERMAX
      INTEGER               :: IXMAX,IYMAX,IZMAX, IXMIN,IYMIN,IZMIN
      INTEGER               :: NSTK,IX,IY,IZ
      INTEGER               :: ICOMM,MYPID,MPIERR
      INTEGER * 8           :: NPIX8
      CHARACTER (LEN=1)     :: NULL = CHAR(0)

      REAL, ALLOCATABLE     :: FIM(:,:,:)


      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C     OPEN INPUT FILE
      NILMAX = NIMAX
      FOUROK = .TRUE.
      CALL OPFILES(0,LUN1,LUNDOC,LUNXM,  
     &             .TRUE.,FILNAM,NLET, 'O',
     &             ITYPE,NX,NY,NZ,NSTK,
     &             NULL,
     &             FOUROK, ILIST,NILMAX,
     &             NDUM,NGOT,IMG1, IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      !write(6,'(a,8i5)')' In  NSTK,ngot1,img1:',NSTK,ngot1,img1

      ALLOCATE(FIM(NX,NY,NZ), STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         CALL ERRT(46,'QSTATLOC, FIM...',NX*NY*NZ)
         RETURN
      ENDIF

      NPIX8 = NX    * NY
      NPIX8 = NPIX8 * NZ           ! DO NOT SIMPLIFY, COMPILER BUG

      IMAMIT  = IMAMI              ! FROM CMBLOCK

      NINDX  = 1
      DO                           ! LOOP OVER ALL IMAGES/STACKS

        IF (IMAMI .NE. 1) THEN    ! IMAMI IS FROM CMBLOCK
           CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV)
        ENDIF

        FVALMAX = FIM(1,1,1)
        IXMAX   = 1
        IYMAX   = 1
        IZMAX   = 1
        NLOCMAX = 1

        FVALMIN = FIM(1,1,1)
        IXMIN   = 1
        IYMIN   = 1
        IZMIN   = 1
        NLOCMIN = 1
         
C       LOAD IMAGE/VOLUME
        CALL REDVOL(LUN1,NX,NY,1,NZ,FIM,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        DO IZ = 1,NZ
           DO IY = 1,NY
             DO IX = 1,NX
                FVAL = (FIM(IX,IY,IZ))

                IF (FVAL > FVALMAX) THEN
                   FVALMAX = FVAL
                   IXMAX   = IX
                   IYMAX   = IY
                   IZMAX   = IZ
                   NLOCMAX = 1
               
                ELSEIF (FVAL == FVALMAX) THEN
                   NLOCMAX = NLOCMAX + 1
                ENDIF

                IF (FVAL < FVALMIN) THEN
                   FVALMIN = FVAL
                   IXMIN   = IX
                   IYMIN   = IY
                   IZMIN   = IZ
                   NLOCMIN = 1

                ELSEIF (FVAL == FVALMIN) THEN
                   NLOCMIN = NLOCMIN + 1
                ENDIF

             ENDDO
           ENDDO
        ENDDO

        PERMIN = 100.0 * FLOAT(NLOCMIN) / FLOAT(NPIX8)  
        PERMAX = 100.0 * FLOAT(NLOCMAX) / FLOAT(NPIX8)  

        IF (NZ > 1) THEN
           WRITE(NOUT,96) NPIX8 
96         FORMAT(/,'  TOTAL VOXELS: ',I0)
           WRITE(NOUT,90) FVALMIN,IXMIN,IYMIN,IZMIN,NLOCMIN,PERMIN 
90         FORMAT('  MIN:', 1PG10.3,' AT: (',I0,', 'I0,', ',I0,')',
     &            '  OCCURANCES: ',I0,'  %-MIN OCCURANCES:',0PF6.2)

           WRITE(NOUT,91) FVALMAX,IXMAX,IYMAX,IZMAX,NLOCMAX,PERMAX 
91         FORMAT('  MAX:', 1PG10.3,' AT: (',I0,', 'I0,', ',I0,')',
     &            '  OCCURANCES: ',I0,'  %-MAX OCCURANCES:',0PF6.2)
        ELSE
           WRITE(NOUT,97) NPIX8 
97         FORMAT(/,'  TOTAL PIXELS: ',I0)

           WRITE(NOUT,92) FVALMIN,IXMIN,IYMIN,NLOCMIN,PERMIN 
92         FORMAT('  MIN:', 1PG10.3,' AT: (',I0,', 'I0,')',
     &            '  OCCURANCES: ',I0,'  %-MIN OCCURANCES:',0PF6.2)

           WRITE(NOUT,93) FVALMAX,IXMAX,IYMAX,NLOCMAX,PERMAX 
93         FORMAT('  MAX:', 1PG10.3,' AT: (',I0,', 'I0,')',
     &            '  OCCURANCES: ',I0,'  %-MAX OCCURANCES:',0PF6.2)
        ENDIF
        WRITE(NOUT,*) ' '
 

C       OPEN NEXT INPUT FILE, UPDATE NINDX 
        CALL NEXTFILE(NINDX,   ILIST, 
     &                 FOUROK,  LUNXM,
     &                 NGOT,    NSTK,  
     &                 LUN1,    0,  
     &                 FILNAM,  'O',
     &                 IMG1,    IRTFLG)

        IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK
        IMAMIT = IMAMI               ! IMAMI IS FROM CMBLOCK

      ENDDO

      CALL REG_SET_NSEL(1,4,FVALMIN,FLOAT(IXMIN),FLOAT(IYMIN),
     &                              FLOAT(IZMIN), 0.0, IRTFLG)
      CALL REG_SET_NSEL(5,4,FVALMAX,FLOAT(IXMAX),FLOAT(IYMAX),
     &                              FLOAT(IZMAX), 0.0, IRTFLG)

9999  CLOSE(LUN1)
      CLOSE(LUNXM)

      END

