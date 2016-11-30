C++*********************************************************************
C    QALI.F
C            3D ORIENTATION SEARCH                  2/14/94
C            ROTATION MATRIX CORRECTED
C            USED OPFILE                            DEC 99 ArDean Leith
C            OPFILEC                                FEB 03 ArDean Leith
C            PROMPTS                                JUN 13 ArDean Leith
C            NULLIFY(XPO) FOR IFORT 17              NOV 16 ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C     QALI(MODE)
C
C     PARAMETERS  MODE    ROTATION AROUND ARBITRARY POINT MODE='A'
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C **********************************************************************

        SUBROUTINE QALI(MODE)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 

        COMMON /POINT/ XPO,YPO
        REAL,  POINTER          :: XPO(:,:,:)
        REAL,  POINTER          :: YPO(:,:,:)

        CHARACTER(LEN=MAXNAM)   :: FILNAM
        CHARACTER * 1           :: MODE
        CHARACTER * 1           :: NULL = CHAR(0)
        REAL                    :: PIT(4)

        DOUBLE PRECISION        :: AA,AB
        COMMON  /QNORMA/   AA,AB
        COMMON  /DIMSPEC/  R


        DATA  LUN1,LUN2/77,78/

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'O',IFORM,NSAM,NROW,NSLICE,
     &               MAXIM,'REFERENCE VOLUME~',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0)  GOTO 9999

        IF (IFORM.NE.3 .OR. 
     &     (.NOT.(NSAM == NROW .AND. NROW == NSLICE))) THEN
           CALL  ERRT(101,'MUST BE CUBIC VOLUME ',NE)
           GOTO 9999
        ENDIF

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN2,'O',IFORM,NSM,NRW,NSL,
     &               MAXIM,'EXPERIMENTAL VOLUME~',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL SIZCHK(NULL,NSAM,NROW,NSLICE,0, NSM,NRW,NSL,0, IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL  RDPRI1S(NRAD,NOT_USED,'MASK RADIUS',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        R = NRAD

        IF (MODE == 'A')  THEN
           NZ = -999999
           CALL RDPRI3S(NX,NY,NZ,NOT_USED,
     &             'ROTATION CENTER  X,Y,Z',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (NZ == -999999) THEN
             CALL  RDPRI1S(NZ,NOT_USED,'ROTATION CENTER Z',IRTFLG)
             IF (IRTFLG .NE. 0) GOTO 9999
           ENDIF
        ENDIF

        PIT(3) = -999999
        CALL RDPRM3S(PIT(1),PIT(2),PIT(3),NOT_USED,
     &             'INITIAL EULERIAN ANGLES  PHI,THETA,PSI',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (PIT(3) == -999999) THEN
           CALL  RDPRM1S(PIT(3),NOT_USED,
     &           'INITIAL EULERIAN ANGLE; PSI',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF
 
           NSM = NSAM
           NRM = NROW
           NLM = NSLICE

           IF (MODE == 'A')  THEN
              KLX = -NX+1
              KNX = NSAM-NX
              KLY = -NY+1
              KNY = NROW-NY
              KLZ = -NZ+1
              KNZ = NSLICE-NZ

           ELSE
              KLX = -NSAM/2
              IF (MOD(NSAM,2) == 0)  THEN
                 KNX = NSAM/2-1
              ELSE
                 KNX = -KLX
              ENDIF

              KLY = -NROW/2
              IF (MOD(NROW,2) == 0)  THEN
                 KNY = NROW/2-1
              ELSE
                 KNY = -KLY
              ENDIF

              KLZ = -NSLICE/2
              IF (MOD(NSLICE,2) == 0)  THEN
                 KNZ = NSLICE/2-1
              ELSE
                 KNZ = -KLZ
              ENDIF
           ENDIF

           ALLOCATE (XPO(KLX:KNX,KLY:KNY,KLZ:KNZ), 
     &               YPO(KLX:KNX,KLY:KNY,KLZ:KNZ), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
              MWANT = (KNX-KLX) * (KLY-KNY) *(KLZ-KNZ) 
              CALL ERRT(46,'QALI, XPO...',MWANT)
              GOTO 9998
           ENDIF

        CALL READV(LUN1,XPO,NSAM,NROW,NSAM,NROW,NSLICE)

        CALL READV(LUN2,YPO,NSAM,NROW,NSAM,NROW,NSLICE)

        WRITE(NOUT,90) 
90      FORMAT(/,
     &   '  Iteration   Phi     Theta       Psi    Distance         R')

        CALL QNRF(XPO,YPO,KLX,KNX,KLY,KNY,KLZ,KNZ,R,AA,AB)
       
        CALL UQU(PIT)  

        WRITE(NOUT,91) 
91      FORMAT(/,'  Final   Phi       Theta       Psi        R')

        WRITE(NOUT,92)  (PIT(L),L=1,4)
92      FORMAT(8X,F9.4, 1X,F9.4, 1X,F9.4, 1X,F10.7)

        CALL REG_SET_NSEL(1,4,PIT(1),PIT(2),PIT(3),PIT(4),0.0,IRTFLG)

        WRITE (NDAT,93)
93      FORMAT (/,' ',30('-'),' END OF COMPUTATION ',30('-'),/)

9998    IF (ASSOCIATED(XPO)) DEALLOCATE (XPO)
	IF (ASSOCIATED(YPO)) DEALLOCATE (YPO)

        NULLIFY(XPO)  ! FOR IFORT 17  nov 2016
        NULLIFY(YPO)

9999    CLOSE(LUN1)
        CLOSE(LUN2)

        END

