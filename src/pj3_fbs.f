C++*********************************************************************
C                                                                      *
C PJ3_FBS.F   NEW                               DEC 2011 G KISHCHENKO  *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
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
C                                                                      *
C PJ3_FBS                                                              *
C                                                                      *
C PURPOSE:  COMPUTES A SINGLE PROJECTION OF A 3D VOLUME ACCORDING TO   *
C           THREE EULERIAN ANGLES.                                     *
C           USES 3D FBS INTERPOLATION.                                 *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PJ3_FBS()
 
        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC' 
        
        CHARACTER(LEN=MAXNAM) :: FINPIC
        REAL, ALLOCATABLE     :: VOLIN(:,:,:)
        REAL, ALLOCATABLE     :: XYZ(:,:,:)
        REAL, ALLOCATABLE     :: X1 (:,:,:)
        REAL, ALLOCATABLE     :: Y1 (:,:,:)
        REAL, ALLOCATABLE     :: Z1 (:,:,:)
        REAL, ALLOCATABLE     :: XY2(:,:,:)
        REAL, ALLOCATABLE     :: XZ2(:,:,:)
        REAL, ALLOCATABLE     :: YZ2(:,:,:)
        REAL, ALLOCATABLE     :: PRJOUT(:,:)

        INTEGER, PARAMETER    :: INPIC = 60
        INTEGER, PARAMETER    :: IOPIC = 61

        INTEGER               :: NX,NY,NZ,IRTFLG,NE,NXP,NYP,NOT_USED
        INTEGER               :: NXLD,MWANT,K,J,I,maxim,itype,nzp
        REAL                  :: PSI,PHI,THETA,TVAL

        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FINPIC,INPIC,'O',ITYPE,
     &               NX,NY,NZ,
     &               MAXIM,'INPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (ITYPE .NE. 3)  THEN
           CALL  ERRT(101,'THIS OPERATION ONLY WORKS ON VOLUMES',NE)
           GOTO 9999
        ENDIF

        NXP = NX
        NYP = 0
        CALL RDPRIS(NXP,NYP,NOT_USED,
     &              'PROJECTION DIMENSIONS IN X,Y',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        IF (NYP.EQ.0) NYP = NINT(NXP*(FLOAT(NY)/FLOAT(NX)))
        IF (NYP.EQ.0) THEN
           CALL  ERRT(102,'ILLEGAL Y DIMENSION',NYP)
           GOTO 9999
        ENDIF

        MAXIM = 0
        NZP   = 1
        ITYPE = 1
        CALL OPFILEC(0,.TRUE.,FINPIC,IOPIC,'U',ITYPE,
     &               NXP,NYP,NZP,
     &               MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999 

        PSI = HUGE(PSI)
        CALL RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &              'PHI, THETA, & PSI',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999 

        TVAL = HUGE(PSI)
        IF (PSI .EQ. TVAL) THEN 
           CALL RDPRM1S(PSI,NOT_USED,'PSI',IRTFLG)
           IF (IRTFLG.NE.0) GOTO 9999 
        ENDIF

        NXLD = NX + 2 - MOD(NX,2)

        ALLOCATE (VOLIN (NX,  NY,NZ),
     &            PRJOUT(NXP, NYP),
     &            XYZ   (NXLD,NY,NZ),
     &            X1    (NXLD,NY,NZ),
     &            Y1    (NXLD,NY,NZ),
     &            Z1    (NXLD,NY,NZ),
     &            XY2   (NXLD,NY,NZ),
     &            XZ2   (NXLD,NY,NZ),
     &            YZ2   (NXLD,NY,NZ),STAT=IRTFLG)

        IF (IRTFLG.NE.0) THEN
           MWANT = NX*NY*NZ + NXP*NYP + 7*NXLD*NY*NZ
           CALL ERRT(46,'PJ3_FBS; VOLIN,...',MWANT)
           GOTO 9999
        ENDIF

        CALL READV(INPIC,VOLIN,NX,NY,NX,NY,NZ)  ! NO PAD

C       COPY VOLIN INTO XYZ WITH ROW PADDING
        DO K=1,NZ
           DO J=1,NY
              DO I=1,NX
                 XYZ(I,J,K) = VOLIN(I,J,K)
              ENDDO

              DO I = NX+1,NXLD
                 XYZ(I,J,K) = 0
              ENDDO
           ENDDO
        ENDDO

C       WRITE(6,*), ' NXLD = ............',  NXLD
C       WRITE(6,*), ' NX = ..............',  NX
C       WRITE(6,*), ' NY = ..............',  NY
C       WRITE(6,*), ' NZ = ..............',  NZ
C       WRITE(6,*), ' PSI = .............',  PSI
C       WRITE(6,*), ' THETA = ...........',  THETA
C       WRITE(6,*), ' PHI = .............',  PHI

C       CALCULATE PROJECTION DERIVATIVES USING 3D FFT
        CALL FBS3_PREP(XYZ, NXLD, NX, NY, NZ,
     &                  X1, Y1, Z1, XY2, XZ2, YZ2)

        CALL PJ33_FBS(VOLIN,PRJOUT, NXLD,NX,NY,NZ,NXP,NYP,
     &                PSI,THETA,PHI,
     &                XYZ, X1, Y1, Z1, XY2, XZ2, YZ2)

        NZP = 1
        CALL WRITEV(IOPIC,PRJOUT, NXP,NYP,NXP,NYP,NZP) ! NOT PADDED

9999    IF (ALLOCATED(VOLIN)) DEALLOCATE (VOLIN)
        IF (ALLOCATED(XYZ))   DEALLOCATE (XYZ)
        IF (ALLOCATED(X1))    DEALLOCATE(X1)
        IF (ALLOCATED(Y1))    DEALLOCATE(Y1)
        IF (ALLOCATED(Z1))    DEALLOCATE(Z1)
        IF (ALLOCATED(XY2))   DEALLOCATE(XY2)
        IF (ALLOCATED(XZ2))   DEALLOCATE(XZ2)
        IF (ALLOCATED(YZ2))   DEALLOCATE(YZ2)
        IF (ALLOCATED(PRJOUT))DEALLOCATE(PRJOUT)

        CLOSE(IOPIC)
        CLOSE(INPIC)

        END


C       ************ SUBROUTINE  PJ33_FBS ****************************

        SUBROUTINE  PJ33_FBS(VOLIN,PRJOUT, NXLD,NX,NY,NZ,NXP,NYP,
     &                  PSI,THETA,PHI,
     &                  XYZ, X1, Y1, Z1, XY2, XZ2, YZ2)

        IMPLICIT NONE

        REAL              :: VOLIN (NX, NY,NZ) 
        REAL              :: PRJOUT(NXP,NYP)
        INTEGER           :: NX,NY,NZ,NXP,NYP
        REAL              :: PSI,THETA,PHI
        REAL              :: XYZ(NXLD,NY,NZ)
        REAL              :: X1 (NXLD,NY,NZ)
        REAL              :: Y1 (NXLD,NY,NZ)
        REAL              :: Z1 (NXLD,NY,NZ)
        REAL              :: XY2(NXLD,NY,NZ)
        REAL              :: XZ2(NXLD,NY,NZ)
        REAL              :: YZ2(NXLD,NY,NZ)

        REAL              :: DM(9)
        DOUBLE PRECISION  :: CPHI,SPHI,CTHE,STHE,CPSI,SPSI
        DOUBLE PRECISION  :: RAD_TO_DGR

        INTEGER           :: I,J,K,I2,J2
        INTEGER           :: IOX, IOY, IOZ
        INTEGER           :: N, NXLD,LDPX,LDPY,LDPZ,MPPX,MPPY

        REAL              :: XB, YB, ZB
        REAL              :: DM1,DM2,DM3

        REAL              :: fbs3

        DOUBLE PRECISION,PARAMETER :: QUADPI = 3.141592653589793238462643383279
        DOUBLE PRECISION,PARAMETER :: DGR_TO_RAD = QUADPI/180

        CPHI = DCOS(DBLE(PHI)  * DGR_TO_RAD)
        SPHI = DSIN(DBLE(PHI)  * DGR_TO_RAD)
        CTHE = DCOS(DBLE(THETA)* DGR_TO_RAD)
        STHE = DSIN(DBLE(THETA)* DGR_TO_RAD)
        CPSI = DCOS(DBLE(PSI)  * DGR_TO_RAD)
        SPSI = DSIN(DBLE(PSI)  * DGR_TO_RAD)

        DM(1) =  CPHI*CTHE*CPSI-SPHI*SPSI
        DM(2) =  SPHI*CTHE*CPSI+CPHI*SPSI
        DM(3) = -STHE*CPSI
        DM(4) = -CPHI*CTHE*SPSI-SPHI*CPSI
        DM(5) = -SPHI*CTHE*SPSI+CPHI*CPSI
        DM(6) =  STHE*SPSI

        DM(7) = STHE*CPHI
        DM(8) = STHE*SPHI
        DM(9) = CTHE

        DM1   = DM(1)
        DM2   = DM(2)
        DM3   = DM(3)

C       ZERO THE WHOLE PRJOUT ARRAY
        PRJOUT = 0.0

        LDPX = NX/2+1
        LDPY = NY/2+1
        LDPZ = NZ/2+1

        MPPX = NXP/2+1
        MPPY = NYP/2+1

        DO K=1,NZ
           DO  J=1,NY
               XB  = (1-LDPX)*DM(1) + (J-LDPY)*DM(4) +
     &               (K-LDPZ)*DM(7) + LDPX

               YB  = (1-LDPX)*DM(2) + (J-LDPY)*DM(5) +
     &               (K-LDPZ)*DM(8) + LDPY

               ZB  = (1-LDPX)*DM(3) + (J-LDPY)*DM(6) +
     &               (K-LDPZ)*DM(9) + LDPZ

c               J2  = J - LDPY + MPPY

c               IOX = IFIX(XB)
c               IOY = IFIX(YB)
c               IOZ = IFIX(ZB)

c               IF (J2  >= 1 .AND. J2  <= NYP .AND.
c     &             IOX >= 1 .AND. IOX < NX   .AND.
c     &             IOY >= 1 .AND. IOY < NY   .AND.
c     &             IOZ >= 1 .AND. IOZ < NZ)  THEN

c                   DO  I=1,NX
c                      I2 = I - LDPX + MPPX

c                      IF (I2 >= 1 .AND. I2 <= NXP) THEN


               J2 = J-LDPY+MPPY
               IF (.NOT.(J2.LT.1 .OR. J2.GT.NYP))  THEN
                   DO  I=1,NX
                     IOX    = IFIX(XB)
                     IF(.NOT.(IOX.LT.1 .OR. IOX.GE.NX))  THEN
                        IOY    = IFIX(YB)
                        IF (.NOT.(IOY.LT.1 .OR. IOY.GE.NY))  THEN
                            IOZ    = IFIX(ZB)
                            IF (.NOT.(IOZ.LT.1 .OR. IOZ.GE.NZ)) THEN
                               I2=I-LDPX+MPPX
                               IF (.NOT.(I2.LT.1 .OR. I2.GT.NXP))THEN



                          PRJOUT(I2,J2) = PRJOUT(I2,J2) + 
     &                              FBS3(XB,YB,ZB,
     &                              NXLD, NX, NY, NZ,
     &                              VOLIN,NX, XYZ,
     &                              X1, Y1, Z1,
     &                              XY2,XZ2,YZ2)

                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF

                      XB = XB + DM1
                      YB = YB + DM2
                      ZB = ZB + DM3
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         END

