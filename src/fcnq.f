C ++********************************************************************
C                                                                      *
C                                                                      *
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
C  FCNQ(P)                                                                    *
C                                                                      *
C  PURPOSE: sloppily written using commons!                                                            *
C                                                                      *
C  PARAMETERS:  XPO,YPO                          (IMPORTED IN COMMON)                                                        *
C               AA,AB                            (IMPORTED IN COMMON)                                                        *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        REAL FUNCTION FCNQ(P)

        INCLUDE 'CMBLOCK.INC'

        REAL              :: P(3)
        INTEGER           :: IM(3)

        DOUBLE PRECISION  :: AA,AB
        POINTER           :: XPO(:,:,:), YPO(:,:,:)

        COMMON  /QNORMA/  AA,AB
        COMMON  /DIMSPEC/ R
        COMMON  /POINT/   XPO,YPO
        COMMON  /ITERU/   ITER

        DOUBLE PRECISION  :: CHI,AVD,RM(3,3),QR(3),DX,DY,DZ

        DATA  PI/3.1415926/

        KLX  = LBOUND(XPO, DIM =1)
        KNX  = UBOUND(XPO, DIM =1)
        KLY  = LBOUND(XPO, DIM =2)
        KNY  = UBOUND(XPO, DIM =2)
        KLZ  = LBOUND(XPO, DIM =3)
        KNZ  = UBOUND(XPO, DIM =3)

        ITER = ITER + 1


        PHI   = P(1) * 180.0 / PI
        THETA = P(2) * 180.0 / PI
        PSI   = P(3) * 180.0 / PI
        CHI   = 0.0
        RR    = R * R

        CALL BLDR(RM,PSI,THETA,PHI)

c$omp parallel do private(ix,iy,iz,rz,ry,rt,qr,a1,a2,a3,a4,a5,a6,a7,a8,
c$omp&   iox,ioy,ioz,dx,dy,dz,avd),reduction(+:chi)
        DO IZ=KLZ,KNZ
           RZ=IZ*IZ

           DO IY=KLY,KNY
              RY=IY*IY+RZ

              DO IX=KLX,KNX
                 RT=IX*IX+RY

                 IF (RT .LE. RR) THEN  

C                   DO  3  I3=1,3
C                   QR(I3)=0.0
C                   DO  3  I2=1,3
C3                  QR(I3)=QR(I3)+RM(I2,I3)*IM(I2)

                    QR(1) = RM(1,1)*IX+RM(2,1)*IY+RM(3,1)*IZ
                    QR(2) = RM(1,2)*IX+RM(2,2)*IY+RM(3,2)*IZ
                    QR(3) = RM(1,3)*IX+RM(2,3)*IY+RM(3,3)*IZ	            
                    IOX   = QR(1)+FLOAT(1-KLX)
                    DX    = QR(1)+FLOAT(1-KLX)-IOX
                    DX    = DMAX1(DX,1.0D-5)
                    IOX   = IOX+KLX-1
                    IOY   = QR(2)+FLOAT(1-KLY)
                    DY    = QR(2)+FLOAT(1-KLY)-IOY
                    DY    = DMAX1(DY,1.0D-5)
                    IOY   = IOY+KLY-1
                    IOZ   = QR(3)+FLOAT(1-KLZ)
                    DZ    = QR(3)+FLOAT(1-KLZ)-IOZ
                    DZ    = DMAX1(DZ,1.0D-5)
                    IOZ   = IOZ+KLZ-1

C FASTER VERSION:
                    A1 = YPO(IOX,IOY,IOZ)
                    A2 = YPO(IOX+1,IOY,IOZ) - YPO(IOX,IOY,IOZ)
                    A3 = YPO(IOX,IOY+1,IOZ) - YPO(IOX,IOY,IOZ)
                    A4 = YPO(IOX,IOY,IOZ+1) - YPO(IOX,IOY,IOZ)
       A5 = YPO(IOX,IOY,IOZ) - YPO(IOX+1,IOY,IOZ) - YPO(IOX,IOY+1,IOZ)
     &   + YPO(IOX+1,IOY+1,IOZ)
       A6 = YPO(IOX,IOY,IOZ) - YPO(IOX+1,IOY,IOZ) - YPO(IOX,IOY,IOZ+1)
     &   + YPO(IOX+1,IOY,IOZ+1)
       A7 = YPO(IOX,IOY,IOZ) - YPO(IOX,IOY+1,IOZ) - YPO(IOX,IOY,IOZ+1)
     &   + YPO(IOX,IOY+1,IOZ+1)
       A8 = YPO(IOX+1,IOY,IOZ) + YPO(IOX,IOY+1,IOZ)+ YPO(IOX,IOY,IOZ+1)
     & - YPO(IOX,IOY,IOZ)- YPO(IOX+1,IOY+1,IOZ) - YPO(IOX+1,IOY,IOZ+1)
     &   - YPO(IOX,IOY+1,IOZ+1) + YPO(IOX+1,IOY+1,IOZ+1)
       AVD = A1 + DZ*(A4 + A6*DX + (A7 + A8*DX)*DY) + A3*DY
     &   + DX*(A2 + A5*DY)

C**********************************************************
C                   CHI=CHI+(XPO(IX,IY,IZ)-AVD)*(XPO(IX,IY,IZ)-AVD)

                    CHI = CHI + XPO(IX,IY,IZ) * AVD
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        CHI  = (CHI-AB) / AA

        FCNQ = 1.0 - CHI

        WRITE(NOUT,90) ITER,(P(L)*180.0/PI,L=1,3), FCNQ, CHI
90      FORMAT(' ',I6,1x,3(1X,F9.4),2x,G13.7,1x,F8.6)

        END

