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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C
C IMAGE_PROCESSING_ROUTINE                                             *
C
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        subroutine  docors3_s(bcke,bcn,n,nmat,ipcube,nn,ala,ANGS)

        INCLUDE 'PAR.INC'
C	     PAR includes INTEGER LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1,LDP,NM,LDPNM

        dimension    bcke(nmat),bcn(n,n,n),ipcube(5,nn),angs(3)
        DIMENSION  IM(3)
        DOUBLE PRECISION  RM(3,3),QR(3),THETA,PHI,PSI,DX,DY,DZ
C       EQUIVALENCE  (IM(1),IX),(IM(2),IY),(IM(3),IZ)
        psi=angs(1)/45.0*datan(1.0d0)
        theta=angs(2)/45.0*datan(1.0d0)
        phi=angs(3)/45.0*datan(1.0d0)
        RM(1,1)=DCOS(THETA)*DCOS(PHI)*DCOS(PSI)-DSIN(PHI)*DSIN(PSI)
        RM(2,1)=-DCOS(THETA)*DCOS(PHI)*DSIN(PSI)-DSIN(PHI)*DCOS(PSI)
        RM(3,1)=DSIN(THETA)*DCOS(PHI)
        RM(1,2)=DCOS(THETA)*DSIN(PHI)*DCOS(PSI)+DCOS(PHI)*DSIN(PSI)
        RM(2,2)=-DCOS(THETA)*DSIN(PHI)*DSIN(PSI)+DCOS(PHI)*DCOS(PSI)
        RM(3,2)=DSIN(THETA)*DSIN(PHI)
        RM(1,3)=-DSIN(THETA)*DCOS(PSI)
        RM(2,3)=DSIN(THETA)*DSIN(PSI)
        RM(3,3)=DCOS(THETA)
C
c$omp parallel do private(i,j,ix,iy,iz,qr,qt,a1,a2,a3,a4,a5,
c$omp& a6,a7,a8,iox,ioy,ioz,dx,dy,dz,ixx)
        do   i=1,nn
           ixx=ipcube(3,i)-1-LDP
           iy=ipcube(4,i)-LDP
           iz=ipcube(5,i)-LDP
           do   j=ipcube(1,i),ipcube(2,i)
              ix=ixx+(j+1-ipcube(1,i))
C
C       DO  3  I3=1,3
C       QR(I3)=0.0
C       DO  3  I2=1,3
C3      QR(I3)=QR(I3)+RM(I2,I3)*IM(I2)
              QR(1)=RM(1,1)*IX+RM(2,1)*IY+RM(3,1)*IZ
              QR(2)=RM(1,2)*IX+RM(2,2)*IY+RM(3,2)*IZ
              QR(3)=RM(1,3)*IX+RM(2,3)*IY+RM(3,3)*IZ
C
              IOX=QR(1)+FLOAT(LDPNM)
              DX=QR(1)+LDPNM-IOX
              DX=DMAX1(DX,1.0D-5)
              IOY=QR(2)+FLOAT(LDPNM)
              DY=QR(2)+LDPNM-IOY
              DY=DMAX1(DY,1.0D-5)
              IOZ=QR(3)+FLOAT(LDPNM)
              DZ=QR(3)+LDPNM-IOZ
              DZ=DMAX1(DZ,1.0D-5)
C
c             QT=
c     &         +(1-DX)*(1-DY)*(1-DZ)*bcn(IOX,IOY,IOZ)
c     &         +   DX *(1-DY)*(1-DZ)*bcn(IOX+1,IOY,IOZ)
c     &         +(1-DX)*   DY *(1-DZ)*bcn(IOX,IOY+1,IOZ)
c     &         +(1-DX)*(1-DY)*   DZ *bcn(IOX,IOY,IOZ+1)
c     &         +   DX *   DY *(1-DZ)*bcn(IOX+1,IOY+1,IOZ)
c     &         +   DX *(1-DY)*   DZ *bcn(IOX+1,IOY,IOZ+1)
c     &         +(1-DX)*   DY *   DZ *bcn(IOX,IOY+1,IOZ+1)
c     &         +   DX *   DY *   DZ *bcn(IOX+1,IOY+1,IOZ+1)
C faster version :
c
        a1 = bcn(IOX,IOY,IOZ) 
        a2 = bcn(IOX+1,IOY,IOZ) - bcn(IOX,IOY,IOZ) 
        a3 = bcn(IOX,IOY+1,IOZ) - bcn(IOX,IOY,IOZ) 
        a4 = bcn(IOX,IOY,IOZ+1) - bcn(IOX,IOY,IOZ)
        a5 = bcn(IOX,IOY,IOZ)-bcn(IOX+1,IOY,IOZ)-bcn(IOX,IOY+1,IOZ)
     &     + bcn(IOX+1,IOY+1,IOZ)
        a6 = bcn(IOX,IOY,IOZ)-bcn(IOX+1,IOY,IOZ)-bcn(IOX,IOY,IOZ+1)
     &     + bcn(IOX+1,IOY,IOZ+1)
        a7 = bcn(IOX,IOY,IOZ)-bcn(IOX,IOY+1,IOZ)-bcn(IOX,IOY,IOZ+1)
     &     + bcn(IOX,IOY+1,IOZ+1)
        a8 = bcn(IOX+1,IOY,IOZ)+bcn(IOX,IOY+1,IOZ)+bcn(IOX,IOY,IOZ+1)
     &   - bcn(IOX,IOY,IOZ)-bcn(IOX+1,IOY+1,IOZ)-bcn(IOX+1,IOY,IOZ+1)
     &     - bcn(IOX,IOY+1,IOZ+1) + bcn(IOX+1,IOY+1,IOZ+1)
        QT= a1 + dz*(a4 + a6*dx + (a7 + a8*dx)*dy) + a3*dy
     &         + dx*(a2 + a5*dy)
C**********************************************************
c
              bcke(j)=bcke(j)+ala*QT
           enddo
        enddo
        end

