
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
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE  PRJC2(CUBE,LTC,DM,IM,B,IPCUBE,NN,NSAM)

      INCLUDE 'PAR.INC'
C     Includes INTEGER LDPX,LDPY,LDPZ,LDPNMX,LDPNMY,NZ1,LDP,NM,LDPNM

      DIMENSION DM(9,IM)
      DIMENSION  CUBE(LTC),B(NSAM,IM)
      INTEGER  IPCUBE(5,NN)

cc$omp parallel do  private(i,j)
      DO    i=1,im
         DO    j=1,nsam
         B(j,i)=0.0
         ENDDO
      ENDDO
c$omp parallel do  private(i),schedule(static)
      DO    I=1,IM
         CALL  PRJC0(CUBE,LTC,DM(1,I),B(1,i),NSAM,IPCUBE,NN)
      ENDDO
      END
