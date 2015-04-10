
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



      SUBROUTINE MRKUR3(Y,N,Y0,IAUS,IMAX)


C***** AUSDRUCKEN DES GRAPHEN VON Y, MAXIMALE HOEHE IMAX=100******

C      OCTOBER 1980 M.RADERMACHER

      CHARACTER*1 A(110),BLK,ST,SI,PLU
      DIMENSION Y(N)
      COMMON /UNITS/LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT
      DATA BLK,ST,SI,PLU/' ','*','I','+'/

      LY0=0
      IF(IAUS.GT.0) WRITE(NOUT,100) Y
100   FORMAT(' *KURVE*'/20((1X,5E11.4)/))
      WRITE(NOUT,102) N,Y0
102   FORMAT(' ',///1X,'N=',I5,'  Y0=',E11.4)
      DO  K=1,IMAX
       A(K)=PLU
      ENDDO
      F=0.0
	DF=0.5/FLOAT(N)
      WRITE(NOUT,3) F,(A(LL),LL=1,IMAX)
	F=-DF
      IF(Y0.GE.0.) GOTO 4
      LY0=-INT(Y0)
4     DO1 K=1,N
      DO2 KK=1,IMAX
    2 A(KK)=BLK
      IF(LY0.GE.1.AND.LY0.LE.IMAX) A(LY0)=SI
      L=INT(Y(K)-Y0)
	F=F+DF
      IF(L.GT.IMAX) GOTO1
      IF(L.LT.1) GOTO 1
      A(L)=ST
   1  WRITE(NOUT,3) F,(A(LL),LL=1,IMAX)
      DO  K=1,IMAX
       A(K)=PLU
      ENDDO
      F=0.5
      WRITE(NOUT,3) F,(A(LL),LL=1,IMAX)
C**3  FORMAT('-',F3.2,' ',<IMAX>A1)
   3  FORMAT('-',F3.2,' ',100A1)
      END
