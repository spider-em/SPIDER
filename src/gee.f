C ++********************************************************************
C                                                                      *
C  GEE                                                                 *
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
C  GEE                                                                *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

          DOUBLE PRECISION  FUNCTION  GEE(K,N)

          COMMON TDM(2000),
     &	  PI2,AD,DEVI,X,Y,GRID,DES,WT,ALPHA,IEXT,NFCNS,NGRID 
          DIMENSION  IEXT(66),AD(66),ALPHA(66),X(66),Y(66)
          DIMENSION  DES(1056),GRID(1056),WT(1056)
          DOUBLE PRECISION  AD,DEVI,X,Y
          DOUBLE PRECISION  PI2,P,C,D,XF

          P=0.0 
          XF=GRID(K)
          XF=DCOS(PI2*XF) 
          D=0.0 
          DO    J=1,N
             C=XF-X(J) 
             C=AD(J)/C 
             D=D+C 
             P=P+C*Y(J)
	  ENDDO
          GEE=P/D 
          END 
