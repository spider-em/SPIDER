C ++********************************************************************
C                                                                      *
C  FITT                                                              *
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
C  FITT                                                              *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE FITT(X,Y,NDATA,SIG2,MWT,A,B,SIGA,SIGB)

         DIMENSION X(NDATA),Y(NDATA),SIG2(NDATA)

        SX=0.
        SY=0.
        ST2=0.
        B=0.
        IF(MWT.NE.0) THEN
        SS=0.
        DO  I=1,NDATA
          WT=1./(SIG2(I)**2)
          SS=SS+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
	ENDDO
        ELSE
          DO  I=1,NDATA
          SX=SX+X(I)
          SY=SY+Y(I)
	  ENDDO
        SS=FLOAT(NDATA)
        ENDIF
        SXOSS=SX/SS
        IF(MWT.NE.0) THEN
        DO  I=1,NDATA
          T=(X(I)-SXOSS)/SIG2(I)
          ST2=ST2+T*T
          B=B+T*Y(I)/SIG2(I)
	ENDDO
        ELSE
        DO I=1,NDATA
          T=X(I)-SXOSS
          ST2=ST2+T*T
          B=B+T*Y(I)
	ENDDO
        ENDIF
        B=B/ST2
        A=(SY-SX*B)/SS
        SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
        SIGB=SQRT(1./ST2)
        END
