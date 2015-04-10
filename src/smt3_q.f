C ++********************************************************************
C                                                                      *
C   SMT3_Q                                                             *
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
C  SMT3_Q(T,ALA,X,Y,NSAM,NROW,NSLICE,IPCUBE,NN)                        *
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

        SUBROUTINE  SMT3_Q(T,ALA,X,Y,NSAM,NROW,NSLICE,IPCUBE,NN)

        DIMENSION   X(NSAM,NROW,NSLICE),Y(NSAM*NROW*NSLICE)
        DIMENSION  IPCUBE(5,NN)

        QT=1.0-6*T*ALA
        TL=T*ALA
        M=0
        DO KN=1,NN
           J=IPCUBE(4,KN)
           K=IPCUBE(5,KN)
           DO I=IPCUBE(3,KN),IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN)
              M=M+1
              Y(M)=QT*X(I,J,K)+TL*(X(I-1,J,K)+X(I+1,J,K)+
     &                      X(I,J-1,K)+X(I,J+1,K)+
     &                      X(I,J,K-1)+X(I,J,K+1))
           ENDDO
        ENDDO

        N=M+1
        DO KN=NN,1,-1
           J=IPCUBE(4,KN)
           K=IPCUBE(5,KN)
           DO I=IPCUBE(3,KN)+IPCUBE(2,KN)-IPCUBE(1,KN),IPCUBE(3,KN),-1
              N=N-1
              X(I,J,K)=Y(N)
           ENDDO
        ENDDO

C       PUT ZEROS OUTSIDE ...

        NT=1
        DO    I=1,NSAM*NROW*NSLICE
           IF(NT.GT.NN)  THEN
              Y(I)=0.0
           ELSEIF(I.LT.IPCUBE(1,NT))  THEN
              Y(I)=0.0
           ELSEIF(I.EQ.IPCUBE(2,NT))  THEN
              NT=NT+1
           ENDIF
        ENDDO
        END
