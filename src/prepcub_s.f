
C ++********************************************************************
C                                                                      *
C  PREPCUB_S.F          COMMON PAR REMOVED       DEC 2010 ARDEAN LEITH *                                                                  *
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
C  PREPCUB_S(N,NN,IPCUBE,RI,MD,LDP)                                    *
C                                                                      *
C  PURPOSE: MAKES A LIST OF VOXEL LOCATIONS ON EACH LINE IN THE        *
C           VOLUME WHICH ARE WITHIN A SPECIFED RADIUS SPHERE IN        *
C           VOLUME.  IF (MD= FALSE) JUST RETURNS NUMBER OF LINES IN    *
C           VOLUME WHICH NEED A VOXEL LIST. VOLUME IS A CUBE!          *
C                                                                      *
C  PARAMETERS:   N         VOLUME DIMENSION                       SENT *
C                NN        NO. OF ROWS IN VOXEL LIST       SENT OR RET *
C                IPCUBE    VOXEL LIST                              RET *
C                RI        RADIUS                                 SENT *
C                MD        IF FALSE = ONLY WANT NN RETURN         SENT *
C                LDP       CENTER OF VOLUME                       SENT *
C                                                                      *
C  VARIABLES:                                                          *
C     IPCUBE:    1 - BEGINNING VOXEL ON LINE                           *
C                2 - ENDING VOXEL ON LINE                              *
C                3 - IX FOR VOXEL                                      *
C                4 - IY FOR VOXEL                                      *
C                5 - IZ FOR VOXEL                                      *
C                                                                      *
C NOT SAME AS PREPCUB!!!! al                                           *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE PREPCUB_S(N,NN,IPCUBE,RI,MD,LDP)

        IMPLICIT NONE

        INTEGER  :: N,NN
        INTEGER  :: IPCUBE(5,NN)
        REAL     :: RI
        LOGICAL  :: MD
        INTEGER  :: LDP

        LOGICAL  :: FIRST
        REAL     :: R,T,XX,YY,RC
        INTEGER  :: NMAT,I1,I2,I3

        R    = RI * RI
        NMAT = 0
        NN   = 0

        DO I1=1,N
           T  = I1 - LDP
           XX = T * T

           DO I2=1,N
              T     = I2 - LDP
              YY    = T * T + XX
              FIRST = .TRUE.

              DO I3=1,N
                 NMAT = NMAT + 1
                 T    = I3 - LDP
                 RC   = T * T + YY

                 IF (FIRST)  THEN
C                   FIRST PIXEL ON THIS LINE,
                    IF (RC > R)  CYCLE

                    FIRST = .FALSE.
                    NN    = NN + 1

                    IF (MD) THEN
                       IPCUBE(1,NN) = NMAT
                       IPCUBE(2,NN) = NMAT
                       IPCUBE(3,NN) = I3
                       IPCUBE(4,NN) = I2
                       IPCUBE(5,NN) = I1 
                   ENDIF

                 ELSEIF (MD .AND. RC <= R) THEN 
                    IPCUBE(2,NN) = NMAT
                 ENDIF
              ENDDO
           ENDDO
        ENDDO

        END



