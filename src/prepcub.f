
C ++********************************************************************
C                                                                      *
C  PREPCUB.F          REWRITTEN                 SEP 2003 PAWEL PENCZYK                                                                   *
C                     LDPX,LDPY,LDPZ            FEB 2005 ArDean Leith                                                        *
C                     RENAMED FROM PRPCUB_Q_N   NOV 2011 ArDean Leith                                                        *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C  PREPCUB(NSAM,NROW,NSLICE,NN,IPCUBE,RI,MD,LDPX,LDPY,LDPZ)            *
C                                                                      *
C  PURPOSE: MAKES A LIST OF VOXELS ON EACH LINE IN THE                 *
C           VOLUME  WHICH ARE WITHIN A SPECIFED RADIUS SPHERE IN       *
C           VOLUME.  IF (MD= FALSE) JUST RETURNS NUMBER OF LINES IN    *
C           VOLUME WHICH NEED A VOXEL LIST.                            *                               *
C                                                                      *
C  PARAMETERS:   NSAM..    VOLUME DIMENSION                       SENT *
C                NN        NO. OF ROWS IN VOXEL LIST       SENT OR RET *
C                IPCUBE    VOXEL LIST                              RET *
C                RI        RADIUS                                 SENT *
C                MD        IF FALSE = ONLY WANT NN RETURN         SENT *
C                LDPX      CENTER OF VOLUME                       SENT *
C                LDPY      CENTER OF VOLUME                       SENT *
C                LDPY      CENTER OF VOLUME                       SENT *
C                                                                      *
C     IPCUBE:    1 - BEGINNING VOXEL ON LINE                           *
C                2 - ENDING VOXEL ON LINE                              *
C                3 - IX FOR VOXEL                                      *
C                4 - IY FOR VOXEL                                      *
C                5 - IZ FOR VOXEL                                      *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE PREPCUB(NSAM,NROW,NSLICE,
     &                     NN,IPCUBE,RI,MD,LDPX,LDPY,LDPZ)

        IMPLICIT NONE

        INTEGER         :: NSAM,NROW,NSLICE,NN
        INTEGER         :: IPCUBE(5,NN)
        REAL            :: RI
        LOGICAL         :: MD
        INTEGER         :: LDPX,LDPY,LDPZ

        LOGICAL         :: FIRST
        INTEGER         :: I1,I2,I3
        REAL            :: R,T,XX,YY,RC

        R  = RI * RI

        NN = 0                  ! RETURNED

        DO I1=1,NSLICE
           T  = I1 - LDPZ
           XX = T*T

           DO I2=1,NROW
              T     = I2 - LDPY
              YY    = T * T + XX
              FIRST = .TRUE.

              DO I3=1,NSAM
                 T  = I3 - LDPX
                 RC = T * T + YY

                 IF (FIRST)  THEN
C                   FIRST PIXEL ON THIS LINE,
                    IF (RC .GT. R)  CYCLE
                    FIRST = .FALSE.

                    NN = NN + 1
                    IF (MD) THEN
                       IPCUBE(1,NN) = I3
                       IPCUBE(2,NN) = I3
                       IPCUBE(3,NN) = I3-LDPX
                       IPCUBE(4,NN) = I2-LDPY
                       IPCUBE(5,NN) = I1-LDPZ
                   ENDIF
                 ELSE
C                   SECOND OR LATER PIXEL ON THIS LINE,
                    IF (MD) THEN
                       IF (RC .LE. R)  IPCUBE(2,NN) = I3
                    ENDIF
                ENDIF
              ENDDO
           ENDDO
        ENDDO

        END
