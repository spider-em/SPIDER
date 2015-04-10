
C ++********************************************************************
C                                                                      *
C   BPCQP         OPFILEC                         FEB 03 ARDEAN LEITH
C                 PARTITION                       MAY 03 ARDEAN LEITH
C                 PARTITION                       SEP 03 ARDEAN LEITH
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
C
C  BPCQP(PROJ,W,NNNN,NSAM,NROW,CB,NX3D,NY3D,NZC,ILIST,DM,SS,
C        NANG,NANGP,SNR,FINPAT,FINPIC)                                                                    *
C                                                                      *
C  PURPOSE:    VOLUME IN MEMORY                                                        *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE BPCQP(PROJ,W,NNNN,NSAM,NROW,CB,NX3D,NY3D,NZC,
     &           ILIST,ILISTP,DM,SS,NANG,NANGP,SNR,FINPAT,FINPIC,INPROJ)

        DIMENSION   PROJ(NNNN,NROW),CB(NX3D,NY3D,NZC),W(NNNN/2,NROW)
        DIMENSION   ILISTP(NANGP),DM(9,NANG),SS(6,NANG),ILIST(NANG)

        LOGICAL          :: PARTITION
        CHARACTER(LEN=*) :: FINPIC,FINPAT

        NLET  = LEN(FINPAT)

c$omp   parallel do private(i,j,k)
        DO K=1,NZC
           DO J=1,NY3D
              DO I=1,NX3D
                 CB(I,J,K) = 0.0
	      ENDDO
	   ENDDO
	ENDDO

        PARTITION = (NANG .NE. NANGP)

        DO K=1,NANGP
           IANGNOW = K
           IF (PARTITION) IANGNOW = ILISTP(K)

           CALL FILGET(FINPAT,FINPIC,NLET,ILISTP(K),IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,LSAM,LROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          LOAD THE PROJECTION, DO NOT USE REDVOL
           DO J=1,NROW
              CALL REDLIN(INPROJ,PROJ(1,J),NSAM,J)
           ENDDO
           CLOSE(INPROJ)

           IF (SNR .GT. 0.0) THEN
              CALL WTF(ILIST,PARTITION,PROJ,W,NNNN,NSAM,NROW,SS,
     &                 NANG,SNR,IANGNOW)
           ELSEIF (SNR .LT. 0.0) THEN
              CALL WTM(ILIST,PARTITION,PROJ,W,NNNN,NSAM,NROW,SS,
     &                 NANG,-SNR,IANGNOW)
           ENDIF

           CALL BPCQ(PROJ,NNNN,NSAM,NROW,CB,NX3D,NY3D,NZC,DM(1,IANGNOW))
	ENDDO

        END
