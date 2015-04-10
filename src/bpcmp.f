
C ++********************************************************************
C                                                                      *
C BPCMP                                                               *
C                 PARTITION                       MAY 03 ARDEAN LEITH
C                 PARTITION BUG                   SEP 03 ARDEAN LEITH
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
C  BPCMP(PROJ,W,NNNN,NSAM,NROW,LPRJ,CB,NX3D,NY3D,NZC,ILIST,
C           DM,SS,NANG,SNR,IOPIC,FINPAT,FINPIC)  
C                                                                      *
C  PURPOSE:  PROJECTIONS IN MEMORY                                                          *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE BPCMP(PROJ,W,NNNN,NSAM,NROW,LPRJ,CB,NX3D,NY3D,NZC,
     &          ILIST,DM,SS,NANG,SNR,IOPIC,FINPAT,FINPIC,INPROJ)

        DIMENSION :: PROJ(NNNN,NROW,LPRJ),CB(NX3D,NY3D),W(NNNN/2,NROW)
        DIMENSION :: DM(9,NANG),SS(6,NANG),ILIST(NANG)

        CHARACTER(LEN=*) ::  FINPIC,FINPAT
        LOGICAL          ::  FIRST

        FIRST = .TRUE.
        NLET  = LEN(FINPAT)

        DO K=1,NANG
           IT = ILIST(K)
           CALL FILGET(FINPAT,FINPIC,NLET,IT,IRTFLG)
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,INPROJ,'O',IFORM,LSAM,LROW,NSL,
     &               MAXIM,' ',.FALSE.,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          LOAD THE PROJECTION, DO NOT USE REDVOL
           LP = MOD(K-1,LPRJ)+1
           DO J=1,NROW
              CALL REDLIN(INPROJ,PROJ(1,J,LP),NSAM,J)
	   ENDDO
           CLOSE(INPROJ)

           IF (SNR .GT. 0.0)  THEN
              CALL WTF(IDUM,.FALSE.,PROJ(1,1,LP),W,NNNN,NSAM,NROW,SS,
     &                 NANG, SNR,K)
	   ELSEIF (SNR .LT. 0.0)  THEN
              CALL WTM(IDUM,.FALSE.,PROJ(1,1,LP),W,NNNN,NSAM,NROW,SS,
     &                 NANG, -SNR,K)
	   ENDIF

           IF (LP.EQ.LPRJ .OR. K.EQ.NANG)  THEN
              CALL BPCM(PROJ,NNNN,NSAM,NROW,LP,CB,NX3D,NY3D,NZC,
     &                  DM(1,K),IOPIC,FIRST)
              FIRST = .FALSE.

           ENDIF
	ENDDO

        END
