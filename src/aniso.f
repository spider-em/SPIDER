C ++********************************************************************
C                                                                      
C ANISO      CREATED                          APRIL 2002  ArDean Leith   
C            STACKS SUPPORT                   OCT   2002  ArDean Leith  
C            STACKS/SERIES SUPPORT            FEB   2020  ArDean Leith 
C                                                                  
C **********************************************************************
C=*                                                                    *
C=* Author: ArDean Leith                                               *                                                            *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2010  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  ANISO_PAR(ANS,ITER,HT,SIGMA,FLAMBDA,IRTFLG)
C  ANISO_RUN(ANS,LUN1,LUN2,NX,NY,NZ,ITER,BUF,
C            FMIN1,FMAX1,HT,SIGMA,FLAMBDA,W, IRTFLG)
C
C  PARAMETERS: LUN1,LUN2  IO UNITS                            (INPUT)
C              NX,NY,NZ   DIMENSIONS                          (INPUT)
C              IRTFLG     ERROR FLAG                          (OUTPUT)
C
C  PURPOSE: ALTER IMAGE/VOLUME CONTRAST USING ANISOTROPIC
C           DIFFUSION
C         
C  CALLER: UTIL_11
C                                                             
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

C       -------------------- ANISO_PAR -------------------------------


	SUBROUTINE ANISO_PAR(ANS,ITER,HT,SIGMA,FLAMBDA,IRTFLG)

        IMPLICIT NONE

        CHARACTER (LEN=1)  :: ANS
        INTEGER            :: ITER
        REAL               :: HT,SIGMA,FLAMBDA,W
        INTEGER            :: IRTFLG

        CHARACTER (LEN=1)  :: NULL = CHAR(0)
        INTEGER            :: NC,NOT_USED


        CALL RDPRMC(ANS,NC,.TRUE.,'CPF, MCD, OR HEG?',NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        ITER = 10
        IF (ANS == 'H') ITER = 60
        CALL RDPRI1S(ITER,NOT_USED,'ITERATIONS',IRTFLG)
	IF (IRTFLG .NE. 0) RETURN

        IF (ANS == 'H') THEN
C          USE HEGERL & FRANGAKIS FORMULATION

C          HT IS TIME STEP
           HT = 0.01
           CALL RDPRM1S(HT,NOT_USED,'TIME STEP(0...0.25)',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

C          FLAMBDA IS A CONTRAST PARAMETER,  SIGMA IS A NOISE SCALE
           SIGMA   = 1.0
           FLAMBDA = 0.01
           CALL RDPRM2S(SIGMA,FLAMBDA,NOT_USED,'SIGMA & LAMBDA',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
C          MEAN CURVATURE OR EED FORMULATION
           IF (ANS == 'M') THEN
C             MEAN CURVATURE FORMULATION
              W = 1.0
              CALL RDPRM1S(W,NOT_USED,'WEIGHTING FACTOR',IRTFLG)
	      IF (IRTFLG .NE. 0) RETURN
           ENDIF

         ENDIF
 
        END


C       -------------------- ANISO_RUN -------------------------------

	SUBROUTINE ANISO_RUN(ANS,LUN1,LUN2,NX,NY,NZ,ITER,BUF,
     &                       FMIN1,FMAX1,HT,SIGMA,FLAMBDA,W, IRTFLG)

        IMPLICIT NONE

	INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER (LEN=1)  :: ANS
        INTEGER            :: LUN1,LUN2,NX,NY,NZ,ITER
        REAL               :: BUF(NX*NY*NZ)
        REAL               :: FMIN1,FMAX1,HT,SIGMA,FLAMBDA,W
        INTEGER            :: IRTFLG

        CHARACTER (LEN=1)  :: NULL = CHAR(0)
        REAL               :: FCON,RANGE
        INTEGER            :: ILOC,I


        IF (ANS == 'H') THEN
C          USE HEGERL & FRANGAKIS FORMULATION

           IF (NZ > 3) THEN
C             VOLUME OR VOLUME STACK 
	      CALL ANISOF3(LUN1,LUN2,NX,NY,NZ,ITER,HT,
     &                     FLAMBDA,SIGMA,IRTFLG)
           ELSE
C             2D IMAGE OR IMAGE STACK
	      CALL ANISOF(LUN1,LUN2,NX,NY,NZ,ITER,HT,
     &                    FLAMBDA,SIGMA,IRTFLG)
            ENDIF
            IF (IRTFLG .NE. 0) RETURN

         ELSE
C           USE  CPF OR MCD FORMULATION
c            write(3,*)' In aniso fmin1...:',fmin1,fmax1

C           LOAD VOLUME
            CALL REDVOL(LUN1,NX,NY,1,NZ,BUF,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
       
C           NORMALIZE VOLUME OVER 0....1
            FCON = 1.0 / (FMAX1 - FMIN1) 
            BUF  = (BUF - FMIN1) * FCON
  
C           GO THRU VOLUME SLICE-BY-SLICE 
            ILOC = 1
            DO I = 1,NZ
               IF (ANS == 'M') THEN
C                 USE MEAN CURVATURE FORMULATION
                  CALL ANISOE_M(W,BUF(ILOC),NX,NY,ITER,IRTFLG)

               ELSE
C                 USE CORNER PRESERVING FORMULATION
                  CALL ANISOE(BUF(ILOC),NX,NY,ITER,IRTFLG)
               ENDIF 
               IF (IRTFLG .NE. 0) RETURN
               ILOC = ILOC + NX * NY
            ENDDO

C           UN-NORMALIZE VOLUME
            RANGE = (FMAX1 - FMIN1) 
            BUF   = (RANGE * BUF) + FMIN1

C           OUTPUT VOLUME
            CALL WRTVOL(LUN2,NX,NY,1,NZ,BUF,IRTFLG)
         ENDIF

        END








