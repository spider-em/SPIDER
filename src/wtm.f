
C ++********************************************************************
C                                                                      *
C WTM.F          ADDED ISELECT PARAMETER         SEP 03 ARDEAN LEITH   *
C                PARTITION BUG                   SEP 03 ARDEAN LEITH   *
C                                                                      *
C **********************************************************************
C=*                                                                    *
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
C                                                                      *
C  WTM(ISELECT,PROJ,W,NNNN,NSAM,NROW,SS,NANG,DIAMETER,K)               *
C                                                                      *
C  PURPOSE:     WEIGHTING FUNCTION IN 2D ACCORDING TO MVH              *
C                                                                      *
C  PARAMETERS:  ISELECT,    SELECTED PROJECTIONS LIST           (SENT) *
C               PROJ                                       (SENT/RET.) *
C               W                                               (WORK) *
C               NNNN                                            (SENT) *
C               NSAM,NROW                                       (SENT) *
C               SS                                              (SENT) *
C               NANG                                            (SENT) *
C               DIAMETER                                        (SENT) *
C               K         CURRENT PROJ/ANGLE NUMBER             (SENT) *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE WTM(ISELECT,USESELECT,PROJ,W,NNNN,NSAM,NROW,SS,NANG,
     &                 DIAMETER,NUMP)

        REAL    :: SS(6,NANG),PROJ(NNNN,NROW),W(NNNN/2,NROW)
	REAL    :: CC(3),VV(3),CP(2),VP(2),RI(3,3)
        INTEGER :: ISELECT(NANG)
        LOGICAL :: USESELECT

	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))
	PARAMETER (RAD_TO_DGR = (180.0/QUADPI))

C        R(1,1)=CPHI*CTHETA*CPSI-SPHI*SPSI
C        R(2,1)=-CPHI*CTHETA*SPSI-SPHI*CPSI
C        R(3,1)=CPHI*STHETA
C        R(1,2)=SPHI*CTHETA*CPSI+CPHI*SPSI
C        R(2,2)=-SPHI*CTHETA*SPSI+CPHI*CPSI
C        R(3,2)=SPHI*STHETA
C        R(1,3)=-STHETA*CPSI
C        R(2,3)=STHETA*SPSI
C        R(3,3)=CTHETA

         RI(1,1)=SS(1,NUMP)*SS(3,NUMP)*SS(5,NUMP)-SS(2,NUMP)*SS(6,NUMP)
         RI(2,1)=-SS(1,NUMP)*SS(3,NUMP)*SS(6,NUMP)-SS(2,NUMP)*SS(5,NUMP)
         RI(3,1)=SS(1,NUMP)*SS(4,NUMP)
         RI(1,2)=SS(2,NUMP)*SS(3,NUMP)*SS(5,NUMP)+SS(1,NUMP)*SS(6,NUMP)
         RI(2,2)=-SS(2,NUMP)*SS(3,NUMP)*SS(6,NUMP)+SS(1,NUMP)*SS(5,NUMP)
         RI(3,2)=SS(2,NUMP)*SS(4,NUMP)
         RI(1,3)=-SS(4,NUMP)*SS(5,NUMP)
         RI(2,3)=SS(4,NUMP)*SS(6,NUMP)
         RI(3,3)=SS(3,NUMP)

	NR2=NROW/2

	THICK=NSAM/DIAMETER/2.0
c$omp   parallel do private(i,j)
	DO J=1,NROW
           DO I=1,NNNN/2
C             SET W TO 1 - PROJECTION WITH ITSELF
              W(I,J) = 1.0
           ENDDO
	ENDDO

	DO  LT=1,NANG

        L = LT
        IF (USESELECT) L = ISELECT(LT)

	IF (L .NE. NUMP)  THEN
	  CC(1)=SS(2,L)*SS(4,L)*SS(3,NUMP)-SS(3,L)*SS(2,NUMP)*SS(4,NUMP)
	  CC(2)=SS(3,L)*SS(1,NUMP)*SS(4,NUMP)-SS(1,L)*SS(4,L)*SS(3,NUMP)
	  CC(3)=SS(1,L)*SS(4,L)*SS(2,NUMP)*SS(4,NUMP)-
     &		SS(2,L)*SS(4,L)*SS(1,NUMP)*SS(4,NUMP)

	  CCN=AMAX1(AMIN1(SQRT(CC(1)**2+CC(2)**2+CC(3)**2),1.0),-1.0)
	  ALPHA=RAD_TO_DGR*ASIN(CCN)
	  IF (ALPHA.GT.180.0) ALPHA=ALPHA-180.0
	  IF (ALPHA.GT.90.0) ALPHA=180.0-ALPHA
	  IF (ALPHA.LT.1.0E-6)  THEN
c$omp        parallel do private(i,j)
	     DO J=1,NROW
	        DO I=1,NNNN/2
                   W(I,J)=W(I,J)+1.0
                ENDDO
             ENDDO
	  ELSE
             FM=THICK/(ABS(SIN(ALPHA*DGR_TO_RAD)))
             CC   = CC/CCN
             VV(1)= SS(2,L)*SS(4,L)*CC(3)-SS(3,L)*CC(2)
             VV(2)= SS(3,L)*CC(1)-SS(1,L)*SS(4,L)*CC(3)
             VV(3)= SS(1,L)*SS(4,L)*CC(2)-SS(2,L)*SS(4,L)*CC(1)
             CP   = 0.0
             VP   = 0.0
             DO LL=1,2
                DO M=1,3
                   CP(LL) = CP(LL)+RI(LL,M)*CC(M)
                   VP(LL) = VP(LL)+RI(LL,M)*VV(M)
                ENDDO
             ENDDO
             TMP = CP(1)*VP(2)-CP(2)*VP(1)

C            PREVENT TMP TO BE TOO SMALL, SIGN IS IRRELEVANT
             TMP = AMAX1(1.0E-4,ABS(TMP))

c$omp        parallel do private(i,j,jy,fv,rt)
             DO  J=1,NROW
                JY = (J-1)
                IF (JY.GT.NR2)  JY=JY-NROW
                DO  I=1,NNNN/2
                   FV     = ABS((JY*CP(1)-(I-1)*CP(2))/TMP)
                   RT     = 1.0-FV/FM
                   W(I,J) = W(I,J)+AMAX1(RT,0.0)
                ENDDO
             ENDDO
          ENDIF
	ENDIF

        ENDDO

	INV = +1
        CALL FMRS_2(PROJ,NSAM,NROW,INV)
        IF (INV .EQ. 0)THEN
           CALL ERRT(38,'WTM',NE)
           RETURN
        ENDIF

c$omp   parallel do private(i,j,ww)
        DO  J=1,NROW
           DO  I=1,NNNN,2
              WW          = 1.0/W((I+1)/2,J)
              PROJ(I,J)   = PROJ(I,J)*WW
              PROJ(I+1,J) = PROJ(I+1,J)*WW
           ENDDO
	ENDDO

        INV = -1
        CALL FMRS_2(PROJ,NSAM,NROW,INV)
        END
