C++*********************************************************************
C
C HIANG.F
C
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
C
C  HIANG(ANG,NANG,DM,LB,LO)
C
C  PURPOSE:     COMPRESS ANGLES AND CONVERT TO DTM FORMAT
C
C  PARAMETERS:   ANG                                          SENT/RET
C                NANG                                         SENT
C                DM                                           RET
C                LB                                           SENT/RET
C                LO                                           RET
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE HIANG(ANG,NANG,DM,LB,LO)

        DIMENSION         :: ANG(3,NANG),LB(NANG),DM(9,NANG)
        DOUBLE PRECISION  :: CPHI,SPHI,CTHE,STHE,CPSI,SPSI
        DOUBLE PRECISION  :: QUADPI,DGR_TO_RAD

        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (DGR_TO_RAD = (QUADPI/180))

C       ANG(1 - PHI, ANG(2 - THETA, ANG(3 - PSI
        DO I=1,NANG
           LB(I )= 1
           IF (ANG(2,I) .GE. 90.0)  THEN
              ANG(2,I) = 180.0 - ANG(2,I)
              ANG(1,I) = ANG(1,I) + 180.0
              IF (ANG(1,I) .GE. 360.0)  ANG(1,I) = ANG(1,I)-360.0
           ENDIF
        ENDDO

        DO I=1,NANG
           JJJ      = INT(100.0 * ANG(1,I))
           ANG(1,I) = FLOAT(JJJ) / 100.0
        ENDDO

        DO I=1,NANG-1
           IF (LB(I) .NE. 0)  THEN
              DO  J=I+1,NANG
                 IF (LB(I) .NE. 0)  THEN
                    IF (ANG(1,I) .EQ. ANG(1,J) .AND. 
     &                  ANG(2,I) .EQ. ANG(2,J)) THEN
                       LB(I) = LB(I) + 1
                       LB(J) = 0
                    ENDIF
                 ENDIF
              ENDDO
           ENDIF
        ENDDO

        LO = 0
        DO I=1,NANG
           IF (LB(I) .NE. 0) THEN
              LO        = LO+1
              ANG(1,LO) = ANG(1,I)
              ANG(2,LO) = ANG(2,I)

C -           PSI SET TO ZERO - IRRELEVANT
              ANG(3,LO) = 0.0

              LB(LO)    = LB(I)
           ENDIF
        ENDDO
C       NANG=LO

        DO I=1,LO
           CPHI    = DCOS(DBLE(ANG(1,I))*DGR_TO_RAD)
           SPHI    = DSIN(DBLE(ANG(1,I))*DGR_TO_RAD)
           CTHE    = DCOS(DBLE(ANG(2,I))*DGR_TO_RAD)
           STHE    = DSIN(DBLE(ANG(2,I))*DGR_TO_RAD)
           CPSI    = 1.0D0
           SPSI    = 0.0D0

           DM(1,I) =  CPHI*CTHE*CPSI-SPHI*SPSI
           DM(2,I) =  SPHI*CTHE*CPSI+CPHI*SPSI
           DM(3,I) = -STHE*CPSI
           DM(4,I) = -CPHI*CTHE*SPSI-SPHI*CPSI
           DM(5,I) = -SPHI*CTHE*SPSI+CPHI*CPSI
           DM(6,I) =  STHE*SPSI
           DM(7,I) =  STHE*CPHI
           DM(8,I) =  STHE*SPHI
           DM(9,I) =  CTHE     
        ENDDO

        END
