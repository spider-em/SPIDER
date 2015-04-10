
C ++********************************************************************
C                                                                      *
C  MRCALERR                                                                    *
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
C  MRCALERR(XYPTS,RFPTS,ERTOT,ERVW,ERPT,PTACTIVE,NUMPTS,NTVW,NTPT)                                                                  *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE MRCALERR(XYPTS,RFPTS,ERTOT,ERVW,ERPT,
     &                 PTACTIVE,NUMPTS,NTVW,NTPT)

C       CALCULATE ALL POSSIBLE ERRORS

        PARAMETER (LV=300)
        PARAMETER (LS=256)

        LOGICAL           PTACTIVE(LS,LV)
        INTEGER           NUMPTS(LV)
        DIMENSION         XYPTS(2,LS,LV),RFPTS(2,LS,LV)
        DIMENSION         ERVW(LV),ERPT(LS),NERPT(LS)
	DOUBLE PRECISION  ERTVR,ERTOT

       ERTOT = 0.0D0
C      ERTVR = 0.0D0

      NTE = 0
      DO IPT=1,NTPT
         ERPT(IPT)  = 0.0
         NERPT(IPT) = 0
      ENDDO

      DO IVIEW=1,NTVW
          ERVW(IVIEW) = 0.0

          DO IPT=1,NTPT
            IF (PTACTIVE(IPT,IVIEW)) THEN
              ERNOW = (XYPTS(1,IPT,IVIEW) - RFPTS(1,IPT,IVIEW))**2 +
     &                (XYPTS(2,IPT,IVIEW) - RFPTS(2,IPT,IVIEW))**2
              ERVW(IVIEW) = ERVW(IVIEW) + ERNOW
              ERPT(IPT)   = ERPT(IPT)   + ERNOW
              NERPT(IPT)  = NERPT(IPT)  + 1
              ERTOT       = ERTOT + ERNOW
C             ERTVR       = ERTVR + ERNOW * ERNOW
              NTE         = NTE   + 1
            ENDIF
         ENDDO
         ERVW(IVIEW) = SQRT(ERVW(IVIEW) / NUMPTS(IVIEW))
      ENDDO

      DO  IPT=1,NTPT
         ERPT(IPT) = SQRT(ERPT(IPT) / NERPT(IPT))
      ENDDO

      ERTOT = SQRT(ERTOT/NTE)

      END

C     UNUSED BELOW -------------------------------------------
C     ERTVR=DSQRT((ERTVR-ERTOT*ERTOT*NTE)/(NTE-1))
C     PRINT *,'  AVERAGE AND SIGMA OF ERROR PER POINT ',ERTOT,ERTVR
C     DETECT POINTS WITH ERROR OVER 3*SIGMA
C      DO 700 IVIEW=1,NTVW
C          DO 700 IPT=1,NTPT
C            IF(PTACTIVE(IPT,IVIEW)) THEN
C              ERNOW= SQRT(
C     &               (XYPTS(1,IPT,IVIEW)-RFPTS(1,IPT,IVIEW))**2 +
C     &               (XYPTS(2,IPT,IVIEW)-RFPTS(2,IPT,IVIEW))**2    )
C		IF(ERNOW.GT.3*ERTVR)  THEN
C		PRINT *,' POINT #',IPT,' VIEW #',IVIEW,' ERROR=',ERNOW
C		ENDIF
C	    ENDIF
C700	CONTINUE

