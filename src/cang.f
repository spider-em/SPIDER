
C ++********************************************************************
C                                                                      *
C  CANG         ADDED DOSS                         JUL 03 ARDEAN LEITH *
C                                                                      *
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
C  CANG(PHI,THETA,PSI,DOSS,SS,DM)                                      *
C                                                                      *
C  PURPOSE:    CONVERTS EULER ANGLES TO ROTATION MATRIX                *
C                                                                      *
C  PARAMETERS: PHI,THETA,PSI     ANGLES                           SENT *
C              DOSS              FLAG TO RETURN SS ANGLES         SENT *
C              SS                ANGLES                           RET  *
C              DM                ROTATION MATRIX                  RET  *
C                                                                      *
C  CALLED BY:  buildm.f, bp3f.f, bprp.f, builds.f, bpcg.f              *
C                                                                      *
C  SEE ALSO:   bldr.f, hiang,f                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE CANG(PHI,THETA,PSI,DOSS,SS,DM)

        REAL              :: PHI,THETA,PSI
        LOGICAL           :: DOSS
        DOUBLE PRECISION  :: CPHI,SPHI,CTHE,STHE,CPSI,SPSI
        REAL              :: DM(9),SS(6)

	DOUBLE PRECISION  :: QUADPI, DGR_TO_RAD
	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))

        CPHI = DCOS(DBLE(PHI)   * DGR_TO_RAD)
        SPHI = DSIN(DBLE(PHI)   * DGR_TO_RAD)
        CTHE = DCOS(DBLE(THETA) * DGR_TO_RAD)
        STHE = DSIN(DBLE(THETA) * DGR_TO_RAD)
        CPSI = DCOS(DBLE(PSI)   * DGR_TO_RAD)
        SPSI = DSIN(DBLE(PSI)   * DGR_TO_RAD)

        IF (DOSS) THEN
C          WANT TO RETURN SS 
	   SS(1) = SNGL(CPHI)
	   SS(2) = SNGL(SPHI)
	   SS(3) = SNGL(CTHE)
	   SS(4) = SNGL(STHE)
	   SS(5) = SNGL(CPSI)
	   SS(6) = SNGL(SPSI)
        ENDIF

        DM(1) =  CPHI*CTHE*CPSI - SPHI*SPSI
        DM(2) =  SPHI*CTHE*CPSI + CPHI*SPSI
        DM(3) = -STHE*CPSI
        DM(4) = -CPHI*CTHE*SPSI - SPHI*CPSI
        DM(5) = -SPHI*CTHE*SPSI + CPHI*CPSI
        DM(6) =  STHE*SPSI
        DM(7) =  STHE*CPHI
        DM(8) =  STHE*SPHI
        DM(9) =  CTHE

        END

#ifdef NEVER
        SUBROUTINE BLDR(R,PSI,THETA,PHI)

        DOUBLE  PRECISION :: R(3,3)
        REAL              :: PSI,THETA,PHI 

        DOUBLE  PRECISION :: CPHI,SPHI,CTHETA,STHETA,CPSI,SPSI
	DOUBLE PRECISION  :: QUADPI,DGR_TO_RAD

	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))

        CPHI   = DCOS(DBLE(PHI)   * DGR_TO_RAD)
        SPHI   = DSIN(DBLE(PHI)   * DGR_TO_RAD)
        CTHETA = DCOS(DBLE(THETA) * DGR_TO_RAD)
        STHETA = DSIN(DBLE(THETA) * DGR_TO_RAD)
        CPSI   = DCOS(DBLE(PSI)   * DGR_TO_RAD)
        SPSI   = DSIN(DBLE(PSI)   * DGR_TO_RAD)

        R(1,1) =  CPHI*CTHETA*CPSI - SPHI*SPSI
        R(2,1) = -CPHI*CTHETA*SPSI - SPHI*CPSI
        R(3,1) =  CPHI*STHETA
        R(1,2) =  SPHI*CTHETA*CPSI + CPHI*SPSI
        R(2,2) = -SPHI*CTHETA*SPSI + CPHI*CPSI
        R(3,2) =  SPHI*STHETA
        R(1,3) = -STHETA*CPSI
        R(2,3) =  STHETA*SPSI
        R(3,3) =  CTHETA

        END
#endif
