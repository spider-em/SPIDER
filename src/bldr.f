C++*********************************************************************
C
C $$ BLDR.FOR
C
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
C
C IMAGE_PROCESSING_ROUTINE
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************
         SUBROUTINE  BLDR(R,PSI,THETA,PHI)
C
         DOUBLE  PRECISION  R(3,3),CPHI,SPHI,CTHETA,STHETA,CPSI,SPSI
	 DOUBLE PRECISION  QUADPI,DGR_TO_RAD
	 PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	 PARAMETER (DGR_TO_RAD = (QUADPI/180))
C
         CPHI=DCOS(DBLE(PHI)*DGR_TO_RAD)
         SPHI=DSIN(DBLE(PHI)*DGR_TO_RAD)
         CTHETA=DCOS(DBLE(THETA)*DGR_TO_RAD)
         STHETA=DSIN(DBLE(THETA)*DGR_TO_RAD)
         CPSI=DCOS(DBLE(PSI)*DGR_TO_RAD)
         SPSI=DSIN(DBLE(PSI)*DGR_TO_RAD)
         R(1,1)=CPHI*CTHETA*CPSI-SPHI*SPSI
         R(2,1)=-CPHI*CTHETA*SPSI-SPHI*CPSI
         R(3,1)=CPHI*STHETA
         R(1,2)=SPHI*CTHETA*CPSI+CPHI*SPSI
         R(2,2)=-SPHI*CTHETA*SPSI+CPHI*CPSI
         R(3,2)=SPHI*STHETA
         R(1,3)=-STHETA*CPSI
         R(2,3)=STHETA*SPSI
         R(3,3)=CTHETA
         END
