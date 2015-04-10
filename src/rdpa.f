
C ++********************************************************************
C                                                                      *
C  RDPA                                                                     *
C                                                                      *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C  RDPA(NANG,ANG,DM)                                                                     *
C                                                                      *
C  PURPOSE:  FILLS DM WITH ANGLES MATRIX                                                          *
C                                                                      *
C  PARAMETERS:  NANG      NUMBER OF ANGLES                    (SENT)   *
C               ANG       THETA                               (SENT)
C               DM        ANGLES MATRIX                       (RET.)
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE  RDPA(NANG,ANG,DM)

        INCLUDE 'CMBLOCK.INC' 

	DIMENSION  ANG(NANG),DM(9,NANG)
	DOUBLE PRECISION  QUADPI,DGR_TO_RAD
	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))

	DO K=1,NANG

           WRITE(NOUT,333) K,ANG(K)
333        FORMAT('  Projection: ',I5,' THETA: ',F7.3)

           DM(1,K) = DCOS(DBLE(ANG(K))  * DGR_TO_RAD)
           DM(2,K) = 0.0
           DM(3,K) = -DSIN(DBLE(ANG(K)) * DGR_TO_RAD)
           DM(4,K) = 0.0
           DM(5,K) = 1.0
           DM(6,K) = 0.0
           DM(7,K) = -DM(3,K)
           DM(8,K) = 0.0
           DM(9,K) = DM(1,K)
	ENDDO

	END
