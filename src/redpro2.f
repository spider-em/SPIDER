
C ++********************************************************************
C                                                                      *
C   REDPRO2                                                            *
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
C  REDPRO2 (NSAM,NROWL,NROWH,NANG,PROJ,
C	   ANG,LTB,LTBN,ILIST,IPCUBE,NN,DM,RI,ABA,NOUT)                                                                    *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE  REDPRO2 (NSAM,NROWL,NROWH,NANG,PROJ,
     &	                ANG,LTB,LTBN,ILIST,IPCUBE,NN,DM,RI,ABA,NOUT)

	DIMENSION         PROJ(NSAM,NANG),ANG(NANG)
	DIMENSION         ILIST(NANG),IPCUBE(5,NN),DM(9,NANG)
	DOUBLE PRECISION  ABA,ABIN
	DOUBLE PRECISION  QUADPI,DGR_TO_RAD

	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))

	DATA      IOFF/6/

	ABA  = 0.0D0
	KLP  = 0
	ABIN = 0.0D0
	KLIN = 0
	AMI  = 1.0E23
	AMA  = -AMI

	DO K=1,NANG
           DO K2=NROWL,NROWH
              CALL REDLIN(IOFF+K,PROJ,NSAM,K2)
              CALL ASTCYL(PROJ,NSAM,RI,ABA,KLP,ABIN,KLIN,AMI,AMA)
           ENDDO

           WRITE(NOUT,333)  K,ANG(K)
333	   FORMAT('  PROJECTION:',I7,' THETA:',F6.1)

           DM(1,K) = DCOS(DBLE(ANG(K))*DGR_TO_RAD)
           DM(2,K) = 0.0
           DM(3,K) = -DSIN(DBLE(ANG(K))*DGR_TO_RAD)
           DM(4,K) = 0.0
           DM(5,K) = 1.0
           DM(6,K) = 0.0
           DM(7,K) = -DM(3,K)
           DM(8,K) = 0.0
           DM(9,K) = DM(1,K)
	ENDDO
c
	WRITE(NOUT,2044) KLIN
2044	FORMAT(/,'  TOTAL NUMBER OF POINTS IN PROJECTIONS: ',I0)

	KLIN = KLIN + KLP
	ABIN = (ABIN+ABA)/KLIN
	ABA  = ABA/KLP
	LTB  = NSAM*NANG
	LTBN = NSAM*NANG

C       PRINT STATISTICS
	WRITE(NOUT,1111)  KLIN,AMI,AMA,ABIN,ABA
1111	FORMAT('  NUMBER OF POINTS: ',I0,
     &	   /,'  MIN:',1PE10.3,'  MAX:',1PE10.3,'  AVERAGE:',1PE10.3,
     &     /,'  AVERAGE OUTSIDE THE WINDOW IS SUBTRACTED ',1PE10.3)

	END
