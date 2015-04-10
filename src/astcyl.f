
C++*********************************************************************
C
C  ASTCYL                               
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
C   ASTCYL
C
C--*********************************************************************

	SUBROUTINE  ASTCYL(X,NSAM,RI,ABA,KLP,ABIN,KLIN,AMI,AMA)

	DIMENSION         X(NSAM)
	DOUBLE PRECISION  ABA,ABIN

	LR=RI
	NCS=NSAM/2+1
	NCS1=MAX0(NCS-LR,1)
	NCS2=MIN0(NCS+LR,NSAM)

	DO   I=1,NSAM
	   AMA=AMAX1(AMA,X(I))
	   AMI=AMIN1(AMI,X(I))
	   IF (I.GT.NCS1 .AND. I.LT.NCS2) THEN
	      ABIN=ABIN+DBLE(X(I))
	      KLIN=KLIN+1
	   ELSE
	      ABA=ABA+DBLE(X(I))
	      KLP=KLP+1
           ENDIF
	ENDDO

	END
