C ++********************************************************************
C                                                                      *
C                                                                      *
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C        0         2         3         4         5         6         7 *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

	SUBROUTINE  OPSF3_P(LUN,ROT,NF,Q,NFIL)

	DIMENSION  Q(NFIL)
	DOUBLE PRECISION  ROT(29,29,29)

	DO    I=NF,1,-1
	   DO    K=1,NF
	      NR1=NF-I+1+(K-1)*NFIL
	      NR2=I+NF-1+(K-1)*NFIL
	      NR3=NF-I+1+(NFIL-K)*NFIL
	      NR4=I+NF-1+(NFIL-K)*NFIL
	      DO    J=1,NF
	         Q(J)=ROT(NF-J+1,I,NF-K+1)
	         Q(NFIL-J+1)=ROT(NF-J+1,I,NF-K+1)
	      ENDDO
	      CALL  WRTLIN(LUN,Q,NFIL,NR1)
	      CALL  WRTLIN(LUN,Q,NFIL,NR2)
	      CALL  WRTLIN(LUN,Q,NFIL,NR3)
	   CALL  WRTLIN(LUN,Q,NFIL,NR4)
	   ENDDO
	ENDDO
	END
