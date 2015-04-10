
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

	SUBROUTINE  RCNV3_P(LUN1,LUN2,B,Q,NSAM,NROW,NSLICE,PSF,NSPRD)

	DIMENSION  PSF(NSPRD,NSPRD,NSPRD),Q(NSAM,NROW,NSPRD),B(NSAM)
#ifdef SP_MP
#else
	DOUBLE PRECISION  Z
#endif
	KM(K)=MOD(K-1,NSPRD)+1
C
	NQ=NSPRD/2
	NB=NQ+1
	NSE=NSAM-NQ
	NRE=NROW-NQ
	NLE=NSLICE-NQ
C
	DO    K=1,NSPRD-1
	 CALL  RDSL_P(LUN1,Q(1,1,K),NSAM,NROW,K)
	ENDDO
	DO    K=1,NQ
	 CALL  WRSL_P(LUN2,Q(1,1,K),NSAM,NROW,K)
	ENDDO
C
	DO   K=NB,NLE
	 CALL  RDSL_P(LUN1,Q(1,1,KM(K+NQ)),NSAM,NROW,K+NQ)
	 KKTT=KM(K)
	 K1=K-NQ
	 K2=K+NQ 
	 K3=K1-1
	  DO    J=1,NQ
	   CALL  WRTLIN(LUN2,Q(1,J,KKTT),NSAM,(K-1)*NROW+J)
	  ENDDO
	DO    J=NB,NRE
	J1=J-NQ
	J2=J+NQ
	J3=J1-1
	DO    I=1,NQ
	 B(I)=Q(I,J,KKTT)
	ENDDO
#ifdef SP_MP
c$omp parallel do private(I)
	 DO    I=NB,NSE
	  B(I)=RCNV3_PS(PSF,Q,NSAM,NROW,NSPRD,NQ,I,J1,J2,J3,K1,K2,K3)
	 ENDDO
#else
	 DO    I=NB,NSE
	  Z=0.0
	  I1=I-NQ
	  I2=I+NQ
C
	   DO    KT=K1,K2
	    KTT=KM(KT)
	    K4=KT-K3
	     DO    JT=J1,J2
	      J4=JT-J3
	       DO    IT=I1,I2
	        Z=Z+Q(IT,JT,KTT)*DBLE(PSF(IT-I1+1,J4,K4))
	       ENDDO
	     ENDDO
	   ENDDO
C
	  B(I)=Z
	 ENDDO
#endif
C
	  DO    I=NSE+1,NSAM
	   B(I)=Q(I,J,KKTT)
	  ENDDO
	 CALL  WRTLIN(LUN2,B,NSAM,(K-1)*NROW+J)
	 ENDDO
C
	 DO    J=NRE+1,NROW
	  CALL  WRTLIN(LUN2,Q(1,J,KKTT),NSAM,(K-1)*NROW+J)
	 ENDDO
	ENDDO
C
	DO    K=NLE+1,NSLICE
	 CALL  WRSL_P(LUN2,Q(1,1,KM(K)),NSAM,NROW,K)
	ENDDO
	END
C
#ifdef SP_MP
	FUNCTION  RCNV3_PS(PSF,Q,NSAM,NROW,NSPRD,NQ,I,J1,J2,J3,K1,K2,K3)
	DIMENSION  PSF(NSPRD,NSPRD,NSPRD),Q(NSAM,NROW,NSPRD)
	DOUBLE PRECISION  Z
	KM(K)=MOD(K-1,NSPRD)+1
	  Z=0.0
	  I1=I-NQ
	  I2=I+NQ
C
	   DO    KT=K1,K2
	    KTT=KM(KT)
	    K4=KT-K3
	     DO    JT=J1,J2
	      J4=JT-J3
	       DO    IT=I1,I2
	        Z=Z+Q(IT,JT,KTT)*DBLE(PSF(IT-I1+1,J4,K4))
	       ENDDO
	     ENDDO
	   ENDDO
	 RCNV3_PS=Z
	END
#endif
