
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

        SUBROUTINE  RCNV2_P(LUN1,LUN2,B,Q,NSAM,NROW,PSF,NQ)

        DIMENSION  PSF(-NQ:NQ,-NQ:NQ),Q(NSAM,-NQ:NQ),B(NSAM)
#ifdef SP_MP
#else
        DOUBLE PRECISION  Z
#endif
C
        KM(K)=MOD(K+3*NQ,2*NQ+1)-NQ
        KQ(K)=MOD(K+NROW,NROW)+1
C
        DO    J=-NQ,NQ-1
        CALL  REDLIN(LUN1,Q(1,J),NSAM,KQ(J))
        ENDDO
C
        DO    J=1,NROW
        CALL  REDLIN(LUN1,Q(1,KM(J+NQ)),NSAM,KQ(J+NQ-1))
#ifdef SP_MP
c$omp parallel do private(I)
        DO    I=1,NSAM
         B(I)=RCNV2_PS(PSF,Q,NSAM,NQ,I,J)
        ENDDO
#else
        DO    I=1,NSAM
         Z=0.0D0
C
          DO    JT=-NQ,NQ
           JTM=KM(JT+J)
            DO    IT=-NQ,NQ
             Z=Z+Q(MOD(IT+I+NSAM-1,NSAM)+1,JTM)*DBLE(PSF(IT,JT))
            ENDDO
          ENDDO
C
         B(I)=Z
        ENDDO
#endif
        CALL  WRTLIN(LUN2,B,NSAM,J)
        ENDDO
        END
C
#ifdef SP_MP
        FUNCTION RCNV2_PS(PSF,Q,NSAM,NQ,I,J)
        DIMENSION  PSF(-NQ:NQ,-NQ:NQ),Q(NSAM,-NQ:NQ)
        DOUBLE PRECISION  Z
        KM(K)=MOD(K+3*NQ,2*NQ+1)-NQ
C
        Z=0.0D0
C
          DO    JT=-NQ,NQ
           JTM=KM(JT+J)
            DO    IT=-NQ,NQ
             Z=Z+Q(MOD(IT+I+NSAM-1,NSAM)+1,JTM)*DBLE(PSF(IT,JT))
            ENDDO
          ENDDO
C
         RCNV2_PS=Z
        END
#endif
