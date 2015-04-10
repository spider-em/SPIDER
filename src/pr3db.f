C++*********************************************************************
C
C PR3DB.F              COSMETIC                  FEB 2011 ARDEAN LEITH  
C                      CYCLE                     FEB 2012 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C PURPOSE: CALCULATE THE 3-D PHASE RESIDUE OUTSIDE MISSING CONE, 
C          PR OF FOURIER RINGS(RADIUS, DIRECTION RELATIVE TO Z) AND 
C          OF FOURIER SHELLS(RADIUS) IS CALCULATED. 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE PR3DB(A,B,PR,AMP,CSUM1,LR,CSUM,CSUM2,
     &     AVSUM,LSD,NX,NY,NZ,DSCALE,NSCALE,SCALE1,SSANG,
     &     INC,Y1,WI,SER)

        REAL            :: A(LSD,NY,NZ), B(LSD,NY,NZ)
        REAL            :: PR(NSCALE,INC), AMP(NSCALE,INC)
        REAL            :: AVSUM(NSCALE,INC)
        REAL            :: CSUM1(INC),CSUM2(INC)
        INTEGER         :: LR(INC)
        REAL            :: CSUM(INC)
        CHARACTER*1     :: SER

        REAL, PARAMETER :: QUADPI     = 3.1415926535897932384626
        REAL, PARAMETER :: RAD_TO_DGR = (180.0/QUADPI)
        REAL, PARAMETER :: DGR_TO_RAD = (QUADPI/180)

        ZANG = (90.0 - SSANG) * DGR_TO_RAD
        ND2  = NX / 2
        NR2  = NY / 2

C       FOR 2D CASE SET NS2 TO ONE
        NS2 = MAX(1,NZ/2)

        DO  L=1,INC
           LR(L)    = 0
           CSUM(L)  = 0.0
           CSUM1(L) = 0.0
           CSUM2(L) = 0.0

           DO  NSC=1,NSCALE
              PR(NSC,L)    = 0.0
              AVSUM(NSC,L) = 0.0
              AMP(NSC,L)   = 0.0
           ENDDO
        ENDDO

        DO K=1,NZ
           IIK = (K - 1)
           IF (IIK > NS2)  IIK = IIK-NZ
           PK = (FLOAT(IIK) / FLOAT(NS2))**2

           DO J=1,NY
              IIJ = (J-1)
              IF (IIJ > NR2)  IIJ = IIJ - NY
              PJ = (FLOAT(IIJ)/FLOAT(NR2))**2

              DO I=1,LSD,2
                III=(I-1)/2

C               SKIP HERMITIAN RELATED VALUES
                IF (III > 0 .OR. (IIK >= 0 .AND. 
     &             (IIJ >=0 .OR. IIK.NE.0))) THEN

                PII = (FLOAT(III) / FLOAT(ND2))**2
                R   = SQRT(PII+PJ+PK) * 0.5

		IF (R  >  0.0)  THEN
                   R1 = SQRT(PII+PK)*0.5          ! r1=(x^2+z^2)^0.5
                   !if(r1 > 0.05) cycle

                   IF (SER == 'C')  THEN
                      IF (ACOS(AMIN1(1.0,SQRT(PK)/R)) < ZANG) CYCLE

                   ELSEIF (SER == 'W')  THEN
                      IF (R1 == 0.0)  THEN
                         FI = ACOS(1.0)
                      ELSE
                         FI = ACOS(AMIN1(1.0,SQRT(PK)/R1))
                      ENDIF
                      IF (FI  <  ZANG) CYCLE
                   ENDIF

                   L = NINT(R*2*(INC-1) )+ 1
                   IF (L  >  INC) CYCLE

                   LR(L) = LR(L) + 2
                   IF (A(I,J,K) .NE. 0.0)THEN
                      PHA = ATAN2(A(I+1,J,K),A(I,J,K)) * RAD_TO_DGR
                   ELSE
                      PHA = 0
                   ENDIF 

                   IF (B(I,J,K) .NE. 0.0)THEN
                      PHB = ATAN2(B(I+1,J,K),B(I,J,K)) * RAD_TO_DGR
                   ELSE
                      PHB = 0.0
                   ENDIF 

                   DPH = PHA - PHB
                   IF (DPH  >  180.0)  DPH = 360.0-DPH
                   IF (DPH  <  -180.0) DPH = 360.0+DPH

                   QA = SQRT(A(I,J,K)**2 + A(I+1,J,K)**2)
                   QB = SQRT(B(I,J,K)**2 + B(I+1,J,K)**2)

                   CSUM(L)=
     &                CSUM(L)+A(I,J,K)*B(I,J,K)+A(I+1,J,K)*B(I+1,J,K)

                   CSUM1(L) = CSUM1(L) + QA * QA
                   CSUM2(L) = CSUM2(L) + QB * QB

C                  SCALES
                   DPH = DPH * DPH

                   DO NSC=1,NSCALE
                      SCALE        = SCALE1 + (NSC-1) * DSCALE
                      QBS          = QB * SCALE
                      AVSUM(NSC,L) = AVSUM(NSC,L) + ABS(QA - QBS)
                      PR(NSC,L)    = PR(NSC,L)  + ((QA+QBS)/2.) * DPH
                      AMP(NSC,L)   = AMP(NSC,L) + ((QA+QBS)/2.)
                   ENDDO
		ENDIF
                ENDIF

             ENDDO
          ENDDO
        ENDDO

        END
