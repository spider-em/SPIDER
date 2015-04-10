C++*********************************************************************
C
C  BFACT.F     ORIGINAL AUTHOR:  F.J. Asturias              SEPT. 2001
C
C************************************************************************
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
C
C  BFACT(LUN1,LUN2)
C
C  PARAMETERS:
C        LUN1,LUN2         LOGICAL UNIT NUMBERS 
C
C  NOTE: THIS SUBROUTINE IS USED BY 'FF' WHERE
C        'FF' (ffilts.f)            SETS X1 TO: (NX  / 2)**2 BUT
C        'FQ' (four_fq.f or fq_q.f) SETS X1 TO: (NXF / 2)**2 WHERE
C             NX  IS X DIMENSION OF POSSIBLY PADDED IMAGE
C             NXF IS SLIGHTLY LARGER DUE TO MIXED RADIX FOURIER PAD
C        SO THEY GIVE SLIGHTLY DIFFERENT RESULTS.  I SUSPECT THAT
C        'FF' IS ACTUALLY CORRECT?
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE BFACT(LUN1,LUN2,NX,NY,NZ,NXO)

         INCLUDE 'CMBLOCK.INC'

         COMMON        B(1)
         COMPLEX       CB(1)
         EQUIVALENCE   (CB,B)

         WRITE(NOUT,*)
     &     ' NORMALIZES AMPLITUDES BY A "B" TEMPERATURE FACTOR'
         WRITE(NOUT,*)' AMP = AMP*D(EXP(Bs**2))'

         CALL RDPRM1S(PARM_B,     NOT_USED,'B FACTOR',   IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CALL RDPRM1S(PARM_D,     NOT_USED,'D CONSTANT', IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         CALL RDPRM1S(PARM_CUTOFF,NOT_USED,'FQ CUTOFF',  IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         FQCSQ = PARM_CUTOFF**2

         NS2   = NX/2
         NR2   = NY/2
         NL2   = NZ/2

         X1    = FLOAT(NXO/2)**2
         Y1    = FLOAT(NR2)**2

         IF (NZ == 1) THEN
            Z1 = 1.0
         ELSE
            Z1 = FLOAT(NL2)**2
         ENDIF

         DO K=1,NZ
            IZ = K-1
            IF (IZ > NL2) IZ = IZ - NZ

            DO J=1,NY
               IY = J - 1
               IF (IY > NR2) IY = IY - NY
               NR = J + (K - 1) * NY

               CALL REDLIN(LUN1,B,NX,NR )

               DO I=1,NS2
                  IX = I - 1

C                 B-FACTOR NORMALIZATION
                  FSQ = (FLOAT(IX*IX) / X1 + 
     &                   FLOAT(IY*IY) / Y1 + 
     &                   FLOAT(IZ*IZ) / Z1) * 0.25

                  CORR = FSQ * PARM_B

                  IF (FSQ < FQCSQ)  THEN
                     CB(I) = CB(I) * PARM_D * EXP(CORR)
                  ELSE
                     CB(I) = CB(I)
                  ENDIF
               ENDDO

               CALL  WRTLIN(LUN2,B,NX,NR)
            ENDDO
         ENDDO

         END
