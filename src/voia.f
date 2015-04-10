C++*********************************************************************
C
C VOIA.F                          NEW               AUG 04 ARDEAN LEITH
C                                 USED surfrot.c
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
C  VOIA
C
C  PURPOSE:  DETERMINES INCLUDED ANGLE BETWEEEN 2 VECTORS
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************


       SUBROUTINE VOIA

       INCLUDE 'CMBLOCK.INC' 

       REAL, DIMENSION(3) :: VEC1,VEC2

       CALL  RDPRM3S(VEC1(1),VEC1(2),VEC1(3),NOT_USED,
     &                   'FIRST VECTOR  - X, Y & Z',IRTFLG)
       IF (IRTFLG .NE. 0)  RETURN

       CALL  RDPRM3S(VEC2(1),VEC2(2),VEC2(3),NOT_USED,
     &                   'SECOND VECTOR - X, Y & Z',IRTFLG)
       IF (IRTFLG .NE. 0)  RETURN

       CALL INCANGV(VEC1,VEC2,COSANG,ANG,IRTFLG)

       WRITE(6,90) ANG 
90     FORMAT('  INCLUDED ANGLE: ',F7.2)

       CALL REG_SET_NSEL(1,1, ANG,0.0,0.0,0.0,0.0, IRTFLG)
 
       END



        SUBROUTINE INCANGV(VEC1,VEC2,COSANG,ANG,IRTFLG)
  
        REAL                   :: VEC1(3),VEC2(3)

        REAL, PARAMETER        :: QUADPI     = 3.14159265358979323846 
        REAL, PARAMETER        :: DGR_TO_RAD = (QUADPI/180)

        REAL, PARAMETER        :: FLTZER     = 10E-30
        
        DIS1 = VEC1(1)*VEC1(1) + VEC1(2)*VEC1(2) + VEC1(3)*VEC1(3)
        DIS2 = VEC2(1)*VEC2(1) + VEC2(2)*VEC2(2) + VEC2(3)*VEC2(3)

C       FIND ANGLE BETWEEN LINES  
        DIS13 = DIS1 * DIS2

        IF (ABS(DIS13) .LT. FLTZER) THEN
C          REJECT DUPLICATED POINTS OR GET DIVISION BY ZERO!
           COSANG = 1.0
           ANG    = 0.0
           IRTFLG = 1

        ELSE 

C          COSANG =  DOT(VEC1,VEC2) / SQRT(DIS13)
           COSANG = (VEC1(1)*VEC2(1) + 
     &               VEC1(2)*VEC2(2) + 
     &               VEC1(3)*VEC2(3))  / SQRT(DIS13)

#if defined (SP_GFORTRAN)
           ANG =  ACOS(COSANG) / DGR_TO_RAD
#else
           ANG =  ACOSD(COSANG)
#endif

           IRTFLG = 0
        ENDIF

        END

