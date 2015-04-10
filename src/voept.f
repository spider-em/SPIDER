
C++*********************************************************************
C
C VOEPT.F                         NEW                AUG 04 ARDEAN LEITH
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
C  VOEPT
C
C  PURPOSE:  FINDS POINT ON EULER ANGLE VECTOR
C
C  NOTES:    I KNOW THE ORDER IS FLOOEY, SOMEWHERE THERE IS A DIFFERENCE
C            IN NAMING CONVENTION
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

       SUBROUTINE VOEPT

       INCLUDE 'CMBLOCK.INC' 

       XI    = 0
       YI    = 0
       ZI    = 1
       CALL  RDPRM3S(XI,YI,ZI,NOT_USED,
     &              'INITIAL POINT - X, Y & Z',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL  RDPRM3S(PHI,THETA,PSI,NOT_USED,
     &              'EULER ANGLES - PHI, THETA & PSI',IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL ANGD(PSI,THETA,PHI, XI,YI,ZI, X,Y,Z, IRTFLG)

       WRITE(6,90) PHI,THETA,PSI, XI,YI,ZI, X,Y,Z 
90     FORMAT('  PHI: ',F5.1,' THETA: ',F5.1,' PSI: ',F5.1,' : ',
     &        '(',F8.3,',',F8.3,',',F8.3,')',' --> ',
     &        '(',F8.3,',',F8.3,',',F8.3,')')

       CALL REG_SET_NSEL(1,3, X,Y,Z,0.0,0.0, IRTFLG)
 
       END



      SUBROUTINE ANGD(PHI1,THETA1,PSI1, XI,YI,ZI, X,Y,Z, IRTFLG)

      PARAMETER (QUADPI = 3.1415926535897932384626)
      PARAMETER (DGR_TO_RAD = (QUADPI/180))

      PHI   = PHI1   * DGR_TO_RAD
      THETA = THETA1 * DGR_TO_RAD
      PSI   = PSI1   * DGR_TO_RAD 

      RM11  =  COS(THETA)*COS(PHI)*COS(PSI)-SIN(PHI)*SIN(PSI)
      RM21  = -COS(THETA)*COS(PHI)*SIN(PSI)-SIN(PHI)*COS(PSI)
      RM31  =  SIN(THETA)*COS(PHI)

      RM12  =  COS(THETA)*SIN(PHI)*COS(PSI)+COS(PHI)*SIN(PSI)
      RM22  = -COS(THETA)*SIN(PHI)*SIN(PSI)+COS(PHI)*COS(PSI)
      RM32  =  SIN(THETA)*SIN(PHI)

      RM13  = -SIN(THETA)*COS(PSI)
      RM23  =  SIN(THETA)*SIN(PSI)
      RM33  =  COS(THETA)

      X     = RM11 * XI + RM21 * YI + RM31 * ZI   
      Y     = RM12 * XI + RM22 * YI + RM32 * ZI   
      Z     = RM13 * XI + RM23 * YI + RM33 * ZI   

CC      WRITE(6,90) X,Y,Z, VEC1(1),VEC1(2),VEC1(3)
cc90    FORMAT('(',F4.1,','F4.1,','F4.1,')....(',F4.1,','F4.1,','F4.1,')')


      END




   
