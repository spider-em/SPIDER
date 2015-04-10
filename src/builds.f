
C++*********************************************************************
C
C  BUILDS.F
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
C  BUILDS
C
C  PURPOSE: BULID ROTATION MATRICES FROM THREE EULERIAN ANGLES
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE  BUILDS(DS,NANG,ANGSYM,IRTFLG)

        INCLUDE    'CMBLOCK.INC' 

        DIMENSION   DS(9,NANG),ANGSYM(4,NANG)

         DO K=1,NANG

C          READ ANGLES FROM THE DOCUMENT FILE.
C          ORDER IN THE DOCUMENT FILE IS PSI, THETA, PHI AND ANGLES 
C          ARE IN DEGREES! IN ANG ARRAY IT IS OTHER WAY AROUND

           ICOUNT = ANGSYM(1,K)
           IF (ICOUNT .LE. 0) THEN
C             MISSING KEY
              CALL ERRT(102,'MISSING SYMMETRY',K)
              IRTFLG = 1
              RETURN
           ENDIF

           CALL CANG(ANGSYM(4,K),ANGSYM(3,K),ANGSYM(2,K),
     &               .FALSE.,SSDUM,DS(1,K))

           WRITE(NOUT,333)  K,(ANGSYM(J,K),J=2,4)
333        FORMAT(' SYMMETRY #',I5,
     &         '; PSI=',F6.1,' THETA=',F6.1,' PHI=',F6.1)
        ENDDO

        IRTFLG = 0

        END


