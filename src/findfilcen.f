C++*******************************************************************
C
C FINDFILCEN.F  NEW                               DEC 14 ARDEAN LEITH
C
C **********************************************************************
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
C  PURPOSE: FIND SPIDER FILE CENTER USING SPIDER'S CENTER CONVENTION
C
C  FINDFILCEN(ITYPE,NX,NY,NZ,IXC,IYC,IZC)
C
C  PARAMETERS: ITYPE         FILE TYPE
C	       NX,NY,NZ      FILE DIMENSIONS
C	       IXC,IYC,IZ    FILE CENTER
C	       IRTFLG        ERROR FLAG
C
C--*******************************************************************
        
       SUBROUTINE FINDFILCEN(ITYPE,NX,NY,NZ,IXC,IYC,IZC,IRTFLG)

       IMPLICIT NONE

       INTEGER  :: ITYPE,NX,NY,NZ,IXC,IYC,IZC,IRTFLG

       INTEGER  :: IDUM

       IRTFLG = 0

       IF (ITYPE == 1 .OR. ITYPE == 3) THEN
   
          IXC = NX / 2 + 1

          IYC = 0
          IF (NY > 0) IYC = NY / 2 + 1 

          IZC = 0
          IF (NZ > 0) IZC = NZ / 2 + 1 

          IRTFLG = 0

       ELSE

          CALL ERRT(102,'PGM ERROR, FOURIER INPUT NOT SUPPORTED',IDUM)
          IRTFLG = 1
       ENDIF

       END

