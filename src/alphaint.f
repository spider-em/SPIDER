C++*********************************************************************
C
C  AF.F                               
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
C   AF(MAXDIM,LUN1,LUN2,LUN3,FILNAM,FILNAMO,NSAM,NROW,NSLICE,
C         MAXIM,IRTFLG)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*******************************************************************      DOUBLE PRECISION FUNCTION  ALPHAINT(F,N1,N2)
 
      DOUBLE PRECISION FUNCTION  ALPHAINT(F,N1,N2)

      IMPLICIT REAL*8  (A-H,O-Z)
     
      INTEGER*4 LUN51

      DATA LUN51/28/

      XN1=N1
      XN2=N2
      XM=XN2/(XN2+XN1*F)
      CALL  INVBT(XM,XN2/2.0,XN1/2.0,IFAULT,HEJ)
      IF  (IFAULT.NE.0)  THEN
         WRITE(LUN51,1)  IFAULT
 1       FORMAT(//'  E R R O R   I N   A L P H A I N T  ',I1)
         CLOSE(LUN51)
	 ALPHAINT=0.0
         STOP
      ENDIF
      ALPHAINT=HEJ
      END
