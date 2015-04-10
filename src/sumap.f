
C++*********************************************************************
C
C $$ SUMAP.FOR
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
C
C IMAGE_PROCESSING_ROUTINE
C
C        1         2         3         4         5         6         7
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************


        SUBROUTINE SUMAP(
     &  ANGOLD,XSHOLD,YSHOLD,ANGNEW,XSHNEW,YSHNEW,ANGSUM,XSHSUM,YSHSUM)
        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
        PARAMETER (DGR_TO_RAD = (QUADPI/180))

        C=COS(ANGNEW*DGR_TO_RAD)
        S=-SIN(ANGNEW*DGR_TO_RAD)
        XSHSUM=XSHNEW+XSHOLD*C-YSHOLD*S
        YSHSUM=YSHNEW+XSHOLD*S+YSHOLD*C
        ANGSUM=AMOD(ANGNEW+ANGOLD,360.0)
        END

