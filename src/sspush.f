
C++*********************************************************************
C
C    SSPUSH.F                  REWRITTEN FOR FILE BASED INPUT JAN 98 al
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
C    SSPUSH(LUNOUT,DATA,NDATA,IRTFLG)
C
C    PURPOSE:  TO WRITE THE CONTOUR DATA TO LUNSS
C
C    PARAMETERS:        
C        DATA(3,NDATA)  CO-ORDINATES OF CONTOURS                 (SENT)
C        NDATA          NO. OF POINTS                            (SENT)
C        islice         slice number                             (SENT)
C        IRTFLG         ERROR RETURN FLAG.(ZERO IS NORMAL)       (RET.)
C
C    VARIABLES:
C        ZCOO           Z- COORDINATE FOR THE CONTOUR (IN COMMON)
C        LUNSS          LOGICAL UNIT NUMBER  ASSIGNED TO STORAGE FILE.
C
C--********************************************************************

	SUBROUTINE SSPUSH(LUNOUT,DATA,NDATA,ISLICE,IRTFLG)

        INCLUDE 'CMBLOCK.INC'

        CHARACTER *4 CLASSNAME

        PARAMETER   (NSIZE = 2000)
        DIMENSION   DATA(3,NSIZE)

        COMMON /POLY/    MINPTS,ISLICET,ICNT

C       SET BAD ERROR RETURN
        IRTFLG = 1
       
        ISEC      = ISLICE
        IPLAN     = 1
        CLASSNAME = 'CONT'
        ICNT      = ICNT + 1  
        ZCOO      = ISEC
        I3D       = 2

        WRITE(LUNOUT,90)ISEC,IPLAN,CLASSNAME,ICNT,NDATA,I3D,ZCOO
   90   FORMAT(2I10,6X,A,3I10,F10.3)

        DO I = 1,NDATA
           WRITE(LUNOUT,*) DATA(1,I),DATA(2,I),ZCOO
        ENDDO

        IRTFLG = 0

        RETURN
        END

