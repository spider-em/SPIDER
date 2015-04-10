
C ++********************************************************************
C
C POLLSQ.F                             REDUCED ORDER MAY 99 ARDEAN LEITH
C                                      USED REG_SET AUG 00 AL
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
C  PROGRAM TO FIT A POLYNOMIAL TO NUMBERS IN A DOCUMENT FILE READ
C  BEFORE CALLING THIS PROGRAM WITH UNSDAL INTO DBUF.
C  MAXIMUM ORDER OF POLYNOM IS 10  10/14/87 MR.
C
C  REDUCED ORDER TO 5 DUE TO INABILITY OF ALGORITHM TO COPE WITH
C  HIGHER ORDER DATA IN SOME CIRCUMSTANCES MAY 99 al
C                                                                      
C***********************************************************************

	SUBROUTINE POLLSQ(DBUF,MAXREG,MAXKEY,NDOC)

        INCLUDE 'CMBLOCK.INC'

        COMMON     BUF(1)
        CHARACTER  YN,NULL
        DIMENSION  REG(3)
        DIMENSION  DBUF(MAXREG,MAXKEY),C(11)

        NULL = CHAR(0)

        CALL RDPRMI(KEY1,KEY2,NOT_USED,'FIRST, LAST KEY NUMBER')

	CALL RDPRMI(IABS,IORD,NOT_USED,
     &             'COLUMN #S FOR ABSCISSA, ORDINATE (0=KEY)')

        CALL RDPRM2(SCALEX,SCALEY,NOT_USED,'SCALEX, SCALEY (DEF=1)')
        IF (SCALEX .LE. 0.0) SCALEX = 1.0
        IF (SCALEY .LE. 0.0) SCALEY = 1.0

        NORD = 1
        CALL RDPRI1S(NORD,NOT_USED,
     &              'ORDER OF POLYNOMIAL (1-5, DEFAULT= 1)',IRTFLG)
        IF (NORD .LE. 0) NORD = 1
        IF (NORD .GT. 5) THEN
           CALL ERRT(31,'POLLSQ',IDUM)
           GOTO 9999
        ENDIF

C       INCREASE IABS & IORD BY 1 BECAUSE DBUF(1,*) CONTAINS KEY, 
C       FOLLOWED BY REGISTERS 
        IABS = IABS + 1
        IORD = IORD + 1

C       STORE COORDINATES INTO ARRAYS. X IN BUF(1)... , 
C                                      Y IN BUF(IOFF+1)...

        ILINE = 0                          
        IDO   = 1  
        IOFF  = KEY2 - KEY1 + 1
        DO  I=KEY1,KEY2        
           BUF(IDO)      = DBUF(IABS,I) * SCALEX
           BUF(IDO+IOFF) = DBUF(IORD,I) * SCALEY
           IDO           = IDO + 1
	ENDDO

C       WRITE OUT INPUT VALUES
        WRITE(NOUT,101) (BUF(L),BUF(IOFF+L),L=1,IOFF)
101     FORMAT(' ',2E12.4)

C       POLFIT DOES THE LEAST SQUARE FIT, INCLUDING THE CREATION OF THE
C       MATRIX SYSTEM. RESULT IS RETURNED IN C(1)... C(NORD+1). 
C       IOFF IS DIMENSION  OF ARRAYS CONTAINING X AND Y
 
        CALL POLFIT(BUF(1),BUF(IOFF+1),NORD,IOFF,C)
 
        IF (NORD+1 .LT. 11) THEN
C          ZERO UNUSED C VALUES
           DO K = NORD+2,11
              C(K) = 0.0
           ENDDO
        ENDIF

C       STORE RESULT INTO REGISTERS (MAXIMUM 6 I.E. POLYNOM ORDER 5)
        ILAST = MAX(6,NORD+1)
        CALL REG_SET_NSELA(ILAST,C,.TRUE.,IRTFLG)
        
C       WRITE OUT RESULTS
        WRITE(NOUT,100) (C(I),I=1,NORD+1) 
100     FORMAT('  RESULTING POLYNOM:',/,1X,              
     &   E12.4,' + '
     &  ,E12.4,' *X + '
     &  ,E12.4,' *X**2 + '
     &  ,E12.4,' *X**3 + '
     &  ,E12.4,' *X**4 + '
     &  ,E12.4,' *X**5 + ',/)

C       1X
C    &  ,E12.4,' *X**6 + '
C    &  ,E12.4,' *X**7 + '
C    &  ,E12.4,' *X**8 + '
C    &  ,E12.4,' *X**9 + '
C    &  ,E12.4,' *X**10 + ')
 
C       CURRENTLY A DOCUMENT FILE CAN BE CREATED, THAT CAN BE PLOTTED LATER.
        CALL RDPRMC(YN,NCHAR,.TRUE.,
     &     'CREATE OUTPUT DOCUMENT FILE (Y/N)',NULL,IRTFLG)
        IF (YN .NE. 'Y') RETURN

        CALL RDPRM2(XA,XE,NOT_USED,'X FROM, TO')
        CALL RDPRM(XIN,NOT_USED,'IN INTERVALS OF')

        NPOINT = (XE - XA) / XIN + 1.0
        X      = XA
        DO  I=1,NPOINT
           POLY = C(1)
           DO  K=1,NORD
              POLY = POLY + C(K+1) * X ** K
	   ENDDO
 
C          OUTPUT INTO DOCUMENT FILE:
           REG(1) = I
           REG(2) = X
           REG(3) = POLY
           CALL SAVD(NDOC,REG,3,IRTFLG)
           X = X + XIN
	ENDDO
       
9999    CALL SAVDC
        CLOSE(NDOC)

        RETURN
        END

