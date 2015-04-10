
C ++********************************************************************
C                                                                      *
C TRINTER.F        M.RADERMACHER, JUNE 1984
C                                                                      *
C                                                                      *
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
C                                                                      *
C  TRINTER                                                                    *
C                                                                      *
C  PURPOSE: DO A BILINEAR INTERPOLATION USING ONLY 
C           THREE POINTS OF A FACET.
C           IN PARTS BASED ON SUBROUTINE GENINT, AUTHOR: R.SMITH 1983
C------------------------------------------------------------------------------
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************


        SUBROUTINE TRINTER(LUNO,LUNN,
     $              OLDIMX,OLDIMY,OLDIMZ,NEWDIMX,NEWDIMY,NEWDIMZ)

        COMMON BUF(1)

        INTEGER OLDIMX,OLDIMY,OLDIMZ,ZCOLD,XCOLD,YCOLD,YLINE,XLINE
        INTEGER YLINEO,XLINEO,YLINE2,XLINE2

        NS1=1
        NS2=OLDIMX+1
        NS3=NS2+OLDIMX

	!WRITE(6,100)LUNO,LUNN,OLDIMX,OLDIMY,OLDIMZ,NEWDIMX,NEWDIMY,
C                    NEWDIMZ,NS1,NS2,NS3
100     FORMAT(' ',11I5)

        FACZ=FLOAT(OLDIMZ)/FLOAT(NEWDIMZ)
        FACX=FLOAT(OLDIMX)/FLOAT(NEWDIMX)
        FACY=FLOAT(OLDIMY)/FLOAT(NEWDIMY)

        DO  I=1,NEWDIMZ
        ZCOO=FLOAT(I)*FACZ
        ZCOLD=INT(ZCOO)
        ZDIF=AMOD(ZCOO,1.)
        ZLINE=(ZCOLD-1)*OLDIMY
        DO  K=1,NEWDIMY
        YLINEO=YLINE
        YLINEO2=YLINE2
        YCOO=FLOAT(K)*FACY
        YCOLD=INT(YCOO)
        YLINE=ZLINE+YCOLD
        YLINE2=YLINE+1
        IF(YLINE2.GT.(ZLINE+OLDIMY)) YLINE2=YLINE
        IF(YLINE.EQ.ZLINE) YLINE=YLINE2
        YDIF=AMOD(YCOO,1.)
        IF(YLINE.NE.YLINEO) THEN 
            
          IF(YLINE.EQ.YLINEO2) THEN
            
            NS1P=NS1
            NS1=NS2
            NS2=NS1P
            CALL REDLIN(LUNO,BUF(NS2),OLDIMX,YLINE2) 
          ELSE
            CALL REDLIN(LUNO,BUF(NS1),OLDIMX,YLINE)
            CALL REDLIN(LUNO,BUF(NS2),OLDIMX,YLINE2)        
          ENDIF
        ENDIF

        DO  L=1,NEWDIMX
        XCOO=FLOAT(L)*FACX
        XCOLD=INT(XCOO)       
C       FIRST SOLVE THE TWODIMENSIONAL PROBLEM:
        XDIF=AMOD(XCOO,1.)
        IC1=NS1+XCOLD
        IC2=IC1+1
        IC3=NS2+XCOLD
        IC4=IC3+1
        Z1=BUF(IC1)
        Z2=BUF(IC2)
        Z3=BUF(IC3)
        Z4=BUF(IC4)
        Z5=Z2-Z1
        Z6=Z4-Z3
        Z7=Z4-Z2
        Z8=Z1-Z3
        IF(ABS(Z1-Z4)-ABS(Z2-Z3).LT.0) GOTO 21
        IF(XDIF+YDIF.GT.1) GOTO 22
        Z=Z1+XDIF*Z5-YDIF*Z8
        GOTO 200

22      Z=Z4+(XDIF-1.)*Z6+(YDIF-1)*Z7
        GOTO 200

21      IF(XDIF.GT.YDIF) GOTO 23
        Z=Z3+Z8*(1-YDIF)+Z6*XDIF
        GOTO 200

23      Z=Z2+Z5*(XDIF-1.)+Z7*YDIF
200     BUF(NS3+L)=Z
	ENDDO

C       WRITE(6,100)K
        CALL WRTLIN(LUNN,BUF(NS3),NEWDIMX,K)

	ENDDO
	ENDDO

        END
 
 
