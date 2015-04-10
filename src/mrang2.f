
C ++********************************************************************
C                                                                      *
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
C  MRANG2(RPT,VPT,IVIEW,ANG,PTACTIVE,NTPT)                                                                    *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

      SUBROUTINE MRANG2(RPT,VPT,IVIEW,ANG,PTACTIVE,NTPT)

      PARAMETER (LV=300)
      PARAMETER (LS=256)

      LOGICAL           PTACTIVE(LS,LV)
      DIMENSION         RPT(2,LS),VPT(2,LS)

      DOUBLE PRECISION  SUM1,SUM2

      SUM1 = 0.0
      SUM2 = 0.0

      DO  IPT=1,NTPT
        IF (PTACTIVE(IPT,IVIEW))  THEN
          SUM1 = SUM1 + RPT(1,IPT) * DBLE(VPT(1,IPT))
     &	              + RPT(2,IPT) * DBLE(VPT(2,IPT))
          SUM2 = SUM2 + RPT(1,IPT) * DBLE(VPT(2,IPT))
     &	              - RPT(2,IPT) * DBLE(VPT(1,IPT))
	ENDIF
      ENDDO

C     FOLLOWING LINE DOESN'T HAVE TO BE IN DOUBLE PRECISION,
C     BUT IN SINGLE PRECISION COMPILER MAKES OPTIMIZATION ERROR.
      ANG = DATAN2(SUM2,SUM1)

      END
