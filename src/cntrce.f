
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

C++*************************************************************************
C
C      CNTRCE.FOR 
C
C      LAST UPDATE        11/20/89 al
C      PREVIOUS UPDATES   03/2/88  al        01/08/78      22/11/73
C                         4/6/87   al        3/25/85 al
C **************************************************************************
C
C    CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS)
C
C    PURPOSE:      CONTUR SUBROUTINE TO FOLLOW CURVE
C
C    PARAMETERS:   AM       2D ARRAY FOR THIS IMAGE
C                  IRRX     WORKING ARRAY
C                  X,Y      COORDINATES OF POINTS ON THIS CONTOUR
C                  NMAX     DIMENSION OF X,Y ARRAYS
C                  LUN      LOGICAL UNIT FOR CONTOUR FILE
C                  MULTIZ   LOGICAL FLAG FOR MULTIPLE Z LEVELS
C                  MAXPTS   MAX. NUMBER OF POINTS ON A CONTOUR
C
C    CALLED BY:    CNSCAN
C
C    CALLS:        CNCALC    CNSTUFF
C
C--********************************************************************

      SUBROUTINE CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,
     &                  MAXIRR,IRTFLGR)


C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

C-------- START OF EM-PLOTT-COMMON-------------------------------------
C     INTEGERS
      COMMON /CONT1/ ICALL, IDIDIT, IDONE, IDX, IDY, ILINE, INTT,
     &               IRCD, ISS, ISTART, ISUM1, ISUM2, ISUM3, IT, IV, 
     &               IXX1, IXX2, IXX3, IX, IY, JSUM1, JSUM2, JSUM3, JT,
     &               LEVEL, LW, M, MF, MI, MT, N, NDIV, NF, NI, NT, NW

C     FLOATING POINT
      COMMON /CONT2/ APDIV, APDIVX, CV, DL, PY, RA, RC, RS, SKALE, THE,
     &               SX, SY, DENSL

C     ARRAYS
      COMMON /CONT3/ INCX(3), IORGX(3), INX(8),
     &               INY(8),  IPT(3,3), IMAP(12), NG(3), NP(3)

      COMMON /CONT4/ CTRI(6),FCTR(6),CTRDEL(6),ICNDSH(6),ICNCEL

C--------END OF EM-PLOTT-COMMON----------------------------------------

      COMMON /UNITS/  LUNC,NIN,NOUT,NECHO,IFOUND,NPROC,NDAT

      INTEGER        IRRX(MAXIRR)
      DIMENSION      AM(KAM),X(NMAX),Y(NMAX)
      LOGICAL        MULTIZ

      PARAMETER      (NSIZE=2000)

      IRTFLGR = 0

      PY     = 0.0
      JT     = 0
      N      = 0
      IT     = 0
      IX0    = IX
      IY0    = IY
      ISX    = IDX + 2
      ISY    = IDY + 2
      IS     = IPT(ISX,ISY)
      IS0    = IS
      IF (IS0 .GT. 8) IS0 = IS0 - 8

C ****************** MAIN CONTOUR PIXEL VISITION LOOP ***************

    5 CALL CNCALC(AM,KAM,X,Y,NMAX,IRTFLGD)
      IF (N .GE. NSIZE) THEN
         CALL CNSTUFF(LUN,X,Y,N,MULTIZ,MAXPTS,1,IRTFLG)
         X(1) = X(N)
         Y(1) = Y(N)
         N    = 1
      ENDIF
              
C**      IF (IT+JT-1 .GT. 0)  THEN
      IF (IT+JT-1 .GT. 0 .AND. N .GT. 1)  THEN
C        EXCHANGE LAST TWO POINTS IN CONTOUR
         XS     = X(N-1)
         YS     = Y(N-1)
         X(N-1) = X(N)
         Y(N-1) = Y(N)
         X(N)   = XS
         Y(N)   = YS
      ENDIF

      IS  = IS + 1
      JT  = IT

9     IF (IS .GE. 9) IS = IS - 8
      IDX = INX(IS)
      IDY = INY(IS)
      IX2 = IX + IDX
      IY2 = IY + IDY

C     CHECK TO SEE IF BACK AT STARTING PIXEL OF THIS CONTOUR
      IF (ISS .GT. 0   .AND.
     &     IX .EQ. IX0 .AND. IY .EQ. IY0 .AND. IS .EQ. IS0) THEN
         CALL CNCALC(AM,KAM,X,Y,NMAX,IRTFLGD)
         IPEN = 6
         GO TO 73
      ENDIF

C     CHECK TO SEE IF END OF IMAGE SCANNING AREA
      IF ((IX2-MI+1) .EQ. 0 .OR. IX2 .GT. MF .OR. (IY2-NI+1) .EQ. 0 .OR.
     &    IY2 .GT. NF) THEN
          IPEN = 5
          GOTO 73
      ENDIF 
      K = (IY2-1)*M+IX2

      IF (CV .GT. AM(K)) GOTO 5

      IF (IDX**2+IDY**2-1 .NE. 0) THEN
         K1  = (IY-1)*M + IX
         K2  = K1 + IDX
         K3  = (IY2-1)*M + IX
         DCP = (AM(K1) + AM(K2) + AM(K3) + AM(K )) / 4.0
         IF (DCP .LT. CV) GOTO 5
         IF (INX(IS-1) .NE. 0) THEN
            IX  = IX + IDX
            IDX = -IDX
            PY  = 2.0
            CALL CNCALC(AM,KAM,X,Y,NMAX,IRTFLGD)
            IF (N .GE. NSIZE) THEN
              CALL CNSTUFF(LUN,X,Y,N,MULTIZ,MAXPTS,2,IRTFLG)
              X(1) = X(N)
              Y(1) = Y(N)
              N    = 1
            ENDIF
            IX  = IX + IDX
         ELSE
            IY  = IY + IDY
            IDY = -IDY
            PY  = 2.0
            CALL CNCALC(AM,KAM,X,Y,NMAX,IRTFLGD)
            IF (N .GE. NSIZE) THEN
              CALL CNSTUFF(LUN,X,Y,N,MULTIZ,MAXPTS,3,IRTFLG)
              X(1) = X(N)
              Y(1) = Y(N)
              N    = 1
            ENDIF
            IY  = IY + IDY
         ENDIF
      ENDIF

    6 K = (IY-1)*M + IX - 1
      IF (AM(K) .LT. CV) THEN
        IF (NPP .GE. MAXIRR) THEN
           WRITE(NOUT,*) ' *** CNTRCE, IRRX OVERFLOW,NPP:',NPP
           IRTFLGR = 1
           IPEN    = 8
           GOTO 73
        ENDIF

C       RECORD PIXEL AS VISITED IN IRRX
        NW  = NW + 1
        NPP = NW
        IRRX(NPP) = K+1
      ENDIF

      IS = IS + 5
      IX = IX2
      IY = IY2
      GO TO 9

C ************************** END OF MAIN PIXEL VISITING LOOP **********

 73   CALL CNSTUFF(LUN,X,Y,N,MULTIZ,MAXPTS,IPEN,IRTFLG)
       
      RETURN
      END
