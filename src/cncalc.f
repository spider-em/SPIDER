
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
C
C  CNCALC.F 
C
C **********************************************************************
C
C  LAST UPDATE    11/20/89 al         
C  PREVIOUS UPDATES    01/08/78 WRS     13/12/74
C
C **********************************************************************
C
C   CNCALC(AM,KAM,X,Y,NMAX)
C
C   PURPOSE:  CONTUR SUBROUTINE TO CALCULATE PLOT POINTS
C
C   PARAMETERS:   AM    2-D ARRAY FOR THIS SLICE OF IMAGE
C                 KAM   NUMBER OF PIXELS IN AM
C                 X,Y   ARRAYS FOR COORDINATES FOR THIS CONTOUR
C                 NMAX  DIMENSION OF X AND Y
C  
C   CALLED BY:    CNTRCE
C
C   CALLS:        NONE
C
C--********************************************************************

      SUBROUTINE CNCALC(AM,KAM,X,Y,NMAX,IRTFLG)

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

      INCLUDE 'CMLIMIT.INC'
      COMMON  /UNITS/LUNC,NIN,NOUT

      DIMENSION        AM(KAM),X(NMAX),Y(NMAX)

      PARAMETER        (NSIZE = 2000)
      DIMENSION        WORK(3,NSIZE),WORK1(3,NSIZE)
      EQUIVALENCE      (BUF(1),WORK(1,1)), (BUF(6400),WORK1(1,1))
      COMMON /IOBUF/   BUF(NBUFSIZ)

      DATA      FMAX/1.E-10/

      K1  = (IY-1)*M+IX
      IT  = 0
      N   = N+1

      IF (IDX**2 + IDY**2 -1 .EQ. 0)  THEN
        IF (IDX .EQ. 0) THEN
            X(N) = IX
            Z   = IY
            DY  = IDY
            K2  = K1+IDY*M
            DVS = AM(K1) - AM(K2)
            IF (ABS(DVS) .LT. FMAX) DVS = FMAX
            Y(N) = ((AM(K1)-CV)/DVS)*DY+Z
         ELSE
            Y(N) = IY
            W    = IX
            DX   = IDX
            K2   = K1+IDX
            DVS  = AM(K1) - AM(K2)
            IF (ABS(DVS) .LT. FMAX) DVS = FMAX
            X(N) = ((AM(K1)-CV)/DVS)*DX+W
         ENDIF

      ELSE
         W   = IX
         Z   = IY
         DX  = IDX
         DY  = IDY
         K2  = K1+IDX
         K3  = K1+IDY*M
         K4  = K3+IDX
         DCP = (AM(K1)+AM(K2)+AM(K3)+AM(K4))/4.0

         IF (PY .EQ. 2.0 .OR. DCP .LE. CV) THEN
           AL   = AM(K1)-DCP
           IF (ABS(AL) .LT. FMAX) AL = FMAX
           V    = .5*(AL+DCP -CV)/AL
           X(N) = V*DX+W
           Y(N) = V*DY+Z
           PY   = 0.0
         ELSE

           IT   = 1
           AL   = AM(K4)-DCP
           IF (ABS(AL) .LT. FMAX) AL = FMAX
           V    = .5*(AL+DCP-CV)/AL
           X(N) = -V*DX+W + DX
           Y(N) = -V*DY+Z  + DY
         ENDIF
      ENDIF

      IRTFLG = 0

      IF (N .GE. (NMAX-1)) THEN 
C       ARRAY ALMOST FULL, TRY TO THIN THE CONTOUR
        DO  I = 1,N
          WORK(1,I) = X(I)
          WORK(2,I) = Y(I)
        ENDDO

        N0     = N
C       ANGLE FOR RETAINING
        RETANG = 178.0
C       NUMBER OF POINTS TO TRY TO REDUCE TO
        MAXPTS = 1400
C       NUMBER OF ITERATIONS
        MAXTRY = 1
        CALL RATHIN(WORK,N,0.0333,RETANG,MAXPTS,MAXTRY,WORK1,IRTFLG)

        WRITE(NOUT,*) 'CONTOUR LENGTH REDUCED FROM: ',N0,
     &                ' TO: ',N,' POINTS.'

        DO  I = 1,N
           X(I) = WORK(1,I) 
           Y(I) = WORK(2,I) 
        ENDDO

        IRTFLG = 1
      ENDIF

      RETURN
      END
