
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
C                                                                      *
C                                                                      *
C  PURPOSE:                                                            *
C                                                                      *
C CNSCAN.FOR 
C
C  LAST UPDATE  22/11/73    MRS         WRS:CONTUR.SCAN
C               4/6/87      AL          MULTIZ ADDED
C
C **********************************************************************
C
C    CNSCAN(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS)
C
C    PURPOSE:  CONTUR SUBROUTINE TO CALCULATE PLOT POINTS
C
C    PARAMETERS:
C       AM         2-D MATRIX TO BE CONTOURED. 
C       IRRX       WORKING ARRAY
C       X,Y        ARRAYS FOR CONTOUR POINTS
C       NMAX       LENGTH OF X AND Y ARRAYS
C       LUN        LOGICAL OUTPUT UNIT FOR CONTOUR FILE
C       MULITIZ    FLAG FOR MULTIPLE Z LEVELS
C       MAXPTS     MAX. NO. OF POINTS WANTED ON A STORED CONTOUR
C
C    VARIABLES:
C       MT, NT     AM ARRAYS  X AND Y DIMENSIONS
C       CV         THE CONTOUR LEVEL.
C       RA         RATIO OF THE LENGTH OF ONE DIVISION IN THE Y DIRECTION TO
C                  ONE DIVISION IN X.
C       THE        COS OF THE ANGLE BETWEEN THE X AND Y AXIS
C
C   CALLED BY:    CNTUR
C
C   CALLS:        CNTRCE
C
C--********************************************************************

      SUBROUTINE CNSCAN(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,MAXIRR)

 

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

      COMMON /UNITS/LUNDOC,NIN,NOUT

      INTEGER   IRRX(MAXIRR)
      DIMENSION AM(KAM),X(1),Y(1)
      LOGICAL   MULTIZ

C     RESETS NW COUNTER FOR IRRX TAGS
      NW = 0

      IF (ISTART .EQ. 0)  THEN
C        SET UP OFFSET POINTERS, ETC.
         ISTART   =1
         RC       = THE*RA
         RS       = SQRT(1.0-THE**2)*RA
         IPT(1,1) = 8
         IPT(1,2) = 1
         IPT(1,3) = 2
         IPT(2,1) = 7
         IPT(2,3) = 3
         IPT(3,1) = 6
         IPT(3,2) = 5
         IPT(3,3) = 4
         INX(1)   = -1
         INX(2)   = -1
         INX(3)   = 0
         INX(4)   = 1
         INX(5)   = 1
         INX(6)   = 1
         INX(7)   = 0
         INX(8)   = -1
         INY(1)   = 0
         INY(2)   = 1
         INY(3)   = +1
         INY(4)   = +1
         INY(5)   = 0
         INY(6)   = -1
         INY(7)   = -1
         INY(8)   = -1
      ENDIF

C*********FEB 92 AL
C**   3 DO 58 J=1,NMAX
C**   58 IRRX(J) = 0
C*******

      ISS = 0

C     SCAN TOP ROW
    2 M5 = MI+(NI-1)*M
      M6 = MF-1+(NI-1)*M

      DO  I=M5,M6
         IF (AM(I) .LT. CV .AND. AM(I+1) .GE. CV) THEN
            IY  = NI
            IX  = I+1-M*(NI-1)
            IDX = -1
            IDY = 0
            CALL CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,
     &                  MAXIRR,IRTFLG)
         ENDIF
      ENDDO

C     SCAN RIGHT EDGE
      N5 = (NI-1)*M+MF
      N6 = (NF-2)*M+MF
      DO  I=N5,N6,M
         IF (AM(I) .LT. CV .AND. AM(I+M) .GE. CV) THEN
            IX  = MF
            IY  = (I+M-1)/M+1
            IDX = 0
            IDY = -1
            CALL CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,
     &                  MAXIRR,IRTFLG)
         ENDIF
      ENDDO

C     SCAN BOTTOM ROW
   22 MT3 = M*(NF-1)+MF+1
      MT1 = MF-MI
      DO  I=1,MT1
         MT2 = MT3-I
         IF (AM(MT2).LT. CV .AND. AM(MT2-1) .GE. CV) THEN
           IX  = MF-I
           IY  = NF
           IDX = 1
           IDY = 0
           CALL CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,
     &                 MAXIRR,IRTFLG)
         ENDIF
      ENDDO

C     SCAN LEFT EDGE
      NT3 = (NF-1)*M+MI
      NT1 = (NF-NI)-1
      DO  I=1,NT1
         NT2 = NT3-M*(I-1)
         IF (AM(NT2) .LT. CV .AND. AM(NT2-M) .GE. CV) THEN 
           IX  = MI
           IY  = (NT2-M-1)/M+1
           IDX = 0
           IDY = 1
           CALL CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,
     &                  MAXIRR,IRTFLG)
         ENDIF
      ENDDO
      ISS=1

C     SCAN MIDDLE OF MATRIX
      NT5 = NI+1
      NT6 = NF-1
      MT6 = MF-1
C     SEARCH  EACH ROW
      DO 10 J = NT5,NT6
         K3=(J-1)*M
C        SEARCH EACH COLUMN
         DO 10 I=MI,MT6
C           FIND PIXEL POINTER
            K=K3 +I
            IF (AM(K) .LT. CV .AND. AM(K+1) .GE. CV) THEN
C              PIXEL IS AT THE CONTOUR LEVEL

    7          IF (NW .NE. 0) THEN
C                 ALREADY HAVE SOME CONTOURS ON THIS LEVEL
                  NPP = NW
                  KTEMP = K + 1

C                 CHECK TO SEE IF THIS PIXEL IS TAGGED IN IRRX YET
                  DO ID = 1,NPP
                    IF (IRRX(ID) .EQ. KTEMP)  GOTO 10 
                  ENDDO             
               ENDIF

C              THIS PIXEL IS ON A NEW CONTOUR AT THIS LEVEL
               IX  = I+1
               IY  = J
               IDX = -1
               IDY = 0
               CALL CNTRCE(AM,KAM,IRRX,X,Y,NMAX,LUN,MULTIZ,MAXPTS,
     &                     MAXIRR,IRTFLG)

               IF (IRTFLG .NE. 0) THEN
C                 OVERFLOW OF IRRX (ALREADY CONTOURED FLAGS) OR
C                 OVERFLOW OF X,Y ARRAYS.
                  WRITE(NOUT,*) 'CONTOURING ABORTED ON THIS LEVEL'
                  RETURN
               ENDIF           
            ENDIF
   10 CONTINUE

      RETURN
      END
