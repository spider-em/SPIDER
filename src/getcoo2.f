C++*********************************************************************
C
C 
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

      SUBROUTINE GETCOO2 ( ITOT, JTOT, JACT, NFAC, LUV, IDROW, 
     &                     S, D, PJ, PJA, SOM, SOMA, COORD, CO, U, W, 
     &                     LSAV, LPIX, LIMA)

 

C     I DO NOT KNOW IF SAVE IS NEEDED FEB 99 al
      SAVE

      REAL    S(JACT, JACT), D(JACT), CO(JTOT), COORD(JTOT, NFAC),
     &  PJ(JTOT), PJA(JACT), U(JTOT), W(NFAC)
      INTEGER IDROW(ITOT), LUV(JTOT)
      DATA IAI / 1 /

C.......... COORDINATES AND DISTANCES TO ORIGIN FOR THE ROWS

      CALL REW ( LSAV, 1)
      CALL REW ( LPIX, 1)

C.......... TRANSFORM COLUMN WEIGHTS FROM ABSOLUTE TO RELATIVE

      DO  J = 1,JTOT
          PJ(J) = PJ(J) / SOM
          CO(J) = 0.0
      ENDDO
      DO  JA = 1,JACT
          PJA(JA) = PJA(JA) / SOMA
      ENDDO
      DO   I = 1,ITOT
          READ (LSAV)  (U(J),J=1,JTOT), PIT, IDPIXL
	  PIA = 0.0
	  DO  J=1,JTOT
	      IF (LUV(J) .EQ. 1) PIA = PIA + U(J)
   	  ENDDO
          DO   K = 1,NFAC
              W(K)  = 0.0
	      JA    = 0
              DO 20  J = 1,JTOT
		  IF (LUV(J) .NE. 1) GOTO 20
		  JA = JA + 1
		  W(K)  = W(K) + (U(J) * S(JA,K)) / PIA
   20         CONTINUE
  	  ENDDO
          DOR =  0.0
          DO 40  J = 1,JTOT
	      IF (LUV(J) .NE. 1) GOTO 40
	      DOR = DOR + (U(J) / PIT - PJ(J)) ** 2 / PJ(J)
   40     CONTINUE
          IF (DOR .LT. 1.0E-10) DOR = 1.0E-10
          PI    = PIA / SOMA
          DO   J = 1,JTOT
	      CO(J) = CO(J) + (U(J)/(SOM*PJ(J)) - PI)**2 / PI
   	  ENDDO
   60     WRITE (LPIX)  (W(K),K=1,NFAC), PI, DOR, IDPIXL, IAI
      ENDDO
C.......... COMPUTE COORDINATES OF THE ACTIVE COLUMNS
      DO   K = 1,NFAC
          IF (D(K) .LT. 1.0E-9) D(K) = 1.0E-9
C.......... TRANSFORM EIGENVALUES INTO THEIR SQUARE ROOT   
          D(K)  = SQRT(D(K))
	  JA = 0
          DO 80  J = 1,JTOT
	      IF (LUV(J) .NE. 1) GOTO 80
	      JA = JA + 1
	      COORD(J,K) = S(JA,K) * D(K)
   80     CONTINUE
      ENDDO
C.......... COMPUTE COORDINATES OF THE INACTIVE COLUMNS
C           TRANSFORM BACK TO ABSOLUTE WEIGHTS OF THE COLUMNS
      DO   J = 1,JTOT
	  PJ(J) = SOM * PJ(J)
      ENDDO
      DO 120  J = 1,JTOT
	  IF (LUV(J) .EQ. 1) GOTO 120
          DO   K  = 1,NFAC
	      COORD(J,K)= 0.0
	  ENDDO
  120 CONTINUE
      CALL REW ( LSAV, 1)
      CALL REW ( LPIX, 1)
      DO   I = 1,ITOT
          READ (LSAV)  (U(J),J=1,JTOT), PIA, IDPIXL
          READ (LPIX)  (W(K),K=1,NFAC), PIA, DOR, IDPIXL, LAI
          DO 135  J = 1,JTOT
	      IF (LUV(J) .EQ. 1) GOTO 135
              DO   K = 1,NFAC
                  COORD(J,K) = COORD(J,K) + (U(J)*W(K)) / (PJ(J)*D(K))
	      ENDDO
  135     CONTINUE
      ENDDO
C.......... GO BACK TO THE EIGENVALUES
      DO   K = 1,NFAC
       D(K)  = D(K)*D(K)
      ENDDO
C.......... WRITE THE COORDINATES OF THE COLUMNS (IMAGES)
      DO   J = 1,JTOT
	PI=PJ(J)/SOMA
      WRITE(LIMA) (COORD(J,K),K=1,NFAC), PI, CO(J), IDROW(J), LUV(J)
      ENDDO 
      RETURN
      END
