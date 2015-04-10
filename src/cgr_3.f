
C++*********************************************************************
C
C  CGR_3.F
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
C  CALCULATES CENTER OF GRAVITY INSIDE ELLIPSOIDE AROUND HIGHEST PEAK
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE CGR_3(LUN,BUF,NSAM,NROW,NSLICE,
     &                    ELIPX,ELIPY,ELIPZ,NPC,
     &                    RXCEN,RYCEN,RZCEN,RSUM)

         DIMENSION  BUF(NSAM),NPC(3)

         XCOOR=NPC(1)
         YCOOR=NPC(2)
         ZCOOR=NPC(3)

         XMAXI = XCOOR + ELIPX
         XMINI = XCOOR - ELIPX
         YMAXI = YCOOR + ELIPY
         YMINI = YCOOR - ELIPY
         ZMAXI = ZCOOR + ELIPZ
         ZMINI = ZCOOR - ELIPZ

         RSUM  = 0

         RXCEN = 0
         RYCEN = 0
         RZCEN = 0

         DO I=ZMINI,ZMAXI
            RRAUZI = ((FLOAT(I)-ZCOOR)/ELIPZ)**2
            IR     = MOD(I-1+NSLICE,NSLICE)+1

            DO II=YMINI,YMAXI
               IIR = MOD(II-1+NROW,NROW)+1
               CALL REDLIN(LUN,BUF,NSAM,(IR-1)*NROW+IIR)

               RRAUXI = ((FLOAT(II)-YCOOR)/ELIPY)**2 + RRAUZI

               DO III=XMINI,XMAXI
                  RELIPS = ((FLOAT(III)-XCOOR)/ELIPX)**2 + RRAUXI

                  IF (RELIPS .LE. 1.0) THEN
                     IIIR = MOD(III-1+NSAM,NSAM) + 1
                     IF (BUF(IIIR) .LT. 0.0)  THEN
                        RSUM = 0.0
                        RETURN
                     ENDIF
                     RXCEN = BUF(IIIR)*III + RXCEN
                     RYCEN = BUF(IIIR)*II  + RYCEN
                     RZCEN = BUF(IIIR)*I   + RZCEN
                     RSUM  = BUF(IIIR)+RSUM
                  ENDIF
	       ENDDO
	    ENDDO
	 ENDDO
         IF (RSUM .EQ. 0.0) RETURN

         RXCEN = RXCEN/RSUM
         RYCEN = RYCEN/RSUM
         RZCEN = RZCEN/RSUM

         END
