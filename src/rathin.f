
C **********************************************************************
C
C RATHIN.FOR  -- CREATED OCT 87
C                   
C **********************************************************************
C *  AUTHOR: ArDean Leith                                                  *
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
C        RATHIN(DATA,NDATA,SCALE,MAXPTS,DATAT,IRTFLG)
C
C        PURPOSE:   THINS OUT POINTS ON A CONTOUR UNTIL IT REACHES MAXPTS POINTS
C
C        CALLED BY:  SSUPER
C
C        CALLS:   
C
C        PARAMETERS:     DATA   = ARRAY FOR COORD. INPUT
C                        NDATA  = NUMBER OF POINTS ON CNT FILE LOOP
C                        SCALE  = SCALE FACTOR FOR COORD.
C                        RETANGT= ANGLE BETWEEN LINE SEGS.
C                        MAXPTS = MAX LENGTH DESIRED FOR THINNED CONTOUR
C                        MAXTRYS= MAX. THINNING ITERATIONS 
C                        DATAT  = TEMPORARY ARRAY FOR COORD. OUTPUT
C                                 (PASSED TO USE COMMON ARRAY SPACE ONLY)
C                        IRTFLG = ERROR FLAG
C                                 (1 IS NOT THINNED TO MAXPTS)
C
C--********************************************************************


        SUBROUTINE RATHIN(DATA,NDATA,SCALE,RETANGT,MAXPTS,MAXTRYS,
     &                    DATAT,IRTFLG)
        COMMON /UNITS/ LUNDOC,NIN,NOUT
        
        PARAMETER (NSIZE  = 2000)
C..     MAXIMUM OF 2000 POINTS / CONTOUR

        DIMENSION      DATA(3,NSIZE), DATAT(3,NSIZE)
	PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
	PARAMETER (DGR_TO_RAD = (QUADPI/180))

        DATA MININT/-32768/,MAXINT/32767/,FLTMIN/-10E15/,FLTZER/10E-30/,
     &     FLTMAX/10E30/
        
C       MIN LENGTH LINE TO BE RETAINED (UNLESS ITS DELETION GIVES A
C       NEW LINE LONGER THAN DELMAX)

        DELMAX = 0.1
        DELMAX = DELMAX*DELMAX/(SCALE*SCALE)


C       MAXIMUM ANGLE TO BE RETAINED
CCCC    RETANG = 178 USED TO USE THIS VALUE
        RETANG = RETANGT

        NTRY   = 0

1001    RETCOS = COS(RETANG*DGR_TO_RAD)
        RETCOSQ = RETCOS * RETCOS

           X1   = DATA(1,1)
           Y1   = DATA(2,1)
           X3   = DATA(1,2)
           Y3   = DATA(2,2)
           NK   = 1
           NTRY = NTRY + 1

           DATAT(1,1) = X1
           DATAT(2,1) = Y1

           DO 1244 N3 = 3,NDATA
            X2 = X3
            Y2 = Y3
            X3 = DATA(1,N3)
            Y3 = DATA(2,N3)

C           LENGTHS OF LINES
            DIS13 = ((X3-X1) * (X3-X1)) + ((Y3-Y1) * (Y3-Y1))
C           KEEP THE POINT IF ITS DELETION WOULD GIVE A LINE LONGER 
C           THAN THE REQUESTED MAX LENGTH
            IF (DIS13 .GT. DELMAX) GOTO 1230

            X21 = X1 - X2
            Y21 = Y1 - Y2
            X23 = X3 - X2
            Y23 = Y3 - Y2
            DIS23 = X23*X23 + Y23*Y23
C           REJECT ANY POINTS WHICH ARE CLOSER THAN RETMIN

C           FIND ANGLE BETWEEN LINES  
            DIS12   = X21*X21 + Y21*Y21
            DIS1223 = DIS12*DIS23

C           REJECT DUPLICATED POINTS OR GET DIVISION BY ZERO!
            IF (DIS1223 .LT. FLTZER) GOTO 1244  

C***            COST = (X21*X23+Y23*Y21)/SQRT(DIS1223)
C***            IF (COST .LT. RETCOS) GOTO 1244

            COSTSQ = (X21*X23+Y23*Y21)**2 / DIS1223
            IF (COSTSQ .LT. RETCOSQ) GOTO 1244

C           KEEP THE POINT
 1230       NK = NK + 1
            DATAT(1,NK) = X2
            DATAT(2,NK) = Y2
 
            X1 = X2
            Y1 = Y2           

 1244     CONTINUE                    


C         KEEP LAST POINT ON THE CONTOUR ALWAYS
          NK = NK + 1
          DATAT(1,NK) = X3
          DATAT(2,NK) = Y3
        
          IF (NK .GT. MAXPTS) THEN
C            CONTOUR IS LONGER THAN DESIRED MAXIMUM
             RETANG = RETANG - 1

C            CHECK NUMBER OF ITERATIONS
             IF (NTRY .LT. MAXTRYS) GOTO 1001
             IRTFLG = 1
             WRITE(NOUT,97) NTRY, NK 
97           FORMAT(
     &         ' THINNED:',I3,' TRIES BUT CONTOUR LENGTH STILL: ',I4)
          ELSE
             IRTFLG = 0
          ENDIF

          DO I = 1,NK
             DATA(1,I) = DATAT(1,I)
             DATA(2,I) = DATAT(2,I)
          ENDDO

          NDATA = NK

          RETURN
          END
