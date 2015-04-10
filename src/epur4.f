
C ++********************************************************************
C                                                                      *
C EPUR4     REMOVED FROM HPLAN.FOR FILE        APRIL 89 al                    *
C           ID(IDIM,NID)                       12/7/93 ML                                  *             
C           UPDATED SOME                       AUG 1999 ARDEAN LEITH                           *
C           SIMPLIFIED LOGIC FOR OUTLIERS      DEC 2005 ARDEAN LEITH                                                           *
C           KLIC REMOVED (BUGGY)               JAN 2014 ARDEAN LEITH  
C                                                         *
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C   EPUR4(IDIM,NPTS, X,Y,ID, MOD, PEX, KP,KLIC,IRTFLG,NDAT) 
C 
C   POINTS WHICH FALL BEYOND PEX STANDARD DEVIATIONS COLLAPSED ONTO         
C   THE FRAME OF THE GRAPH.
C                                                     
C   OUTPUT:                 
C   IF THERE ARE MORE THAN 2048 POINTS ON THE FRAME THE GRAPH IS ABORTED         
C   IF MOD=1 LABELS HAVE A1 FORMAT - IF MOD=4 LABELS HAVE A7 FORMAT 
C            
C   WARNING: X(*) AND Y(*) ARE DESTROYED IF KP .NE. 0                           
C
C   CALLED BY: HISMAP,  HISMAP4
C
C **********************************************************************

        SUBROUTINE EPUR4(IDIM,NPTS,X,Y,ID,MOD,PEX,KP,IDUMC,IRTFLG,NDAT)

        CHARACTER(LEN=7) :: ID(IDIM)
        REAL             :: X(IDIM), Y(IDIM)

        IRTFLG  = 1
        IF (PEX <= 0.0) RETURN

        NPTSP1  = NPTS + 1
        IF (NPTSP1 > IDIM) RETURN
                               
        X(NPTSP1) = 0.0
        Y(NPTSP1) = 0.0

        NPTSP1    = NPTS
        SX        = 0.0
        SY        = 0.0
        DO I = 1,NPTS                                                      
           SX = SX + X(I) * X(I)
           SY = SY + Y(I) * Y(I)
        END DO

        SX    = SQRT(SX/FLOAT(NPTS))
        SY    = SQRT(SY/FLOAT(NPTS))
        PX    = PEX * SX 
        PY    = PEX * SY
        KP    = 0 

        DO  I = 1, NPTS                                                      
          IF (ABS(X(I)) > PX .OR. ABS(Y(I)) > PY) THEN
C           POINT IS BEYOND SD BOUNDARY REQUESTED 
            KP = KP + 1 
            IF (MOD .NE. 1 .AND. KP == 1) THEN
C              HEADER FOR OUTLIER LISTING
               WRITE(NDAT,90) PEX 
90            FORMAT (//,'  WARNING, FOLLOWING POINTS BEYOND: ',F5.1,
     &          ' SD FROM ORIGIN PLACED ON EDGES OF MAP!',/)

               WRITE(NDAT,91)
91             FORMAT('  ** POINT ********* LOCATION ',12('*'))
            ENDIF

            IF (KP > 2048) THEN     !OVERFLOW TRAP
               CALL ERRT(101,'MORE THAN 2048 POINTS ON FRAME',IDUM)
               RETURN
            ENDIF
              
            !KLIC(KP) = I  ! not used anymore al

            IF (MOD .NE. 1) WRITE(NDAT,92) ID(I),X(I),Y(I)
92          FORMAT ('     ',A7,'  (',F12.5,',',F12.5,')')

C            LOCATIONS SET TO: FRAME
             IF (ABS(X(I)) > PX)  X(I) = SIGN(PX, X(I))
             IF (ABS(Y(I)) > PY)  Y(I) = SIGN(PY, Y(I))
          ENDIF
        ENDDO

        IF (KP .NE. 0) WRITE(NDAT,93)
93      FORMAT('  ',40('*'),//)

        IRTFLG  = 0

        END 
                                                                    
