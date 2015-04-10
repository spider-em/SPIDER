
C++*********************************************************************
C
C  SURFTOVOL.FOR -- CREATED JUL 96 al
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
C  SURFTOVOL(LUNVOL,NSAM,NROW,NSLICE,NX,NY,NZ,TINY)
C
C  PURPOSE:  CONVERTS PROGRAMMED SURFACE FUNCTION TO VOLUME. 
C            VOLUME IS CONTINUOUS FROM 0.0 TO 1.0, EXCEPT FOR 
C            SURFACE 2D WHICH IS CONTINUOUS FROM 0.4...1.0.   
C
C  PARAMETERS:   LUNVOL             LUN FOR OUTPUT              (SENT)
C                NSAM,NROW,NSLICE   VOLUME SIZE                 (SENT)
C                FNXYZ              REPEATS IN X, Y, & Z        (SENT)
C                IORDER             ORDER FOR  X, Y, & Z        (SENT)
C                ANS                SURFACE TYPE (3 AVAIL.)     (SENT)
C                TINY               ANYTHING < THIS IS SURFACE  (SENT)
C
C        0         2         3         4         5         6         7     
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE SURFTOVOL(LUNVOL,ISIZE,FNXYZ,IORDER,ANS,TINY)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       COMMON       BUF(NBUFSIZ)

       CHARACTER    *(*) ANS
       INTEGER      ISIZE(3),IORDER(3)
       REAL         FNXYZ(3)
       LOGICAL      GEE,TWOD

C      X, Y, & Z AXIS ARE INTERCHANGEABLE
       NSAM   = ISIZE(IORDER(1))
       NROW   = ISIZE(IORDER(2))
       NSLICE = ISIZE(IORDER(3))

       FNX    = FNXYZ(IORDER(1))
       FNY    = FNXYZ(IORDER(2))
       FNZ    = FNXYZ(IORDER(3))

C      TYPE D AND TYPE G SURFACES AVAILABLE
       GEE   = ANS(1:1) .EQ. 'G'  .OR. ANS(1:1) .EQ. 'g'
       TWOD  = ANS(1:2) .EQ. 'D2' .OR. ANS(1:2) .EQ. 'd2' 

       TWOPI = 2 * 3.14159
       IREC  = 0
      
       DO ISLICE = 1,NSLICE
          Z = FLOAT(ISLICE) / (NSLICE / FNZ) 

          DO IROW = 1,NROW

             Y     = FLOAT(IROW) / (NROW / FNY) 
             FXCON = 1.0 / (NROW  / FNX)

             IF (GEE) THEN
C               G SURFACE SCALED SO THAT DESIRED SURFACE IS NEAR 1.0 
                DO ISAM = 1,NSAM
                   X = FLOAT(ISAM) * FXCON 

                   VAL = SIN(TWOPI * X) * COS(TWOPI * Y) +
     &                   SIN(TWOPI * Y) * COS(TWOPI * Z) +
     &                   COS(TWOPI * X) * SIN(TWOPI * Z) 
      
                   BUF(ISAM) = 1.0 - ABS(VAL)
                ENDDO
             ELSEIF (.NOT. TWOD) THEN
C               D SURFACE SCALED SO THAT DESIRED SURFACE IS NEAR 1.0 
                DO ISAM = 1,NSAM
                   X = FLOAT(ISAM) * FXCON 

                   VAL = COS(TWOPI * X - TWOPI * Y) * COS(TWOPI * Z) +
     &                   SIN(TWOPI * X + TWOPI * Y) * SIN(TWOPI * Z)
      
                   BUF(ISAM) = 1.0 - ABS(VAL)
                ENDDO

             ELSEIF (TWOD) THEN
C               LEAVE RAW DATA IN VOLUME
                DO ISAM = 1,NSAM
                   X = FLOAT(ISAM) * FXCON 

                   BUF(ISAM) =
     &                   COS(TWOPI * X - TWOPI * Y) * COS(TWOPI * Z) +
     &                   SIN(TWOPI * X + TWOPI * Y) * SIN(TWOPI * Z)
      
                ENDDO
             ENDIF
             IREC = IREC + 1
             CALL WRTLIN(LUNVOL,BUF,NSAM,IREC)
          ENDDO
       ENDDo
    
       CONTINUE
       CLOSE(LUNVOL)

       RETURN
       END
