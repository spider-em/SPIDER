C++*********************************************************************
C
C  IMTOVOL.F   -- CREATED JUL 96 al
C                 USED AUTOMATIC ARRAYS APRIL 00 ARDEAN LEITH
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
C  IMTOVOL(LUNIM,NSAM,NROW,NSLICE,LUNVOL,FMINT,FMAXT, MAXDIM)
C
C  PURPOSE:  CONVERTS IMAGES TO VOLUME ASSUMING THAT THE IMAGE
C            VALUE IS DEPTH.  VOLUME IS BINARY CONTAINING 0 AND 1
C
C        0         2         3         4         5         6         7     
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE IMTOVOL(LUNIM,NSAM,NROW,NSLICE,LUNVOL,FMINT,FMAXT,
     &                    MAXDIM)

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC'

       REAL, DIMENSION(NSAM*NROW*NSLICE) :: VOLBUF

       COMMON /IOBUF/ BUF(NBUFSIZ)

C      ZERO THE WHOLE VOLUME
       VOLBUF = 0.0

C      FIND CONVERSION FACTOR FROM IMAGE MIN..MAX  TO DEPTH (SLICE)
       FACT = (FLOAT(NSLICE - 1) / (FMAXT - FMINT))      

C      FIND PIXELS PER SLICE
       IPPSL = NSAM * NROW

       NSLICEM1 = NSLICE - 1

       DO IROW = 1,NROW
C         READ IMAGE ROW
          CALL REDLIN(LUNIM,BUF,NSAM,IROW)
          ICON1 = (IROW - 1) * NSAM

C         CONVERT IMAGE VALUE TO DEPTH (SLICE) FOR EACH PIXEL
          DO ISAM = 1,NSAM

C            INVERT THE SLICES SO TOP IS UP
             ISLICEM1 = NSLICEM1 - (BUF(ISAM) - FMINT) * FACT
             ICON2    = ICON1 + ISAM

C            SET THE PIXEL AT THAT LOCATION TO 1.0 FOR ALL SLICES BELOW
             DO ISLICE = ISLICEM1, NSLICEM1           
                VOLBUF(ISLICE * IPPSL + ICON2) = 1.0
             ENDDO
          ENDDO
       ENDDo
    
C      SAVE THE VOLUME BUFFER  
       ILOC = 1 
       DO IREC = 1, NROW * NSLICE
          CALL WRTLIN(LUNVOL,VOLBUF(ILOC),NSAM,IREC)
          ILOC = ILOC + NSAM
       ENDDO

999    CONTINUE
       CLOSE(LUNVOL)
       CLOSE(LUNIM)

       RETURN
       END
