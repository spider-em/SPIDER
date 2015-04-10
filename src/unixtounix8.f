
C++*********************************************************************
C
C UNIXTOUNIX8.F  WRITTEN                             JUL 93 ARDEAN LEITH
C                REWRITTEN                           FEB 99 ARDEAN LEITH
C                NATIVE BYTE_ORDER                   JUL 09 ARDEAN LEITH
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
C    UNIXTOUNIX8(LUNO,LUNN,NMIN,NMAX,FMINT,FMAXT,
C               NSAM,NROW,NSLICE,NORMAL,IPAD,IRTFLG)
C
C    PURPOSE:       CONVERT A UNIX 32 BIT IMAGE FILE TO 8 BIT RAW
C
C    PARAMETERS:
C        LUNO       INPUT IO UNIT                                 (SENT)
C        LUNN       OUTPUT IO UNIT                                (SENT)
C        NMIN,NMAX  NORMALIZATION INTERVAL (USUALLY 0..255)       (SENT)
C        FMIN,FMAX  IMAGE RANGE                                   (SENT)
C        NORMAL     LOGICAL VARIABLE FOR NORMALIZATION WANTED     (SENT)
C        IPAD       PADDING AMT                                   (SENT)
C        IRTFLG     ERROR RETURN FLAG. (0 IS NORMAL)              (RET.)
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

	SUBROUTINE UNIXTOUNIX8(LUNO,LUNN,NMIN,NMAX,FMINT,FMAXT,
     &             NSAM,NROW,NSLICE,NORMAL,IPAD,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'
        COMMON /IOBUF/ BUFO(NBUFSIZ)

        INTEGER   * 1  :: LBUF
        INTEGER * 2    :: I2VAL
	LOGICAL        :: NORMAL

        IRTFLG = 1

        IF (NORMAL) THEN
           FN    = (NMAX - NMIN) / (FMAXT - FMINT)
           FNCON = NMIN - FN * FMINT

        ELSEIF (FMAXT .GT. 255) THEN
C          CAN NOT USE THIS FOR 8 BIT OUTPUT
           WRITE(NOUT,*) ' FMAXT:',FMAXT,
     &              'TOO LARGE FOR 8 BIT UN-NORMALIZED OUTPUT'
           CALL ERRT(100,'UNIXTOUNIX8',IDUM)
           RETURN

        ELSEIF (FMINT .LT. 0) THEN
C          CAN NOT USE THIS FOR 8 BIT OUTPUT
           WRITE(NOUT,*) ' FMINT:',FMINT,
     &              'TOO SMALL FOR 8 BIT UN-NORMALIZED OUTPUT'
           CALL ERRT(100,'UNIXTOUNIX8',IDUM)
           RETURN
        ELSE
           FN    = 1.0
           FNCON = 0.0
        ENDIF
  
        IRECOUT = 0
        DO  I = 1,NROW*NSLICE
C          READ EACH RECORD OF 32 BIT UNIX SPIDER INPUT FILE   
           CALL REDLIN(LUNO,BUFO,NSAM,I)

C          CONVERT FLOATING POINT NUMBERS TO -128...127 LOGICAL *1
           DO J=1,NSAM
	     I2VAL         = BUFO(J) * FN + FNCON
             IRECOUT       = IRECOUT + 1
	     LBUF          = I2VAL
             WRITE(LUNN,90,REC=IRECOUT,IOSTAT=IRTFLGT) LBUF
90           FORMAT(A1)
             IF (IRTFLGT .NE. 0) THEN
                WRITE(NOUT,*) '*** ERROR WRITING PIXEL:',IRECOUT
                CALL ERRT(101,'UNIXTOUNIX8',IDUM)
                RETURN
             ENDIF

           ENDDO

           IF (IPAD .GT. 0) THEN
C            ADD PADDING AFTER EACH LINE
             DO L = 1,IPAD
                LBUF    = 0
                IRECOUT = IRECOUT + 1
                WRITE(LUNN,90,REC=IRECOUT,IOSTAT=IRTFLGT) LBUF
                WRITE(LUNN,REC=IRECOUT,IOSTAT=IRTFLGT) LBUF
                IF (IRTFLGT .NE. 0) THEN
                   WRITE(NOUT,*) '*** ERROR WRITING PIXEL:',IRECOUT
                   CALL ERRT(101,'UNIXTOUNIX8',IDUM)
                   RETURN
                ENDIF
             ENDDO
           ENDIF   
        ENDDO 

        WRITE(NOUT,93) IRECOUT
93      FORMAT(' RAW IMAGE SIZE:', I10,' BYTES'/)

        IRTFLG = 0

        RETURN
	END
 

c***************************
c           it1 = lval1
c           it2 = lval2
c           write(6,*) 'i2val,j,bufo(j):',i2val,j,bufo(j),it1,it2
c           iend = iend + 1
c           if (iend .gt. 20) stop
c******************************
 
