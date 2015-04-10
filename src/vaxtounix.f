
#ifdef SP_SUN4 
 
C   THIS ROUTINE NOT AVAILABLE ON SUN
 
       SUBROUTINE vaxtounix(FILOLD,FILNEW,LUNO,LUNN,IRTFLG)
 
       COMMON /UNITS/LUNC,NIN,NOUT
 
       WRITE(NOUT,*) 'DUMMY CALL: vaxtounix'
       RETURN
       END
#else  

C++*********************************************************************
C
C  VAXTOUNIX.FOR  -- ADAPTED FROM GETVAX
C                    CONVERTED TO RUN ON UNIX           JULY 93 al
C                    TOTREC REMOVED                     DEC. 10 al
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
C    VAXTOUNIX(FILOLD,FILNEW,LUNO,LUNN,IRTFLG)
C
C    PURPOSE:       TO CONVERT A VAX SPIDER IMAGE FILE TO UNIX FORMAT
C                   RUNNING ON A UNIX MACHINE
C
C    PARAMETERS:
C        FILOLD     VAX FORMAT FILE NAME                       (SENT)
C        LUNO       LOGICAL UNIT NUMBER  FILOLD.               (SENT)
C        FILNEW     UNIX  FORMAT FILE NAME                     (SENT)
C        LUNN       LOGICAL UNIT NUMBER  FILNEW                (SENT)
C        IRTFLG     ERROR RETURN FLAG. (0 IS NORMAL)           (RET.)
C
C    VARIABLES:   IFORM  = FILE TYPE SPECIFIER. 
C                        = +3    3-D IMAGE
C	                 = +1    2-D IMAGE
C                        = -1    2-D FOURIER TRANSFORM
C                        = -3    3-D FOURIER TRANSFORM
C                        = +8    8 BIT BLACK AND WHITE IMAGE
C                        = 11    8 BIT COLOR IMAGE
C
C--********************************************************************

	SUBROUTINE VAXTOUNIX(FILOLD,FILNEW,LUNO,LUNN,IRTFLG)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOERR/  IERR
        COMMON /IOBUF/ BUFO(NBUFSIZ)

        CHARACTER *(*) FILOLD,FILNEW
        CHARACTER *1   NULL

        NULL   = CHAR(0)

        IRTFLG = 1

C       RETRIEVE ALL LABEL INFO (LABEL CONTAINS >=256 * 4 BYTES)
        CALL OPAUXFILE(.TRUE.,FILOLD,DATEXC,LUNO,1024,'O',
     &                  'VAX/VMS INPUT FILE',.TRUE.,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       READ FIRST 1024 BYTES OF VAX INPUT FILE (FIRST RECORD)
        CALL REDLIN8(LUNO,BUFO,1024,1,IRTFLGT)
        IF (IRTFLGT .NE. 0) THEN
           CALL ERRT(101,'*** CAN NOT READ VAX HEADER',NE)
           RETURN
        ENDIF

C       CONVERT HEADER TO UNIX FLOATING POINT NUMBERS           
        DO I=1,211
           CALL VAX32U(BUFO(I))
        ENDDO

        NSLICE = BUFO(1)
        IF (NSLICE .EQ. 0) THEN
C          FIX UP OLD LABELS NSLICE
           NSLICE  = 1
           BUFO(1) = NSLICE

        ELSEIF (NSLICE .GT. 0) THEN
C          OLD SPIDER SHORT LABEL FILE
           WRITE(NOUT,*) '*** SHORT LABEL CONVERSION NOT IMPLEMENTED'
           CALL ERRT(2,'VAXTOUNIX',NE)
           RETURN
        ELSE
            NSLICE = -NSLICE
        ENDIF

        NROW     = BUFO(2)

        IFORM    = BUFO(5)
        IF (BUFO(5) .LT. 0.0) IFORM = BUFO(5) - 0.5
        IF (IFORM .NE. 1 .AND. IFORM .NE. 3) THEN
          WRITE(NOUT,*)
     &         'WARNING: NOT TESTED FOR THIS FORMAT: ',IFORM
        ENDIF

C       PIXELS IN EACH RECORD IS IN BUFO(12) = NSAM
        NSAM     = BUFO(12)

        IF (NSAM .LE. 0 .OR. NSAM .GT. 8000) THEN
           WRITE(NOUT,*) 'NSAM IN HEADER LOCATION 12:', NSAM
           CALL RDPRI1S(NSAM,NOT_USED,'NSAM',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ENDIF

        LENBYT   = NSAM * 4
        BUFO(23) = LENBYT

        LABREC   = 1024 / LENBYT
        IF (MOD(1024,LENBYT) .NE. 0) LABREC = LABREC + 1
        BUFO(13) = LABREC

        LABBYT   = LABREC * LENBYT
        BUFO(22) = LABBYT

C       CLOSE THE VAX FILE AND REOPEN WITH CORRECT RECORD LENGTH
        CLOSE(LUNO)

        CALL OPAUXFILE(.FALSE.,FILOLD,NULL,LUNO,NSAM*4,'O',
     &                  ' ',.TRUE.,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       PLACE HEADER IN NEW UNIX FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNEW,LUNN,'N',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'OUTPUT',.FALSE.,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       FIND NUMBER OF IMAGE RECORDS EXPECTED IN FILE (IEND)
c       CALL TOTREC(NSAM,NROW,NSLICE,IFORM,IEND) removed dec 2010

        IF (IFORM .EQ. -1) THEN
C          2D FOURIER FILE
           IENDI = (NSAM + 2) * NROW
           IEND  = IENDI / NSAM
           IF (MOD(IENDI,NSAM) .NE. 0) IEND = IEND + 1

        ELSE
C          OTHER FILE TYPE, ASSUME JUST IMAGE RECORDS
           IEND = NROW * NSLICE
        ENDIF

        IERR = 0
        DO I = 1,IEND

C           READ EACH RECORD OF VAX FILE
            CALL REDLIN(LUNO,BUFO,NSAM,I+LABREC)

            IF (IERR .NE. 0) THEN
              CALL ERRT(12,'VAXTOUNIX',NE)
              RETURN
            ENDIF

C           CONVERT FLOATING POINT NUMBERS TO UNIX FORMAT
            DO J=1,NSAM
              CALL VAX32U(BUFO(J))
            ENDDO

C           WRITE RECORD TO UNIX FILE
            CALL WRTLIN(LUNN,BUFO,NSAM,I)

        ENDDO

C       REOPEN TO ECHO DATE AND TIME
        CLOSE(LUNN)
        MAXIM = 0
        WRITE(NOUT,*) ' '
        CALL OPFILEC(0,.FALSE.,FILNEW,LUNN,'O',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,' ',.TRUE.,IRTFLGT)
        IF (IRTFLGT .NE. 0) RETURN

C       SET FLAG FOR NORMAL RETURN	
        IRTFLG = 0

        RETURN
        END
 
#endif
 
