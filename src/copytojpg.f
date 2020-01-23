
C++*********************************************************************
C
C  COPYTOJPG.F  -- NEW                             APR 13 ArDean Leith
C                  ECHO                            APR 16 ArDean Leith
C                  WORKS ON VOLUME SLICES NOW      JAN 18 ArDean Leith
C                  WORKS ON MRC IMAGES NOW         DEC 19 ArDean Leith
C                  EXTRACTED DISP()                JAN 20 ArDean Leith
C **********************************************************************
C=* AUTHOR: A. LEITH                                                   *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2020  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
C=*                                                                    *
C=* SPIDER is free software; you can redistribute it and/or            *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* SPIDER is distributed in the hope that it will be useful,          *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* merchantability or fitness for a particular purpose.  See the      *
C=* GNU General Public License (www.gnu.org/licenses) for details.     *
C=*                                                                    *
C **********************************************************************
C
C   COPYTOJPG(LUNO,FILOLD,LUNT,NX,NY,NZ,VERBOSET,IDELAY))
C
C   PURPOSE:      CONVERT SPIDER IMAGE FILE TO JPEG FORMAT USING
C                 IMAGEMAGICK,  CAN ALSO DISPLAY IMAGE
C
C   CALLED BY:    COPY1, UTIL1
C
C--*********************************************************************

        SUBROUTINE COPYTOJPG(LUNO,LUNT,FILNEW,NX,NY,NZ,VERBOSET,IDELAY)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER               :: LUNO,LUNT,NX,NY,NZ
        LOGICAL               :: VERBOSET
        CHARACTER(LEN=MAXNAM) :: FILNEW
        INTEGER               :: IDELAY

        INTEGER               :: NRECS,ISLICE,NOT_USED,ISKIP
        INTEGER               :: IRTFLG

        REAL                  :: BUF(NX)
        CHARACTER *240        :: COMLIN
        CHARACTER *16         :: STRNX, STRNY
        INTEGER               :: NC2, NC3
        INTEGER               :: NDUM,LENREC,IERR,NLETN
        INTEGER               :: IROW,IDUM
        REAL                  :: FLOW,FHI,RANGENEW
        REAL                  :: RANGEOLD,CON,CON2
        CHARACTER * 1         :: NULL     = CHAR(0)
        LOGICAL               :: ASKNAM   = .FALSE.
        LOGICAL               :: CALLERRT = .TRUE.
        LOGICAL               :: erri2

        CHARACTER(LEN=17)     :: FILTMP = 'JUNK_FOR_JPG.gray' 
 
        INTEGER               :: lnblnkn

        ISLICE = 1
        IF (NZ > 1) THEN
C          THIS IS A VOLUME
           CALL RDPRI1S(ISLICE,NOT_USED,'SLICE',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
           IF (ERRI2(ISLICE,IDUM, 1, 1,NZ, IDUM,IDUM)) RETURN

        ENDIF

C       OPEN TEMP FILE FOR NORMALIZED DATA
C       I DO NOT SEE ANY WAY THAT IMAGEMAGICK CAN HANDLE FLOATING DATA
C       THAT HAS NEGATIVE VALUES WITH ANY OPTION I HAVE TRIED!!
       
        LENREC = NX * 4
        CALL OPAUXFILE(ASKNAM,FILTMP,NULL,LUNT,LENREC,'U',
     &                 ' ',CALLERRT,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (IMAMI == 0) THEN
C          GET INPUT IMAGE STATISTICS FIRST
           CALL NORM3(LUNO,NX,NY,NZ,FMAX,FMIN,AV)
        ENDIF
        IF (FMAX == FMIN) THEN
           CALL ERRT(101,'BLANK FILE SKIPPED',NDUM)
           RETURN
        ENDIF
  
        FLOW      = 0.0   ! RANGE NEEDED BY IMAGEMAGICK
        FHI       = 1.0

        RANGENEW  = FHI  - FLOW
        RANGEOLD  = FMAX - FMIN
        CON       = RANGENEW / RANGEOLD
        CON2      = FLOW - CON * FMIN

        ISKIP    = NY * (ISLICE - 1)
        DO  IROW = 1,NY
           CALL REDLIN(LUNO,BUF,NX,IROW + ISKIP)

           BUF = CON2 + CON * BUF    ! ARRAY OP

C          WRITE TO TEMP FILE AS POSITIVE FLOATING POINT
           WRITE(LUNT,REC=IROW,IOSTAT=IERR) BUF
        ENDDO
        CLOSE(LUNO)
        CLOSE(LUNT)

        IF (FILNEW == NULL) THEN
C          GET NAME FOR JPEG FILE
        
           CALL FILERD(FILNEW,NLETN,'jpg','JPEG OUTPUT~9',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
        ELSE
           NLETN = lnblnkn(FILNEW)
        ENDIF
         
        CALL INTTOCHAR(NX, STRNX, NC2,1)
        CALL INTTOCHAR(NY, STRNY, NC3,1)

        IF (VERBOSET) THEN
           WRITE(COMLIN,8005) STRNX(:NC2),STRNY(:NC3),
     &                        FILTMP,FILNEW(1:NLETN)
8005       FORMAT( 
     &     ' convert -verbose -depth 32 -size ',A,'x',A, 
     &     ' -define quantum:format=floating-point ',
     &     '-define quantum:scale=65536.0 -endian msb ',
     &     A,' ',A )
        ELSE
           WRITE(COMLIN,8006) STRNX(:NC2),STRNY(:NC3),
     &                        FILTMP,FILNEW(1:NLETN)
8006       FORMAT( 
     &     ' convert -depth 32 -size ',A,'x',A, 
     &     ' -define quantum:format=floating-point ',
     &     '-define quantum:scale=65536.0 -endian msb ',
     &     A,' ',A, ' >& /dev/null' )

        ENDIF

        !write(6,*) COMLIN
 
C       ECHO COMLIN
        CALL CSVMS(COMLIN,.TRUE.,IRTFLG)

C       A DELAY HERE. (KLUDGY)
        CALL sleep(IDELAY)

        CLOSE(LUNT,STATUS='DELETE',IOSTAT=IRTFLG)

        END


