
C++*********************************************************************
C
C  COPYTOJPG.F  -- NEW                             APR 13 ArDean Leith                    ARDEAN LEITH
C                  ECHO                            APR 16 ArDean Leith
C **********************************************************************
C=* AUTHOR: A. LEITH                                                   *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
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
C=* merchantability or fitness for a particular purpose.  See the      *
C=* GNU General Public License (www.gnu.org/licenses) for details.     *
C=*                                                                    *
C **********************************************************************
C
C   COPYTOJPG(LUNO,FILOLD,LUNT,NX,NY,NZ)
C
C   PURPOSE:      CONVERT SPIDER IMAGE FILE TO JPEG FORMAT USING
C                 IMAGEMAGICK
C
C   CALLED BY:    COPY1
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

        INTEGER               :: NRECS
        INTEGER               :: IRTFLG

        REAL                  :: BUF(NX)
        CHARACTER *240        :: COMLIN
        CHARACTER *16         :: STRNX, STRNY
        INTEGER               :: NC2, NC3
        INTEGER               :: NDUM,LENREC,IERR,NLETN
        INTEGER               :: IROW
        REAL                  :: FLOW,FHI,RANGENEW
        REAL                  :: RANGEOLD,CON,CON2
        CHARACTER * 1         :: NULL     = CHAR(0)
        LOGICAL               :: ASKNAM   = .FALSE.
        LOGICAL               :: CALLERRT = .TRUE.

        CHARACTER(LEN=17)     :: FILTMP = 'JUNK_FOR_JPG.gray' 
 
        INTEGER               :: lnblnkn

        IF (NZ .NE. 1) THEN
           CALL ERRT(101,'DOES NOT WORK ON VOLUMES',NDUM)
           RETURN
        ENDIF

C       OPEN TEMP FILE FOR NORMALIZED DATA
C       I DO NOT SEE ANY WAY THAT IMAGEMAGICK CAN HANDLE FLOATING DATA
C       THAT HAS NEGATIVE VALUES WITH ANY OPTION I HAVE TRIED
       
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

        DO  IROW = 1,NY
           CALL REDLIN(LUNO,BUF,NX,IROW)

           BUF = CON2 + CON * BUF    ! ARRAY OP

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


C **********************************************************************
C
C  DISP
C
C  PURPOSE:  COPIES SPIDER IMAGE TO JPEG USING IMAGEMAGICK THEN USES  
C            SYSTEM COMMAND TO DISPLAY IMAGE. 
C
C--*******************************************************************

        SUBROUTINE DISP()

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FILOLD,FILNEW
        CHARACTER(LEN=MAXNAM) :: OPTIONS
        CHARACTER(LEN=160)    :: COMLIN
 
        CHARACTER(LEN=1)      :: NULL = CHAR(0)
        LOGICAL               :: VERBOSET,WANTOUT
        INTEGER               :: NX,NY,NZ,MAXIM,ITYPE,IRTFLG,NLET,NLETC
        INTEGER               :: ICOMM,MYPID,MPIERR, lnblnkn,NLETO

        INTEGER, PARAMETER    :: LUN1   = 14 
        INTEGER, PARAMETER    :: LUN2   = 15 
        INTEGER, PARAMETER    :: IDELAY = 3
          
        CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C       OPEN INPUT FILE, WHOLE STACK NOT ALLOWED
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILOLD,LUN1,'O',ITYPE,
     &               NX,NY,NZ,MAXIM,'SPIDER INPUT',
     &               .FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (IMAMI .NE. 1) 
     &      CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV)

C       JPEG OUTPUT FILE NAME
        NLET   = lnblnkn(FILOLD)
        FILNEW = FILOLD(:NLET) // '.jpg'
        NLET   = NLET + 4

C       CONVERT SPIDER IMAGE FILE INTO JPG FILE
        VERBOSET = .FALSE.
        CALL COPYTOJPG(LUN1,LUN2,FILNEW,NX,NY,NZ, VERBOSET,IDELAY)

        IRTFLG = -999   ! KEEP LOWERCASE
        CALL RDPRMC(OPTIONS,NLETO,.TRUE.,
     &          'IMAGEMAGICK DISPLAY OPTIONS (or <CR>) ',NULL,IRTFLG)

        WRITE(COMLIN,90) OPTIONS(1:NLETO),FILNEW(1:NLET)
90      FORMAT( ' display ', A,' ',A, ' &' )

C       DO NOT ECHO COMLIN
        CALL CSVMS(COMLIN,.FALSE.,IRTFLG)

        END
