C ++********************************************************************
C                                                                      *
C RDCLUS.F                 ADDED FACTOR CHOICE    AUG 00  ArDean Leith  *
C                          ALLOWED 50  FACTORS    MAY 01  ArDean Leith  *
C                          REWRITE FOR NEW FORMAT NOV 03  ArDean Leith  *
C                          NFACT = NFAC BUG       JUN 09  ArDean Leith
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
c
C RDCLUS
C
C PURPOSE:  TRANSFER IMAGE COORDINATES ALONG FACTORS FROM 
C           CORRESPONDENCE ANALYSIS COORDINATES FILE (_IMC00) INTO A
C           DOCUMENT FILE. THIS DOCUMENT FILE CAN BE USED IN WEB 
C           'CORR MAP' OPERATION.                                                                  *
C***********************************************************************

	SUBROUTINE RDCLUS

	INCLUDE 'CMBLOCK.INC' 
	INCLUDE 'CMLIMIT.INC' 

C       MAXNAM, INUMBR & NIMAX IS FROM: CMLIMIT.INC
        CHARACTER(LEN=MAXNAM) :: FILPRE,FILNMC,COMMENT
	REAL,ALLOCATABLE      :: FACLIST(:),DLIST(:)
        CHARACTER(LEN=1)      :: NULL

	DATA  LUNI,NDOC/70,71/

        NULL = CHAR(0)

        NILMAX = NIMAX                    ! FROM CMLIMIT
        ALLOCATE(FACLIST(NIMAX),
     &           DLIST(NIMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'RDCLUS; FACLIST....',2*NIMAX)
           RETURN
        ENDIF


        CALL FILERD(FILPRE,NLET,NULL,
     &              'CORAN/PCA FILE PREFIX (e.g. CORAN_)~',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       OPEN IMAGE COORDINATE FILE
        FILNMC = FILPRE(1:NLET) // '_IMC'// NULL
        CALL OPAUXFILE(.FALSE.,FILNMC,DATEXC,LUNI,0,
     &                       'O', ' ',.TRUE.,IRTFLG)
        IF (IRTFLG .NE.0) RETURN

C       READ HEADER 
C       NITEM -- NUMBER OF IMAGES
C       NFAC  -- NUMBER OF FACTORS USED IN THE ANALYSIS
C       NSAM, NROW -- NUMBER OF SAMPLES AND ROWS IN THE UNMASKED IMAGE

        READ(LUNI,*) NITEM, NFAC, NSAM, NROW, NDUM, KIND_PCA

        WRITE(NOUT,*) ' NUMBER OF FACTORS AVAILABLE: ',NFAC

        NMAX = NIMAX
        NFACT = NFAC
        CALL RDPRAI(INUMBR,NMAX,NFACT,1,NFAC,'FACTOR NUMBERS',
     &              NULL,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

	DO  I=1,NITEM
C          WRITE(LUNI,90) (BLW(K),K=1,NFAC),WEIGHTI(J),CO(J),FIM,ACT

           READ(LUNI,*) (DLIST(K),K=1,NFAC),RWGT,DIO,FACLIST(1),ACTIV
C          KEY        = IMAGE NUMBER

           DO J = 1,NFACT
              FACLIST(J+1) = DLIST(INUMBR(J))
           ENDDO
           CALL SAVD(NDOC,FACLIST,NFACT+1,IRTFLG)
	ENDDO

        COMMENT = ' '
        WRITE(COMMENT,95)(INUMBR(I),I=1,NFACT)
95      FORMAT(' FACTORS: ',9(I2,', '),' ...')

        CALL LUNDOCPUTCOM(NDOC,COMMENT(1:77),IRTFLG)

	CALL SAVDC
	CLOSE(NDOC)
        CLOSE(LUNI)
        IF (ALLOCATED(FACLIST)) DEALLOCATE(FACLIST)
        IF (ALLOCATED(DLIST))   DEALLOCATE(DLIST)

        RETURN
	END
