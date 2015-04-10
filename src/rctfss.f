C++*********************************************************************
C
C RCTFSS.F                 ADDED 2D OPERATION       NOV 00    HAI  XIAO
C                          INPIC BUG                APR 02 ARDEAN LEITH
C                          OPFILEC                  FEB 03 ARDEAN LEITH
C                          ILIST FIXED              MAY 09 ARDEAN LEITH
C                          PUT MULFC3 HERE          NOV 10 ARDEAN LEITH
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2012  Health Research Inc.,                         *
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
C  RCTFSS(LUN1,LUN2)
C
C  PARAMETERS   
C
C  PURPOSE: 2D and 3D CTF CORRECTION USING WIENER FILTERING APPROACH
C
C           F(K)=SUM(H I(K)FI(K))/(SUM(HI(K)^2)+1/SNR), 
C           WHERE: HI(K) IS CTF FUNCTION FOR THE ITH IMAGE;
C           FI(K) IS THE ITH IMAGE; SNR IS SIGNAL-TO-NOISE RATIO.
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE RCTFSS(LUN1,LUN2)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        CHARACTER(LEN=MAXNAM) :: FINPAT,FILNAM,FINPIC
        INTEGER, ALLOCATABLE  :: ILIST(:)
        COMPLEX, ALLOCATABLE  :: VOLIN(:,:,:) 
        COMPLEX, ALLOCATABLE  :: VOLSUM(:,:,:)
        COMPLEX, ALLOCATABLE  :: CTF2(:,:,:)
        CHARACTER(LEN=1)      :: NULL
        LOGICAL               :: FIRST

        CALL SET_MPI(ICOMM,MYPID,MPIERR)  ! SET MYPID

        NULL = CHAR(0)

        NILMAX = NIMAX
        ALLOCATE(ILIST(NILMAX),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL  ERRT(46,'TF CTS, ILIST',NILMAX)
           RETURN
        ENDIF

C       ASK FOR DATA FILE
        CALL FILELIST(.TRUE.,LUN2,FINPAT,NLET,ILIST,NILMAX,NIMA,
     &      'TEMPLATE FOR IMAGE/VOLUME FILES',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        CALL FILERD(FILNAM,NLEP,NULL,
     &       'TEMPLATE FOR IMAGE/VOLUME CTF',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        SNR = 100
        CALL RDPRM1S(SNR,NOT_USED,'SNR',IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999
        IF (SNR <= 0.0)  THEN
           CALL ERRT(101,'SNR MUST BE > 0',IER)
           GOTO 9999
        ENDIF

        S = 1.0 / SNR

        DO  K=1,NIMA

C          READ CURRENT IMAGE/VOLUME FILE
           CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,LUN1,'O',IFORM,
     &               NSAM1,NROW1,NSLICE1,
     &               MAXIM,'DUMMY',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0)  GOTO 9999

C          PROCESS THE VOLUME

           INUMDEM = 0  ! INUMDEM RECORDS 3D OR 2D INUMDEM: 3=3D 1=2D

           IF (IFORM==-11)  THEN
              LRCL    = NSAM1
              LS      = NSAM1
              NSAM1   = NSAM1-1
              INUMDEM = 1
           ELSEIF (IFORM==-12)  THEN
              LRCL    = NSAM1
              LS      = NSAM1
              NSAM1   = NSAM1-2
              INUMDEM = 1
           ELSEIF (IFORM == 1) THEN
              LS      =NSAM1+2-MOD(NSAM1,2)
              LRCL    =NSAM1
              INUMDEM = 1
           ELSEIF (IFORM == -21)  THEN
              LRCL    = NSAM1
              LS      = NSAM1
              NSAM1   = NSAM1-1
              INUMDEM = 3
           ELSEIF (IFORM == -22)  THEN
              LRCL    = NSAM1
              LS      = NSAM1
              NSAM1   = NSAM1-2
              INUMDEM = 3
           ELSEIF (IFORM == 3) THEN
              LS      = NSAM1+2-MOD(NSAM1,2)
              LRCL    = NSAM1
              INUMDEM = 3
           ELSE
              CALL ERRT(2,'TF CTS',NE)
              GOTO 9999
           ENDIF

C          GET THE DIMENSION INFORMATION FOR THE FIRST TIME
           IF (K == 1) THEN
              NSAM   = NSAM1
              NROW   = NROW1
              NSLICE = NSLICE1

C             VOLUMES NEEDED

              MWANT = (LS/2) * NROW * NSLICE
              ALLOCATE(VOLIN( LS/2,NROW,NSLICE),
     &                 VOLSUM(LS/2,NROW,NSLICE),
     &                 CTF2(  LS/2,NROW,NSLICE),  STAT=IRTFLG)
              IF (IRTFLG .NE. 0) THEN
                  CALL ERRT(46,'RCTFSS, VOLIN...',MWANT*3)
                  RETURN
              ENDIF
   
C             SET OUTPUT VOLUME AND CTF^2 TO ZERO
              VOLSUM = 0.0
              CTF2   = 0.0
           ELSE
C
C             VERIFY THE DIMENSIONS
              IF (NSAM.NE.NSAM1.OR.NROW.NE.NROW1.OR.NSLICE.NE.NSLICE1)
     &           THEN
                 CALL ERRT(101,'INCONSISTENT DIMENSIONS',NE)
                 GOTO 9999
              ENDIF
           ENDIF

C          READ INPUT FILE
           DO L=1,NSLICE
              DO J=1,NROW
                 CALL REDLIN(LUN1,VOLIN(1,J,L),LRCL,J+(L-1)*NROW)
              ENDDO
           ENDDO
           IF (MYPID <= 0) CLOSE(LUN1)

           IF ((IFORM == 1) .OR. (IFORM == 3)) THEN
C             FOURIER TRANSFORM
              INV = +1
              CALL FMRS_3(VOLIN,NSAM,NROW,NSLICE,INV)

              IF (INV == 0)  THEN
                 CALL ERRT(101,'FFT ERROR',NE)
                 GOTO 9999
              ENDIF
           ENDIF

C          GET CTF VOLUME
           CALL FILGET(FILNAM,FINPIC,NLEP,ILIST(K),IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           MAXIM = 0
           CALL OPFILEC(0,.FALSE.,FINPIC,LUN1,'O',IFORM,
     &             NSAM1,NROW1,NSLICE1,
     &             MAXIM,'DUMMY',.TRUE.,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           LRCL1 = NSAM1
C          ONLY 2D FOURIER

           INUMCTF = 0  !INUMCTF RECORDS 3D OR 2D: 3=3D 1=2D

           IF    (IFORM == -11) THEN
              NSAM1   = NSAM1-1
              INUMCTF = 1
           ELSEIF (IFORM == -12) THEN
              NSAM1   = NSAM1-2
              INUMCTF = 1
           ELSEIF (IFORM == -21) THEN
              NSAM1   = NSAM1-1
              INUMCTF = 3
           ELSEIF (IFORM == -22) THEN
              NSAM1   = NSAM1-2
              INUMCTF = 3
           ELSE
              CALL  ERRT(101,'CTF FILE MUST BE IN FOURIER FORMAT',NE)
              GOTO 9999
           ENDIF

           IF (INUMCTF .NE. INUMDEM) THEN
              CALL  ERRT(101,'CTF FILE MUST BE IN FOURIER FORMAT',NE)
              GOTO 9999
           ENDIF

c          VERIFY DIMENSIONS OF THE CTF
           IF (NSAM.NE.NSAM1.OR.NROW.NE.NROW1.OR.NSLICE.NE.NSLICE1)THEN
              CALL ERRT(101,'INCONSISTENT DIMENSIONS',NE)
              GOTO 9999
           ENDIF

C          MUTLIPLY INPUT VOLUME STORED IN THE MEMORY VOLIN, PLACE OUTPUT
C          IN VOLSUM AND SQUARED CTF IN CTF2.
           CALL MULFC3(LUN1,VOLIN,VOLSUM,CTF2,LRCL1/2,NROW,NSLICE)

        ENDDO              ! END OF DO-LOOP OVER ALL VOLUMES

C       DIVIDE BY THE SUM OF H^2+1/SNR
        VOLIN = VOLSUM / (REAL(CTF2) + S)   ! ARRAY OPERATION

        INV = -1
        CALL FMRS_3(VOLIN,NSAM,NROW,NSLICE,INV)

        IFORM = INUMCTF
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'N',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'OUTPUT',.FALSE.,IRTFLG)
        IF (IRTFLG.NE.0) GOTO 9999

C       WRITE THE OUTPUT
        DO L=1,NSLICE
           DO J=1,NROW
              CALL WRTLIN(LUN1,VOLIN(1,J,L),NSAM,J+(L-1)*NROW)
           ENDDO
        ENDDO


9999    IF (MYPID <= 0) CLOSE(LUN1)

        IF (ALLOCATED(CTF2))   DEALLOCATE(CTF2)
        IF (ALLOCATED(VOLSUM)) DEALLOCATE(VOLSUM)
        IF (ALLOCATED(VOLIN))  DEALLOCATE(VOLIN)
        IF (ALLOCATED(ILIST))  DEALLOCATE(ILIST)
        END



C++*********************************************************************
C
C MULFC3.F                         
C
C **********************************************************************
C
C  MULFC3(LUN1,VOLIN,VOLSUM, H2, NN2,NROW,NSLICE)
C
C  PARAMETERS:  LUN1      INPUT UNIT FOR CTF FILE              (INPUT)
C               VOLIN     FFT OF INPUT IMAGE                   (INPUT)
C               VOLSUM    VOLSUM + VOLIN * CONJG(B1)          (OUTPUT)
C               H2        H2 + B1^2                           (OUTPUT)
C               NN2       ROW LENGTH                           (INPUT)
C               NROW      IMAGE ROWS                           (INPUT)        
C               NSLICE    VOLUME SLICES                        (INPUT) 
C 
C  PURPOSE: USED IN CTF CORRECTION USING WIENER FILTERING APPROACH
C           READS CTF FILE FOR THIS FOR CURRENTLY OPEN VOLUME.   
C           ACCUMULATES: VOLIN    * CONJ(CTF) IN: VOLSUM
C           ACCUMULATES: REAL CTF * CONJ(CTF) IN: H2
C            
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE MULFC3(LUN1,VOLIN,VOLSUM, H2, NN2,NROW,NSLICE)

        COMPLEX  :: VOLIN(NN2,NROW,NSLICE)
        COMPLEX  :: VOLSUM(NN2,NROW,NSLICE)

        COMPLEX  :: H2(NN2,NROW,NSLICE)

        COMPLEX  :: B1(NN2)

        DO K=1,NSLICE
          DO J=1,NROW
            NR1 = J+(K-1)*NROW

C           READS CTF FOR THIS VOLUME      
            CALL REDLIN(LUN1,B1,2*NN2,NR1)

            DO I=1,NN2
               VOLSUM(I,J,K) = VOLSUM(I,J,K) + 
     &                         VOLIN(I,J,K) * CONJG(B1(I))

               H2(I,J,K) = H2(I,J,K)  + 
     &                     CMPLX(REAL(B1(I) * CONJG(B1(I))), 0.0)
            ENDDO
          ENDDO
        ENDDO

        END

 

