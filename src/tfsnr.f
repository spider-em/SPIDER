C **********************************************************************
C *  TFSNR.F
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2002, L. Joyeux & P. A. Penczek                       *
C=*                                                                    *
C=* University of Texas - Houston Medical School                       *
C=*                                                                    *
C=* Email:  pawel.a.penczek@uth.tmc.edu                                *
C=*                                                                    *
C=* This program is free software; you can redistribute it and/or      *
C=* modify it under the terms of the GNU General Public License as     *
C=* published by the Free Software Foundation; either version 2 of the *
C=* License, or (at your option) any later version.                    *
C=*                                                                    *
C=* This program is distributed in the hope that it will be useful,    *
C=* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
C=* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
C=* General Public License for more details.                           *
C=*                                                                    *
C=* You should have received a copy of the GNU General Public License  *
C=* along with this program; if not, write to the                      *
C=* Free Software Foundation, Inc.,                                    *
C=* 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.      *
C=*                                                                    *
C=**********************************************************************
C                         
C                  OPFILEC                         FEB  03 ARDEAN LEITH
C                  PROMPTS                         NOV  10 ARDEAN LEITH
C
C **********************************************************************

      SUBROUTINE  TFSNR

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
      INCLUDE 'F90ALLOC.INC'

      INTEGER, PARAMETER             :: NILMAX=998

      CHARACTER(LEN=MAXNAM)          :: FINPAT,FILNAM,FINPIC
      CHARACTER(LEN=MAXNAM)          :: FILPAT1,FILNM1,FILNMO

      INTEGER, ALLOCATABLE           :: ILIST(:)
      COMPLEX, ALLOCATABLE           :: VOLIN(:,:,:),VOLSUM(:,:,:)
      COMPLEX, ALLOCATABLE           :: CTF2(:,:,:)
      CHARACTER(LEN=1)               :: NULL
      LOGICAL                        :: FIRST
      REAL,  POINTER                 :: RPOINTER(:,:)

C     AUTOMATIC ARRAY
      INTEGER,PARAMETER              :: MAXSNR=998
      INTEGER                        :: IMAR(MAXSNR)
      
      DATA    LUN1/10/

      ALLOCATE(ILIST(NILMAX),STAT=IRTFLG)
      IF (IRTFLG.NE.0) THEN
         CALL  ERRT(46,'TF CTF, ILIST',NILMAX)
         GOTO 9999
      ENDIF

C     ASK FOR DATA FILE
      CALL FILELIST(.TRUE.,INPIC,FINPAT,NLET,ILIST,NILMAX,NIMA,
     &    'TEMPLATE FOR IMAGE FILES',IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      NULL = CHAR(0)

      CALL  FILERD(FILNAM,NLEP,NULL,'TEMPLATE FOR CTF',IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      CALL  FILERD(FILPAT1,NLET1,NULL,'SNR TEMPLATE',IRTFLG)
      IF (IRTFLG .NE. 0) GOTO 9999

      CALL RDPRM(WI, NOT_USED, 'RING WIDTH')

      CALL RDPRM(SNRF, NOT_USED,'FACTOR APPLIED ON SNR (0.5,1,...)')

      DO  K=1,NIMA

         CALL FILGET(FINPAT,FINPIC,NLET,ILIST(K),IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999

         MAXIM = 0
         CALL OPFILEC(0,.FALSE.,FINPIC,LUN1,'O',IFORM,NSAM1,NROW1,
     &             NSLICE1,MAXIM,'DUMMY',.TRUE.,IRTFLG)
         IF (IRTFLG.NE.0)  GOTO 9999
         
C         PROCESS THE VOLUME

         INUMDEM = 0
CC       USE INUMDEM TO RECORD 3D OR 2D
CC             INUMDEM = 3    3D
CC                     1    2D

         IF(IFORM.EQ.-11)  THEN
            LRCL=NSAM1
            LS=NSAM1
            NSAM1=NSAM1-1
              INUMDEM = 1
         ELSEIF(IFORM.EQ.-12)  THEN
            LRCL=NSAM1
            LS=NSAM1
            NSAM1=NSAM1-2
              INUMDEM = 1
         ELSEIF(IFORM.EQ.1) THEN
            LS=NSAM1+2-MOD(NSAM1,2)
            LRCL=NSAM1
              INUMDEM = 1
         ELSEIF(IFORM.EQ.-21)  THEN
            LRCL=NSAM1
            LS=NSAM1
            NSAM1=NSAM1-1
              INUMDEM = 3
         ELSEIF(IFORM.EQ.-22)  THEN
            LRCL=NSAM1
            LS=NSAM1
            NSAM1=NSAM1-2
              INUMDEM = 3
         ELSEIF(IFORM.EQ.3) THEN
            LS=NSAM1+2-MOD(NSAM1,2)
            LRCL=NSAM1
              INUMDEM = 3
         ELSE
            CLOSE(LUN1)
            CALL  ERRT(2,'TF CTF',NE)
            GOTO 9999
         ENDIF

C        GET THE DIMENSION INFORMATION FOR THE FIRST TIME
         IF(K.EQ.1) THEN
            NSAM=NSAM1
            NROW=NROW1
            NSLICE=NSLICE1
C           VOLUMES NEEDED

            ALLOCATE(VOLIN(LS/2,NROW,NSLICE),STAT=IRTFLG)
            IF(IRTFLG.NE.0) CALL  ERRT(46,'TF CTF, VOLIN',IER)

            ALLOCATE(VOLSUM(LS/2,NROW,NSLICE),STAT=IRTFLG)
            IF(IRTFLG.NE.0) CALL  ERRT(46,'TF CTF, VOLSUM',IER)

            ALLOCATE(CTF2(LS/2,NROW,NSLICE),STAT=IRTFLG)
            IF(IRTFLG.NE.0) CALL  ERRT(46,'TF CTF, CTF2',IER)

C           SET OUTPUT VOLUME AND CTF^2 TO ZERO
            VOLSUM=0.0
            CTF2=0.0
         ELSE

C           VERIFY THE DIMENSIONS
            IF(NSAM.NE.NSAM1.OR.NROW.NE.NROW1.OR.NSLICE.NE.NSLICE1) THEN
               CLOSE(LUN1)
               CALL  ERRT(1,'TF CTF',NE)
               GOTO 9999
            ENDIF
         ENDIF

C        READ INPUT FILE
         DO L=1,NSLICE
            DO J=1,NROW
               CALL  REDLIN(LUN1,VOLIN(1,J,L),LRCL,J+(L-1)*NROW)
            ENDDO
         ENDDO
         CLOSE(LUN1)


         IF ((IFORM.EQ.1).OR.(IFORM.EQ.3)) THEN
C           FOURIER TRANSFORM
            INV=+1
            CALL FMRS_3(VOLIN,NSAM,NROW,NSLICE,INV)
            IF (INV.EQ.0)  THEN
               CALL  ERRT(38,'TF CTF',NE)
               GOTO 9999
            ENDIF
         ENDIF


C        GET CTF
         CALL FILGET(FILNAM,FINPIC,NLEP,ILIST(K),IRTFLG)
         IF (IRTFLG .NE. 0) GOTO 9999
         
         MAXIM=0
         CALL OPFILEC(0,.FALSE.,FINPIC,LUN1,'O',IFORM,
     &          NSAM1,NROW1,NSLICE1,MAXIM,'DUMMY',.TRUE.,IRTFLG)
         IF (IRTFLG.NE.0) GOTO 9999

         LRCL1=NSAM1
C        ONLY 2D FOURIER

           INUMCTF = 0
CC       USE INUMCTF TO RECORD 3D OR 2D
CC             INUMCTF = 3    3D
CC                     1    2D

         IF (IFORM.EQ.-11) THEN
            NSAM1=NSAM1-1
              INUMCTF = 1
         ELSEIF (IFORM.EQ.-12) THEN
            NSAM1=NSAM1-2
              INUMCTF = 1
         ELSEIF (IFORM.EQ.-21) THEN
            NSAM1=NSAM1-1
              INUMCTF = 3
         ELSEIF (IFORM.EQ.-22) THEN
            NSAM1=NSAM1-2
              INUMCTF = 3
         ELSE
            CLOSE(LUN1)
            CALL  ERRT(2,'TF CTF',NE)
            GOTO 9999
         ENDIF

         IF (INUMCTF.NE.INUMDEM) THEN
             CLOSE(LUN1)
             CALL ERRT(2,'TF CTF',NE)
             GOTO 9999
         ENDIF

c        VERIFY DIMENSIONS OF THE CTF
         IF (NSAM.NE.NSAM1.OR.NROW.NE.NROW1.OR.NSLICE.NE.NSLICE1) THEN
            CLOSE(LUN1)
            CALL ERRT(1,'TF CTF',NE)
            GOTO 9999
         ENDIF
C         READ the snr
           CALL FILGET(FILPAT1, FILNM1, NLET1,ILIST(K), IRTFLG)
           CALL FILNAMANDEXT(FILNM1,DATEXC,FILNMO,NFN,.TRUE.,IRTFLG)
         MAXXT = 0
         MAXYT = 0
           CALL GETDOCDAT("Nothing", .FALSE., FILNMO, LUN2, .TRUE.,
     F         MAXXT, MAXYT, RPOINTER, IRTFLG)

C        MUTLIPLY INPUT VOLUME STORED IN THE MEMORY volin, PLACE OUTPUT
C        IN VOLSUM AND SQUARED CTF IN CTF2.
         CALL MULFC3_1(LUN1,VOLIN,VOLSUM,CTF2,LRCL1,NROW,NSLICE,
     &               RPOINTER,MAXXT,MAXYT,WI,SNRF)
         
         DEALLOCATE(RPOINTER)
         CLOSE(LUN1)
C        END OF DO-LOOP OVER ALL VOLUMES
      END DO ! K
   
C        DIVIDE BY THE SUM OF H^2+1/SNR
      VOLIN=VOLSUM/(REAL(CTF2)+1)
c     CALL  DIVCR3(B(KAD2),B(KAD2),B(KAD3),S,LRCL1/2,NROW,NSLICE)

      INV=-1
      CALL FMRS_3(VOLIN,NSAM,NROW,NSLICE,INV)

      IFORM = INUMCTF
      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FILNAM,LUN1,'N',IFORM,NSAM,NROW,NSLICE,
     &              MAXIM,'OUTPUT',.FALSE.,IRTFLG)
      IF(IRTFLG.NE.0) GOTO 9999

C     WRITE THE OUTPUT
      DO L=1,NSLICE
         DO J=1,NROW
            CALL  WRTLIN(LUN1,VOLIN(1,J,L),NSAM,J+(L-1)*NROW)
         ENDDO
      ENDDO
      CLOSE(LUN1)

9999  IF(ALLOCATED(CTF2)) DEALLOCATE(CTF2)
      IF(ALLOCATED(VOLSUM)) DEALLOCATE(VOLSUM)
      IF(ALLOCATED(VOLIN)) DEALLOCATE(VOLIN)
      IF(ALLOCATED(ILIST)) DEALLOCATE(ILIST)

      END
#if 1
CCCCCCCCCCCCCCCCCCCCCC   
   
      SUBROUTINE MULFC3_1(LUN1,A1,OUT,H2,NSAM,NROW,NSLICE,RPOINTER,
     F            MAXXT,MAXYT,WI,SNRF)
      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'F90ALLOC.INC'
      COMPLEX  A1(NSAM/2,NROW,NSLICE),OUT(NSAM/2,NROW,NSLICE)
      COMPLEX  H2(NSAM/2,NROW,NSLICE)
      COMPLEX  B1(NSAM/2)
      REAL      RPOINTER(MAXXT,MAXYT)
      INTEGER   RADIUS
C
      NS2 = NSLICE / 2
      NR2 = NROW / 2
      NN2 = NSAM / 2
      Y1 =  FLOAT(MAX0(NSAM,NROW,NSLICE))
        RPOINTER(2,:) = RPOINTER(2,:)*SNRF
      DO    K=1,NSLICE
         IZ = K-1
         IF (IZ .GT. NS2)  IZ = IZ - NSLICE
         PK = (FLOAT(IZ)/FLOAT(NS2))**2
         DO    J=1,NROW
            IY = J-1
            IF (IY .GT. NR2) IY = IY - NROW
               NR1=J+(K-1)*NROW
             CALL  REDLIN(LUN1,B1,2*NN2,NR1)
            PJ=(FLOAT(IY)/FLOAT(NR2))**2
            
            DO    I=1,NSAM,2
               III    = (I-1)/2
               PIII    = (FLOAT(III)/FLOAT(NN2))**2
               RR    = SQRT(PIII+PJ+PK)*0.5
               IR    = MIN0(MAX0(NINT(RR*Y1/WI)+1,1),MAXYT)
               IF(IR>MAXYT) THEN
                  WRITE(NOUT,*) 'ERROR:',IR,' is larger than: ', MAXYT
                  CALL ABORT
               END IF
                 FR     = RPOINTER(2,IR)
               OUT(III,J,K) = OUT(III,J,K)
     &               +A1(III,J,K)*CONJG(B1(III))*FR
                H2(III,J,K) = H2(III,J,K)
     &               +CMPLX(REAL(B1(III)*CONJG(B1(III))),0.0)*FR
            END DO
          END DO
      END DO
      
      
      END
#endif
