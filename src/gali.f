C++*********************************************************************
C
C    GALI.F           DYNAMIC MEMORY ALLOCATION SEPT 2000 BIMAL RATH 
C                     REPLACED 'INCORE'         SEPT 2001 ARDEAN LEITH
C                     OPFILEC                   FEB  2003 ARDEAN LEITH
C                     GALI_P INSERTED           MAR  2012 ARDEAN LEITH
C                     BAD LUNIN                 DEC  2013 ARDEAN LEITH
C                     LEN=MAXNAM                JUL  2014 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2017, Health Research Inc.,                         *
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
C   GALI
C
C   CALLS:  GALI_P 
C           ALROSI_Q 
C           ALROSF_Q 
C           RTQS_Q 
C           CROSRNG_Q
C           SHFI_Q 
C           BLOB_Q 
C           FINDMX_Q 
C           ALRQ_Q
C           ALPRBS_Q 
C           FOURING_Q 
C           RTQ_Q 
C           CENT_Q 
C           OUTIM_Q 
C           FMRS_1 
C           FMRS_1D 
C           QUADRI_Q
C
C   NOTE: SOMETIMES FAILS, IS BUGGY, SHOULD BE REWRITTEN al
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE GALI

        INCLUDE 'CMBLOCK.INC' 
        INCLUDE 'CMLIMIT.INC'
 
        INTEGER, ALLOCATABLE   :: NUMR(:,:)
        REAL,    ALLOCATABLE   :: BLOB(:,:)

        CHARACTER(LEN=MAXNAM)  :: FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC
        COMMON /FISPEC/
     &     FINPAT,FINPIC,DOCFIL,OUTIMA,OUTDOC,NLET,NLETI,NLIMA,NLDOC

        COMMON /MXR/     MAXRIN
C       MXR in: ang.f, gali.f, hali_p.f

        LOGICAL                :: FOUROK = .FALSE.
        CHARACTER*1            :: MODE
        LOGICAL                :: USEBLOB

        CHARACTER(LEN=1)       :: NULL = CHAR(0)

        INTEGER,PARAMETER      :: INPIC  = 78
        INTEGER,PARAMETER      :: LUNDOC = 80   
        INTEGER,PARAMETER      :: LUNXM  = 0  ! SELFILE NOT ALLOWED
        INTEGER                :: IMGNUM

C       OPEN INPUT IMAGE(S)
        MAXIM = 0
        NILMAX = NIMAX
        IMGNUM = 0    ! Passing uninitialized variable may result in 'INVALID IMAGE NUMBER' error
        CALL OPFILES(0,INPIC,LUNDOC,LUNXM,
     &             .TRUE.,FINPAT,NLET, 'O',
     &             IFORM,NX,NY,NZ,MAXIM,
     &             'INPUT FILE TEMPLATE (E.G. PIC****)~',
     &             FOUROK,INUMBR,NILMAX, 
     &             NDUM,NIMA,IMGNUM, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN
        CLOSE(INPIC)


C       NIMA - TOTAL NUMBER OF IMAGES
        IF (NIMA > 0)  THEN
           WRITE(NOUT,2001) NIMA
2001       FORMAT('  Number of images: ',I5)
        ELSE
           CALL ERRT(101,'NO IMAGES',NDUM)
           GOTO 9999
        ENDIF

C       FOURIER SIZE
        LSD = NX + 2 - MOD(NX,2)

        CALL RDPRI1S(NSI,NOT_USED,'EXPECTED SIZE OF THE OBJECT',IRTFLG)
        IF (NSI >= NX)  THEN
           WRITE(NOUT,*)' OBJECT SIZE CANNOT BE LARGER THAN IMAGE SIZE'
           NSI = NX - 2
           WRITE(NOUT,*) ' OBJECT SIZE SET TO: ',NSI
        ENDIF
        NSI = MAX(1,(NX-NSI)/2)

        CALL RDPRMI(MR,NR,NOT_USED,'FIRST AND LAST RING RADIUS')

        IF (MR.LE.0 .OR. NR .GE. MIN((NX/2),(NY/2))) THEN
           CALL ERRT(31,'AP SR',NE)
           GOTO 9999
        ENDIF

        IF (NR .GT. NX/2-1)  THEN
           NR = NX/2-1
           WRITE(NOUT,*) ' LAST RING RADIUS LIMITED TO: ',NR
        ENDIF
        MR    = MAX(1,MIN(NR,MR))
        ISKIP = 1
        MODE  = 'F'

C       FIND TOTAL NUMBER OF RINGS
        NRING = 0
        DO I=MR,NR,ISKIP
           NRING = NRING + 1
        ENDDO

        ALLOCATE (NUMR(3,NRING), STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           CALL ERRT(46,'GALI; NUMR',3*NRING)
           GOTO 9999
        ENDIF

C       FILL RINGS POINTER
        NRING = 0
        DO I=MR,NR,ISKIP
           NRING         = NRING + 1
           NUMR(1,NRING) = I
        ENDDO

C       FILL NUMR
        CALL ALPRBS_Q(NUMR,NRING,LCIRC,MODE)

        MAXRIN = NUMR(3,NRING) - 2  ! WHY -2 ??


C       CHECK WHETHER USER HAS OWN REFERENCE FILE TO CENTER THE AVERAGE.
C       DISP = 'Z' WILL NOT CALL ERRT IN OPFIL IF NOT EXISTING
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,OUTIMA,INPIC,'Z',IFORM,NX,NY,NZ,
     &    MAXIM,'IMAGE TO BE USED TO CENTER THE AVERAGE',.FALSE.,IRTFLG)

        IF (IRTFLG .NE. 0)  THEN
           WRITE(NOUT,*) ' NO IMAGE GIVEN, DEFAULT PROCEDURE USED'
           USEBLOB = .TRUE.
        ELSE
           USEBLOB = .FALSE.
        ENDIF

C       TEMPLATES FOR OUTPUT FILES
        CALL  FILERD(OUTIMA,NLIMA,NULL,
     &        'TEMPLATE FOR AVERAGE FILES',IRTFLG)

        CALL  FILERD(OUTDOC,NLDOC,NULL,
     &        'TEMPLATE FOR ALIGNMENT DOC FILES',IRTFLG)
	
        JACUP = 0  ! UNKNOWN PURPOSE

        ALLOCATE(BLOB(LSD,NY),STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL  ERRT(46,'AP SR; BLOB',LSD*NY)
           GOTO 9999
        ENDIF

        IF (.NOT.USEBLOB)  THEN 
           CALL READV(INPIC,BLOB,LSD,NY,NX,NY,NZ)
           CLOSE(INPIC)
        ENDIF

        CALL GALI_P(LSD,NX,NY,NSI,INUMBR,NIMA,MODE,JACUP,
     &           LCIRC,NUMR,NRING,MAXRIN,BLOB,USEBLOB,NR,NOUT,ITER)
 
        CALL REG_SET_NSEL(1,1,FLOAT(ITER),0.0,0.0,0.0,0.0,IRTFLG)

        !!WRITE (NOUT,2600)
2600    FORMAT (/'  ',12('-'),' END OF COMPUTATION ',12('-')/)

9999    IF (ALLOCATED(BLOB))  DEALLOCATE(BLOB)
        IF (ALLOCATED(NUMR))  DEALLOCATE(NUMR)
                
        END



C++*********************************************************************
C
C GALI_P.F            DYNAMIC MEMORY ALLOCATION SEP 2000 BIMAL RATH    
C                     REPLACED 'INCORE'         SEP 2001 ARDEAN LEITH
C                     ADDED ITER PARAM          MAR 2012 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2017  Health Research Inc.,                         *
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
C GALI_P
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE GALI_P(LSD,NX,NY,NSI,ILIST,NIMA,MODE,JACUP,
     &             LCIRC,NUMR,NRING,MAXRIN,BLOB,USEBLOB,NR,NOUT,NITER)

        REAL, ALLOCATABLE, DIMENSION(:,:) :: A, B, C, D, E
        REAL, ALLOCATABLE, DIMENSION(:,:) :: REFER, REFERN, REFERTMP
        REAL, ALLOCATABLE, DIMENSION(:,:) :: PARA
        REAL, ALLOCATABLE, DIMENSION(:)   :: A_CIRC, REFER_CIRC
        REAL, ALLOCATABLE, DIMENSION(:)   :: QIMAGES

        INTEGER    :: ILIST(NIMA)
        REAL       :: BLOB(LSD,NY)
        LOGICAL    :: CHANGE, NOCHANGE, USEBLOB, IMAGE(NIMA), INCORE
        INTEGER    :: NUMR(3,NRING)

C       TEMP IS AN AUTOMATIC ARRAY
        DOUBLE PRECISION  :: TEMP(MAXRIN+2,2),FNRM,DNRM

        CHARACTER*1       :: MODE

        DATA  ZERO/0.0/
        DATA  INPIC/77/

        NSNR = LSD * NY
        
        ALLOCATE(A(LSD,NY), B(LSD,NY),C(LSD,NY), D(LSD,NY),
     &           E(LSD,NY), REFER(LSD,NY), REFERN(LSD,NY),
     &           REFERTMP(LSD,NY), PARA(3,NIMA), A_CIRC(LCIRC),
     &           REFER_CIRC(LCIRC),    STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           MWANT = 8*LSD*NY + 3*NIMA + 2*LCIRC
           CALL  ERRT(46,'GALI_P; A,B,C,D,....',MWANT)
           GOTO 9999
        ENDIF

        INCORE = .FALSE.
        ISIZE  = NX * NY * NIMA
        ALLOCATE(QIMAGES(ISIZE), STAT=IRTFLG)
        IF (IRTFLG == 0) THEN
           WRITE(NOUT,*) ' -- IMAGES IN_CORE VERSION -- '
C          LOAD IMAGE DATA INTO QIMAGES TO SPEED THINGS UP
           DO IMI = 1,NIMA
              CALL GETIMA_Q(QIMAGES(NX*NY*(IMI-1)+1),NX,NX,NY,
     &                      IMI,ILIST,NIMA,QIMAGES,INPIC,INCORE)
           ENDDO
           INCORE = .TRUE.
        ELSE
           WRITE(NOUT,*) ' -- IMAGES FROM DISK VERSION -- '
        ENDIF

C       PREPARE 'BLOB'  TO CENTER INPUT IMAGES ....
        IF (USEBLOB) CALL BLOB_Q(BLOB,LSD,NX,NY,NR)
	
        INS = 1
        CALL FMRS_2(BLOB,NX,NY,INS)
        IF (INS == 0)  THEN
           CALL  ERRT(38,'AP SR',NE)
           GOTO 9999
        ENDIF

C       BUILD FIRST AVERAGE

C       IMAGE  IS USED HERE FOR THE RANDOM SELECTION OF IMAGES
C       INITIALIZE ALL IMAGE() TO .FALSE.

        IMAGE = .FALSE.
        CALL RANDOM_NUMBER(CIID)
        IMI         = MIN0(NIMA,MAX0(1,INT(CIID*NIMA+0.5)))

        IMAGE(IMI)  = .TRUE.
        PARA(1,IMI) = 0.0
        PARA(2,IMI) = 0.0
        PARA(3,IMI) = 0.0

C       GET IMAGE DATA
        CALL GETIMA_Q(A,LSD,NX,NY,IMI,ILIST,NIMA,
     &                QIMAGES,INPIC,INCORE)

C       COPIES IMAGE DATA FROM A INTO B (PARALLEL)
        CALL COP(A,B,NSNR)

C       FFT ON IMAGE DATA IN ARRAY B
        INS = 1
        CALL FMRS_2(B,NX,NY,INS)

C       BLOB - BLOB; B-CURRENT IMAGE; C-CCF.

C       CROSS CORRELATE BLOB WITH B
        LSC = NX+2-MOD(NX,2)
        CALL CCRS_2(BLOB,B,C, LSC,NX,NY)

C       FIND PEAK IN CROSS CORRELATION
        CALL FINDMX_Q(C,LSD,NX,NY,NSI,CMX1,SX1,SY1)

C       SHIFT A INTO REFERN
        ISX1        = SX1
        ISY1        = SY1
        PARA(2,IMI) = ISX1
        PARA(3,IMI) = ISY1
        CALL SHFI_Q(A,REFERN,LSD,NX,NY,ISX1,ISY1)

        !PRINT  *,' FIRST IMAGE WITH BLOB',IMI,ISX1,ISY1

C       COPY REFERN INTO REFER
        CALL COP(REFERN,REFER,NSNR)

C       INTERPOLATE REFER INTO POLAR COORDINATES IN REFER_CIRC
        CALL ALRQ_Q(REFER,LSD,NX,NY,NUMR,REFER_CIRC,LCIRC,
     &                NRING,MODE,IPIC)

C       FFT OF POLAR RINGS
        CALL FOURING_Q(REFER_CIRC,LCIRC,NUMR,NRING,TEMP,MODE)

        !write(6,*)' In gali_p, refer1:', refer(61,34)

C       FORWARD FFT ON REFER
        INS = 1
        CALL FMRS_2(REFER,NX,NY,INS)


C       GO THROUGH ALL THE REMAINING IMAGES.
        CHANGE = .FALSE.
        DO KTN = 2,NIMA
           DO
              CALL RANDOM_NUMBER(CIID)
              M   = MIN0(NIMA,MAX0(1,INT(CIID*(NIMA-KTN+1)+0.5)))
              IMI = 0
              DO  I=1,NIMA
                 IF (.NOT.IMAGE(I))  THEN
                    IMI = IMI + 1
                    IF (IMI == M)  THEN
                       IMI = I
                       GOTO  801
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO

801        IMAGE(IMI) = .TRUE.
C          GET IMAGE DATA
           CALL GETIMA_Q(A,LSD,NX,NY,IMI,ILIST,NIMA,
     &                   QIMAGES,INPIC,INCORE)

C          COPY IMAGE DATA ARRAY A INTO ARRAY B
           CALL COP(A,B,NSNR)

C          FORWARD FFT ON IMAGE DATA IN ARRAY B
           INS = 1
           CALL FMRS_2(B,NX,NY,INS)

C          BLOB - BLOB; B-CURRENT IMAGE; C-CCF.

C          CROSS CORRELATE BLOB WITH B INTO C
           LSC = NX+2-MOD(NX,2)
           CALL CCRS_2(BLOB,B,C, LSC,NX,NY)

C          FIND PEAK IN CROSS CORRELATION IMAGE C

           CALL FINDMX_Q(C,LSD,NX,NY,NSI,CMX1,SX1,SY1)

           ISX1        = SX1
           ISY1        = SY1
           PARA(1,IMI) = 0.0
           PARA(2,IMI) = ISX1
           PARA(3,IMI) = ISY1

C          SHIFT IMAGE A INTO B
           CALL SHFI_Q(A,B,LSD,NX,NY,ISX1,ISY1)

           !PRINT  *,' NEXT IMAGE WITH BLOB',IMI, ISX1,ISY1

C          RUN ROTATION-SHIFT ALIGNMENT
C          INPUT:   THREE REAL IMAGES, A, B AND REFERN
C          OUTPUT:  B - REAL ALIGNED IMAGE, PARA - UDATED PARAMETERS
C          SCRATCH: D
C          NOCHANGE = .TRUE. - NO CHANGE IN THE POSITION OF THE IMAGE
C
C          CONTENTS:
C          BLOB   : IMAGE FOR CENTERING AVERAGE *
C          A      : IMAGE DATA
C          B      : COPY OF IMAGE DATA --> FFT OF IMAGE DATA 
C                   THEN SHIFTED IMAGE A
C          D      : ??
C          C      : CROSS CORRELATION OF BLOB WITH B
C          REFER  : FFT OF SHIFTED IMAGE
C          REFERN : SHIFTED IMAGE 

           CALL ALROSI_Q(A,B,D,C,REFER,LSD,NX,NY,NSI,
     &          PARA(1,IMI),NOCHANGE,
     &          A_CIRC,REFER_CIRC,LCIRC,JACUP,NUMR,NRING,MAXRIN,
     &          TEMP,MODE,KTN)
           IF (.NOT. NOCHANGE)  CHANGE = .TRUE.

C          WEIGHTED ADD SHIFTED IMAGE B TO THE SHIFTED REFERENCE REFERN
C          CALL UPDTF  (REFERN,B,NSNR,KTN)
           CALL UPDTF_R(REFERN,B,NX,NY,LSD,KTN)

        ENDDO

C       CORRECT CENTER

        IF (USEBLOB)  THEN
           CALL CENT_Q(REFERN,LSD,NX,NY,XS,YS)
           XS = -(XS-NX/2-1)
           YS = -(YS-NY/2-1)
        ELSE
C          BLOB - BLOB; REFERN-CURRENT AVERAGE; D-CCF.

C          COPY REFERN INTO REFER
           CALL COP(REFERN,REFER,NSNR)

C          FORWARD FFT ON REFER
           INS = +1
           CALL FMRS_2(REFER,NX,NY,INS)

C          CROSS CORRELATE BLOB WITH REFER INTO D
           LSC = NX+2-MOD(NX,2)
           CALL CCRS_2(BLOB,REFER,D, LSC,NX,NY)

C          FIND PEAK IN D
           CALL FINDMX_Q(D,LSD,NX,NY,NSI,CMX1,XS,YS)
        ENDIF

C       PRINT *, '  CENTER CORRECTION  ',XS,YS

C       ROTATE AND SHIFT REFERN INTO REFER
        CALL RTQS_Q(REFERN,REFER,LSD,NX,NY,ZERO,XS,YS)
        DO IMI=1,NIMA
           PARA(2,IMI) = PARA(2,IMI)+XS
           PARA(3,IMI) = PARA(3,IMI)+YS
        ENDDO

C       SAVE IMAGE FROM REFER
        CALL OUTIM_Q(REFER,LSD,NX,NY,1)

        CALL OUTPR(PARA,NIMA,1)

        TEMP(1,1) = FNRM(REFER,NSNR)

        I = 1
        WRITE(NOUT,5091)  I,TEMP(1,1)
5091    FORMAT('  ITERATION:',I5,'   SUM OF SQUARES:',1PD12.5)

        NITER = I

C       IF NO CHANGES TERMINATE
        IF (.NOT. CHANGE) GOTO 9999

C       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       REFINE THE ALIGNMENT

        ITER = 0
        DNRM = 0.0D0

C       ITERATE THE REFINEMENT
        DO
           CHANGE = .FALSE.
           ITER   = ITER + 1

c$omp      parallel do private(j,i)
           DO J=1,NY
              DO I=1,LSD
                 REFERN(I,J) = 0.0
              ENDDO
           ENDDO
c$omp      end parallel do
	   
c$omp      parallel do private(imi)
           DO IMI=1,NIMA
              IMAGE(IMI) = .FALSE.
           ENDDO
c$omp      end parallel do

           DO KTN=1,NIMA
              DO
                 CALL RANDOM_NUMBER(CIID) 
                 M   = MIN0(NIMA,MAX0(1,INT(CIID*(NIMA-KTN+1)+0.5)))
                 IMI = 0
                 DO I=1,NIMA
                    IF (.NOT. IMAGE(I))  THEN
                       IMI = IMI+1
                       IF (IMI == M)  THEN
                          IMI = I
                          GOTO  901
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO

901           IMAGE(IMI) = .TRUE.

              !PRINT *,' IMAGE  #  ',IMI
C             GET IMAGE DATA IN ARRAY A
              CALL GETIMA_Q(A,LSD,NX,NY,IMI,ILIST,NIMA,
     &                      QIMAGES,INPIC,INCORE)

              CALL RTQS_Q(A,C,LSD,NX,NY,PARA(1,IMI),
     &                    PARA(2,IMI),PARA(3,IMI))

C             SUBTRACT THE IMAGE FROM THE REFERENCE
              CALL SUBAF(REFER,C,REFERTMP,NSNR,NIMA)

C             COPY REFERTMP INTO E
              CALL COP(REFERTMP,E,NSNR)

C             RUN ROTATION-SHIFT ALIGNMENT
C             INPUT:  TWO REAL IMAGES, B AND REFER
C             OUTPUT: C - REAL ALIGNED IMAGE, PARA - UDATED PARAMETERS
C             NOCHANGE=.TRUE. - NO CHANGE IN THE POSITION OF THE IMAGE

              CALL  ALROSF_Q(A,C,D,B,REFERTMP,LSD,NX,NY,NSI,
     &                       PARA(1,IMI),NOCHANGE,
     &                       A_CIRC,REFER_CIRC,LCIRC,JACUP,
     &                       NUMR,NRING,MAXRIN,TEMP,MODE)

              IF (.NOT. NOCHANGE)  THEN
                 CHANGE = .TRUE.

C                ADD IMAGE TO THE REFERENCE
                 CALL UPDF(E,C,REFER,NSNR,NIMA)
              ENDIF
C              CALL UPDTF(REFERN,C,NSNR,KTN)
               CALL UPDTF_R(REFERN,C,NX,NY,LSD,KTN)
           ENDDO

C          CORRECT CENTER

           IF (USEBLOB)  THEN
              CALL CENT_Q(REFERN,LSD,NX,NY,XS,YS)
              XS = -(XS-NX/2-1)
              YS = -(YS-NY/2-1)
           ELSE

C             BLOB - BLOB; REFERN-CURRENT AVERAGE; D-CCF.

C             COPY REFERN INTO REFER
              CALL COP(REFERN,REFER,NSNR)

C             FORWARD FFT ON REFER
              INS = +1
              CALL FMRS_2(REFER,NX,NY,INS)

C             CROSS CORRELATE BLOB WITH REFER INTO D
              LSC = NX+2-MOD(NX,2)
              CALL CCRS_2(BLOB,REFER,D, LSC,NX,NY)

C             FIND PEAK IN D
              CALL FINDMX_Q(D,LSD,NX,NY,NSI,CMX1,XS,YS)
           ENDIF

C          ROTATE AND SHIFT REFERN INTO REFER
           CALL RTQS_Q(REFERN,REFER,LSD,NX,NY,ZERO,XS,YS)

c$omp      parallel do private(imi)
           DO IMI=1,NIMA
              PARA(2,IMI) = PARA(2,IMI) + XS
              PARA(3,IMI) = PARA(3,IMI) + YS
           ENDDO
c$omp      end parallel do

C          SAVE REFER IMAGE
           CALL OUTIM_Q(REFER,LSD,NX,NY,ITER+1)

           CALL OUTPR(PARA,NIMA,ITER+1)

           TEMP(1,1) = FNRM(REFER,NSNR)
           WRITE(NOUT,5091)  ITER+1,TEMP(1,1)
           NITER = ITER + 1

           IF (NITER > 100) THEN
              CALL ERRT(101,
     &         'POSSIBLE ENDLESS LOOP AFTER, ITERATIONS:',NITER)
           GOTO 9999
        ENDIF

           IF (TEMP(1,1) >= DNRM)  THEN
              DNRM = TEMP(1,1)
              IF (.NOT.CHANGE)  GOTO 9999
           ELSE
              GOTO 9999
           ENDIF

C          END OF INFINITE ITERATIONS
        ENDDO

9999    IF (ALLOCATED(A))          DEALLOCATE(A)
        IF (ALLOCATED(B))          DEALLOCATE(B)
        IF (ALLOCATED(C))          DEALLOCATE(C)
        IF (ALLOCATED(D))          DEALLOCATE(D)
        IF (ALLOCATED(E))          DEALLOCATE(E)
        IF (ALLOCATED(REFER))      DEALLOCATE(REFER)
        IF (ALLOCATED(REFERN))     DEALLOCATE(REFERN)
        IF (ALLOCATED(REFERTMP))   DEALLOCATE(REFERTMP)
        IF (ALLOCATED(PARA))       DEALLOCATE(PARA)
        IF (ALLOCATED(A_CIRC))     DEALLOCATE(A_CIRC)
        IF (ALLOCATED(REFER_CIRC)) DEALLOCATE(REFER_CIRC)
        IF (ALLOCATED(QIMAGES))    DEALLOCATE(QIMAGES)

        END
