
C++*********************************************************************
C                                                                      *
C CCRS.F         NEW FROM CCRS_* FOR SPEEDUP       FEB 08 ARDEAN LEITH *
C                OUTPUT IN X                       APR 09 ARDEAN LEITH *
C                NSLCIE = 1                        NOV 11 ARDEAN LEITH *
C                                                                      *
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
C  CCRS(X,Y, LSC,NSAM,NROW,NSLICE, SPIDER_SIGN,SPIDER_SCALE, IRTFLG)
C
C  PURPOSE: CROSS-CORRELATE TWO INPUTS.
C
C  PARAMETERS: X             SAMPLE ROW/IMAGE/VOLUME 1            SENT 
C                            CORRELATED IMAGE ROW/IMAGE/VOLUME     RET.
C              Y             SAMPLE ROW/IMAGE/VOLUME 2            SENT
C              LSC           NSAM+2-MOD(NSAM,2))                   SENT
C              NSAM..        DIMENSIONS                            SENT
C              SPIDER_SIGN   CHANGE SIGN OF IMAGINARY TO SPIDER    SENT
C              SPIDER_SCALE  SCALE THE OUTPUT FLAG                 SENT
C              IRTFLG        ERROR FLAG                            RET.  
C
C  NOTE:   SHOULD CONTAIN ROW/IMAGE/VOLUME WITH A ROW 
C          LENGTH: NSAM+2-MOD(NSAM,2)).   FOR BEST ACCURACY IMAGES
C          SHOULD BE PADDED TO 2X ORIGINAL SIZE, THEN ROW LENGTH WILL
C          BE: 2*NSAM+2 (WHERE NSAM IS ORIGINAL IMAGE ROW LENGTH)
C
C--*********************************************************************

        SUBROUTINE CCRS(X,Y, LSC,NSAM,NROW,NSLICE,
     &                  SPIDER_SIGN,SPIDER_SCALE, IRTFLG)

        USE TYPE_KINDS
        INTEGER(KIND=I_8)    :: IPLAN = 0     !STRUCTURE POINTER 

        REAL, INTENT(INOUT)  :: X(LSC, NROW,NSLICE)
        REAL, INTENT(INOUT)  :: Y(LSC, NROW,NSLICE)
 
        INTEGER, INTENT(IN)  :: LSC,NSAM,NROW,NSLICE
        LOGICAL, INTENT(IN)  :: SPIDER_SIGN
        LOGICAL, INTENT(IN)  :: SPIDER_SCALE

        INTEGER, INTENT(OUT) :: IRTFLG

        REAL, PARAMETER      :: QUADPI = 3.141592653589793
        REAL, PARAMETER      :: PI2=2*QUADPI
 
        ITMPX = NSAM  / 2
        ITMPY = NROW  / 2
        ITMPZ = NSLICE / 2

        !write(6,*) ' itmpx,itmpy,itmpz:', itmpx,itmpy,itmpz 
        !write(6,*) ' itmpx,itmpy,itmpz:', itmpx,itmpy,itmpz 

        IF ( MOD(NSAM,2)   .EQ. 0  .AND. 
     &       MOD(NROW,2)   .EQ. 0  .AND. 
     &      (MOD(NSLICE,2) .EQ. 0  .OR.  NSLICE == 1)) THEN

C           EVEN NROW & NSLICE DIMENSIONS
C           CAN SKIP SIN & COS TESTS AS THEY ARE INTEGRALS OF PI
C           SIN(i*PI) IS ALWAYS ZERO. COS(i*PI) is either 1 or -1    

           NPISX = 2.0 * FLOAT(ITMPX) / FLOAT(NSAM)
           NPISY = 2.0 * FLOAT(ITMPY) / FLOAT(NROW)
           NPISZ = 2.0 * FLOAT(ITMPZ) / FLOAT(NSLICE)

           !write(6,*) ' npisx,npisy,npisz:',npisx,npisy,npisz 
   
           DO K=1,NSLICE
              IZ  = K - 1
              IF (IZ .GT. (NSLICE/2))  IZ = IZ - NSLICE
              NPISK = NPISZ * IZ

C$omp         parallel do private(j,iy,npisj,i,npis,icos,tmpr,tmpi)
              DO J=1,NROW
                 IY = J - 1
                 IF (IY .GT. (NROW/2)) IY = IY - NROW
                 NPISJ = NPISY * IY + NPISK

                 DO I=1,LSC,2
                    NPIS = NPISX * (I-1) / 2 + NPISJ
                    ICOS = 1
                    IF (MOD(NPIS,2) .NE. 0) ICOS = -1

                    !if (k .eq. nslice .and. j .eq. nrow) then
                    !write(6,*) ' npisk,npisj,npis:', npisk,npisj,npis,icos 
                    !endif

                    TMPR       =  X(I,J,K)   * Y(I,J,K)  + 
     &                            X(I+1,J,K) * Y(I+1,J,K)
                    TMPI       =  X(I+1,J,K) * Y(I,J,K)  - 
     &                            X(I,J,K)   * Y(I+1,J,K)
                    
                    X(I,J,K)   = TMPR * ICOS
                    X(I+1,J,K) = TMPI * ICOS
                 ENDDO
              ENDDO
           ENDDO
        ELSE
C          CAN NOT SKIP SIN & COS TESTS AS THEY ARE NOT INTEGRALS OF PI

           SX = PI2 * FLOAT(ITMPX) / FLOAT(NSAM)
           SY = PI2 * FLOAT(ITMPY) / FLOAT(NROW)
           SZ = PI2 * FLOAT(ITMPZ) / FLOAT(NSLICE)

           DO K=1,NSLICE
              IZ  = K - 1
              IF (IZ .GT. (NSLICE/2))  IZ = IZ - NSLICE
              ARGZ = SZ * IZ

              !write(6,*) 'iz,argz,k:',iz,argz,k
C$omp         parallel do private(j,iy,argy,i,arg,tmpr,tmpi)
              DO J=1,NROW
                 IY = J - 1
                 IF (IY .GT. (NROW/2)) IY = IY - NROW
                 ARGY = SY * IY + ARGZ
                 ! write(6,*) 'iy,argz,k:',iz,argz,k

                 DO I=1,LSC,2
                    ARG        = SX * (I-1) / 2 + ARGY

                    TMPR       =  X(I,J,K)   * Y(I,J,K)  + 
     &                            X(I+1,J,K) * Y(I+1,J,K)
                    TMPI       =  X(I+1,J,K) * Y(I,J,K)  - 
     &                            X(I,J,K)   * Y(I+1,J,K)

                    X(I,J,K)   = TMPR * COS(ARG) - TMPI * SIN(ARG)
                    X(I+1,J,K) = TMPI * COS(ARG) + TMPR * SIN(ARG)
                 ENDDO
              ENDDO
           ENDDO
        ENDIF

C       FOURIER OUTPUT, TRANSFORM IT BACK TO REAL
        INV = -1
        CALL FMRS(X, NSAM,NROW,NSLICE, IPLAN, 
     &            SPIDER_SIGN, SPIDER_SCALE, INV,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(101,FOURIER TRANSFORM FAILED,NE)
            RETURN
        ENDIF

        END



#ifdef NEVER
C       -------------------------- CCRS_SINEPLAN ---------------------
C       CURRENTLY UNUSED

        SUBROUTINE CCRS_SINEPLAN(X,Y, LSC,NSAM,NROW,NSLICE,
     &                           SPIDER_SIGN,SPIDER_SCALE, IRTFLG)

#ifdef SP_LIBFFTW3
        USE TYPE_KINDS
        INTEGER(KIND=I_8)    :: IPLAN = 0     !STRUCTURE POINTER 
#endif

        REAL, INTENT(IN)        :: X(LSC, NROW,NSLICE)
        REAL, INTENT(IN)        :: Y(LSC, NROW,NSLICE)
 
        INTEGER, INTENT(IN)     :: NSAM,NROW,NSLICE
        LOGICAL, INTENT(IN)     :: SPIDER_SIGN
        LOGICAL, INTENT(IN)     :: SPIDER_SCALE
        INTEGER, INTENT(OUT)    :: IRTFLG

        REAL, ALLOCATABLE, SAVE :: FSIN(:,:,:)

        REAL, PARAMETER         :: QUADPI = 3.141592653589793
        REAL, PARAMETER         :: PI2=2*QUADPI

        INTEGER, SAVE :: NSAMO=0,  NROWO=0,  NSLICEO=0

        IF (NSAM  .NE. NSAMO .OR. 
     &      NROW   .NE. NROWO .OR. 
     &      NSLICE .NE. NSLICEO) THEN
C           CREATE NEW SINE PLAN
           IF (ALLOCATED(FSIN)) DEALLOCATE(FSIN)

           write(6,*) ' Creating new sine plan '

           ALLOCATE(FSIN(LSC,NROW,NSLICE), STAT=IRTFLG)            
           IF (IRTFLG .NE. 0) THEN
              MWANT = LSC * NROW * NSLICE
              CALL ERRT(46,'FSIN',MWANT)
              RETURN
           ENDIF

           ITMPX = NSAM / 2
           SX    = PI2 * FLOAT(ITMPX) / FLOAT(NSAM)
           ITMPY = NROW / 2
           SY    = PI2 * FLOAT(ITMPY) / FLOAT(NROW)
           ITMPZ = NSLICE / 2
           SZ    = PI2 * FLOAT(ITMPZ) / FLOAT(NSLICE)

           DO K=1,NSLICE
              IZ  = K - 1
              IF (IZ .GT. (NSLICE/2))  IZ = IZ - NSLICE
              ARGZ = SZ * IZ

C$omp         parallel do private(i,j,iy,argy,arg)
              DO J=1,NROW
                 IY = J - 1
                 IF (IY .GT. (NROW/2)) IY = IY - NROW
                 ARGY = SY * IY + ARGZ

                 DO I=1,LSC,2
                    ARG           = SX * (I-1) / 2 + ARGY
                    FSIN(I,J,K)   = SIN(ARG)
                    FSIN(I+1,J,K) = COS(ARG)

#ifdef DEBUGD
                    it = mod((i-1),4)
                    temp = 1.0
                    if (it .eq. 0) temp = -1.0
                    
                    if ( j .eq. 50 .and. mod(i,1) .eq. 0) then
                    write(6,*)' i,arg,sin(arg),cos(arg),1-sin: ',
     &                          i,arg,sin(arg),cos(arg),temp  
                    endif
#endif
                  ENDDO
              ENDDO
           ENDDO
           NSAMO   = NSAM
           NROWO   = NROW
           NSLICEO = NSLICE
        ENDIF

        DO K=1,NSLICE
C$omp      parallel do private(i,j,sint,cost,tmpr,tmpi)
           DO J=1,NROW
              DO I=1,LSC,2
                 TMPR       =  X(I,J,K)   * Y(I,J,K)  + 
     &                         X(I+1,J,K) * Y(I+1,J,K)

                 TMPI       =  X(I+1,J,K) * Y(I,J,K)  - 
     &                         X(I,J,K)   * Y(I+1,J,K)

                 SINT       = FSIN((I),J,K)
                 COST       = FSIN((I+1),J,K)

                 X(I,J,K)   = TMPR * COST - TMPI * SINT
                 X(I+1,J,K) = TMPI * COST + TMPR * SINT
              ENDDO
           ENDDO
        ENDDO

C       FOURIER OUTPUT, TRANSFORM IT TO REAL
        INV = -1
        CALL FMRS(X, NSAM,NROW,NSLICE, IPLAN, 
     &            SPIDER_SIGN, SPIDER_SCALE, INV,IRTFLG)
        IF (IRTFLG .NE. 0) THEN
            CALL ERRT(101,FOURIER TRANSFORM FAILED,NE)
            RETURN
        ENDIF

        END
#endif








#ifdef DEBUGF
C       -------------------------- CCRS_1------------------------------

        SUBROUTINE CCRS_1(X,Y,O, LSC,NSAM,NROW,NSLICE,IRTFLG)

        REAL  X(LSC,NROW,NSLICE)
        REAL  Y(LSC,NROW,NSLICE)
        REAL  O(LSC,NROW,NSLICE)
 
        REAL, PARAMETER :: QUADPI = 3.1415926535897932384
        REAL, PARAMETER :: PI2=2*QUADPI

        ITMP = NSAM / 2
        SX   = PI2 * FLOAT(ITMP)/FLOAT(NSAM)
        ITMP = NROW / 2
        SY   = PI2 * FLOAT(ITMP )/ FLOAT(NROW)

           k = 1
C$omp      parallel do private(i,j,iy,argy,arg,tmpr,tmpi)
           DO J=1,NROW
              IY = J - 1
              IF (IY .GT. (NROW/2)) IY = IY - NROW
              ARGY = SY * IY 

              DO I=1,LSC,2
                 ARG  = SX * (I-1)/2 + ARGY

                 TMPR = X(I,J,1)   * Y(I,J,1)  + 
     &                  X(I+1,J,1) * Y(I+1,J,1)

                 TMPI = X(I+1,J,1) * Y(I,J,1)  - 
     &                  X(I,J,1)   * Y(I+1,J,1)

                 O(I,J,1)   = TMPR * COS(ARG) - TMPI * SIN(ARG)
                 O(I+1,J,1) = TMPI * COS(ARG) + TMPR * SIN(ARG)

                 if(j.eq.1.and.i.eq.1) then
                    sinta = SIN(ARG)
                    costa = COS(ARG)
                    write(6,*)' sinta:',sinta,' arg: ',arg
                    write(6,*)' costa:',costa
                    write(6,*)' O:    ',O(I,J,1),O(I+1,J,1)
                 endif
                 if(j.eq.20.and.i.eq.21) then
                   sinta = SIN(ARG)
                    costa = COS(ARG)
                    write(6,*)' sinta:',sinta,' arg: ',arg
                    write(6,*)' costa:',costa
                    write(6,*)' O:    ',O(I,J,1),O(I+1,J,1)
                 endif


              ENDDO
           ENDDO
        END
#endif

#ifdef NEVER
                    IF (SPIDER_SIGN) THEN
                       TMPR       =  X(I,J,K)   * Y(I,J,K)  + 
     &                               X(I+1,J,K) * Y(I+1,J,K)
                       TMPI       =  X(I+1,J,K) * Y(I,J,K)  - 
     &                               X(I,J,K)   * Y(I+1,J,K)
                    ELSE
                       TMPR       =  X(I,J,K)   * Y(I,J,K)  - 
     &                               X(I+1,J,K) * Y(I+1,J,K)
                       TMPI       =  - X(I+1,J,K) * Y(I,J,K)  = 
     &                                 X(I,J,K)   * Y(I+1,J,K)
                    ENDIF
 

IX,IY, bufi(ix,iy,1), XSHNEW,YSHNEW
 2*102,  176132.344,  -6.638765335E-4,  9.420156479E-3
  
  PEAKV, XSHNEW,YSHNEW
 175961.156,  6.638765335E-4,  -9.420156479E-3
 Testing Shift    --- SX :        0.00  SY :       -0.01

shud be_________________________________
 sinta: 0.E+0  arg:  0.E+0
  costa: 1.
  O:     0.E+0,  0.E+0
  sinta: 1.041018095E-6  arg:  91.1061859
  costa: -1.
  O:     -6.68330336,  -0.267273873
   
 IX,IY, bufi(ix,iy,1), XSHNEW,YSHNEW
 2*102,  176132.344,  -6.638765335E-4,  9.420156479E-3
  
  PEAKV, XSHNEW,YSHNEW
 175961.156,  6.638765335E-4,  -9.420156479E-3
 Testing Shift    --- SX :        0.00  SY :       -0.01
 T

#endif
