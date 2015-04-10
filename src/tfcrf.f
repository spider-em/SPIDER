C **********************************************************************
C *  TFCRF
C=**********************************************************************
C=* From: SPIDER - MODULAR IMAGE PROCESSING SYSTEM                     *
C=* Copyright (C)2001, P. A. Penczek                                   *
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
C  PURPOSE:  DETERMINE DEFOCUS AND AMPL CONTR. RATIO BASED ON 
C            CROSS-RESOLUTION CURVE
C
C  PARAMETERS: NONE    
C
C **********************************************************************

        SUBROUTINE  TFCRF

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'F90ALLOC.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL, DIMENSION(:,:), POINTER   :: PANG
        REAL, ALLOCATABLE, DIMENSION(:) :: FREQ,CCF,OCCF
        REAL ::                            LAMBDA

        CHARACTER(LEN=MAXNAM)           :: FILE

        PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)

        CALL RDPRM2(P,CS,NOT_USED,
     &      'PIXEL SIZE[A] & SPHERICAL ABERRATION CS[MM]')

        IF (CS < 0.0001)    CS = 0.0001

C       CONVERT CS TO [A]
        CS = CS*1.0E07

        CALL RDPRM2(DZ,LAMBDA,NOT_USED,
     &    'DEFOCUS(ANGSTROEMS), WAVELENGTH LAMBDA(ANGSTROMS)[A]')

        CALL RDPRM(WGH,NOT_USED,'AMPL. CONTRAST RATIO [0-1]')

C       RETRIEVE ARRAY WITH ANGLES DATA IN IT
        MAXXT = 0
        MAXYT = 0
        CALL GETDOCDAT('CROSS-RESOLUTION DOC',.TRUE.,FILE,
     &                       77,.TRUE.,MAXXT,MAXYT,PANG,IRTFLG)
 
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(4,'TF CRF ',NE)
           RETURN
        ENDIF    
   
        ALLOCATE (FREQ(MAXYT),CCF(MAXYT),OCCF(MAXYT), STAT=IRTFLG)
        IF (IRTFLG.NE.0) THEN 
           CALL ERRT(46,'TF CRF; FREQ, CCF & OCCF',IER)
           RETURN
        ENDIF

C       RETAIN SIGNS OF THE CCF
	L=0
        DO  I=1,MAXYT
	  IF(PANG(1,I) .GT. 0.0)  THEN
	     L       = L+1
             FREQ(L) = PANG(2,I)/P
             CCF(L)  = SIGN(1.0,PANG(4,I))
             OCCF(L) = PANG(4,I)**2
	  ENDIF
        ENDDO

        DEALLOCATE(PANG)

C       DETERMINE DEFOCUS AND PHASE BY MATCHING SIGNS
        AA   = 0.5*QUADPI*CS*LAMBDA**3
        BB   = QUADPI*LAMBDA
        CZRM = 0.0
        
C       WAS ESTIMATION OF THE AMPLITUDE CONTRAST REQUESTED?
        IF (WGH .EQ. 0.0)  THEN
           LAN = 199
        ELSE
           LAN = 0
        ENDIF

C       WAS ESTIMATION OF THE DEFOCUS REQUESTED?
        IF (DZ .EQ. 0.0)  THEN
           MDZN = 50000
        ELSE
           MDZN=0
        ENDIF

        DO  LA=0,LAN
           IF(LA.EQ.0)  THEN
              WGHT=WGH
           ELSE
              WGHT=REAL(LA)/200
           ENDIF
C          GET OFFSET A
           A=ATAN(WGHT/(1.0-WGHT))
        
           DO  MDZ=0,MDZN,10
              IF(MDZ.EQ.0)  THEN
                 DZT=DZ
              ELSE
                 DZT=MDZ 
              ENDIF
C
              CZR=0.0
              DO  I=1,L
                AK=FREQ(I)
                CZR=CZR+OCCF(I)*
     &          ABS(CCF(I)-SIGN(1.0,SIN((AA*AK**2-BB*DZT)*AK**2-A)))/2.0
              ENDDO
              IF (CZR.GT.CZRM)  THEN
                 CZRM=CZR
                 DEFOCUS=DZT
                 WGHM=WGHT
                 DEFOCUSN=DZT
                 WGHN=WGHT
              ELSEIF (CZR.EQ.CZRM)  THEN
C                CATCH THE UPPER LIMIT OF BOTH VALUES
                 DEFOCUSN=DZT
                 WGHN=WGHT
              ENDIF
          ENDDO
        ENDDO

C        CONTROL PRINTOUT, NOT NEEDED
C        A=ATAN(WGHM/(1.0-WGHM))
C        RAT=SQRT((1.0-WGHM)**2+WGHM**2)
C          DO  I=1,L
C         AK=FREQ(I)
C        PRINT  541,I,CCF(I)
C     &  ,SIGN(1.0,SIN((AA*AK**2-BB*DEFOCUS)*AK**2-A))
C     &  ,SIN((AA*AK**2-BB*DEFOCUS)*AK**2-A)*RAT
C     & ,(1-WGHM)*SIN((AA*AK**2-BB*DEFOCUS)*AK**2)-
C     &  WGHM*COS((AA*AK**2-BB*DEFOCUS)*AK**2)
C 541	FORMAT(1X,I4,4F10.6)
C        ENDDO

        DEALLOCATE(FREQ,CCF,OCCF)

C       RETURN IN REGISTERS DEFOCUS, WGHM, DEFOCUSN, WGHN, AND CZRM
        CALL REG_SET_NSEL(1,5,DEFOCUS,WGHM,DEFOCUSN,WGHN,CZRM,IRTFLG)

        WRITE(NOUT,654)  DEFOCUS,WGHM,DEFOCUSN,WGHN,CZRM
654     FORMAT(2('DEFOCUS: ',F10.2,' AMPL. CONTRAST RATIO: ',F4.2,/),
     &          ' ERROR: ',F14.6)

        END
