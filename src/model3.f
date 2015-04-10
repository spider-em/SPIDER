C++*********************************************************************
C
C MODEL3.F                         DOCUMENT FILE OPTION FEB   88 JF
C                                  FILENAMES LENGTHENED DEC   88 AL
C                                  CHAR. VARIABLES      AUG   89 AL
C                                  REWRITTEN SOME       APR   97 AL
C                                  HELIX BUGS FIXED     AUG   02 AL
C                                  RDPRAF REMOVED       DEC   05 AL
C                                  NORMAL GAUSSIAN      APR   09 AL
C                                  CYL BUG              SEP   09 AL
C                                  G2..                 JAN   12 AL
C
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
C    MODEL3(LUN1,NDOC,FILNAM,NX,NY,NZ)
C
C    PURPOSE:   PREPARES 3-D MODEL DENSITY DISTRIBUTIONS
C
C    PARAMETERS:       LUN1        LOGICAL UNIT NUMBER OF VOLUME FILE
C                      NDOC        LOGICAL UNIT NUMBER OF DOCUMENT FILE
C                      FILNAM      NAME OF VOLUME FILE
C                      NX,NY,NZ    DIMENSIONS OF FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--********************************************************************

	SUBROUTINE MODEL3(LUN1,NDOC,FILNAM,NX,NY,NZ)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        INTEGER          :: LUN1,NDOC,NX,NY,NZ
        CHARACTER(LEN=*) :: FILNAM

	PARAMETER (MAXSIN=50)
	PARAMETER (MAXSPH=300)

	INTEGER          :: NCX(MAXSIN),NCY(MAXSIN),NCZ(MAXSIN)
	REAL             :: CX(MAXSIN),ZZ(MAXSPH),YY(MAXSPH),XX(MAXSPH)
	REAL             :: FL2(MAXSPH),FL22(MAXSPH)
	REAL             :: RAD1(MAXSPH), RAD2(MAXSPH), DENS(MAXSPH)

        REAL                  :: A0(NX)
        REAL                  :: FWA(4),ADUM(3)
	LOGICAL               :: LDOC,GETOLD,GFLAG
        CHARACTER(LEN=MAXNAM) :: DOCNAM
        CHARACTER(LEN=4)      :: ANS
        CHARACTER(LEN=1)      :: GA,YN

        CHARACTER(LEN=1)      :: NULL = CHAR(0)

	REAL, PARAMETER       :: QUADPI = 3.1415926535897932
	REAL, PARAMETER       :: TWOPI=2*QUADPI
      
	IF (COPT == 'I' .OR. VERBOSE) WRITE(NOUT,100)
 100    FORMAT(
     &    ' .MENU: B    -- BLANK, CONSTANT DENSITY VOLUME'/
     &    '        C    -- CYLINDERS'/
     &    '        G    -- GAUSSIAN SPHERE/ELLIPSOID'/
     &    '        G1   -- 1ST ORDER GAUSSIAN SPHERE/ELLIPSOID (0...1)'/
     &    '        G2   -- 2ND ORDER GAUSSIAN SPHERE/ELLIPSOID (0...1)'/
     &    '        G3   -- 3RD ORDER GAUSSIAN SPHERE/ELLIPSOID (0...1)'/
     &    '        H    -- HELIX OF SPHERES'/
     &    '        HA   -- HELIX OF SPHERES, ADD DENSITIES'/
     &    '        NUM  -- LINE NUMBERS'/
     &    '        R    -- RANDOM DENSITY PATTERN'/
     &    '        T    -- TWO 3D SINE WAVES'/
     &    '        S    -- 3D SINE WAVES'/
     &    '        SP   -- SPHERES'/
     &    '        SPA  -- SPHERES, ADD DENSITIES'/
     &    '        SPV  -- SPHERES, VARIABLE DENSITIES'/
     &    '        W    -- DENSITY WEDGE'/)

       CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &      'B/C/G/G1/G2/G3/H/HA/NUM/R/T/S/SP/SPA/SPV/W', NULL,IRTFLG)

      SELECT CASE (ANS(:NCHAR))

      CASE ('T','S')
C     TEST & SINE WAVES ********************************** TEST & SINE

      IF (ANS(1:1)  ==  'T') THEN
  200    NS     = 1
         CX(1)  = 1.0
         NCX(1) = 2
         NCY(1) = 2
         NCZ(1) = 2
      ELSE

         NS = 3
         CALL RDPRI1S(NS,IDUM,'NUMBER OF SINE WAVES',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         IF (NS > MAXSIN) THEN
            WRITE(NOUT,400) MAXSIN
  400       FORMAT(' *** RESTRICTED TO ',I3,' SINE WAVES')
            NS = MAXSIN
         ENDIF

         AV = 0.0
         DO I = 1,NS
            CALL RDPRA('AMPLITUDE, SPATIAL FREQ. IN (X,Y,Z)',
     &                 4,0,.FALSE.,FWA,NVAL,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            CX(I)  = FWA(1)
            NCX(I) = FWA(2)
            NCY(I) = FWA(3)
            NCZ(I) = FWA(4)
            AV     = AV + CX(I)
         ENDDO

         DO I = 1,NS
            CX(I) = CX(I) * 2.0 / AV
            WRITE(NOUT, 900)I,CX(I),NCX(I),NCY(I),NCZ(I)
  900       FORMAT(1X,I5,F5.2,3I6)
         ENDDO
      ENDIF

      DO L = 1,NZ
         PHASEZ = FLOAT(L-1) * TWOPI / FLOAT(NZ) 
                
         DO I=1,NY
           PHASE = FLOAT(I-1) * TWOPI / FLOAT(NY)

           A0    = 0.
           DO J=1,NX
             DO K = 1,NS
             A0(J) = CX(K) * SIN(FLOAT(J-1) * TWOPI * FLOAT(NCX(K)) /
     &                FLOAT(NX) +
     &            PHASE * FLOAT(NCY(K)) + PHASEZ * FLOAT(NCZ(K)))+ A0(J)
             ENDDO
           ENDDO
           CALL WRTLIN(LUN1,A0,NX,I+(L-1)*NY)
        ENDDO
      ENDDO

C     WEDGE **************************************************** WEDGE
      CASE ('W')

      SCF = 2.0 / FLOAT(NY+NX+NZ)

      DO  L=1,NZ
         FI = FLOAT(L) * SCF

         DO I=1,NY
            FI = FI + SCF * FLOAT(I)

            DO J = 1,NX
              A0(J) = FI + SCF * FLOAT(J)
 	    ENDDO
            CALL WRTLIN(LUN1,A0,NX,(L-1)*NY+I)
 	 ENDDO
      ENDDO

C     SPHERES ************************************************ SPHERES

      CASE ('SP','SPA','SPV')
        GETOLD = .FALSE.
        CALL RDPRMC(YN,NCHAR,.TRUE.,
     &        'GET COORDINATES FROM DOCUMENT FILE (Y/N)',
     &        NULL,IRTFLG)
        IF (IRTFLG == -1) RETURN

        LDOC = .FALSE.
        XOFF = 0.
        YOFF = 0.
        ZOFF = 0.

        IF (YN .NE. 'N') THEN
C          USE DOC FILE FOR COORDINATES
           LDOC  = .TRUE.
           NOPEN = 0
           CALL FILERD(DOCNAM,NLET,NULL,'DOCUMENT',IRTFLG)
           IF (IRTFLG  ==  -1) RETURN

           NNSPH = 5
           IRET  = 1
           CALL RDPRM2S(NNSPH,IREG,NOT_USED,
     &       'NUMBER OF SPHERES, STARTING COL. FOR X,Y,Z',
     &        IRTFLG)

           CALL RDPRM3S(XOFF,YOFF,ZOFF,NOT_USED,
     &                  'X, Y, & Z OFFSETS',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
	ENDIF


        GETOLD = .FALSE.
        IF (ANS(1:3)  ==  'SPV') THEN
C          VARIABLE SPHERE DENSITIES READ FROM DOC FILE
           IDEN = 1
           IF (ANS(4:4)  ==  'A') GETOLD = .TRUE.
        ELSE
C          ALL SPHERES WILL HAVE SAME DENSITIES
           IDEN = 0

c          DEFAULT VALUE 2.0 
           AMQ = 2.0         
           CALL RDPRM1S(AMQ, NOT_USED, 
     &        'DENSITY INSIDE SPHERES (<CR> = 2.0)',
     &        IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (AMQ  ==  0.0) AMQ = 2.0     
        ENDIF

        NSPH = 0
C       THE DEFAULT LOCATION OF THE ALL THE SPHERES ARE CENTER.

        DO  L = 1, MAXSPH
           RAD1(L) = 0.0
           RAD2(L) = 0.0
           XX(L)   = NX/2 + 1
           YY(L)   = NY/2 + 1
           ZZ(L)   = NZ/2 + 1
        ENDDO
		
        ROUT = MIN(NX,NY,NZ) / 2 - 5
        RIN  = 0.0
 	CALL RDPRM2S(ROUT,RIN,NOT_USED,
     &       'OUTER AND INNER RADII OF SPHERES',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

C       CHANGED 10/12/88 MR
C       FIND HOW MANY RADII HAVE TO BE READ FROM DOCUMENT FILE:
	NSPH = 1
        IRAD = 0
        IF (ROUT < 0) IRAD = 1
        IF (RIN  < 0) IRAD = 2

	IF (.NOT. LDOC) 
     &      WRITE(NOUT,*) ' TYPE: 0, 0, 0  TO STOP INPUT'

2200	RAD2(NSPH)  = ROUT
	RAD1(NSPH)  = RIN

C       'SPA' OPTION -- ADD MASSES WHERE SPHERES OVERLAP
C       'SP'  OPTION -- LET SPHERES INTERPENETRATE

	SWITCH = 0.0
	IF (ANS(3:3) == 'A') SWITCH = 1.0

	IF (LDOC) THEN
           CALL UNSAV(DOCNAM,NOPEN,NDOC,NSPH,CX,
     &                IREG+2+IRAD+IDEN,
     &                LERR,1)  
c          SEQUENTIAL READ SPECIFIED FOR DOC FILE
           NOPEN   = 1
           ADUM(1) = CX(IREG)
           ADUM(2) = CX(IREG+1)
           ADUM(3) = CX(IREG+2)
           IF (IRAD > 0) RAD2(NSPH) = CX(IREG+3)
           IF (IRAD > 1) RAD1(NSPH) = CX(IREG+4)
           IF (IDEN > 0) DENS(NSPH) = CX(IREG+3+IRAD)
           WRITE(NOUT,2400)NSPH,(ADUM(K),K=1,3)
	ELSE
	   ADUM(1) = NX/2 +1
	   ADUM(2) = NY/2 +1
	   ADUM(3) = NZ/2 +1
           CALL RDPRM3S(ADUM(1),ADUM(2),ADUM(3),NOT_USED,
     &       'CENTER COORDINATES X,Y,Z (or <CR> FOR IMAGE CENTER)',
     &       IRTFLG)

C*         ALTERED SO DOES NOT ACCEPT LAST 0.0,0.0,0,0 al june 88
           IF (ADUM(1) == 0. .AND. 
     &         ADUM(2) == 0. .AND.
     &         ADUM(3) == 0.) THEN
               NSPH = NSPH - 1
               IF (NSPH  ==  0) RETURN
               GOTO 3000
           ENDIF

           WRITE(NOUT,2400) NSPH,(ADUM(K),K=1,3)
 2400      FORMAT(I4,3F12.6)
	ENDIF

	XX(NSPH) = ADUM(1) + XOFF
	YY(NSPH) = ADUM(2) + YOFF
	ZZ(NSPH) = ADUM(3) + ZOFF
 
        IF (ZZ(NSPH)+RAD2(NSPH)  >  NZ .OR.
     &      ZZ(NSPH)-RAD2(NSPH)  <   0 .OR.
     &      YY(NSPH)+RAD2(NSPH)  >  NY .OR.
     &      YY(NSPH)-RAD2(NSPH)  <   0 .OR.
     &      XX(NSPH)+RAD2(NSPH)  >  NX .OR.
     &      XX(NSPH)-RAD2(NSPH)  <   0) WRITE(NOUT,2450) NSPH
 2450   FORMAT(' *** SPHERE ',I3,' WILL BE TRUNCATED')

	IF (LDOC .AND. NSPH  ==  NNSPH) GOTO 3000  
     
	IF (NSPH  >= MAXSPH) THEN
	   WRITE(NOUT,2500)MAXSPH
 2500	   FORMAT(' ** NUMBER OF SPHERES LIMITED TO',I3)
	   GOTO 3000
	ENDIF

	NSPH = NSPH + 1
	GOTO 2200


      CASE ('HA','H')
C       SPHERES ARRANGED ON HELIX ****************************** HELIX

        GETOLD = .FALSE.
        IDEN   = 0
        AMQ    = 2.0
        SWITCH = 0.0
        CALL RDPRM1S(AMQ, NOT_USED, 
     &      'DENSITY VALUE FOR SPHERES (DEFAULT = 2.0)',
     &       IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        RADT = 3
        RHEL = 3
	CALL RDPRM2S(RADT,RHEL,NOT_USED,
     &            'SPHERE RADIUS, HELIX RADIUS',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        NSPH  = 100
        NTURN = 5
	CALL RDPRI2S(NSPH,NTURN,NOT_USED,
     &               'NO. OF SPHERES, NO. OF TURNS',IRTFLG)
	IF (NSPH  >=  MAXSPH) THEN
	   CALL ERRT(102,'NUMBER OF SPHERES LIMITED TO',MAXSPH)
	   RETURN
	ENDIF

	DO N=1,NSPH
           RAD1(N) = 0.0
           RAD2(N) = RADT
           XX(N)   = NX/2+1    + RHEL * COS(TWOPI * FLOAT(N-1) * 
     &               FLOAT(NTURN) / FLOAT(NSPH))
           ZZ(N)   = NZ/2+1  + RHEL * SIN(TWOPI * FLOAT(N-1) *
     &               FLOAT(NTURN) / FLOAT(NSPH))
           YY(N)   = FLOAT(N-1)  * FLOAT(NY) / FLOAT(NSPH)
        ENDDO

 3000 IF (NSPH  ==  0) NSPH = 1
      DO L1 = 1, NZ
         DO  N = 1, NSPH
           FL2(N) = (L1-ZZ(N))**2
         ENDDO

         DO L2 = 1, NY
           DO N  =1, NSPH
              FL22(N) = FL2(N) + (L2-YY(N))**2
           ENDDO

           IF (GETOLD) THEN
C             WANT TO ADD NEW SPERES TO EXISTING FILE
              CALL REDLIN(LUN1,A0,NX,(L1-1)*NY+L2)
           ENDIF

           DO L3 = 1, NX
             IF (.NOT. GETOLD) A0(L3) = 0.0

             DO N=1, NSPH
               RAD41 = RAD1(N)**2
               RAD44 = RAD2(N)**2
               IF (IDEN  >  0) AMQ = DENS(N)
               XYZ   = FL22(N) + (L3-XX(N))**2
               IF (XYZ  >= RAD41 .AND.  XYZ  <=  RAD44) THEN
C                NAIK GENERAL AMQ NOT 2
                 A0(L3) = AMQ + A0(L3) * SWITCH    
               ENDIF
             ENDDO
           ENDDO
           CALL WRTLIN(LUN1,A0,NX,(L1-1)*NY+L2)
         ENDDO
      ENDDO

      CASE ('R','RAN')
C     RANDOM DISTRIBUTION ************************************* RANDOM

      CALL RDPRMC(GA,NCHAR,.TRUE.,
     &            'GAUSSIAN DISTRIBUTED? (Y/N)',
     &            NULL,IRTFLG)
      IF (IRTFLG  ==  -1) RETURN
      
      GFLAG = (GA  ==  'Y')
      IF (GFLAG) THEN
	CALL RDPRM2s(XM,SD,NOT_USED,
     &     'MEAN AND STANDARD DEVIATION OF GAUSSIAN DIST.',
     &     IRTFLG)
      ENDIF

      DO L = 1,NZ
        DO I = 1,NY
          IF (GFLAG) THEN
             DO K = 1,NX
                A0(K) = RANN(XM,SD)
             ENDDO
          ELSE
	     CALL  RANDOM_NUMBER(HARVEST=A0)
          ENDIF
          CALL WRTLIN(LUN1,A0,NX,I+(L-1)*NY)
        ENDDO
      ENDDO


C     BLANK (CONSTANT DENSITY) ********************************* BLANK
      CASE ('B','BL')

      BVAL = 0.0
      CALL RDPRM1S(BVAL,NOT_USED,'DENSITY',IRTFLG)
      IF (IRTFLG  .NE. 0) RETURN

C     CREATE BACKGROUND VOLUME
      A0 = BVAL

      DO L=1,NZ                                  
         DO I=1,NY                                 
            CALL WRTLIN(LUN1,A0,NX,I+(L-1)*NY)   
         ENDDO
      ENDDO



      CASE ('C','CYL')
C     CYLINDERS (INTERPENETRATING) ************************* CYLINDERS

C     CLEAR VOLUME FIRST (SET VOLUME BACKGROUND)

      BVAL = 0.0
      CALL RDPRM1S(BVAL,NOT_USED,
     &             'BACKGROUND DENSITY',IRTFLG)
      IF (IRTFLG .NE.0) RETURN

C     CREATE BACKGROUND VOLUME
      A0 = BVAL
      DO L=1,NZ                                  
         DO I=1,NY                                 
            CALL WRTLIN(LUN1,A0,NX,I+(L-1)*NY)   
         ENDDO
      ENDDO

 5000 CALL RDPRMC(YN,NA,.TRUE.,
     &   'CYLINDER AXIS; X, Y, Z (or Q TO END CYLINDER ENTRY)',
     &      NULL,IRTFLG)
      IF (YN .NE. 'X' .AND. 
     &    YN .NE. 'Y' .AND. 
     &    YN .NE. 'Z') RETURN 

      RADC = 5
      HT   = ((MIN(NX,NY,NZ))/2) - 3
      CALL RDPRM2S(RADC,HT,NOT_USED,'RADIUS, HEIGHT',IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      IF (RADC == 0.0 .OR.  HT == 0.0) RETURN

      HT = HT / 2.0

      XH = (NX/2) + 1
      YH = (NY/2) + 1
      CALL RDPRM2S(XH,YH,NOT_USED,
     &     'CENTER COORDINATES X,Y (or <CR> FOR IMAGE CENTER)',
     &     IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      ZH    = (NZ/2) + 1
      VALUE = 1.0
      CALL RDPRM2S(ZH,VALUE,NOT_USED,
     &             'CENTER COORDINATE Z, DENSITY',
     &             IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     ****************** CYLINDER AXIS ALONG Y ***********

      IF (YN == 'Y') THEN 
         NSL1  = ZH - RADC - 0.5
         NSL2  = ZH + RADC + 0.5
         IF (NSL1  <   1) NSL1 = 1
         IF (NSL2  >  NZ) NSL2 = NZ
         NX1   = XH - RADC - 0.5
         NX2   = XH + RADC + 0.5
         IF (NX1  <   1) NX1 = 1
         IF (NX2  >  NX) NX2 = NX
         NY1   = YH - HT - 0.5
         NY2   = YH + HT + 0.5
         IF (NY1  <   1) NY1 = 1
         IF (NY2  >  NY) NY2 = NY

         RAD44 = RADC**2

         DO L=NSL1,NSL2
           I1 = (L - 1) * NY
           Z2 = (FLOAT(L) - ZH) **2
           DO I=NY1,NY2
             CALL REDLIN(LUN1,A0,NX,I+I1)
             DO K=NX1,NX2
                IF ((FLOAT(K)-XH)**2 + Z2 <= RAD44) A0(K) = VALUE                 
             ENDDO
             CALL WRTLIN(LUN1,A0,NX,I1+I)
           ENDDO
         ENDDO
      ENDIF
 
C     ****************** CYLINDER AXIS ALONG Z ***********
      IF (YN  ==  'Z') THEN 

         NSL1  = ZH - HT - 0.5
         NSL2  = ZH + HT + 0.5
         IF (NSL1 < 1)      NSL1 = 1
         IF (NSL2 > NZ) NSL2 = NZ
         NX1   = XH - RADC - 0.5
         NX2   = XH + RADC + 0.5
         IF (NX1 <  1)    NX1=  1
         IF (NX2 >  NX) NX2 = NX
         NY1   = YH - RADC - 0.5
         NY2   = YH + RADC + 0.5
         IF( NY1 <  1)    NY1 = 1
         IF (NY2 >  NY) NY2 = NY

         RAD44 = RADC**2

         DO L=NSL1,NSL2
           I1 = (L-1)*NY

           DO I=NY1,NY2
             CALL REDLIN(LUN1,A0,NX,I+I1)
             Y2 = (FLOAT(I)-YH)**2

             DO  K=NX1,NX2
                IF ((FLOAT(K)-XH)**2+Y2 <= RAD44) A0(K) = VALUE 
             ENDDO                
             CALL WRTLIN(LUN1,A0,NX,I1+I)
           ENDDO
         ENDDO
      ENDIF

C     ****************** CYLINDER AXIS ALONG X ***********
      IF (YN  ==  'X') THEN 
         NSL1  = ZH - RADC - 0.5
         NSL2  = ZH + RADC + 0.5
         IF (NSL1  <   1) NSL1=1
         IF (NSL2  >  NZ) NSL2=NZ
         NX1   = XH - HT - 0.5
         NX2   = XH + HT + 0.5
         IF (NX1  <   1) NX1 = 1
         IF (NX2  >  NX) NX2 = NX
         NY1   = YH - RADC - 0.5
         NY2   = YH + RADC + 0.5
         IF (NY1 <  1) NY1=1
         IF (NY2 > NY) NY2=NY
         RAD44 = RADC**2

         DO L=NSL1,NSL2
            I1 = (L-1) * NY
            Z2 = (FLOAT(L)-ZH)**2

            DO I=NY1,NY2
               Y2 = (FLOAT(I)-YH)**2

               IF (Z2+Y2  <=  RAD44) THEN
                  CALL REDLIN(LUN1,A0,NX,I+I1)
                  DO  K=NX1,NX2
                     A0(K) = VALUE                 
                  ENDDO
                CALL WRTLIN(LUN1,A0,NX,I1+I)
              ENDIF
           ENDDO
         ENDDO
      ENDIF
      GOTO 5000

      CASE ('NUM','N')
C     PUT NUMBERS INTO VOLUME:********************************* NUMBERS

      IPOS = 1
      CALL RDPRI1S(IPOS,NOT_USED,'POSITION IN LINE',IRTFLG)

      SCALE = 1
      CALL RDPRM1S(SCALE,NOT_USED,'SCALING FACTOR',IRTFLG)

      SUM = 0.0
      DO L=1,NZ
        DO I=1,NY
          A0       = 0.0
          A0(IPOS) = FLOAT((L-1) * NY +I ) * SCALE
          SUM      = SUM + A0(IPOS)

          CALL WRTLIN(LUN1,A0,NX,I+(L-1)*NY)
        ENDDO
      ENDDO

      CALL RDPRMC(ANS,NCHAR,.TRUE.,
     &            'NUMBER PIXELS IN LINE? (Y/N)',
     &             NULL,IRTFLG)

      IF (ANS(1:1)  ==  'Y') THEN
        CALL RDPRI2S(ISLIC,ILINE,NOT_USED,
     &              'SLICE NUMBER, LINE NUMBER',IRTFLG)
        if (irtflg .ne. 0) then
           write(6,*) 'Non-zero irtflg:',irtflg,ISLIC,ILINE
           return
        endif

        CALL RDPRM1S(SCALE2,NOT_USED,
     &               'SCALING FACTOR',IRTFLG)
        if (irtflg .ne. 0) then
           write(6,*) 'Non-zero irtflg:',irtflg,scale2
           return
        endif

        ILL = (ISLIC-1) * NY + ILINE
        CALL REDLIN(LUN1,A0,NX,ILL)
        DO  I=1,NX
           A0(I) = FLOAT(I) * SCALE2
           SUM   = SUM + A0(I)
        ENDDO

        SUM = SUM - FLOAT(IPOS)*SCALE
        CALL WRTLIN(LUN1,A0,NX,ILL)
      ENDIF

      CASE ('G','G1','G2','G3')
C     ******************************************* 3D GAUSSIAN FUNCTION

C      CREATES A 3D GAUSSIAN DENSITY DISTRIBUTION NORMALIZED 
C      SUCH THAT THE SUM OF ALL THE VOXEL DENSITIES EQUAL 1.0.

        XCEN = (NX/2) + 1
        YCEN = (NY/2) + 1
        ZCEN = (NZ/2) + 1
        CALL RDPRM3S(XCEN,YCEN,ZCEN,NOT_USED,
     &    'CENTER COORDINATES X,Y,Z (or <CR> FOR IMAGE CENTER)',
     &    IRTFLG)
        if (irtflg .ne. 0) then
           write(6,*) 'Non-zero irtflg:',irtflg,xcen
        endif

        STDX = (NX/2) - 1
        STDY = (NY/2) - 1
        STDZ = (NZ/2) - 1
        CALL RDPRM3S(STDX,STDY,STDZ,NOT_USED,
     &              'RADII IN  X,Y,Z ( = STD. DEV.)',IRTFLG)

        WRITE(NOUT,90) XCEN,YCEN,ZCEN, STDX,STDY,STDZ
90      FORMAT('  CENTER:',G8.2,',',G8.2,',',G8.2,
     &         '  RADII:', G8.2,',',G8.2,',',G8.2) 

C       SET THE ORDER FOR SUPERGAUSSIAN 
        NORDER = 1
        IF (ANS(2:2) == '2') NORDER = 2
        IF (ANS(2:2) == '3') NORDER = 3

	GNM = 1.0 / STDX / STDY / STDZ / (2*QUADPI)**1.5
        IF (NCHAR > 1 ) GNM = 1.0    ! G2....

	TNM = ALOG(1.0 / TINY(GNM))

        STDXSQ = STDX**2
        STDYSQ = STDY**2
        STDZSQ = STDZ**2

        DO K = 1,NZ
           DO J = 1,NY
              DO I = 1,NX
	        EEE = 0.5 * ((I - XCEN) **2 / STDXSQ +
     &		             (J - YCEN) **2 / STDYSQ +
     &		             (K - ZCEN) **2 / STDZSQ)
	        IF (EEE  >= TNM) THEN
	           A0(I) = 0.0
	        ELSE  
	           EEE   = 0.5 * (2*EEE)**NORDER
                   A0(I) = GNM * EXP(-EEE)
	        ENDIF
             ENDDO

             CALL WRTLIN(LUN1,A0,NX,J+(K-1)*NY)
           ENDDO
        ENDDO

      CASE DEFAULT
         CALL ERRT(101,'UNIDENTIFIED OPTION',NE)

      END SELECT

      END

