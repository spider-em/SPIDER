C**************************************************************************
C
C TRAFSIM.F                   BIMAL RATH Initial revision      9/26/05 
C
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
C TRAFSIM.F
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

         SUBROUTINE TRAFSIM(LUN1)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC' 
 
         CHARACTER(LEN=MAXNAM)   ::  FILNAM
         COMMON /COMMUN/FILNAM
         REAL, ALLOCATABLE, DIMENSION(:,:) :: ATOMLOC
         REAL, ALLOCATABLE, DIMENSION(:,:) :: QB
         
         REAL          LAMBDA,KAPPA,KM       
         CHARACTER    NULL 
        
         CHARACTER(LEN=MAXNAM):: PDBFILE, RECLIN
         CHARACTER(LEN=3)::      ATOM*3
         CHARACTER               TYPE_IMAGE, CHANGE_INT  

         PARAMETER (QUADPI = 3.141592653589793238462643383279502884197)
         PARAMETER (KPP = 11.09035491943359380000)
         PARAMETER (FP = -9.8696050643920898400000)
         
C       TYPES OF ATOMS SUPPORTED. TO ADD MORE ATOMS, CHANGE VALUE OF
C       NTYPES_ATOM TO DESIRED NUMBER AND ENTER VALUES FOR BIG_A AND
C       BIG_B FROM JOACHIM HAASE'S 1970 PAPER IN ZEITSCHRIFT FUR 
C       NATURFORSCHUNG,PAGE 1219-1235 (THOMAS-FERMI-DIRAC-FALL 
C       POTENTIAL) ADD A LINE e.g.: ELSEIF (ATOM(2:2) .EQ. 'U') THEN
C       VATOM=92.0 REFLECTING THE ATOMIC SYMBOL AND ATOMIC NUMBER OF
C       THE ADDED ATOMS AT THE APPROPRIATE PLACE IN THE CODES BELOW.
                 
         PARAMETER (NTYPES_ATOM = 9)


C       EACH BIG_A COLUMN HAS: ATOMIC #, A0,0; A0,1; A0,2; A1,0; A1,1; A1,2;
C       A2,0; A2,1; A2,2;A3,0; A3,1; A3,2; A4,0; A4,1; A4,2
 
         DOUBLE PRECISION, DIMENSION(16,NTYPES_ATOM) :: BIG_A =
     &    RESHAPE(( /6.0, 0.1009E+01, 0.1991E-02, -0.1725E-05, 
     &     -.3722, 0.1129E-04, -0.4992E-07, 
     &     0.1223E-01,-0.5984E-06, 0.2483E-08,
     &     -0.2109E-03,0.8397E-08, -0.3356E-10, 
     &     0.1381E-05, -0.7332E-10, 0.3155E-12,
     &     7.0, 0.1123E+01, 0.2012E-02, -0.1819E-05, 
     &     -0.3687, 0.1524E-04, -0.6472E-07,
     &     0.1209E-01, -0.8274E-06, 0.3336E-08, 
     &     -0.2084E-03, 0.1239E-07, -0.4726E-10, 
     &     0.1365E-05, -0.1032E-09, 0.4138E-12, 
     &     8.0, 0.1219E+01, 0.2031E-02, -0.1892E-05,
     &     -0.3656, 0.2102E-04, -0.9057E-07,
     &     0.1197E-01, -0.1176E-05, 0.4927E-08,
     &     -0.2062E-03, 0.1906E-07, -0.7797E-10,
     &     0.1351E-05, -0.1534E-10, 0.6454E-12,
     &     15.0, 0.1646E+01, 0.2222E-02, -0.2692E-05,
     &     -0.3511, 0.7217E-04, -0.3063E-06,
     &     0.1144E-01,-0.4186E-05, 0.1783E-07,
     &     -0.1965E-03, 0.7491E-07, -0.3154E-09,
     &     0.1287E-05, -0.5395E-09, 0.2270E-11, 
     &     16.0, 0.1687E+01, 0.2253E-02, -0.2822E-05,
     &     -0.3498, 0.8137E-04, -0.3451E-06,
     &     0.1140E-01, -0.4721E-05, 0.2013E-07,
     &     -0.1957E-03, 0.8504E-07, -0.3592E-09,
     &     0.1282E-05, -0.6070E-09, 0.2560E-11,
     &     20.0, 0.1825E+01, 0.2379E-02, -0.3349E-05,
     &     -0.3459, 0.1231E-03, -0.5174E-06,
     &     0.1128E-01, -0.7112E-05, 0.3007E-07,
     &     -0.1938E-03, 0.1308E-06, -0.5485E-09,
     &     0.1269E-05, -0.9123E-09, 0.3804E-11,            
     &     26.0, 0.1977E+01, 0.2573E-02, -0.4135E-05,
     &     -0.3429, 0.1925E-03, -0.8038E-06,
     &     0.1124E-01, -0.1106E-04, 0.4653E-07,
     &     -0.1933E-03, 0.2090E-06, -0.8773E-09,
     &     0.1266E-05, -0.1434E-08, 0.5991E-11,          
     &     53.0, 0.2337E+01, 0.3340E-02, -0.7153E-05,
     &     -0.3458, 0.4716E-03, -0.1868E-05,
     &     0.1176E-01, -0.2742E-04, 0.1105E-06,
     &     -0.2079E-03, 0.5787E-06, -0.2374E-08,
     &     0.1370E-05, -0.4159E-08, 0.1727E-10,     
     &     92.0, 0.2560E+01, 0.4140E-02, -0.1030E-04,
     &     -0.3566, 0.6479E-03, -0.2367E-05,
     &     0.1252E-01, -0.3584E-04, 0.1301E-06,
     &     -0.2224E-03, 0.7265E-06, -0.2614E-08,
     &     0.1484E-05, -0.5299E-08, 0.1924E-10/),  
     &     (/16,NTYPES_ATOM/)) 

C       EACH BIG_B COLUMN HAS: ATOMIC #, B0,0; B0,1; B0,2; B1,0; B1,1; B1,2; 
C       B2,0; B2,1; B2,2; B3,0; B3,1; B3,2; B4,0; B4,1; B4,2

         DOUBLE PRECISION, DIMENSION(16,NTYPES_ATOM) :: BIG_B =
     &    RESHAPE(( /6.0, 0.1931E-01, -0.1632E-03, 0.6395E-06,
     &     0.3922E-01, -0.3323E-03, 0.1304E-05,
     &     -0.8743E-03, 0.6744E-05, -0.2599E-07,
     &     0.1236E-04, -0.1041E-06, 0.4072E-09,
     &     -0.6837E-07, 0.5822E-09, -0.2286E-11,
     &     7.0, 0.2274E-01, -0.1912E-03, 0.7485E-06,
     &     0.4436E-01, -0.3751E-03, 0.1471E-05,
     &     -0.9709E-03, 0.7561E-05, -0.2927E-07,
     &     0.1344E-04, -0.1125E-06, 0.4385E-09,
     &     -0.7846E-07, 0.6219E-09, -0.2433E-11,
     &     8.0, 0.2607E-01, -0.2170E-03, 0.8462E-06, 
     &     0.4931E-01, -0.4165E-03, 0.1632E-06,
     &     -0.1061E-02, -0.8327E-05, -0.3235E-07,
     &     0.1441E-04, -0.1204E-06, 0.4692E-09,
     &     -0.7793E-07, 0.6601E-09, -0.2584E-11,
     &     15.0, 0.4786E-01, -0.3707E-03, 0.1414E-05,
     &     0.7954E-01, -0.6645E-03, 0.2596E-05,
     &     -0.1537E-02, 0.1223E-04, -0.4753E-07,
     &     0.1860E-04, -0.1504E-06, 0.5784E-09,
     &     -0.9274E-07, 0.7628E-09, -0.2945E-11, 
     &     16.0, 0.5068E-01, -0.3876E-03, 0.1472E-05,
     &     0.8335E-01, -0.6952E-03, 0.2715E-05,
     &     -0.1589E-02, 0.1264E-04, -0.4912E-07,
     &     0.1893E-04, -0.1522E-06, 0.5843E-09,
     &     -0.9320E-07, 0.7612E-09, -0.2931E-11,
     &     20.0, 0.6099E-01, -0.4417E-03, 0.1650E-05,
     &     0.9767E-01, -0.8088E-03, 0.3148E-05,
     &     -0.1771E-02, 0.1401E-04, -0.5417E-07,
     &     0.1968E-04, -0.1554E-06, 0.5894E-09,
     &     -0.9260E-07, 0.7276E-09, -0.2748E-11,        
     &     26.0, 0.7356E-01, -0.4791E-03, 0.1724E-05,
     &     0.1167, -0.9532E-03, 0.3701E-05,
     &     -0.1985E-02, 0.1558E-04, -0.6005E-07,
     &     0.2041E-04, -0.1544E-06, 0.5797E-09,
     &     -0.8631E-07, 0.6352E-09, -0.2333E-11,           
     &     53.0, 0.9503E-01, -0.2194E-03, 0.3460E-06,
     &     0.1796, -0.1396E-02, 0.5369E-05,
     &     -0.2645E-02, 0.2227E-04, -0.8921E-07,
     &     0.2211E-04, -0.2145E-06, 0.9156E-09,
     &     -0.5834E-07, 0.6844E-09, -0.3190E-11,      
     &     92.0, 0.8521E-01, 0.3492E-03, -0.1949E-05,
     &     0.2360, -0.1733E-02, 0.6579E-05,
     &     -0.3231E-02, 0.2956E-04, -0.1190E-06,
     &     0.2778E-04, -0.3886E-06, 0.1734E-08,
     &     -0.9099E-07, 0.2260E-08, -0.1117E-10 /),
     &     (/16,NTYPES_ATOM/)) 


C       EACH SMALL_A COLUMN HAS: ATOMIC  #, a0, a1, a2, a3, a4        
C       EACH SMALL_B COLUMN HAS: ATOMIC  #, b0, b1, b2, b3, b4        

         DOUBLE PRECISION, DIMENSION(6,NTYPES_ATOM) :: SMALL_A
         DOUBLE PRECISION, DIMENSION(6,NTYPES_ATOM) :: SMALL_B

         CALL FILERD(FILNAM,NLET,NULL,'OUTPUT',IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

C       MICROSCOPE OPERATING VOLTAGE IN KV       
         CALL RDPRM(VOLTK,NOT_USED,'OPERATING VOLTAGE [KV]')             
         
         DO I = 1,NTYPES_ATOM
            DO J = 1,6
               IF (J .EQ. 1) THEN
                  SMALL_A(J,I) = BIG_A(J,I)
               ELSE 
                  K = (J-1)*3 - 1
                  SMALL_A(J,I) = BIG_A(K,I) + BIG_A(K+1,I) * 
     &                         VOLTK + BIG_A(K+2,I) * VOLTK * VOLTK
     
               ENDIF      
            ENDDO
          
         ENDDO

         DO I = 1,NTYPES_ATOM
            DO J = 1,6
               IF (J .EQ. 1) THEN
                  SMALL_B(J,I) = BIG_B(J,I)
               ELSE 
                  K = (J-1)*3 - 1
                  SMALL_B(J,I) = BIG_B(K,I) + BIG_B(K+1,I) * 
     &                        VOLTK + BIG_B(K+2,I) * VOLTK * VOLTK
     
               ENDIF      
            ENDDO
          
         ENDDO   

C       OPEN PDB FILE (WHICH IS READABLE ASCII) AS FORMATTED SEQ.    
         
         LENREC = 0
         NULL = CHAR(0)
         IRTFLG = 0 
         NATOM = 0       
         XMAX  = -1.0E23
         XMIN  = 1.0E23
         YMAX  = -1.0E23
         YMIN  = 1.0E23
         ZMAX  = -1.0E23
         ZMIN  = 1.0E23
         UX    = 0.0
         UY    = 0.0
         UZ    = 0.0     
         
         CALL OPAUXFILE(.TRUE.,PDBFILE,NULL,LUN1,LENREC,'O',
     &                 'PDB INPUT',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
         
C       FIND THE NUMBER OF ATOMS IN THE PDB FILE. HELPS IN ALLOCATING MEMORY 
C       FOR ARRAY ATOMLOC
15       READ(LUN1,20) RECLIN
20       FORMAT(A80)

         IF (RECLIN(1:4) .EQ. 'ATOM') THEN
            READ (RECLIN,50) HEAD,ATOM,X,Y,Z,OCCUPANCY,TEMPERATURE
50          FORMAT(A6,5X,1X,A3,2X,3X,1X,1X,4X,1X,3X,3F8.3,2F6.2)        
C          WRITE(NOUT,*)'ATOM #',NATOM,'  ',ATOM,' AT ',X,Y,Z,
C     &                    OCCUPANCY,TEMPERATURE          
            NATOM=NATOM+1
            
C         COORDINATE SYSTEM IS DIFFERENT IN O AND SPIDER 
            W = X
            X = Y
            Y = W
            Z = -Z

            XMAX=AMAX1(XMAX,X)
            XMIN=AMIN1(XMIN,X)
            YMAX=AMAX1(YMAX,Y)
            YMIN=AMIN1(YMIN,Y)
            ZMAX=AMAX1(ZMAX,Z)
            ZMIN=AMIN1(ZMIN,Z)
            UX=UX+X
            UY=UY+Y
            UZ=UZ+Z
               
            GOTO 15
         ELSEIF (RECLIN(1:3) .NE. 'END') THEN
            GOTO 15
         ENDIF


         IF (NATOM .LE. 0)  THEN
            CALL ERRT(101,'NO ATOMS IN THE PDB FILE',NE)
            GOTO 9999
         ENDIF

C       ATOMLOC HAS 5 COLUMNS. THEY HOLD, ATOM NUMBER (1 - NATOM), ATOMIC NUMBER,
C       X, Y, Z POSITIONS STARTING FROM FIRST COLUMN.
         ALLOCATE (ATOMLOC(NATOM,5), STAT=IRTFLG)
         IF (IRTFLG .NE. 0) THEN 
            CALL ERRT(46,'TRAFSIM, ATOMLOC',IER)
            GOTO 9999
         ENDIF


C       SET CELL DIMENSION

        SX=(XMAX-XMIN)+3.0
        SY=(YMAX-YMIN)+3.0
        SZ=(ZMAX-ZMIN)+3.0
        WRITE(NOUT,*)' CELL SIZE: ',SX,SY,SZ

C       INPUT SAMPLE SIZE
       
        CALL RDPRM(PIXEL,NOT_USED,'VOXEL SIZE [A]')

        UX=UX/NATOM
        UY=UY/NATOM
        UZ=UZ/NATOM
          
        NX = NINT(SX/PIXEL)+1
        NY = NINT(SY/PIXEL)+1
        NZ = NINT(SZ/PIXEL)+1
        WRITE(NOUT,506) NX,NY,NZ
506     FORMAT(' Minimum size of the volume: ',3I6)

        CALL RDPRI3S(NX,NY,NZ,NOT_USED,'VOLUME NX, NX & NZ',IRTFLG)

        I=NX/2+1
        XOFFSET=PIXEL*I
        I=NY/2+1
        YOFFSET=PIXEL*I
        I=NZ/2+1
        ZOFFSET=PIXEL*I
        
C       CHECK IF THE VOLUME IS LARGE ENOUGH TO CONTAIN THE OBJECT

        IF (XMIN-UX .LT. (-1.*(XOFFSET)) .OR. YMIN-UY .LT. 
     &    (-1.*(YOFFSET)) .OR. ZMIN-UZ .LT. 
     &                       (-1.*(ZOFFSET))) THEN

            CALL ERRT(101,'VOLUME TOO SMALL TO CONTAIN OBJECT',NE)
            GOTO 9999
        ENDIF
        
        TMPX=FLOAT(NX-1)*PIXEL-(XOFFSET)
        TMPY=FLOAT(NY-1)*PIXEL-(YOFFSET)
        TMPZ=FLOAT(NZ-1)*PIXEL-(ZOFFSET)        
        
        IF (XMAX-UX .GT. TMPX .OR. YMAX-UY .GT. TMPY
     &              .OR. ZMAX-UZ .GT. TMPZ) THEN

           CALL ERRT(101,'VOLUME TOO SMALL TO CONTAIN OBJECT.',NE)
           GOTO 9999
        ENDIF

            
         REWIND   LUN1
        
         IATOMNUM = 0
         
40       CONTINUE        
         READ(LUN1,20) RECLIN
         IF (RECLIN(1:4) .EQ. 'ATOM' ) THEN
            READ(RECLIN,50)HEAD,ATOM,X,Y,Z,OCCUPANCY,TEMPERATURE

            IATOMNUM = IATOMNUM + 1


C          COORDINATE SYSTEM IS DIFFERENT IN O AND SPIDER 
            W = X
            X = Y
            Y = W
            Z = -Z
                    
            X = X-UX
            Y = Y-UY
            Z = Z-UZ
    
            X = IFIX(X+XOFFSET)     
            Y = IFIX(Y+YOFFSET)     
            Z = IFIX(Z+ZOFFSET)
                    
            IF (ATOM(2:2) .EQ. 'H') THEN
               VATOM=1.0
            ELSEIF ((ATOM(2:2) .EQ. 'C') .AND.(ATOM(3:3) .NE. 'a')) THEN
               VATOM=6.0
            ELSEIF (ATOM(2:2) .EQ. 'N') THEN
               VATOM=7.0
            ELSEIF (ATOM(2:2) .EQ. 'O') THEN
               VATOM=8.0
            ELSEIF (ATOM(2:2) .EQ. 'S') THEN
               VATOM=16.0
            ELSEIF (ATOM(2:2) .EQ. 'P') THEN
               VATOM=15.0
            ELSEIF ((ATOM(2:2) .EQ. 'C') .AND. (ATOM(3:3) .EQ. 'a'))THEN
               VATOM=20.0                             
            ELSEIF ((ATOM(2:2) .EQ. 'F') .AND. (ATOM(3:3) .EQ. 'e'))THEN
               VATOM=26.0              
            ELSEIF (ATOM(2:2) .EQ. 'J') THEN
               VATOM=53.0              
            ELSEIF (ATOM(2:2) .EQ. 'U') THEN
               VATOM=92.0                      
            ELSE
               WRITE(NOUT, *) ' SPECIAL ATOM ENCOUNTERED IN: '
               WRITE(NOUT,*) RECLIN 
               GOTO 9999
            ENDIF       
            ATOMLOC(IATOMNUM, 1) = FLOAT(IATOMNUM)
            ATOMLOC(IATOMNUM, 2) = VATOM
            ATOMLOC(IATOMNUM, 3) = X                
            ATOMLOC(IATOMNUM, 4) = Y
            ATOMLOC(IATOMNUM, 5) = Z            
         
            GOTO 40
         ELSEIF (RECLIN(1:3) .NE. 'END') THEN
            GOTO 40
         ENDIF

        CALL RDPRM(CS,NOT_USED,'SPHERICAL ABERRATION CS [MM]')
           IF (CS < 0.0001)    CS = 0.0001

         CALL RDPRM(DZ,NOT_USED,'DEFOCUS[A]')
         LAMBDA = (12.3)/(SQRT((VOLTK*1000 + VOLTK*VOLTK)))

C       MAXIMUM SPATIAL FREQUENCY KM = 1/2.0*PIXEL SIZE
         KM = 1/(2.0 * PIXEL)

         CALL RDPRM2(Q,DS,NOT_USED,
     &        'SOURCE SIZE[A-1], DEFOCUS SPREAD [A]')

         CALL RDPRM2(DZA,AZZ,NOT_USED,'ASTIGMATISM [A], AZIMUTH [DEG]')
         
         CALL RDPRM(ENV,NOT_USED,'GAUSSIAN ENVELOPE HALFWIDTH [1/A]')
         ENV = 1./ENV**2
         
         CALL RDPRMC(TYPE_IMAGE,NCHAR,.TRUE.,
     &     'IMAGE CREATED BY  (A)AMPLITUDE ONLY, (P)PHASE ONLY, OR (B)
     &BOTH AMPLITUDE AND PHASE',NULL,IRT)
         
     
         IF (TYPE_IMAGE .NE. 'A' .AND. TYPE_IMAGE .NE. 'P'  
     &                         .AND. TYPE_IMAGE .NE. 'B' ) THEN
            WRITE(NOUT,*)'INVALID INPUT PARAMETER FOR 
     &     "IMAGE CREATION OPTION" '
            GOTO 9999
         ENDIF

         CALL RDPRMC(CHANGE_INT,NCHAR,.TRUE.,
     & 'IMAGE INTENSITY TO BE MULTIPLIED BY,  (E)ENVELOPE FUNCTION ONLY,
     &(B)BOTH ENVELOPE FUNCTION AND GAUSSIAN FUNCTION, (N)NONE'
     & ,NULL,IRT)

         
         IF (CHANGE_INT .NE. 'E' .AND. CHANGE_INT .NE. 'B'  
     &                         .AND. CHANGE_INT .NE. 'N' ) THEN
            WRITE(NOUT,*)'INVALID INPUT PARAMETER FOR 
     &                            "INTENSITY MULTIFICATION" '
            GOTO 9999
         ENDIF

         NSAM = NX
         NROW = NY


         IF (MOD(NSAM,2) .EQ. 0)  THEN
            IFORM = 1
            LSM   = NSAM+2
         ELSE
            IFORM = 1
            LSM   = NSAM+1
         ENDIF

         IYC    = NROW/2+1

         NSLICE = 1
         MAXIM  = 0
         IRTFLG = 0
         LUN = 11

         NXX = NX
         NYY = NY       
         CALL OPFILEC(0,.FALSE.,FILNAM,LUN,'U',IFORM,NXX,NYY,NSLICE,

     &                MAXIM,' ',.TRUE.,IRTFLG)
         IF (IRTFLG .NE. 0) RETURN

         ALLOCATE (QB(LSM,NROW), STAT=IRTFLG)
        
         IF (IRTFLG.NE.0) THEN
            CALL ERRT(45,'TF SIM',NE)      
            GOTO 9998
         ENDIF  


         SCX = 2.0 / NSAM
         SCY = 2.0 / NROW

         CS  = CS*1.E7
 
         DO  K=1,NROW
            KY = K-1
            IF (K.GT.IYC) KY = KY-NROW
    
            DO  I=1,LSM,2
               KX = (I-1)/2

               AK = KM * SQRT((KX*SCX)**2 + (KY*SCY)**2)
               
C            CALCULATE GAMMA (QQT)
               F1=1./SQRT(CS*LAMBDA)
               F2=SQRT(SQRT(CS*LAMBDA**3))
               AKK=AK*F2

               IF (KX.NE.0) THEN
                  AZ = ATAN2(FLOAT(KY),FLOAT(KX)) + QUADPI/2.
               ELSE
                  AZ =  QUADPI/2.
               ENDIF
               
               AZR = AZZ*(QUADPI/180.)
               DZZ = DZ+DZA/2*SIN(2*(AZ-AZR))

               DZ1=DZZ*F1
               QQT=2.*QUADPI*(.25*AKK**4-.5*DZ1*AKK**2)

C            CALCULATE PRODUCT OF ENVELOPE FUNCTION AND GAUSSIAN FALL OFF
               Q1 = (Q*F2)**2
               ENV1=ENV/F2**2  
                                       
               DS1=DS*DS*F1
               KAPPA=DS1*FP/KPP
               PDD=AKK**3-DZ1*AKK
               CH=EXP(AKK**4*KAPPA)
               EF_GFO=(EXP(FP*Q1*PDD**2)*CH)*EXP(-ENV1*AKK**2)
               EF=(EXP(FP*Q1*PDD**2)*CH) 

C             THE AK USED FOR PHI AND FMOD IS MODIFIED HERE TO AKR AS 
C             PER HAASE'S DEFINITION OF WAVE VECTOR.
C             THE KX AND KY REPRESENTING UNMODIFIED AK ARE USED FOR CALCULATING
C             OR'S AND OI'S (SCALAR PRODUCT K.R = KX.X + KY.Y)AS PER J. FRANK'S
C             "THREE DIMENSIONAL ELECTRON MICROSCOPY OF MACROMOLECULES" BOOK 
C             DEFINITION FOR OR AND OI.   

                AKR = AK*2.0* QUADPI
                OR_REAL = 0.0
                OR_IMG  = 0.0
                OI_REAL = 0.0
                OI_IMG  = 0.0
                            
               DO IJ = 1, NATOM
               
C             LOOK AT THE ATOMLOC(IK,2) TO FIND WHAT ATOM IS IT
C             LOOK AT SMALL_A AND SMALL_B TABLE TO LOCATE THE ATOM
C             USE FORMULA TO CALCULATE "FMOD" AND "PHI"
C             CALCULATE F_REAL AND F_IMG
C             O_REAL = O_REAL + FORMULA FOR O_REAL
C             O_IMG = O_IMG + FORMULA FOR O_IMG      
  
               
                  KANUM = ATOMLOC(IJ,2)
 
C               NO ENTRY FOR HYDROGEN. LIGHT ATOM SO CAN BE NEGLECTED             
                  IF (KANUM .NE. 1) THEN                  
                   
                     DO LA = 1, NTYPES_ATOM
                        IF(SMALL_A(1,LA) .EQ. KANUM) GOTO 25                                 
                     ENDDO
25                   CONTINUE   

                     DO LB = 1, NTYPES_ATOM
                        IF(SMALL_B(1,LB) .EQ. KANUM) GOTO 30                                 
                     ENDDO
30                   CONTINUE   

 

                         
C                  AK HAS UNITS 1/ANGSTROM; ARGUMENTS OF EXP() SHOULD BE A REAL NUMBER !          
                     FMOD = EXP(SMALL_A(2,LA) + SMALL_A(3,LA)*AKR +       
     &                 SMALL_A(4,LA)*AKR*AKR + 
     &                      SMALL_A(5,LA)*AKR*AKR*AKR +
     &                            SMALL_A(6,LA)*AKR*AKR*AKR*AKR)
     
C                  PHI SHOULD BE IN RADIAN !                 
                     PHI = SMALL_B(2,LB) + SMALL_B(3,LB)*AKR +    
     &                 SMALL_B(4,LB)*AKR*AKR + 
     &                     SMALL_B(5,LB)*AKR*AKR*AKR +
     &                         SMALL_B(6,LB)*AKR*AKR*AKR*AKR
 
 
                     F_REAL = FMOD  * COS(PHI)
                     F_IMG  = FMOD  * SIN(PHI)
                                     
                     THETA = 2.0*QUADPI*(KM*SCX*KX*ATOMLOC(IJ,3) + 
     &                                      KM*SCY*KY*ATOMLOC(IJ,4))
                      
                     OR_REAL = OR_REAL + F_REAL * COS(THETA)
C                  RESULTS IN ROTATION OF THE IMAGE BY 180 DEGREES      
C                  OR_IMG  = OR_IMG  - F_REAL * SIN(THETA)
                     OR_IMG  = OR_IMG  + F_REAL * SIN(THETA)


                     OI_REAL = OI_REAL + F_IMG * COS(THETA)
C                  RESULTS IN ROTATION OF THE IMAGE BY 180 DEGREES                   
C                  OI_IMG  = OI_IMG  - F_IMG * SIN(THETA)                 
                     OI_IMG  = OI_IMG  + F_IMG * SIN(THETA)               
 
                  ENDIF                        
               ENDDO
              
C            USES BOTH PHASE AND AMPLITUDE
               IF (TYPE_IMAGE .EQ. 'B') THEN          
                  TI_REAL = SIN(QQT) * OR_REAL - COS(QQT) * OI_REAL
                  TI_IMG  = SIN(QQT) * OR_IMG  - COS(QQT) * OI_IMG
C               REVERSE SIGN TO HAVE SAME DENSITY CONVENTION AS SPIDER
                  TI_REAL = -TI_REAL
                  TI_IMG = -TI_IMG                        
C             USES PHASE ONLY                             
               ELSE IF (TYPE_IMAGE .EQ. 'P')  THEN        
                  TI_REAL = SIN(QQT) * OR_REAL
                  TI_IMG = SIN(QQT) * OR_IMG
C               REVERSE SIGN TO HAVE SAME DENSITY CONVENTION AS SPIDER
                  TI_REAL = -TI_REAL
                  TI_IMG = -TI_IMG                        
C             USES AMPLITUDE ONLY                
               ELSE   
                  TI_REAL = -COS(QQT) * OI_REAL
                  TI_IMG = -COS(QQT) * OI_IMG
C               REVERSE SIGN TO HAVE SAME DENSITY CONVENTION AS SPIDER
                  TI_REAL = -TI_REAL
                  TI_IMG = -TI_IMG                
               ENDIF


               IF (CHANGE_INT .EQ. 'B')  THEN           
C                MULTIPLY INTNSITY BY BOTH ENVELOPE FUNCTION AND GAUSSIAN FALL OFF                     
                  TI_REAL = TI_REAL * EF_GFO    
                  TI_IMG =  TI_IMG * EF_GFO      
               ELSE IF (CHANGE_INT .EQ. 'E')  THEN 
C                MULTIPLY INTNSITY BY ONLY THE ENVELOPE FUNCTION                       
                  TI_REAL = TI_REAL * EF
                  TI_IMG =  TI_IMG * EF                
               ENDIF
               
               QB(I,K)   = TI_REAL
               QB(I+1,K) = TI_IMG   
C             TO PRINT OUT THE GAMMA (QQT) USE THE FOLLOWING AND USE SPIDER
C             SPIDER COMMAND "PW" TO SEE THE CTF          
C              QB(I,K)   = qqt
C              QB(I+1,K) = 0                         
            ENDDO        
            
         ENDDO

C        DO INVERSE FOURIER TRANSFORM
         INV = -1
         NSLICE = 1
         CALL FMRS_2(QB,NSAM,NROW,INV)
         IF (INV .EQ. 0)THEN
            CALL ERRT(38,'TRAFSIM.F',NE)
            CLOSE(LUN)
            GOTO 9998
         ELSE  
                                     
         CALL WRITEV(LUN,QB,LSM,NROW,NSAM,NROW,NSLICE) 
          
         ENDIF                 
                            
9998     CLOSE(LUN)      
9999     CLOSE(LUN1)

         END


 
         
         
        
