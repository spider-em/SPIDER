C++*********************************************************************
C
C SPEAKC.F               REVISED                  5/21/86  MR
C                        RENAMED FROM SPEAKG, I2  8/18/87  MR   
C                        USED REG_GET             AUG 00   ARDEAN LEITH
C                        PARAMETERS               FEB 2001 ARDEAN LEITH
C                        NPTR MOVED AWAY FROM BUF SEP 2003 ARDEAN LEITH
C                        IPAR BUG                 SEP 2007 ARDEAN LEITH
C                        LUNDOCWRTDAT REWRITE     DEC 2012 ARDEAN LEITH
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
C  PURPOSE
C
C       THIS SUBROUTINE SEARCHES FOR THE ML HIGHEST PEAKS IN THE (REAL)
C       IMAGE FILNAM AND PRINTS OUT POSITIONS AND VALUES OF THESE PEAKS.
C       HAS SUB-PIXEL ACCURACY USING BOTH CG AND POLYGONAL WEIGHTING
C       EXCLUSIVELY
C
C     SPEAKC(FILNAM,LUN,NX,NY,MAXDIM,OPTO,NDOCOUT,ML,NOR)
C         FILNAM     FILE NAME
C         LUN        LOGICAL UNIT NUMBER OF IMAGE
C         NX,NY  DIMENSIONS OF IMAGE
C         MAXDIM     MAXIMUM BUFFER SPACE AVAILABLE
C         OPTO       OUTPUT OPTION
C         OPTO=' '   DEFAULT: NO DOCUMENT OUTPUT
C         OPTO='X'   I.E.,REGISTER LIST FOLLOWING):
C                    OUTPUT OF POSITION & VALUE OF PEAK IN REGISTERS
C         OPTO='D'   DOCUMENT OUTPUT: NUMBER,POSITION, AND VALUE
C                    OF PEAK ARE WRITTEN INTO A DOCUMENT FILE
C         NDOCOUT    I/O UNIT NUMBER FOR DOC FILE
C         ML         NUMBER OF PEAKS WANTED
C         NOR        NEW ORIGIN WANTED
C    
C     REGISTER POSITIONS
C          1= CG X-SHIFT
C          2= CG Y-SHIFT
C          3= ABSOLUTE PEAK HEIGHT
C          4= RATIO PEAK HEIGHT
C          5= SUB-PIXEL NON-CG FLOATING PT. X-SHIFT
C          6= SUB-PIXEL NON-CG FLOATING PT. Y-SHIFT
C          7= VALUE OF  NON-CG EXTREMUM OF SUB-PIXEL PARABOLOID
C
C     DOC FILE POSITIONS FOR 'PK DC'
C          1= CG X-SHIFT
C          2= CG Y-SHIFT
C          3= ABSOLUTE PEAK HEIGHT
C          4= RATIO PEAK HEIGHT
C          5= NEGATIVITY WARNING FLAG
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE SPEAKC(FILNAM,LUN,NX,NY,MAXDIM,OPTO,NDOCOUT,ML,NOR)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         COMMON /IOBUF/ BUF(NBUFSIZ)
         COMMON NPTR(64)

         CHARACTER (LEN=*)      :: FILNAM
         INTEGER                :: LUN,NX,NY,MAXDIM
         CHARACTER              :: OPTO
         INTEGER                :: NDOCOUT,ML,NOR

         INTEGER, DIMENSION(ML) :: IPOS
         REAL,    DIMENSION(ML) :: XILIST,YILIST, XLIST,YLIST
         REAL,    DIMENSION(ML) :: ALIST,RLIST
         REAL,    DIMENSION(9) ::  RSQ

         REAL                   :: DLIST(6)
         CHARACTER              :: POSE,NULL

         CHARACTER (LEN=80)     :: COMMEN
         CHARACTER (LEN=MAXNAM) :: DOCOUT
         LOGICAL                :: ASKNAM,ADDEXT,ISOLD
         LOGICAL                :: WRTCOM,APPEND,NEWFILE

         CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

         NULL = CHAR(0)

         IF (MAXDIM < 3*NX+200) THEN
            CALL ERRT(102,'INSUFFICIENT BUFFER SPACE',MAXDIM)
            RETURN
         ENDIF

         DO  K = 1,ML
            IPOS(K)  = 0
            XILIST(K) = 0.0
            YILIST(K) = 0.0
            RLIST(K) = 0.0
            ALIST(K) = -HUGE(XVAL)
         ENDDO

         NTAB  = 1
         NXCTR = NX/2+1
         NYCTR = NY/2+1
         IF (NOR .NE. 0) THEN

            CALL RDPRIS(NXCTR,NYCTR,NOT_USED,
     &                  'ORIGIN COORDINATES',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            CALL RDPRI1S(NTAB,NOT_USED,
     &             'PEAK NUMBER FOR RATIO',IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            IF (NTAB <= 0 .OR. NTAB > ML) THEN
               CALL ERRT(102,'NOT CONTAINED IN TABLE',NTAB)
               RETURN
            ENDIF
          ENDIF

30       CALL RDPRM2(CG,CG2,NOT_USED,
     &       'ELLIPSE AXES (X,Y) FOR CGR CALCULATION')
         IF (CG2 == 0) CG2 = CG
         CG2SQ = CG2*CG2
         CGSQ  = CG*CG
         CALL RDPRMC(POSE,NC,.TRUE.,'POSITIVITY ENFORCED? (Y/N)',
     &         NULL,IRTFLG)

         CALL RDPRM(FNEIGH,NOT_USED,'NEIGHBORHOOD DISTANCE')

         CALL RDPRMI(IEWX,IEWY,NOT_USED,'EDGE EXCLUSION WIDTH X, Y')
         IF (IEWY <= 0) IEWY = IEWX

C        MAKE SURE THAT THE EDGE EXCLUSION IS AT LEAST AS LARGE AS THE
C        CENTER OF GRAVITY AXES:
         IF (IEWX < CG)  IEWX = CG
         IF (IEWY < CG2) IEWY = CG2

         WRITE(NOUT,2001) IEWX,IEWY
2001     FORMAT('  EDGE EXCLUSION WIDTH USED: X: ',I4,'  Y: ',I4)

         IESAM = NX - IEWX
         IEROW = NY - IEWY

C        CALCULATE HOW MANY RECORDS HAVE TO BE READ IN TO COVER THE COMPLETE
C        FIELD USED IN THE CENTER OF GRAVITY CALCULATIONS:
         NREC = CG2 + 0.5

C        MAKE DIAMETER AN ODD NUMBER (=TOTAL NUMBER OF RECORDS TO
C        BE AVAILABLE AT ANY TIME)
         NREC = NREC * 2 + 1
         IF  (CG == 0.0) NREC=3

         NMAX   = 0
         NEG    = 0

C        INITIALIZE MREAD
         CALL MREAD(0,BUF,NX,NREC,NPTR)

         NREC2 = NREC/2+1
         NRECH = NREC/2

         DO  I = NREC2,NY-NREC2+1
            CALL MREAD(LUN,BUF,NX,NREC,NPTR)

C           GET BEGINNING LOCATION OF THREE LINES
            I1A = (NPTR(NREC2-1)-1)*NX
            I2A = (NPTR(NREC2)-1)*NX
            I3A = (NPTR(NREC2+1)-1)*NX

            DO  K = NREC2,NX-NREC2+1
               A = BUF(K+I2A)
C              CHECK IF A IS SMALLER THAN THE 8 SURROUNDING POINTS:
               IF (A <= BUF(K-1+I2A)) CYCLE
               IF (A <= BUF(K-1+I1A)) CYCLE
               IF (A <= BUF(K+I1A))   CYCLE
               IF (A <= BUF(K+1+I1A)) CYCLE
               IF (A <= BUF(K+1+I2A)) CYCLE
               IF (A <= BUF(K+1+I3A)) CYCLE
               IF (A <= BUF(K+I3A))   CYCLE
               IF (A <= BUF(K-1+I3A)) CYCLE

C              MAKE SURE THAT PEAK NOT NEAR EDGE (DEFINED BY IEW)
               IF (K <= IEWX .OR. K >= IESAM) CYCLE
               IF (I <= IEWY .OR. I >= IEROW) CYCLE

C              NEW LOCAL PEAK FOUND
               NMAX = NMAX + 1
c              write(0,*) ' nmax: ',nmax,'   ml: ',ml,'  A: ',A
c              write(0,*) '  '

C              CHECK IF PEAK IS LARGER THAN THE ML PEAKS FOUND PREVIOUSLY
C              AND PUT IT INTO THE CORRECT RANK POSITION

               DO L = 1,ML                     !L IS PEAK NUMBER

                 IF (A <= ALIST(L)) CYCLE   ! CYCLE IF SMALLER

                 IF (L .NE. ML) THEN
C                   NOT LAST PEAK NUMBER, REORDER THE ML PEAKS
                    L1 = L+1                       
                    DO J = L1,ML
                       JO        = ML-J+L1
                       JN        = ML-J+L1-1
                       ALIST(JO) = ALIST(JN)
                       XILIST(JO) = XILIST(JN)
                       YILIST(JO) = YILIST(JN)
                       IPOS(JO)  = IPOS(JN)
                    ENDDO     ! END OF: DO J = L1,ML
                 ENDIF        ! END OF: IF (L .NE. ML) THEN

C                STORE PEAK HEIGHT AND COORDINATES
                 IF (L == 1) THEN
                     IPAR = I
                     KPAR = K
                 ENDIF

                 ALIST(L) = A      
                 IPOS(L)  = 0
                 XILIST(L) = I
                 YILIST(L) = K

C                CENTER OF GRAVITY SEARCH
                 IF (CG == 0.0) EXIT
                 CNY   = 0
                 CNX   = 0
                 CGR   = 0
                 PCORR = 0.

                 IF (POSE == 'Y') THEN
C                   POSITIVITY OF AREA IS ENFORCED, DETERMINE AREA'S MINIMUM:
                    AMIN = 1.E10
                    DO ICG=1,NREC
                       IND1 = (NPTR(ICG)-1)*NX
                       YSQ  = (ICG-NREC2)**2
                       DO KCG=K-NRECH,K+NRECH
                          XSQ  = (KCG-K)**2
                          CRIT = XSQ/CGSQ+YSQ/CG2SQ
                          IF (CRIT <= 1) THEN
                             IND2 = IND1+KCG
                             XXX  = BUF(IND2)
                             IF (AMIN > XXX) AMIN = XXX
                          ENDIF
                       ENDDO    ! END OF: DO KCG=K-NRECH,K+NRECH
                    ENDDO       ! END OF: DO ICG=1,NREC
                    PCORR = AMIN
                 ENDIF          ! END OF: IF (POSE == 'Y') THEN


C                CALCULATE CENTER OF GRAVITY:
                 DO ICG=1,NREC
                    IND1 = (NPTR(ICG)-1) * NX
                    YSQ  = (ICG-NREC2) ** 2

                    DO  KCG=K-NRECH,K+NRECH
                       XSQ  = (KCG-K) **2
                       CRIT = XSQ / CGSQ + YSQ / CG2SQ

                       IF (CRIT > 1) CYCLE

C                      POINT INSIDE ELLIPTICAL AREA FOR CGR DETERMINATION:
                       IND2 = IND1 + KCG

                       ADD = BUF(IND2) - PCORR
                       IF (ADD < 0) THEN
                          IPOS(L) = 1
                          NEG     = NEG + 1
CCCC                      WRITE(NOUT,1999) ADD,ICG,KCG
1999                      FORMAT(' *** NEGATIVE VALUE ',G12.4,
     &                      '  FOUND AT PIXEL ',2I5,'  CGR MEANINGLESS')
                       ENDIF

                       CGR = CGR + ADD
                       CNY = CNY + FLOAT(ICG-NREC2)*ADD
                       CNX = CNX + FLOAT(KCG-K)*ADD
                    ENDDO    ! END OF: DO KCG=K-NRECH,K+NRECH 
                 ENDDO       ! END OF: DO ICG=1,NREC

                 IF (IPOS(L) == 0) THEN
C                   MEANINGFULL CGR
                    XILIST(L) = CNY / CGR + XILIST(L)
                    YILIST(L) = CNX / CGR + YILIST(L)
                    EXIT
                 ENDIF
              ENDDO          ! END OF: DO L = 1,ML
           ENDDO             ! END OF: DO K = NREC2,NX-NREC2+1
         ENDDO

         IF (NMAX == 0)  THEN
            IF (MYPID <= 0) WRITE(NDAT,*) '*** NO PEAKS FOUND'
            IF (MYPID <= 0 .AND. NOUT .NE. NDAT) 
     &                        WRITE(NOUT,*) '*** NO PEAKS FOUND'

            CALL REG_SET_NSEL(1,5,0.0, 0.0, 0.0, 0.0, 0.0,IRTFLG)
            CALL REG_SET_NSEL(6,2,0.0, 0.0, 0.0, 0.0, 0.0,IRTFLG)
            RETURN
         ENDIF

         IF (OPTO == 'D') THEN
            ADDEXT  = .TRUE.
            ASKNAM  = .TRUE.
            ISOLD   = .FALSE.
            APPEND  = .FALSE.
            WRTCOM  = .TRUE.
            CALL OPENDOC(DOCOUT,ADDEXT,NLET,
     &                  NDOCOUT,NICDOCOUT,ASKNAM,'OUTPUT DOCUMENT',
     &                  ISOLD,APPEND,WRTCOM,NEWFILE,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN
C
            COMMEN = 
     &       '            XO            YO        HEIGHT' //
     &       '         RATIO         ERRORS'
C             123456789 123456789 123456789 123456789 123456789 123456
            CALL LUNDOCPUTCOM(NDOCOUT,COMMEN(1:80),IRTFLG)
         ENDIF

               DLIST(1) = XLIST(L)
               DLIST(2) = YLIST(L)
               DLIST(3) = ALIST(L)
               DLIST(4) = RLIST(L)
               DLIST(5) = IPOS(L)


         IF (MYPID <= 0) THEN
            IF (NEG > 0) THEN
               WRITE(NOUT,304)NEG
304            FORMAT(/,'  WARNING: ',I6,' NEGATIVE PIXEL VALUES FOUND',
     &                 ' DURING CALCULATION OF CGR',/)
            ENDIF

            WRITE(NDAT,299)
            IF (NDAT .NE. NOUT) WRITE(NOUT,299)
299         FORMAT('    NO      IX      IY      XO      YO   HEIGHT',
     &   '     RATIO',5X,' 1 IF NEGATIVE'/5X,58X,'VALUE IN CGR AREA')
         ENDIF

         MLIST = MIN(ML,NMAX)

C        FOR EACH PEAK, MAKE SURE IT DOES NOT FALL WITHIN RADIUS=FNEIGH
C        OF A HIGHER-RANKING ONE

         FN2 = FNEIGH * FNEIGH
C        CHANGED 2/8/88 MR

         DO L=MLIST,2,-1
            FK = YILIST(L)
            FI = XILIST(L)
            DO LI=1,L-1
C              IF OUTIDE EXCLUSION ZONE, CYCLE
               IF ((FK-YILIST(LI))**2+(FI-XILIST(LI))**2 > FN2) CYCLE
               ALIST(L) = -HUGE(XVAL)
            ENDDO
         ENDDO
         NPKCNT = 0

C        NOW SELECT PEAKS THAT COMPLY WITH CONDITIONS 
         DO L=1,MLIST
            IF (ALIST(L) == -HUGE(XVAL)) THEN
C              PEAK HAS BEEN EXCLUDED
               IF (L .NE. NTAB) CYCLE
               NTAB = NTAB + 1

               IF (MYPID <= 0) WRITE(NOUT,319) NTAB
               IF (MYPID <= 0 .AND. NOUT .NE. NDAT) 
     &            WRITE(NOUT,319) NTAB
319            FORMAT(' NUMBER OF REF. PEAKS INCREASED BY ONE TO: ',I5)
               CYCLE
            ENDIF

            NPKCNT    = NPKCNT + 1
            XLIST(L) = YILIST(L) - NXCTR
            YLIST(L) = XILIST(L) - NYCTR
            RLIST(L)  = ALIST(L) / ALIST(NTAB)

            IF (MYPID <= 0) THEN
C              WRITE TO TERMINAL AND RESULTS FILE
               WRITE(NDAT,301) NPKCNT,YILIST(L),XILIST(L),
     &             XLIST(L),YLIST(L),ALIST(L),RLIST(L),IPOS(L)

301            FORMAT(I6, 4F8.2, G12.3, F10.5,5X,I2)

               IF (NOUT .NE. NDAT) THEN
                  WRITE(NOUT,301) NPKCNT,YILIST(L),XILIST(L),
     &               XLIST(L),YLIST(L), ALIST(L),RLIST(L),IPOS(L)
               ENDIF
            ENDIF     ! END OF: IF (MYPID <= 0) THEN

            IF (OPTO == 'D') THEN
C              WRITE PEAK TO DOCUMENT FILE
               DLIST(1) = XLIST(L)
               DLIST(2) = YLIST(L)
               DLIST(3) = ALIST(L)
               DLIST(4) = RLIST(L)
               DLIST(5) = IPOS(L)

C              PUSH DLIST INTO OUTPUT DOC. FILE
               CALL LUNDOCWRTDAT(NICDOCOUT,L,DLIST,5,IRTFLG)


               !!CALL SAVD(NDOC,DLIST,5,IRTFLG)
            ENDIF
         ENDDO        ! END OF: DO L=1,MLIST

         CLOSE(NDOCOUT)

C        SET RESULTS IN COMMAND LINE REGISTERS 1..4
         CALL REG_SET_NSEL(1,4,XLIST(NTAB),YLIST(NTAB),
     &                         ALIST(NTAB), RLIST(NTAB),0.0,IRTFLG)

C        9/25/81 PARABOLIC FIT TO THE 3X3 NEIGHBORHOOD OF THE PEAK POINT
C        PROGRAM SECTION SENT BY M.VAN HEEL, MODIFIED FOR SPIDER USE. JF

C        APPLIED ONLY TO HIGHEST RANKING PEAK

         KL = KPAR
         DO  I=1,3
            if (Ipar > NY .or. ipar <= 0) then
               write(6,*) ' IPAR: ',IPAR,' > NY:',NY
               write(6,*) ' I: ',I
               CALL ERRT(101,'BAD ROW REFERENCED',NE)
            endif

            IROW = IPAR + I - 2    ! GET ROW NUMBER
            IF (IROW < 1)    IROW = IROW + NY
            IF (IROW > NY) IROW = IROW - NY

C           READ ORIGINAL DATA
            CALL REDLIN(LUN,BUF,NX,IROW)
            I1 = (I-1)*3
            DO  K=1,3
               ISAM = KL+K-2
               IF (ISAM <1)     ISAM=ISAM+NX
               IF (ISAM > NX) ISAM=ISAM-NX

C              VALUES FOR NINE LOCATIONS AROUND THE PEAK & PEAK
               RSQ(I1+K)= BUF(ISAM)
            ENDDO
         ENDDO

C        FIND XSH & YSH (RANGE -1.0 .... 1.0)
         CALL PARABL(RSQ,XSH,YSH,PEAKV)

C        SET SUB-PIXEL RESULTS IN COMMAND LINE REGISTERS 5..7
         CALL REG_SET_NSEL(5,3,XSH+KPAR-NXCTR, YSH+IPAR-NYCTR,
     &                       PEAKV, 0,0,IRTFLG)

         IF (MYPID <= 0) THEN
            IF (VERBOSE) WRITE(NOUT,*) ' '
            WRITE(NDAT,302)XSH,YSH,PEAKV
            IF (NDAT .NE. NOUT) 
     &         WRITE(NOUT,302)XSH,YSH,PEAKV
302            FORMAT('  SUB-PIXEL OFFSET OF HIGHEST PEAK: (',
     &               F5.2,', ',F5.2,')  HEIGHT: ',G15.7)

            IDIFF = ML - NPKCNT
            WRITE(NOUT,401) ML,NPKCNT,IDIFF
401         FORMAT('  PEAKS SPECIFIED: ',I6,'  PEAKS PASSED: ',I6,
     &             '  PEAKS NEAR EDGE EXCLUDED: ',I6)
            IF (MYPID .LE. 0) WRITE(NOUT,*) ' '

         ENDIF

         END
