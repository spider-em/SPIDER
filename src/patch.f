
C++*********************************************************************
C
C  PATCH.F      USED REG_SET                        AUG 00 ARDEAN LEITH
C               ADDED 3D to 'PA' & 'IN'             FEB 01 ARDEAN LEITH
C               0,0 BUG                             FEB 03 ARDEAN LEITH
C               1,1 'PD' SPEEDUP                    NOV 03 ARDEAN LEITH
C               1,1 BACKGROUND BUG                  JAN 07 ARDEAN LEITH
C               BORDER CALC BUG                     FEB 11 ARDEAN LEITH
C               WRITE  PAD INTO REGISTER ALWA       FEB 11 ARDEAN LEITH      
C               IREC BUG WHEN OUTSIDE               OCT 11 ARDEAN LEITH      
C               FMIN BUG on 'IN S'                  JAN 14 ARDEAN LEITH      
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2014  Health Research Inc.,                         *
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
C   PATCH(LUN2,LUN1,NSAM2,NROW2,NSAM1,NROW1,LOCX,LOCY,IN,AVS,OPT)
C
C   PURPOSE: PATCH, PAD OR INSERT IMAGE
C
C      THIS SUBROUTINE TAKES IMAGE STORED ON LUN1 AND ADDS IT ONTO OR
C      INSERTS IT INTO IMAGE STORED ON LUN2 AT A SPECIFIED LOCATION.
C
C   PARAMETERS:
C        LUN2          SMALL INPUT IMAGE                     (SENT)
C        LUN1          BIG INPUT IMAGE (= OUTPUT IMAGE)      (SENT)
C        NSAM1,NROW1   SIZE OF BIG IMAGE 1                   (SENT)
C        NSLIC1                                              (SENT)
C        NSAM2,NROW2   SIZE OF IMAGE 2 (NSAM2 <= NSAM1)      (SENT)
C        NSLIC2                        (NROW2 <= NROW1)      (SENT)
C        LOCX,LOCY   UL COORDINATES OF PATCH IN BIG IMAGE  (SENT)
C        IN            0   PATCH                       'PA'  (SENT)
C                      1   INSERT                      'IN'
C                      2   CREATE IMAGE WITH VALUE AVS 'PD'
C                          AND IMAGE FROM LUN2 INSERTED.       
C 		       3   SAME AS 2 WITH CIRCULAR PAD  'PD'
C        AVS           PADDING VALUE                         (SENT)
C        OPT           FCHAR(4:4)                            (SENT)
C        FMIN2,FMAX2   MIN & MAX FOR SMALL INPUT FILE        (SENT)
C
C--********************************************************************

      SUBROUTINE PATCH(LUN2,LUN1,NSAM2,NROW2,NSLIC2,
     &                           NSAM1,NROW1,NSLIC1,
     &                           LOCX,LOCY,LOCZ,
     &                           IN,AVS,OPT,
     &                           FMIN2,FMAX2,USEBORDER)

      INCLUDE 'CMBLOCK.INC'    
      INCLUDE 'CMLIMIT.INC'

      COMMON /IOBUF/ BUF(NBUFSIZ)

      INTEGER            :: IN
      CHARACTER          :: OPT
      REAL               :: AVS,FMIN2,FMAX2 
      LOGICAL            :: USEBORDER

C     AUTOMATIC ARRAYS
      REAL               :: BUF1(NSAM1)
      REAL               :: BUF2(NSAM1)

      DOUBLE PRECISION   :: RRRBOR

      CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

C     SKIP DIRECTLY TO CIRCULAR PADDING IF DESIRED --------------
      IF (IN == 3) GOTO 300


      NSAMP  = NSAM1 + 1

C     TRUNCATE INSERTION AT BOUNDARIES
      KS = LOCX
      IF (KS <= 0)    KS = 1
      KE = KS + NSAM2 - 1
      IF (KE > NSAM1) KE = NSAM1

      NS = LOCY
      IF (NS <= 0)    NS = 1
      NE = NS + NROW2 - 1
      IF (NE > NROW1) NE = NROW1

      IZGO = LOCZ        
      IF (IZGO <= 0)     IZGO = 1
      IZEND = IZGO + NSLIC2 - 1 
      IF (IZEND > NSLIC1) IZEND = NSLIC1

C     FOR UNCHANGED CODE FOR PADDING ETC.
      KKS = IZGO        
      KKE = IZEND 

      IF (IN == 0) THEN
C        FOR PATCH OPERATION --------------------------------------'PA'

         DO ISLICE = LOCZ,NSLIC1
           IF (ISLICE < 1) CYCLE
           ISLICE2 = ISLICE - LOCZ + 1
           IF (ISLICE2 > NSLIC2) CYCLE

           DO IROW = LOCY, NROW1
              IF (IROW < 1) CYCLE
              IROW2 = IROW - LOCY + 1
              IF (IROW2 > NROW2) CYCLE

              IREC1 = (ISLICE-1) *NROW1 + IROW
              IREC2 = (ISLICE2-1)*NROW2 + IROW2
 
              CALL REDLIN(LUN1,BUF1,NSAM1,IREC1)
              CALL REDLIN(LUN2,BUF2,NSAM2,IREC2)

              DO K = KS,KE
C                 ADD THE TWO PIXELS
                  BUF1(K) = BUF1(K) + BUF2(K-KS+1)
              ENDDO

              CALL WRTLIN(LUN1,BUF1,NSAM1,IREC1)
            ENDDO
         ENDDO
         RETURN

      ELSEIF (IN == 1) THEN
C       FOR INSERT OPERATION --------------------------------------'IN'

60	IF (OPT == 'S') THEN 
           SCALE  = (FMAX - FMIN) / (FMAX2 - FMIN2)
           SCALE2 = FMIN - FMIN2 * SCALE

           !write(6,*) ' '
           !write(6,*) ' fmin, fmin2: ', fmin,fmin2
           !write(6,*) ' fmax, fmax2: ', fmax,fmax2
           !write(6,*) ' scale,scale2:',scale,scale2

         ELSEIF (OPT == 'M') THEN
	   FMINS = FMIN
	   FMAXS = FMAX
	   SCALE = 2.0 / (FMAXS-FMINS)
         ENDIF
        !write(6,*) ' islice2,irow2,irec2:',islice2,irow2,irec2
        !write(6,*) ' izgo,izend, ns,ne:',izgo,izend, ns,ne


         DO ISLICE = LOCZ,NSLIC1
           IF (ISLICE < 1) CYCLE
           ISLICE2 = ISLICE - LOCZ + 1
           IF (ISLICE2 > NSLIC2) CYCLE

           DO IROW = LOCY, NROW1
              IF (IROW < 1) CYCLE
              IROW2 = IROW - LOCY + 1
              IF (IROW2 > NROW2) CYCLE

              IREC1 = (ISLICE-1) *NROW1 + IROW
              IREC2 = (ISLICE2-1)*NROW2 + IROW2
 
              !write(6,*) ' islice2,irow2,irec2:',islice2,irow2,irec2

              CALL REDLIN(LUN1,BUF1,NSAM1,IREC1)
              CALL REDLIN(LUN2,BUF2,NSAM2,IREC2)

              IF (OPT .NE. 'S' .AND. OPT .NE. 'M') THEN
C                SIMPLE 'IN'
                 DO K = KS,KE
C                   REPLACE THE PIXEL
                    BUF1(K) = BUF2(K-KS+1)
                 ENDDO

              ELSE
C                'IN S'  OR   'IN M' ('IN M' IS UNDOCUMENTED!)

  	         DO K = KS,KE
C 	            BUF(K) = FMIN + (BUF(NSAMP+K-KS) - FMIN2) * 
C    &                       (FMAX - FMIN) / (FMAX2 - FMIN2)
C                   SCALE = (FMAX - FMIN) / (FMAX2 - FMIN2)
C 	            BUF(K) = FMIN + (BUF(NSAMP+K-KS) - FMIN2) * SCALE
C 	            BUF(K) = FMIN + BUF(NSAMP+K-KS) * SCALE - FMIN2 * SCALE
C                   SCALE2 = FMIN - FMIN2 * SCALE

C                   REPLACE THE PIXEL
                    BUF1(K) = BUF2(K-KS+1) * SCALE + SCALE2
                 ENDDO
              ENDIF

              CALL WRTLIN(LUN1,BUF1,NSAM1,IREC1)
           ENDDO
        ENDDO
        
        RETURN
      ENDIF

C     FOR PADDING ----------------------------------------------- "PD' 
C     MODIFIED FOR 3D PADDING (11/11/86. JMC)
C     FIRST FILL WITH AVS ALL THE PLANES BETWEEN 1 AND (LOCZ-1)
C     AVS MAY BE THE AVERAGE IN THE OUTER TWO-PLANES OF THE VOLUME

C     IT APPEARS THAT AVS WAS ALSO USED TO SIGNAL (AVS .LT. 99999). 
C     BAD BUG IF ACTUAL AVS WAS > 99999 Nov 03 al

      IF (USEBORDER) THEN 
C        CALCULATE THE MEAN (=AVS) OF THE 2 OUTSIDE VOXELS OF THE FILE 
C        (IT IS DONE THIS WAY TO AVOID PASSING A NEW VARIABLE IN A 
C         PROGRAM THAT IS CALLED SOMEWHERE ELSE ALSO).

         KKKCOU = 0
         RRRBOR = 0

         IF (NSLIC2 .NE. 1) THEN
C           3D VOLUME

C                      END-SLICES          END-COLUMNS
            KKKCOU = (2*NSAM2)*NROW2*2 + (4*NROW2)*(NSLIC2-4)
C                 END-SLICES       END-ROWS            END-COLUMNS
            KKK=NSAM2*NROW2*4 + (NSLIC2-4)*NSAM2*4 +(4*NROW2)*(NSLIC2-4)

            DO  KKSLIC=1,2           ! LOOP OVER END 2 SLICES
               IREC1 = (KKSLIC-1)*NROW2 
               IRECN = (NSLIC2-KKSLIC)*NROW2 
               DO  KKROW=1,NROW2
                  CALL REDLIN(LUN2,BUF,         NSAM2,IREC1+KKROW)
                  CALL REDLIN(LUN2,BUF(NSAM2+1),NSAM2,IRECN+KKROW)

                  RRRBOR = RRRBOR + SUM(BUF(1:2*NSAM2))
	       ENDDO
	    ENDDO

            DO  KKSLIC=3,NSLIC2-2   ! LOOP OVER NON-END SLICES
               IREC1 = (KKSLIC-1)*NROW2 
               IRECN = (NSLIC2-KKSLIC)*NROW2 

               DO  KKROW=1,2        ! LOOP OVER END 2 ROWS
                  CALL REDLIN(LUN2,BUF,         NSAM2,IREC1+KKROW) 
                  CALL REDLIN(LUN2,BUF(NSAM2+1),NSAM2,IRECN+KKROW)

                  RRRBOR = RRRBOR + SUM(BUF(1:2*NSAM2))
               ENDDO

               DO  KKROW=1,NROW2    ! LOOP OVER END COLUMNS
                  CALL REDLIN(LUN2,BUF,         NSAM2,IREC1+KKROW)
 
                  RRRBOR = RRRBOR+BUF(1)+BUF(2)+BUF(NSAM2-1)+BUF(NSAM2)
	       ENDDO
	    ENDDO

          ELSE
C           2D IMAGE NOT A VOLUME
C                 END-ROWS    END-COLUMNS
            KKK = NSAM2*4  +  4*(NROW2-4)

            DO  KKROW=1,2           ! LOOP OVER END ROWS     
               CALL REDLIN(LUN2,BUF,NSAM2,         KKROW)
               CALL REDLIN(LUN2,BUF(NSAM2+1),NSAM2,NROW2-(KKROW-1))

               RRRBOR = RRRBOR + SUM(BUF(1:2*NSAM2))
	    ENDDO

            DO  KKROW=3,(NROW2-2)  ! LOOP OVER NON-END ROWS 
               CALL REDLIN(LUN2,BUF,NSAM2,KKROW)
               RRRBOR = RRRBOR + BUF(1)+BUF(2)+BUF(NSAM2)+BUF(NSAM2-1)
	    ENDDO
         ENDIF

C        COMPUTE AVERAGE
         AVST = RRRBOR / KKK

         IF (MYPID <= 0) WRITE (NOUT,3555) AVST,KKK
3555     FORMAT('  OUTER AVERAGE: ',1PG15.5,'  VOXELS:',I8)

      ELSE
 
         AVST = AVS
      ENDIF     ! END OF: IF (USEBORDER) 

C     WRITE THE PAD INTO OPTIONAL REGISTER  
      CALL REG_SET_NSEL(1,1, AVST,0.0,0.0,0.0,0.0,IRTFLG)
 
      IF (LOCX == 1 .AND. LOCY == 1 .AND. LOCZ == 1 .AND.
     &    NSAM2  <= NSAM1  .AND. NROW2 <= NROW1 .AND. 
     &    NSLIC2 <= NSLIC1 .AND. IN == 2) THEN
C         FAST & SIMPLE PAD AROUND UPPER LEFT CORNER

C         SET PAD TO AVST
          IF (NSAM2 .LT. NSAM1) BUF(NSAM2+1:NSAM1) = AVST

          IREC1 = 0
          DO ISLICE2 = 1,NSLIC2
             IGO2 = (ISLICE2 -1) * NROW2 + 1
             DO IREC2 = IGO2,IGO2+NROW2-1
                CALL REDLIN(LUN2,BUF,NSAM2,IREC2)
                IREC1 = IREC1 + 1
                CALL WRTLIN(LUN1,BUF,NSAM1,IREC1)
             ENDDO

             IF (NROW2 .LT. NROW1) THEN
C               STILL HAVE SOME PAD ROWS LEFT FOR OUTPUT
                BUF(1:NSAM1) = AVST

                DO IREC = 1,NROW1-NROW2 
                   IREC1 = IREC1 + 1
                   CALL WRTLIN(LUN1,BUF,NSAM1,IREC1)
                ENDDO
             ENDIF
          ENDDO

          IF (NSLIC2 .LT. NSLIC1) THEN
C            STILL HAVE SOME PAD SLICES LEFT FOR OUTPUT
             BUF(1:NSAM1) = AVST

             DO IREC = NROW1 * NSLIC2 + 1, NROW1 * NSLIC1
                CALL WRTLIN(LUN1,BUF,NSAM1,IREC)
             ENDDO
          ENDIF
          RETURN
      ENDIF

      NSAMP  = NSAM1 + 1
      NS     = LOCY
      NE     = NS + NROW2 - 1

      IF (KKS .NE. 1) THEN
         KKS1 = KKS-1
         DO  JJ=1,NSAM1
             BUF(JJ) = AVST
         ENDDO
         DO  KK=1,KKS1
            DO  K=1,NROW1
               CALL WRTLIN(LUN1,BUF,NSAM1,K+(NROW1*(KK-1)))
	    ENDDO
         ENDDO
      ENDIF

C     THE PLANES BETWEEN LOCZ AND LOCZ WILL BE REALLY PADDED

      DO 130 KK=KKS, KKE

         IF (NS.NE.1) THEN
            NS1=NS-1
            DO  I=1,NS1
               DO  K=1,NSAM1
                   BUF(K) = AVST
	       ENDDO
               CALL WRTLIN(LUN1,BUF,NSAM1,I+(NROW1*(KK-1)))
	    ENDDO
         ENDIF

         DO  I = NS,NE
            IF (KS==1) GOTO 160
            KS1 = KS - 1
            DO  K = 1,KS1
              BUF(K) = AVST
	    ENDDO
160         CALL REDLIN
     &         (LUN2,BUF(NSAMP),NSAM2,(I-NS+1)+(NROW2*(KK-KKS)))
            DO  K = KS,KE
               BUF(K) = BUF(NSAMP+K-KS)
	    ENDDO
            IF (KE. EQ. NSAM1) GOTO 179
            KE1 = KE + 1
            DO  K = KE1,NSAM1
               BUF(K) = AVST
	    ENDDO
179         CALL WRTLIN(LUN1,BUF,NSAM1,I+(NROW1*(KK-1)))
	 ENDDO

         IF (NE==NROW1) GOTO 130
         IE1 = NE + 1
         DO  I = IE1,NROW1
            DO  K = 1,NSAM1
               BUF(K) = AVST
	    ENDDO
         CALL WRTLIN(LUN1,BUF,NSAM1,I+(NROW1*(KK-1)))
	 ENDDO

130   CONTINUE

C     FILL THE REST OF THE PLANES BETWEEN LOCZ AND NSLIC1 WITH AVST 
       IF (KKE.NE.NSLIC1) THEN

          KKE1=KKE+1
          DO  I=1,NSAM1
             BUF(I) = AVST
          ENDDO

          DO  KK=KKE1,NSLIC1
            DO  K=1,NROW1
               CALL WRTLIN(LUN1,BUF,NSAM1,K+(NROW1*(KK-1)))
	    ENDDO
         ENDDO
      ENDIF
      RETURN


300   NSAMP  = NSAM1 + 1
      NS     = LOCY
      NE     = NS + NROW2 - 1


C INDEX ASSIGNMENTS FOR CIRCULAR PADDING --------------------------- PD
C
C	     KS2        KE2              KS1        KE1
C	NS2  I----------------------------------------I
C	     I          .                .
C	     I          .                .
C	     I          .                .
C	NE2  I...........                .............
C	     I
C	     I
C	     I
C	NS1  I...........
C	     I          .
C	     I          .
C	NE1  I------------------------------------------I
C
C
      NE2 = 0
      NS2 = 0
      IF (NS.LT.0) NS = NS+NROW1
      NS1 = NS
      NE1 = NS+NROW2-1
      NPASS = 1
      IF (NE1 <= NROW1) GOTO 320
      NS2 = 1
      NE2 = NE1-NROW1
      NE1 = NROW1
320   KS = LOCX
      KE2 = 1
      KS2 = 0
      IF (KS.LT.0) KS = KS+NSAM1
      KS1 = KS
      KE1 = KS+NSAM2-1
      IF (KE1<=NSAM1) GOTO 3205
      KS2 = 1
      KE2 = KE1 - NSAM1
      KE1 = NSAM1
3205  ISTART = NS2
      IEND = NE2
      IFILLS = NE2+1
      IFILLE = NS1-1
      IRECOF = NROW1-LOCY+1
      NPASS = 1
321   IF (ISTART==IEND) GOTO 350
      DO 340 I = ISTART,IEND
         IRECIN = I-ISTART+IRECOF+1
         CALL REDLIN(LUN2,BUF(NSAMP),NSAM2,IRECIN)
         IF (KS2==0) GOTO 330
         NOFFS = NSAMP+KE2-2*KS2
         DO  K = KS2,KE2
           BUF(K) = BUF(K+NOFFS+1)
	 ENDDO

         DO  K = KE2+1,KS1-1
           BUF(K) = AVST
	 ENDDO
330      NOFFS = NSAMP-KS1
         DO  K = KS1,KE1
           BUF(K) = BUF(K+NOFFS)
	 ENDDO

         IF (KE1==NSAM1) GOTO 340
         DO  K = KE1+1,NSAM1
           BUF(K) = AVST
	 ENDDO

340   CALL WRTLIN(LUN1,BUF,NSAM1,I)

350   IF (IFILLS==IFILLE) GOTO 390
      DO  I = IFILLS,IFILLE
         DO  K = 1,NSAM1
           BUF(K) = AVST
	 ENDDO
      CALL WRTLIN(LUN1,BUF,NSAM1,I)
      ENDDO

390   GOTO (400,410),NPASS
400   NPASS  = 2
      ISTART = NS1
      IEND   = NE1
      IFILLS = NE1 +1
      IFILLE = NROW1
      IRECOF = 0
      GOTO 321

410   RETURN
      END



