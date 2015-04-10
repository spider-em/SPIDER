C++*********************************************************************
C
C  NOISE.F   FIXED UNWORKING CODE                NOV 2014 ARDEAN LEITH
C
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
C NOISE:  CALCULATE BACKGROUND NOISE OF POWER SPECTRUM AND SUBTRACT
C         IT USING LEAST-SQUARE METHOD TO FIT POINTS INTO GAUSSIAN 
C         PROFILE USING STEEPEST DESCENT METHOD.
C         F(XI,A) = A1 * EXP(-(XI/A2)**2) + A3
C
C NOTE:   COULD NOT HAVE WORKED SINCE 1999 DUE TO LUN2 BUG  al NOV 2014
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        SUBROUTINE NOISE(LUN1,LUN2)

        IMPLICIT NONE

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        REAL :: BUF
        COMMON /IOBUF/ BUF(NBUFSIZ)   ! FROM CMLIMITS.INC

        INTEGER               :: LUN1,LUN2

        INTEGER, PARAMETER    :: MAXMINS = 120
        INTEGER, PARAMETER    :: MAXX  = NBUFSIZ
        INTEGER ,PARAMETER    :: LUNDOC  = 81

        REAL                  :: Y1(MAXX),Y2(MAXX)
        INTEGER               :: NUMLIS(MAXMINS)
        REAL                  :: X(MAXMINS),RKFR(MAXMINS),RMINS(MAXMINS)
        REAL                  :: SFR(MAXMINS)   ! UNUSED HERE
        REAL                  :: KM,KS
        CHARACTER(LEN=MAXNAM) :: OUTNAME,IMFILE,DOCOUT
        INTEGER               :: NE,NLUNDOC,IDUM,NUMT,IP

        INTEGER               :: MAXIM,NX,NY,NZ,IRTFLG,INUM,N,NUM,I,J
        INTEGER               :: NSTEP,NOT_USED,NLET
        REAL                  :: SPMAX,A1,A2,A3,P,X1,X0,DA1,DA2,DA3
        REAL                  :: DXA1,DXA2,XMIN,DXA3,SUM,D

        CHARACTER(LEN=30)     :: COMMEN
        LOGICAL               :: ASKNAM,NEWFILE,ADDEXT,ISOLD
        LOGICAL               :: WRTCOM,APPEND
        

C       OPEN 1D POWER SPECTRUM IMAGE FILE
        MAXIM = 0
        CALL OPFILEC(0,.TRUE.,IMFILE,LUN1,'O',IFORM,NX,NY,NZ,
     &             MAXIM,'POWER SPECTRUM 1D IMAGE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        WRITE(NOUT,10 ) NX,NY
10      FORMAT('  IMAGE DIMENSIONS:', I0,' x ',I0)

C       FIND MINIMA IN POWER SPECTRUM
        INUM = 1
        CALL DEFO003(INUM,N,RKFR,SFR,Y2,NX,SPMAX,LUN1,BUF,
     &               MAXMINS,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        IF (N <= 0) THEN
           CALL ERRT(101,'INVALID INPUT, NO MINIMA FOUND',I)
           RETURN
        ENDIF

C       REREAD INPUT 1D IMAGE 
        CALL REDLIN(LUN1,Y2,NX,1)

        CALL FLUSHRESULTS

        NUMT = N
        CALL RDPRAI(NUMLIS,MAXMINS,NUMT, IDUM,IDUM,
     &              'LIST OF MINIMA TO BE USED','F',IRTFLG)
        IF (IRTFLG .NE. 0) RETURN
        IF (NUMT <= 0) THEN
           DO I = 1,N
              NUMLIS(I) = I
           ENDDO
           NUMT = N
        ENDIF


        KM = SPMAX
        KS = KM / FLOAT(NX)
       
        NUM = 0
        DO I=1,NUMT

           IF (NUMLIS(I) < 1) THEN
              CALL ERRT(102,'MINIMA BELOW VALID RANGE',NUMLIS(I))
              RETURN
           ELSEIF (NUMLIS(I) > N) THEN
              WRITE(NOUT,'(A,I0)') 
     &              'INVALID MINIMA DISCARDED: ',NUMLIS(I)
              CYCLE
           ENDIF

           RMINS(I) = RKFR(NUMLIS(I))

           Y1(I) = (1.-(RMINS(I)-INT(RMINS(I)))) * BUF(INT(RMINS(I))) +
     &                 (RMINS(I)-INT(RMINS(I)))  * BUF(INT(RMINS(I))+1)
           X(I)  = RMINS(I) * KS
 
           !write(6,*) ' x:',x(i),'  y1:',y1(i), rmins(i)
           !write(6,*) ' x:',x(i),'  y1:',y1(i), buf(int(rmins(i)))
           NUM = NUM + 1
        ENDDO

C       OPEN OUTPUT 1D IMAGE
        MAXIM = 0
        CALL  OPFILEC(LUN1,.TRUE.,OUTNAME,LUN2,'U',IFORM, NX,1,1,
     &                MAXIM,'DENOISED PROFILE',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        ADDEXT  = .TRUE.
        ASKNAM  = .TRUE.
        ISOLD   = .FALSE.
        APPEND  = .FALSE.
        WRTCOM  = .TRUE.
        CALL OPENDOC(DOCOUT,ADDEXT,NLET,
     &               LUNDOC,NLUNDOC,ASKNAM,'OUTPUT DOCUMENT',
     &               ISOLD,APPEND,WRTCOM,NEWFILE,IRTFLG)

        IF (IRTFLG == -1) THEN
C          DO NOT WANT OUTPUT DOC FILE
           NLUNDOC = 0
        ENDIF
        IF (IRTFLG .NE. 0) GOTO 9999

C                 123456789 123456789 123456789 1234567890
        COMMEN = '   DENOISED VALUES'
        IF (NLUNDOC > 0) CALL LUNDOCPUTCOM(NLUNDOC,COMMEN(1:19),IRTFLG)


           
C       FOR  SINGLE MINIMUM, HACK TO AVOID JUMP INTO LOOP --------
        IF (NUM == 1) THEN

           CALL RDPRM1S(A2,NOT_USED,'A2 VALUE [A-1]',IRTFLG)
           IF (IRTFLG .NE. 0) RETURN

           IF (A2 == 0) THEN
              CALL ERRT(101,'INVALID INPUT GIVES DIVISION BY ZERO',I)
              GOTO 9999
           ENDIF

           A3 = 0
           A1 = Y1(1) / (EXP(-(X(1) / A2)**2))
           P  = A2 / KS

           GOTO 99

        ENDIF


C       HAS MORE THAN ONE MINIMA IN LIST ------------------------

C       SET INITIAL VALUE FOR A1,A2 AND A3

        A3 = 0             ! FOR NUM=2 ONLY
        IF (NUM > 2) A3 = Y1(NUM) * 0.8
       
        A2 = SQRT( (X(2)**2 - X(1) **2) / LOG(Y1(1) / Y1(2)) )
        A1 = (Y1(1) - A3) / EXP(-(X(1) / A2)**2)

c       write(6,*) ' Num:', num
c       write(6,*) ' Top:', (x(2)**2 - x(1)**2)
c       write(6,*) ' Bot:', log(y1(1)/y1(2))
c       write(6,*) ' X(1), x(2) :',x(1), x(2)
c       write(6,*) ' Y1(1), y1(2) :',y1(1), y1(2)
c       write(6,*) ' Y1(1)/y1(2)  :',y1(1)/y1(2)
c       write(6,*) ' Initial a1 :',a1
c       write(6,*) ' Initial a2 :',a2
c       write(6,*) ' Initial a3 :',a3
       
C       SET INITIAL VALUE OF X**2
        DO  I=1,NUM
          Y2(I) = A1 * EXP(-(X(I)/A2)**2) + A3
        ENDDO

        X0 = 0
        DO I=1,NUM
           X0 = X0 + (Y1(I) - Y2(I))**2
        ENDDO

        NSTEP = 0
999     DA1   = 0.001 * A1
        DA2   = 0.001 * A2
        DA3   = 0.001 * A3

C       CALCULATE dX**2 / dA1
        DO  I=1,NUM
           Y2(I) = (A1 + 0.1 * DA1) * EXP(-(X(I)/A2)**2) + A3
        ENDDO

        X1 = 0
        DO  I=1,NUM
           X1 = X1 + (Y1(I) - Y2(I))**2
        ENDDO

        DXA1 = (X1-X0) / (0.1*DA1)
c       write(6,*) ' dxa1:',dxa1,'=',x1,x0,da1

C       CALCULATE DERIVATE dX**2/dA2
        DO  I=1,NUM
           Y2(I) = A1 * EXP(-(X(I) / (A2 + 0.1 * DA2))**2) + A3
        ENDDO

        X1 = 0
        DO I=1,NUM
           X1 = X1 + (Y1(I) - Y2(I))**2
        ENDDO

        DXA2 = (X1-X0) / (0.1 * DA2)
        

C       CALCULATE DERIVATE dX**2 / dA3
        DO  I=1,NUM
           Y2(I) = A1 * EXP(-(X(I)/A2)**2) + (A3 + 0.1 * DA3)
        ENDDO

        X1 = 0
        DO I=1,NUM
           X1 = X1 + (Y1(I) - Y2(I)) **2
        ENDDO

        IF (NUM < 3) THEN
           DXA3 = 0
        ELSE
           DXA3 = (X1-X0) / (0.1 * DA3)
        ENDIF

c       write(6,*) ' for sum:',dxa1,da1
c       write(6,*) ' for sum:',dxa2,da2
c       write(6,*) ' for sum:',dxa3,da3

        SUM = SQRT((DXA1*DA1)**2 + (DXA2*DA2)**2 + (DXA3*DA3)**2)

c       write(6,*) ' sum:',sum

        A1 = A1 - DXA1 * DA1**2 / SUM
        A2 = A2 - DXA2 * DA2**2 / SUM
        A3 = A3 - DXA3 * DA3**2 / SUM
        
C       CRITERIA FOR ITERATION  ...... CALCULATE THE Y2
        DO  I=1,NUM
           Y2(I) = A1 * EXP(-(X(I )/ A2)**2) + A3
        ENDDO

        X1 = 0
        DO  I=1,NUM
           X1 = X1 + (Y1(I) - Y2(I))**2
        ENDDO
        D = X1 - X0

C       SOLVE GOTO 899 BY REWRITING THE IF () THEN .. LOOP . ML 7/5/95
        IF (D <= 0) THEN
           X0 = X1
           P  = A2 / KS
C          WRITE(NOUT,*)'A1=',A1,' A2=',A2,'(A-1) ',P,'POINTS',' A3',A3
           NSTEP = NSTEP + 1
           GOTO 999
        ENDIF

        A1 = A1 + 0.5 * DXA1 * DA1**2 / SUM
        A2 = A2 + 0.5 * DXA2 * DA2**2 / SUM
        A3 = A3 + 0.5 * DXA3 * DA3**2 / SUM
        P  = A2 / KS



99      IP = P          ! TRUNCATE -----------------------------------

        WRITE(NOUT,*) ' '
        WRITE(NOUT,*) ' LOWER ENVELOPE EQUATION PARAMETERS:'
        WRITE(NOUT,91) A1, A2, A3, IP
91      FORMAT(  '  A1:',     ES11.4,' ',  
     &           '  A2:',     ES11.4,' ',
     &           '  A3:',     ES11.4,' (A-1) ',
     &           '  POINTS:', I0,/)

C       ADD ADDITIONAL VALUE, SO ALL THE MINIMI ARE ABOVE ZERO, 
C       AT THE OTHER PLACES VALUES BELOW ZERO ARE SET TO ZERO.

        XMIN = 0
        DO I=1,NX

           X1    = FLOAT(I) * KS
           Y1(I) = A1 * EXP(-(X1/A2)**2) + A3
           Y2(I) = Y2(I) - Y1(I)

           DO  J=1,NUM
              IF ( I == INT(RMINS(J))-1 .OR. 
     &             I == INT(RMINS(J))     .OR.
     &             I == NINT(RMINS(J))+1 ) 
     &           XMIN = MIN( Y2(I), XMIN )
           ENDDO
        ENDDO


        DO I=1,NX
           Y2(I) = Y2(I) - XMIN
           IF (Y2(I) < 0) Y2(I) = 0

C          PUSH Y2(I) INTO OUTPUT DOC. FILE
           IF (NLUNDOC > 0) 
     &        CALL LUNDOCWRTDAT(NLUNDOC,I,Y2(I),1,IRTFLG)

        ENDDO

        CALL WRTLIN(LUN2,Y2,NX,1)
  
        CALL REG_SET_NSEL(1,1,FLOAT(NUM),0,0,0,0, IRTFLG)
          
9999    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUNDOC)

        END

