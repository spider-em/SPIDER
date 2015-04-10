
#ifdef SP_SUN4 
 
C   THIS ROUTINE FAILED TO COMPILE ON SUN
 
       SUBROUTINE tvwn3(MAXDIM)
 
       COMMON /UNITS/LUNC,NIN,NOUT
 
       WRITE(NOUT,*) 'DUMMY CALL: tvwn3'
       RETURN
       END
 
#else
C++*********************************************************************
C
C TVWN3.F                           ADAPTED FROM TVWN2.F  -- SEPT 92 al
C                                   USED OPFILE NOV 00 ARDEAN LEITH
C                                   MAXKEY SEPT 02 ARDEAN LEITH
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
C TVWN3()
C
C PURPOSE:     TV WINDOW WITHOUT THE TV.  RECIPROCAL LATTICE SEARCH 
C              AND REFINEMENT
C       
C PARAMETERS:  NONE
C
C CALLER:      DRIV2  --> TVWN# --> LATCEN --> SOLV2D
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

        SUBROUTINE TVWN3(MAXDIM)

        INCLUDE 'CMBLOCK.INC'
        INCLUDE 'CMLIMIT.INC'

        COMMON /IOBUF/ BUFF(NBUFSIZ)

C       PEAK ARRAY SIZE
        PARAMETER (NREFLM = 200)
	COMPLEX   C
	INTEGER   HIND,KIND
	REAL      FINT,KX,KY
	COMMON    C(NREFLM),HIND(NREFLM),KIND(NREFLM),
     &            FINT(NREFLM),KX(NREFLM),KY(NREFLM),AREA(1)
C       WARNING:  MAXF ALSO USES THIS UNLABELED COMMON

        PARAMETER (MAXKEY=9999) 
        PARAMETER (MAXREG = 7)
        COMMON /DOC_BUF/ DBUF(MAXREG,MAXKEY*2)

        CHARACTER(LEN=MAXNAM)  ::  FILNAM
	COMPLEX        CC
        DIMENSION      PLIST(7)

        DATA LUNPOW,LUNFOU,LUNFILT,LUNDOC/11,12,13,14/
	DATA PI/3.14159/

        NLIST  = 7
        WGT    = 1.0
        CX     = 0.0
        CY     = 0.0

	IDM    = 0
	IF (FCHAR(5:5) .EQ. 'D') IDM = 1

	WRITE(NDAT,47)
47      FORMAT(' --- WT : AUTOMATIC LATTICE INDEXING. VERSION 11/6/96 ',
     &         '---      WEIGHTING USED: AMP'/) 

C       OPEN POWER SPECTRUM FILE -----------------------------------
48      MAXIM = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNPOW,'O',IFORM,
     &        NSAM1,NROW1,NSLICE1,MAXIM,'POWER SPECTRUM',.FALSE.,IRTFLG)
        IF (IRTFLG .NE. 0) RETURN

        WRITE(NDAT,49)
49      FORMAT(' xxxxxxxxxxxxxxxx  POWER SPECTRUM xxxxxxxxxxx'/)

C       OPEN FOURIER FILE -------------------------------------------
        WRITE(NOUT,*) ' '
        MAXIM = 0
	CALL OPFILEC(0,.TRUE.,FILNAM,LUNFOU,'O',IFORM,
     &      NSAM3,NROW3,NSLICE3, MAXIM,'FOURIER',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 48

        IF (IFORM .GE. 0) THEN
C          NOT A FOURIER FILE
           CALL ERRT(2,'TVWN3',NE)
           GOTO 918
        ENDIF

        WRITE(NDAT,61)
61      FORMAT(' xxxxxxxxxxxx FOURIER TRANSFORM xxxxxxxxxxxxx'/)
        WRITE(NDAT,62)
62      FORMAT(6X,'  AMPL     PHASE   ENERGY      WEIGHT'/)
	WRITE(NDAT,9132)
9132    FORMAT(6X,' REFL.NO.     H        K     KX COMP.   KY COMP.')

        CRIT = 1.0
        R    = 100.0

C       COPY F00 FOR USE IN FILTER FILE
        CALL REDLIN(LUNFOU,BUFF,NSAM3,1)
	F00 = BUFF(1)

        WRITE(NOUT,*) ' '

C       OPEN REFLECTION CONTAINING DOC FILE --------------------------
306     CALL FILERD(FILNAM,NLET,DATEXC,'REFLECTIONS DOC',IRTFLG)

        IKEY  = 1
        ICALL = 0
        CALL UNSDAL(FILNAM,ICALL,LUNDOC,IKEY,PLIST,NLIST,
     &              DBUF,NREFLM,MAXREG,NKEY,LERR)
        IF (LERR .EQ. 2) GOTO 306
             
        NTOT = 0
        DO I = 1,NKEY
           IF (DBUF(1,I) .GT. 0) THEN
C             FOUND A VALID KEY FOR REFLECTION
              NTOT = NTOT + 1
              IF (NTOT .GT. NREFLM) THEN
                 CALL ERRT(6,'TVWN3',NE)
                 GOTO 916
              ENDIF
              HIND(NTOT) = DBUF(3,I)
              KIND(NTOT) = DBUF(4,I)
              KX(NTOT)   = DBUF(5,I)
              KY(NTOT)   = DBUF(6,I)

              WRITE(NOUT,*) NTOT,HIND(NTOT),KIND(NTOT),KX(NTOT),KY(NTOT)
           ENDIF
        ENDDO


        CALL RDPRMI(IXB,IXB1,NOT_USED,'WINDOW SIZES (SEARCH,FILT)')
	IF (IXB1.EQ.0) IXB1 = IXB
	IYB    = IXB
        NSAMD2 = NSAM1 / 2
        NROWD2 = NROW1 / 2

C       REFINE EACH OF THE REFLECTIONS

        WRITE(NOUT,*) ' '
        WRITE(NOUT,*) 
     &     ' COORDINATES:   INPUT      TRANSFORMED       EXACT'

        NREFL = 0

500     NREFL  = NREFL + 1
C        KXPOS = KY(NREFL) - NROWD2
C        KYPOS = KX(NREFL) - NSAMD2
         KXPOS = KX(NREFL) - NSAMD2
         KYPOS = KY(NREFL) - NROWD2

C        RECORD COORDINATES IN INTEGER FORM FOR OUTPUT
         KXI   = KX(NREFL)
         KYI   = KY(NREFL)
         KXT   = KXPOS
         KYT   = KYPOS

         CALL MAXF(LUNPOW,LUNFOU,NSAM1,NROW1,NSAM3,NROW3,
     &       IXB,IYB,KXPOS,KYPOS,C(NREFL),FINT(NREFL),IDM,R,CRIT)

         WRITE(NOUT,9001) KXI,KYI,KXT,KYT,KXPOS,KYPOS
9001     FORMAT(13X,3('(',I5,',',I5,')',3X))

         KX(NREFL) = KXPOS
         KY(NREFL) = KYPOS
         WGT = CABS(C(NREFL))

         IF (FINT(NREFL) .GE. 0.0) THEN
C          VALID REFLECTION

           IF (NREFL .EQ. 1) CALL LATCEN(IERR,0,IDUM,IDUM,
     &                    FDUM,FDUM,FDUM,FDUM,FDUM,FDUM,FDUM,FDUM,FDUM)

           CALL LATCEN(IERR,1,HIND(NREFL),KIND(NREFL),KX(NREFL),
     &                 KY(NREFL),WGT,CX,CY,FDUM,FDUM,FDUM,FDUM)
           CC  = C(NREFL)
           AMP = CABS(CC)
           PH  = ATAN2(AIMAG(CC),REAL(CC))*180./PI
           WRITE(NDAT,9111)NREFL,HIND(NREFL),KIND(NREFL),KX(NREFL),
     &                     KY(NREFL),AMP,PH,FINT(NREFL),WGT

         ELSEIF (NREFL .LT. NTOT) THEN
C          DISCARD THIS REFLECTION FROM FURTHER CONSIDERATION
           DO I = NREFL,NTOT-1
             KX(I)   = KX(I+1)
             KY(I)   = KY(I+1)
             HIND(I) = HIND(I+1)
             KIND(I) = KIND(I+1)
             FINT(I) = FINT(I+1)
             C(I)    = C(I+1)
             NTOT    = NTOT - 1
           ENDDO
           NREFL = NREFL - 1

         ELSE
C          DISCARD LAST REFLECTION OF THE SERIES
           NTOT  = NTOT - 1
           NREFL = NTOT
         ENDIF

C        LOOP THROUGH ALL REFLECTIONS
         IF (NREFL .LT. NTOT) GOTO 500


        CALL LATCEN(IERR,2,NDUM,NDUM,DUM,DUM,WGT,CX,CY,AX,AY,BX,BY)
	IF (IERR .NE. 0) GOTO 916

C       GENERATE ALL LATTICE POSITIONS WITHIN RESOLUTION RANGE
        CALL RDPRMI(NRESOL,NDUM,NOT_USED,
     &     'MAXIMUM RESOLUTION DISTANCE') 
        IF (NRESOL .EQ. 0) GOTO 918
        RES2 = FLOAT(NRESOL)**2

C       DETERMINE LIMITING INDEX POSITIONS
        NRUN = 1
        NACC = NREFL
C       NRUN = 1 IS THE INTERACTIVE SEARCH, NRUN = 2 THE REFINEMENT RUN

C       ***********************************************************
C       CALCULATE AND PRINT LENGTHS AND ANGLES OF UNIT VECTORS
C       ***********************************************************

9121    DISTA2 = AX**2 + AY**2
        DISTB2 = BX**2 + BY**2
        DISTA  = SQRT(DISTA2)
        DISTB  = SQRT(DISTB2)
        IF (AY .NE. 0.) PHIA = ATAN2(AX,AY)*180./PI
        IF (AY .EQ. 0.) PHIA = SIGN(90.,AX)
        IF (BY .NE. 0.) PHIB = ATAN2(BX,BY)*180./PI
        IF (BY .EQ. 0.) PHIB = SIGN(90.,BX)

        WRITE(NOUT,912) NACC
912     FORMAT(' UNIT VECTORS BASED ON ',I3,' REFLECTIONS'//12X,
     &         'LENGTH     ANGLE     KX COMP.   KY COMP.'/)

        IF (NDAT .NE. NOUT) WRITE(NDAT,912) NACC
        WRITE(NOUT,913) DISTA,PHIA,AX,AY
913     FORMAT(8X,4F10.2)

        IF (NDAT .NE. NOUT) WRITE(NDAT,913) DISTA,PHIA,AX,AY
        WRITE(NOUT,913) DISTB,PHIB,BX,BY
        IF (NDAT .NE. NOUT) WRITE(NDAT,913) DISTB,PHIB,BX,BY
        DPHI = PHIB - PHIA
        WRITE(NOUT,9131) CX,CY,DPHI
        IF (NDAT .NE. NOUT) WRITE(NDAT,9131) CX,CY,DPHI
9131    FORMAT('0ORIGIN AT (',F8.2,',',F8.2,')   DIFF.ANGLE ',F8.2)
        KA1  = 0
        IF (DISTA .EQ. 0.0 .OR. DISTB .EQ. 0.0) THEN
           WRITE(NOUT,*) ' *** ERROR; DIVISION BY ZERO IN TVWN3'
           WRITE(NOUT,*) ' DISTA,DISTB: ',DISTA,DISTB
           GOTO 916
        ENDIF

        KA2   = FLOAT(NRESOL)/DISTA
        KB1L  = -FLOAT(NRESOL)/DISTB
        KB2   = -KB1L
        NREFL = 1

C       ***********************************************************
C       NOW CREATE ALL LATTICE POSITIONS WITHIN RESOLUTION RANGE
        IF (NRUN .EQ. 1) WRITE(NDAT,9051)
9051    FORMAT('0EXTRAPOLATED RECIPROCAL LATTICE'/)
        IF (NRUN .EQ. 2) WRITE(NDAT,9052)
9052    FORMAT('0REFINED RECIPROCAL LATTICE'/)
C       ***********************************************************

        WRITE(NDAT,9132)
        DO KA = KA1,KA2
          KB1 = KB1L
          DO KB = KB1,KB2
            HIND(NREFL) = KA
            KIND(NREFL) = KB
            KX(NREFL)   = KA*AX + KB*BX
            KY(NREFL)   = KA*AY + KB*BY
            RES = KX(NREFL)**2+KY(NREFL)**2
            IF (RES.GT.RES2 .OR. (KA.EQ.0.AND.KB.EQ.0))  CYCLE

            WRITE(NDAT,9111)NREFL,KA,KB,KX(NREFL),KY(NREFL)
            IF (NRUN.EQ.2 .AND. FINT(NREFL) .LT. 0.0) WRITE(NDAT,9098)
9098        FORMAT('+',57X,'*')
9111        FORMAT(1X,3I10,2F10.2,7X,G12.4,F8.2,2G12.4)
            NREFL = NREFL+1
            IF (NREFL .GT. NREFLM) THEN
               CALL ERRT(6,'TVWN3',NE)
               GOTO 916
            ENDIF
9105      ENDDO
        ENDDO


        IF (NRUN .EQ. 2) WRITE(NDAT,9097)
9097    FORMAT(' REFLECTIONS MARKED WITH "*" WERE OMITTED IN THE',
     1         ' REFINEMENT RUN'/)
        NREFL = NREFL-1

        IF (NRUN .NE. 1) GOTO 9145

C       CONTINUE HERE FOR REFINEMENT (NRUN=1); SKIP THIS SECTION
C       AFTER REFINEMENT FINISHED
        CALL LATCEN(IERR,0,IDUM,IDUM,
     &            FDUM,FDUM,FDUM,FDUM,FDUM,FDUM,FDUM,FDUM,FDUM)

C       *********************************************************
C       REDEFINE WINDOW SIZES IXB,IYB IF NECESSARY
C       *********************************************************

        NRUN = 2
        WRITE(NDAT,9133)
9133    FORMAT('0REFLECTION MAXIMA FOUND'/)
        WRITE(NDAT,9132)
        WRITE(NDAT,62)

C       POWER SPECTRUM CENTERED AT (NSAM1/2+1,NROW1/2)

C       *********************************************************
C       HERE NEW TOP LEFT COOS OF WINDOW ARE CALCULATED
C        *********************************************************

        NACC = 0
        DO I = 1,NREFL
	  KXPOS = KX(I)+0.5
	  IF (KX(I).LT.0.)KXPOS = KX(I)-0.5
	  KYPOS = KY(I)+0.5
	  IF (KY(I) .LT. 0.0) KYPOS = KY(I)-0.5
          CALL MAXF(LUNPOW,LUNFOU,NSAM1,NROW1,NSAM3,NROW3,
     &          IXB,IYB,KXPOS,KYPOS,C(I),FINT(I),IDM,R,CRIT)
	  KX(I) = KXPOS
	  KY(I) = KYPOS
          AMP   = CABS(C(I))
          WGT   = AMP
          PH    = ATAN2(AIMAG(C(I)),REAL(C(I)))
          WRITE(NDAT,9111)I,HIND(I),KIND(I),KX(I),
     &                    KY(I),AMP,PH,FINT(I),WGT
          CALL LATCEN(IERR,1,HIND(I),KIND(I),KX(I),KY(I),
     1                WGT,CX,CY,FDUM,FDUM,FDUM,FDUM)
        ENDDO

        CALL LATCEN(IERR,2,NDUM,NDUM,DUM,DUM,WGT,CX,CY,AX,AY,BX,BY)
	IF (IERR .NE. 0) GOTO 918

C       GO TO CALCULATION AND PRINTOUT OF LENGTHS AND ANGLES
        GOTO 9121


9145    IF (NRESOL .EQ. 0) GOTO 920

C       *********************************************************
C       NOW PUT REFLECTIONS INTO DATA FILE (FILTER FILE)
C       *********************************************************

        IFORM  = 1
        IDUMMY = 1
        MAXIM  = 0
        CALL OPFILEC(0,.TRUE.,FILNAM,LUNFILT,'U',IFORM,20,NREFL,IDUMMY,
     &                   MAXIM,'FILTER',.TRUE.,IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 916

        WRITE(NDAT,9146)
9146    FORMAT(' ',15X,' FILTER FILE'/)

C       IACC COUNTS ACCEPTED REFLECTIONS

        DO I = 1,NREFL

C         THIS CALL TO MAXF GETS CORRECTLY ASSIGNED FOU AMPLITUDES
C         (THEY MAY BE OUT OF STEP IN THE REFLECTION COUNTING BECAUSE
C         THE REFINEMENT MAY THROW SOME REFLECTIONS OUT OR ENTER NEW ONES)

          KXPOS = KX(I)+0.5
          IF (KX(I) .LT. 0.0) KXPOS=KX(I)-0.5
          KYPOS = KY(I) + 0.5
          IF (KY(I) .LT. 0.0) KYPOS=KY(I)-0.5
          CALL MAXF(LUNPOW,LUNFOU,NSAM1,NROW1,NSAM3,NROW3,
     1       IXB1,IXB1,KXPOS,KYPOS,C(I),FINT(I),IDM,R,CRIT)
          BUFF(1)  = I
          BUFF(4)  = HIND(I)
          BUFF(5)  = KIND(I)
          BUFF(2)  = KX(I)
          BUFF(3)  = KY(I)
9156      BUFF(6)  = 1.0
          BUFF(7)  = 0.0
          BUFF(8)  = REAL(C(I))
          BUFF(9)  = AIMAG(C(I))
          BUFF(10) = CABS(C(I))
          BUFF(11) = F00
          BUFF(12) = DISTA
          BUFF(13) = DISTB
          BUFF(14) = DPHI
          CALL WRTLIN(LUNFILT,BUFF,20,I)
        ENDDO

916     CLOSE(LUNDOC)
917     CLOSE(LUNFILT)
918     CLOSE(LUNFOU)
919     CLOSE(LUNPOW)

920     RETURN
        END
#endif
