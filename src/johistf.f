
C++*********************************************************************
C
C JOHISTF.F                 NEW APRIL  2003  ArDean Leith
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
C    JOHISTF(LUN1,LUN2,NSAM,NROW,NSLICE,NBINSA,NBINSP,IRTFLG)
C
C    PURPOSE:    MUTUAL SHARED INFO FOR 2 FOURIER INPUT IMAGES.
C
C    PARAMETERS:  LUN1           IO UNIT NUMBER OF IMAGE FILE
C                 LUN2           IO UNIT NUMBER OF IMAGE FILE
C                 NSAM,NROW      DIMENSIONS OF IMAGE              
C                 NSLICE         DIMENSIONS OF IMAGE
C                 NBINSA,NBINSP  BINS IN HISTOGRAM                (SENT)
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE JOHISTF(LUN1,LUN2,NSAM,NROW,NSLICE,
     &                   NBINSA,NBINSP,IRTFLG)

      INCLUDE 'CMBLOCK.INC'

#ifndef SP_32
      INTEGER *8 NPIX
      INTEGER *8 FREQA1(NBINSA),FREQA2(NBINSA)
      INTEGER *8 FREQP1(NBINSP),FREQP2(NBINSP)

      INTEGER* 8 FREQA(NBINSA,NBINSp) 
      INTEGER* 8 FREQP(NBINSa,NBINSP)
#else
      INTEGER *4 NPIX
      INTEGER *4 FREQA1(NBINSA),FREQA2(NBINSA)
      INTEGER *4 FREQP1(NBINSP),FREQP2(NBINSP)

      INTEGER* 4 FREQA(NBINSA,NBINSp) 
      INTEGER* 4 FREQP(NBINSa,NBINSP)
#endif

      REAL, DIMENSION(NSAM*NROW/2,2) :: BUF1,BUF2

      DOUBLE PRECISION :: FNPIX

C     GET IMAGE VALUES
      CALL REDVOL(LUN1,NSAM,NROW,1,NSLICE,BUF1,IRTFLG)
      CALL REDVOL(LUN2,NSAM,NROW,1,NSLICE,BUF2,IRTFLG)

      NPIX   = NSAM * NROW / 2 * NSLICE
      FNPIX  = 1.0 / FLOAT(NPIX)
      ANGF   = 180.0 / 3.14159

      FMINA1 = HUGE(FMINA1)
      FMAXA1 = -FMINA1

      FMINA2 = FMINA1
      FMAXA2 = FMAXA1
      FMINP1 = FMINA1
      FMAXP1 = FMAXA1
      FMINP2 = FMINA1
      FMAXP2 = FMAXA1

C     GET AMPLITUDE MIN/MAX FROM IMAGE VALUES
      DO  IPIX=1,NPIX

C        FIND AMPLITUDE & PHASE OF FIRST FILE
         BR   = BUF1(IPIX,1)
         BI   = BUF1(IPIX,2)
         AM   = SQRT(BI**2 + BR**2)
         PH   = 0.0
         IF (BI .NE. 0. .OR. BR .NE. 0.) PH = ATAN2(BI,BR)*ANGF

         FMINA1  = MIN(FMINA1,AM)
         FMAXA1  = MAX(FMAXA1,AM)

         FMINP1  = MIN(FMINP1,PH)
         FMAXP1  = MAX(FMAXP1,PH)

         BUF1(IPIX,1) = AM
         BUF1(IPIX,2) = PH

C        FIND AMPLITUDE & PHASE OF SECOND FILE
         BR   = BUF2(IPIX,1)
         BI   = BUF2(IPIX,2)
         AM   = SQRT(BI**2 + BR**2)
         PH   = 0.0
         IF (BI .NE. 0. .OR. BR .NE. 0.) PH = ATAN2(BI,BR)*ANGF

         FMINA2  = MIN(FMINA2,AM)
         FMAXA2  = MAX(FMAXA2,AM)

         FMINP2  = MIN(FMINP2,PH)
         FMAXP2  = MAX(FMAXP2,PH)

         BUF2(IPIX,1) = AM
         BUF2(IPIX,2) = PH
      ENDDO  

      WRITE(NOUT,*) ' '
      WRITE(NOUT,90) NPIX,FMINA1,FMAXA1,FMINP1,FMAXP1,
     &                    FMINA2,FMAXA2,FMINP2,FMAXP2

90    FORMAT(
     & ' PIXELS:            ',I11,/,   
     & ' FIRST FILE AMP. RANGE:  ',1PG11.4,'   .........    ',1PG11.4,/,
     & ' FIRST FILE PHASE RANGE: ',1PG11.4,'   .........    ',1PG11.4,/,
     & ' 2ND   FILE AMP. RANGE:  ',1PG11.4,'   .........    ',1PG11.4,/,
     & ' 2ND   FILE PHASE RANGE: ',1PG11.4,'   .........    ',1PG11.4)


C     ZERO THE HISTOGRAM FREQUENCIES
      FREQA  = 0 
      FREQP  = 0 

      FREQA1 = 0 
      FREQA2 = 0
 
      FREQP1 = 0 
      FREQP2 = 0 

C     FIND HISTOGRAM BIN SELECTION FACTORS
      HDIFFA1 = FMAXA1 - FMINA1
      FFA1    = (NBINSA - 1.0) / HDIFFA1

      HDIFFA2 = FMAXA2 - FMINA2
      FFA2    = (NBINSA - 1.0) / HDIFFA2

      HDIFFP1 = FMAXP1 - FMINP1
      FFP1    = (NBINSP - 1.0) / HDIFFP1

      HDIFFP2 = FMAXP2 - FMINP2
      FFP2    = (NBINSP - 1.0) / HDIFFP2

C     GET HISTOGRAMS FROM IMAGE VALUES
      DO  IPIX=1,NPIX

C        FIND BIN1 NUMBER
         BVAL   = BUF1(IPIX,1)
         JBINA1 = INT((BVAL - FMINA1) * FFA1) + 1.5

C        FIND BIN2 NUMBER
         BVAL   = BUF2(IPIX,1)
         JBINA2 = INT((BVAL - FMINA2) * FFA2) + 1.5

         IF (JBINA1.GE.1 .AND. JBINA1.LE.NBINSA  .AND.
     &       JBINA2.GE.1 .AND. JBINA2.LE.NBINSA) THEN
C           WITHIN HISTOGRAM RANGE
            FREQA1(JBINA1)       = FREQA1(JBINA1) + 1
            FREQA2(JBINA2)       = FREQA2(JBINA2) + 1
            FREQA(JBINA1,JBINA2) = FREQA(JBINA1,JBINA2) + 1
         ENDIF
 
C        FOR PHASES
C        FIND BIN1 NUMBER
         BVAL   = BUF1(IPIX,2)
         JBINP1 = INT((BVAL - FMINP1) * FFP1) + 1.5

C        FIND BIN2 NUMBER
         BVAL   = BUF2(IPIX,2)
         JBINP2 = INT((BVAL - FMINP2) * FFP2) + 1.5

         IF (JBINP1.GE.1 .AND. JBINP1.LE.NBINSP  .AND.
     &       JBINP2.GE.1 .AND. JBINP2.LE.NBINSP) THEN
C           WITHIN HISTOGRAM RANGE
            FREQP1(JBINP1)       = FREQP1(JBINP1) + 1
            FREQP2(JBINP2)       = FREQP2(JBINP2) + 1
            FREQP(JBINP1,JBINP2) = FREQP(JBINP1,JBINP2) + 1
         ENDIF
      ENDDO

C     FIND ENTROPY
      ENTROPYAA = 0.0
      ENTROPYPP = 0.0

      DO IBIN = 1,NBINSA
         DO JBIN=1,NBINSP
            FTAA = FREQA(IBIN,JBIN) * FNPIX
            IF (FTAA .GT. 0.0) ENTROPYAA = ENTROPYAA - FTAA * LOG(FTAA)

            FTPP = FREQP(IBIN,JBIN) * FNPIX
            IF (FTPP .GT. 0.0) ENTROPYPP = ENTROPYPP - FTPP * LOG(FTPP)
         ENDDO 
      ENDDO

      ENTROPYA1 = 0.0
      DO IBIN = 1,NBINSA
         FTA1 = FREQA1(IBIN) * FNPIX
         IF (FTA1 .GT. 0.0) ENTROPYA1 = ENTROPYA1 - FTA1 * LOG(FTA1) 
      ENDDO


      ENTROPYA2 = 0.0
      DO IBIN = 1,NBINSA
         FTA2 = FREQA2(IBIN) * FNPIX
         IF (FTA2 .GT. 0.0) ENTROPYA2 = ENTROPYA2 - FTA2 * LOG(FTA2) 
      ENDDO

      ENTROPYP1 = 0.0
      DO IBIN = 1,NBINSP
         FTP1 = FREQP1(IBIN) * FNPIX
         IF (FTP1 .GT. 0.0) ENTROPYP1 = ENTROPYP1 - FTP1 * LOG(FTP1) 
      ENDDO

      ENTROPYP2 = 0.0
      DO IBIN = 1,NBINSP
         FTP2 = FREQP2(IBIN) * FNPIX
         IF (FTP2 .GT. 0.0) ENTROPYP2 = ENTROPYP2 - FTP2 * LOG(FTP2) 
      ENDDO

      FMSIA = ENTROPYA1 + ENTROPYA2 - ENTROPYAA

      FMSIP = ENTROPYP1 + ENTROPYP2 - ENTROPYPP

      FMSI = FMSIA + FMSIP

      WRITE(NOUT,91) NBINSA,NBINSP,
     &               ENTROPYA1,ENTROPYA2,ENTROPYAA,
     &               ENTROPYP1,ENTROPYP2,ENTROPYPP,
     &               FMSIA,FMSIP,FMSI

91    FORMAT(
     &    ' AMPLITUDE BINS:      ',I11,    '  PHASE BINS:  ',I11,//,
     &    ' FIRST FILE AMP. ENTROPY:      ',1PG11.4,/,
     &    ' 2ND   FILE AMP ENTROPY:       ',1PG11.4,/,
     &    ' JOINT AMPLITUDE ENTROPY:      ',1PG11.4,/,
     &    ' FIRST FILE PHASE ENTROPY:     ',1PG11.4,/,
     &    ' 2ND   FILE PHASE ENTROPY:     ',1PG11.4,/,
     &    ' JOINT PHASE ENTROPY:          ',1PG11.4,/,
     &    ' AMPLITUDE MUTUAL SHARED INFO: ',1PG11.4,/,
     &    ' PHASE MUTUAL SHARED INFO:     ',1PG11.4,/,
     &    ' JOINT MUTUAL SHARED INFO:     ',1PG11.4/)

      CALL REG_GET_USED(NSEL_USED)

      IF (NSEL_USED .GT. 0) THEN
        CALL REG_SET_NSEL(1,3,FMSIA,FMSIP,FMSI, 0.0, 0.0,IRTFLG)
      ENDIF

      RETURN
      END
  
