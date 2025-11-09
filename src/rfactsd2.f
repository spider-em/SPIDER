
C ++********************************************************************
C                                                                      
C RFACTSD2                                                             
C                  REGISTER OUTPUT ADDED         JAN 2005 ARDEAN LEITH   
C                  VERBOSE                       FEB 2006 ARDEAN LEITH  
C                  DOC FILE HEADER               NOV 2009 ARDEAN LEITH  
C                  DOC FILE HEADER               JUN 2011 ARDEAN LEITH  
C                  MINIMUM WARNING               JUN 2011 ARDEAN LEITH  
C                  FSC                           FEB 2012 ARDEAN LEITH 
C                  FSCCUT                        SEP 2012 ARDEAN LEITH
C                  WANTSQRTS                     MAY 2014 ARDEAN LEITH
C                  REPORTS MIN FSC NOT LAST      MAY 2014 ARDEAN LEITH
C                  RESOLUTION FOR GPL            APR 2016 ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2016  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
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
C RFACTSD2(PR,AMP,CSUM1,LR,CSUM,CSUM2,AVSUM, 
C          NSCALE,INC,WI,FACT,NOUT,TWOD,FSCOP,WANTSQRTS)
C
C PURPOSE: PUTS FRC/FSC STATISTICS IN DOC FILE
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE RFACTSD2(PR,AMP,CSUM1,LR,CSUM,CSUM2,AVSUM,
     &                       NSCALE,INC,WI,FACT,TWOD,
     &                       LUNDOC,FSCOP,FMAXSPFREQ,LUNGP,FSCCUT,
     &                       WANTSQRTS)


         INCLUDE 'CMLIMIT.INC' 
         INCLUDE 'CMBLOCK.INC'
 
         REAL     :: PR(NSCALE,INC),AMP(NSCALE,INC),CSUM1(INC)
         INTEGER  :: LR(INC)
         REAL     :: CSUM(INC),CSUM2(INC),AVSUM(NSCALE,INC)
         REAL     :: RVAL(INC),FSCVAL(INC)
         INTEGER  :: NSCALE,INC
         REAL     :: WI,FACT
         LOGICAL  :: TWOD
         INTEGER  :: LUNDOC
         LOGICAL  :: FSCOP
         REAL     :: FMAXSPFREQ,PIXSIZ
         INTEGER  :: LUNGP,ILOCMIN
         REAL     :: FSCCUT,FSCMIN
         LOGICAL  :: WANTSQRTS

         CHARACTER (LEN=MAXNAM)   :: DOCNAM,GPLOTFILE,PROMPT
         CHARACTER (LEN=2*MAXNAM) :: MSG

         INTEGER  :: ICOMM,MYPID,MPIERR
         LOGICAL  :: IFOUNDIT,NEWFILE

         LOGICAL  :: WANTGPLOT = .FALSE.

         INTEGER, PARAMETER :: NLIST = 6
         REAL     :: DLIST(NLIST)

         CALL SET_MPI(ICOMM,MYPID,MPIERR)

         CALL REG_GET_USED(NSEL_USED)

         XPREV  = 0
         DLIST  = HUGE(FSCLAST)

         WIP    = 1.0 / WI  ! TO PIXELS
         PIXLEN = 0.5 / FMAXSPFREQ

         IF (TWOD) THEN
C           FRC: WRITE FRC RESULT FILE OUTPUT

            CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOC,NICDOC,.TRUE.,
     &                   'FRC OUTPUT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,
     &                   NEWFILE,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            IF (FSCOP) THEN
C              OPEN FORMATTED, SEQUENTIAL FILE FOR GNUPLOT COMMANDS
               CALL OPAUXFILE(.TRUE.,GPLOTFILE,DATEXC,LUNGP,0,'N',
     &                          'GNUPLOT',.TRUE.,IRTFLG)
               NLET      = lnblnkn(GPLOTFILE)
          !write(6,*) 'gplot:',nlet,irtflg,gplotfile

               WANTGPLOT =  (NLET > 0) .AND. (IRTFLG == 0)  
               IF (IRTFLG > 0) RETURN
               MSG = '       NORM-FREQ      ANGST         FRC' 
               NMSG = 44


            ELSE
               MSG = '       NORM-FREQ      DPH           FRC'//
     &               '          FRCCRIT        PIXELS'    
C                             10        20        30        40
C                     1234567890123456789012345678901234567890 
               NMSG = 72  
   
            ENDIF

            IF (VERBOSE .AND. MYPID <= 0) THEN
               WRITE(NOUT,*)' '
               IF (FSCOP) THEN
                  WRITE(NOUT,90) WIP
90                FORMAT('  FRC,   RING WIDTH (PIXELS):',G10.3)
                  WRITE(NOUT,*) ' '

                  WRITE(NOUT,91)
91                FORMAT(10X,'|NORM-FREQ|     |ANGST|         |FSC|') 
               ELSE

                  WRITE(NOUT,92) WI
92                FORMAT('  PHASE RESIDUE & FRC, ',
     &                   '  RING WIDTH (RECIPROCAL UNITS):',G10.3)
                  WRITE(NOUT,*) ' '
                  WRITE(NOUT,93)
93                FORMAT(10X,'|NORM-FREQ|      |DPH|          |FRC|',
     &                       '     |FRCCRIT|     |PIXELS|')
               ENDIF
            ENDIF
         ELSE
C           FSC: WRITE RESULT FILE

            CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOC,NICDOC,.TRUE.,
     &                'FSC OUTPUT DOCUMENT',.FALSE.,.FALSE.,.TRUE.,
     &                NEWFILE,IRTFLG)
            IF (IRTFLG .NE. 0) RETURN

            IF (FSCOP) THEN
C              OPEN FORMATTED, SEQUENTIAL FILE FOR GNUPLOT COMMANDS
               CALL OPAUXFILE(.TRUE.,GPLOTFILE,DATEXC,LUNGP,0,'N',
     &                          'GNUPLOT',.TRUE.,IRTFLG)
               NLET      = lnblnkn(GPLOTFILE)
               WANTGPLOT =  (NLET > 0) .AND. (IRTFLG == 0)  
               IF (IRTFLG > 0) RETURN
C                     12345678901234567890123456789012345678901234567890  
               MSG = '       NORM-FREQ      ANGST         FSC' 
               NMSG = 44  

               IF (WANTSQRTS) THEN
                  MSG = '       NORM-FREQ      ANGST         FSC'//
     &                  '      |SQRT(FSC)| |SQRT(2FSC/(FSC+1))|'    
C                                10        20        30        40
C                        1234567890123456789012345678901234567890 
                  NMSG = 79
               ENDIF

            ELSE
               MSG = '       NORM-FREQ      DPH           FSC' //
     &               '          FSCCRIT        VOXELS'    
C                            10        20        30        40
C                    1234567890123456789012345678901234567890 
               NMSG = 72  
            ENDIF

            IF (VERBOSE .AND. MYPID <= 0) THEN
               WRITE(NOUT,*)' '
               IF (FSCOP) THEN
                  WRITE(NOUT,94) WIP
94                FORMAT('  FSC,  SHELL WIDTH (VOXELS):',G10.3)
                  WRITE(NOUT,*) ' '
                  WRITE(NOUT,95)
95                FORMAT(10X,'|NORM-FREQ|     |ANGST|         |FSC|')      
               ELSE
                  WRITE(NOUT,96) WI
96                FORMAT('  PHASE RESIDUE & FSC,',
     &                   '   SHELL WIDTH (RECIPROCAL UNITS):',G10.3)
                  WRITE(NOUT,*) ' '
                  WRITE(NOUT,97)
97                FORMAT(10X,'|NORM-FREQ|      |DPH|          |FSC|',
     &                       '     |FSCCRIT|     |VOXELS|')
                ENDIF
            ENDIF
         ENDIF

         CALL LUNDOCPUTCOM(LUNDOC,MSG(:NMSG),IRTFLG)

         DO  L=1,INC
            NVOX = LR(L)

            IF (NVOX .NE. 0) THEN
               DLIST(1) = L         ! NUMBER

               FSCMIN   = HUGE(FSCMIN)      ! MIN OF FSC CURVE
               SPFLAST  = DLIST(2)  ! WAS HUGE
               DLIST(2) = FLOAT(L-1) / FLOAT(INC-1)*0.5  ! NORM FREQ

               IF (FSCOP) THEN
                  DLIST(3) = 0.5 / (FMAXSPFREQ * DLIST(2)) ! RESOL
                  RVAL(L)  = DLIST(3) 
               ENDIF
               DLIST(5) = MIN(1.0, FACT/SQRT(FLOAT(NVOX))) ! FSCCRIT
               DLIST(6) = NVOX                           ! # VOXELS

               RFMIN    = -HUGE(RFMIN)
               NSCM     = 1
               IFOUNDIT = .FALSE.

               DO NSC=1,NSCALE
		 IF (AMP(NSC,L) > TINY(RFMIN))  THEN
                    RFM = AVSUM(NSC,L) / AMAX1(1.0,AMP(NSC,L))
                    IF (RFM  <  RFMIN) THEN
                       NSCM     = NSC
                       RFMIN    = RFM
                       IFOUNDIT = .TRUE.
                    ENDIF
		 ENDIF
               ENDDO

C              NSCM IS THE NUMBER OF THE ELEMENT IN EACH ARRAY WITH THE
C              CORRECT SCALING. SCALE IS THE CORRECT SCALING.

               BK1 = AMP(NSCM,L)
               BK2 = PR(NSCM,L)
               IF (BK1 > TINY(BK3)) THEN
		  IF (.NOT. FSCOP) DLIST(3) = SQRT(BK2/BK1)
               ELSE
		  IF (.NOT. FSCOP) DLIST(3) = 0.0
               ENDIF
               BK3     = CSUM2(L)
               BK4     = CSUM1(L)

               FSCLAST = DLIST(4)  ! PREVIOUS FSC (STARTS AS HUGE)

               IF (BK3 > TINY(BK3) .AND. 
     &             BK4 > TINY(BK3)) THEN
		  DLIST(4) = CSUM(L) / SQRT(BK4 * BK3)   ! FSC
               ELSE
		  DLIST(4) = 0.0
               ENDIF
               IF (FSCOP) THEN
                  FSCVAL(L)  = DLIST(4) 
               ENDIF

               IF (WANTSQRTS) THEN
C                 WRITE FSC FILTER REGISTER COLUMNS
                  FSCZ   = MAX(DLIST(4), 0.0)
                  IF (FSCZ == 0.0) THEN
                     DLIST(5) = 0.0
		     DLIST(6) = 0.0
                  ELSE
		     DLIST(5) = SQRT(FSCZ)
		     DLIST(6) = SQRT( 2*FSCZ/ (FSCZ+1) )
                  ENDIF
               ENDIF


C              WRITE TO DOC FILE
               NLEN = 5
               IF (FSCOP) NLEN = 3
               IF (FSCOP .AND. WANTSQRTS) NLEN = 5
               CALL LUNDOCWRTDAT(LUNDOC,L,DLIST(2),NLEN,IRTFLG)

               IF (VERBOSE .AND. MYPID <= 0) THEN
                  IF (FSCOP) THEN
                     WRITE(NOUT,98) L,(DLIST(K),K=2,4)
                  ELSEIF (IFOUNDIT) THEN
                     WRITE(NOUT,98) L,(DLIST(K),K=2,5),NVOX
98                   FORMAT (1X,I4,4(2X,F12.5),4X,I6)
                  ELSE
                     WRITE(NOUT,99) L,(DLIST(K),K=2,5),NVOX
99                   FORMAT (1X,I4,4(2X,F12.5),4X,I6,'  LACKS MINIMUM!')
                  ENDIF
               ENDIF

               IF (DLIST(4) < FSCMIN) THEN
C                 MIN VALUE ON FSC CURVE
                  FSCMIN = DLIST(4)
                  ILOCMIN = L
               ENDIF

               IF (L  >=  3 .AND. 
     &            FSCLAST  >= FSCCUT .AND.
     &            DLIST(4) <  FSCCUT) THEN

C                 CROSSED FSCCUT GOING DOWN
                  XPREV   = L - 1     ! LAST INDEX ABOVE CUTOFF

                  FSCPREV = FSCLAST
                  FSCNOW  = DLIST(4)

                  SPFPREV = SPFLAST
                  SPFNOW  = DLIST(2)
               ENDIF
            ENDIF
         ENDDO

         IF (VERBOSE .AND. MYPID <= 0) WRITE(NOUT,*) ' '

         IF (FSCOP .OR. NSEL_USED > 0) THEN
C           RESOLUTION NEEDED
            IF (XPREV > 0) THEN
               FINTERP   = (FSCCUT - FSCPREV) / (FSCNOW - FSCPREV)  

               XINTERP   = XPREV   + FINTERP * (1)
               SPFINTERP = SPFPREV + FINTERP * (SPFNOW - SPFPREV)          
            ELSE
               XINTERP   = ILOCMIN
               SPFINTERP = FSCMIN    ! NOT INTERPOLATED?
            ENDIF

            RESOL = 0.0
            IF (FSCOP) THEN
               RESOL = 0.5 / (FMAXSPFREQ * SPFINTERP )
               !write(6,*) 'resol:',RESOL,FMAXSPFREQ,SPFINTERP,xinterp 
            ENDIF
         ENDIF

         IF (NSEL_USED > 0) THEN
C           OUTPUT TO SPIDER'S REGISTERS NEEDED
            CALL REG_SET_NSEL(1,3,XINTERP,SPFINTERP,RESOL,
     &                            0.0,0.0,IRTFLG)
         ENDIF

         IF (WANTGPLOT) THEN
C           WRITE GNUPLOT FILE OUTPUT

            WRITE(LUNGP,'(A)') 'set xlabel "Angstroms"' 
            WRITE(LUNGP,191)   'set title " At FSC: 0.5  Resolution:', 
     &                          RESOL,' Angstroms"'
191         FORMAT(A, F7.2, A)
            WRITE(LUNGP,'(A)') 'set yrange [0:1.0]' 
            NPIXLEND2 = NINT(INC * 0.5)

            WRITE(LUNGP,192) 'set xrange [0:',NPIXLEND2,'] reverse' 
192         FORMAT(A, I6, A)

            !WRITE(LUNGP,'(A)''SET TERM POSTSCRIPT'
            !WRITE(LUNGP,'(A)''save "outfile"'

             WRITE(LUNGP,'(A)') 'plot 0.5 ,  "-" using 2:3 with line'

            DO  L=1,INC
               IF (LR(L) .NE. 0) THEN
                  WRITE(LUNGP,195) L,RVAL(L),FSCVAL(L)
195               FORMAT(' ',I5,'  ',F7.2,' ',F10.3)
               ENDIF
           ENDDO

         ENDIF

         END
         !write(6,*) 'float(l-1):',FLOAT(L-1)
         !write(6,*) 'float(inc-1):',FLOAT(inc-1)
         !write(6,*) 'normfreq:',normfreq
