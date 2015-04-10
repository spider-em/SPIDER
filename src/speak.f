
C++*********************************************************************
C SPEAK.F     PARAMETERS                        FEB 2001 ARDEAN LEITH
C             FORMATTING                        MAY 2001 ARDEAN LEITH
C             NDAT/NOUT BUG                     MAY 2007 ARDEAN LEITH
C             REWRITE                           NOV 2012 ARDEAN LEITH
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
C  PURPOSE:  SEARCHES FOR THE ML HIGHEST PEAKS IN THE (REAL) IMAGE
C            FILNAM AND PRINTS OUT POSITIONS AND VALUES OF THESE PEAKS.
C
C  SPEAK(FILNAM,LUN,NX,NY,MAXDIM,OPT,NDOC)
c
C  PARAMETERS:    
C        FILNAM       FILE NAME
C        LUN          LOGICAL UNIT NUMBER OF IMAGE
C        NX,NY        DIMENSIONS OF IMAGE
C        MAXDIM       MAXIMUM BUFFER SPACE AVAILABLE
C        OPT          OUTPUT OPTION
C        OPT=' '      DEFAULT:NO DOCUMENT OUTPUT
C        OPT='X'      I.E.,REGISTER LIST FOLLOWING):NO DOCUMENT OUTPUT
C                       BUT OUTPUT OF POSITION & VALUE OF PEAK IN REGISTERS)
C        OPT='D'      DOCUMENT OUTPUT:NUMBER,POSITION, AND VALUE
C                       OF PEAK ARE WRITTEN INTO A DOCUMENT FILE
C        NDOCOUT      I/O UNIT NUMBER FOR DOCUMENT FILE
C
C        REGISTER POSITIONS 1= INTEGER X-SHIFT
C          2= INTEGER Y-SHIFT
C          3= ABSOLUTE PEAK HEIGHT
C
C          5= FLOATING PT. X-SHIFT
C          6= FLOATING PT. Y-SHIFT
C          7= VALUE OF EXTREMUM OF PARABOLOID
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

         SUBROUTINE SPEAK(FILNAM,LUN,NX,NY,MAXDIM,OPT,NDOCOUT,ML,NOR)

         INCLUDE 'CMBLOCK.INC'
         INCLUDE 'CMLIMIT.INC'

         COMMON  BUF(1)

         CHARACTER (LEN=*)      :: FILNAM
         INTEGER                :: LUN,NX,NY,MAXDIM
         CHARACTER              :: OPT
         INTEGER                :: NDOCOUT,ML,NOR


         REAL                   :: DLIST(5),RSQ(9)
         REAL                   :: ALIST(ML),RLIST(ML)
         INTEGER                :: ILIST(ML),KLIST(ML)
         INTEGER                :: KXLIST(ML),IXLIST(ML)
         CHARACTER (LEN=60)     :: COMMEN
         CHARACTER (LEN=MAXNAM) :: DOCOUT
         LOGICAL                :: ASKNAM,ADDEXT,ISOLD
         LOGICAL                :: WRTCOM,APPEND,NEWFILE

         CALL SET_MPI(ICOMM,MYPID,MPIERR) ! SETS ICOMM AND MYPID

         IF (MAXDIM < 3 * NX) THEN
            CALL ERRT(102,'INSUFFICIENT BUFFER SPACE',MAXDIM)
            RETURN
         ENDIF

         DO K = 1,ML
            ILIST(K) = 0
            KLIST(K) = 0
            RLIST(K) = 0.0
            ALIST(K) = -HUGE(THRESH)
         ENDDO

         NTAB     = 1
         NXCTR    = NX/2+1
         NYCTR    = NY/2+1

         IF  (NOR .NE. 0) THEN
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


30       THRESH = -HUGE(THRESH)
         NMAX   = 0
         NX1    = NX-1
         NY1    = NY -1
         I1     = 0
         I2     = 1
         I3     = 2
         CALL REDLIN(LUN,BUF,NX,1)
         CALL REDLIN(LUN,BUF(1+NX),NX,2)

         DO  I = 3,NY
            I1  = MOD(I1,3)+1
            I2  = MOD(I2,3)+1
            I3  = MOD(I3,3)+1
            I1A = (I1-1)*NX
            I2A = (I2-1)*NX
            I3A = (I3-1)*NX

            CALL REDLIN(LUN,BUF(1+I3A),NX,I)

            DO  K = 2,NX1
               A = BUF(K+I2A)

C              IGNORE PIXEL IF LOWER THAN LOWEST PIXEL ON PEAK LIST 
               IF (A <= THRESH)       CYCLE

C              IGNORE PIXEL IF LESS OR EQUAL TO ANY OF 8 NEIGHBORS
               IF (A <= BUF(K-1+I2A)) CYCLE
               IF (A <= BUF(K-1+I1A)) CYCLE
               IF (A <= BUF(K+I1A))   CYCLE
               IF (A <= BUF(K+1+I1A)) CYCLE
               IF (A <= BUF(K+1+I2A)) CYCLE
               IF (A <= BUF(K+1+I3A)) CYCLE
               IF (A <= BUF(K+I3A))   CYCLE
               IF (A <= BUF(K-1+I3A)) CYCLE

               NMAX = NMAX + 1

               DO L = 1,ML
                  IF (A  <=  ALIST(L)) CYCLE
                     IF (L .NE. ML) THEN
                     L1 = L + 1

                     DO J = L1,ML
                        JO        = ML-J+L1
                        JN        = ML-J+L1-1
                        ALIST(JO) = ALIST(JN)
                        ILIST(JO) = ILIST(JN)
                       KLIST(JO)  = KLIST(JN)
                     ENDDO 
                  ENDIF

                  ALIST(L) = A
                  ILIST(L) = I-1
                  KLIST(L) = K
                  IF (NMAX > ML) THRESH = ALIST(ML)
                  EXIT
               ENDDO
            ENDDO
	 ENDDO

         IF (NMAX == 0)  THEN
            IF (MYPID <= 0) WRITE(NDAT,*) ' NO PEAK FOUND'
            IF (NDAT .NE. NOUT .AND. MYPID <= 0) 
     &         WRITE(NOUT,*) ' NO PEAK FOUND'
            CALL REG_SET_NSEL(1, 5,0.0, 0.0, 0.0, 0.0, 0.0,IRTFLG)
            CALL REG_SET_NSEL(6, 2,0.0, 0.0, 0.0, 0.0, 0.0,IRTFLG)
            RETURN
         ENDIF

        IF (OPT == 'D') THEN
           ADDEXT  = .TRUE.
           ASKNAM  = .TRUE.
           ISOLD   = .FALSE.
           APPEND  = .FALSE.
           WRTCOM  = .TRUE.
           CALL OPENDOC(DOCOUT,ADDEXT,NLET,
     &                  NDOCOUT,NICDOCOUT,ASKNAM,'OUTPUT DOCUMENT',
     &                  ISOLD,APPEND,WRTCOM,NEWFILE,IRTFLG)
           IF (IRTFLG .NE. 0) RETURN
C                123456789 123456789 123456789 123456789 123456789 123
           COMMEN = 
     &      '            XO            YO        HEIGHT         RATIO'
C            123456789 123456789 123456789 123456789 123456789 123456
           CALL LUNDOCPUTCOM(NDOCOUT,COMMEN(1:56),IRTFLG)
         ENDIF

         IF (MYPID <= 0) THEN
            IF (VERBOSE) WRITE(NOUT,*) ' '
            WRITE(NDAT,299)
299         FORMAT(
     &      '    NO    IX    IY     XO    YO      HEIGHT         RATIO')
            IF (NDAT .NE. NOUT) WRITE(NOUT,299)
            IF (VERBOSE) WRITE(NOUT,*) ' '
         ENDIF

         MLIST = MIN(ML,NMAX)

         DO L = 1,MLIST
            KXLIST(L) = KLIST(L) - NXCTR
            IXLIST(L) = ILIST(L) - NYCTR
            RLIST(L)  = ALIST(L) / ALIST(NTAB)

            IF (OPT == 'D') THEN
               DLIST(1) = KXLIST(L)
               DLIST(2) = IXLIST(L)
               DLIST(3) = ALIST(L)
               DLIST(4) = RLIST(L)

C              PUSH DLIST INTO OUTPUT DOC. FILE
               CALL LUNDOCWRTDAT(NICDOCOUT,L,DLIST,4,IRTFLG)

cc               CALL SAVD(NDOC,DLIST,5,IRTFLG)
            ENDIF
            IF (NDAT .NE. NOUT .AND. MYPID <= 0)
     &         WRITE(NOUT,301)L,KLIST(L),ILIST(L),KXLIST(L),IXLIST(L),
     &                        ALIST(L),RLIST(L)
301         FORMAT(3I6,' ',2I6,G16.7,2X,F8.5)
C dgm 301         FORMAT(3I6,' ',2I6G16.7,2X,F8.5)

            IF (MYPID <= 0)  
     &         WRITE(NDAT,301) L,KLIST(L),ILIST(L),KXLIST(L),
     &                         IXLIST(L), ALIST(L),RLIST(L)
         ENDDO

         CLOSE(NDOCOUT)

         CALL REG_SET_NSEL(1,4,FLOAT(KXLIST(1)),FLOAT(IXLIST(1)),
     &                        ALIST(1),RLIST(1),0.0,IRTFLG)

C        9/25/81 PARABOLIC FIT TO 3X3 NEIGHBORHOOD OF THE PEAK POINT
C        PROGRAM SECTION SENT BY M.VAN HEEL, MODIFIED FOR SPIDER. JF

         KL = KLIST(1)

         DO I=1,3
            IROW = ILIST(1)+I-2
            IF (IROW < 1)  IROW = IROW + NY
            IF (IROW > NY) IROW = IROW - NY

            CALL REDLIN(LUN,BUF,NX,IROW)
            I1 = (I-1) * 3

            DO K=1,3
               ISAM = KL + K - 2
               IF (ISAM < 1)  ISAM = ISAM + NX
               IF (ISAM > NX) ISAM = ISAM - NX

               RSQ(I1+K) = BUF(ISAM)
            ENDDO
         ENDDO

         CALL PARABL(RSQ,XSH,YSH,PEAKV)

C        HACK FOR ALMOST BINARY PEAK ERROR
         IF (PEAKV < MAXVAL(RSQ)) PEAKV = MAXVAL(RSQ)

         IF (MYPID <= 0) THEN
            IF (VERBOSE) WRITE(NOUT,*) ' '
            WRITE(NDAT,302)XSH,YSH,PEAKV
            IF (NDAT .NE. NOUT) 
     &         WRITE(NOUT,302)XSH,YSH,PEAKV
302            FORMAT('  SUB-PIXEL OFFSET OF HIGHEST PEAK: (',
     &               F5.2,', ',F5.2,')  HEIGHT: ',G15.7)

            IF (MYPID .LE. 0) WRITE(NOUT,*) ' '
         ENDIF

         XT = XSH + KXLIST(1)
         YT = YSH + IXLIST(1)
 
         CALL REG_SET_NSEL(5,3,XT,YT,PEAKV, 0.0,0.0,IRTFLG)

         END

