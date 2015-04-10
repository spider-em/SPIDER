
C **********************************************************************
C
C  MAPDIST.FOR  -- CREATED JAN 91
C                  OPFILEC                         FEB  03 ARDEAN LEITH
C                  GETDOCDAT                       SEP  09 ARDEAN LEITH
C **********************************************************************
C * AUTHOR: ArDean Leith                                               *
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
C    MAPDIST(MAPIMAG,DOCRES,IRTFLG) 
C
C    PURPOSE:       READS A DOC. FILE CONTAINING CLUSTER NUMBER (KEY),
C                   NUMBER OF VOXELS (REG 1) IN THE CLUSTER, AND COM. 
C                   DETERMINES DISTANCE BETWEEN ALL COMS.  IF DISTANCE
C                   BETWEEN PRESENT CLUSTER AND TEST CLUSTER IS LESS
C                   THAN THRESHOLD, SETS TEST CLUSTER TO PRESENT CLUSTER.
C                   REMAPS IMAGE TO NEW CLUSTER NUMBERS.  
C
C      PARAMETERS:  MAPIMAG          MAP INTO NEW FILE  
C                   DOCRES           OUTPUT NEW DOC FILE
C                   IRTFLG           ERROR FLAG
C
C      CALLED BY:   UTIL6
C
C      CALLS:       MAKTAB   
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE MAPDIST(MAPIMAG,DOCRES,IRTFLG)                                                                                                

       INCLUDE 'CMBLOCK.INC'
       INCLUDE 'CMLIMIT.INC' 
       INCLUDE 'F90ALLOC.INC'

       PARAMETER  (MAXLUT=10000)       ! ARBITRARY SIZE
       PARAMETER  (NEQMAX=16000)       ! ARBITRARY SIZE

       REAL                           :: TABLE(MAXLUT)         
       INTEGER                        :: IEQUIV(2,NEQMAX)     

       CHARACTER(LEN=MAXNAM)          :: DOCNAM,IMFILE,OUTFILE
       LOGICAL                        :: MAPIMAG,DEBUGGING,DOCRES
 
       REAL                           :: BUF(NBUFSIZ)  ! FROM CMLIMIT
       REAL, POINTER                  :: DBUF(:,:)
       REAL                           :: DLIST(5)

       INTEGER, PARAMETER             :: LUND   = 81   
       INTEGER, PARAMETER             :: LUNIM  = 21    
       INTEGER, PARAMETER             :: LUNOUT = 22     
       INTEGER, PARAMETER             :: LUND2  = 82     

       DEBUGGING = .FALSE.

       WRITE(NOUT,90)
90     FORMAT(/,
     & '  Reads a document file from the <EC STAT> output',/,
     & '  and merges all clusters having center of mass less than   ',/,
     & '  threshold distance away into a new overall cluster.'/)

C      RETRIEVE DBUF ARRAY 
       MAXXT = 5
       MAXYT = 0
       CALL GETDOCDAT('CENTER-OF-MASS DOC',.TRUE.,DOCNAM,LUND,
     &                 .TRUE.,MAXXT, NVAL,DBUF,IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 999

       WRITE(NOUT,*) ' READ DATA FOR :',NVAL,' CLUSTERS.'
       IF (NVAL > MAXLUT) THEN
C         TABLE OVERFLOW
          MAXL = MAXLUT
          WRITE(NOUT,*)' *** IN MAPDIST, TABLE LIMIT:',MAXL,' GOT:',NVAL
          CALL ERRT(100,'MAPDIST',NE)
          GOTO 999
       ENDIF

       CALL RDPRM2S(SEPMIN,FMINCLU,NOT_USED,
     &    'THRESHOLD DISTANCE & INITIAL CLUSTER',IRTFLG)
       IF (IRTFLG .NE. 0) GOTO 999
        
       MINCLU = FMINCLU

C      MINIMUM SEPARATION FOR ANY CELL PAIR (IN UM) SQ
       SEPMINSQ = SEPMIN * SEPMIN

C      INITIALIZE NUMBER OF EQUIVAENT CLUSTERS
       NEQUIV   = 0

       DO I = MINCLU,NVAL-1

C        FIND IDENTITY OF PRESENT CLUSTER
         ICLU1 = DBUF(1,I) 
         N1    = DBUF(2,I) 
         X1    = DBUF(3,I)
         Y1    = DBUF(4,I) 
         Z1    = DBUF(5,I) 

         IF (N1 <= 0) THEN
C           IMPOSSIBLE NUMBER OF VALUES IN THIS CLUSTER
            WRITE(NOUT,'(A,I8,A,I9,A)')
     &                            '  *** SKIPPING CLUSTER:',I,
     &                            ' WITH:',N1, ' OCCUPANTS'
            CYCLE
         ENDIF

         DO J = I+1,NVAL

C           FIND LOCATION OF TEST CLUSTER

            ICLU2 = DBUF(1,J)
            N2    = DBUF(2,J)
            X2    = DBUF(3,J) 
            Y2    = DBUF(4,J) 
            Z2    = DBUF(5,J) 

            IF (N2 <= 0) THEN
C               IMPOSSIBLE NUMBER OF VALUES IN CLUSTER
               WRITE(NOUT,'(A,I8,A,I9,A)')
     &                            '  *** Skipping cluster:',J,
     &                            ' WITH:',N2, ' occupants'
                CYCLE
            ENDIF

C           FIND SEPARATION BETWEEN TWO ORIGINAL CLUSTERS
            DXSQ = (X1 - X2)**2      
            DYSQ = (Y1 - Y2)**2
            DZSQ = (Z1 - Z2)**2

            DSQ  = DXSQ + DYSQ + DZSQ
    
            IF (VERBOSE) THEN

               DIST = SQRT(DSQ)
               WRITE(NOUT, '(A,I7,A,I7, A,F8.2)') 
     &            '  Clusters:' ,I,' -->',J,
     &            '   Distance:', DIST

            ENDIF

            IF (DSQ <= SEPMINSQ) THEN
C              PUT IN SAME OVERALL CLUSTER
               NEQUIV = NEQUIV + 1
               IEQUIV(1,NEQUIV) = I
               IEQUIV(2,NEQUIV) = J
            ENDIF

          ENDDO
       ENDDO

       WRITE(NOUT,94) NEQUIV
94     FORMAT(/,'  Number of equivalences: ',I8)

       CALL MAKTAB(IEQUIV,NEQUIV,TABLE,NVAL,NTAB, IRTFLG)

       IF (MAPIMAG) THEN
C         USE TABLE TO MAP OLD IMAGE TO NEW IMAGE VALUES

          MAXIM = 0
          CALL OPFILEC(0,.TRUE.,IMFILE,LUNIM,'O',IFORM,NX,NY,NZ,
     &       MAXIM,'CLUSTER INPUT',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 999

          MAXIM = 0
          CALL OPFILEC(LUNIM,.TRUE.,OUTFILE,LUNOUT,'U',IFORM,
     &        NX,NY,NZ,MAXIM,'CLUSTER OUTPUT',.FALSE.,IRTFLG)
          IF (IRTFLG .NE. 0) GOTO 999

          NREC2 = NY * NZ
          CALL MAPIM(LUNIM,LUNOUT,NX,1,NREC2,TABLE,NDUM,BUF,IRTFLG)

       ELSEIF (DOCRES) THEN
C         SAVE RESULTS IN A DOC FILE
          DO I = 1,NVAL
            IF (TABLE(I) > 0) THEN
               DLIST(1) = I
               DLIST(2) = TABLE(I)
               DLIST(3) = DBUF(3,I)
               DLIST(4) = DBUF(4,I)
               DLIST(5) = DBUF(5,I)
               CALL SAVD(LUND2,DLIST,5,IRTFLG)

               IF (DEBUGGING) THEN
C                 WRITE STATISTICS ON STANDARD OUTPUT
                  WRITE(NOUT,901) I-1,INUM,(DLIST(J),J=3,5)
901               FORMAT(2I10,3G12.3)
               ENDIF
            ENDIF
         ENDDO
         CLOSE(LUND2)                 
       ENDIF
      
       IRTFLG = 0 

 999   CLOSE(LUNIM)
       CLOSE(LUNOUT)
       CLOSE(LUND)
       IF (ASSOCIATED(DBUF)) DEALLOCATE(DBUF)
         
       END

