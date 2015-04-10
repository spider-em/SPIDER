C++*********************************************************************
C
C  WARP.F           PORTED FROM WARP.C           DEC 03   ARDEAN LEITH
C                   CONVERTED TO 3D              MAR 04   ARDEAN LEITH
C                   CLOSE LUNOUT                 SEP 08   ARDEAN LEITH
C                   IX,IU OVERFLOW TRAP          OCT 09   ARDEAN LEITH
C                   X WARP BUG                   MAY 11   ARDEAN LEITH
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011  Health Research Inc.,                         *
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
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C 
C  Three pass mesh warping based on algorithm in Wolberg  
C  input image in has height NROW and width NSAM          
C  xs, ys contain coor. of source mesh                  
C  xd, yd contain coor. of destination mesh             
C  their height and width dimensions are th and tw.    
C  The output is stored in out. Due to frozen edge     
C  assumptions, out has same dimensions as in.           
C   
C   
C   xs,xd  (tw wide x th high)        First pass (phase one): Create tables
C     |                               ts and ti for x-intercepts of
C    col to row                       vertical splines in s and i.
C     |
C     v
C   xrow1 (th)         indx(NROW)
C     |                  |
C     v                  v
C   map1, map2 ( NROW  )
C     |
C    row to col
C     |
C     v
C   ts,ti ( tw x NROW)                 First pass (phase two): warp x using 
C     |                               ts and ti.
C     |
C    row
C     |
C     v
C   x1,x2 (tw)         indx(NSAM)
C     |                  |
C     v                  v
C   map1, map2 ( NSAM  )
C     |
C     |               src (NSAM)
C     |                  |
C     rows              row
C     |                  |
C     v                  v
C   temp (NSAM x NROW)
C   
C    ----------------------------------------------------------
C   
C   xs,xd  (tw wide x th high)        Second pass (phase one): Create tables
C     |                               ts and ti for y-intercepts of
C    row                              horiz splines in ti and td.
C     |
C     v         indx(NROW)
C   x1,y1(tw)       |
C     |             |
C     v             v
C   ti,td ( NSAM x th )                Second pass (phase two): warp y using 
C     |                               ti and td.
C    col
C     |
C     v
C   xrow1,yrow1 (th)     indx(NSAM)
C     |                  |
C     v                  v
C   map1( NSAM  )
C     |
C     |               temp (NSAM x NROW)
C     |                  |
C     row               col
C     |                  |
C     |                  v
C     |               mine1 (NROW)
C     |                  |
C     v                  v
C   mine2 (NROW)
C     |
C    col to row
C     |
C     v
C   output (NSAM x NROW)
C   
C   
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

      SUBROUTINE WARP3(LUNDOCT,LUNIN,LUNOUT)

      INCLUDE 'CMBLOCK.INC'
      INCLUDE 'CMLIMIT.INC'
      INCLUDE 'F90ALLOC.INC'

      CHARACTER(LEN=MAXNAM)             :: DOCNAM,DOCNAMPE,FILNAM
      CHARACTER(LEN=1)                  :: NULL
      LOGICAL                           :: NEWFILE
      REAL, DIMENSION(7)                :: DLIST

      REAL, DIMENSION(:,:), POINTER     :: DBUF
      REAL, ALLOCATABLE, DIMENSION(:)   :: BUFIN,BUFOUT

      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XS,YS,ZS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XD,YD,ZD
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: X ,Y ,Z 
      REAL, ALLOCATABLE, DIMENSION(:)     :: XV,YV,ZV

      NULL   = CHAR(0)
      IRTFLG = 1

      CALL OPENDOC(DOCNAM,.TRUE.,NLET,LUNDOCT,LUNDOC,.TRUE.,
     &    'WARP COORDINATES DOC.',.TRUE.,.FALSE.,.FALSE.,NEWFILE,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     WANT A COMMENTED KEY
      ICOUNT = 3
      CALL LUNDOCGETCOM(LUNDOC,-1,DLIST,ICOUNT,.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
      CLOSE(LUNDOC)

      IW = DLIST(1) + 2
      IH = DLIST(2) + 2
      ID = DLIST(3) + 2
      
      WRITE(NOUT,90) iw,ih,id
90    FORMAT('  Warping mesh dimensions: ',I3,' x ',I3,' x ',I3,/)

C     MERGE DOCNAM WITH DATEXC
      CALL FILNAMANDEXT(DOCNAM,DATEXC,DOCNAMPE,NLETPE,.TRUE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     READ THE SOURCE  MESH    
C     MAXX IS 1 + NUM OF REGISTERS SINCE DBUF CONTAINS KEY ALSO
      MAXX = 7
      MAXY = IW * IH * ID
      CALL GETDOCDAT(' ',.FALSE., DOCNAMPE,LUNDOC,.FALSE.,MAXX, MAXY,
     &               DBUF,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     OPEN INPUT FILE, NO FOURIER INPUT ALLOWED 
      MAXIM = 0
      CALL OPFILEC(0,.TRUE.,FILNAM,LUNIN,'O',ITYPE,NSAM,NROW,NSLICE,
     &             MAXIM,'INPUT',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     ALLOCATE MEMORY FOR INPUT AND OUTPUT IMAGES
      ALLOCATE(BUFIN(NSAM*NROW*NSLICE),
     &         BUFOUT(NSAM*NROW*NSLICE),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         MWANT = NSAM*NROW*NSLICE
         CALL ERRT(46,'BUFIN & BUFOUT',MWANT)
         RETURN
      ENDIF

C     FILL THE BUFIN ARRAY WITH FIRST FILE
      CALL REDVOL(LUNIN,NSAM,NROW,1,NSLICE,BUFIN,IRTFLG)
      CLOSE(LUNIN)

C     ALLOCATE MEMORY FOR INPUT AND OUTPUT MESHES
      ALLOCATE(XS(IW,IH,ID), YS(IW,IH,ID), ZS(IW,IH,ID),
     &         XD(IW,IH,ID), YD(IW,IH,ID), ZD(IW,IH,ID),
     &         X (IW,IH,ID), Y (IW,IH,ID), Z (IW,IH,ID),
     &         XV (IW),      YV(IH),       ZV(ID),
     &         STAT=IRTFLG)
      IF (IRTFLG .NE. 0) THEN
         MWANT = 9 *IW*IH*ID + IW + IH + ID
         CALL ERRT(46,'XS....',MWANT)
         RETURN
      ENDIF

C     RETRIEVE XS,YS.. & XD,YD..  COOR. FROM DBUF  

C     CORE VALUES
      IKEY = 1 
      DO K=2,ID-1
         DO J=2,IH-1
            DO I=2,IW-1
               XS(I, J, K) = DBUF(2 , IKEY)
               YS(I, J, K) = DBUF(3 , IKEY)
               ZS(I, J, K) = DBUF(4 , IKEY)

               XD(I, J, K) = DBUF(2 , IKEY) + DBUF(5 , IKEY) 
               YD(I, J, K) = DBUF(3 , IKEY) + DBUF(6 , IKEY)
               ZD(I, J, K) = DBUF(4 , IKEY) + DBUF(7 , IKEY)

               IKEY = IKEY + 1
            ENDDO
         ENDDO
      ENDDO

C     WALL VALUES
      XS(:,  1 , : ) = XS(:   , 2    , :)
      XS(:,  IH, : ) = XS(:   , IH-1 , :)

      XS(1,  : , : ) = XS(2   , :    , :)
      XS(IW, : , : ) = XS(IW-1, :    , :)

      XS(:,  : , 1 ) = XS(:   , :    , 2)
      XS(:,  : , ID) = XS(:   , :    , ID-1)


      YS(:,  1 , : ) = YS(:   , 2    , :)
      YS(:,  IH, : ) = YS(:   , IH-1 , :)

      YS(1,  : , : ) = YS(2   , :    , :)
      YS(IW, : , : ) = YS(IW-1, :    , :)

      YS(:,  : , 1 ) = YS(:   , :    , 2)
      YS(:,  : , ID) = YS(:   , :    , ID-1)

      ZS(1,  : , : ) = ZS(2   , :    , :)
      ZS(IW, : , : ) = ZS(IW-1, :    , :)

      ZS(:,  1 , : ) = ZS(:   , 2    , :)
      ZS(:,  IH, : ) = ZS(:   , IH-1 , :)

      ZS(:,  : , 1 ) = ZS(:   , :    , 2)
      ZS(:,  : , ID) = ZS(:   , :    , ID-1)

C     WALL VALUES
      XD(:,  1 , : ) = XD(:   , 2    , :)
      XD(:,  IH, : ) = XD(:   , IH-1 , :)

      XD(1,  : , : ) = XD(2   , :    , :)
      XD(IW, : , : ) = XD(IW-1, :    , :)

      XD(:,  : , 1 ) = XD(:   , :    , 2)
      XD(:,  : , ID) = XD(:   , :    , ID-1)


      YD(:,  1 , : ) = YD(:   , 2    , :)
      YD(:,  IH, : ) = YD(:   , IH-1 , :)

      YD(1,  : , : ) = YD(2   , :    , :)
      YD(IW, : , : ) = YD(IW-1, :    , :)

      YD(:,  : , 1 ) = YD(:   , :    , 2)
      YD(:,  : , ID) = YD(:   , :    , ID-1)

      ZD(1,  : , : ) = ZD(2   , :    , :)
      ZD(IW, : , : ) = ZD(IW-1, :    , :)

      ZD(:,  1 , : ) = ZD(:   , 2    , :)
      ZD(:,  IH, : ) = ZD(:   , IH-1 , :)

      ZD(:,  : , 1 ) = ZD(:   , :    , 2)
      ZD(:,  : , ID) = ZD(:   , :    , ID-1)

C     SET UP BOUNDARY VALUES
      XS(1,  :, :)  = 1.0
      XS(IW, :, :)  = NSAM

      YS(:, 1,  :) = 1.0
      YS(:, IH, :) = NROW

      ZS(:, :,  1) = 1.0
      ZS(:, :, ID) = NSLICE

      IF (NSLICE .EQ. 1) THEN
C        WARPING IMAGE NOT VOLUME al
         ZS(:, :,  1) = 0.0
         ZS(:, :, ID) = 2
      ENDIF

C     SET UP BOUNDARY VALUES
      XD(1,  :, :)  = 1.0
      XD(IW, :, :)  = NSAM

      YD(:, 1,  :) = 1.0
      YD(:, IH, :) = NROW

      ZD(:, :,  1) = 1.0
      ZD(:, :, ID) = NSLICE

      X = XD - XS
      Y = YD - YS
      Z = ZD - ZS


#ifdef DEBUG
      write(6,*) "---- x input mesh ---- "
      WRITE(6,91) XS


      write(6,*) "---- y input mesh ---- "
      WRITE(6,91) YS
      write(6,*) "---- z input mesh ---- "
      WRITE(6,91) ZS

      write(6,*) "---- x output mesh ---- "
      WRITE(6,91) XD
      write(6,*) "---- y output mesh ---- "
      WRITE(6,91) YD
      write(6,*) "---- z output mesh ---- "
      WRITE(6,91) ZD

      write(6,*) "---- x offset mesh ---- "
      WRITE(6,91) X
      write(6,*) "---- y offset mesh ---- "
      WRITE(6,91) Y
      write(6,*) "---- z offset mesh ---- "
      WRITE(6,91) Z
#endif

91    FORMAT(7(7(F7.2,' ')/)/)

      XV = XS(:,1,1)  
      YV = YS(1,:,1)
      ZV = ZS(1,1,:)

C     DO THE WARPING
      CALL WARP_MESH(XV,YV,ZV, X,Y,Z, IW,IH,ID, NSAM,NROW,NSLICE,
     &               BUFIN,BUFOUT,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     OUTPUT THE WARPED VOLUME
      !write(6,*) "OUTPUTTING WARPED VOL ",nsam,nrow,nslice,itype

      MAXIM = 0
      ITYPE = 3
      IF (NSLICE .LE. 1) ITYPE = 1
      CALL OPFILEC(0,.TRUE.,FILNAM,LUNOUT,'U',ITYPE,NSAM,NROW,NSLICE,
     &             MAXIM,'OUTPUT',.FALSE.,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

      IF (NSLICE .GT. 1) THEN
         CALL WRTVOL(LUNOUT,NSAM,NROW,1,NSLICE,BUFOUT,IRTFLG)
      ELSE
C        BUFFERS ARE FLIPPED
         CALL WRTVOL(LUNOUT,NSAM,NROW,1,NSLICE,BUFIN,IRTFLG)
      ENDIF

      CLOSE(LUNIN)
      CLOSE(LUNOUT)

      END


C     ******************* WARP_MESH ********************************* 

      SUBROUTINE WARP_MESH(XV,YV,ZV, X,Y,Z, IW,IH,ID, NSAM,NROW,NSLICE, 
     &                     BUFIN,BUFOUT, IRTFLG)

      INCLUDE 'CMBLOCK.INC'

      REAL, DIMENSION(*)                :: XV
      REAL, DIMENSION(*)                :: YV
      REAL, DIMENSION(*)                :: ZV
      REAL, DIMENSION(IW*IH*ID)         :: X, Y, Z
      REAL, DIMENSION(NSAM,NROW,NSLICE) :: BUFIN,BUFOUT

      !REAL, DIMENSION(ID,IH,IW)         :: XT
       REAL, DIMENSION(IH,ID,IW)         :: XT
       REAL, DIMENSION(IW,ID,IH)         :: YT 

C     3-PASS MESH WARPING BASED ON ALGORITM IN SMYTHE (90).    
C     INPUT IMAGE IS: NSAM x NROW x NSLICE.               
C     MESH SIZE IS IW x IH x ID
C     OUTPUT IS STORED IN BUFOUT. BUFOUT MUST HAVE SAME DIMENSION 
C     AS BUFIN 

      WRITE(NOUT,*) ' WARPING IN: X'
      ILOC = 0

      DO IY = 1,ID

         DO IX = 1,IH

            DO IZ = 1,IW
               ILOC         = ILOC + 1
               XT(IX,IY,IZ) = X(ILOC)
            ENDDO
         ENDDO
      ENDDO

#ifdef DEBUG
      write(6,*) "---- x warping mesh ---- ",id,ih,iw
      WRITE(6,91) XT
91    FORMAT(7(7(F7.2,' ')/))
#endif

      CALL WARP_PASS(1, YV,ZV,XT, IH,ID,IW, NROW,NSLICE,NSAM, 
     &               NSAM,NROW,NSLICE, BUFIN,BUFOUT, IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

c-----------------------------
c     OUTPUT THE WARPED VOLUME
      MAXIM = 0
      ITYPE = 3
      IF (NSLICE .LE. 1) ITYPE = 1
      CALL OPFILEC(0,.false.,'jnk-warpedx',51,'U',ITYPE,
     &             NSAM,NROW,NSLICE,
     &             MAXIM,' ',.FALSE.,IRTFLG)
      CALL WRTVOL(51,NSAM,NROW,1,NSLICE,BUFOUT,IRTFLG)
      close(51)
c-----------------------------

      WRITE(NOUT,*) ' WARPING IN: Y'
      ILOC = 0
      DO IY = 1,ID
         DO IZ = 1,IH
            DO IX = 1,IW
               ILOC         = ILOC + 1
               YT(IX,IY,IZ) = Y(ILOC)
            ENDDO
         ENDDO
      ENDDO

c-----------------------------
      CALL WARP_PASS(2, XV,ZV,YT, IW,ID,IH, NSAM,NSLICE,NROW, 
     &               NSAM,NROW,NSLICE, BUFOUT,BUFIN, IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

c     OUTPUT THE WARPED VOLUME
      MAXIM = 0
      ITYPE = 3
      IF (NSLICE .LE. 1) ITYPE = 1
      CALL OPFILEC(0,.false.,'jnk-warpedy',51,'U',ITYPE,
     &             NSAM,NROW,NSLICE,
     &             MAXIM,' ',.FALSE.,IRTFLG)
      CALL WRTVOL(51,NSAM,NROW,1,NSLICE,BUFIN,IRTFLG)
      close(51)

c-----------------------------
      IF (NSLICE .GT. 1) THEN
         WRITE(NOUT,*) ' WARPING IN: Z'

         CALL WARP_PASS(3, XV,YV,Z, IW,IH,ID, NSAM,NROW,NSLICE, 
     &                  NSAM,NROW,NSLICE, BUFIN,BUFOUT, IRTFLG)
         IF (IRTFLG .NE. 0) RETURN
      ENDIF

      END


C     ******************** WARP_PASS **********************************

      SUBROUTINE WARP_PASS(IPASS, XV,YV,Z, IW,IH,ID, NSAM,NROW,NSLICE, 
     &                     NSAMO,NROWO,NSLICEO, BUFIN,BUFOUT, IRTFLG)


      REAL, DIMENSION(*)                   :: XV
      REAL, DIMENSION(*)                   :: YV
      REAL, DIMENSION(IW,IH,ID)            :: Z
      REAL, DIMENSION(NSAMO,NROWO,NSLICEO) :: BUFIN,BUFOUT

C     AUTOMATIC ARRAYS
      REAL, DIMENSION(NSAM*NROW)           :: XI,YI 
      REAL, DIMENSION(NSAM*NROW,ID)        :: ZI
      REAL, DIMENSION(NSLICE)              :: Z2,ZVALS,ZIVALS,ZOVALS

      REAL, DIMENSION(ID)                  :: YT

C     WK    = THREE DIMENSIONAL ARRAY USED AS A WORK AREA.
      REAL, DIMENSION(3*IW*IH)             :: WK

C      (PHASE ONE): CREATE TABLES  AND  FOR X-INTERCEPTS OF 
C     VERTICAL SPLINES IN S AND I. TABLES HAVE IW COLUMNS OF HEIGHT NSLICE   
C     INDICES USED TO SAMPLE VERTICAL SPLINES  
      NDP = IW * IH

      NIP = NSAM * NROW
      IX  = NSAM
      IY  = 0
      DO I = 1,NIP
         IX = IX + 1 
         IF (IX .GT. NSAM) THEN
            IX = 1 
            IY = IY + 1 
         ENDIF
         XI(I) = IX
         YI(I) = IY
      ENDDO

      MD = 1

      DO IZ = 1,ID

C        RGBI3P  PERFORMS INTERPOLATION  OF A BIVARIATE FUNCTION,
C        Z(X,Y), ON A RECTANGULAR GRID IN THE X-Y PLANE. IT IS 
C        BASED ON A REVISED AKIMA METHOD.

         CALL RGBI3P(MD,IW,IH,XV, YV, Z(1,1,IZ),
     &               NIP, XI,YI, ZI(1,IZ), IRTFLG, WK)

         IF (IRTFLG .NE. 0) THEN
            CALL ERRT(102,'WARPING IN RGBI3P',IRTFLG)
            RETURN
         ENDIF
      ENDDO  ! END OF DO IZ = 1,ID

C     (PHASE TWO): WARP Z USING: XI,YI,ZI -------------

C     SET INDICES USED TO SAMPLE SURFACES
      DO IZ = 1,NSLICE
         Z2(IZ) = IZ
      ENDDO

      ILOC = 0
      DO IY=1,NROW        !     VISIT EACH ROW  OF VOLUME  

C (I+1)not allowed    parallel do private(ix,iloc,yt,z2,zvals,zivals,zovals)
         DO IX=1,NSAM     !     VISIT EACH COL. OF VOLUME  

            ILOC = ILOC + 1

C           LOAD YT WITH  DELTA Z FOR EACH XY PLANE OF CONTROL PLANES
            YT   = ZI(ILOC,:)

C           FIT SPLINE TO XY-INTERCEPTS; RESAMPLE OVER ALL SLICES. 
C           RETURNS ZVALS: DELTA Z FOR EACH SLICE

            CALL ISPLINE_GEN(XV, YT, ID, 
     &                       Z2, ZVALS, NSLICE, IRTFLG) 
            IF (IRTFLG .NE. 0) RETURN

C           CHANGE ZVALS TO ACTUAL Z FOR EACH SLICE
            ZVALS = Z2 - ZVALS - 1.0

            IF (IPASS .EQ. 1) THEN

C              LOAD ZIVALS WITH UNWARPED SOURCE VALUES FOR THIS COL.
               ZIVALS = BUFIN(:,IX,IY)

C              RESAMPLE SOURCE COL. (ZIVALS) BASED ON ZVALS 
C              RETURNS ZOVALS: DELTA Z FOR EACH SLICE
               CALL RESAMPLE_GEN(ZVALS, ZIVALS, ZOVALS, NSLICE)

C              STORE THE WARPED VALUES IN BUFOUT
               BUFOUT(:,IX,IY) = ZOVALS

            ELSE IF (IPASS .EQ. 2) THEN

C              LOAD ZIVALS WITH UNWARPED SOURCE VALUES FOR THIS COL.
               ZIVALS = BUFIN(IX,:,IY)

C              RESAMPLE SOURCE ROW BASED ON ZVALS 
               CALL RESAMPLE_GEN(ZVALS, ZIVALS, ZOVALS, NSLICE)

C              STORE THE WARPED VALUES IN BUFOUT
               BUFOUT(IX,:,IY) = ZOVALS

            ELSEIF (IPASS .EQ. 3) THEN

C              LOAD ZIVALS WITH UNWARPED SOURCE VALUES FOR THIS COL.
               ZIVALS = BUFIN(IX,IY,:)

C              RESAMPLE SOURCE ROW BASED ON ZVALS 
               CALL RESAMPLE_GEN(ZVALS, ZIVALS, ZOVALS, NSLICE)

C              STORE THE WARPED VALUES IN BUFOUT
               BUFOUT(IX,IY,:) = ZOVALS
            ENDIF

         ENDDO
      ENDDO

      IRTFLG = 0
      END
      

C     ******************* ISPLINE_GEN ********************************

      SUBROUTINE ISPLINE_GEN(X1,Y1,LEN1, 
     &                       X2,Y2,LEN2, IRTFLG)

C     INTERPOLATES CUBIC SPLINE FUNCTION FOR IRREGULARLY-SPACED POINTS 
C     INPUT: X1,Y1 IS LIST OF IRREGULAR DATA POINTS (LEN1 ENTRIES) 
    
C     OUTPUT Y2 IS CUBIC SPLINE SAMPLED ACCORDING TO ID (LEN2 ENTRIES) 
C     ASSUME THAT X1,Y1 ENTRIES ARE MONOTONICALLY INCREASING           

      REAL, DIMENSION(0:LEN1-1)             :: X1,Y1
      REAL, DIMENSION(0:LEN2-1)             :: X2,Y2

C     AUTOMATIC ARRAY
      DOUBLE PRECISION, DIMENSION(0:LEN1-1) :: YD

      IRTFLG = 1

C     COMPUTE 1ST DERIV. AT EACH POINT -> YD 
      CALL GETYD_GEN(X1,Y1,YD,LEN1,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     ERROR CHECKING 
      IF (X2(0) .LT. X1(0) .OR. X2(LEN2-1) .GT. X1(LEN1-1)) THEN
          CALL ERRT(101,'OUT OF RANGE IN ISPLINE_GEN3',NE)
          RETURN
      ENDIF

C     P1 IS LEFT  ENDPT. OF INTERVAL,  P2 IS RESAMPLING POSITION 
C     P3 IS RIGHT ENDPT. OF INTERVAL                            
C     J IS INPUT INDEX FOR CURRENT INTERVAL                     

      P3 = X2(0) - 1              ! FORCE COEFF. INIT.
 
      J = 0
      DO I=0,LEN2-1
C        CHECK IF IN NEW INTERVAL 
         P2 = X2(I)
         IF (P2 .GT. P3) THEN
C           FIND THE INTERVAL WHICH CONTAINS P2 
            DO WHILE (J .LT. LEN1 .AND. P2 .GT. X1(J) )
              J = J + 1 
            ENDDO

            IF (P2 .LT. X1(J)) J = J - 1
            P1 = X1(J)         ! UPDATE LEFT ENDPOINT 
            P3 = X1(J+1)       ! UPDATE RIGHT ENDPOINT 

C           COMPUTE SPLINE COEFF.
            DX = 1.0 / (X1(J+1) - X1(J))
            DY = (Y1(J+1) - Y1(J)) * DX
            A0 = Y1(J)
            A1 = YD(J)
            A2 = DX * (3.0 * DY - 2.0 * YD(J) - YD(J+1))
            A3 = DX * DX * (-2.0 * DY + YD(J) + YD(J+1))
         ENDIF

C        USE HORNER'S RULE TO CALCULATE CUBIC POLYNOMIAL
         X     = P2 - P1
         Y2(I) = ((A3 * X + A2) * X + A1) * X + A0
      ENDDO

      IRTFLG = 0

      END


C      *************** RESAMPLE_GEN ***********************************

       SUBROUTINE RESAMPLE_GEN(F, BUFIN, BUFOUT, LENOUT)

       DOUBLE PRECISION                       :: ACC,INTENSITY,INSFAC
       DOUBLE PRECISION                       :: INSEG,OUTSEG
       DOUBLE PRECISION,DIMENSION(0:LENOUT-1) :: INPOS,INPOST
       REAL,DIMENSION(0:LENOUT-1)             :: F,BUFIN,BUFOUT

C      PRECOMPUTE INPUT INDEX FOR EACH OUTPUT PIXEL 
      
       IU = 0
       DO IX=0,LENOUT-1

CC        DO WHILE (F(IU+1) .LT. IX) 
          DO WHILE (IU .LE. (LENOUT -3) .AND. F(IU+1) .LT. IX) 
             IU = IU + 1
          ENDDO
          if (iu .ge. (lenout-1)) THEN
             iuold = iu
             iu    = lenout-2
             write(6,*)'iu fixed:',iuold,'...',iu
          endif

          INPOS(IX) = IU +  (FLOAT(IX) - F(IU)) / (F(IU+1)- F(IU))

          if (inpos(ix) .gt. lenout) THEN
             told      = inpos(ix)
             inpos(ix) = lenout-1
c             write(6,*)'inpos(',ix,') fixed:',told,'...',inpos(ix)

          elseif (inpos(ix) .lt. 0) THEN
             told      = inpos(ix)
             inpos(ix) = 0
c            write(6,*)'inpos(',ix,') fixed:',told,'...',inpos(ix)
         endif
            
          if (inpos(IX) .gt. lenout .or. inpos(ix) .lt. 0) then
              write(6,*) 'bad inpos() for ix:',ix,inpos(ix)
              write(6,*) 'iu,fiu,fiu+1:',iu,f(iU),f(iu+1)
              stop
          endif
       ENDDO

       INSEG  = 1.0        ! ENTIRE INPUT PIXEL IS AVAILABLE             
       OUTSEG = INPOS(1)   ! # INPUT PIXELS THAT MAP ONTO 1 OUTPUT PIXEL 
       INSFAC = OUTSEG     ! INVERSE SCALE FACTOR                        
       ACC    = 0.0        ! CLEAR ACCUMULATOR                           

C      COMPUTE ALL OUTPUT PIXELS 
       IU = 0
       IX = 0
       DO WHILE(IX .LT. LENOUT)
C         USE LINEAR INTERPOLATION FOR RECONSTRUCTION
          IF (IU .LT. (LENOUT -1)) THEN           ! al dec 04???
             INTENSITY = (INSEG * BUFIN(IU)) + 
     &                   ((1 - INSEG) *BUFIN(IU + 1))
          ELSE
             INTENSITY = BUFIN(LENOUT-1)          ! al dec 04???
          ENDIF

          IF (INSEG .LT. OUTSEG) THEN
C            INPUT PIXEL IS ENTIRELY CONSUMED BEFORE OUTPUT PIXEL
             ACC    = ACC  + (INTENSITY * INSEG) ! ACCUMULATE WEIGHTED CONTRIBUTION  
             OUTSEG = OUTSEG - INSEG             ! INSEG PORTION HAS BEEN FILLED      
             INSEG  = 1.0                        ! NEW INPUT PIXEL WILL BE AVAILABLE  
             IU     = IU + 1                     ! INDEX INTO NEXT INPUT PIXEL 
          ELSE
C            INPUT PIXEL IS NOT ENTIRELY CONSUMED FIRST 
             ACC = ACC  + (INTENSITY * OUTSEG)    ! ACCUM. WEIGHTED CONTRIB.                
             if (insfac .ne. 0) then
                BUFOUT(IX) = ACC / INSFAC         ! INIT OUTPUT WITH NORMALIZED ACCUMULATOR 
             else
                BUFOUT(IX) = BUFIN(IX)            ! al apr 04???
             endif

             ACC        = 0.0                     ! RESET ACCUMULATOR FOR NEXT OUTPUT PIX   
             INSEG      = INSEG - OUTSEG          ! OUTSEG PORTION OF INPUT HAS BEEN USED   
             IX         = IX + 1                  ! INDEX INTO NEXT OUTPUT PIXEL             
             if (ix .lt. (lenout -1)) then        ! al oct09
                INSFAC  = INPOS(IX+1) - INPOS(IX) ! INIT SPATIALLY VARYING INSFAC 
             else
                insfac  = 0                       ! INIT SPATIALLY VARYING INSFAC al oct09
             endif
             OUTSEG     = INSFAC                  ! INIT SPATIALLY VARYING SIZFAC 
          ENDIF
      ENDDO

      END


C     ***************** GETYD_GEN *************************************

      SUBROUTINE GETYD_GEN(X, Y, YD, LEN, IRTFLG)

      REAL, DIMENSION(0:LEN-1)             :: X,Y
      DOUBLE PRECISION, DIMENSION(0:LEN-1) :: YD
      DOUBLE PRECISION                     :: H0, H1, R0, R1

C     AUTOMATIC ARRAYS FOR TRIDIAGONAL BANDS A,B, & C 
      DOUBLE PRECISION, DIMENSION(0:LEN-1) :: A,B,C

C     INIT FIRST ROW DATA 
      H0    = X(1) - X(0)
      H1    = X(2) - X(1)
      R0    = (Y(1) - Y(0)) / H0
      R1    = (Y(2) - Y(1)) / H1
      B(0)  = H1 * (H0 + H1)
      C(0)  = (H0 + H1) **2 
      YD(0) = R0 * (3 * H0 * H1 + 2 * H1 **2 ) + R1 * H0**2

C     INITIALIZE  TRIDIAG BANDS A,B, C, AND COL. VECTOR YD
C     YD WILL LATER BE USED TO RETURN THE DERIVATIVES     
      DO I=1,LEN-2
         H0    =  X(I)   - X(I-1)
         H1    =  X(I+1) - X(I)
         R0    = (Y(I)   - Y(I-1)) / H0
         R1    = (Y(I+1) - Y(I))   / H1
         A(I)  = H1
         B(I)  = 2 * (H0 + H1)
         C(I)  = H0
         YD(I) = 3 * (R0 * H1 + R1 * H0)
      ENDDO

C     LAST ROW 
      A(I)  = (H0 + H1) **2
      B(I)  = H0 * (H0 + H1)
      YD(I) = R0 * H1**2 + R1 * (3 * H0 * H1 + 2 * H0**2)

C     SOLVE FOR THE TRIDIAGONAL  MATRIX: YD = YD * INV(TRIDIAG MATRIX) 
      CALL TRIDIAG_GEN(A, B, C, YD, LEN,IRTFLG)

      END

C     ********************** TRIDIAG_GEN ***************************

      SUBROUTINE TRIDIAG_GEN(A, BB, C, D , LEN, IRTFLG)

      DOUBLE PRECISION, DIMENSION(0:LEN-1) :: A,BB,C,D

      DOUBLE PRECISION                     :: B

C     AUTOMATIC ARRAYS
      DOUBLE PRECISION, DIMENSION(0:LEN-1) :: F

C     GAUSS ELIMINATION; FORWARD SUBSTITUTIION
      B    = BB(0)
      D(0) = D(0) / B

      DO I = 1,LEN-1
         F(I) = C(I-1) / B
         B    = BB(I) - A(I) * F(I)
         IF (B .EQ. 0.0) THEN
            CALL ERRT(101,'TRIDIAG_GEN: DIVIDE BY ZERO',NE)
            IRTFLG = 1
            RETURN
         ENDIF
         D(I) = (D(I) - D(I-1) * A(I)) / B
      ENDDO

C     BACKSUBSTITUTION 
      DO I=LEN-2,0,-1
         D(I) = D(I) - (D(I+1) * F(I+1))
      ENDDO
      IRTFLG = 0

      END

