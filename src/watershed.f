C ++********************************************************************
C                                                                      *
C   WATERSHED           K WAS WRONG                MAR 01 ARDEAN LEITH *
C                                                                      *
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
C                                                                      *
C  WATERSHED(LUN1,LUN2,NX,NY,NZ,FMINT)                         *
C                                                                      *
C  PURPOSE: WATERSHED                                                  *
C                                                                      *
C  PARAMETERS:                                                         *
C                                                                      *
C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************

        SUBROUTINE WATERSHED(LUN1,LUN2,NX,NY,NZ,FMINT)

        INCLUDE 'CMBLOCK.INC'
 
        REAL, ALLOCATABLE   :: VOL(:),VOLOUT(:)

        NPIX = NX * NY * NZ

        ALLOCATE (VOL(NPIX),VOLOUT(NPIX),STAT=IRTFLG)

        IF (IRTFLG .NE. 0) THEN 
           MWANT = NPIX * 2
           CALL ERRT(46,'WATERSHED; VOL..',MWANT)
           RETURN
        ENDIF

        IF (NZ == 1)  THEN
C          FOR IMAGE

C          READ INPUT VOLUME
           CALL REDVOL(LUN1,NX,NY,1,1,VOL,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999 

           CALL WATER2(VOL,VOLOUT,NX,NY,FMINT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999 

C          CREATE OUTPUT VOLUME
           CALL WRTVOL(LUN2,NX,NY,1,1,VOLOUT,IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999 


        ELSE
C          FOR VOLUME
           CALL ERRT(101,'NOT IMPLEMENT FOR VOLUMES',NE)
           GOTO 9999

        ENDIF

9999    IF (ALLOCATED(VOL))    DEALLOCATE(VOL)
        IF (ALLOCATED(VOLOUT)) DEALLOCATE(VOLOUT)
       
        END

C       ----------------- WATER2 -------------------------------------

        SUBROUTINE WATER2(XIMG,XOUT,NX,NY,FMINT,IRTFLG)

        REAL              :: XIMG(NX*NY),XOUT(NX*NY)

        REAL, ALLOCATABLE :: XSORT(:),RLOC(:)
        REAL, ALLOCATABLE :: NEXTPIX(:)

        INCLUDE 'CMBLOCK.INC'
 
        IRTFLG     = 1
        NPIX       = NX * NY

        NQUESIZE = NPIX

        ALLOCATE (XSORT(NPIX),RLOC(NPIX),NEXTPIX(NQUESIZE), 
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN 
           MWANT = NPIX * 2 + NQUESIZE
           CALL ERRT(46,'WATER2, RLOC,... ',MWANT)
           RETURN
        ENDIF

        DO IPIX = 1,NPIX
           XSORT(IPIX) = XIMG(IPIX)
           RLOC(IPIX)  = IPIX
        ENDDO

C       SORT THE PIXELS (SORT IS A SINGLETON SORT), RLOC IS REORDERD
        CALL SORT(XSORT,XOUT,RLOC,NPIX)

C       SET CURRENT SEGMENTATION LEVEL
        RLEVEL     = 1.0
        NLEVEL     = 0

C       INITIALIZE ARRAY XOUT
        XOUT   = -1.0

        DO IPIX = 1,NPIX             ! LOOP OVER ALL SORTED PIXELS
           ILOC = RLOC(IPIX)         ! STARTING PIXEL LOCATION

C          UPDATE LEVEL
           IF (NLEVEL > 0) THEN
              RLEVEL = RLEVEL + 1.0
              NLEVEL = 0

              DO I=1,NPIX
                 IF(XOUT(I) == 0) XOUT(I) = -1
              ENDDO
           ENDIF

C          PUT STARTING PIXEL IN QUE
           IQUE          = 0
           NQUE          = 1
           NEXTPIX(NQUE) = ILOC
C          write(6,*)' ----------- stack( ',nstack,'): ',ILOC

           DO WHILE(IQUE < NQUE) 

C             GET NEXT PIXEL FROM QUE
              IQUE  = IQUE + 1
              ISLOC = NEXTPIX(IQUE)
c             write(6,*)' from stack( ',istack,') got: ',ISLOC

              IF (XOUT(ISLOC) <= 0.0) THEN
C                PIXEL NOT VISITED OR SET YET 

C                GIVE UNVISITED PIXEL CURRENT WATERSHED LEVEL
                 XOUT(ISLOC) = RLEVEL

c                write(6,*)' XOUT(',ISLOC,'): ',RLEVEL

C                INCREMENT NUMBER OF PIXELS IN AT THIS LEVEL
                 NLEVEL = NLEVEL + 1

C                PLACE PIXEL'S UNVISITED HIGHER NEIGHBORS IN QUE

C                UPPER
                 NEIGH = ISLOC - NX
                 IF (NEIGH > 0 .AND. 
     &                XOUT(NEIGH) < 0.0 .AND.
     &               (XIMG(NEIGH) >=  XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,' QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                LOWER
                 NEIGH = ISLOC + NX
                 IF (NEIGH <= NPIX .AND. 
     &                XOUT(NEIGH) < 0.0 .AND.
     &               (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                LEFT
                 NEIGH = ISLOC - 1
                 IF (MOD(NEIGH,NX) > 0 .AND. 
     &               XOUT(NEIGH) < 0 .AND.
     &               (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                RIGHT
                 NEIGH = ISLOC + 1 
                 IF (MOD(ISLOC,NX) > 0 .AND.
     &               XOUT(NEIGH) < 0 .AND.
     &              (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                UPPER-LEFT
                 NEIGH = ISLOC - NX - 1 
                 IF (NEIGH > 0 .AND.
     &               MOD(NEIGH,NX) > 0 .AND.
     &               XOUT(NEIGH) < 0 .AND.
     &              (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                LOWER-LEFT
                 NEIGH = ISLOC + NX - 1 
                 IF (NEIGH <= NPIX .AND.
     &               MOD(NEIGH,NX) > 0 .AND.
     &               XOUT(NEIGH) < 0 .AND.
     &              (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                UPPER-RIGHT
                 NEIGH = ISLOC - NX + 1 
                 IF (NEIGH > 0 .AND.
     &               MOD(ISLOC,NX) > 0 .AND.
     &               XOUT(NEIGH)   < 0 .AND.
     &              (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

C                LOWER-RIGHT
                 NEIGH = ISLOC + NX + 1 
                 IF (NEIGH <= NPIX .AND.
     &               MOD(ISLOC,NX) > 0 .AND.
     &               XOUT(NEIGH) < 0 .AND.
     &              (XIMG(NEIGH) >= XIMG(ISLOC))) THEN
                     NQUE = NQUE + 1
                     IF (NQUE > NQUESIZE) THEN
                        CALL ERRT(102,'QUE OVERFLOW',NQUE)
                        RETURN
                     ENDIF
                     NEXTPIX(NQUE) = NEIGH
                     XOUT(NEIGH)   = 0.0
                 ENDIF

              ENDIF  ! END OF: IF (XOUT(ISLOC) <= 0  .AND......
           ENDDO     ! END OF: DO WHILE (IQUE <= NQUE)

           NSHEDS = RLEVEL
           IF (NLEVEL > 0) THEN
              WRITE(NOUT,90)NSHEDS,NLEVEL
 90           FORMAT('  WATERSHED: ',I0,'  HAS: ',I0,' PIXELS')
           ENDIF

        ENDDO        ! END OF: IPIX = 1,NPIX  LOOP OVER SORTED PIXELS
       
        NSHEDS = RLEVEL
        IF (NLEVEL <= 0) NSHEDS = NSHEDS -1
        WRITE(NOUT,*)' ' 
        WRITE(NOUT,*)' WATERSHEDS: ',NSHEDS 

        IRTFLG = 0
        END
