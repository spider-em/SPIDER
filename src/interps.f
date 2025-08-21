C++*********************************************************************
C
C INTERPS.F         FROM UTIL2                    ARDEAN LEITH 11/24/10
C                   OPFILES                       ARDEAN LEITH  1/11/11
C 	            IP FS --- FOURIER SPLINE      ARDEAN LEITH  5/23/11
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2011	  Health Research Inc.,                        *
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
C  INTERPS(MAXDIM)
C
C  PURPOSE: CALLER FOR INTERPOLATION OPERATIONS
C	    IP         BILINEAR                                       
C	    IP T       TRIANGULAR 
C	    IP FT      FOURIER    
C 	    IP FS      FOURIER SPLINE 
C
C23456789012345678901234567890123456789012345678901234567890123456789012
C--*********************************************************************

        SUBROUTINE INTERPS(MAXDIM)

        IMPLICIT NONE

	INCLUDE 'CMBLOCK.INC'
	INCLUDE 'CMLIMIT.INC'

        INTEGER                 :: MAXDIM

        CHARACTER(LEN=MAXNAM)   :: FILNAM1,FILNAM2
        CHARACTER(LEN=MAXNAM)   :: PROMPT
	CHARACTER(LEN=1)        :: NULL = CHAR(0)

        INTEGER, ALLOCATABLE    :: ILIST1(:),ILIST2(:)
        INTEGER                 :: NSHUD

        REAL, ALLOCATABLE       :: X1F(:), X2F(:) ! ACTUAL MULTI-DIM.     

	INTEGER, PARAMETER      :: LUN1   = 21
	INTEGER, PARAMETER      :: LUN2   = 22
	INTEGER, PARAMETER      :: LUNDOC = 70
	INTEGER, PARAMETER      :: LUNXM1 = 71
	INTEGER, PARAMETER      :: LUNXM2 = 72

        INTEGER                 :: NILMAX,MAXIM1,NLET,IFORM1,IFORM2
        INTEGER                 :: NX1,NX2,NY1,NY2,NDUM
        INTEGER                 :: NZ1,NZ2,NGOT1,NGOT2,IMG1
        INTEGER                 :: NOT_USED,IRTFLG,MAXIM2,NLET2,IMG2
        INTEGER                 :: NINDX1,NINDX2
        INTEGER                 :: NE,LSD,LSDN,MWANT1,MWANT2,IDUM

        REAL                    :: SCALEX,SCALEY,SCALEZ,SCALET

        NILMAX = NIMAXPLUS            ! FROM CMLIMIT
        ALLOCATE(ILIST1(NILMAX),
     &           ILIST2(NILMAX),
     &           STAT=IRTFLG)
        IF (IRTFLG .NE. 0) THEN
           CALL ERRT(46,'INTERPS; ILIST....',2*NILMAX)
           RETURN
        ENDIF

C       OPEN FIRST INPUT FILE
        MAXIM1  = 0
        IMG1 = 0    ! Passing uninitialized variable may result in 'INVALID IMAGE NUMBER' error
        CALL OPFILES(0,LUN1,LUNDOC,LUNXM1,  
     &               .TRUE.,FILNAM1,NLET, 'O',
     &               IFORM1,NX1,NY1,NZ1,MAXIM1,
     &               NULL,
     &              .FALSE., ILIST1,NILMAX, 
     &               NDUM,NGOT1,IMG1, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN
        !write(6,'(a,8i5)') ' Input maxim1,ngot1,img1:',maxim1,ngot1,img1
        !write(6,'(a,8i5)')' xxxInput maxim1,ngot1,img1:',maxim1,ngot1,img1

       
        IF (FCHAR(4:4) == 'T' .AND. NZ1 > 1) THEN
C          TRIANGULAR INTERPOLATION
           CALL ERRT(101,'OPERATION DOES NOT WORK ON VOLUMES',NE)
           GOTO 9999
        ENDIF

C	GET OUTPUT FILE NAME
C                 123456789 123456789 123456789 123456789 1234567890
        PROMPT = 'OUTPUT FILE NAME OR TEMPLATE (E.G. STK@****)'
        CALL FILELIST(.TRUE.,LUNDOC,FILNAM2,NLET2,
     &                ILIST2,NILMAX,NGOT2,PROMPT(1:44),IRTFLG)
        IF (IRTFLG .NE. 0) GOTO 9999

        NX2 = 0
        NY2 = 0
        IF (NZ1 > 1) THEN
           NZ2 = -1
           CALL RDPRI3S(NX2,NY2,NZ2,NOT_USED,
     &                'OUTPUT DIMENSIONS, NX, NY, & NZ',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999

           IF (NZ2 <= -1) THEN
              CALL RDPRI1S(NZ2,NOT_USED,'NZ',IRTFLG)
              IF (IRTFLG .NE. 0) GOTO 9999
           ENDIF

           IF (NZ2 == 0) THEN
              NZ2 = (FLOAT(NX2) / FLOAT(NX1)) * FLOAT(NZ2)
           ENDIF
        ELSE
           NZ2 = 1
           CALL RDPRIS(NX2,NY2,NOT_USED,
     &                'OUTPUT DIMENSIONS, NX & NY',IRTFLG)
           IF (IRTFLG .NE. 0) GOTO 9999
        ENDIF

C	FOR RECTANGULAR IMAGES, THE USER IS ALLOWED TO ENTER ONLY
C	ONE DIMENSION, THE OTHER DIM. IS COMPUTED TO KEEP THE SAME 
C	RELATION AS THE ONE BETWEEN THE DIMS. OF INPUT.

        IF (NY2 .EQ. 0) THEN
           NY2 = (FLOAT(NX2) / FLOAT(NX1)) * FLOAT(NY1)
        ENDIF

      
C       OPEN FIRST OUTPUT FILE NOW
        IMG2 = IMG1
        CALL OPFILES(LUN1,LUN2,LUNDOC,LUNXM2, 
     &             .FALSE.,FILNAM2,NLET,'N',
     &             IFORM1,NX2,NY2,NZ2,MAXIM2,
     &             FILNAM2(1:NLET2),
     &             .FALSE., ILIST2,NILMAX, 
     &             NDUM,NGOT2,IMG2, IRTFLG) 
        IF (IRTFLG .NE. 0) RETURN
        !write(6,'(a,8i5)')' xxoutput maxim1,ngot1,img1:',maxim1,ngot1,img1

c        write(6,'(A,4i6)') 
c     &        ' Output maxim2,ngot2,img2:',maxim2,ngot2,img2

        SCALEX = FLOAT(NX1) / FLOAT(NX2)
        SCALEY = FLOAT(NY1) / FLOAT(NY2)
        SCALEZ = FLOAT(NZ1) / FLOAT(NZ2)

C       IMAGE FILE HEADER FOR PIXSIZ HAS CHANGED
        SCALET = SCALEX
        IF (SCALEY .NE. SCALEX) SCALET = 0.0   ! X NOT SAME AS Y
        IF (NZ1 > 1 .AND. 
     &      SCALEZ .NE. SCALEX) SCALET = 0.0 ! Z NOT SAME AS Y
 
        LSD    = NX1+2-MOD(NX1,2)
        LSDN   = NX2+2-MOD(NX2,2)

        IF (FCHAR(4:4) .EQ. 'F') THEN
           MWANT1 = LSD *NY1*NZ1 
           MWANT2 = LSDN*NY2*NZ2 ! FOR FBS THIS ADDS UNUSED PAD

           ALLOCATE (X1F(MWANT1),
     &               X2F(MWANT2), STAT=IRTFLG)
           IF (IRTFLG .NE. 0) THEN 
               CALL ERRT(46,'INTERPS, X1F & X2F',MWANT1+MWANT2)
               GOTO 9999
           ENDIF 
        ENDIF 


        NINDX1 = 1
        NINDX2 = 1
        DO 
           IF (FCHAR(4:5) == 'FT') THEN
C             FOURIER INTERPOLATION -------------------------

              WRITE(NOUT,*) ' Using fourier interpolation'
              CALL READV(LUN1,X1F, LSD,NY1, NX1,NY1,NZ1)

              IF (NZ1 == 1) THEN
                 CALL FINT(X1F,X2F,NX1,NY1,NX2,NY2,LSD,LSDN)
              ELSE
C                VOLUME
                 CALL FINT3(X1F,X2F,NX1,NY1,NZ1,NX2,NY2,
     &                      NZ2,LSD,LSDN)
              ENDIF 

              CALL WRITEV(LUN2,X2F, LSDN,NY2, NX2,NY2,NZ2)

           ELSEIF (FCHAR(4:5) == 'FS') THEN
C             FOURIER BASED SPLINE INTERPOLATION -------------

              WRITE(NOUT,*) ' Using fourier based spline interpolation'
              CALL READV(LUN1,X1F, LSD,NY1, NX1,NY1,NZ1) ! PAD FOR FFT

              IF (NZ1 .EQ. 1) THEN
C                X2F IS NOT FOURIER PADDED
                 CALL INTERP_FBS(X1F, X2F, 
     &                           LSD,NX1,NY1,
     &                               NX2,NY2,  IRTFLG)
              ELSE
                 CALL INTERP_FBS3(X1F,X2F,NX1,NY1,NZ1,
     &                                    NX2,NY2,NZ2,  LSD)
              ENDIF 

C             X2F IS NOT FOURIER PADDED!
              CALL WRITEV(LUN2,X2F, NX2,NY2, NX2,NY2,NZ2)

           ELSEIF (FCHAR(4:4) == 'T') THEN
C             TRIANGULAR INTERPOLATION ----------------

              WRITE(NOUT,*) ' Using trianglular interpolation'

              CALL TRINTER(LUN1,LUN2,NX1,NY1,1,NX2,NY2,1)


           ELSE
C             BI-LINEAR INTERPOLATION ---------------------
              WRITE(NOUT,*) ' Using bilinear interpolation'

              IF (NZ1 .EQ. 1)  THEN
                 CALL INTERP(LUN1,LUN2,NX1,NY1,NX2,NY2,MAXDIM)
              ELSE
                 CALL INTERP3(LUN1,LUN2,NX1,NY1,NZ1,
     &                        NX2,NY2,NZ2,IDUM)
              ENDIF
           ENDIF

C          UPDATE THE INCORE HEADER VALUE & FILE HEADER FOR PIXSIZ
           CALL SETPRMS(LUN2, SCALET,IRTFLG)

C          OPEN NEXT SET OF I/O FILES, INCREMENTS NWANT 

c          write(6,'(A,4i6)') ' Next ------------------'
c          write(6,'(A,4i6)') 
c     &        ' Calling Nextfiles, nindx1,nindx2:',nindx1,nindx2

           CALL NEXTFILES(NINDX1,NINDX2,  ILIST1,ILIST2, 
     &                    .FALSE., LUNXM1,LUNXM2,
     &                    NGOT1,NGOT2,    MAXIM1,MAXIM2,  
     &                    LUN1,LUN1,LUN2, FILNAM1,FILNAM2,
     &                    IMG1,IMG2, IRTFLG) 
!          write(6,'(A,5i6)') 
!     &    ' Nextfiles img?,nindx?,eflg:',img1,img2,nindx1,nindx2,irtflg
           IF (IRTFLG .NE. 0) EXIT      ! ERROR / END OF INPUT STACK

        ENDDO

9999    CLOSE(LUN1)
        CLOSE(LUN2)
        CLOSE(LUNDOC)
        CLOSE(LUNXM1)
        CLOSE(LUNXM2)

        IF (ALLOCATED(ILIST1)) DEALLOCATE(ILIST1)
        IF (ALLOCATED(ILIST2)) DEALLOCATE(ILIST2)
        IF (ALLOCATED(X1F))    DEALLOCATE(X1F)
        IF (ALLOCATED(X2F))    DEALLOCATE(X2F)

        RETURN
        END
