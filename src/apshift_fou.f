
C++*********************************************************************
C                                                                      *
C APSHIFT_FOU.F     NEW                          AUG 2011 ARDEAN LEITH *
C                   QUADRI_COR FOR SPEED         NOV 2011 ARDEAN LEITH *
C                   FBS2                         DEC 2011 ARDEAN LEITH *
C                   FBS2  NXP                    DEC 2011 ARDEAN LEITH *
C                   DO_FFT_R SET                 APR 2012 ARDEAN LEITH *
C                                                                      *
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
C  APSHIFT_FOU(LUNIN, REFPAT,IMGREF,  NX,NY, NXP,NYP, 
C              EXPBUF,AVI,SIGI, ISHRANGE,
C              ROTANG,XSH,YSH,MIRRORIT,PEAKV,IRTFLG)
C                                                                       
C  PURPOSE:  ROTATE EXPBUF IMAGE BY ROTANG, MIRROR IF NECESSARY, 
C            CROSS-CORRELATE WITH REFERENCE IMAGE, FIND CC PEAK.
C                                                                       
C  PARAMETERS: LUNIN        IO UNIT                               SENT
C              EXPPAT       EXP. IMAGE TEMPLATE                   SENT
C              IMGEXP       EXP. IMAGE NUMBER                     SENT
C              REFPAT       REF. IMAGE TEMPLATE                   SENT
C              IMGREF       REF. IMAGE NUMBER (UNUSED)            SENT
C              NX,NY        ACTUAL INPUT IMAGE DIM.               SENT
C              EXPBUF       EXP. IMAGE BUFFER                     SENT  
C              AVI,SIGI     EXP. IMAGE STATISTICS                 SENT
C              AVR          REF. IMAGE AVERAGE (UNUSED)           SENT
C              SIGR         REF. IMAGE STATISTICS                 SENT
C              LMASK                                              SENT
C              BUFPADI      EXP IMAGE FFT                         WORK
C              BUFPADR      REF IMAGE FFT                         SENT
C              F0,X1,Y1,XY2 FBS ARRAYS                            WORK
C              ISHRANGE..   POSSIBLE IMAGE SHIFTS                 SENT
C              ROTANG       INPLANE ROTATION ANGLE                SENT
C              XSH,YSH      SHIFTS                                RET.
C              MIRRORIT     LOGICAL FLAG THAT REF. NEEDS MIRROR   SENT
C              PEAKV        PEAK HEIGHT                           RET. 
C              IRTFLG       ERROR FLAG                            RET.  
C                                                                       
C  NOTE:   NXP INCLUDES THE EXTRA SPACE FOR FOURIER XFORM
C                                                                       
C--*********************************************************************

        SUBROUTINE APSHIFT_FOU(LUNIN,  IMGREF,
     &                     NX,NY,      NXLD, 
     &                     EXPBUF,     AVI,SIGI, AVR,SIGR, LMASK,
     &                     BUFPADI,    BUFPADR,
     &                     F0,         X1,Y1,XY2, 
     &                     ISHRANGEX,  ISHRANGEY,
     &                     ROTANG,     XSH,YSH,
     &                     MIRRORIT,   PEAKV, IRTFLG)

        IMPLICIT NONE

        INCLUDE 'CMLIMIT.INC'
        INCLUDE 'CMBLOCK.INC'

        INTEGER,INTENT(IN)     :: LUNIN
        INTEGER,INTENT(IN)     :: IMGREF
        INTEGER,INTENT(IN)     :: NX,NY, NXLD
        REAL                   :: EXPBUF(NX,NY)
        REAL,   INTENT(IN)     :: AVI,SIGI
        REAL,   INTENT(IN)     :: avr,SIGR
        LOGICAL,INTENT(IN)     :: LMASK(NXLD,NY)
        INTEGER,INTENT(IN)     :: ISHRANGEX,ISHRANGEY
        REAL,INTENT(OUT)       :: BUFPADI(NXLD,NY)  
        REAL,INTENT(IN)        :: BUFPADR(NXLD,NY)    
        REAL                   :: F0(NXLD,NY), XY2(NXLD,NY)
        REAL                   :: X1(NXLD,NY), Y1 (NXLD,NY)

        REAL,   INTENT(INOUT)  :: ROTANG
        REAL,   INTENT(OUT)    :: XSH,YSH
        LOGICAL,INTENT(IN)     :: MIRRORIT
        REAL,   INTENT(OUT)    :: PEAKV
        INTEGER,INTENT(OUT)    :: IRTFLG

        CHARACTER(LEN=MAXNAM)  :: FILNAM
        INTEGER * 8            :: IPLAN = 0     !STRUCTURE POINTER 

        INTEGER                :: NLET,NSLICE,MAXIM,IYCEN,IXCEN
        INTEGER                :: IREP,IY,IGO,IX,ITYPE,NXP1
        INTEGER                :: NE,IBEST,IODD,MWANT,INV,IRAD,IRADSQ
        integer                :: NXt,NYt,itypet
                   
        REAL                   :: PADVAL,RNY1,RNX1,RNY2,RNX2
        REAL                   :: COD,SID,YI,YCOD,YSID,XOLD,YOLD
        REAL                   :: COST(2),SINT(2)
        REAL                   :: XI,ZDUM

        REAL                   :: fbs2,quadri,quadri_cor

        LOGICAL                :: DO_FFT_R     = .FALSE.  
        LOGICAL                :: DO_FFT_I     = .TRUE.
        LOGICAL                :: SPIDER_SCALE = .FALSE.
        LOGICAL                :: SPIDER_SIGN  = .FALSE.
        LOGICAL                :: NORMIT       = .TRUE.
        LOGICAL                :: SKIP_PEAK    = .FALSE.
        INTEGER                :: ILOCS(2)

        REAL                   :: PEAKVT(2), ROTANGT(2)
        REAL                   :: XSHT(2),   YSHT(2)
        DOUBLE PRECISION       :: DAVT(2)    ! UNUSED
        DOUBLE PRECISION       :: DSIGT(2)
        INTEGER                :: ILOC

        REAL, PARAMETER        :: QUADPI = 3.1415926535897932384626
        REAL, PARAMETER        :: DGR_TO_RAD = QUADPI / 180          

        real :: temp


         IF (USE_FBS_INTERP) THEN
C           PREPARE DATA FOR FSBI
            write(6,*) ' Using fbs_interp not implemented'
            stop

            F0(1:NX, 1:NY) = EXPBUF(1:NX, 1:NY)  ! PADDING
            CALL FBS2_PREP(F0, X1,Y1,XY2, NXLD,NX,NY, IRTFLG)
            IF (IRTFLG .NE. 0) GOTO 9999
         ENDIF

         IRAD       = (MIN(NX,NY)) / 2 - 1
         IRADSQ     = IRAD**2

         NXP1        =  NX + 1

C        SPIDER IMAGE CENTER
         IYCEN       =  NY/2+1
         IXCEN       =  NX/2+1

C        IMAGE DIMENSIONS AROUND ORIGIN
         RNY1         = -NY/2
         RNX1         = -NX/2

         RNY2         = -RNY1
         RNX2         = -RNX1
         IF (MOD(NX,2) == 0) RNX2 = RNX2 - 1.0
         IF (MOD(NY,2) == 0) RNY2 = RNY2 - 1.0

         ROTANGT(1) = ROTANG
         ROTANGT(2) = ROTANG + 180.0

C        CREATE TRANSFORMATION MATRIX

         SINT(1)     = SIN(ROTANG * DGR_TO_RAD)
         !COST(1)    = COS(ROTANG * DGR_TO_RAD)
         !temp       = SQRT(1 - SINT(1)**2)  !faster??
         IF (ROTANG > 90 ) THEN
             COST(1) = -SQRT(1 - SINT(1)**2)  
         ELSE
             COST(1) =  SQRT(1 - SINT(1)**2)  
         ENDIF 

         SINT(2)     = -SINT(1)
         COST(2)     = -COST(1)

#ifdef NEVER
         SINT(1)     = SIN(ROTANGT(1) * DGR_TO_RAD)
         COST(1)     = COS(ROTANGT(1) * DGR_TO_RAD)
         SINT(2)     = SIN(ROTANGT(2) * DGR_TO_RAD)
         COST(2)     = COS(ROTANGT(2) * DGR_TO_RAD)
         write(6,*)'  ang,sin,cos 1:', ROTANGT(1),sint(1),COST(1)
         write(6,*)'  ang,sin,cos 2:', ROTANGT(2),sint(2),COST(2)
#endif

         IF (     (1 - IYCEN) < RNY1) THEN
            CALL ERRT(102,'PGM ERR, BAD (1 - IYCEN)',(1 - IYCEN))
            STOP
         ELSEIF ( (1 - IYCEN) > RNY2) THEN
            CALL ERRT(102,'PGM ERR, BAD (1 - IYCEN)',(1 - IYCEN))
            STOP
         ELSEIF ( (1 - IXCEN) < RNX1) THEN
            CALL ERRT(102,'PGM ERR, BAD (1 - IXCEN)',(1 - IXCEN))
            STOP
         ELSEIF ( (1 - IXCEN) > RNX2) THEN
            CALL ERRT(102,'PGM ERR, BAD (1 - IXCEN)',(1 - IXCEN))
            STOP
         ENDIF
         !call chkfile('jnkexpbufa',66,1,nx,ny,1, expbuf,irtflg)

         DO IREP = 1,2

            !bufpadi(nx+1:nxld,1:ny) = 0.0   ! ZERO PAD
            !call chkcmplx('rings',bufpadr,nxld*NY, nxld,4,NY*nxld)

C           ROTATE EXP IMAGE BY ROTANG & PAD FOR FFT  ------ ROTATE

C           CURRENT TRANSFORMATION MATRIX
            COD  = COST(IREP)
            SID  = SINT(IREP)

C           OMP ON OUTER LOOP MUCH FASTER!
c$omp       parallel do private(iy,yi,ycod,ysid,
c$omp&                          ix,xi,xold,yold,iloc)
            DO IY=1,NY

               YI   = IY - IYCEN

               YCOD =  YI * COD + IYCEN
               YSID = -YI * SID + IXCEN

               DO IX=1,NX
                  XI   = IX - IXCEN
 
                  XOLD = XI * COD  + YSID
                  YOLD = XI * SID  + YCOD

                  ILOC = IX
                  IF (MIRRORIT) ILOC = NXP1 - IX

                 !write(6,*)'  ILOC,IY zero:', ILOC,IY,xold,yold 

                  IF ( XOLD < 2 .OR. XOLD >= NX .OR. 
     &                 YOLD < 2 .OR. YOLD >= NY) THEN
C                      OUTSIDE OF UNROTATED IMAGE
                     BUFPADI(ILOC,IY) = AVI    !AVI
                     !write(6,*)'  ILOC,IY zero:', ILOC,IY,xold,yold 

                     CYCLE
                  ELSEIF (LMASK(ILOC,IY)) THEN
                     BUFPADI(ILOC,IY) = AVI    !AVI
                     CYCLE
 
                  ENDIF
                  BUFPADI(ILOC,IY) = 
     &                  QUADRI_COR(XOLD,YOLD,NX,NY,EXPBUF)
	       ENDDO
	    ENDDO

C           SINCE IMAGE IS ROTATED AND SHIFTED STATS MAY CHANGE.
C           SLOWER IF NORM CODE IS EXTRACTED AND PUT IN ABOVE LOOP!!
            CALL NORMVALSP(BUFPADI,NX,NY,1,
     &                     NXLD,NY,1, 
     &                     DAVT(IREP),DSIGT(IREP),.TRUE.)
 
            !if (irep == 1)  then
            !bufpadi(nx+1:nxld,1:ny) = 0.0   ! ZERO PAD
            !call chkfile('jnkrtsq1',66,1,nxld,ny,1, bufpadi,irtflg)
            !call chkfile('jnkpadr', 66,1,nxld,ny,1, bufpadr,irtflg)
            !else
            !bufpadi(nx+1:nxld,1:ny) = 0.0   ! ZERO PAD
            !call chkfile('jnkrtsq2',66,1,nxld,ny,1, bufpadi,irtflg)
            !endif

            !write(6,*)' Calling apcc_new, irep: ',irep
C           FFT CROSS CORRELATION BUFPADI vs BUFPADR --------------- CC
            CALL APCC_NEW(NXLD,  NX,NY,1, 
     &                BUFPADI,   BUFPADR,
     &                DO_FFT_I,  DO_FFT_R,
     &                SKIP_PEAK, NORMIT,SPIDER_SIGN, 
     &                ISHRANGEX, ISHRANGEY,0,
     &                XSHT(IREP),YSHT(IREP),ZDUM, 
     &                PEAKVT(IREP),IRTFLG)
            IF (IRTFLG .NE. 0)  THEN
               CALL ERRT(101,'APSHIFT_FOU CC ERROR',NE)
               RETURN
            ENDIF

            !if (irep == 1)  then
            !call chkfile('jnkcc1',66,1,nxld,ny,1, bufpadi,irtflg)
            !else
            !call chkfile('jnkcc2',66,1,nxld,ny,1, bufpadi,irtflg)
            !endif

        ENDDO

C       NORMALIZATION 
C       NO SCALING IN FMRS OR IN CCRS, IS SCALED BY 1/NX*NY in APCC
 
        PEAKVT(1) = PEAKVT(1) / FLOAT(NX*NY-1) / DSIGT(1) / SIGR
        PEAKVT(2) = PEAKVT(2) / FLOAT(NX*NY-1) / DSIGT(2) / SIGR

        IBEST = 1
        IF (PEAKVT(2) > PEAKVT(1)) IBEST = 2

        PEAKV  = PEAKVT  (IBEST)
        ROTANG = ROTANGT(IBEST)
        XSH    = XSHT (IBEST)
        YSH    = YSHT (IBEST)
      
        !write(6,*)' 1:',rotangt(1),xsht(1),ysht(1),peakvt(1)
        !write(6,*)' 2:',rotangt(2),xsht(2),ysht(2),peakvt(2)

C       ADJUST FOR MIRROR OF IMAGE NOT MIRROR OF REF. IMAGE
        IF (MIRRORIT) THEN              ! BUG FIX MAR 2011
           IODD = MOD(NX,2)
           IF (IODD == 0) THEN
              XSH = 1 - XSH
           ELSE
              XSH = -XSH
           ENDIF
        ENDIF            
    
9999    CONTINUE

        END


C     ------------------- QUADRI_COR ---------------------------------

      REAL FUNCTION QUADRI_COR(X,Y, NX,NY, FDATA)

C     MUST HAVE ALREADY HANDLED  CORNER VALUES IN CALLER AND 
C     ENSURED THAT X & Y ARE INSIDE IMAGE!!
    
      IMPLICIT NONE

      REAL      :: X,Y           ! INPUT LOCATION
      INTEGER   :: NX,NY         ! IMAGE DIMENSIONS
      REAL      :: FDATA(NX,NY)  ! INPUT IMAGE

      INTEGER   :: I,J, IM1,IP1,JM1,JP1
      REAL      :: DX0,DY0,DXB,DYB,F0,C1,C2,C3,C4,C5

      I   = IFIX(X)
      J   = IFIX(Y)

      DX0 = X - I          ! SLOWER IF MOVED DOWN!!
      DY0 = Y - J
      DXB = DX0 - 1.0
      DYB = DY0 - 1.0

      IM1 = I-1            ! USE SPEEDS UP!
      IP1 = I+1
      JM1 = J-1
      JP1 = J+1

      F0  = FDATA(I,J)

      C1  = FDATA(IP1,J) - F0

      C2  = (C1 - F0 + FDATA(IM1,J)) * 0.5

      C3  = FDATA(I,JP1) - F0 

      C4  = (C3 - F0 + FDATA(I,JM1)) * 0.5 

      IF (DX0 < 0) THEN
         IF (DY0 < 0) THEN   ! -1 -1
            C5 =  (FDATA(IM1,JM1) - F0 + C1 - 2 * (C2 + C4) + C3)

         ELSE                ! -1 +1
            C5 =  (FDATA(IM1,JP1) + F0 - C1 + 2 * C2 + C3)  

         ENDIF
      ELSE
         IF (DY0 < 0) THEN   ! +1 -1
            C5 =  (FDATA(IP1,JM1) + F0 + C1 - C3 + 2 * C4) 
 
         ELSE                ! +1 +1
            C5 =  (FDATA(IP1,JP1) - F0 - C1 - C3)
         ENDIF
      ENDIF

      QUADRI_COR = F0 + 
     &              DX0 * (C1 + DXB * C2 + DY0 * C5) + 
     &              DY0 * (C3 + DYB * C4)

      END











#ifdef DEBUGNEVER
               IF (USE_FBS_INTERP) THEN

c$omp            parallel do private(IX,xi,xold,yold)
                 DO IX=1,NX

                    !IRADNOW = IY**2 + IX**2
                    !IF (IRADNOW > IRADSQ) THEN
C                   !   OUTSIDE OF ROTATIONAL MASK
                    !   BUFPADI(IX,IY)
                    !ENDIF
                    XI = IX - IXCEN 

                    IF (XI .LT. RNX1) XI = MIN(RNX2+XI-RNX1+1.0,RNX2)
                    IF (XI .GT. RNX2) XI = MAX(RNX1+XI-RNX2-1.0,RNX1)

                    YOLD = XI * SID  + YCOD
                    XOLD = XI * COD  + YSID

                    IF (.NOT. MIRRORIT) THEN
C                      NO MIRROR 
                       BUFPADI(IX,IY) = FBS2(XOLD,YOLD,
     &                      NXLD,NX,NY,EXPBUF,NX,X1,Y1,XY2,.FALSE.)
                    ELSE
C                      MIRROR THE IMAGE IN PLACE  
                       BUFPADI(NXP1-IX,IY) = FBS2(XOLD,YOLD,
     &                      NXLD,NX,NY, EXPBUF,NX, X1,Y1,XY2, .FALSE.)
                    ENDIF
	         ENDDO

              ELSE
#endif


