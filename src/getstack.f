
C++*********************************************************************
C
C GETSTACK.F  -- CREATED 10 SEP 02 ArDean Leith
C                USED GETOLDSTACK, GETNEWSTACK    APR 99  ArDean Leith
C                INDEXED STACK                    JAN 02  ArDean Leith
C                GETNEWSTACK PARAM.               FEB 03  ArDean Leith
C                MRC SUPPORT                      SEP 19  ArDean Leith
C
C **********************************************************************
C=*                                                                    *
C=* This file is part of:   SPIDER - Modular Image Processing System.  *
C=* SPIDER System Authors:  Joachim Frank & ArDean Leith               *
C=* Copyright 1985-2019  Health Research Inc.,                         *
C=* Riverview Center, 150 Broadway, Suite 560, Menands, NY 12204.      *
C=* Email: spider@health.ny.gov                                        *
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
C  GETSTACK(LUN1,LUN2,IMGNUM,MAXIM,SAYIT,LOADIT,BUF,NORMIT,IRTFLG)
C
C  PURPOSE:  LOCATE NEXT IMAGE/VOL. IN INPUT & OUTPUT STACK & 
C            OPTIONALLY LOAD IT INTO: BUF.
C 
C  PARAMETERS:  LUN1      LUN FOR INPUT IMAGE & COPYING FMIN...  (SENT)            (SENT)
C               LUN2      LUN FOR OUTPUT IMAGE                   (SENT) 
C               IMGNUM    LAST IMAGE                        (SENT/RET.)
C                 USE -3 ON INPUT FOR NON STACK LOOP END         (SENT)
C               MAXIM     STACK FLAG                             (SENT)
C               SAYIT     REPORT FILE OPENING INFO               (SENT)
C               LOADIT    LOAD IMAGE/VOL. IN BUF                 (SENT)
C               BUF       IMAGE/VOLUME AREA                      (RET.)
C               NORMIT    FIND IMAGE/VOL. MIN,MAX....            (SENT)
C               IRTFLG    ERROR FLAG                             (RET.)
C
C  NOTE:      THE THREE SUBROUTINES: GETSTACK, GETOLDSTACK, AND
C             GETNEWSTACK WERE WRITTEN TO HANDLE WHOLE STACK (@STACK)
C             OPERATIONS. THEY ARE ONLY USED IN THE OPERATIONS:
C             'AF' (UNDOCUMENTED), 'CP TO STK','DE','CE AD', 
C             'CE' (MANY?), 'AR SCA',
C             'AD D2', 'AD DF', 'AD DR', 'AD F', 'AD R', 'DIV 2', 
C             'MUL 2', 'MUL O', 'MU O', & 'SU B2'.
C             I DECIDED THAT WHOLE STACK OPERATION SUPPORT WAS LESS
C             IMPORTANT THAN SUPPORTING OPERATIONS ON SELECTED 
C             STACKED IMAGES SO I STOPPED CONVERTING OPERATIONS TO
C             HANDLING WHOLE STACKS WITHOUT SELECTED STACK SUPPORT.
C             THE NEW PARADIGM USES OPTILES AND NEXTFILE WHICH
C             HANDLES BOTH SELECTED AND WHOLE STACKS. 
c
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************

       SUBROUTINE GETSTACK(LUN1,LUN2,IMGNUM,MAXIM,
     &                     SAYIT,LOADIT,BUF,NORMIT,IRTFLG)

       IMPLICIT NONE
       INCLUDE    'CMBLOCK.INC' 

       INTEGER  :: LUN1,LUN2,IMGNUM,MAXIM,IRTFLG
       LOGICAL  :: WANTNEXT,MUSTGET,SAYIT,LOADIT,NORMIT
       REAL     :: BUF(*)

       INTEGER  :: NX,NY,NZ
       LOGICAL  :: IS_BARE,IS_MRC
       LOGICAL  :: COPYSTAT = .FALSE.  ! COPY FMIN.. FROM INPUT FILE

       WANTNEXT = .TRUE.
       MUSTGET  = .FALSE.

C      GET CURRENT NX,NY,NZ
       CALL LUNGETSIZE(LUN1,NX,NY,NZ,IRTFLG)   ! MRC OK
       IF (IRTFLG .NE. 0) RETURN

       CALL LUNGETIS_MRC(LUN1,IS_MRC,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       CALL LUNGETISBARE(LUN1,IS_BARE,IRTFLG)
       IF (IRTFLG .NE. 0) RETURN

       !write(6,*)' In getstack - maxim,bare,mrc: ',maxim,is_bare,is_mrc

       IF (MAXIM < 0  .OR. (.NOT. IS_BARE .AND. IS_MRC) ) THEN
C         NOT A WHOLE STACK OPERATION
          IMGNUM = MAXIM 

       ELSE
C         NEGATIVE IMGNUM SET TO ZERO
          IMGNUM = MAX(IMGNUM+1,1)

C         GET INPUT IMAGE FROM INPUT STACK            
          CALL GETOLDSTACK(LUN1,NX,IMGNUM,
     &                     WANTNEXT,MUSTGET,SAYIT,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

          !write(6,*) ' In getstack, imgnum 1: ',imgnum,irtflg

C         CREATE OUTPUT IMAGE IN OUTPUT STACK
C         COPYSTAT   COPY FMIN.. FROM INPUT FILE                
          CALL GETNEWSTACK(LUN1,LUN2,COPYSTAT,NX,IMGNUM,IRTFLG)
          IF (IRTFLG .NE. 0) RETURN

       ENDIF

C      LOAD CURRENT IMAMI..  IN COMMON AREA
       CALL LUNGETSTAT(LUN1,IMAMI,FMIN,FMAX,AV,SIG,IRTFLG) ! MRC OK

       !write(6,*)' In getstack - imami, max,min: ',imami,fmax,fmin

       IF (IMAMI .NE. 1 .AND. NORMIT) THEN
          CALL NORM3(LUN1,NX,NY,NZ,FMAX,FMIN,AV)
          IF (IRTFLG .NE. 0) RETURN
       ENDIF

       IF (LOADIT) THEN
C         LOAD IMAGE/VOL. IN BUF
          CALL REDVOL(LUN1,NX,NY,1,NZ,BUF,IRTFLG)
       ENDIF

       END
