

C++*********************************************************************
C
C  NEXTFILES.F  NEW                              12/15/06 ArDean Leith
C               OVERUN OUTPUT LIST = -99          1/15/12 ArDean Leith
C               SPLIT OUT FROM OPFILES           12/13/25 ArDean Leith
C               BARE STACK BUGS FIXED             2/18/26 ArDean Leith
C               OVERUN OUTPUT LIST BUG            2/18/26 ArDean Leith
C
C **********************************************************************
C 
C  CONTAINS BOTH : NEXTFILE & NEXTFILES
C
C  NEXTFILES(NINDX1, NINDX2, ILIST1,ILIST2, 
C            FOUROK,LUNXM1,LUNXM2,
C            NLIST1,NLIST2,   NSTK1,NSTK2,   
C            LUN1,LUNCP,LUN2, FILPAT1,FILPAT2,
C            ISTK1,ISTK2, IRTFLG) 
C
C PURPOSE:  GETS NEXT INPUT OR OUTPUT FILE FOR BARE STACK INPUT,
C           A NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST, 
C           OR A  STACKED IMAGE WITH/WITHOUT LIST
C           STACKS MUST BE OPENED WITH OPFILES!!!
C
C PARAMETERS: NINDX1,NINDX2    LIST INDICES                  (SENT/RET)
C             ILIST1,ILIST2    IMAGE NUMBER LISTS                (SENT)
C             FOUROK           FOURIER INPUT IS OK               (SENT)
C             LUNXM1,LUNXM2    LUN FOR SELFILE INPUT             (SENT)
C             NLIST1,NLIST2    NUMBER OF IMAGES                  (SENT)
C             NSTK1,NSTK2      MAX IMAGE IN STACK                (SENT)
C             LUN1             LUN FOR INPUT  (0 = NO FILE IN)   (SENT)
C             LUNCP            LUN FOR OUTPUT HEADER COPY        (SENT)
C             LUN2             LUN FOR OUTPUT (0 = NO FILE OUT)  (SENT)
C             FILPAT,FILPAT2   FILE NAME PATTERNS                (SENT)
C             ISTK1,ISTK2      IMAGE NUMBERS                 (SENT/RET)
C             IRTFLG           ERROR (0 IS OK, -1 IS END STACK)   (RET)
C
C CALLED BY:  bp32f.f, bp3f.f, copyd.f, copyfrommrc.f, 
C             copytomrc.f, copytostk.f, interps.f, rotqss.f,           
C             softmask.f,  util_1010.f, util_1011.f, util_1110.f,
C             util_11.f,        
C
C CALL TREE:
C     ONLY USED AFTER OPFILES HAS BEEN CALLED 
C
C       NEXTFILES 
C          |      
C          |      Old image (see opfiles.f)
C          ` ---> GETOLDIMG    
C          |      
C          |      New image (see opfiles.f)                         
C          ` ---> GETNEWIMG
C
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--********************************************************************* 
      SUBROUTINE NEXTFILES(NINDX1, NINDX2, ILIST1,ILIST2, 
     &                     FOUROK,LUNXM1,LUNXM2,
     &                     NLIST1,NLIST2,   NSTK1,NSTK2,   
     &                     LUN1,LUNCP,LUN2, FILPAT1,FILPAT2,
     &                     ISTK1,ISTK2, IRTFLG) 
 
      IMPLICIT NONE
      INCLUDE 'CMLIMIT.INC' 
 
      INTEGER           :: NINDX1,NINDX2
      INTEGER           :: ILIST1(NLIST1),ILIST2(NLIST2)
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1,LUNXM2
      INTEGER           :: NLIST1,NLIST2
      INTEGER           :: NSTK1,NSTK2,LUN1,LUNCP,LUN2
      CHARACTER(LEN=*)  :: FILPAT1,FILPAT2
      INTEGER           :: ISTK1,ISTK2,ISTK2_BAD,IRTFLG

      INTEGER           :: NWANT1,NWANT2, it
      INTEGER           :: NLETP1,NLETP2
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: GOTAST1, GOTAST2
      LOGICAL           :: GOTAT1,  GOTAT2
      LOGICAL           :: IS_BARE1,IS_BARE2    
      LOGICAL           :: IS_ONE_IMG1,IS_ONE_IMG2

      CHARACTER(LEN=MAXNAM) :: FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT 
      LOGICAL               :: IS_MRC, IS_MRCS,GOT_FILE2
      INTEGER               :: LIST_VAL



      INTEGER           :: lnblnkn       ! FUNCTION

C     INCREMENT CURRENT LIST POINTERS

      NINDX1 = NINDX1 + 1
      NINDX2 = NINDX2 + 1

      NLETP1 = lnblnkn(FILPAT1)
      NLETP2 = lnblnkn(FILPAT2)

      CALL LUNGETISBARE(LUN1,IS_BARE1,IRTFLG)

cc#define SP_DBUGIO
#if defined (SP_DBUGIO)
      write(3,*)' '
      write(3,*)' In Nextfiles-0; lun1:     ', lun1
      write(3,*)' In Nextfiles-0; lun2:     ', lun2
      write(3,*)' In Nextfiles-0; is_bare1:           ', is_bare1
      write(3,*)' In Nextfiles-0; nlist1:   ', nlist1
      write(3,*)' In Nextfiles-0; nlist2:   ', nlist2

      write(3,*)' In Nextfiles-0; filpat1:  ', filpat1(:nletp1)
      write(3,*)' In Nextfiles-0; filpat2:  ', filpat2(:nletp2)
      write(3,*)' '
#endif

      IF (.NOT. IS_BARE1 .AND.  (NLIST1 <= 0 .AND.  NLIST2 <= 0)) THEN
C         LISTS FINISHED 
          IRTFLG = -1

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In Nextfiles-0; irtflg:   ', irtflg
          write(3,*)' In Nextfiles-0; Returning at list end'
#endif
          RETURN
      ENDIF

#if defined (SP_DBUGIO)
      write(3,*)' In Nextfiles-0; nindx1,nindx2: ', nindx1,nindx2
      write(3,*)' In Nextfiles-0; nstk1,nstk2:   ', nstk1,nstk2
      write(3,*)' In Nextfiles-0; istk1,istk2:   ', istk1, istk2
      write(3,*)' '
#endif

         
      IF (LUN1 > 0) THEN  
C        OPEN NEXT INPUT IMG (MUST BE OLD & EXISTING, MAY BE IN STACK)
 
         GOTAT1  = (INDEX(FILPAT1,'@') > 0)
         GOTAST1 = (INDEX(FILPAT1,'*') > 0)

C        IS THIS A BARE STACK OPERATION?  (OK FOR SPIDER & MRC)
         CALL LUNGETISBARE(LUN1,IS_BARE1,IRTFLG)

C        DETERMINE IF A SINGLE IMG, (NO LIST, NOT BARE STACK) 
         IS_ONE_IMG1 = (GOTAT1 .AND. (.NOT. GOTAST1) .AND.
     &                               (.NOT. IS_BARE1))                


         IF (ISTK1 == -1 .AND. LUNXM1 > 0 ) THEN
C           XMIPP SELFILE LISTED IMAGE -----------------------------

            NWANT1 = -1

         ELSEIF (IS_ONE_IMG1 ) THEN
C           NON-STACK IMG/FILE (NO LIST) FILE ALREADY OPEN ----------

#if defined (SP_DBUGIO)
            write(3,*)' In Nextfiles-O1; filpat1:  ', filpat1(:nletp1)
            write(3,*)' In Nextfiles-O1; Can skip opening'
            write(3,*)' '

#endif
            GOTO 99

         ELSEIF ( IS_BARE1 ) THEN
C           MRC OR SPIDER BARE STACK INPUT (NO LIST) ---------------

            IF (ISTK1 < 0) ISTK1 = 1
            NWANT1 = ISTK1 + 1

#if defined (SP_DBUGIO)
          write(3,*)' '
          write(3,*)' In Nextfiles-OBS1; istk1:   ',  istk1
          write(3,*)' In Nextfiles-OBS1; nwant1:  ', nwant1
#endif

            IF (NWANT1 > NSTK1) THEN
C              FINISHED THE WHOLE STACK
               IRTFLG = -1
               RETURN
            ENDIF

         ELSEIF (.NOT. GOTAT1 .AND. GOTAST1 ) THEN
C           NON-STACKED IMG WITH TEMPLATED LIST --------------------

c            IF (NLIST1 > 0 .AND. NINDX1 > NLIST1) THEN
            IF (NINDX1 > NLIST1) THEN
C              OVERUN INPUT LIST
               IRTFLG = -1
               RETURN
            ENDIF

C           OPEN NEXT INPUT FILE 
            NWANT1 = ILIST1(NINDX1)

#if defined (SP_DBUGIO)
            write(3,*)' '
            write(3,*)' In Nextfiles-OTL1; nstk1:     ', nstk1 
            write(3,*)' In Nextfiles-OTL1; nindx1:    ', nindx1
            write(3,*)' In Nextfiles-OTL1; nlist1:    ', nlist1
            write(3,*)' In Nextfiles-OTL1; nwant1:    ', nwant1
#endif

         ELSEIF (NLIST1 == 0 .AND. NSTK1 == 0  .AND. ISTK1 == 0) THEN

C           NO INPUT LISTS, INPUT NOT IN STACK. NOT A FILELIST ----

#if defined (SP_DBUGIO)
            write(3,*)' In nextfiles-1; No Lists. Finished  '
#endif
            IRTFLG = -1
            RETURN

         ELSEIF (NSTK1 == -2  .OR.
     &           NSTK1 == -1  .OR.
     &           NSTK1  >    0 ) THEN
C           NON-STACKED IMAGE WITH/WITHOUT TEMPLATED LIST ----------
C           STACKED     IMAGE WITH/WITHOUT LIST         

#if defined (SP_DBUGIO)
            write(3,*)' '
            write(3,*)' In Nextfiles_ONL; nstk1 :      ', nstk1 
            write(3,*)' In Nextfiles_ONL; nindx1:      ', nindx1
            write(3,*)' In Nextfiles_ONL; nlist1:      ', nlist1
#endif

            IF (NLIST1 > 0 .AND.  NINDX1 > NLIST1) THEN
C              OVERUN INPUT LIST
               IRTFLG = -1
               RETURN
            ENDIF

C           OPEN NEXT INPUT FILE 
            NWANT1 = ILIST1(NINDX1)

         ENDIF


#if defined (SP_DBUGIO)
         write(3,*)' '
         write(3,*)' In Nextfiles-OLD; nwant1:   ', nwant1 
         write(3,*)' In Nextfiles-OLD; istk1:    ', istk1
         write(3,*)' In Nextfiles-OLD; Calling getoldimg  '
         write(3,*)' '
#endif

         CALL GETOLDIMG(LUN1,LUNXM1,FILPAT1,NWANT1,SAYIT, 
     &                  FOUROK,ISTK1,IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)'  '
         write(3,*)' In Nextfiles-OLD; After getoldimg  '
         write(3,*)' In Nextfiles-OLD; irtflg:   ', irtflg 
         write(3,*)' In Nextfiles-OLD; nwant1:   ', nwant1 
         write(3,*)' In Nextfiles-OLD; istk1:    ', istk1
         write(3,*)' In Nextfiles-OLD; nindx1:   ', nindx1
         write(3,*)' In Nextfiles-OLD; nindx2:   ', nindx2
         write(3,*)' In Nextfiles-OLD; is_bare1:          ', is_bare1
         write(3,*)'  '
#endif

         IF (IRTFLG < 0)    RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

         IF (IS_BARE1) THEN
C           INPUT FROM A BARE STACK 
            NINDX1 = ISTK1
            NINDX2 = ISTK1
         ENDIF
      ENDIF

  99  CONTINUE

      IF (LUN2 > 0) THEN  ! ----------------------------------------- 

C        OPEN NEXT NEW OUTPUT IMG, (IS A NEW FILE, IF NON-STACKED IMG)

         GOTAT2  = (INDEX(FILPAT2,'@') > 0)
         GOTAST2 = (INDEX(FILPAT2,'*') > 0)

         CALL GET_FILNAM_INFO(FILPAT2, LIST_VAL,
     &                        FIL_NOAT,FIL_DIRS,FIL_BASE,FIL_EXT,
     &                        IS_MRC, IS_MRCS, IS_BARE2, GOT_FILE2, 
     &                        ISTK2_BAD,   IRTFLG)



C        IS THIS IS A BARE STACK OPERATION?  (OK FOR SPIDER & MRC)
ccc         CALL LUNGETISBARE(LUN2,IS_BARE2,IRTFLG)

C        DETERMINE IF A SINGLE IMG, (NO LIST, NOT BARE STACK) 
         IS_ONE_IMG2 = (GOTAT2 .AND. (.NOT. GOTAST2) .AND.
     &                               (.NOT. IS_BARE2))                

#if defined (SP_DBUGIO)
         write(3,*) '  '
         write(3,*)' In Nextfiles-N0; is_bare2:            ', is_bare2
        !write(3,*)' In nextfiles-N0; ilist2(1): ', ilist2(1)
         write(3,*)' In nextfiles-N0; nindx2:    ', nindx2
         write(3,*)' In nextfiles-N0; nstk2:     ', nstk2
         write(3,*)' In nextfiles-N0: istk2:     ', istk2
         write(3,*)' In nextfiles-N0: gotat2:              ', gotat2
         write(3,*)' In nextfiles-N0: gotast2:             ', gotast2
         write(3,*)' In nextfiles-N0: is_one_img2:         ',is_one_img2
         write(3,*)' '
#endif

         IF (ISTK2 == -1 .AND. LUNXM2 > 0  ) THEN
C           XMIPP SELFILE LISTED IMAGE ----------------------

            NWANT2 = -1
#if defined (SP_DBUGIO)
            write(3,*)' In nextfiles-N1; nwant2:  ', nwant2
#endif


         ELSEIF (IS_ONE_IMG2 ) THEN

C           NON-STACK IMG/FILE (NO LIST) ------- FILE ALREADY OPEN

#if defined (SP_DBUGIO)
            write(3,*)' In Nextfiles-N2; filpat2:  ', filpat2(:nletp2)
            write(3,*)' In nextfiles-N2; Can skip opening'
            write(3,*)' '
#endif
            GOTO 999


         ELSEIF (LUN1 > 0 .AND. IS_BARE1 .AND. ISTK1 > 1 ) THEN
C           MRC OR SPIDER BARE STACK INPUT (NO LIST) -------
C           USE SAME OUTPUT IMG NUMBER AS INPUT

            NWANT2 = ISTK1                !!!???????????

#if defined (SP_DBUGIO)
            write(3,*)' In nextfiles-N3; istk2:   ', istk2
            write(3,*)' In nextfiles-N3; nwant2:  ', nwant2
#endif

         ELSEIF (IS_BARE2 ) THEN
C           MRC OR SPIDER BARE STACK OUTPUT (NO LIST) --------

            IF (ISTK2 < 0) ISTK2 = 1
            NWANT2 = ISTK2 + 1

            IF (LUN1 > 0 .AND. IS_BARE1) THEN
C              BARE STACK INPUT , USE SAME  OUTPUT IMAGE NUMBER
               NWANT2 = ISTK2 + 1
            ENDIF

#if defined (SP_DBUGIO)
            write(3,*)' In nextfiles-N4; istk2:     ', istk2
            write(3,*)' In nextfiles-N4; nwant2:    ', nwant2
#endif


         ELSEIF (NSTK2 == -2 .OR. NSTK2 == -1 .OR. NSTK2 > 0 ) THEN

C           NON-STACK IMG WITH/WITHOUT TEMPLATE LIST -------

            IF (NLIST2 > 0 .AND.  NINDX2 > NLIST2) THEN
C               OVERUN OUTPUT LIST
                IRTFLG = -99

                RETURN
            ENDIF

C           OPEN NEXT OUTPUT FILE 
            NWANT2 = ILIST2(NINDX2)

#if defined (SP_DBUGIO)
            write(3,*)' In nextfiles-N5; nstk2:   ', nstk2
            write(3,*)' In nextfiles-N5; nwant2:  ', nwant2
            write(3,*)' In nextfiles-N5; nindx2:  ', nindx2
            write(3,*)' In nextfiles-N5; nlist2:  ', nlist2
#endif

         ELSEIF (GOTAST2 .AND. .NOT. GOTAT2) THEN

C           NON-STACK IMG WITH TEMPLATE LIST --------------

            IF (NLIST2 > 0 .AND.  NINDX2 > NLIST2) THEN
C               OVERUN OUTPUT LIST
                IRTFLG = -99

                RETURN
            ENDIF

C           OPEN NEXT OUTPUT FILE 
            NWANT2 = ILIST2(NINDX2)

#if defined (SP_DBUGIO)
            write(3,*)' In nextfiles-N6; nindx2: ', nindx2
            write(3,*)' In nextfiles-N6; nwant2: ', nwant2
#endif
         ENDIF



#if defined (SP_DBUGIO)
         write(3,*)'   '
         write(3,*)' In nextfiles-99; is_bare2: ', is_bare2
         write(3,*)' In nextfiles-99; nstk2:    ', nstk2
         write(3,*)' In nextfiles-99; nindx2:   ', nindx2
         write(3,*)' In nextfiles-99; nlist2:   ', nlist2
         write(3,*)' In nextfiles-99; nwant2:   ', nwant2
         write(3,*)' In nextfiles-99; istk2:    ', istk2
         write(3,*)' In nextfiles-99; Calling getnewimg  '
#endif

         CALL GETNEWIMG(LUNCP,LUN2,LUNXM2,FILPAT2,NWANT2,
     &                  SAYIT,ISTK2,IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)' Gotnewimg; nwant2: ', nwant2
         write(3,*)' Gotnewimg; istk2:  ', istk2
#endif
       ENDIF

999   CONTINUE

#if defined (SP_DBUGIO)
      write(3,*)' Leaving nextfiles; filpat1:  ', filpat1(:nletp1)
      write(3,*)' Leaving nextfiles; filpat2:  ', filpat2(:nletp2)

      write(3,*)' Leaving nextfiles, nlist1,nlist2: ', nlist1,nlist2
      write(3,*)' Leaving nextfiles; nwant1,nwant2: ', nwant1,nwant2
      write(3,*)' Leaving nextfiles; nstk1,nstk2:   ', nstk1, nstk2
      write(3,*)' Leaving nextfiles; istk1,istk2:   ', istk1, istk2
      write(3,*)'     '
#endif

      END




C++********************************************************************C++   NEEDS WORK BELOW      !!!!!!!!!!!!!!!!!
C++
C++********************************************************************







C++*********************************************************************
C
C NEXTFILE.F    NEW                              12/15/06 ArDean Leith
C               OVERUN OUTPUT LIST = -99          1/15/12 ArDean Leith
C               MRC SUPPORT                       8/22/19 ArDean Leith
C
C **********************************************************************
C
C NEXTFILE(NINDX1, ILIST1, FOUROK, LUNXM1, NLIST1,   NSTK1,   
C          LUN1,   LUNCP,   FILPAT1,  DISP, ISTK1, IRTFLG)  
C
C PURPOSE:  GETS NEXT INPUT OR OUTPUT FILE FOR BARE STACK INPUT,
C           A NON STACKED IMAGE WITH/WITHOUT TEMPLATED LIST, 
C           OR A  STACKED IMAGE WITH/WITHOUT LIST
C           STACKS MUST BE OPENED WITH OPFILES!!!
C
C PARAMETERS: NINDX1         LIST INDEX                     (SENT/RET)
C             ILIST1         IMAGE NUMBER LIST                  (SENT)
C             FOUROK         FOURIER INPUT IS OK                (SENT)
C             LUNXM1         LUN FOR SELFILE INPUT              (SENT)
C             NLIST1         NUMBER OF IMAGES IN LIST           (SENT)
C             NSTK1          HIGHEST IMAGE IN STACK             (SENT)
C             LUN1           LUN FOR I/0                        (SENT)
C             LUNCP          LUN FOR OUTPUT HEADER COPY         (SENT)
C             FILPAT1        FILE NAME PATTERN                  (SENT)
C             DISP           IMAGE EXISTANCE                    (SENT)
C             ISTK1          IMAGE NUMBER                  (SENT/RET.)
C             IRTFLG         ERROR (0 IS OK, -1 IS END STACK)   (RET.)
C
C CALLED BY:  adds_n.f,    apscc.f,     bp3f.f,    bpcg.f,  bprp3.f, 
C             copytostk.f, copytovol.f, opfiles.f, qstat.f,
C             qstatloc.f,  ssnrb.f,     util2sup.f 
c
C CALL TREE:
C     ONLY USED AFTER OPFILES HAS BEEN CALLED:
C
C       NEXTFILE 
C          |      
C          |      Old image (see opfiles.f)
C          ` ---> GETOLDIMG    
C          |        |
C          |        |  Templated simple image:  IMG*** 
C          |        ` ---> FILGET --> OPFILEC  
C          |        | 
C          |        |  Templated stacked image: STK@*     
C          |        ` ---> LUN***
C          |        | 
C          |        |  Whole image stack:       STK@  
C          |        ` ---> LUNS***
C          |        | 
C          |        |  Templated stacked MRC:  *@MRC or *@.MRCS
C          |        ` ---> GETOLDIMG_MRC --> OPFILEC  
C          |      
C          |      
C          |      New image (see opfiles.f)                         
C          ` ---> GETNEWIMG
C                   |
C                   |  Templated simple image:  IMG*** 
C                   ` ---> FILGET --> OPFILEC  
C                   |             
C                   |  Templated stacked image: STK@*     
C                   ` ---> LUN***
C                   | 
C                   |  Whole image stack:       STK@  
C                   ` ---> LUNS***  
C                   | 
C                   |  Templated stacked MRC:  *@MRC or *@.MRCS
C                   ` ---> GETNEWIMG_MRC --> OPFILEC  
C
C23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
C--*********************************************************************  
      SUBROUTINE NEXTFILE(NINDX1,  ILIST1, 
     &                    FOUROK,  LUNXM1,
     &                    NLIST1,  NSTK1,   
     &                    LUN1,    LUNCP, 
     &                    FILPAT1, DISP,
     &                    ISTK1,   IRTFLG) 
 
      IMPLICIT NONE

      INTEGER           :: NINDX1 
      INTEGER           :: ILIST1(NLIST1) 
      LOGICAL           :: FOUROK
      INTEGER           :: LUNXM1 
      INTEGER           :: NLIST1 
      INTEGER           :: NSTK1, LUN1,LUNCP 
      CHARACTER(LEN=*)  :: FILPAT1
      CHARACTER(LEN=1)  :: DISP
      INTEGER           :: ISTK1
      INTEGER           :: IRTFLG

      INTEGER           :: NWANT1, NLETP1, it
      LOGICAL           :: SAYIT = .TRUE.
      LOGICAL           :: IS_BARE1,IS_MRC1
      LOGICAL           :: GOTAST1,GOTAT1

      INTEGER           :: lnblnkn     ! FUNCTION

      NINDX1 = NINDX1 + 1

      NLETP1 = lnblnkn(FILPAT1)

      GOTAT1  = (INDEX(FILPAT1,'@') > 0)
      GOTAST1 = (INDEX(FILPAT1,'*') > 0)

C     IS THIS IS A BARE STACK OPERATION?
      CALL LUNGETISBARE(LUN1,IS_BARE1,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN

C     IS THIS IS A MRC FILE?
      CALL LUNGETIS_MRC(LUN1,IS_MRC1,IRTFLG)
      IF (IRTFLG .NE. 0) RETURN
 
#if defined (SP_DBUGIO)
      write(3,*)'   '
      write(3,*)' In nextfile, ----------------------: '
      write(3,*)' In nextfile, filpat1:      ', filpat1(:nletp1)
      write(3,*)' In nextfile, gotast1:      ', gotast1
      write(3,*)' In nextfile, isbare1:      ', is_bare1
      write(3,*)' In nextfile, ismrc1:       ', is_mrc1
      write(3,*)' In nextfile, nindx1:       ', nindx1
      write(3,*)' In nextfile, nlist1:       ', nlist1 
      write(3,*)' In nextfile, nstk1:        ', nstk1
      write(3,*)' In nextfile, disp:         ', disp
      write(3,*)' In nextfile, istk1,irtflg: ', istk1,irtflg
      !write(3,*)' In nextfile, ilist1(1):   ', ilist1(1)
#endif

      IF ( IS_MRC1 ) THEN  ! DEC 2025
C        MRC FILE IN USE 

         CALL NEXTFILE_MRC(NINDX1,  ILIST1, 
     &                    FOUROK,  LUNXM1,
     &                    NLIST1,  NSTK1,   
     &                    LUN1,    LUNCP, 
     &                    FILPAT1, DISP,
     &                    ISTK1,   IRTFLG) 
         RETURN
      ENDIF

          
      IF (ISTK1 == -1  .AND. LUNXM1 > 0 ) THEN
C        XMIPP SELFILE LISTED IMAGE
         NWANT1 = -1

      ELSEIF ( IS_BARE1 ) THEN
C        SPIDER BARE STACK INPUT   (NO LIST)

         IF (ISTK1 < 0) ISTK1 = 1
         NWANT1 = ISTK1 + 1

         IF (NWANT1 > NSTK1) THEN
C           FINISHED THE WHOLE STACK
            IRTFLG = -1
            RETURN
         ENDIF

      ELSEIF (.NOT. GOTAT1.AND. GOTAST1 ) THEN
C        NON STACKED IMAGE WITH TEMPLATED LIST

         IF (NLIST1 > 0 .AND.  NINDX1 > NLIST1) THEN
C           OVERUN INPUT LIST
            IRTFLG = -1
            RETURN
         ENDIF

C        OPEN NEXT INPUT FILE 
         NWANT1 = ILIST1(NINDX1)

      ELSEIF (NSTK1 == -2  .OR.
     &        NSTK1 == -1  .OR.
     &        NSTK1  >    0 ) THEN
C        NON STACKED IMAGE WITHOUT TEMPLATED LIST
C        STACKED     IMAGE WITH/WITHOUT LIST

         IF (NLIST1 > 0 .AND. NINDX1 > NLIST1) THEN
C           OVERAN I/O LIST
            IRTFLG = -1
            RETURN

         ELSEIF ( NINDX1 >= NLIST1) THEN
C           OVERAN I/O LIST
            IRTFLG = -1
            RETURN
         ENDIF

C        OPEN NEXT I/O FILE 
         NWANT1 = ILIST1(NINDX1)

      ENDIF

      !write(3,*)' In nextfile, nindx1,nlist1,nstk1,istk1: ',
      !     &                   nindx1,nlist1,nstk1,istk1
      !write(3,*)'  In nextfile, calling get, nwant1,istk1,nstk1:',
      !&                                      nwant1,istk1,nstk1

      IF (DISP == 'O' .OR. DISP == 'B' .OR. 
     &    DISP == 'Z' .OR. DISP == 'E') THEN 
  
C        OPEN NEXT INPUT FILE ------------------------------

#if defined (SP_DBUGIO_NEVER)
         write(3,*) ' Call gotoldimg, nwant1: ',nwant1,istk1
#endif

         CALL GETOLDIMG(LUN1,LUNXM1,FILPAT1,NWANT1,SAYIT, 
     &                  FOUROK,ISTK1,IRTFLG)

#if defined (SP_DBUGIO)
         write(3,*)' In nextfile; nwant1:   ', nwant1
         write(3,*)' In nextfile; istk1:    ', istk1
         write(3,*)' In nextfile; irtflg:   ', irtflg
         write(3,*)'     '
#endif

         IF (IRTFLG    < 0) RETURN    ! END OF WHOLE-STACK
         IF (IRTFLG .NE. 0) RETURN    ! ERROR

      ELSE   


C        OPEN NEXT OUTPUT FILE ------------------------------ 
#if defined (SP_DBUGIO_NEVER)
         !write(3,*) ' Call gotnewimg, nwant1: ',nwant1,istk1
#endif

         CALL GETNEWIMG(LUNCP,LUN1,LUNXM1,FILPAT1,NWANT1,
     &                  SAYIT,ISTK1,IRTFLG)

#if defined (SP_DBUGIO_NEVER)
          write(3,*) ' Gotnewimg, nwant1,istk1,nstk1,irtflg:',
     &                            nwant1,istk1,nstk1,irtflg
#endif
         IF (IRTFLG .NE. 0) RETURN   ! ERROR

      ENDIF

      IF (IS_BARE1) THEN
C        BARE STACK, OUTPUT IMAGE HAS SAME # AS INPUT ALWAYS
         NINDX1 = ISTK1
      ENDIF

#if defined (SP_DBUGIO)
      write(3,*)' Leaving nextfile, is_bare1: ', is_bare1
      write(3,*)' Leaving nextfile, nlist1:   ', nlist1
      write(3,*)' Leaving nextfile; nwant1:   ', nwant1
      write(3,*)' Leaving nextfile; nstk1:    ', nstk1
      write(3,*)' Leaving nextfile; istk1:    ', istk1
      write(3,*)'     '
#endif

      END













