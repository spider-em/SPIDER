
C      F90 EXPLICIT INTERFACES NEEDED BY SUBROUTINES THAT CALL
C      GETDOCDAT

C       --------------------------------------

        INTERFACE ! ---------------------
        SUBROUTINE GETDOCDAT(PROMPT,ASKNAME,DOCNAM,LUNDOC,GETSIZE,
     &                       MAXXT,MAXYT,RPOINTER,IRTFLG)

        CHARACTER *(*), INTENT(IN) ::                 PROMPT
        LOGICAL, INTENT(IN) ::                        ASKNAME
        CHARACTER *(*), INTENT(INOUT) ::              DOCNAM
        INTEGER, INTENT(IN) ::                        LUNDOC
        LOGICAL, INTENT(IN) ::                        GETSIZE
        INTEGER, INTENT(INOUT) ::                     MAXXT
        INTEGER, INTENT(INOUT) ::                     MAXYT
        REAL, DIMENSION(:,:), POINTER  ::             RPOINTER
        INTEGER, INTENT(OUT) ::                       IRTFLG
        END SUBROUTINE GETDOCDAT
        END INTERFACE !--------------------


C      ---------------------------------------------

