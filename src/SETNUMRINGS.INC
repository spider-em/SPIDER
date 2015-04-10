C       ----------- SETNUMRINGS ------------------------------------

C      SETNUMRINGS EXPLICIT INTERFACE NEEDED BY SUBROUTINES THAT CALL

        INTERFACE  

        SUBROUTINE SETNUMRINGS(MR,NR,ISKIP,MODE, IRAY,
     &                         NUMR,NRING,LCIRC,
     &                         IRTFLG)

        INTEGER, INTENT(IN)     :: MR,NR,ISKIP ! RING START, END, SKIP
        INTEGER, INTENT(IN)     :: IRAY        ! RAY SKIP
        CHARACTER*1, INTENT(IN) :: MODE        ! FULL OR HALF CIRCLE
        INTEGER, ALLOCATABLE    :: NUMR(:,:)   ! RING SPEC. ARRAY
        INTEGER, INTENT(OUT)    :: NRING       ! # OF RINGS
        INTEGER, INTENT(OUT)    :: LCIRC       ! TOTAL LENGTH OF CIRCLE ARRAYS
        INTEGER, INTENT(OUT)    :: IRTFLG      ! ERROR FLAG

        END SUBROUTINE SETNUMRINGS

        END INTERFACE 

