      INTERFACE      ! -------------------------------------- INTERFACE
      ! C FUNCTION TO INITIALIZE 1D FFT  
      INTEGER FUNCTION cufftPlan1d(IFFTPLAN, NVAL, ITYPE, NBATCH)
     &           BIND(C,NAME='cufftPlan1d')
        USE ISO_C_BINDING
        INTEGER               :: IFFTPLAN   ! Pointer to cufftHandle object
        INTEGER(C_INT), VALUE :: NVAL       ! Transform size 
        INTEGER(C_INT), VALUE :: ITYPE      ! Transform data type
        INTEGER(C_INT), VALUE :: NBATCH     ! # of transforms of size = nx
      END FUNCTION cufftPlan1d
 

      ! C FUNCTION TO EXECUTE 1D FFT
      INTEGER FUNCTION cufftExecC2R(IFFTPLAN, FIN, FOUT) 
     &           BIND(C,NAME='cufftExecC2R')
        USE ISO_C_BINDING
        INTEGER(C_INT), VALUE    :: IFFTPLAN  ! Pointer to cufftHandle object
        COMPLEX(C_FLOAT), DEVICE :: FIN       ! Complex input data
        COMPLEX(C_FLOAT), DEVICE :: FOUT      ! Real output data
      END FUNCTION cufftExecC2R
  
                               
      ! C FUNCTION TO DESTROY 1D FFT PLAN
      INTEGER FUNCTION cufftDestroy(IFFTPLAN)
     &           BIND(C,NAME='cufftDestroy')
         USE ISO_C_BINDING
         INTEGER(C_INT) :: IFFTPLAN          
      END FUNCTION cufftDestroy


      ! CUBLAS FUNCTION TO FIND INDEX OF MAX IN BUFF
      INTEGER FUNCTION cublasIsamax(N, CBUF, INCX)
     &           BIND(C,NAME='cublasIsamax')
         USE ISO_C_BINDING
         INTEGER(C_INT),  VALUE   :: N
         COMPLEX(C_FLOAT),DEVICE  :: CBUF
         INTEGER(C_INT),  VALUE   :: INCX
      END FUNCTION cublasIsamax

      END INTERFACE                        ! ------------ END INTERFACE
