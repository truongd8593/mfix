      MODULE parallel


      Use param
      Use param1


!   Definitions for parallel programming
!
!                      Chunks of arrays used by each processor
      INTEGER          Chunk_size
      PARAMETER (chunk_size=8)

      LOGICAL          USE_DOLOOP
      PARAMETER( USE_DOLOOP=.TRUE. )


      END MODULE parallel                                                                        
