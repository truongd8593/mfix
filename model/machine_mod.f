      MODULE machine


!      Use param
!      Use param1


!              record length used in open statement for unformatted,
!              direct access file, with 512 bytes per record
      INTEGER  OPEN_N1
!
!              number of DOUBLE PRECISION words in 512 bytes
      INTEGER  NWORDS_DP
!
!              number of REAL words in 512 bytes
      INTEGER  NWORDS_R
!
!              number of INTEGER words in 512 bytes
      INTEGER  NWORDS_I
!
      LOGICAL :: JUST_FLUSH = .TRUE.


      END MODULE machine                                                                         
