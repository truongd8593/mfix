C              record length used in open statement for unformatted,
C              direct access file, with 512 bytes per record
      INTEGER  OPEN_N1
C
C              number of DOUBLE PRECISION words in 512 bytes
      INTEGER  NWORDS_DP
C
C              number of REAL words in 512 bytes
      INTEGER  NWORDS_R
C
C              number of INTEGER words in 512 bytes
      INTEGER  NWORDS_I
C
      COMMON /MACHINE/ OPEN_N1, NWORDS_DP, NWORDS_R, NWORDS_I
