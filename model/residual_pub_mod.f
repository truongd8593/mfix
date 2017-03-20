! -*- f90 -*-
MODULE residual_pub

      Use param, only: dim_n, dim_m
      use residual


!                      Residual sum within a group of equations
      LOGICAL          :: GROUP_RESID
   CONTAINS

      FUNCTION GET_RESID_STRING(INDEX)
         IMPLICIT NONE
         CHARACTER(LEN=4) :: GET_RESID_STRING
         INTEGER, INTENT(IN) :: INDEX

         GET_RESID_STRING = RESID_STRING(INDEX)

      END FUNCTION GET_RESID_STRING

      FUNCTION GET_RESID_GRP_STRING(INDEX)
         IMPLICIT NONE
         CHARACTER(LEN=8) :: GET_RESID_GRP_STRING
         INTEGER, INTENT(IN) :: INDEX

         GET_RESID_GRP_STRING = RESID_GRP_STRING(INDEX)

      END FUNCTION GET_RESID_GRP_STRING

      FUNCTION GET_RESID(INDEX)
         IMPLICIT NONE
         DOUBLE PRECISION :: GET_RESID
         INTEGER, INTENT(IN) :: INDEX
         INTEGER :: RI, RI2

         IF (INDEX > SIZE(RESID_INDEX,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR INDEX ",INDEX
            GET_RESID = 0.0
            RETURN
         ENDIF
         RI = RESID_INDEX(INDEX,1)
         RI2 = RESID_INDEX(INDEX,2)
         IF (RI > SIZE(RESID,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR RESID_INDEX 1 ",RI
            GET_RESID = 0.0
            RETURN
         ENDIF
         IF (RI2 > SIZE(RESID,1)) THEN
            ! PRINT *,__FILE__," INVALID VALUE FOR RESID_INDEX 2 ",RI2
            GET_RESID = 0.0
            RETURN
         ENDIF
         GET_RESID = RESID(RI,RI2)

      END FUNCTION GET_RESID

      FUNCTION GET_RESID_GRP(INDEX)
         IMPLICIT NONE
         DOUBLE PRECISION :: GET_RESID_GRP
         INTEGER, INTENT(IN) :: INDEX

          GET_RESID_GRP = RESID_GRP(INDEX)

      END FUNCTION GET_RESID_GRP

   END MODULE residual_pub