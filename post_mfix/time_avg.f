!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TIME_AVG (TAVG,NTAVG,NPOINTS,VAR_INDEX,M_USE)                C
!  Purpose: Update the TIME-AVERAGED sum                               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 03-MAR-92  C
!  Reviewer:                                                           C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: EP_g, P_g, P_star, U_g, V_g, W_g, U_s, V_s    C
!                        W_s, ROP_s, T_g, T_s1, T_s2                   C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: L                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE TIME_AVG(TAVG,NTAVG,NPOINTS,VAR_INDEX,M_USE)
!
!
      Use param
      Use param1
      Use fldvar
      Use geometry
      IMPLICIT NONE
!
!     passed arguments
!
      DOUBLE PRECISION TAVG(DIMENSION_3,*)
      INTEGER          NPOINTS
      INTEGER          VAR_INDEX
      INTEGER          NTAVG
      INTEGER          M_USE
!
!     local variables
!
      INTEGER L
!
      DO L = 1,IJKMAX2
        IF (VAR_INDEX.EQ.01) TAVG(L,NTAVG) = TAVG(L,NTAVG) + EP_g(L)
        IF (VAR_INDEX.EQ.02) TAVG(L,NTAVG) = TAVG(L,NTAVG) + P_g(L)
        IF (VAR_INDEX.EQ.03) TAVG(L,NTAVG) = TAVG(L,NTAVG) + P_star(L)
        IF (VAR_INDEX.EQ.04) TAVG(L,NTAVG) = TAVG(L,NTAVG) + U_g(L)
        IF (VAR_INDEX.EQ.05) TAVG(L,NTAVG) = TAVG(L,NTAVG) + V_g(L)
        IF (VAR_INDEX.EQ.06) TAVG(L,NTAVG) = TAVG(L,NTAVG) + W_g(L)
        IF (VAR_INDEX.EQ.07) TAVG(L,NTAVG) = TAVG(L,NTAVG) +&
                                                     U_s(L,M_USE)
        IF (VAR_INDEX.EQ.08) TAVG(L,NTAVG) = TAVG(L,NTAVG) +&
                                                     V_s(L,M_USE)
        IF (VAR_INDEX.EQ.09) TAVG(L,NTAVG) = TAVG(L,NTAVG) +&
                                                     W_s(L,M_USE)
        IF (VAR_INDEX.EQ.10) TAVG(L,NTAVG) = TAVG(L,NTAVG) +&
                                                     ROP_s(L,M_USE)
        IF (VAR_INDEX.EQ.11) TAVG(L,NTAVG) = TAVG(L,NTAVG) + T_g(L)
        IF (VAR_INDEX.EQ.12) TAVG(L,NTAVG) = TAVG(L,NTAVG)& 
                                             + T_s(L,M_USE)
        IF (VAR_INDEX.EQ.13) TAVG(L,NTAVG) = TAVG(L,NTAVG) &
                                             + T_s(L,M_USE)
      END DO
      NPOINTS = NPOINTS + 1
!
      RETURN
      END
