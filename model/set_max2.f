!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_MAX2                                               C
!  Purpose: calculate IMAX1,IMAX2,JMAX1,JMAX2,KMAX1,KMAX2,IJMAX2       C
!                     IJKMAX2                                          C
!                                                                      C
!  Author: P. Nicoletti                               Date: 04-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX, JMAX, KMAX, NO_I, NO_J, NO_K            C
!  Variables modified: IMAX1, IMAX2, JMAX1, JMAX2, KMAX1, KMAX2        C
!                      IJMAX2, IJKMAX2, IMIN1, JMIN1, KMIN1            C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_MAX2 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      IF (NO_I) THEN 
         IMIN1 = 1 
         IMAX1 = 1 
         IMAX2 = 1 
      ELSE 
         IMIN1 = 2 
         IMAX1 = IMAX + 1 
         IMAX2 = IMAX + 2 
      ENDIF 
!
      IF (NO_J) THEN 
         JMIN1 = 1 
         JMAX1 = 1 
         JMAX2 = 1 
      ELSE 
         JMIN1 = 2 
         JMAX1 = JMAX + 1 
         JMAX2 = JMAX + 2 
      ENDIF 
!
      IF (NO_K) THEN 
         KMIN1 = 1 
         KMAX1 = 1 
         KMAX2 = 1 
      ELSE 
         KMIN1 = 2 
         KMAX1 = KMAX + 1 
         KMAX2 = KMAX + 2 
      ENDIF 
!
      IJMAX2 = IMAX2*JMAX2 
      IJKMAX2 = IMAX2*JMAX2*KMAX2 
      IF (DO_K) THEN 
         IJKMIN1 = IJMAX2 + 1 
         IJKMAX1 = IJKMAX2 - IJMAX2 
      ELSE 
         IJKMIN1 = IMAX2 + 1 
         IJKMAX1 = IJKMAX2 - IMAX2 
      ENDIF 
!
      RETURN  
      END SUBROUTINE SET_MAX2 
