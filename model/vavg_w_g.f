!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: VAVG_W_g                                               C
!  Purpose: Calculate volume averaged W_g                              C
!                                                                      C
!  Author: M. Syamlal                                 Date: 28-APR-94  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      DOUBLE PRECISION FUNCTION VAVG_W_G () 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run 
      USE parallel 
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      USE compar       
      USE mpi_utility  
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
! 
!                      Indices 
      INTEGER          IJK, I_DUM 
! 
!                      Integral of W_g*EP_g for entire volume 
      DOUBLE PRECISION SUM_W_g 
! 
!                      Total volume of computational cells 
      DOUBLE PRECISION SUM_VOL
!
!                      To check when velocity becomes NaN to auto_restart MFIX
      CHARACTER *80 notnumber 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Integrate the velocity values for the whole domain,
!
      SUM_W_G = ZERO 
      SUM_VOL = ZERO 

!!$omp   parallel do private(IJK) reduction(+:SUM_VOL,SUM_W_G)
      DO IJK = IJKSTART3, IJKEND3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK), J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN 
            SUM_VOL = SUM_VOL + VOL_W(IJK) 
            SUM_W_G = SUM_W_G + W_G(IJK)*EP_G(IJK)*VOL_W(IJK) 
         ENDIF 
      END DO 

      CALL GLOBAL_ALL_SUM(SUM_VOL)
      CALL GLOBAL_ALL_SUM(SUM_W_G)
      VAVG_W_G = SUM_W_G/SUM_VOL
      WRITE(notnumber,*) VAVG_W_G
! Check for NaN's
! See if velocity (a real number) contains a letter "n" or symbol "?"
! in which case it's a NaN (Not a Number)
!
      IF(INDEX(notnumber,'?') > 0 .OR.     &
         INDEX(notnumber,'n') > 0 .OR.     &
         INDEX(notnumber,'N') > 0 ) THEN
        write(*,*) VAVG_W_G,  ' NaN being caught in VAVG_W_G.f '
        AUTOMATIC_RESTART = .TRUE.
      ENDIF 
!
      RETURN  
      END FUNCTION VAVG_W_G 



      DOUBLE PRECISION FUNCTION VAVG_Flux_W_G () 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE fldvar
      USE bc
      USE geometry
      USE physprop
      USE indices
      USE compar       
      USE mpi_utility  
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
! 
!                      Indices 
      INTEGER          IJK 
! 
!                      Integral of W_g*EP_g for entire volume 
      DOUBLE PRECISION SUM_W_g 
! 
!                      Total volume of computational cells 
      DOUBLE PRECISION SUM_AREA 
! 
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!  Integrate the velocity values for the whole domain,
!
      SUM_W_G = ZERO 
      SUM_AREA = ZERO 

!!$omp   parallel do private(IJK) reduction(+:SUM_AREA,SUM_W_G)
      DO IJK = IJKSTART3, IJKEND3
      IF(.NOT.IS_ON_myPE_wobnd(I_OF(IJK), J_OF(IJK), K_OF(IJK))) CYCLE
         IF (FLUID_AT(IJK)) THEN 
	   IF(W_g(IJK) > ZERO) THEN
             SUM_W_G = SUM_W_G + W_G(IJK)*ROP_G(IJK)*AXY(IJK) 
	   ELSE
             SUM_W_G = SUM_W_G + W_G(IJK)*ROP_G(TOP_OF(IJK))*AXY(IJK) 
	   ENDIF
           SUM_AREA = SUM_AREA + AXY(IJK)
         ENDIF 
      END DO 

      CALL GLOBAL_ALL_SUM(SUM_AREA)
      CALL GLOBAL_ALL_SUM(SUM_W_G)
      VAVG_Flux_W_G = SUM_W_G/SUM_AREA 
!
      RETURN  
      END FUNCTION VAVG_Flux_W_G 
