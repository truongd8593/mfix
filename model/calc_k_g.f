!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_K_g(IER)                                          C
!  Purpose: Calculate the effective conductivity of fluid phase        C
!                                                                      C
!  Author:M. Syamlal                                  Date: 24-APR-96  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
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
      SUBROUTINE CALC_K_G(IER) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!     Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE physprop
      USE fldvar
      USE geometry
      USE indices
      USE constant
      USE compar
      USE sendrecv
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      Error index
      INTEGER          IER
!
!                      Indices
      INTEGER          IJK
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      IF (K_G0 /= UNDEFINED) RETURN  

!!$omp parallel do private(ijk) &
!!$omp& schedule(dynamic,chunk_size)
      DO IJK = ijkstart3, ijkend3      
         IF (FLUID_AT(IJK)) THEN 
!           Gas conductivity (air) in cal/s.cm.K
!           Bird, Stewart, and Lightfoot (1960) -- Temperature dependence from formula
!           8.3-12 on p. 255 and conductivity value at 300 K from p. 263
            K_G(IJK) = 6.02D-5*SQRT(T_G(IJK)/300.) 
         ELSE 
            K_G(IJK) = ZERO 
         ENDIF 
      END DO 

      CALL send_recv(K_G, 2)     
      
      RETURN  
      END SUBROUTINE CALC_K_G 

!// Comments on the modifications for DMP version implementation      
!// 001 Include header file and common declarations for parallelization 
!// 350 Changed do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
!// 400 Added sendrecv module and send_recv calls for COMMunication
