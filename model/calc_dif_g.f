!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_DIF_g(IER)
!  Purpose: Calculate the effective diffusivity of fluid phase        C
!                                                                      C
!  Author:M. Syamlal                                  Date: 13-FEB-98  C
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
      SUBROUTINE CALC_DIF_G(IER) 
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
      USE compar   !//d
      USE sendrecv  !// 400
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
      INTEGER          IJK, N
!
!-----------------------------------------------
      INCLUDE 'function.inc'
!
!
      IF (DIF_G0 /= UNDEFINED) RETURN  
!!$omp  parallel do private(ijk) &
!!$omp& schedule(dynamic,chunk_size)

      DO N = 1, NMAX(0) 
!// 350 1112 MTP changed do loop limits 1,ijkmax2 ==> ijkstart3, ijkend3
         DO IJK = IJKSTART3, IJKEND3 	 
            IF (.NOT.WALL_AT(IJK)) THEN 
!           Gas diffusion coefficient
!           Bird, Stewart, and Lightfoot (1960) -- CO2--N2 at 298.2 K
               DIF_G(IJK,N) = ROP_G(IJK)*0.165 
            ELSE 
               DIF_G(IJK,N) = ZERO 
            ENDIF 
         END DO 
      END DO 
      
!//S 1113 try to move this COMM to the end of transport_prop to do all COMMs
!//       at certain locations, provided that no data dependency in between.

!// 400 1113 MTP communicate boundaries
      CALL SEND_RECV(DIF_G, idbg)     
   
      RETURN  
      END SUBROUTINE CALC_DIF_G 
