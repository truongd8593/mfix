!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_GEOMETRY1                                          C
!  Purpose: Calculate cell volumes and face areas                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 1-MAY-96   C
!  Reviewer:                                                           C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: COORDINATES, IMAX2, DT, DX, JMAX2, DY, KMAX2, C
!                        DZ,                                           C
!                                                                      C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE SET_GEOMETRY1 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!  Include param.inc file to specify parameter values
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE parallel 
      USE run
      USE geometry
      USE indices
      USE compar
      USE funits !//AIKEPARDBG
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
      INTEGER          I, J, K, IP, JP, KP, IJK 
! 
!-----------------------------------------------
!
!
!                      Indices
!
!//AIKEPARDBG
!       do ijk=ijkstart3,ijkend3
!	 write(UNIT_LOG,"(' I_OF(',I4,') = ',I5, &
!	              & '  J= ',I5,'  K= ',I5)") &
!		      & ,IJK,I_OF(ijk),J_OF(ijk),K_OF(ijk)  !//AIKEPARDBG
!	 write(UNIT_LOG,"(' IP1(',I4,') = ',I5, &
!	              & '  JP1= ',I5,'  KP1= ',I5)") I_OF(ijk),IP1(I_OF(ijk)), &
!		      &  JP1(J_OF(IJK)),KP1(K_OF(IJK))  !//AIKEPARDBG
!       end do
!      call mfix_exit(myPE) !//AIKEPARDBG       

!!$omp  parallel do private( I, J, K, IP, JP, KP, IJK)  &
!!$omp  schedule(dynamic,chunk_size)
!
!//? Make sure all references to KP are used in 1D arrays otherwise need to
!//? communicate the calculated values also check prequisites.


!// 350 1025 change do loop limits: 1,ijkmax2-> ijkstart3, ijkend3
      DO IJK = ijkstart3, ijkend3
      
!
         I = I_OF(IJK) 
         IP = IP1(I) 
         J = J_OF(IJK) 
         JP = JP1(J) 
         K = K_OF(IJK) 
         KP = KP1(K) 

!// 200 1108 Check if any of i,j,k is point the 2nd ghost layer
!//          KP is k+1 which is the problematic due data dependency,
!//          however, only used in DZ() which is known globally on all PEs	 
	 
!
         VOL(IJK) = DX(I)*DY(J)*(X(I)*DZ(K)) 
         VOL_U(IJK) = HALF*(DX(I)+DX(IP))*DY(J)*(HALF*(X(I)+X(IP))*DZ(K)) 
         VOL_V(IJK) = DX(I)*HALF*(DY(J)+DY(JP))*(X(I)*DZ(K)) 
         VOL_W(IJK) = DX(I)*DY(J)*(X(I)*HALF*(DZ(K)+DZ(KP))) 
!
         AYZ(IJK) = DY(J)*(X_E(I)*DZ(K)) 
         AYZ_U(IJK) = DY(J)*(X(IP)*DZ(K)) 
         AYZ_V(IJK) = HALF*(DY(J)+DY(JP))*(X_E(I)*DZ(K)) 
         AYZ_W(IJK) = DY(J)*(X_E(I)*HALF*(DZ(K)+DZ(KP))) 
!
         AXY(IJK) = DX(I)*DY(J) 
         AXY_U(IJK) = HALF*(DX(I)+DX(IP))*DY(J) 
         AXY_V(IJK) = DX(I)*HALF*(DY(J)+DY(JP)) 
         AXY_W(IJK) = AXY(IJK) 
!
         AXZ(IJK) = DX(I)*(X(I)*DZ(K)) 
         AXZ_U(IJK) = HALF*(DX(I)+DX(IP))*(HALF*(X(I)+X(IP))*DZ(K)) 
         AXZ_V(IJK) = AXZ(IJK) 
         AXZ_W(IJK) = DX(I)*(X(I)*HALF*(DZ(K)+DZ(KP))) 
	 

      END DO 
      RETURN  
      END SUBROUTINE SET_GEOMETRY1 
