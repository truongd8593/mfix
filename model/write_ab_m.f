!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: Write_Ab_m(A_m, b_m, IJKMAX2, M, IER)                  C                     C
!  Purpose: Write the sparse matrix coefficients and the               C
!           source vector.                                             C
!                                                                      C
!                                                                      C
!  Author: M. Syamlal                                 Date: 16-MAY-96  C
!  Reviewer:                                          Date:            C
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
      SUBROUTINE WRITE_AB_M(A_M, B_M, IJKMAX2A, M, IER)    ! pnicol
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
      USE matrix 
!//d      USE funits 
      USE compar        !//d
      USE mpi_utility   !//d
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
!                      Local index
      INTEGER          L
!
!                      Phase index
      INTEGER          M
!
!                      cell index
      INTEGER          IJK
!
!                      Maximum dimension
      INTEGER          IJKMAX2A  ! pnicol
!
!                      Septadiagonal matrix A_m
      DOUBLE PRECISION A_m(DIMENSION_3, -3:3, 0:DIMENSION_M)
!
!                      Source vector
      DOUBLE PRECISION b_m(DIMENSION_3, 0:DIMENSION_M)

      double precision, allocatable :: array1(:) , array2(:)   !//
      double precision              :: am(-3:3)                !//
!
!-----------------------------------------------
!
      if (myPE.ne.PE_IO)  return    !//d only put in ROOT log file

      allocate (array1(ijkmax2))    !//d
      allocate (array2(ijkmax2))    !//d

      CALL START_LOG 
      WRITE (UNIT_LOG,*) ' A_m and B_m arrays below are in the '
      WRITE (UNIT_LOG,*) ' mfix INTERNAL order'
      WRITE (UNIT_LOG,*) ' '
      WRITE (UNIT_LOG, '(A,A)') &
         '  IJK  b         s         w         p         e       ', &
         '  n         t         Source' 


      call gather(b_m(:,M),array2,root)        !//d
      DO IJK = 1, IJKMAX2 

!// start
         do L = -3,3
            call gather(a_m(:,L,M),array1,root)
            am(l) = a_m(ijk,L,M)
         end do
!//         WRITE (UNIT_LOG, '(I5, 8(1X,G9.2))') IJK, (A_M(IJK,L,M),L=-3,3), B_M(&
!//            IJK,M) 
         WRITE (UNIT_LOG, '(I5, 8(1X,G9.2))') IJK, (AM(L),L=-3,3), array2(IJK) 

!// end
      END DO 
      CALL END_LOG 

      deallocate (array1)    !//
      deallocate (array2)    !//

      RETURN  
      END SUBROUTINE WRITE_AB_M 
