!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OUT_ARRAY_KC (ARRAY, K)                                C
!  Purpose: print out a 2D (constant k-plane) array to standard output C
!           (character)                                                C
!                                                                      C
!  Author: P.Nicoletti                                Date: 02-DEC-91  C
!  Reviewer: W. Rogers, M. Syamlal, S. Venkatesan     Date: 31-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: IMAX2, JMAX2                                  C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: NCOL, NTAB, LL1, LL2, LL3, IFORM1, IFORM2, IJK, IJ2
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE OUT_ARRAY_KC(ARRAY, K_PLANE) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE fldvar
      USE physprop
      USE indices
      USE funits 
      USE compar        !//d
      USE mpi_utility   !//EFD to use gather
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!                      2D array to print out
      CHARACTER*3      ARRAY(*)
!
!                      K plane
!//EFD use K_plane instead of K
      INTEGER,intent(in) ::          K_PLANE

!
! local variables
!
!
!                      A line of characters to print
      CHARACTER*132    LINE
!
!                      number of columns to print out across the page
      INTEGER          NCOL
!
!                      number of tables the 2D array must be split into
!                      for printing
      INTEGER          NTAB
!
!                      loop indices
      INTEGER          LL1, LL2, LL3, LL4
!
!                      start and end 'I' for current table
      INTEGER          IFORM1 , IFORM2
!
!                      start 'IJ' and end 'IJ' for a given 'J' to print out
      INTEGER          IJK , IJ2
!    
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!//EFD create global versions of flag arrays
      INTEGER, DIMENSION(:), allocatable :: FLAG_E_G, FLAG_N_G, FLAG_T_G, FLAG_TEMP
      integer :: i,j,k, ijk_io, ijk_gl, flag_size

!//SP Global Functions
!                      Logical function to identify IP at East of the cell
      LOGICAL          IP_AT_E_G
!
!                      Logical function to identify IP at North of the cell
      LOGICAL          IP_AT_N_G
!
!                      Logical function to identify IP at Top of the cell
      LOGICAL          IP_AT_T_G
!
!                      Logical function to identify IP at East of the cell
      LOGICAL          SIP_AT_E_G
!
!                      Logical function to identify IP at North of the cell
      LOGICAL          SIP_AT_N_G
!
!                      Logical function to identify IP at Top of the cell
      LOGICAL          SIP_AT_T_G
!


      IP_AT_E_G(IJK)     = FLAG_E_G(IJK) .LT. 1000
      IP_AT_N_G(IJK)     = FLAG_N_G(IJK) .LT. 1000
      IP_AT_T_G(IJK)     = FLAG_T_G(IJK) .LT. 1000
      SIP_AT_E_G(IJK)    = (FLAG_E_G(IJK) .LT. 2000)
      SIP_AT_N_G(IJK)    = (FLAG_N_G(IJK) .LT. 2000)
      SIP_AT_T_G(IJK)    = (FLAG_T_G(IJK) .LT. 2000)

      INCLUDE 'function.inc'

!//EFD allocate storage for temporary flag arrays

      flag_size = 1
      if (myPE.eq.root) then
          flag_size = ijkmax3
      endif

      allocate( flag_e_g(flag_size) )
      allocate( flag_n_g(flag_size) )
      allocate( flag_t_g(flag_size) )
      allocate( flag_temp(flag_size) )

      flag_e_g(:) = undefined_i
      flag_n_g(:) = undefined_i
      flag_t_g(:) = undefined_i
      flag_temp(:) = undefined_i

      

      call gather( flag_e, flag_e_g )
      call gather( flag_n, flag_n_g )
      call gather( flag_t, flag_t_g )


!//EFD
!//     reorder to conform with funijk_io
  if (myPE.eq.PE_IO) then
 
      flag_temp(:) = flag_e_g(:)

      do k=kmin2,kmax2
      do j=jmin2,jmax2
      do i=imin2,imax2
          ijk_io = funijk_io(i,j,k)
          ijk_gl = funijk_gl(i,j,k) 
          
          flag_e_g(ijk_io) = flag_temp( ijk_gl )
      enddo
      enddo
      enddo


	  
      flag_temp(:) = flag_n_g(:)

      do k=kmin2,kmax2
      do j=jmin2,jmax2
      do i=imin2,imax2
          ijk_io = funijk_io(i,j,k)
          ijk_gl = funijk_gl(i,j,k)
          
          flag_n_g(ijk_io) = flag_temp( ijk_gl )
      enddo
      enddo
      enddo


      flag_temp(:) = flag_t_g(:)

      do k=kmin2,kmax2
      do j=jmin2,jmax2
      do i=imin2,imax2
          ijk_io = funijk_io(i,j,k)
          ijk_gl = funijk_gl(i,j,k)
          
          flag_t_g(ijk_io) = flag_temp( ijk_gl )
      enddo
      enddo
      enddo

      
      
      


!
! NOTE:  IF NCOL IS CHANGED TO A NUMBER GREATER THAN 30, THEN THE "30"
!        IN FORMATS 5050 AND 5100 MUST BE CHANGED TO THAT NUMBER.
!
      NCOL = 30 
      NTAB = IMAX2/NCOL + 1 
      IF (MOD(IMAX2,NCOL) == 0) NTAB = NTAB - 1 
!
      DO LL1 = 1, NTAB 
         IFORM1 = 1 + NCOL*(LL1 - 1) 
         IFORM2 = NCOL*LL1 
         IFORM2 = MIN(IFORM2,IMAX2) 
         WRITE (UNIT_OUT, 5050) (LL3,LL3=IFORM1,IFORM2) 
         DO LL2 = JMAX2, 1, -1 
            IJK = funijk_io(IFORM1,LL2,1) 
            IJ2 = funijk_io(IFORM2,LL2,1) 
!efd
            WRITE (LINE, 5100) LL2, (ARRAY(LL3),LL3=IJK,IJ2) 
!            WRITE (LINE, 5100) LL2, (ARRAY(funijk_io(LL3,LL2,1)),LL3=IFORM1,IFORM2)

            IJK = funijk_io(IFORM1,LL2,K) 
            IJ2 = funijk_io(IFORM2,LL2,K) 
            LL4 = 12 
            DO LL3 = IJK, IJ2 
               LL4 = LL4 + 4 
!
!//SP Global Functions
               IF (IP_AT_E_G(LL3)) THEN 
                  LINE(LL4:LL4) = 'E' 
               ELSE IF (SIP_AT_E_G(LL3)) THEN 
                  LINE(LL4:LL4) = 'e' 
               ENDIF 
!
!//SP Global Functions
               IF (IP_AT_N_G(LL3)) THEN 
                  LINE(LL4:LL4) = 'N' 
               ELSE IF (SIP_AT_N_G(LL3)) THEN 
                  LINE(LL4:LL4) = 'n' 
               ENDIF 
!
!//SP Global Functions
               IF (IP_AT_T_G(LL3)) THEN 
                  LINE(LL4:LL4) = 'T' 
               ELSE IF (SIP_AT_T_G(LL3)) THEN 
                  LINE(LL4:LL4) = 't' 
               ENDIF 
            END DO 
           WRITE (UNIT_OUT, '(A)') LINE(1:LL4) 
         END DO 
      END DO 
    endif


!//EFD deallocate storage of temporary flag arrays

    deallocate( flag_e_g )
    deallocate( flag_n_g )
    deallocate( flag_t_g )
    deallocate( flag_temp )





 5050 FORMAT(3X,'J',3X,'I=',3X,30(I3,1X)) 
 5100 FORMAT(1X,I3,8X,30(A3,1X)) 
      RETURN  
      END SUBROUTINE OUT_ARRAY_KC 
