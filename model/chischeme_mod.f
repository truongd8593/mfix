MODULE ChiScheme
!  Purpose: Determine factors for Chi-scheme of Darwish and Moukalled  
!           to ensure consistency of equation sets (e.g., species mass 
!           fractions add up to 1)                                                 
!                                                                      
!  Author: M. Syamlal                                 Date: 1-AUG-03   
!  Reviewer:                                          Date:            
!                                                                      
!                                                                      
!  Literature/Document References: Darwish, M. and Moukalled, F., 
!   "The Chi-shemes: a new consistent high-resolution formulation based on
!    the normalized variable methodology," Comput. Methods Appl. Mech. Engrg.,
!    192, 1711-1730 (2003)
!                                                                      

! To initiate Chi-Scheme
!   call set_chi( ...) 
!
! and to terminate Chi-Scheme
!   call unset_chi(ier)
!

      USE param 
      USE param1 
      USE run
      USE geometry
      USE indices
      USE compar
      USE sendrecv
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_e
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_n
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_t

      
   LOGICAL :: ChiScheme_allocated = .false.
   LOGICAL :: Chi_flag = .false.
      
   CONTAINS
      SUBROUTINE Set_Chi(DISCR, PHI, Nmax, U, V, W, IER) 
!
!                      Error index
      INTEGER          IER
!
!                      discretization method
      INTEGER          DISCR
!
!                      Second dimension of Phi array 
      INTEGER          NMax
!
!                      convected quantity
      DOUBLE PRECISION PHI(DIMENSION_3, Nmax)
!
!                      Velocity components
      DOUBLE PRECISION U(DIMENSION_3), V(DIMENSION_3), W(DIMENSION_3)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_e_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_n_temp
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Chi_t_temp
!
!                      index
      INTEGER          IJK, N

      if(.not.ChiScheme_allocated)then
        Allocate( Chi_e(DIMENSION_3) , &
                  Chi_n(DIMENSION_3) , &
                  Chi_t(DIMENSION_3)     )
        ChiScheme_allocated = .true.
      endif

      if(Chi_flag)then
        ! Error: Chi-Scheme is already active.  This routine cannot be called
        ! again before unsetting the flag
	Print *, 'Module ChiScheme: Cannot call Set_Chi again, before Unset_chi'
	call Mfix_Exit(0)    
      
      else
      ! Set Chi_flag to indicate that future calls to calc_Xsi will use
      ! the Chi-Scheme for discretization
        Chi_e = large_number
        Chi_n = large_number
        Chi_t = large_number
        Chi_flag = .true.
      Endif

      Allocate( Chi_e_temp(DIMENSION_3) , &
                Chi_n_temp(DIMENSION_3) , &
                Chi_t_temp(DIMENSION_3)  )
      
!  Start Chi calculations
      DO N = 1, Nmax

        CALL CALC_CHI(DISCR, PHI(1,N), U, V, W, Chi_e_temp, Chi_n_temp, Chi_t_temp,0) 

!!!$omp    parallel do private(IJK)
        DO IJK = ijkstart3, ijkend3
          Chi_e(IJK) = MIN(Chi_e(IJK), Chi_e_temp(IJK))
          Chi_n(IJK) = MIN(Chi_n(IJK), Chi_n_temp(IJK))
          Chi_t(IJK) = MIN(Chi_t(IJK), Chi_t_temp(IJK))
	ENDDO
	
      ENDDO

!  End Chi calculations

      call send_recv(CHI_E,2)
      call send_recv(CHI_N,2)
      call send_recv(CHI_T,2)

      Deallocate( Chi_e_temp , &
                  Chi_n_temp , &
                  Chi_t_temp  )


      RETURN  
      END SUBROUTINE Set_Chi
      
      
      SUBROUTINE Unset_Chi(IER) 
        INTEGER          IER
        Chi_flag = .false.
      RETURN  
      END SUBROUTINE Unset_Chi 
      
END MODULE CHIScheme                                                      

