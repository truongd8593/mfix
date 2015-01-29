!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is     C
!           user-definable.  The user may insert code in this routine  C 
!           or call appropriate user defined subroutines.              C 
!           This routine is not called from an IJK loop, hence         C
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
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
      SUBROUTINE USR3 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!


        CALL calculate_de_norms
      

      RETURN 

      END SUBROUTINE USR3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_de_norms	                               !
!  Purpose: Calculates discretization error norms for the problem      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jan 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_de_norms
      Use param1, only    : zero      
      Use compar, only    : ijkstart3, ijkend3, myPE    
      Use functions, only : i_of, j_of, k_of
      Use functions, only : IS_ON_myPE_owns
      Use usr, only       : de_p_g, de_u_g, de_v_g      
      Use usr, only       : p_g_ex, u_g_ex, v_g_ex
      Use usr, only       : lnorms_p_g, lnorms_u_g, lnorms_v_g  
      Use fldvar, only    : p_g, u_g, v_g 
      Use geometry, only  : imax, imin1, imax1
      Use geometry, only  : jmax, jmin1, jmax1
      Use funits, only    : newunit
      Use mpi_utility
      IMPLICIT NONE

! number of data points for norm calculation
        integer   :: var_size

! looping variables        
        integer   :: ijk, i, j, k

! file unit        
        integer   :: f1


        call calculate_exact_solution_channel

! scalar variable (pressure)
        de_p_g = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle

          de_p_g(ijk) = p_g(ijk) - p_g_ex(ijk)
        end do

        var_size = imax*jmax
        call calculate_lnorms(var_size, de_p_g, lnorms_p_g)

! vector variable (x)
        de_u_g = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.ge.imax1)) cycle  ! note: >=
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle

          de_u_g(ijk) = u_g(ijk) - u_g_ex(ijk)
        end do

        var_size = (imax-1)*jmax
        call calculate_lnorms(var_size, de_u_g, lnorms_u_g)


! vector variable (x)
        de_v_g = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.ge.jmax1)) cycle  ! note: >=

          de_v_g(ijk) = v_g(ijk) - v_g_ex(ijk)
        end do

        var_size = imax*(jmax-1)
        call calculate_lnorms(var_size, de_v_g, lnorms_v_g)


! Output DE norms data to a file
        if(myPE == PE_IO) then
          open(unit=newunit(f1), file="de_norms.dat", status='unknown')
          write(f1,*) "# DE Norms for Horizontal Channel Flow:"
          write(f1,*) "# imax= ",imax, " jmax=", jmax
          write(f1,*) "# 1st line: L1 Norms, 2nd line: L2 Norms, &
                       &3rd line: Linf Norms"
          write(f1,*) "# Columns: P_g : U_g : V_g"
          write(f1,*) lnorms_p_g(1), lnorms_u_g(1), lnorms_v_g(1)
          write(f1,*) lnorms_p_g(2), lnorms_u_g(2), lnorms_v_g(2)
          write(f1,*) lnorms_p_g(3), lnorms_u_g(3), lnorms_v_g(3)
          close(f1)
        end if
        
      RETURN

      END SUBROUTINE calculate_de_norms


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_lnorms	                               !
!  Purpose: Calculates L1, L2 and Lifinity norms for any input         !
!  variable of defined size                                            ! 
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jan 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_lnorms(var_size, var, lnorms)
      Use compar, only      : ijkstart3, ijkend3      
      Use param, only       : dimension_3
      Use param1, only      : zero
      Use mpi_utility, only : global_sum, global_max
      IMPLICIT NONE

! total number of data points for norm calculation
        integer, intent(in)                             :: var_size

! variable for which norms are to be calculated      
        double precision, dimension(dimension_3), intent(in)   :: var 

! the three norms: L1, L2 and Linfinity        
        double precision, dimension(3), intent(out)     :: lnorms

! sum or max of lnorms (over all processors
        double precision, dimension(3)                  :: lnorms_all

! looping indices        
        integer         :: ijk


! calculate L1, L2 and Linfinity norms        
        lnorms(1:3) = zero
        do ijk = ijkstart3, ijkend3
            lnorms(1) = lnorms(1) + abs(var(ijk))
            lnorms(2) = lnorms(2) + abs(var(ijk))**2
            lnorms(3) = max(lnorms(3), abs(var(ijk)))
        end do

! save global sum in lnorms_all variables        
        call global_sum(lnorms(1), lnorms_all(1))
        call global_sum(lnorms(2), lnorms_all(2))
        call global_max(lnorms(3), lnorms_all(3))

! put final result in lnorms         
        lnorms(1) = lnorms_all(1)/real(var_size)
        lnorms(2) = sqrt(lnorms_all(2)/real(var_size))
        lnorms(3) = lnorms_all(3)


      RETURN

      END SUBROUTINE calculate_lnorms


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_exact_solution_channel                       !
!  Purpose: Calculates exact solution for the channel problem          !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jan 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_exact_solution_channel
      Use compar, only    : ijkstart3, ijkend3      
      Use functions, only : i_of, j_of, k_of
      Use geometry, only  : dx, dy, dz, imax1
      Use usr, only       : xtr, ytr, ztr
      Use usr, only       : p_g_ex, u_g_ex, v_g_ex
      Use param1, only    : zero, half
      IMPLICIT NONE

! user-defined channel problem parameters 
! (must match the input in mfix.dat)
! (better to set these up via an input file in future) 
        double precision, parameter :: mu_usr     = 0.001d0   ! viscosity
        double precision, parameter :: dpdx_usr   = 1200.d0   ! dp/dx
        double precision, parameter :: h_usr      = 0.01d0    ! height
        double precision, parameter :: l_usr      = 0.2d0     ! length
        double precision, parameter :: p0_usr     = 101325.d0 ! ref. pressure
      
! looping indices      
        integer :: ijk, i, j, k

        integer :: ii, jj, kk

! temporary variables        
        double precision  :: xt, yt, zt

        
! generate grid locations for exact solution calculation 
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          xt = zero - dx(1)
          yt = zero - dy(1) 
          zt = zero - dz(1)

          do ii = 1, i
            xt = xt + dx(ii)
          end do

          do jj = 1, j
            yt = yt + dy(jj)
          end do

          do kk = 1, k
            zt = zt + dz(kk)
          end do

          xtr(ijk) = xt
          ytr(ijk) = yt
          ztr(ijk) = zt 
        end do

! set exact solution
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)

! scalar variable (pressure)
          xt = xtr(ijk) - dx(i)*half
          p_g_ex(ijk) = p0_usr + &
              dpdx_usr*((l_usr - half*dx(imax1)) - xt) 

! vector variable (x)
          yt = ytr(ijk) - dy(j)*half         
          u_g_ex(ijk) = half/mu_usr*dpdx_usr*yt*(h_usr-yt) 

! vector variable (y)          
          v_g_ex(ijk) = zero
        end do


      RETURN

      END SUBROUTINE calculate_exact_solution_channel
