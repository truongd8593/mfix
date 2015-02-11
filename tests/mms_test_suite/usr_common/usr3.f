!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is
!           user-definable.  The user may insert code in this routine
!           or call appropriate user defined subroutines.
!           This routine is not called from an IJK loop, hence
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
      Use mms, only         : deallocate_mms_vars
      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
      logical               :: tecplot_output = .FALSE.
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!

      
      Call calculate_de_norms

      !if(tecplot_output) Call write_tecplot_data 

      Call deallocate_mms_vars


      RETURN
      END SUBROUTINE USR3



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: calculate_de_norms	                               !
!  Purpose: Calculates discretization error norms for the problem      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Feb 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      ! 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_de_norms
      Use param, only     : dimension_3
      Use param1, only    : zero      
      Use compar, only    : ijkstart3, ijkend3, istart1_all, iend1_all,&
                            jstart1_all, jend1_all, &
                            kstart1_all, kend1_all, &
                            ijkstart3_all, ijkend3_all, &
                            istart, iend, jstart, jend, kstart, kend
      Use functions, only : funijk, funijk_gl                          
      Use functions, only : IS_ON_myPE_owns
      Use indices, only   : i_of, j_of, k_of
      Use usr             ! all variables
      Use mms             ! all variables
      Use fldvar, only    : p_g, u_g, v_g, w_g, u_s, v_s, w_s, &
                            t_g, t_s, ep_g, rop_s, theta_m
      Use geometry, only  : imax, imin1, imax1
      Use geometry, only  : jmax, jmin1, jmax1
      Use geometry, only  : kmax, kmin1, kmax1     
      Use funits, only    : newunit
      Use compar, only    : myPE, PE_IO, numPEs
      Use mpi_utility, only : bcast
      IMPLICIT NONE

! number of data points for norm calculation
        integer   :: var_size

! looping variables        
        integer   :: ijk, i, j, k, proc

! file unit        
        integer   :: f1

! pressure shift value        
        double precision  :: delta_p_g        

! owner of ijk_sh
        integer   :: ijk_sh_owner


! allocate de_ variables and lnorms_ variables
        allocate(de_ep_g(dimension_3))
        allocate(de_p_g(dimension_3))
        allocate(de_u_g(dimension_3))
        allocate(de_v_g(dimension_3))
        allocate(de_w_g(dimension_3))
        allocate(de_t_g(dimension_3))
        allocate(de_rop_s(dimension_3))
        allocate(de_u_s(dimension_3))
        allocate(de_v_s(dimension_3))
        allocate(de_w_s(dimension_3))
        allocate(de_t_s(dimension_3))
        allocate(de_theta_m(dimension_3))

        allocate(lnorms_ep_g(3))
        allocate(lnorms_p_g(3))
        allocate(lnorms_u_g(3))
        allocate(lnorms_v_g(3))
        allocate(lnorms_w_g(3))
        allocate(lnorms_t_g(3))
        allocate(lnorms_rop_s(3))
        allocate(lnorms_u_s(3))
        allocate(lnorms_v_s(3))
        allocate(lnorms_w_s(3))
        allocate(lnorms_t_s(3))
        allocate(lnorms_theta_m(3))

!! evaluate pressure shift value
!        do proc = 0, numPEs-1
!            
!          i = i_of(ijk_sh)
!          j = j_of(ijk_sh)
!          k = k_of(ijk_sh)
!
!          write(*,*) proc, i, istart1_all(proc), iend1_all(proc)
!          write(*,*) proc, j, jstart1_all(proc), jend1_all(proc)
!          write(*,*) proc, k, kstart1_all(proc), kend1_all(proc)
!
!! *** ask jmusser
!! is there a function that check whether an ijk is owned by a
!! processor 'proc'?
!! istart_all for all 'proc' holds the same value? is this a problem?
!          if (.not.((i.ge.istart1_all(proc)).and.(i.le.iend1_all(proc))))&
!            cycle
!          if (.not.((j.ge.jstart1_all(proc)).and.(j.le.jend1_all(proc))))&
!            cycle
!          if (.not.((k.ge.kstart1_all(proc)).and.(k.le.kend1_all(proc))))&
!            cycle
!
!! 'proc' has ijk_sh
!          ijk_sh_owner = proc
!
!! if 'proc' is also the current process then set delta_p_g else set some
!! random value
!          if(IS_ON_myPE_owns(i,j,k)) then
!            delta_p_g = P_G(ijk_sh) - MMS_P_G(ijk_sh)      
!! only ijk_sh_owner broadcasts the value of delta_p_g
!            call bcast(delta_p_g)
!          else
!            delta_p_g = zero
!          end if            
!
!        end do

        do ijk = ijkstart3, ijkend3

          i = i_of(ijk_sh)
          j = j_of(ijk_sh)
          k = k_of(ijk_sh)
          
          if(IS_ON_myPE_owns(i,j,k)) then
            write(*,*) "myPE= ", myPE, ijk_sh, i, j, k
            delta_p_g = P_G(ijk_sh) - MMS_P_G(ijk_sh)      
! only ijk_sh_owner broadcasts the value of delta_p_g
            call bcast(delta_p_g)
            exit
          end if

        end do


        write(*,*) "myPE, dpg", myPE, delta_p_g

! scalar variables
        de_ep_g = zero
        de_p_g = zero
        de_t_g = zero
        de_rop_s = zero
        de_t_s = zero
        de_theta_m = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle
          if((k.lt.kmin1).or.(k.gt.kmax1)) cycle

          de_ep_g(ijk)    = ep_g(ijk)     - mms_ep_g(ijk)
          de_p_g(ijk)     = p_g(ijk)      - mms_p_g(ijk) - delta_p_g
          de_t_g(ijk)     = t_g(ijk)      - mms_t_g(ijk)
          de_rop_s(ijk)   = rop_s(ijk,1)  - mms_rop_s(ijk)
          de_t_s(ijk)     = t_s(ijk,1)    - mms_t_s(ijk)
          de_theta_m(ijk) = theta_m(ijk,1)- mms_theta_m(ijk)

        end do

        var_size = imax*jmax*kmax
        call calculate_lnorms(var_size, de_ep_g, lnorms_ep_g)
        call calculate_lnorms(var_size, de_p_g, lnorms_p_g)
        call calculate_lnorms(var_size, de_t_g, lnorms_t_g)
        call calculate_lnorms(var_size, de_rop_s, lnorms_rop_s)
        call calculate_lnorms(var_size, de_t_s, lnorms_t_s)
        call calculate_lnorms(var_size, de_theta_m, lnorms_theta_m)

! vector variables (x)
        de_u_g = zero
        de_u_s = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.ge.imax1)) cycle  ! note: >=
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle
          if((k.lt.kmin1).or.(k.gt.kmax1)) cycle

          de_u_g(ijk) = u_g(ijk) - mms_u_g(ijk)
          de_u_s(ijk) = u_s(ijk,1) - mms_u_s(ijk)
        end do

        var_size = (imax-1)*jmax*kmax
        call calculate_lnorms(var_size, de_u_g, lnorms_u_g)
        call calculate_lnorms(var_size, de_u_s, lnorms_u_s)

! vector variables (y)
        de_v_g = zero
        de_v_s = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.ge.jmax1)) cycle  ! note: >=
          if((k.lt.kmin1).or.(k.gt.kmax1)) cycle

          de_v_g(ijk) = v_g(ijk) - mms_v_g(ijk)
          de_v_s(ijk) = v_s(ijk,1) - mms_v_s(ijk)
        end do

        var_size = imax*(jmax-1)*kmax
        call calculate_lnorms(var_size, de_v_g, lnorms_v_g)
        call calculate_lnorms(var_size, de_v_s, lnorms_v_s)

! vector variables (z)
        de_w_g = zero
        de_w_s = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle        

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle
          if((k.lt.kmin1).or.(k.ge.kmax1)) cycle  ! note: >=

          de_w_g(ijk) = w_g(ijk) - mms_w_g(ijk)
          de_w_s(ijk) = w_s(ijk,1) - mms_w_s(ijk)
        end do

        var_size = imax*jmax*(kmax-1)
        call calculate_lnorms(var_size, de_w_g, lnorms_w_g)
        call calculate_lnorms(var_size, de_w_s, lnorms_w_s)

! Output DE norms data to a file
        if(myPE == PE_IO) then
          open(unit=newunit(f1), file="de_norms.dat", status='unknown')
          write(f1,*) "# DE Norms for MMS test cases:"
          write(f1,*) "# imax= ",imax, " jmax=", jmax, " kmax=", kmax
          write(f1,*) "# 1st line: L1 Norms, 2nd line: L2 Norms, &
                       &3rd line: Linf Norms"
          write(f1,*) "# variables='p_g''u_g''v_g''w_g''u_s''v_s''w_s'&
                      &'t_g''t_s''ep_g''rop_s''theta_m'"
          write(f1,*) lnorms_p_g(1), &
            lnorms_u_g(1), lnorms_v_g(1), lnorms_w_g(1), &
            lnorms_u_s(1), lnorms_v_s(1), lnorms_w_s(1), &
            lnorms_t_g(1), lnorms_t_s(1), &
            lnorms_ep_g(1), lnorms_rop_s(1), &
            lnorms_theta_m(1)
          write(f1,*) lnorms_p_g(2), &
            lnorms_u_g(2), lnorms_v_g(2), lnorms_w_g(2), &
            lnorms_u_s(2), lnorms_v_s(2), lnorms_w_s(2), &
            lnorms_t_g(2), lnorms_t_s(2), &
            lnorms_ep_g(2), lnorms_rop_s(2), &
            lnorms_theta_m(2)
          write(f1,*) lnorms_p_g(3), &
            lnorms_u_g(3), lnorms_v_g(3), lnorms_w_g(3), &
            lnorms_u_s(3), lnorms_v_s(3), lnorms_w_s(3), &
            lnorms_t_g(3), lnorms_t_s(3), &
            lnorms_ep_g(3), lnorms_rop_s(3), &
            lnorms_theta_m(3)
          close(f1)
        end if

! de allocate de_ variables and lnorms_ variables
        deallocate(de_ep_g)
        deallocate(de_p_g)
        deallocate(de_u_g)
        deallocate(de_v_g)
        deallocate(de_w_g)
        deallocate(de_t_g)
        deallocate(de_rop_s)
        deallocate(de_u_s)
        deallocate(de_v_s)
        deallocate(de_w_s)
        deallocate(de_t_s)
        deallocate(de_theta_m)

        deallocate(lnorms_ep_g)
        deallocate(lnorms_p_g)
        deallocate(lnorms_u_g)
        deallocate(lnorms_v_g)
        deallocate(lnorms_w_g)
        deallocate(lnorms_t_g)
        deallocate(lnorms_rop_s)
        deallocate(lnorms_u_s)
        deallocate(lnorms_v_s)
        deallocate(lnorms_w_s)
        deallocate(lnorms_t_s)
        deallocate(lnorms_theta_m)

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


        if(var_size.eq.0) then
          lnorms(:) = zero
          return
        end if        

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
        lnorms(1) = lnorms_all(1)/dble(var_size)
        lnorms(2) = sqrt(lnorms_all(2)/dble(var_size))
        lnorms(3) = lnorms_all(3)


      RETURN

      END SUBROUTINE calculate_lnorms


