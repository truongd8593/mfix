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
      use usr, only         : tecplot_output
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

        if(tecplot_output) CALL write_tecplot_data

        CALL deallocate_usr_variables


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
      Use indices, only   : i_of, j_of, k_of
      Use functions, only : IS_ON_myPE_owns
      Use usr, only       : de_P_g, de_u_g, de_v_g
      Use usr, only       : P_g_ex, u_g_ex, v_g_ex
      Use usr, only       : lnorms_P_g, lnorms_u_g, lnorms_v_g
      Use fldvar, only    : P_g, u_g, v_g
      Use geometry, only  : imax, imin1, imax1
      Use geometry, only  : jmax, jmin1, jmax1
      Use funits, only    : newunit
      Use compar, only    : myPE, PE_IO
      IMPLICIT NONE

! number of data points for norm calculation
        integer   :: var_size

! looping variables
        integer   :: ijk, i, j, k

! file unit
        integer   :: f1


        call calculate_exact_solution_channel

! scalar variable (pressure)
        de_P_g = zero

        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)
          k = k_of(ijk)

          if(.NOT.IS_ON_myPE_owns(i,j,k)) cycle

! to select only internal points, else remains zero
          if((i.lt.imin1).or.(i.gt.imax1)) cycle
          if((j.lt.jmin1).or.(j.gt.jmax1)) cycle

          de_P_g(ijk) = P_g(ijk) - P_g_ex(ijk)
        end do

        var_size = imax*jmax
        call calculate_lnorms(var_size, de_P_g, lnorms_P_g)

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
          write(f1,*) lnorms_P_g(1), lnorms_u_g(1), lnorms_v_g(1)
          write(f1,*) lnorms_P_g(2), lnorms_u_g(2), lnorms_v_g(2)
          write(f1,*) lnorms_P_g(3), lnorms_u_g(3), lnorms_v_g(3)
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
        lnorms(1) = lnorms_all(1)/dble(var_size)
        lnorms(2) = sqrt(lnorms_all(2)/dble(var_size))
        lnorms(3) = lnorms_all(3)


      RETURN

      END SUBROUTINE calculate_lnorms

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      		!
!  Module name: calculate_exact_solution_channel                       		!
!  Purpose: Calculates exact solution for the Couette flow problem with 	!
!           favorable pressure gradient						!
!                                                                      		!
!  Author: Avinash Vaidheeswaran                        Date: July 2016		!
!  email: avinash.vaidheeswaran@netl.doe.gov			     		!
!  Reviewer:                                          Date:            		!
!                                                                      		!
!  Revision Number #                                  Date:            		!
!  Author: #                                                           		!
!  Purpose: #                                                          		!
!                                                                      		!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE calculate_exact_solution_channel
      Use compar, only    : ijkstart3, ijkend3
      Use indices, only   : i_of, j_of, k_of
      Use geometry, only  : dx, dy, dz, imax1
      Use usr, only       : xtr, ytr, ztr
      Use usr, only       : P_g_ex, u_g_ex, v_g_ex
      Use param1, only    : zero, half
      IMPLICIT NONE

! user-defined channel problem parameters
! (must match the input in mfix.dat)
! (better to set these up via an input file in future)
        double precision, parameter :: mu_usr     = 0.0000157d0 ! viscosity
        double precision, parameter :: dpdx_usr   = 0.0d0 ! p gradient
        double precision, parameter :: l_usr      = 0.1d0   ! length
        double precision, parameter :: h_usr      = 0.01d0   ! height
        double precision, parameter :: U_usr      = 10.0d0        ! velocity of N wall

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

          yt = ytr(ijk) - dy(j)*half

! vector variable (x)
          u_g_ex(ijk) = (1./2./mu_usr)*(dpdx_usr)*(yt**2.) &
                        + ((U_usr/h_usr) - (1./2./mu_usr)*(dpdx_usr)*h_usr)*yt

! vector variable (y)
          v_g_ex(ijk) = zero
        end do


      RETURN

      END SUBROUTINE calculate_exact_solution_channel


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: deallocate_usr_variables                               !
!  Purpose: Deallocate allocatable variables defined in usr_mod.f      !
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

      SUBROUTINE deallocate_usr_variables
      Use usr, only     : P_g_ex, u_g_ex, v_g_ex
      Use usr, only     : lnorms_P_g, lnorms_u_g, lnorms_v_g
      Use usr, only     : xtr, ytr, ztr
      Use usr, only     : de_P_g, de_u_g, de_v_g
      IMPLICIT NONE

        ! allocate variables defined in usr_mod.f
        deallocate(P_g_ex)
        deallocate(u_g_ex)
        deallocate(v_g_ex)

        deallocate(lnorms_P_g)
        deallocate(lnorms_u_g)
        deallocate(lnorms_v_g)

        deallocate(xtr)
        deallocate(ytr)
        deallocate(ztr)

        deallocate(de_P_g)
        deallocate(de_u_g)
        deallocate(de_v_g)

      RETURN

      END SUBROUTINE deallocate_usr_variables


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: write_tecplot_data                                     !
!  Purpose: Write data for visualization in Tecplot.                   !
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

      SUBROUTINE write_tecplot_data
      use geometry, only    : ijkmax3
      use geometry, only    : imin2, imin1, imax1, imax2
      use geometry, only    : jmin2, jmin1, jmax1, jmax2
      use geometry, only    : kmin2, kmin1, kmax1, kmax2
      use geometry, only    : imax, jmax, kmax
      use geometry, only    : dx, dy, dz
      use compar, only      : myPE, PE_IO
      use funits, only      : newunit
      use functions, only   : funijk_gl
      use fldvar, only      : P_g, u_g, v_g
      use usr, only         : P_g_ex, u_g_ex, v_g_ex
      use usr, only         : xtr, ytr, ztr
      use usr, only         : tec_output_block, tec_no_k, &
                              x_velocity_profile, P_profile
      use param1, only      : half
      use mpi_utility, only : gather
      IMPLICIT NONE

! indices
        integer             :: i, j, k, ijk
        integer             :: imjk, ijmk, ijkm

! file units
        integer             :: ftcc, fxvp1, fxvp2

! arrays to gather variables
        double precision, allocatable :: arr_Pg(:)
        double precision, allocatable :: arr_Pgex(:)
        double precision, allocatable :: arr_Ug(:)
        double precision, allocatable :: arr_Ugex(:)
        double precision, allocatable :: arr_Vg(:)
        double precision, allocatable :: arr_Vgex(:)
        double precision, allocatable :: arr_xtr(:)
        double precision, allocatable :: arr_ytr(:)
        double precision, allocatable :: arr_ztr(:)

! temporary variables for x,y,z coordinates
        double precision              :: xt, yt, zt

! x,y,z coordinate variables (3D arrays)
        double precision, &
          dimension(imin2:imax1,jmin2:jmax1,kmin2:kmax1)  :: &
                                          x_tmp, y_tmp, z_tmp

! output variables (3D arrays)
        double precision, &
          dimension(imin2:imax2,jmin2:jmax2,kmin2:kmax2)  :: &
                                          Pg_tmp, Pgex_tmp, &
                                          Ug_tmp, Ugex_tmp, &
                                          Vg_tmp, Vgex_tmp



! allocate variables as ijkmax3 size 1D arrays
        allocate(arr_Pg(ijkmax3))
        allocate(arr_Pgex(ijkmax3))
        allocate(arr_Ug(ijkmax3))
        allocate(arr_Ugex(ijkmax3))
        allocate(arr_Vg(ijkmax3))
        allocate(arr_Vgex(ijkmax3))
        allocate(arr_xtr(ijkmax3))
        allocate(arr_ytr(ijkmax3))
        allocate(arr_ztr(ijkmax3))

! gather the MFIX global variables into ijkmax3 sized variables
        call gather(P_g,arr_Pg,PE_IO)
        call gather(P_g_ex,arr_Pgex,PE_IO)
        call gather(U_G,arr_Ug,PE_IO)
        call gather(u_g_ex,arr_Ugex,PE_IO)
        call gather(V_G,arr_Vg,PE_IO)
        call gather(v_g_ex,arr_Vgex,PE_IO)
        call gather(xtr,arr_xtr,PE_IO)
        call gather(ytr,arr_ytr,PE_IO)
        call gather(ztr,arr_ztr,PE_IO)

! write out data from PE_IO processor
        if(myPE==PE_IO) then

! create 3D arrays for cell-centered data

! create 3D arrays for mesh nodes
          do k = kmin2, kmax1
          do j = jmin2, jmax1
          do i = imin2, imax1
             ijk = funijk_gl(i,j,k)
             x_tmp(I,J,K) = arr_xtr(IJK)
             y_tmp(I,J,K) = arr_ytr(IJK)
             z_tmp(I,J,K) = arr_ztr(IJK)
          end do
          end do
          end do

! create 3D arrays for solution variables at cell-centers
          do k = kmin1, kmax1
          do j = jmin1, jmax1
          do i = imin1, imax1

            ijk = funijk_gl(i,j,k)
            imjk = funijk_gl(i-1,j,k)
            ijmk = funijk_gl(i,j-1,k)
            if(tec_no_k) then
              ijkm = funijk_gl(i,j,k)
            else
              ijkm = funijk_gl(i,j,k-1)
            end if

            Pg_tmp(i,j,k) = arr_Pg(ijk)
            Pgex_tmp(i,j,k) = arr_Pgex(ijk)

            Ug_tmp(i,j,k) = half*(arr_Ug(ijk)+arr_Ug(imjk))
            Ugex_tmp(i,j,k) = half*(arr_Ugex(ijk)+arr_Ugex(imjk))

            Vg_tmp(i,j,k) = half*(arr_Vg(ijk)+arr_Vg(ijmk))
            Vgex_tmp(i,j,k) = half*(arr_Vgex(ijk)+arr_Vgex(ijmk))

          end do
          end do
          end do

! write tecplot data (point format) for the x-velocity profile at x=L/2:
! - numerical solution for u_g          
! - exact solution for u_g
! - error for u_g

          if(x_velocity_profile) then

            open(unit=newunit(fxvp1), &
             file="solution_x_velocity_profile.dat", status='unknown')

            write(fxvp1,"(6a)") 'variables = "x""y""z""Ug""Ugex""UgErr"'
            write(fxvp1,*) 'zone T="',0,'" '
            write(fxvp1,*) 'J= ',JMAX

            i = imax/2 + 1
            k = kmin1
            do j = jmin1, jmax1

              ijk = funijk_gl(i,j,k)

              xt = arr_xtr(ijk)
              yt = arr_ytr(ijk) - dy(j)*half
              zt = arr_ztr(ijk) - dz(k)*half

              Ug_tmp(I,J,K) = arr_Ug(IJK)
              Ugex_tmp(I,J,K) = arr_Ugex(IJK)

              write(fxvp1,10) xt, yt, zt, &
                Ug_tmp(I,J,K), &
                Ugex_tmp(I,J,K), &
                Ug_tmp(I,J,K) - Ugex_tmp(I,J,K) 

          10 format(f8.5,2x,f8.5,2x,f8.5,2x,f16.12,2x,f16.12,2x,f19.15)

            end do

          endif ! end of if(x_velocity_profile)

         end if ! end of if(myPE==PE_IO)



      RETURN

      END SUBROUTINE write_tecplot_data
