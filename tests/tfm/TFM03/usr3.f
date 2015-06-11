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


        if(tecplot_output) then
          CALL generate_grid_locations
          CALL write_tecplot_data
        endif
        
        CALL deallocate_usr_variables


      RETURN

      END SUBROUTINE USR3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: generate_grid_locations                                !
!  Purpose: Generate grid node location arrays from dx, dy dz          !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE generate_grid_locations
      Use compar, only    : ijkstart3, ijkend3
      Use indices, only   : i_of, j_of, k_of
      Use geometry, only  : dx, dy, dz
      Use usr, only       : xtr, ytr, ztr
      Use param1, only    : zero
      IMPLICIT NONE

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


      END SUBROUTINE generate_grid_locations

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: read_literature_results                                !
!  Purpose: read results for centerline velocities from literature     !  
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE read_literature_results
      Use compar, only      : ijkstart3, ijkend3
      Use indices, only     : i_of, j_of, k_of
      Use funits, only      : newunit
      Use usr, only         : npts
      Use usr, only         : y_loc, x_loc
      Use usr, only         : u_ghia, v_ghia
      use compar, only      : myPE, PE_IO      
      IMPLICIT NONE

! indices
        integer           :: i

! file units
        integer           :: fug, fvg

! temporary variables
        double precision  :: xtmp, ytmp, utmp, vtmp

! read from PE_IO processor
        if(myPE==PE_IO) then

          ! read vertical centerline data 
          open(unit=newunit(fug), file="ghia_results_u.dat", &
                                  status='old')
          do i = 1,3 
            read(fug,*) ! read and discard first 3 lines
          end do

          do i = 1,npts
            read(fug,*) xtmp, y_loc(npts-i+1), u_ghia(npts-i+1), utmp
          end do

          close(fug)

          ! read horizontal centerline data 
          open(unit=newunit(fvg), file="ghia_results_v.dat", &
                                  status='old')
          do i = 1,3 
            read(fvg,*) ! read and discard first 3 lines
          end do

          do i = 1,npts
            read(fvg,*) x_loc(npts-i+1), ytmp, v_ghia(npts-i+1), vtmp
          end do

          close(fvg)

        end if ! end of if (myPE==PE_IO)


      RETURN

      END SUBROUTINE read_literature_results


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: find_centerline_values                                 !
!  Purpose: Find u_g and v_g at centerlines of the cavity.             !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE find_centerline_values
      use geometry, only    : ijkmax3      
      use geometry, only    : imax, jmax
      use geometry, only    : kmin1
      use geometry, only    : dx, dy, dz
      use functions, only   : funijk_gl
      use compar, only      : myPE, PE_IO      
      use fldvar, only      : u_g, v_g
      use usr, only         : npts
      use usr, only         : y_loc, x_loc
      use usr, only         : u_mfix
      use usr, only         : v_mfix
      use param1, only      : zero, half
      use mpi_utility, only : gather      
      IMPLICIT NONE


! indices
        integer             :: i, j, k, ijk
        integer             :: n

! temporary variables for x, y, z coordinates        
        double precision    :: xt, yt, zt

! arrays to gather variables
        double precision, allocatable :: arr_Ug(:)
        double precision, allocatable :: arr_Vg(:)

        
! allocate variables as ijkmax3 size 1D arrays
        allocate(arr_Ug(ijkmax3))
        allocate(arr_Vg(ijkmax3))

! gather the MFIX global variables into ijkmax3 sized variables
        call gather(U_G,arr_Ug,PE_IO)
        call gather(V_G,arr_Vg,PE_IO)

        if(myPE==PE_IO) then
        
! vertical centerline        
        do n = 1,npts
          xt = half     ! cavity of unit dimensions
          yt = y_loc(n)

          CALL CALC_CELL_INTERSECT(ZERO, xt, DX, IMAX, i)
          CALL CALC_CELL_INTERSECT(ZERO, yt, DY, JMAX, j)
          k = kmin1
     
          ijk = funijk_gl(i,j,k)

          u_mfix(n) = arr_Ug(ijk)

        end do

! horizontal centerline
        do n = 1,npts
          xt = x_loc(n)
          yt = half     ! cavity of unit dimensions

          CALL CALC_CELL_INTERSECT(ZERO, xt, DX, IMAX, i)
          CALL CALC_CELL_INTERSECT(ZERO, yt, DY, JMAX, j)
          k = kmin1
     
          ijk = funijk_gl(i,j,k)

          v_mfix(n) = arr_Vg(ijk)
        end do

      end if ! end of (myPE==PE_IO)

! deallocate local arrays        
      deallocate(arr_Ug)
      deallocate(arr_Vg)


      END SUBROUTINE find_centerline_values

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: write_tecplot_data                                     !
!  Purpose: Write data for visualization in Tecplot.                   !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
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
      use fldvar, only      : p_g, u_g, v_g
      use usr, only         : npts
      use usr, only         : y_loc, x_loc
      use usr, only         : u_mfix, u_ghia
      use usr, only         : v_mfix, v_ghia
      use usr, only         : xtr, ytr, ztr
      use usr, only         : tec_output_block, tec_no_k, &
                              uv_profile
      use param1, only      : half
      use mpi_utility, only : gather
      IMPLICIT NONE

! indices
        integer             :: i, j, k, ijk
        integer             :: imjk, ijmk, ijkm

! file units
        integer             :: ftcc, fup, fvp

! arrays to gather variables
        double precision, allocatable :: arr_Pg(:)
        double precision, allocatable :: arr_Ug(:)
        double precision, allocatable :: arr_Vg(:)
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
                                          Pg_tmp, &
                                          Ug_tmp, &
                                          Vg_tmp



! allocate variables as ijkmax3 size 1D arrays
        allocate(arr_Pg(ijkmax3))
        allocate(arr_Ug(ijkmax3))
        allocate(arr_Vg(ijkmax3))
        allocate(arr_xtr(ijkmax3))
        allocate(arr_ytr(ijkmax3))
        allocate(arr_ztr(ijkmax3))

! gather the MFIX global variables into ijkmax3 sized variables
        call gather(P_G,arr_Pg,PE_IO)
        call gather(U_G,arr_Ug,PE_IO)
        call gather(V_G,arr_Vg,PE_IO)
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

            Ug_tmp(i,j,k) = half*(arr_Ug(ijk)+arr_Ug(imjk))

            Vg_tmp(i,j,k) = half*(arr_Vg(ijk)+arr_Vg(ijmk))

          end do
          end do
          end do

! write true tecplot cell-centered data (block format)
! This is the traditional tecplot cell-centered data. Be wary while
! using contour plots which could do inaccurate interpolation near the
! boundaries. It's better to use tecplot's "primary flood" with this
! option.

          if(tec_output_block) then

            open(unit=newunit(ftcc), file="solution_tec_block.dat", &
                                      status='unknown')

            write(ftcc,"(6a)") 'variables = "x""y""z"&
                &"Pg""Ug""Vg"'
            write(ftcc,*) 'zone T="',0,'" '
            write(ftcc,*) 'I=',IMAX1,' J=',JMAX1,' K=',KMAX1

            write(ftcc,*) 'DATAPACKING=BLOCK'
            write(ftcc,*) "VARLOCATION=([4-6]=CELLCENTERED)"

            ! x coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(ftcc,*) (x_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! y coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(ftcc,*) (y_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! z coordinates
            do k=kmin2,kmax1
            do j=jmin2,jmax1
              write(ftcc,*) (z_tmp(i,j,k),i=imin2,imax1)
            end do
            end do

            ! solution variables
            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Pg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Ug_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            do k=kmin1,kmax1
            do j=jmin1,jmax1
              write(ftcc,*) (Vg_tmp(i,j,k),i=imin1,imax1)
            end do
            end do

            close(ftcc)

          end if ! end of if(tec_output_block)

! write tecplot data (point format) for:
! x-velocity profile at x=L/2 (vertical centerline), and
! y-velocity profile at y=H/2 (horizonontal centerline).
! Also included are results from literature.          

          if(uv_profile) then

            call read_literature_results

            call find_centerline_values

            open(unit=newunit(fup), &
             file="u_profile.dat", status='unknown')

            write(fup,"(3a)") 'variables = "y""Ug (MFIX)""Ghia"'
            write(fup,*) 'zone T="',0,'" '
            write(fup,*) 'I= ', npts

            do i = 1,npts
              write(fup,*) y_loc(i), u_mfix(i), u_ghia(i)  
            end do

            close(fup)

            open(unit=newunit(fvp), &
             file="v_profile.dat", status='unknown')

            write(fvp,"(3a)") 'variables = "x""Vg (MFIX)""Ghia"'
            write(fvp,*) 'zone T="',0,'" '
            write(fvp,*) 'I= ', npts

            do i = 1,npts
              write(fvp,*) x_loc(i), v_mfix(i), v_ghia(i)  
            end do

            close(fvp)

          endif ! end of if(uv_profile)


        end if ! end of if(myPE==PE_IO)

! deallocate local arrays        
        deallocate(arr_Pg)
        deallocate(arr_Ug)
        deallocate(arr_Vg)
        deallocate(arr_xtr)
        deallocate(arr_ytr)
        deallocate(arr_ztr)


      RETURN

      END SUBROUTINE write_tecplot_data


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: deallocate_usr_variables                               !
!  Purpose: Deallocate allocatable variables defined in usr_mod.f      !
!                                                                      !
!  Author: Aniruddha Choudhary                        Date: Jun 2015   !
!  email: anirudd@vt.edu					       !
!  Reviewer:                                          Date:            !
!                                                                      !
!  Revision Number #                                  Date:            !
!  Author: #                                                           !
!  Purpose: #                                                          !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE deallocate_usr_variables
      Use usr, only     : xtr, ytr, ztr
      IMPLICIT NONE

! deallocate variables defined in usr_mod.f
        deallocate(xtr)
        deallocate(ytr)
        deallocate(ztr)


      RETURN

      END SUBROUTINE deallocate_usr_variables
