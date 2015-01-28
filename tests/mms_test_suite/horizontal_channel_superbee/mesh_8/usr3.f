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
      logical     :: tecplot_output = .TRUE.
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!

          CALL calculate_de_norms
      
          if(tecplot_output) then

            !CALL write_tecplot_data            
            continue

          endif


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
      Use functions, only : funijk_gl
      Use geometry, only  : imax, jmax
      Use fldvar, only    : p_g, u_g, v_g  
      Use usr, only       : p_g_ex, u_g_ex, v_g_ex
      Use usr, only       : lnorms_p_g, lnorms_u_g, lnorms_v_g     
      Use funits, only    : newunit
      Use mpi_utility     
      IMPLICIT NONE

! temporary array for discretization error: de = p_g - p_g_ex ; etc.        
        double precision, allocatable         :: de_p_g(:)
        double precision, allocatable         :: de_u_g(:)
        double precision, allocatable         :: de_v_g(:)

! file unit        
        integer   :: f1

! looping variables        
        integer   :: ijk, i, j, k
        integer   :: counter

        
        call calculate_exact_solution_channel

        allocate(de_p_g((imax)*(jmax)))
        allocate(de_u_g((imax-1)*(jmax)))
        allocate(de_v_g((imax)*(jmax-1)))

        k = kmin1

        counter = 0
        do j = 2, jmax
        do i = 2, imax
          counter = counter + 1          
          ijk = funijk_gl(i,j,k)
          de_p_g(counter) = p_g(ijk) - p_g_ex(ijk)
        end do
        end do
        call calculate_lnorms(counter, de_p_g, lnorms_p_g) 

        counter = 0
        do j = 2, jmax
        do i = 2, imax-1
          counter = counter + 1          
          ijk = funijk_gl(i,j,k)
          de_u_g(counter) = u_g(ijk) - u_g_ex(ijk)
        end do
        end do
        call calculate_lnorms(counter, de_u_g, lnorms_u_g) 

        counter = 0
        do j = 2, jmax-1
        do i = 2, imax
          counter = counter + 1          
          ijk = funijk_gl(i,j,k)
          de_v_g(counter) = v_g(ijk) - v_g_ex(ijk)
        end do
        end do
        call calculate_lnorms(counter, de_v_g, lnorms_v_g) 

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
      Use param1, only    : zero
      IMPLICIT NONE

! size of the variable
        integer, intent(in)                               :: var_size

! variable for which norms are to be calculated      
        double precision, dimension(var_size), intent(in)  :: var 

! the three norms: L1, L2 and Linfinity        
        double precision, dimension(3), intent(out)       :: lnorms

! looping indices        
        integer         :: ijk


! calculate L1, L2 and Linfinity norms        
        lnorms(1:3) = zero
        do ijk = 1, var_size
            lnorms(1) = lnorms(1) + abs(var(ijk))
            lnorms(2) = lnorms(2) + abs(var(ijk))**2
            lnorms(3) = max(lnorms(3), abs(var(ijk)))
        end do
        lnorms(1) = lnorms(1)/real(var_size)
        lnorms(2) = sqrt(lnorms(2)/real(var_size))


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
      Use geometry, only  : ijkmax3, dx, dy, dz, imax
      Use compar, only    : ijkstart3, ijkend3      
      Use functions, only : i_of, j_of, k_of
      Use param1, only    : zero, half
      Use usr, only       : p_g_ex, u_g_ex, v_g_ex
      Use usr, only       : xtr, ytr, ztr
      Use sendrecv, only  : send_recv
      IMPLICIT NONE

! looping indices      
        integer :: ijk, i, j, k
        integer :: ii, jj, kk

! temporary variables        
        double precision  :: xt, yt, zt

! user-defined channel problem parameters 
! (must match the input in mfix.dat)
! (better to set these up via an input file in future) 
        double precision, parameter :: mu_usr     = 0.001d0   ! viscosity
        double precision, parameter :: dpdx_usr   = 1200.d0   ! dp/dx
        double precision, parameter :: h_usr      = 0.01d0    ! height
        double precision, parameter :: l_usr      = 0.2d0     ! length
        double precision, parameter :: p0_usr     = 101325.d0 ! ref. pressure

        
! generate grid locations for exact solution calculations 
        !** would this work correctly in DMP?
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

        call send_recv(xtr)
        call send_recv(ytr)
        call send_recv(ztr)

! set exact solution
        do ijk = ijkstart3, ijkend3
          i = i_of(ijk)
          j = j_of(ijk)

! scalar variable
          xt = xtr(ijk) - dx(i)*half
          p_g_ex(ijk) = p0_usr + &
              dpdx_usr*((l_usr - half*dx(imax+1)) - xt) 
          !** dx(imax+1)??

! vector variable (x)
          yt = ytr(ijk) - dy(j)*half         
          u_g_ex(ijk) = half/mu_usr*dpdx_usr*yt*(h_usr-yt) 

! vector variable (y)          
          v_g_ex(ijk) = zero
        end do

        call send_recv(p_g_ex)
        call send_recv(u_g_ex)
        call send_recv(v_g_ex)


      RETURN

      END SUBROUTINE calculate_exact_solution_channel


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: write_tecplot_data	                               C
!  Purpose: (1) Write Tecplot data for visualization of fully          C 
!  developed Laminar flow in a horizontal channel		       C
!  (2) Write tecplot file with the L2 Norms of the discretization      C 
!  error for u_g at face center locations 			       C
!                                                                      C
!  Author: Aniruddha Choudhary                        Date: Jan 2015   C
!  email: anirudd@vt.edu					       C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date:            C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!      subroutine write_tecplot_data
!      use geometry
!      use indices
!      use fldvar
!      use functions
!      implicit none
!
!      double precision    :: three = 3.d0, two = 2.d0
!
!      integer   :: i, j, k, ijk
!      integer   :: imjk, ijmk, imjmk, imjpk, ipjmk, ipjk, ijpk
!
!      logical   :: write_solution_tecplot = .FALSE.      
!
!! x, y coordinate variables
!      double precision, dimension(2:imax+2,2:jmax+2)  :: x_tmp, y_tmp ! array
!      double precision  :: xtmp, ytmp   ! local/temporary node locations
!      double precision  :: xctmp, yctmp ! local/temporary cell-center locations
!      
!! solution variables for cell-centered data visualization	
!      double precision, dimension(2:imax+1,2:jmax+1)  :: p_g_tmp_c, &
!                                              u_g_tmp_c, v_g_tmp_c
!      double precision, dimension(2:imax+1,2:jmax+1)  :: p_g_ex_c, &
!                                              u_g_ex_c, v_g_ex_c
!
!! solution variables for node-located data visualization	
!      double precision, dimension(2:imax+2,2:jmax+2)  :: p_g_tmp_n, &
!                                              u_g_tmp_n, v_g_tmp_n
!      double precision, dimension(2:imax+2,2:jmax+2)  :: p_g_ex_n, & 
!                                              u_g_ex_n, v_g_ex_n
!
!! variables for discretization error (DE) norms calculations
!      double precision  :: p_gtmp, u_gtmp, v_gtmp
!      double precision, dimension(1:3)  :: l1de, l2de, linfde
!      ! 1 = P_g, 2 = U_g, 3 = V_g !
!      
!! user-defined channel problem parameters 
!! (must match the input in mfix.dat)
!! (better to set these up via an input file in future) 
!      double precision, parameter :: mu_usr     = 0.001d0   ! viscosity
!      double precision, parameter :: dpdx_usr   = 1200.d0   ! dp/dx
!      double precision, parameter :: h_usr      = 0.01d0    ! height
!      double precision, parameter :: l_usr      = 0.2d0     ! length
!      double precision, parameter :: p0_usr     = 101325.d0 ! ref. pressure
!     
!! create x_tmp, y_tmp array at node locations 
!      do j = 2, jmax+2
!      do i = 2, imax+2
!          x_tmp(i,j) = l_usr*real(i-2)/real(imax)
!          y_tmp(i,j) = h_usr*real(j-2)/real(jmax)
!      end do
!      end do
!       
!! calculate solution variables at cell-centers
!      k = kmin1  ! because 2D domain
!      
!      do j = 2, jmax+1
!      do i = 2, imax+1
!         ijk  = funijk(i,j,k)
!         imjk = im_of(ijk)
!         ijmk = jm_of(ijk)
!      
!         p_g_tmp_c(i,j) = p_g(ijk)
!         u_g_tmp_c(i,j) = half*(u_g(ijk)+u_g(imjk))
!         v_g_tmp_c(i,j) = half*(v_g(ijk)+v_g(ijmk))
!      
!         xctmp = half*(x_tmp(i+1,j)+x_tmp(i,j))
!         yctmp = half*(y_tmp(i,j+1)+y_tmp(i,j)) 
!
!         p_g_ex_c(i,j) = p0_usr + &
!          dpdx_usr*( (x_tmp(imax+2,j) + x_tmp(imax+1,j))*half - xctmp )
!         u_g_ex_c(i,j) = (half*mu_usr)*dpdx_usr*yctmp*(h_usr-yctmp) 
!         v_g_ex_c(i,j) = zero
!      end do
!      end do
!
!! output solution variables in Tecplot cell-centered (or BLOCK) format
!      if(write_solution_tecplot) then
!
!!*** ask jmusser: better way to find open unit numbers?
!        open(unit=777, file="solution_cellc.dat", status='unknown')
!        write(777,*) 'variables = "x""y""Pg""Ug""Vg""Pgex""Ugex""Vgex"'
!        write(777,*) 'zone T="',0,'" '
!        write(777,*) 'i=',imax+1,' j=',jmax+1
!        write(777,*) 'DATAPACKING=BLOCK'
!        write(777,*) "VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)"
!        
!        do j = 2, jmax+2
!           write(777,*) x_tmp(:,j)
!        end do
!        
!        do j = 2, jmax+2
!           write(777,*) y_tmp(:,j)
!        end do
!        
!        do j = 2, jmax+1
!           write(777,*) p_g_tmp_c(:,j)
!        end do
!        
!        do j = 2, jmax+1
!           write(777,*) u_g_tmp_c(:,j)
!        end do
!        
!        do j = 2, jmax+1
!           write(777,*) v_g_tmp_c(:,j)
!        end do
!        
!        do j = 2, jmax+1
!           write(777,*) p_g_ex_c(:,j)
!        end do
!        
!        do j = 2, jmax+1
!           write(777,*) u_g_ex_c(:,j)
!        end do
!        
!        do J = 2, jmax+1
!           write(777,*) v_g_ex_c(:,j)
!        end do
!        
!        close(777)
!        
!      end if ! end of 'if(write_solution_tecplot)'
!
!! calculate solution variables at nodes 
!      ! get pressure values at nodes
!      ! internal cells
!      do j = 3, jmax+1
!      do i = 3, imax+1
!         ijk = funijk(i,j,k)
!         imjk = im_of(ijk)
!         ijmk = jm_of(ijk)
!         imjmk = im_of(ijmk)
!         p_g_tmp_n(i,j) = &
!          (P_G(ijk)+P_G(imjk)+P_G(ijmk)+P_G(imjmk))*0.25d0
!      end do
!      end do
!      
!      ! top and bottom boundaries
!      j = 2
!      do i = 3, imax+1
!         ijk = funijk(i,j,k)
!         ijpk = jp_of(ijk)
!         imjk = im_of(ijk)
!         imjpk = jp_of(imjk)
!         p_g_tmp_n(i,j) = half*( (three*P_G(ijk) - two*P_G(ijpk)) + &
!                                (three*P_G(imjk) - two*P_G(imjpk)) )
!      end do
!      
!      j = jmax+1
!      do i = 3, imax+1
!         ijk = funijk(i,j,k)
!         ijmk = jm_of(ijk)
!         imjk = im_of(ijk)
!         imjmk = jm_of(imjk)
!         p_g_tmp_n(i,j+1) = half*( (three*P_G(ijk) - two*P_G(ijmk)) + &
!                                (three*P_G(imjk) - two*P_G(imjmk)) )
!      end do
!      
!      ! left and right boundaries
!      i = 2
!      do j = 3, jmax+1
!         ijk = funijk(i,j,k)
!         ipjk = ip_of(ijk)
!         ijmk = jm_of(ijk)
!         ipjmk = ip_of(ijmk)
!         p_g_tmp_n(i,j) = half*( (three*P_G(ijk) - two*P_G(ipjk)) + &
!                                (three*P_G(ijmk) - two*P_G(ipjmk)) )
!      end do
!      
!      i = imax+1
!      do j = 3, jmax+1
!         ijk = funijk(i,j,k)
!         imjk = im_of(ijk)
!         ijmk = jm_of(ijk)
!         imjmk = im_of(ijmk)
!         p_g_tmp_n(i+1,j) = half*( (three*P_G(ijk) - two*P_G(imjk)) + &
!                                (three*P_G(ijmk) - two*P_G(imjmk)) )
!      end do
!      
!      ! corners
!      p_g_tmp_n(2,2) = half*&
!      ( (two*p_g_tmp_n(3,2) - p_g_tmp_n(4,2)) + &
!        (two*p_g_tmp_n(2,3) - p_g_tmp_n(2,4)) )
!      p_g_tmp_n(2,jmax+2) = half*&
!      ( (two*p_g_tmp_n(3,jmax+2) - p_g_tmp_n(4,jmax+2)) + &
!        (two*p_g_tmp_n(2,jmax+1) - p_g_tmp_n(2,jmax)) )
!      p_g_tmp_n(imax+2,2) = half*&
!      ( (two*p_g_tmp_n(imax+1,2) - p_g_tmp_n(imax,2)) + &
!        (two*p_g_tmp_n(imax+2,3) - p_g_tmp_n(imax+2,4)) )
!      p_g_tmp_n(imax+2,jmax+2) = half*&
!      ( (two*p_g_tmp_n(imax+1,jmax+2) - p_g_tmp_n(imax,jmax+2)) + &
!        (two*p_g_tmp_n(imax+2,jmax+1) - p_g_tmp_n(imax+2,jmax)) )
!      
!      ! get U_g and V_g values at nodes
!      do j = 2, jmax+2
!      do i = 2, imax+2
!         ijk = funijk(i,j,k)
!         imjk = im_of(ijk)
!         ijmk = jm_of(ijk)
!         imjmk = im_of(ijmk)       
!      
!         u_g_tmp_n(i,j) = half*(U_G(imjk)+U_G(imjmk))
!         v_g_tmp_n(i,j) = half*(V_G(ijmk)+V_G(imjmk))
!      
!         xtmp = x_tmp(i,j)
!         ytmp = y_tmp(i,j)
!
!         p_g_ex_n(i,j) = p0_usr + dpdx_usr*( (x_tmp(imax+2,j) + & 
!                                  x_tmp(imax+1,j))*half - xtmp )  
!         u_g_ex_n(i,j) = one/(two*mu_usr)*dpdx_usr*ytmp*(h_usr-ytmp)   
!         v_g_ex_n(i,j) = zero
!      end do
!      end do
!
!! output solution variables in Tecplot node-located data format
!      if(write_solution_tecplot) then
!
!        open(unit = 778, file="solution_node.dat", status='unknown')
!        write(778,*) 'variables = "x""y""Pg""Ug""Vg""Pgex""Ugex""Vgex"'
!        write(778,*) 'zone T="',0,'" '
!        write(778,*) 'I=',imax+1,' J=',jmax+1
!        write(778,*) 'DATAPACKING=POINT'
!        
!        do j = 2, jmax+2
!        do i = 2, imax+2
!           write(778,*) x_tmp(i,j), y_tmp(i,j), &
!              p_g_tmp_n(i,j), u_g_tmp_n(i,j), v_g_tmp_n(i,j), &
!              p_g_ex_n(i,j), u_g_ex_n(i,j), v_g_ex_n(i,j) 
!        end do
!        end do
!
!        close(778)
!
!      end if
!
!! calculate discretization error norms
!      k = kmin1
!      
!      ! For pressure !
!      l1de(1) = zero
!      l2de(1) = zero
!      linfde(1) = zero
!      do i = 2, imax+1
!      do j = 2, jmax+1
!         l1de(1) = l1de(1) + abs(p_g_tmp_c(i,j) - p_g_ex_c(i,j))
!         l2de(1) = l2de(1) + abs(p_g_tmp_c(i,j) - p_g_ex_c(i,j))**2
!         linfde(1) = max(linfde(1), abs(p_g_tmp_c(i,j) - p_g_ex_c(i,j)))
!      end do
!      end do
!      l1de(1) = l1de(1)/real(imax*jmax)
!      l2de(1) = sqrt(l2de(1)/real(imax*jmax))
!      
!      ! For U_g !
!      l1de(2) = zero
!      l2de(2) = zero
!      linfde(2) = zero
!      do i = 2, imax+2
!      do j = 2, jmax+1
!         yctmp = half*(y_tmp(i,j) + y_tmp(i,j+1))
!         u_gtmp = one/(two*mu_usr)*dpdx_usr*yctmp*(h_usr-yctmp)   
!      
!         ijk = funijk(i,j,k)
!         imjk = im_of(ijk)
!         l1de(2) = l1de(2) + abs(U_G(imjk) - u_gtmp)
!         l2de(2) = l2de(2) + abs(U_G(imjk) - u_gtmp)**2   
!         linfde(2) = max(linfde(2), abs(U_G(imjk) - u_gtmp))
!      end do
!      end do
!      l1de(2) = l1de(2)/real((imax+1)*jmax)
!      l2de(2) = sqrt(l2de(2)/real((imax+1)*jmax))
!      
!      ! For V_g !
!      l1de(3) = zero
!      l2de(3) = zero
!      linfde(3) = zero
!      do i = 2, imax+1
!      do j = 2, jmax+2
!         v_gtmp = zero   
!      
!         ijk = funijk(i,j,k)
!         ijmk = jm_of(ijk)
!         l1de(3) = l1de(3) + abs(V_G(ijmk) - v_gtmp)
!         l2de(3) = l2de(3) + abs(V_G(ijmk) - v_gtmp)**2   
!         linfde(3) = max(linfde(3), abs(V_G(ijmk) - v_gtmp))
!      end do
!      end do
!      l1de(3) = l1de(3)/real(imax*(jmax+1))
!      l2de(3) = sqrt(l2de(3)/real(imax*(jmax+1)))
!
!! output discretization error norms in tecplot format      
!      open(unit = 779, file="de_norms.dat", status='unknown')
!      write(779,*) "DE Norms for Horizontal Channel Flow:"
!      write(779,*) "imax= ",imax, " jmax=", jmax
!      write(779,*) "1st line: L1 Norms, 2nd line: L2 Norms, &
!                   &3rd line: Linf Norms"
!      write(779,*) "Columns: P_g : U_g : V_g"
!      write(779,*) l1de(1), l1de(2), l1de(3)
!      write(779,*) l2de(1), l2de(2), l2de(3)
!      write(779,*) linfde(1), linfde(2), linfde(3)
!      close(779)
!      
!      end subroutine write_tecplot_data
!!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
