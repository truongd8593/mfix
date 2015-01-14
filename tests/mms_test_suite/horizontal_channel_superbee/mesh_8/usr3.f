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
      USE usr
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

        CALL write_tecplot_data

      RETURN  
      END SUBROUTINE USR3 

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
      subroutine write_tecplot_data
      !USE param
      !USE param1
      !USE parallel
      !USE constant
      !USE run
      !USE toleranc
      use geometry
      use indices
      !USE compar
      !USE sendrecv
      use fldvar
      use functions
      implicit none

      integer, parameter  :: dp = kind(1.d0)
!**** ask jmusser:      
      ! or
      ! integer, parameter  :: dp = selected_real_kind(15,307)
      ! this is guaranteed to have double precision. Variables can be declared
      ! as:
      ! real(dp) = a, b, etc.

      integer  :: i, j, k, ijk
      integer  :: imjk, ijmk, imjmk, imjpk, ipjmk, ipjk, ijkp
      
!**** ask jmusser: should these multi-dim. arrays be made allocatable?

! x, y coordinate variables
      double precision, dimension(2:imax+2,2:jmax+2)  :: x_tmp, y_tmp ! array
      double precision  :: xtmp, ytmp   ! local/temporary node locations
      double precision  :: xctmp, yctmp ! local/temporary cell-center locations
      
! solution variables for cell-centered data visualization	
      double precision, dimension(2:imax+1,2:jmax+1)  :: p_g_tmp_c, u_g_tmp_c, &
                                                         v_g_tmp_c
      double precision, dimension(2:imax+1,2:jmax+1)  :: p_g_ex_c, u_g_ex_c, &
                                                         v_g_ex_c

! solution variables for node-located data visualization	
      double precision, dimension(2:imax+2,2:jmax+2)  :: p_g_tmp_n, u_g_tmp_n, &
                                                         v_g_tmp_n
      double precision, dimension(2:imax+2,2:jmax+2)  :: p_g_ex_n, u_g_ex_n, &
                                                         v_g_ex_n

! variables for discretization error (DE) norms calculations
      double precision  :: p_gtmp, u_gtmp, v_gtmp
      double precision, dimension(1:3)  :: l1de, l2de, linfde
      ! 1 = P_g, 2 = U_g, 3 = V_g !
      
! user-defined channel problem parameters 
! (must match the input in mfix.dat)
! (better to set these up via an input file in future) 
      double precision, parameter :: mu_usr     = 0.001   ! viscosity
      double precision, parameter :: dpdx_usr   = 1200.0  ! dp/dx
      double precision, parameter :: h_usr      = 0.01    ! height
      double precision, parameter :: l_usr      = 0.2     ! length
      double precision, parameter :: p0_usr     = 101325.0  ! ref. pressure
     
! create x_tmp, y_tmp array at node locations 

      do i = 2, imax+2
      do j = 2, jmax+2
          x_tmp(i,j) = l_usr*real(i-2,dp)/real(imax)
          y_tmp(i,j) = h_usr*real(j-2,dp)/real(jmax)
      end do
      end do
       
! calculate solution variables at cell-centers

      k = kmin1  ! because 2D domain
      
      do i = 2, imax+1
      do j = 2, jmax+1
         ijk  = funijk(i,j,k)
         imjk = im_of(ijk)
         ijmk = jm_of(ijk)
      
         p_g_tmp_c(i,j) = p_g(ijk)
         u_g_tmp_c(i,j) = half*(u_g(ijk)+u_h(imjk))
         v_g_tmp_c(i,j) = half*(v_g(ijk)+v_g(ijmk))
      
         xctmp = half*(x_tmp(i+1,j)+x_tmp(i,j))
         yctmp = half*(y_tmp(i,j+1)+y_tmp(i,j)) 

         p_g_ex_c(i,j) = p0_usr + dpdx_usr*( (x_tmp(imax+2,j) + x_tmp(imax+1,j))*half - xctmp )   
         u_g_ex_c(i,j) = one/(two*mu_usr)*dpdx_usr*yctmp*(h_usr-yctmp)   
         v_g_ex_c(i,j) = zero
      end do
      end do

! output solution variables in Tecplot cell-centered (or BLOCK) format      
!*** ask jmusser: better way to find open unit numbers?
      open(unit=777, file="solution_cellc.dat", status='unknown')
      write(777,*) 'variables = "x""y""Pg""Ug""Vg""Pgex""Ugex""Vgex"'
      write(777,*) 'zone T="',0,'" '
      write(777,*) 'i=',imax+1,' j=',jmax+1
      write(777,*) 'DATAPACKING=BLOCK'
      write(777,*) "VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)"
      
      do j = 2, jmax+2
         write(777,*) x_tmp(:,j)
      end do
      
      do j = 2, jmax+2
         write(777,*) y_tmp(:,j)
      end do
      
      do j = 2, jmax+1
         write(777,*) p_g_tmp_c(:,j)
      end do
      
      do j = 2, jmax+1
         write(777,*) u_g_tmp_c(:,j)
      end do
      
      do j = 2, jmax+1
         write(777,*) b_g_tmp_c(:,j)
      end do
      
      do j = 2, jmax+1
         write(777,*) p_g_ex_c(:,j)
      end do
      
      do j = 2, jmax+1
         write(777,*) u_g_ex_c(:,j)
      end do
      
      do J = 2, jmax+1
         write(777,*) v_g_ex_c(:,j)
      end do
      
      close(777)
      
! calculate solution variables at nodes 

      ! get pressure values at nodes
      ! internal cells
      do i = 3, imax+1
      do j = 3, jmax+1
         ijk = funijk(i,j,k)
         imjk = im_of(ijk)
         ijmk = jm_of(ijk)
         imjmk = im_of(ijmk)
         p_g_tmp_n(i,j) = (P_G(ijk)+P_G(imjk)+P_G(ijmk)+P_G(imjmk))*0.25_dp
         !*** ask jmusser: does simply 0.25 result in guaranteed double precision?
         ! better is 0.25_dp
      end do
      end do
      
      ! top and bottom boundaries
      j = 2
      do i = 3, imax+1
         ijk = funijk(i,j,k)
         ijpk = jp_of(ijk)
         imjk = im_of(ijk)
         imjpk = jp_of(imjk)
         p_g_tmp_n(i,j) = half*( (three*P_G(ijk) - two*P_G(ijpk)) + &
                                (three*P_G(imjk) - two*P_G(imjpk)) )
      end do
      
      j = jmax+1
      do i = 3, imax+1
         ijk = funijk(i,j,k)
         ijmk = jm_of(ijk)
         imjk = im_of(ijk)
         imjmk = jm_of(imjk)
         p_g_tmp_n(i,j+1) = half*( (three*P_G(ijk) - two*P_G(ijmk)) + &
                                (three*P_G(imjk) - two*P_G(imjmk)) )
      end do
      
      ! left and right boundaries
      i = 2
      do j = 3, jmax+1
         ijk = funijk(i,j,k)
         ipjk = ip_of(ijk)
         ijmk = jm_of(ijk)
         ipjmk = ip_of(ijmk)
         p_g_tmp_n(i,j) = half*( (three*P_G(ijk) - two*P_G(ipjk)) + &
                                (three*P_G(ijmk) - two*P_G(ipjmk)) )
      end do
      
      i = imax+1
      do j = 3, jmax+1
         ijk = funijk(i,j,k)
         imjk = im_of(ijk)
         ijmk = jm_of(ijk)
         imjmk = im_of(ijmk)
         p_g_tmp_n(i+1,j) = half*( (three*P_G(ijk) - two*P_G(imjk)) + &
                                (three*P_G(ijmk) - two*P_G(imjmk)) )
      end do
      
      ! corners
      p_g_tmp_n(2,2) = half*&
      ( (two*p_g_tmp_n(3,2) - p_g_tmp_n(4,2)) + &
        (two*p_g_tmp_n(2,3) - p_g_tmp_n(2,4)) )
      p_g_tmp_n(2,jmax+2) = half*&
      ( (two*p_g_tmp_n(3,jmax+2) - p_g_tmp_n(4,jmax+2)) + &
        (two*p_g_tmp_n(2,jmax+1) - p_g_tmp_n(2,jmax)) )
      p_g_tmp_n(imax+2,2) = half*&
      ( (two*p_g_tmp_n(imax+1,2) - p_g_tmp_n(imax,2)) + &
        (two*p_g_tmp_n(imax+2,3) - p_g_tmp_n(imax+2,4)) )
      p_g_tmp_n(imax+2,jmax+2) = half*&
      ( (two*p_g_tmp_n(imax+1,jmax+2) - p_g_tmp_n(imax,jmax+2)) + &
        (two*p_g_tmp_n(imax+2,jmax+1) - p_g_tmp_n(imax+2,jmax)) )
      
      ! get U_g and V_g values at nodes
      do i = 2, imax+2
      do j = 2, jmax+2
         ijk = funijk(i,j,k)
         imjk = im_of(ijk)
         ijmk = jm_of(ijk)
         imjmk = im_of(ijmk)       
      
         u_g_tmp_n(i,j) = half*(U_G(imjk)+U_G(imjmk))
         v_g_tmp_n(i,j) = half*(V_G(ijmk)+V_G(imjmk))
      
         xtmp = x_tmp(i,j)
         ytmp = y_tmp(i,j)

         p_g_ex_n(i,j) = p0_usr + dpdx_usr*( (x_tmp(imax+2,j) + & 
                                  x_tmp(imax+1,j))*half - xtmp )  
         u_g_ex_n(i,j) = one/(two*mu_usr)*dpdx_usr*ytmp*(h_usr-ytmp)   
         v_g_ex_n(i,j) = zero
      end do
      end do

! output solution variables in Tecplot node-located data format
      open(unit = 778, file="solution_node.dat", status='unknown')
      write(778,*) 'variables = "x""y""Pg""Ug""Vg""Pgex""Ugex""Vgex"'
      write(778,*) 'zone T="',0,'" '
      write(778,*) 'I=',imax+1,' J=',jmax+1
      write(778,*) 'DATAPACKING=POINT'
      
      do J = 2, jmax+2
      do I = 2, imax+2
         write(778,*) x_tmp(i,j), y_tmp(i,j), p_g_tmp_n(i,j), u_g_tmp_n(i,j), &
            v_g_tmp_n(i,j), p_g_ex_n(i,j), u_g_ex_n(i,j), v_g_ex_n(i,j) 
      end do
      end do
      
      ! calculate and output DE norms
      k = kmin1
      
      ! For pressure !
      l1de(1) = zero
      l2de(1) = zero
      linfde(1) = zero
      do i = 2, imax+1
      do j = 2, jmax+1
         l1de(1) = l1de(1) + abs(p_g_tmp_c(i,j) - p_g_ex_c(i,j))
         l2de(1) = l2de(1) + abs(p_g_tmp_c(i,j) - p_g_ex_c(i,j))**2
         linfde(1) = max(linfde(1), abs(p_g_tmp_c(i,j) - p_g_ex_c(i,j)))
      end do
      end do
      l1de(1) = l1de(1)/real(imax*jmax)
      l2de(1) = sqrt(l2de(1)/real(imax*jmax))
      
      ! For U_g !
      l1de(2) = zero
      l2de(2) = zero
      linfde(2) = zero
      do i = 2, imax+2
      do j = 2, jmax+1
         yctmp = half*(y_tmp(i,j) + y_tmp(i,j+1))
         u_gtmp = one/(two*mu_usr)*dpdx_usr*yctmp*(h_usr-yctmp)   
      
         ijk = funijk(i,j,k)
         imjk = im_of(ijk)
         l1de(2) = l1de(2) + abs(U_G(imjk) - u_g_tmp)
         l2de(2) = l2de(2) + abs(U_G(imjk) - u_g_tmp)**2   
         linfde(2) = max(linfde(2), abs(U_G(imjk) - U_g_tmp))
      end do
      end do
      l1de(2) = l1de(2)/real((imax+1)*jmax,dp)
      l2de(2) = sqrt(l2de(2)/real((imax+1)*jmax),dp)
      
      ! For V_g !
      l1de(3) = zero
      l2de(3) = zero
      linfde(3) = zero
      do i = 2, imax+1
      do j = 2, jmax+2
         v_gtmp = zero   
      
         ijk = funijk(i,j,k)
         ijmk = jm_of(ijk)
         l1de(3) = l1de(3) + abs(V_G(ijmk) - v_g_tmp)
         l2de(3) = l2de(3) + abs(V_G(ijmk) - v_g_tmp)**2   
         linfde(3) = max(linfde(3), abs(V_G(ijmk) - v_g_tmp))
      end do
      end do
      l1de(3) = l1de(3)/real(imax*(jmax+1),2)
      l2de(3) = sqrt(l2de(3)/real(imax*(jmax+1)),2)
      
      open(UNIT = 779, FILE="de_norms.dat", Status='unknown')
      write(779,*) "DE Norms for Horizontal Channel Flow:"
      write(779,*) "imax= ",imax, " jmax=", jmax
      write(779,*) "1st line: L1 Norms, 2nd line: L2 Norms, 3rd line: Linf Norms"
      write(779,*) "Columns: P_g : U_g : V_g"
      write(779,*) l1de(1), l1de(2), l1de(3)
      write(779,*) l2de(1), l2de(2), l2de(3)
      write(779,*) linfde(1), linfde(2), linfde(3)
      close(779)
      
      end subroutine write_tecplot_data
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
