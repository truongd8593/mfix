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

        CALL WRITE_TECPLOT_DATA

      RETURN  
      END SUBROUTINE USR3 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_TECPLOT_DATA	                               C
!  Purpose: (1) Write Tecplot data for visualization of fully developedC
!  Laminar flow in a horizontal channel				       C
!  (2) Write tecplot file containing the L2 Norms of the dicretization C 
!  error for U_g at face center locations 			       C
!                                                                      C
!  Author: Aniruddha Choudhary                        Date: 21-Mar-11  C
!  email: anirudd@vt.edu					       C					
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
       SUBROUTINE WRITE_TECPLOT_DATA
       USE param
       USE param1
       USE parallel
       USE constant
       USE run
       USE toleranc
       USE geometry
       USE indices
       USE compar
       USE sendrecv
       USE fldvar
     
       IMPLICIT NONE
       INTEGER :: I,J,K,IJK
       INTEGER :: IMJK, IJMK, IMJMK, IMJPK, IPJMK, IPJK, IJPK
       
       !! X, Y coordinate variables !!
       REAL*8, DIMENSION(2:IMAX+2,2:JMAX+2) :: X_temp, Y_temp ! X, Y node locations
       REAL*8 :: xctmp, yctmp, xtmp, ytmp ! temporary variables holding x, y locations		
       
       !! Variables for visualization !!	
       ! For cell-centered data visualization !
       REAL*8, DIMENSION(2:IMAX+1,2:JMAX+1) :: P_g_temp_c, U_g_temp_c, V_g_temp_c
       REAL*8, DIMENSION(2:IMAX+1,2:JMAX+1) :: P_g_exact_c, U_g_exact_c, V_g_exact_c
       
       ! For node-located data visualization !
       REAL*8, DIMENSION(2:IMAX+2,2:JMAX+2) :: P_g_temp_n, U_g_temp_n, V_g_temp_n
       REAL*8, DIMENSION(2:IMAX+2,2:JMAX+2) :: P_g_exact_n, U_g_exact_n, V_g_exact_n
       
       !! Variables for DE Norms calculations !!
       REAL*8 :: P_g_temp, U_g_temp, V_g_temp
       REAL*8, DIMENSION(1:3) :: L1DE, L2DE, LinfDE
       ! 1 = P_g, 2 = U_g, 3 = V_g !
       
       !! Channel Problem parameters - (better to set these up via an input file) !!
       DOUBLE PRECISION, PARAMETER :: mu_usr=0.001d0, dpdy_usr=1200.0d0, h_usr=0.2d0, L_usr=0.01d0
       DOUBLE PRECISION, PARAMETER :: P0_usr=101325.0d0
       
       !! Double precision constants !!
       DOUBLE PRECISION, PARAMETER :: FOURTH=0.25d0, TWO=2.0d0, THREE=3.0d0
       
       include "function.inc"
     
       ! Create X_temp, Y_temp at node locations ! 
       DO I = 2, IMAX+2
       DO J = 2, JMAX+2
          !X_temp(I,J) = real(i-2)*dx(2)	! dx of first interior cell
          !Y_temp(I,J) = real(j-2)*dy(2)	! dy of     ---"---
          X_temp(I,J) = L_usr*real(i-2)/real(IMAX)
          Y_temp(I,J) = h_usr*real(j-2)/real(JMAX)
       END DO
       END DO
       
       ! Calculate and output cell-centered data visulatization variables !
       K = KMIN1 ! 2-D domain only
       
       DO I = 2, IMAX+1
       DO J = 2, JMAX+1
          ijk = funijk(i,j,k)
          imjk = im_of(ijk)
          ijmk = jm_of(ijk)
       
          P_g_temp_c(i,j) = P_G(IJK)
          U_g_temp_c(i,j) = HALF*(U_G(ijk)+U_G(imjk))
          V_g_temp_c(i,j) = HALF*(V_G(ijk)+V_G(ijmk))
       
          xctmp = HALF*(X_temp(i+1,j)+X_temp(i,j))
          yctmp = HALF*(Y_temp(i,j+1)+Y_temp(i,j)) 

          P_g_exact_c(i,j) = P0_usr + dpdy_usr*( HALF*(Y_temp(i,JMAX+2)+Y_temp(i,JMAX+1)) - yctmp)   
          U_g_exact_c(i,j) = ZERO          
          V_g_exact_c(i,j) = ONE/(TWO*mu_usr)*dpdy_usr*xctmp*(L_usr-xctmp)   
       END DO
       END DO
       
       OPEN(UNIT = 777, FILE="solution_cellc.dat", Status='unknown')
       Write(777,*) 'variables = "x""y""Pg""Ug""Vg""Pgex""Ugex""Vgex"'
       Write(777,*) 'zone T="',0,'" '
       Write(777,*) 'I=',IMAX+1,' J=',JMAX+1
       Write(777,*) 'DATAPACKING=BLOCK'
       WRITE(777,*) "VARLOCATION=([3,4,5,6,7,8]=CELLCENTERED)"
       
       DO J =2, JMAX+2
          Write(777,*) X_temp(:,j)
       END DO
       
       DO J =2, JMAX+2
          Write(777,*) Y_temp(:,j)
       END DO
       
       DO J =2, JMAX+1
          Write(777,*) P_g_temp_c(:,j)
       END DO
       
       DO J =2, JMAX+1
          Write(777,*) U_g_temp_c(:,j)
       END DO
       
       DO J =2, JMAX+1
          Write(777,*) V_g_temp_c(:,j)
       END DO
       
       DO J =2, JMAX+1
          Write(777,*) P_g_exact_c(:,j)
       END DO
       
       DO J =2, JMAX+1
          Write(777,*) U_g_exact_c(:,j)
       END DO
       
       DO J =2, JMAX+1
          Write(777,*) V_g_exact_c(:,j)
       END DO
       
       CLOSE(777)
       
       
!!       !	!! For Debugging only !!
!!       OPEN(UNIT = 780, FILE="test_out.dat", Status='unknown')
!!       !	DO I = 2, IMAX+1
!!       !	DO J = 2, JMAX+1
!!       i=IMAX+1
!!       j=JMAX+1
!!          ijk = funijk(i,j,k)
!!          Write(780,*) "i= ", i, "j= ", j, "ijk= ", ijk
!!          Write(780,*) "X_temp = ", "Y_temp=", X_temp(i,j), Y_temp(i,j)
!!          Write(780,*) "P_g(ijk) = ", P_g(ijk), "Pgex(ijk) = ", P_g_exact_c(i,j)
!!          Write(780,*) "Xclast = ", (X_temp(IMAX+1,j) + X_temp(IMAX,j))*HALF
!!       
!!       !	END DO
!!       !	END DO
!!       !	OPEN(UNIT = 780, FILE="test_out.dat", Status='unknown')
!!       !	DO I = 2, IMAX+2
!!       !	DO J = 2, JMAX+2
!!       !	   ijk = funijk(i,j,k)
!!       !	   Write(780,*) "i= ", i, "j= ", j, "P_g= ", P_G(IJK)
!!       !	END DO
!!       !	END DO
!!       !
!!       !	DO I = 2, IMAX+2
!!       !	DO J = 1, JMAX+2
!!       !	   ijk = funijk(i,j,k)
!!       !	   imjk = im_of(ijk)
!!       !	   Write(780,*) "i= ", i, "j= ", j, "U_g= ", U_G(IMJK)   
!!       !	END DO
!!       !	END DO
!!       !
!!       !	DO I = 1, IMAX+2
!!       !	DO J = 2, JMAX+2
!!       !	   ijk = funijk(i,j,k)
!!       !	   ijmk = jm_of(ijk)
!!       !	   Write(780,*) "i= ", i, "j= ", j, "V_g= ", V_G(IJMK)   
!!       !	END DO
!!       !	END DO
!!       !
!!       CLOSE(780)
       
       ! Get pressure values at nodes !
       ! Internal cells !
       DO I = 3, IMAX+1
       DO J = 3, JMAX+1
          ijk = funijk(i,j,k)
          imjk = im_of(ijk)
          ijmk = jm_of(ijk)
          imjmk = im_of(ijmk)
          P_g_temp_n(i,j) = (P_G(ijk)+P_G(imjk)+P_G(ijmk)+P_G(imjmk))*FOURTH
       END DO
       END DO
       
       ! Top and Bottom Boundaries !
       J = 2
       DO I = 3, IMAX+1
          ijk = funijk(i,j,k)
          ijpk = jp_of(ijk)
          imjk = im_of(ijk)
          imjpk = jp_of(imjk)
          P_g_temp_n(i,j) = HALF*( (THREE*P_G(ijk) - TWO*P_G(ijpk)) + &
          (THREE*P_G(imjk) - TWO*P_G(imjpk)) )
       END DO
       
       J = JMAX + 1
       DO I = 3, IMAX+1
          ijk = funijk(i,j,k)
          ijmk = jm_of(ijk)
          imjk = im_of(ijk)
          imjmk = jm_of(imjk)
          P_g_temp_n(i,j+1) = HALF*( (THREE*P_G(ijk) - TWO*P_G(ijmk)) + &
          (THREE*P_G(imjk) - TWO*P_G(imjmk)) )
       END DO
       
       ! Left and Right Boundaries !
       I = 2
       DO J = 3, JMAX+1
          ijk = funijk(i,j,k)
          ipjk = ip_of(ijk)
          ijmk = jm_of(ijk)
          ipjmk = ip_of(ijmk)
          P_g_temp_n(i,j) = HALF*( (THREE*P_G(ijk) - TWO*P_G(ipjk)) + &
          (THREE*P_G(ijmk) - TWO*P_G(ipjmk)) )
       END DO
       
       I = IMAX+1
       DO J = 3, JMAX+1
          ijk = funijk(i,j,k)
          imjk = im_of(ijk)
          ijmk = jm_of(ijk)
          imjmk = im_of(ijmk)
          P_g_temp_n(i+1,j) = HALF*( (THREE*P_G(ijk) - TWO*P_G(imjk)) + &
          (THREE*P_G(ijmk) - TWO*P_G(imjmk)) )
       END DO
       
       ! Corners !
       P_g_temp_n(2,2) = HALF*( (TWO*P_g_temp_n(3,2) - P_g_temp_n(4,2)) + &
         (TWO*P_g_temp_n(2,3) - P_g_temp_n(2,4)) )
       P_g_temp_n(2,JMAX+2) = HALF*( (TWO*P_g_temp_n(3,JMAX+2) - P_g_temp_n(4,JMAX+2)) + &
         (TWO*P_g_temp_n(2,JMAX+1) - P_g_temp_n(2,JMAX)) )
       P_g_temp_n(IMAX+2,2) = HALF*( (TWO*P_g_temp_n(IMAX+1,2) - P_g_temp_n(IMAX,2)) + &
         (TWO*P_g_temp_n(IMAX+2,3) - P_g_temp_n(IMAX+2,4)) )
       P_g_temp_n(IMAX+2,JMAX+2) = HALF*( (TWO*P_g_temp_n(IMAX+1,JMAX+2) - P_g_temp_n(IMAX,JMAX+2)) + &
         (TWO*P_g_temp_n(IMAX+2,JMAX+1) - P_g_temp_n(IMAX+2,JMAX)) )
       
       ! U_g and V_g calculations !
       DO I = 2, IMAX+2
       DO J = 2, JMAX+2
          ijk = funijk(i,j,k)
          imjk = im_of(ijk)
          ijmk = jm_of(ijk)
          imjmk = im_of(ijmk)       
       
          U_g_temp_n(i,j) = HALF*(U_G(imjk)+U_G(imjmk))
          V_g_temp_n(i,j) = HALF*(V_G(ijmk)+V_G(imjmk))
       
          xtmp = X_temp(i,j)
          ytmp = Y_temp(i,j)

          P_g_exact_n(i,j) = P0_usr + dpdy_usr*( HALF*(Y_temp(i,JMAX+2)+Y_temp(i,JMAX+1)) - ytmp)   
          U_g_exact_n(i,j) = ZERO          
          V_g_exact_n(i,j) = ONE/(TWO*mu_usr)*dpdy_usr*xtmp*(L_usr-xtmp)

       END DO
       END DO
       
       OPEN(UNIT = 778, FILE="solution_node.dat", Status='unknown')
       Write(778,*) 'variables = "x""y""Pg""Ug""Vg""Pgex""Ugex""Vgex"'
       Write(778,*) 'zone T="',0,'" '
       Write(778,*) 'I=',IMAX+1,' J=',JMAX+1
       Write(778,*) 'DATAPACKING=POINT'
       
       DO J = 2, JMAX+2
       DO I = 2, IMAX+2
          Write(778,*) X_temp(i,j), Y_temp(i,j), P_g_temp_n(i,j), U_g_temp_n(i,j), &
             V_g_temp_n(i,j), P_g_exact_n(i,j), U_g_exact_n(i,j), V_g_exact_n(i,j) 
       END DO
       END DO
       
       !! Calculate and output DE Norms !!
       K = KMIN1
       
       ! For pressure !
       L1DE(1) = ZERO
       L2DE(1) = ZERO
       LinfDE(1) = ZERO
       DO I = 2, IMAX+1
       DO J = 2, JMAX+1
          L1DE(1) = L1DE(1) + abs(P_g_temp_c(i,j) - P_g_exact_c(i,j))
          L2DE(1) = L2DE(1) + abs(P_g_temp_c(i,j) - P_g_exact_c(i,j))**2
          LinfDE(1) = MAX(LinfDE(1), abs(P_g_temp_c(i,j) - P_g_exact_c(i,j)))
       END DO
       END DO
       L1DE(1) = L1DE(1)/REAL(IMAX*JMAX)
       L2DE(1) = SQRT(L2DE(1)/REAL(IMAX*JMAX))
       
       ! For U_g !
       L1DE(2) = ZERO
       L2DE(2) = ZERO
       LinfDE(2) = ZERO
       DO I = 2, IMAX+2
       DO J = 2, JMAX+1
          U_g_temp = ZERO
       
          IJK = funijk(i,j,k)
          imjk = im_of(ijk)
          L1DE(2) = L1DE(2) + abs(U_G(imjk) - U_g_temp)
          L2DE(2) = L2DE(2) + abs(U_G(imjk) - U_g_temp)**2   
          LinfDE(2) = MAX(LinfDE(2), abs(U_G(imjk) - U_g_temp))
       END DO
       END DO
       L1DE(2) = L1DE(2)/REAL((IMAX+1)*JMAX)
       L2DE(2) = SQRT(L2DE(2)/REAL((IMAX+1)*JMAX))
       
       ! For V_g !
       L1DE(3) = ZERO
       L2DE(3) = ZERO
       LinfDE(3) = ZERO
       DO I = 2, IMAX+1
       DO J = 2, JMAX+2
          xctmp = ( X_temp(i,j)+X_temp(i+1,j) )*HALF
          V_g_temp = ONE/(TWO*mu_usr)*dpdy_usr*xctmp*(L_usr-xctmp)   
       
          IJK = funijk(i,j,k)
          ijmk = jm_of(ijk)
          L1DE(3) = L1DE(3) + abs(V_G(ijmk) - V_g_temp)
          L2DE(3) = L2DE(3) + abs(V_G(ijmk) - V_g_temp)**2   
          LinfDE(3) = MAX(LinfDE(3), abs(V_G(ijmk) - V_g_temp))
       END DO
       END DO
       L1DE(3) = L1DE(3)/REAL(IMAX*(JMAX+1))
       L2DE(3) = SQRT(L2DE(3)/REAL(IMAX*(JMAX+1)))
       
       OPEN(UNIT = 779, FILE="de_norms.dat", Status='unknown')
       Write(779,*) "DE Norms for Horizontal Channel Flow:"
       Write(779,*) "IMAX= ",IMAX, " JMAX=", JMAX
       Write(779,*) "1st line: L1 Norms, 2nd line: L2 Norms, 3rd line: Linf Norms"
       Write(779,*) "Columns: P_g : U_g : V_g"
       Write(779,*) L1DE(1), L1DE(2), L1DE(3)
       Write(779,*) L2DE(1), L2DE(2), L2DE(3)
       Write(779,*) LinfDE(1), LinfDE(2), LinfDE(3)
       CLOSE(779)
       
    END SUBROUTINE WRITE_TECPLOT_DATA
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
