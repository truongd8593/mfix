!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: SOLID_DRAG                                              C
!  Purpose: Accounting for the equal and opposite drag force on a      C
!  continuous solids phase due to discrete solid particles by          C
!  introducing the drag as a source term. Face centered.               C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
     
      SUBROUTINE SOLID_DRAG(A_M, B_M, VELDIR, IER)

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE matrix 
      USE geometry
      USE physprop
      USE indices
      USE compar    
      USE discretelement 
      
      IMPLICIT NONE
!-----------------------------------------------
! Dummy Arguments
!-----------------------------------------------
! Septadiagonal matrix A_m       
      DOUBLE PRECISION :: A_M(DIMENSION_3, -3:3, 0:DIMENSION_M)
! Vector b_m       
      DOUBLE PRECISION :: B_M(DIMENSION_3, 0:DIMENSION_M)
! flag/index used to specify which gas phase momentum equation
! (x,y,z) is currently being solved 
      INTEGER :: VELDIR
! Error index
      INTEGER :: IER
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! Indices
      INTEGER :: IJK, I, J, K
! Solids phase indices
      INTEGER :: CM, DM
! Face center values of u_sm (i+1/2), v_sm (j+1/2) and w_sm (k+1/2)
      DOUBLE PRECISION :: USFCM, VSFCM, WSFCM
! Volume of x, y, or z cell on staggered grid      
      DOUBLE PRECISION :: VCELL
! temporary variables for matrix A_M and vector B_M that are 
! used for local calculations
      DOUBLE PRECISION tmp_A, tmp_B
! Averaging factor 
! (=0.25 in 3D and =0.5 in 2D)      
      DOUBLE PRECISION :: AVG_FACTOR
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE '../function.inc'
      INCLUDE '../fun_avg1.inc'
      INCLUDE '../fun_avg2.inc'
!-----------------------------------------------      

      AVG_FACTOR = merge(0.5d0, 0.25D0, DO_K)

      IF(VELDIR.EQ.1) THEN
         DO CM = 1, MMAX
            DO DM = 1, DES_MMAX
               DO IJK = IJKSTART3, IJKEND3
                  IF(FLUID_AT(IJK)) THEN
                     I = I_OF(IJK)
                     J = J_OF(IJK)
                     K = K_OF(IJK)

! currently no difference between interpolated and non-interpolated
! implementation of solid-solid drag                  
                     USFCM = AVG_X(DES_U_S(IJK,DM),DES_U_S(EAST_OF(IJK),DM),I)
                     tmp_A =  - VXF_SDS(IJK,CM,DM)
                     tmp_B =  - VXF_SDS(IJK,CM,DM)*USFCM

                     A_M(IJK,0,CM) = A_M(IJK,0,CM) + tmp_A
                     B_M(IJK,CM) = B_M(IJK,CM) + tmp_B

                  ENDIF   ! end if (fluid_at(ijk))
               ENDDO   ! end do (ijk=ijkstart3,ijkend3)
            ENDDO   ! end do (dm=1,des_mmax)
         ENDDO   ! end do (cm=1,mmax)

      ELSEIF(VELDIR.EQ.2) THEN
         DO CM = 1, MMAX
            DO DM = 1, DES_MMAX
               DO IJK = IJKSTART3, IJKEND3
                  IF(FLUID_AT(IJK)) THEN
                     I = I_OF(IJK)
                     J = J_OF(IJK)
                     K = K_OF(IJK)

! currently no difference between interpolated and non-interpolated
! implementation of solid-solid drag
                     VSFCM = AVG_Y(DES_V_S(IJK,DM),DES_V_S(NORTH_OF(IJK),DM),J)
                     tmp_A =  - VXF_SDS(IJK,CM,DM)
                     tmp_B =  - VXF_SDS(IJK,CM,DM)*VSFCM

                     A_M(IJK,0,CM) = A_M(IJK,0,CM) + tmp_A
                     B_M(IJK,CM) = B_M(IJK,CM) + tmp_B

                  ENDIF   ! end if (fluid_at(ijk))
               ENDDO   ! end do (ijk=ijkstart3,ijkend3)
            ENDDO   ! end do (dm=1,des_mmax)
         ENDDO   ! end do (cm=1,mmax)

      ELSEIF(VELDIR.EQ.3) THEN
         DO CM = 1, MMAX
            DO DM = 1, DES_MMAX
               DO IJK = IJKSTART3, IJKEND3
                  IF(FLUID_AT(IJK)) THEN
                     I = I_OF(IJK)
                     J = J_OF(IJK)
                     K = K_OF(IJK)

! currently no difference between interpolated and non-interpolated
! implementation of solid-solid drag                  
                     WSFCM = AVG_Z(DES_W_S(IJK,DM),DES_W_S(TOP_OF(IJK),DM),K)
                     tmp_A =  - VXF_SDS(IJK,CM,DM)
                     tmp_B =  - VXF_SDS(IJK,CM,DM)*WSFCM

                     A_M(IJK,0,CM) = A_M(IJK,0,CM) + tmp_A
                     B_M(IJK,CM) = B_M(IJK,CM) + tmp_B

                  ENDIF   ! end if (fluid_at(ijk))
               ENDDO   ! end do (ijk=ijkstart3,ijkend3)
            ENDDO   ! end do (dm=1,des_mmax)
         ENDDO   ! end do (cm=1,mmax)

      ENDIF

      RETURN
      END SUBROUTINE SOLID_DRAG
