!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  GET_DEL_H                                             C
!  Purpose: Finds the normal distance and unit vector from any point   C
!  (x0,y0,z0) to a Cut face                                            C
!  This subroutine must be called from a cut-cell                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_DEL_H(IJK,TYPE_OF_CELL,X0,Y0,Z0,Del_H,Nx,Ny,Nz)
    
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
      USE quadric
      USE cutcell
      
      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: X0,Y0,Z0,XREF,YREF,ZREF
      INTEGER :: IJK,I,J,K
      DOUBLE PRECISION :: Del_H,Diagonal
      DOUBLE PRECISION :: Nx,Ny,Nz

      DOUBLE PRECISION :: old_Del_H,oldNx,oldNy,oldNz

           
      include "function.inc"   

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            IF(.NOT.CUT_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' SCALAR CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'             
               CALL MFIX_EXIT(myPE) 
            ENDIF

            Nx = NORMAL_S(IJK,1)
            Ny = NORMAL_S(IJK,2)
            Nz = NORMAL_S(IJK,3)

            Xref = REFP_S(IJK,1)
            Yref = REFP_S(IJK,2)
            Zref = REFP_S(IJK,3)            

         CASE('U_MOMENTUM')

            IF(.NOT.CUT_U_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' U-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'             
               CALL MFIX_EXIT(myPE) 
            ENDIF

            Nx = NORMAL_U(IJK,1)
            Ny = NORMAL_U(IJK,2)
            Nz = NORMAL_U(IJK,3)

            Xref = REFP_U(IJK,1)
            Yref = REFP_U(IJK,2)
            Zref = REFP_U(IJK,3)         

            IF(WALL_U_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF   

         CASE('V_MOMENTUM')

            IF(.NOT.CUT_V_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' V-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'             
               CALL MFIX_EXIT(myPE) 
            ENDIF

            Nx = NORMAL_V(IJK,1)
            Ny = NORMAL_V(IJK,2)
            Nz = NORMAL_V(IJK,3)

            Xref = REFP_V(IJK,1)
            Yref = REFP_V(IJK,2)
            Zref = REFP_V(IJK,3)            

            IF(WALL_V_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN
            ENDIF   


         CASE('W_MOMENTUM')

            IF(.NOT.CUT_W_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
               WRITE(*,*)' W-MOMENTUM CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' MFiX will exit now.'             
               CALL MFIX_EXIT(myPE) 
            ENDIF

            Nx = NORMAL_W(IJK,1)
            Ny = NORMAL_W(IJK,2)
            Nz = NORMAL_W(IJK,3)

            Xref = REFP_W(IJK,1)
            Yref = REFP_W(IJK,2)
            Zref = REFP_W(IJK,3)     

            IF(WALL_W_AT(IJK)) THEN
               Nx = ZERO
               Ny = ZERO
               Nz = ZERO
               DEL_H = UNDEFINED
               RETURN  
            ENDIF

         CASE DEFAULT
            WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H:'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
            WRITE(*,*)'SCALAR' 
            WRITE(*,*)'U_MOMENTUM' 
            WRITE(*,*)'V_MOMENTUM' 
            WRITE(*,*)'W_MOMENTUM' 
            CALL MFIX_EXIT(myPE)
      END SELECT


      DEL_H = Nx * (X0 - Xref) + Ny * (Y0 - Yref) + Nz * (Z0 - Zref)

!======================================================================
! Negative values of DEL_H are not accepted
!======================================================================

       I = I_OF(IJK)  
       J = J_OF(IJK) 
       K = K_OF(IJK) 

       Diagonal = dsqrt(DX(I)**2 + DY(J)**2 + DZ(K)**2)

      IF (DEL_H <= TOL_DELH * Diagonal) THEN

         DEL_H = UNDEFINED
         Nx = ZERO
         Ny = ZERO
         Nz = ZERO

      ENDIF

      RETURN


      END SUBROUTINE GET_DEL_H

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: STORE_CUT_FACE_INFO                                    C
!  Purpose: Compute and store unit normal vector and reference point   C
!           Defining a cut face                                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE STORE_CUT_FACE_INFO(IJK,TYPE_OF_CELL,N_CUT_FACE_NODES,COORD_CUT_FACE_NODES,X_MEAN,Y_MEAN,Z_MEAN)
    
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
      USE quadric
      USE cutcell
      
      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK
      INTEGER :: N_CUT_FACE_NODES
      DOUBLE PRECISION, DIMENSION(15,3) :: COORD_CUT_FACE_NODES
      DOUBLE PRECISION :: X_MEAN,Y_MEAN,Z_MEAN,F_MEAN
      DOUBLE PRECISION, DIMENSION(3) :: N,V
      DOUBLE PRECISION :: NORM

      include "function.inc"


     IF(NO_K) THEN  ! 2D case

        N(1) = COORD_CUT_FACE_NODES(1,2) - COORD_CUT_FACE_NODES(2,2) ! y1-y2
        N(2) = COORD_CUT_FACE_NODES(2,1) - COORD_CUT_FACE_NODES(1,1) ! x2-x1
        N(3) = ZERO

     ELSE


!======================================================================
! Make sure there are a least three points along the plane
!======================================================================

         IF(N_CUT_FACE_NODES < 3) THEN
            WRITE(*,*)' ERROR IN SUBROUTINE STORE_CUT_FACE_INFO:'
            WRITE(*,*)' CUT FACE HAS LESS THAN 3 NODES.'
            WRITE(*,*)' MFIX WILL EXIT NOW.'
            CALL MFIX_EXIT(myPE)
         END IF

!======================================================================
!  Find tentative unit normal vector
!  and reverse sign if necessary 
!  (unit vector must be pointing towards the fluid)
!======================================================================

         CALL CROSS_PRODUCT(COORD_CUT_FACE_NODES(2,:)-COORD_CUT_FACE_NODES(1,:),&
                            COORD_CUT_FACE_NODES(3,:)-COORD_CUT_FACE_NODES(1,:),N)

      ENDIF


      NORM = DSQRT(N(1)**2 + N(2)**2 + N(3)**2)
      N = N / NORM

      V(1) = X_MEAN - COORD_CUT_FACE_NODES(1,1)
      V(2) = Y_MEAN - COORD_CUT_FACE_NODES(1,2)
      V(3) = Z_MEAN - COORD_CUT_FACE_NODES(1,3)


      IF (DOT_PRODUCT(N,V) < ZERO) N = - N

!======================================================================
! Store unit normal vector and reference point       
!======================================================================

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            NORMAL_S(IJK,:) = N
            REFP_S(IJK,:)   = COORD_CUT_FACE_NODES(1,:)

         CASE('U_MOMENTUM')

            NORMAL_U(IJK,:) = N
            REFP_U(IJK,:)   = COORD_CUT_FACE_NODES(1,:)

         CASE('V_MOMENTUM')

            NORMAL_V(IJK,:) = N
            REFP_V(IJK,:)   = COORD_CUT_FACE_NODES(1,:)

         CASE('W_MOMENTUM')

            NORMAL_W(IJK,:) = N
            REFP_W(IJK,:)   = COORD_CUT_FACE_NODES(1,:)

         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: STORE_CUT_FACE_INFO'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
            WRITE(*,*)'SCALAR' 
            WRITE(*,*)'U_MOMENTUM' 
            WRITE(*,*)'V_MOMENTUM' 
            WRITE(*,*)'W_MOMENTUM' 
            CALL MFIX_EXIT(myPE)
      END SELECT
          

      RETURN
      
      END SUBROUTINE STORE_CUT_FACE_INFO




