!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  GET_DEL_H                                             C
!  Purpose: Finds the normal distance and unit vector from a cut face  C
!  to any point (x0,y0,z0)                                             C
!  The unit normal vector points away from the boundary,               C
!  towards the fluid.                                                  C                              
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

       IF(NO_K) THEN
          Diagonal = dsqrt(DX(I)**2 + DY(J)**2 )
       ELSE
          Diagonal = dsqrt(DX(I)**2 + DY(J)**2 + DZ(K)**2)
       ENDIF

      IF (DEL_H <= TOL_DELH * Diagonal) THEN

         DEL_H = UNDEFINED
         Nx = ZERO
         Ny = ZERO
         Nz = ZERO

      ENDIF

      RETURN


      END SUBROUTINE GET_DEL_H

  SUBROUTINE GET_DEL_H_DES(IJK,TYPE_OF_CELL,X0,Y0,Z0,Del_H,Nx,Ny,Nz, allow_neg_dist)
    
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
      LOGICAL :: ALLOW_NEG_DIST

           
      include "function.inc"   

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            IF(.NOT.CUT_CELL_AT(IJK)) THEN
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
               WRITE(*,*)' SCALAR CELL',IJK,' IS NOT A CUT CELL'
               WRITE(*,*)' I, J, K =',I_OF(IJK), J_OF(IJK), K_OF(IJK)
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
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
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
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
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
               WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
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
            WRITE(*,*)' EROR IN SUBROUTINE GET_DEL_H_DES:'
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

       IF(.NOT.ALLOW_NEG_DIST) THEN 
          Diagonal = dsqrt(DX(I)**2 + DY(J)**2 + DZ(K)**2)
          
          IF (DEL_H <= TOL_DELH * Diagonal) THEN
             DEL_H = UNDEFINED
             Nx = ZERO
             Ny = ZERO
             Nz = ZERO

          ENDIF
       ENDIF

      RETURN


      END SUBROUTINE GET_DEL_H_DES
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
      DOUBLE PRECISION, DIMENSION(3) :: N,V,N1,N2
      DOUBLE PRECISION :: NORM,N1_dot_N2

      INTEGER :: NODE, Q_ID

      DOUBLE PRECISION :: F_Q
      LOGICAL:: CLIP_FLAG
    

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



      IF(N_CUT_FACE_NODES > 3) THEN     ! FOR 3D geometry, check normal of plane defined by nodes 1,2, and 4

         N1 = N  ! Keep copy of previous N (nodes 1,2,3)

         CALL CROSS_PRODUCT(COORD_CUT_FACE_NODES(2,:)-COORD_CUT_FACE_NODES(1,:),&
                             COORD_CUT_FACE_NODES(4,:)-COORD_CUT_FACE_NODES(1,:),N2)


         NORM = DSQRT(N2(1)**2 + N2(2)**2 + N2(3)**2)
         N2 = N2 / NORM

         V(1) = X_MEAN - COORD_CUT_FACE_NODES(1,1)
         V(2) = Y_MEAN - COORD_CUT_FACE_NODES(1,2)
         V(3) = Z_MEAN - COORD_CUT_FACE_NODES(1,3)


         IF (DOT_PRODUCT(N2,V) < ZERO) N2 = - N2


      ENDIF

      N1_dot_N2 = DOT_PRODUCT(N1,N2)
      DEBUG_CG(IJK,1)=N1_dot_N2

      IF(N1_dot_N2<0.99) THEN

!         What should be done when the unit vectors are different ?

      ENDIF


!======================================================================
! Store unit normal vector and reference point       
!======================================================================

      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')

            NORMAL_S(IJK,:) = N
            REFP_S(IJK,:)   = COORD_CUT_FACE_NODES(1,:)

            IF(DO_K) CALL TEST_DEL_H(IJK,'SCALAR') ! test for negative del_H

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

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: TEST_DEL_H                                             C
!  Purpose: tests the computation of wall distance                     C
!           If a negative distance is detected, the normal vector      C
!           is inverted                                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 22-Feb-12  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE TEST_DEL_H(IJK,TYPE_OF_CELL)
    
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
      DOUBLE PRECISION:: X_COPY,Y_COPY,Z_COPY
      INTEGER :: IJK
      INTEGER :: NODE,N_N1,N_N2,N
      DOUBLE PRECISION :: Del_H
      DOUBLE PRECISION :: Nx,Ny,Nz

      LOGICAL :: ALLOW_NEG_DIST = .TRUE.  ! forces GET_DEL_H_DES to output negative delh
                                           ! i.e. do not let the subroutine overwrite negative values


! This subroutine tests values of del_H for nodes defining a cut cell.
! Only nodes that are in the fluid region are tested.
! Nodes belonging to the cut face should return zero (or near zero) values and are not tested.
! If a negative del_H is detected, the unit normal vector is reversed.

      IF(NO_K) THEN
         N_N1 = 5
         N_N2 = 8
      ELSE
         N_N1 = 1
         N_N2 = 8
      ENDIF

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      DO NODE = 1,NUMBER_OF_NODES(IJK)
         IF(CONNECTIVITY(IJK,NODE)<=IJKEND3) THEN  ! node does not belong to the cut-face
                                                    ! i.e. is in the fluid region

            DO N = N_N1,N_N2                     ! get node coordinate
               IF(CONNECTIVITY(IJK,NODE) == IJK_OF_NODE(N)) THEN
                  X_COPY = X_NODE(N)
                  Y_COPY = Y_NODE(N)
                  Z_COPY = Z_NODE(N)
                  EXIT
               ENDIF
            ENDDO

!           Compute del_H
            CALL GET_DEL_H_DES(IJK,TYPE_OF_CELL,X_COPY,Y_COPY,Z_COPY,Del_H,Nx,Ny,Nz, ALLOW_NEG_DIST) 


            IF(DEL_H<ZERO) THEN


               IF(PRINT_WARNINGS.AND.MyPE==PE_IO) THEN
                  WRITE(*,*),' Warning: Negative delh detected in scalar cell :',IJK
                  WRITE(*,*) ' Location (X,Y,Z) = ',X_COPY,Y_COPY,Z_COPY
                  WRITE(*,*) ' Reverting unit normal vector.'
               ENDIF

               NORMAL_S(IJK,1) = -NORMAL_S(IJK,1)
               NORMAL_S(IJK,2) = -NORMAL_S(IJK,2)
               NORMAL_S(IJK,3) = -NORMAL_S(IJK,3)


            ENDIF


         ENDIF
      ENDDO

      RETURN

      END SUBROUTINE TEST_DEL_H







!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  GET_DISTANCE_TO_WALL                                  C
!  Purpose: Finds the distance fraom any scalar cell to the closest    C
!  wall                                                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 16-May-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE GET_DISTANCE_TO_WALL
    
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
!      CHARACTER (LEN=*) :: TYPE_OF_CELL
      DOUBLE PRECISION:: X0,Y0,Z0,XREF,YREF,ZREF
      INTEGER :: IJK,I,J,K,N,IJK_CUT
      DOUBLE PRECISION :: Del_H,Diagonal
      DOUBLE PRECISION :: Nx,Ny,Nz

      DOUBLE PRECISION :: D_TO_CUT

      INTEGER :: N_CUT_CELLS
      INTEGER :: LIST_OF_CUT_CELLS(DIMENSION_3)

           
      include "function.inc"   

!======================================================================
!  Get a list of cut cells 
!======================================================================

      N_CUT_CELLS = 0

      DO IJK = IJKSTART3, IJKEND3

         IF(CUT_CELL_AT(IJK)) THEN

         N_CUT_CELLS = N_CUT_CELLS + 1

         LIST_OF_CUT_CELLS(N_CUT_CELLS) = IJK
            

         ENDIF

      ENDDO


!======================================================================
!  Brute force: Loop through all scalar cells
!  compute the distance to each cut cell and take the minimum  
!======================================================================


      DO IJK = IJKSTART3, IJKEND3


         IF(INTERIOR_CELL_AT(IJK)) THEN

            DWALL(IJK) = UNDEFINED

!======================================================================
!  Get coordinates of cell center
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'SCALAR')

            X0 = X_NODE(0)
            Y0 = Y_NODE(0)
            Z0 = Z_NODE(0)

            DO N = 1,N_CUT_CELLS

               IJK_CUT = LIST_OF_CUT_CELLS(N)

               Xref = REFP_S(IJK_CUT,1)
               Yref = REFP_S(IJK_CUT,2)
               Zref = REFP_S(IJK_CUT,3)  

               D_TO_CUT = DSQRT((X0 - Xref)**2 + (Y0 - Yref)**2 + (Z0 - Zref)**2)


               IF(D_TO_CUT<DWALL(IJK)) DWALL(IJK) = D_TO_CUT


            ENDDO


         ENDIF

      END DO



      RETURN


      END SUBROUTINE GET_DISTANCE_TO_WALL








