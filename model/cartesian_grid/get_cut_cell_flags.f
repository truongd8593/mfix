!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_CELL_FLAGS                                  C
!  Purpose: Set flags for scalar cut cells, based on intersection      C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_CELL_FLAGS 
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      Use vtk
      USE polygon
      USE stl

      USE physprop
      USE fldvar
      USE scalars
      USE funits 
      USE rxns


      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,L
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,N_N1,N_N2,Q_ID,Q_ID2
      INTEGER :: MIN_INTERSECTIONS,MAX_INTERSECTIONS,N_N
      LOGICAL :: CLIP_FLAG
      DOUBLE PRECISION :: MIN_VOL,MAX_VOL
      DOUBLE PRECISION :: MIN_AYZ,MAX_AYZ
      DOUBLE PRECISION :: MIN_AXZ,MAX_AXZ
      DOUBLE PRECISION :: MIN_AXY,MAX_AXY
      DOUBLE PRECISION :: F_NODE_02
      INTEGER :: BCID



      include "function.inc"


      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'INTERSECTING GEOMETRY WITH SCALAR CELLS...'
      ENDIF
10    FORMAT(1X,A)

      GLOBAL_VAR_ALLOCATED = .FALSE.
      GRID_INFO_PRINTED_ON_SCREEN = .FALSE.

!  Set arrays for computing indices
!
      
      CALL SET_INCREMENTS 

      DX(IMAX3+1) = DX(IMAX3)
      DY(JMAX3+1) = DY(JMAX3)
      DZ(KMAX3+1) = DZ(KMAX3)

!     x-Location of U-momentum cells for original (uncut grid)
      IF (DO_I) THEN 
         XG_E(1) = ZERO
         DO I = IMIN1, IMAX3
            XG_E(I) = XG_E(I-1) + DX(I) 
         END DO 
      ENDIF
 

!     y-Location of V-momentum cells for original (uncut grid)
      IF (DO_J) THEN 
         YG_N(1) = ZERO
         DO J = JMIN1, JMAX3  
            YG_N(J) = YG_N(J-1) + DY(J) 
         END DO 
      ENDIF


!     z-Location of W-momentum cells for original (uncut grid)
      IF (DO_K) THEN 
         ZG_T(1) = ZERO
         DO K = KMIN1, KMAX3  
            ZG_T(K) = ZG_T(K-1) + DZ(K) 
         END DO 
      ENDIF

      PARTITION = DFLOAT(myPE)       ! ASSIGN processor ID (for vizualisation)

      CALL SET_GEOMETRY1   ! INITIALIZE ALL VOLUMES AND AREAS


!======================================================================
!  Grid coordinates:
!======================================================================

      SMALL_CELL_AT = .FALSE.
      SMALL_CELL_FLAG = 0
      NUMBER_OF_SMALL_CELLS = 0

!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.


      DO IJK = IJKSTART3, IJKEND3

         I = I_OF(IJK) 
         J = J_OF(IJK) 
         K = K_OF(IJK) 

         IF(NO_K) THEN   ! 2D case


            INTERIOR_CELL_AT(IJK) = (     (I >= ISTART1 ).AND.(I <= IEND1 )  & 
                                     .AND.(J >= JSTART1 ).AND.(J <= JEND1 ) )

         ELSE            ! 3D case

            INTERIOR_CELL_AT(IJK) = (     (I >= ISTART1 ).AND.(I <= IEND1 )  & 
                                     .AND.(J >= JSTART1 ).AND.(J <= JEND1 )  &
                                     .AND.(K >= KSTART1 ).AND.(K <= KEND1 ) )

         ENDIF


      END DO


      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
            CALL INTERSECT(IJK,'SCALAR',Xn_int(IJK),Ye_int(IJK),Zt_int(IJK))
         END DO
      ELSE
         CALL CAD_INTERSECT('SCALAR',Xn_int,Ye_int,Zt_int)
      ENDIF
!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

!         IF(INTERIOR_CELL_AT(IJK)) THEN

            CALL CLEAN_INTERSECT(IJK,'SCALAR',Xn_int(IJK),Ye_int(IJK),Zt_int(IJK))

!         ENDIF

      END DO

      call SEND_RECEIVE_1D_LOGICAL(SNAP,2)

       NUMBER_OF_NODES = 0

      DO IJK = IJKSTART3, IJKEND3


         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

         CALL GET_CELL_NODE_COORDINATES(IJK,'SCALAR')

!======================================================================
!  Initialize location of velocity nodes
!======================================================================

         X_U(IJK) = X_NODE(8)
         Y_U(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
         Z_U(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))
  
         X_V(IJK) = HALF * (X_NODE(7) + X_NODE(8))
         Y_V(IJK) = Y_NODE(8)
         Z_V(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

         X_W(IJK) = HALF * (X_NODE(7) + X_NODE(8))
         Y_W(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
         Z_W(IJK) = Z_NODE(8)

         IF(INTERIOR_CELL_AT(IJK)) THEN
!======================================================================
!  Get Connectivity
!======================================================================

            IF(NO_K) THEN
               MIN_INTERSECTIONS = 2
               MAX_INTERSECTIONS = 2
               N_N1 = 5
               N_N2 = 8

            ELSE
               MIN_INTERSECTIONS = 3
               MAX_INTERSECTIONS = 6
               N_N1 = 1
               N_N2 = 8
            ENDIF

            CALL GET_CONNECTIVITY(IJK,'SCALAR',NUMBER_OF_NEW_POINTS,NUMBER_OF_NODES(IJK),CONNECTIVITY,&
            X_NEW_POINT,Y_NEW_POINT,Z_NEW_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_int,Ye_int,Zt_int)

            IF(TOTAL_NUMBER_OF_INTERSECTIONS < MIN_INTERSECTIONS ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('SCALAR',IJK,0,F_NODE(0),CLIP_FLAG,BCID)


               IF(F_NODE(0) < ZERO) THEN
                  IF((FLAG(IJK)>=100).AND.(FLAG(IJK)<=102)) THEN
                     BLOCKED_CELL_AT(IJK) = .TRUE.
                     STANDARD_CELL_AT(IJK) = .FALSE.          ! Blocked cell = wall cell
                  ELSE
                     BLOCKED_CELL_AT(IJK) = .FALSE.
                     STANDARD_CELL_AT(IJK) = .TRUE.           ! Regular fluid cell
                  ENDIF
               ELSE
                  FLAG(IJK) = 100
                  BLOCKED_CELL_AT(IJK) = .TRUE.               ! Blocked fluid cell
                  STANDARD_CELL_AT(IJK) = .FALSE.
                  AXY(IJK) = ZERO
                  AXZ(IJK) = ZERO
                  AYZ(IJK) = ZERO
                  VOL(IJK) = ZERO

                  AXY(BOTTOM_OF(IJK)) = ZERO
                  AXZ(SOUTH_OF(IJK)) = ZERO
                  AYZ(WEST_OF(IJK)) = ZERO


               ENDIF

               IF(NO_K) THEN
                  NUMBER_OF_NODES(IJK) = 4
                  CONNECTIVITY(IJK,1) = IJK_OF_NODE(5)
                  CONNECTIVITY(IJK,2) = IJK_OF_NODE(6)
                  CONNECTIVITY(IJK,3) = IJK_OF_NODE(8)
                  CONNECTIVITY(IJK,4) = IJK_OF_NODE(7)
               ELSE
                  NUMBER_OF_NODES(IJK) = 8
                  DO NODE = N_N1,N_N2
                     CONNECTIVITY(IJK,NODE) = IJK_OF_NODE(NODE)
                  END DO
               ENDIF


            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > MAX_INTERSECTIONS ) THEN
               WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN CELL IJK = ',IJK
               WRITE(*,*) 'MAXIMUM NUMBER OF INTERSECTIONS = ',MAX_INTERSECTIONS
               WRITE(*,*) 'TOTAL NUMBER OF INTERSECTIONS = ',TOTAL_NUMBER_OF_INTERSECTIONS
               WRITE(*,*) 'MFIX WILL EXIT NOW.'
               CALL MFIX_EXIT(MYPE) 

            ELSE                                         ! Cut cell


               CUT_CELL_AT(IJK) = .TRUE.
               BLOCKED_CELL_AT(IJK) = .FALSE.
               STANDARD_CELL_AT(IJK) = .FALSE.


               DO NODE = N_N1,N_N2
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     IF(.NOT.SNAP(IJK_OF_NODE(NODE))) THEN
                        NUMBER_OF_NODES(IJK) = NUMBER_OF_NODES(IJK) + 1
                        CONNECTIVITY(IJK,NUMBER_OF_NODES(IJK)) = IJK_OF_NODE(NODE)
                     ENDIF
                  ENDIF
               END DO

               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'SCALAR',NUMBER_OF_NODES(IJK),CONNECTIVITY,&
               X_NEW_POINT,Y_NEW_POINT,Z_NEW_POINT)

            ENDIF

         ENDIF
      END DO

      RETURN
      END SUBROUTINE SET_3D_CUT_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_U_CELL_FLAGS                                C
!  Purpose: Set flags for U-Momentum cut cells, based on intersection  C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_U_CELL_FLAGS
    
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
      USE polygon
      USE stl
      
      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,N_N1,N_N2,Q_ID,Q_ID2
      INTEGER :: MIN_INTERSECTIONS,MAX_INTERSECTIONS,N_N
      LOGICAL :: CLIP_FLAG
      DOUBLE PRECISION :: MIN_VOL,MAX_VOL
      DOUBLE PRECISION :: MIN_AYZ,MAX_AYZ
      DOUBLE PRECISION :: MIN_AXZ,MAX_AXZ
      DOUBLE PRECISION :: MIN_AXY,MAX_AXY
      DOUBLE PRECISION :: MIN_CUT,MAX_CUT
      DOUBLE PRECISION :: F_NODE_02
      INTEGER :: BCID

      include "function.inc"

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'INTERSECTING GEOMETRY WITH U-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)
!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.
      TOL_SNAP = ZERO


      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
            CALL INTERSECT(IJK,'U_MOMENTUM',Xn_U_int(IJK),Ye_U_int(IJK),Zt_U_int(IJK))
         END DO
      ELSE
         CALL CAD_INTERSECT('U_MOMENTUM',Xn_U_int,Ye_U_int,Zt_U_int)
      ENDIF

!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

         IF(INTERIOR_CELL_AT(IJK)) THEN

            CALL CLEAN_INTERSECT(IJK,'U_MOMENTUM',Xn_U_int(IJK),Ye_U_int(IJK),Zt_U_int(IJK))

         ENDIF

      END DO


      NUMBER_OF_NEW_U_POINTS = 0

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'U_MOMENTUM')

!======================================================================
!  Initialize location of velocity nodes at center of E,W,T faces
!======================================================================

                  X_U_ec(IJK) = X_NODE(8)
                  Y_U_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_U_nc(IJK) = Y_NODE(8)
                  Z_U_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8)) 
                  Y_U_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_tc(IJK) = Z_NODE(8)

!======================================================================
!  Get Connectivity
!======================================================================

            IF(NO_K) THEN
               MIN_INTERSECTIONS = 2
               MAX_INTERSECTIONS = 2
               N_N1 = 5
               N_N2 = 8

            ELSE
               MIN_INTERSECTIONS = 3
               MAX_INTERSECTIONS = 6
               N_N1 = 1
               N_N2 = 8
            ENDIF

            CALL GET_CONNECTIVITY(IJK,'U_MOMENTUM',NUMBER_OF_NEW_U_POINTS,NUMBER_OF_U_NODES(IJK),CONNECTIVITY_U,&
            X_NEW_U_POINT,Y_NEW_U_POINT,Z_NEW_U_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_U_int,Ye_U_int,Zt_U_int)


            IF(TOTAL_NUMBER_OF_INTERSECTIONS < MIN_INTERSECTIONS ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('U_MOMENTUM',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(F_NODE(0) < ZERO) THEN
                  BLOCKED_U_CELL_AT(IJK) = .FALSE.
                  STANDARD_U_CELL_AT(IJK) = .TRUE.          ! Regular fluid cell

                  X_U_ec(IJK) = X_NODE(8)
                  Y_U_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_U_nc(IJK) = Y_NODE(8)
                  Z_U_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_U_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8)) 
                  Y_U_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_U_tc(IJK) = Z_NODE(8)

               ELSE
                  BLOCKED_U_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
                  STANDARD_U_CELL_AT(IJK) = .FALSE.
                  AXY_U(IJK) = ZERO
                  AXZ_U(IJK) = ZERO
                  AYZ_U(IJK) = ZERO
                  VOL_U(IJK) = ZERO
               ENDIF


               IF(NO_K) THEN
                  NUMBER_OF_U_NODES(IJK) = 4
                  CONNECTIVITY_U(IJK,1) = IJK_OF_NODE(5)
                  CONNECTIVITY_U(IJK,2) = IJK_OF_NODE(6)
                  CONNECTIVITY_U(IJK,3) = IJK_OF_NODE(8)
                  CONNECTIVITY_U(IJK,4) = IJK_OF_NODE(7)
               ELSE
                  NUMBER_OF_U_NODES(IJK) = 8
                  DO NODE = N_N1,N_N2
                     CONNECTIVITY_U(IJK,NODE) = IJK_OF_NODE(NODE)
                  END DO
               ENDIF


            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > MAX_INTERSECTIONS ) THEN
               WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN U-CELL IJK = ',IJK
               WRITE(*,*) 'TOTAL NUMBER OF INTERSECTIONS = ',TOTAL_NUMBER_OF_INTERSECTIONS
               WRITE(*,*) 'MFIX WILL EXIT NOW.'
               CALL MFIX_EXIT(MYPE) 

            ELSE                                         ! Cut cell

               CUT_U_CELL_AT(IJK) = .TRUE.
               BLOCKED_U_CELL_AT(IJK) = .FALSE.
               STANDARD_U_CELL_AT(IJK) = .FALSE.


               DO NODE = N_N1,N_N2
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     NUMBER_OF_U_NODES(IJK) = NUMBER_OF_U_NODES(IJK) + 1
                     CONNECTIVITY_U(IJK,NUMBER_OF_U_NODES(IJK)) = IJK_OF_NODE(NODE)
                  ENDIF
               END DO


               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'U_MOMENTUM',NUMBER_OF_U_NODES(IJK),CONNECTIVITY_U,&
               X_NEW_U_POINT,Y_NEW_U_POINT,Z_NEW_U_POINT)

            ENDIF

         ENDIF

      END DO



      DO IJK = IJKSTART3, IJKEND3

         CALL GET_CELL_NODE_COORDINATES(IJK,'U_MOMENTUM')

         DELX_Ue(IJK) = X_NODE(8) - X_U(IJK)
         DELX_Uw(IJK) = X_U(IJK) - X_NODE(1)

         DELY_Un(IJK) = Y_NODE(8) - Y_U(IJK)
         DELY_Us(IJK) = Y_U(IJK) - Y_NODE(1)

         DELZ_Ut(IJK) = Z_NODE(8) - Z_U(IJK)
         DELZ_Ub(IJK) = Z_U(IJK) - Z_NODE(1)

      ENDDO

      RETURN

      
      END SUBROUTINE SET_3D_CUT_U_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_V_CELL_FLAGS                                C
!  Purpose: Set flags for V-Momentum cut cells, based on intersection  C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_V_CELL_FLAGS
    
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
      USE polygon
      USE stl
      
      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,N_N1,N_N2,Q_ID,Q_ID2
      INTEGER :: MIN_INTERSECTIONS,MAX_INTERSECTIONS,N_N
      LOGICAL :: CLIP_FLAG
      DOUBLE PRECISION :: MIN_VOL,MAX_VOL
      DOUBLE PRECISION :: MIN_AYZ,MAX_AYZ
      DOUBLE PRECISION :: MIN_AXZ,MAX_AXZ
      DOUBLE PRECISION :: MIN_AXY,MAX_AXY
      DOUBLE PRECISION :: MIN_CUT,MAX_CUT
      DOUBLE PRECISION :: F_NODE_02
      INTEGER :: BCID

      include "function.inc"

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)'INTERSECTING GEOMETRY WITH V-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)
!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.


      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
            CALL INTERSECT(IJK,'V_MOMENTUM',Xn_V_int(IJK),Ye_V_int(IJK),Zt_V_int(IJK))
         END DO
      ELSE
         CALL CAD_INTERSECT('V_MOMENTUM',Xn_V_int,Ye_V_int,Zt_V_int)
      ENDIF

!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

         IF(INTERIOR_CELL_AT(IJK)) THEN

            CALL CLEAN_INTERSECT(IJK,'V_MOMENTUM',Xn_V_int(IJK),Ye_V_int(IJK),Zt_V_int(IJK))

         ENDIF

      END DO

      NUMBER_OF_NEW_V_POINTS = 0

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'V_MOMENTUM')

!======================================================================
!  Initialize location of velocity nodes at center of E,W,T faces
!======================================================================

            X_V_ec(IJK) = X_NODE(8)
            Y_V_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_V_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_V_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_V_nc(IJK) = Y_NODE(8)
            Z_V_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_V_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8)) 
            Y_V_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_V_tc(IJK) = Z_NODE(8)

!======================================================================
!  Get Connectivity
!======================================================================

            IF(NO_K) THEN
               MIN_INTERSECTIONS = 2
               MAX_INTERSECTIONS = 2
               N_N1 = 5
               N_N2 = 8

            ELSE
               MIN_INTERSECTIONS = 3
               MAX_INTERSECTIONS = 6
               N_N1 = 1
               N_N2 = 8
            ENDIF


            CALL GET_CONNECTIVITY(IJK,'V_MOMENTUM',NUMBER_OF_NEW_V_POINTS,NUMBER_OF_V_NODES(IJK),CONNECTIVITY_V,&
            X_NEW_V_POINT,Y_NEW_V_POINT,Z_NEW_V_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_V_int,Ye_V_int,Zt_V_int)


            IF(TOTAL_NUMBER_OF_INTERSECTIONS < MIN_INTERSECTIONS ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('V_MOMENTUM',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(F_NODE(0) < ZERO) THEN
                  BLOCKED_V_CELL_AT(IJK) = .FALSE.
                  STANDARD_V_CELL_AT(IJK) = .TRUE.          ! Regular fluid cell

                  X_V_ec(IJK) = X_NODE(8)
                  Y_V_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_V_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_V_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_V_nc(IJK) = Y_NODE(8)
                  Z_V_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_V_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8)) 
                  Y_V_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_V_tc(IJK) = Z_NODE(8)


               ELSE
                  BLOCKED_V_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
                  STANDARD_V_CELL_AT(IJK) = .FALSE.
                  AXY_V(IJK) = ZERO
                  AXZ_V(IJK) = ZERO
                  AYZ_V(IJK) = ZERO
                  VOL_V(IJK) = ZERO
               ENDIF

               IF(NO_K) THEN
                  NUMBER_OF_V_NODES(IJK) = 4
                  CONNECTIVITY_V(IJK,1) = IJK_OF_NODE(5)
                  CONNECTIVITY_V(IJK,2) = IJK_OF_NODE(6)
                  CONNECTIVITY_V(IJK,3) = IJK_OF_NODE(8)
                  CONNECTIVITY_V(IJK,4) = IJK_OF_NODE(7)
               ELSE
                  NUMBER_OF_V_NODES(IJK) = 8
                  DO NODE = N_N1,N_N2
                     CONNECTIVITY_V(IJK,NODE) = IJK_OF_NODE(NODE)
                  END DO
               ENDIF



            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > MAX_INTERSECTIONS ) THEN
               WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN V-CELL IJK = ',IJK
               WRITE(*,*) 'MFIX WILL EXIT NOW.'
               CALL MFIX_EXIT(MYPE) 

            ELSE                                         ! Cut cell

               CUT_V_CELL_AT(IJK) = .TRUE.
               BLOCKED_V_CELL_AT(IJK) = .FALSE.
               STANDARD_V_CELL_AT(IJK) = .FALSE.

               DO NODE = N_N1,N_N2
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     NUMBER_OF_V_NODES(IJK) = NUMBER_OF_V_NODES(IJK) + 1
                     CONNECTIVITY_V(IJK,NUMBER_OF_V_NODES(IJK)) = IJK_OF_NODE(NODE)
                  ENDIF
               END DO


               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'V_MOMENTUM',NUMBER_OF_V_NODES(IJK),CONNECTIVITY_V,&
               X_NEW_V_POINT,Y_NEW_V_POINT,Z_NEW_V_POINT)

            ENDIF

         ENDIF
      END DO

      DO IJK = IJKSTART3, IJKEND3

         CALL GET_CELL_NODE_COORDINATES(IJK,'V_MOMENTUM')

         DELX_Ve(IJK) = X_NODE(8) - X_V(IJK)
         DELX_Vw(IJK) = X_V(IJK) - X_NODE(1)
 
         DELY_Vn(IJK) = Y_NODE(8) - Y_V(IJK)
         DELY_Vs(IJK) = Y_V(IJK) - Y_NODE(1)

         DELZ_Vt(IJK) = Z_NODE(8) - Z_V(IJK)
         DELZ_Vb(IJK) = Z_V(IJK) - Z_NODE(1)


      ENDDO


      RETURN

      
      END SUBROUTINE SET_3D_CUT_V_CELL_FLAGS

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_W_CELL_FLAGS                                C
!  Purpose: Set flags for W-Momentum cut cells, based on intersection  C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_W_CELL_FLAGS
    
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
      USE polygon
      USE stl
      
      IMPLICIT NONE
      INTEGER :: IJK,I,J,K
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE,Q_ID,Q_ID2
      LOGICAL :: CLIP_FLAG
      DOUBLE PRECISION :: MIN_VOL,MAX_VOL
      DOUBLE PRECISION :: MIN_AYZ,MAX_AYZ
      DOUBLE PRECISION :: MIN_AXZ,MAX_AXZ
      DOUBLE PRECISION :: MIN_AXY,MAX_AXY
      DOUBLE PRECISION :: MIN_CUT,MAX_CUT
      DOUBLE PRECISION :: F_NODE_02
      INTEGER :: BCID

      include "function.inc"

      IF(MyPE == PE_IO) THEN
         WRITE(*,10)'INTERSECTING GEOMETRY WITH W-MOMENTUM CELLS...'
      ENDIF
10    FORMAT(1X,A)
!======================================================================
!  Intersection between quadric Grid
!======================================================================

      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.
      SNAP = .FALSE.
      
      IF(.NOT.(USE_STL.OR.USE_MSH)) THEN
         DO IJK = IJKSTART3, IJKEND3
            CALL INTERSECT(IJK,'W_MOMENTUM',Xn_W_int(IJK),Ye_W_int(IJK),Zt_W_int(IJK))
         END DO
      ELSE
         CALL CAD_INTERSECT('W_MOMENTUM',Xn_W_int,Ye_W_int,Zt_W_int)
      ENDIF

!======================================================================
!  Clean-up intersection flags in preparaton of small cells removal
!======================================================================

      DO IJK = IJKSTART3, IJKEND3

         IF(INTERIOR_CELL_AT(IJK)) THEN

            CALL CLEAN_INTERSECT(IJK,'W_MOMENTUM',Xn_W_int(IJK),Ye_W_int(IJK),Zt_W_int(IJK))

         ENDIF

      END DO


      NUMBER_OF_NEW_W_POINTS = 0

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

            CALL GET_CELL_NODE_COORDINATES(IJK,'W_MOMENTUM')

!======================================================================
!  Initialize location of velocity nodes at center of E,W,T faces
!======================================================================

            X_W_ec(IJK) = X_NODE(8)
            Y_W_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_W_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_W_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
            Y_W_nc(IJK) = Y_NODE(8)
            Z_W_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

            X_W_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8)) 
            Y_W_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
            Z_W_tc(IJK) = Z_NODE(8)

!======================================================================
!  Get Connectivity
!======================================================================

            CALL GET_CONNECTIVITY(IJK,'W_MOMENTUM',NUMBER_OF_NEW_W_POINTS,NUMBER_OF_W_NODES(IJK),CONNECTIVITY_W,&
            X_NEW_W_POINT,Y_NEW_W_POINT,Z_NEW_W_POINT,TOTAL_NUMBER_OF_INTERSECTIONS,Xn_W_int,Ye_W_int,Zt_W_int)


            IF(TOTAL_NUMBER_OF_INTERSECTIONS < 3 ) THEN   ! Not a cut cell

               Q_ID = 1
               CALL EVAL_F('QUADRIC',X_NODE(0),Y_NODE(0),Z_NODE(0),Q_ID,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('POLYGON',X_NODE(0),Y_NODE(0),Z_NODE(0),N_POLYGON,F_NODE(0),CLIP_FLAG)

               CALL EVAL_F('USR_DEF',X_NODE(0),Y_NODE(0),Z_NODE(0),N_USR_DEF,F_NODE(0),CLIP_FLAG)

               CALL EVAL_STL_FCT_AT('W_MOMENTUM',IJK,0,F_NODE(0),CLIP_FLAG,BCID)

               IF(F_NODE(0) < ZERO) THEN
                  BLOCKED_W_CELL_AT(IJK) = .FALSE.
                  STANDARD_W_CELL_AT(IJK) = .TRUE.          ! Regular fluid cell

                  X_W_ec(IJK) = X_NODE(8)
                  Y_W_ec(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_W_ec(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_W_nc(IJK) = HALF * (X_NODE(7) + X_NODE(8))
                  Y_W_nc(IJK) = Y_NODE(8)
                  Z_W_nc(IJK) = HALF * (Z_NODE(4) + Z_NODE(8))

                  X_W_tc(IJK) = HALF * (X_NODE(7) + X_NODE(8)) 
                  Y_W_tc(IJK) = HALF * (Y_NODE(6) + Y_NODE(8))
                  Z_W_tc(IJK) = Z_NODE(8)


               ELSE
                  BLOCKED_W_CELL_AT(IJK) = .TRUE.           ! Blocked fluid cell
                  STANDARD_W_CELL_AT(IJK) = .FALSE.
                  AXY_W(IJK) = ZERO
                  AXZ_W(IJK) = ZERO
                  AYZ_W(IJK) = ZERO
                  VOL_W(IJK) = ZERO
               ENDIF

               NUMBER_OF_W_NODES(IJK) = 8
               DO NODE = 1,8
                  CONNECTIVITY_W(IJK,NODE) = IJK_OF_NODE(NODE)
               END DO

            ELSE IF(TOTAL_NUMBER_OF_INTERSECTIONS > 6 ) THEN
               WRITE(*,*) 'TOO MANY INTERSECTIONS FOUND IN W-CELL IJK = ',IJK
               WRITE(*,*) 'MFIX WILL EXIT NOW.'
               CALL MFIX_EXIT(MYPE) 

            ELSE                                         ! Cut cell

               CUT_W_CELL_AT(IJK) = .TRUE.

               DO NODE = 1,8
                  IF(F_NODE(NODE) < - TOL_F) THEN
                     NUMBER_OF_W_NODES(IJK) = NUMBER_OF_W_NODES(IJK) + 1
                     CONNECTIVITY_W(IJK,NUMBER_OF_W_NODES(IJK)) = IJK_OF_NODE(NODE)
                  ENDIF
               END DO


               CALL GET_CUT_CELL_VOLUME_AND_AREAS(IJK,'W_MOMENTUM',NUMBER_OF_W_NODES(IJK),CONNECTIVITY_W,&
               X_NEW_W_POINT,Y_NEW_W_POINT,Z_NEW_W_POINT)

            ENDIF

         ENDIF
      END DO

      DO IJK = IJKSTART3, IJKEND3

         CALL GET_CELL_NODE_COORDINATES(IJK,'W_MOMENTUM')

         DELX_We(IJK) = X_NODE(8) - X_W(IJK)
         DELX_Ww(IJK) = X_W(IJK) - X_NODE(1)

         DELY_Wn(IJK) = Y_NODE(8) - Y_W(IJK)
         DELY_Ws(IJK) = Y_W(IJK) - Y_NODE(1)

         DELZ_Wt(IJK) = Z_NODE(8) - Z_W(IJK)
         DELZ_Wb(IJK) = Z_W(IJK) - Z_NODE(1)

      ENDDO



      RETURN

      
      END SUBROUTINE SET_3D_CUT_W_CELL_FLAGS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_3D_CUT_CELL_TREATMENT_FLAGS                        C
!  Purpose: Set flags for scalar cut cells, based on intersection      C
!  of the grid with the quadric(s)                                     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_3D_CUT_CELL_TREATMENT_FLAGS
    
      USE param
      USE param1
      USE parallel
      USE constant
      USE run
      USE toleranc
      USE geometry
      USE indices  
      USE compar
      USE mpi_utility 
      USE sendrecv
      USE quadric
      USE cutcell
      Use vtk


      IMPLICIT NONE
      INTEGER :: IJK,I,J,K,IM,IP,JM,JP,KM,KP
      INTEGER :: IMJK,IPJK,IJMK,IJPK,IJKM,IJKP

      include "function.inc"

      IF(MyPE == PE_IO) THEN
         WRITE(*,*)'SETTING CUT CELL TREATMENT FLAGS...'
      ENDIF
10    FORMAT(1X,A)
!======================================================================
!  Set flags identifying cells requiring cut cell treatment:
!  These are the cut cells and their neighbours
!======================================================================

      CUT_TREATMENT_AT   = .FALSE.
      CUT_U_TREATMENT_AT = .FALSE.
      CUT_V_TREATMENT_AT = .FALSE.
      CUT_W_TREATMENT_AT = .FALSE.
  

      DO IJK = IJKSTART3, IJKEND3

         CALL WRITE_PROGRESS_BAR(IJK,IJKEND3 - IJKSTART3 + 1,'C')

         IF(INTERIOR_CELL_AT(IJK)) THEN

            I = I_OF(IJK) 
            J = J_OF(IJK) 
            K = K_OF(IJK)      

            IP = I + 1
            JP = J + 1
            KP = K + 1

            IM = I - 1 
            JM = J - 1 
            KM = K - 1
    
            IMJK   = FUNIJK(IM,J,K)
            IJMK   = FUNIJK(I,JM,K)
            IJKM   = FUNIJK(I,J,KM)

            IPJK   = FUNIJK(IP,J,K)
            IJPK   = FUNIJK(I,JP,K)
            IJKP   = FUNIJK(I,J,KP)

            CUT_TREATMENT_AT(IJK) = (CUT_CELL_AT(IJK ).OR.   &
                                     CUT_CELL_AT(IMJK).OR.   &
                                     CUT_CELL_AT(IPJK).OR.   &
                                     CUT_CELL_AT(IJMK).OR.   &
                                     CUT_CELL_AT(IJPK))

            CUT_U_TREATMENT_AT(IJK) = (CUT_U_CELL_AT(IJK ).OR.   &
                                       CUT_U_CELL_AT(IMJK).OR.   &
                                       CUT_U_CELL_AT(IPJK).OR.   &
                                       CUT_U_CELL_AT(IJMK).OR.   &
                                       CUT_U_CELL_AT(IJPK))

            CUT_V_TREATMENT_AT(IJK) = (CUT_V_CELL_AT(IJK ).OR.   &
                                       CUT_V_CELL_AT(IMJK).OR.   &
                                       CUT_V_CELL_AT(IPJK).OR.   &
                                       CUT_V_CELL_AT(IJMK).OR.   &
                                       CUT_V_CELL_AT(IJPK))

            IF(DO_K) THEN

               CUT_TREATMENT_AT(IJK) = (CUT_TREATMENT_AT(IJK).OR.   &
                                            CUT_CELL_AT(IJKM).OR.   &
                                            CUT_CELL_AT(IJKP))

               CUT_U_TREATMENT_AT(IJK) = (CUT_U_TREATMENT_AT(IJK).OR.   &
                                              CUT_U_CELL_AT(IJKM).OR.   &
                                              CUT_U_CELL_AT(IJKP))

               CUT_V_TREATMENT_AT(IJK) = (CUT_V_TREATMENT_AT(IJK).OR.   &
                                              CUT_V_CELL_AT(IJKM).OR.   &
                                              CUT_V_CELL_AT(IJKP))

               CUT_W_TREATMENT_AT(IJK) = (CUT_W_CELL_AT(IJK ).OR.   &
                                          CUT_W_CELL_AT(IMJK).OR.   &
                                          CUT_W_CELL_AT(IPJK).OR.   &
                                          CUT_W_CELL_AT(IJMK).OR.   &
                                          CUT_W_CELL_AT(IJPK).OR.   &
                                          CUT_W_CELL_AT(IJKM).OR.   &
                                          CUT_W_CELL_AT(IJKP))


            ENDIF

         ENDIF

      END DO


      IF(CG_SAFE_MODE(1)==1) CUT_TREATMENT_AT   = .FALSE.
      IF(CG_SAFE_MODE(3)==1) CUT_U_TREATMENT_AT = .FALSE.
      IF(CG_SAFE_MODE(4)==1) CUT_V_TREATMENT_AT = .FALSE.
      IF(CG_SAFE_MODE(5)==1) CUT_W_TREATMENT_AT = .FALSE.


      RETURN

      
      END SUBROUTINE SET_3D_CUT_CELL_TREATMENT_FLAGS




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SET_GHOST_CELL_FLAGS                                   C
!  Purpose: Set flags for ghost cell flags                             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE SET_GHOST_CELL_FLAGS
    
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
      INTEGER :: IJK,I,J,K,I23,J23,K23
      INTEGER :: IP,IM,JP,JM,KP,KM
      INTEGER :: IPJK,IMJK,IJPK,IJMK,IJKP,IJKM
      INTEGER :: TOTAL_NUMBER_OF_INTERSECTIONS
      INTEGER :: NODE
      LOGICAL :: CLIP_FLAG
      include "function.inc"


!     EAST BOUNDARY
      I = IEND1
      
      IF(I==IMAX1) THEN

         DO I23 = IEND2,IEND3
            DO K = KSTART3,KEND3
               DO J = JSTART3, JEND3

                  IJK  = FUNIJK(I,J,K) 
                  IPJK = FUNIJK(I23,J,K) 

                  BLOCKED_CELL_AT(IPJK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IPJK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IPJK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IPJK) = BLOCKED_W_CELL_AT(IJK)
 
                  IF(BLOCKED_CELL_AT(IPJK)) FLAG(IPJK) = 100

                  VOL(IPJK)   = VOL(IJK)
                  VOL_U(IPJK) = VOL_U(IJK)
                  VOL_V(IPJK) = VOL_V(IJK)
                  VOL_W(IPJK) = VOL_W(IJK)

                  AYZ(IPJK)   = AYZ(IJK)
                  AYZ_U(IPJK) = AYZ_U(IJK) 
                  AYZ_V(IPJK) = AYZ_V(IJK)
                  AYZ_W(IPJK) = AYZ_W(IJK)

                  AXY(IPJK)   = AXY(IJK)
                  AXY_U(IPJK) = AXY_U(IJK) 
                  AXY_V(IPJK) = AXY_V(IJK)
                  AXY_W(IPJK) = AXY_W(IJK)

                  AXZ(IPJK)   = AXZ(IJK)
                  AXZ_U(IPJK) = AXZ_U(IJK) 
                  AXZ_V(IPJK) = AXZ_V(IJK)
                  AXZ_W(IPJK) = AXZ_W(IJK)

               END DO
            END DO
         END DO

      ENDIF

!     WEST BOUNDARY
      I = ISTART1

      IF(I==IMIN1) THEN

         DO I23 = ISTART3,ISTART2
            DO K = KSTART3,KEND3
               DO J = JSTART3, JEND3

                  IJK  = FUNIJK(I,J,K) 
                  IMJK = FUNIJK(I23,J,K) 

                  BLOCKED_CELL_AT(IMJK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IMJK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IMJK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IMJK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IMJK)) FLAG(IMJK) = 100

                  VOL(IMJK)   = VOL(IJK)
                  VOL_U(IMJK) = VOL_U(IJK)
                  VOL_V(IMJK) = VOL_V(IJK)
                  VOL_W(IMJK) = VOL_W(IJK)

                  AYZ(IMJK)   = AYZ(IJK)
                  AYZ_U(IMJK) = AYZ_U(IJK) 
                  AYZ_V(IMJK) = AYZ_V(IJK)
                  AYZ_W(IMJK) = AYZ_W(IJK)

                  AXY(IMJK)   = AXY(IJK)
                  AXY_U(IMJK) = AXY_U(IJK) 
                  AXY_V(IMJK) = AXY_V(IJK)
                  AXY_W(IMJK) = AXY_W(IJK)

                  AXZ(IMJK)   = AXZ(IJK)
                  AXZ_U(IMJK) = AXZ_U(IJK) 
                  AXZ_V(IMJK) = AXZ_V(IJK)
                  AXZ_W(IMJK) = AXZ_W(IJK)

               END DO
            END DO
         END DO

      ENDIF

!     NORTH BOUNDARY
      J = JEND1

      IF(J==JMAX1) THEN

         DO J23 = JEND2,JEND3
            DO K = KSTART3,KEND3
               DO I = ISTART3, IEND3

                  IJK  = FUNIJK(I,J,K) 
                  IJPK = FUNIJK(I,J23,K) 

                  BLOCKED_CELL_AT(IJPK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IJPK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IJPK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IJPK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IJPK)) FLAG(IJPK) = 100

                  VOL(IJPK)   = VOL(IJK)
                  VOL_U(IJPK) = VOL_U(IJK)
                  VOL_V(IJPK) = VOL_V(IJK)
                  VOL_W(IJPK) = VOL_W(IJK)

                  AYZ(IJPK)   = AYZ(IJK)
                  AYZ_U(IJPK) = AYZ_U(IJK) 
                  AYZ_V(IJPK) = AYZ_V(IJK)
                  AYZ_W(IJPK) = AYZ_W(IJK)

                  AXY(IJPK)   = AXY(IJK)
                  AXY_U(IJPK) = AXY_U(IJK) 
                  AXY_V(IJPK) = AXY_V(IJK)
                  AXY_W(IJPK) = AXY_W(IJK)

                  AXZ(IJPK)   = AXZ(IJK)
                  AXZ_U(IJPK) = AXZ_U(IJK) 
                  AXZ_V(IJPK) = AXZ_V(IJK)
                  AXZ_W(IJPK) = AXZ_W(IJK)

               END DO
            END DO
         END DO
 
      ENDIF
!     SOUTH BOUNDARY
      J = JSTART1

      IF(J==JMIN1) THEN

         DO J23 = JSTART3,JSTART2
            DO K = KSTART3,KEND3
               DO I = ISTART3, IEND3

                  IJK  = FUNIJK(I,J,K) 
                  IJMK = FUNIJK(I,J23,K) 

                  BLOCKED_CELL_AT(IJMK) = BLOCKED_CELL_AT(IJK)
                  BLOCKED_U_CELL_AT(IJMK) = BLOCKED_U_CELL_AT(IJK)
                  BLOCKED_V_CELL_AT(IJMK) = BLOCKED_V_CELL_AT(IJK)
                  BLOCKED_W_CELL_AT(IJMK) = BLOCKED_W_CELL_AT(IJK)

                  IF(BLOCKED_CELL_AT(IJMK)) FLAG(IJMK) = 100

                  VOL(IJMK)   = VOL(IJK)
                  VOL_U(IJMK) = VOL_U(IJK)
                  VOL_V(IJMK) = VOL_V(IJK)
                  VOL_W(IJMK) = VOL_W(IJK)

                  AYZ(IJMK)   = AYZ(IJK)
                  AYZ_U(IJMK) = AYZ_U(IJK) 
                  AYZ_V(IJMK) = AYZ_V(IJK)
                  AYZ_W(IJMK) = AYZ_W(IJK)

                  AXY(IJMK)   = AXY(IJK)
                  AXY_U(IJMK) = AXY_U(IJK) 
                  AXY_V(IJMK) = AXY_V(IJK)
                  AXY_W(IJMK) = AXY_W(IJK)

                  AXZ(IJMK)   = AXZ(IJK)
                  AXZ_U(IJMK) = AXZ_U(IJK) 
                  AXZ_V(IJMK) = AXZ_V(IJK)
                  AXZ_W(IJMK) = AXZ_W(IJK)

               END DO
            END DO
         END DO

      ENDIF

      IF(DO_K) THEN

!        TOP BOUNDARY
         K = KEND1

         IF(K==KMAX1) THEN
         
            DO K23=KEND2,KEND3
               DO J = JSTART3,JEND3
                  DO I = ISTART3, IEND3

                     IJK  = FUNIJK(I,J,K) 
                     IJKP = FUNIJK(I,J,K23) 

                     BLOCKED_CELL_AT(IJKP) = BLOCKED_CELL_AT(IJK)
                     BLOCKED_U_CELL_AT(IJKP) = BLOCKED_U_CELL_AT(IJK)
                     BLOCKED_V_CELL_AT(IJKP) = BLOCKED_V_CELL_AT(IJK)
                     BLOCKED_W_CELL_AT(IJKP) = BLOCKED_W_CELL_AT(IJK)

                     IF(BLOCKED_CELL_AT(IJKP)) FLAG(IJKP) = 100

                     VOL(IJKP)   = VOL(IJK)
                     VOL_U(IJKP) = VOL_U(IJK)
                     VOL_V(IJKP) = VOL_V(IJK)
                     VOL_W(IJKP) = VOL_W(IJK)

                     AYZ(IJKP)   = AYZ(IJK)
                     AYZ_U(IJKP) = AYZ_U(IJK) 
                     AYZ_V(IJKP) = AYZ_V(IJK)
                     AYZ_W(IJKP) = AYZ_W(IJK)

                     AXY(IJKP)   = AXY(IJK)
                     AXY_U(IJKP) = AXY_U(IJK) 
                     AXY_V(IJKP) = AXY_V(IJK)
                     AXY_W(IJKP) = AXY_W(IJK)

                     AXZ(IJKP)   = AXZ(IJK)
                     AXZ_U(IJKP) = AXZ_U(IJK) 
                     AXZ_V(IJKP) = AXZ_V(IJK)
                     AXZ_W(IJKP) = AXZ_W(IJK)

                  END DO
               END DO
            END DO

         ENDIF

!        BOTTOM BOUNDARY
         K = KSTART1

         IF(K==KMIN1) THEN

            DO K23 = KSTART3,KSTART2
               DO J = JSTART3,JEND3
                  DO I = ISTART3, IEND3

                     IJK  = FUNIJK(I,J,K) 
                     IJKM = FUNIJK(I,J,K23) 

                     BLOCKED_CELL_AT(IJKM) = BLOCKED_CELL_AT(IJK)
                     BLOCKED_U_CELL_AT(IJKM) = BLOCKED_U_CELL_AT(IJK)
                     BLOCKED_V_CELL_AT(IJKM) = BLOCKED_V_CELL_AT(IJK)
                     BLOCKED_W_CELL_AT(IJKM) = BLOCKED_W_CELL_AT(IJK)

                     IF(BLOCKED_CELL_AT(IJKM)) FLAG(IJKM) = 100

                     VOL(IJKM)   = VOL(IJK)
                     VOL_U(IJKM) = VOL_U(IJK)
                     VOL_V(IJKM) = VOL_V(IJK)
                     VOL_W(IJKM) = VOL_W(IJK)

                     AYZ(IJKM)   = AYZ(IJK)
                     AYZ_U(IJKM) = AYZ_U(IJK) 
                     AYZ_V(IJKM) = AYZ_V(IJK)
                     AYZ_W(IJKM) = AYZ_W(IJK)

                     AXY(IJKM)   = AXY(IJK)
                     AXY_U(IJKM) = AXY_U(IJK) 
                     AXY_V(IJKM) = AXY_V(IJK)
                     AXY_W(IJKM) = AXY_W(IJK)

                     AXZ(IJKM)   = AXZ(IJK)
                     AXZ_U(IJKM) = AXZ_U(IJK) 
                     AXZ_V(IJKM) = AXZ_V(IJK)
                     AXZ_W(IJKM) = AXZ_W(IJK)

                  END DO
               END DO
            END DO
         
         ENDIF   

      ENDIF

      RETURN

      
      END SUBROUTINE SET_GHOST_CELL_FLAGS

