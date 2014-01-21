!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CUT_CELL_PREPROCESSING                                 C
!  Purpose: Perform the cut-cell preprocessing stage:                  C
!           Identify cut cells, define face areas, and volumes         C
!           Set flags                                                  C
!           Compute Interpolations factors                             C
!           Compute Non-orthogonality Corrections terms                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CUT_CELL_PREPROCESSING
    
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
      USE vtk
      use cdist
      
      USE discretelement
      USE des_stl_functions
      USE fldvar

      IMPLICIT NONE

      INTEGER :: IJK
      INTEGER :: I,J,K, COUNT
      INTEGER :: SAFE_MODE_COUNT
      DOUBLE PRECISION :: CPU_PP_START,CPU_PP_END

      include "function.inc"   

      IF(myPE == PE_IO) THEN

         WRITE(*,10)'============================================================================='
         WRITE(*,10)'  ____          ___________    _   ____      ____    __________   __________  '
         WRITE(*,10)' |    \        /           |  (_)  \   \    /   /   |          | |          | '
         WRITE(*,10)' |     \      /      ______|  ___   \   \  /   /    |    ______| |    ______| '
         WRITE(*,10)' |      \    /      |______  |   |   \   \/   /     |   |        |   |        '   
         WRITE(*,10)' |       \  /              | |   |    \      /  === |   |        |   |  ____  '
         WRITE(*,10)' |   |\   \/   /|    ______| |   |    /      \      |   |        |   | |_   | '
         WRITE(*,10)' |   | \      / |   |        |   |   /   /\   \     |   |______  |   |___|  | '
         WRITE(*,10)' |   |  \    /  |   |        |   |  /   /  \   \    |          | |          | '
         WRITE(*,10)' |___|   \__/   |___|        |___| /___/    \___\   |__________| |__________| '
         WRITE(*,10)'                                                                              '
         WRITE(*,10)'============================================================================='
         WRITE(*,10)'MFIX WITH CARTESIAN GRID IMPLEMENTATION.'
         WRITE(*,10)'STARTING CARTESIAN GRID PRE-PROCESSING...'

      ENDIF

      CALL CPU_TIME (CPU_PP_START)



      CALL OPEN_CUT_CELL_FILES

      CALL ALLOCATE_CUT_CELL_ARRAYS 

      CALL DEFINE_QUADRICS

      CALL SET_3D_CUT_CELL_FLAGS 

      CALL GATHER_DATA

!======================================================================
! Gather Data  and writes surface(s) defined by all cut cells 
!======================================================================

      IF(WRITE_VTK_FILES.AND.(.NOT.BDIST_IO)) THEN
         CALL WRITE_CUT_SURFACE_VTK
      ENDIF



      CALL SET_3D_CUT_U_CELL_FLAGS 
      CALL SET_3D_CUT_V_CELL_FLAGS 
      IF(DO_K) CALL SET_3D_CUT_W_CELL_FLAGS 

      CALL SET_3D_CUT_CELL_TREATMENT_FLAGS

      CALL GET_3D_ALPHA_U_CUT_CELL 
      CALL GET_3D_ALPHA_V_CUT_CELL 
      IF(DO_K) CALL GET_3D_ALPHA_W_CUT_CELL 

      CALL SET_GHOST_CELL_FLAGS

      CALL SET_ODXYZ_U_CUT_CELL
      CALL SET_ODXYZ_V_CUT_CELL
      IF(DO_K) CALL SET_ODXYZ_W_CUT_CELL

      CALL GET_U_MASTER_CELLS
      CALL GET_V_MASTER_CELLS
      IF(DO_K) CALL GET_W_MASTER_CELLS

      CALL SEND_RECEIVE_CUT_CELL_VARIABLES

      CALL PRINT_GRID_STATISTICS

      CALL CG_GET_BC_AREA   

      CALL FLOW_TO_VEL(.FALSE.)

      CALL CG_FLOW_TO_VEL

      CALL CONVERT_CG_MI_TO_PS
 
      CALL GET_DISTANCE_TO_WALL

      IF(DISCRETE_ELEMENT.and.USE_STL) then
         IF(myPE == PE_IO) THEN
            WRITE(*,*)'DISCRETE ELEMENT also detected. Pre-Processing for des now'
         endif

         !now call the pre-procssing for the des in order to 
         !assign facets to grid cells 
         !Set N_facets_des to add any more facets needed by
         !dem and not to contaminate the Eulerian-Eulerian CG stuff 
         N_FACETS_DES = N_FACETS 

         !**********************************************************************************
         !this is a temporary routine to triangulate the 
         !default walls (bounding box of the simulation)
         !It will be users' burden to supply the bounding box 
         !as stl         
         if(DES_CONVERT_BOX_TO_FACETS) call cg_des_convert_to_facets
         !*********************************************************************************
         
         CALL bin_facets_to_grid_des
         
         DO IJK = 1, DIMENSION_3 
            COUNT = LIST_FACET_AT_DES(IJK)%COUNT_FACETS 
            IF(COUNT.eq.0) DEALLOCATE(LIST_FACET_AT_DES(IJK)%FACET_LIST)
         ENDDO

         !CALL DEBUG_WRITE_GRID_FACEINFO
         !CALL DEBUG_write_stl_from_grid_facet(WRITE_FACETS_EACH_CELL=.false.)

         IF(myPE == PE_IO) THEN
            WRITE(*,*)'Done with post processing for DISCRETE ELEMENT'
         endif

      endif  !discrete_element  

      CALL CPU_TIME (CPU_PP_END)

      
      
      IF(myPE == PE_IO) THEN
         WRITE(*,20)'CARTESIAN GRID PRE-PROCESSING COMPLETED IN ',CPU_PP_END - CPU_PP_START, ' SECONDS.'
         WRITE(*,10)'============================================================================='
      ENDIF

      IF(myPE == PE_IO) THEN

         SAFE_MODE_COUNT = SUM(CG_SAFE_MODE)

         IF(SAFE_MODE_COUNT>0) THEN


            WRITE(*,10)'######################################################################'
            WRITE(*,10)'######################################################################'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                              \/                                  ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##  ===>   WARNING: RUNNING CARTESIAN GRID IN SAFE MODE !  <===     ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##  SAFE MODE ACTIVATED FOR :                                       ##'
            IF(CG_SAFE_MODE(1)==1) WRITE(*,10)'##                            - All scalar quantities               ##'
            IF(CG_SAFE_MODE(3)==1) WRITE(*,10)'##                            - X-Velocity (Gas and Solids)         ##'
            IF(CG_SAFE_MODE(4)==1) WRITE(*,10)'##                            - Y-Velocity (Gas and Solids)         ##'
            IF(CG_SAFE_MODE(5)==1) WRITE(*,10)'##                            - Z-Velocity (Gas and Solids)         ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'##                              /\                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                              ||                                  ##'
            WRITE(*,10)'##                                                                  ##'
            WRITE(*,10)'######################################################################'
            WRITE(*,10)'######################################################################'



         ENDIF
      ENDIF


      RETURN

10    FORMAT(1X,A)
20    FORMAT(1X,A,F8.2,A)
      
      END SUBROUTINE CUT_CELL_PREPROCESSING


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EVAL_F                                                 C
!  Purpose: Evaluate the function f(x,y,z) defining the boundary       C
!                                                                      C
!  Author: Jeff Dietiker                               Date: 21-FEB-08 C
!  Reviewer:                                           Date:           C
!                                                                      C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!  Modified: ##                                        Date: ##-###-## C
!  Purpose: ##                                                         C
!                                                                      C
!  Literature/Document References: ##                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
        SUBROUTINE EVAL_F(METHOD,x1,x2,x3,Q,f,CLIP_FLAG)
 
      USE parallel
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell

      USE quadric
  

      IMPLICIT NONE
 
      DOUBLE PRECISION x1,x2,x3
      DOUBLE PRECISION f
      DOUBLE PRECISION, DIMENSION(1,3) :: X_VECTOR,XMT
      DOUBLE PRECISION, DIMENSION(3,1) :: TXMT
      DOUBLE PRECISION, DIMENSION(1,1) :: TEMP_1x1
      INTEGER :: Q,Q_ID,BCID
      LOGICAL :: CLIP_X,CLIP_Y,CLIP_Z,CLIP_FLAG
      CHARACTER (LEN = 7) :: METHOD
      CHARACTER(LEN=9) :: GR

      INTEGER :: GROUP,GS,P

      DOUBLE PRECISION,DIMENSION(DIM_GROUP,0:DIM_QUADRIC) :: F_G

      SELECT CASE(METHOD)

         CASE('QUADRIC')

            F_G = - UNDEFINED

            DO GROUP = 1, N_GROUP
               GS = GROUP_SIZE(GROUP)
               GR = TRIM(GROUP_RELATION(GROUP)) 

               DO P = 1 , GS
                  Q_ID = GROUP_Q(GROUP,P)
                  CALL GET_F_QUADRIC(x1,x2,x3,Q_ID,F_G(GROUP,P),CLIP_FLAG)
               ENDDO
               IF(GR == 'AND') THEN
                  F_G(GROUP,0) = MAXVAL(F_G(GROUP,1:GS))
               ELSEIF(GR == 'OR') THEN
                  F_G(GROUP,0) = MINVAL(F_G(GROUP,1:GS))
               ELSEIF(GR == 'PIECEWISE') THEN
                  CALL REASSSIGN_QUADRIC(x1,x2,x3,GROUP,Q_ID)
!                  CLIP_FLAG=.FALSE.
                  CALL GET_F_QUADRIC(x1,x2,x3,Q_ID,F_G(GROUP,0),CLIP_FLAG)
!                  CLIP_FLAG=.TRUE.
               ENDIF

            ENDDO

            f = F_G(1,0)

            DO GROUP = 2, N_GROUP

               GR = TRIM(RELATION_WITH_PREVIOUS(GROUP)) 

               IF(GR =='AND') THEN
                  f = DMAX1(f,F_G(GROUP,0))
               ELSEIF(GR =='OR') THEN
                  f = DMIN1(f,F_G(GROUP,0))
               ENDIF

            ENDDO


         CASE('POLYGON')

            CALL EVAL_POLY_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)

         CASE('USR_DEF')

            CALL EVAL_USR_FCT(x1,x2,x3,Q,f,CLIP_FLAG)

!         CASE('STL')

!            CALL EVAL_STL_FCT(x1,x2,x3,Q,f,CLIP_FLAG,BCID)

         CASE DEFAULT

            WRITE(*,*)'ERROR IN SUBROUTINE EVAL_F.'
            WRITE(*,*)'UNKNOWN METHOD:',METHOD
            WRITE(*,*)'ACCEPTABLE METHODS:'
            WRITE(*,*)'QUADRIC'
            WRITE(*,*)'POLYGON'
            WRITE(*,*)'USR_DEF'
!            WRITE(*,*)'STL'
            CALL MFIX_EXIT(myPE)

      END SELECT



      RETURN
      END SUBROUTINE EVAL_F     



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: intersect_line                                         C
!  Purpose: Finds the intersection between the quadric surface ,       C
!           and the line (xa,ya,za) and (xb,yb,zb).                    C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE INTERSECT_LINE(METHOD,xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
    
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
      
      IMPLICIT NONE
      DOUBLE PRECISION:: x1,y1,z1,x2,y2,z2,x3,y3,z3
      DOUBLE PRECISION:: xa,ya,za,xb,yb,zb,xc,yc,zc
      INTEGER :: Q_ID,niter
      DOUBLE PRECISION:: x1save,x2save
      DOUBLE PRECISION :: y_plane,z_plane
      DOUBLE PRECISION :: x_intersection
      DOUBLE PRECISION :: f1,f2,f3,fa,fb
      DOUBLE PRECISION :: t1,t2,t3
      DOUBLE PRECISION :: xac,xbc,yac,ybc
      LOGICAL :: CLIP_FLAG,CLIP_FLAG1,CLIP_FLAG2,CLIP_FLAG3,INTERSECT_FLAG
      CHARACTER (LEN=7) ::METHOD

      x1 = xa    ! Initial guesses
      y1 = ya
      z1 = za
      t1 = ZERO

      x2 = xb
      y2 = yb
      z2 = zb
      t2 = ONE

      CALL EVAL_F(METHOD,x1,y1,z1,Q_ID,f1,CLIP_FLAG1)
      CALL EVAL_F(METHOD,x2,y2,z2,Q_ID,f2,CLIP_FLAG2)

!======================================================================
!  The line from (x1,y1,z1) and (x2,y2,z2) is parametrized 
!  from t1 = ZERO to t2 = ONE
!======================================================================

      niter = 1

      CLIP_FLAG = (CLIP_FLAG1).AND.(CLIP_FLAG2)

      if(DABS(f1)<TOL_F) then  ! ignore intersection at corner
        xc = UNDEFINED
        yc = UNDEFINED
        zc = UNDEFINED
        INTERSECT_FLAG = .FALSE.
      elseif(DABS(f2)<TOL_F) then ! ignore intersection at corner
        xc = UNDEFINED
        yc = UNDEFINED
        zc = UNDEFINED
        INTERSECT_FLAG = .FALSE.
      elseif(f1*f2 < ZERO) then
        niter = 0
        f3 = 2.0d0*TOL_F
        do while (   (abs(f3) > TOL_F)   .AND.   (niter<ITERMAX_INT)       )
         
          t3 = t1 - f1*(t2-t1)/(f2-f1)

          x3 = x1 + t3 * (x2 - x1)
          y3 = y1 + t3 * (y2 - y1)
          z3 = z1 + t3 * (z2 - z1)

          CALL EVAL_F(METHOD,x3,y3,z3,Q_ID,f3,CLIP_FLAG3)
          if(f1*f3<0) then
            t2 = t3
            f2 = f3
          else 
            t1 = t3
            f1 = f3
          endif   
          niter = niter + 1

        end do
        if (niter < ITERMAX_INT) then
           xc = x3
           yc = y3
           zc = z3
          INTERSECT_FLAG = .TRUE.
        else
           WRITE(*,*)'   Subroutine intersect_line:'
           WRITE(*,*)   'Unable to find the intersection of quadric:',Q_ID
           WRITE(*,1000)'between (x1,y1,z1)= ', xa,ya,za
           WRITE(*,1000)'   and  (x2,y2,z2)= ', xb,yb,zb
           CALL EVAL_F(METHOD,xa,ya,za,Q_ID,fa,CLIP_FLAG1)
           CALL EVAL_F(METHOD,xb,yb,zb,Q_ID,fb,CLIP_FLAG1)
           WRITE(*,1000)'f(x1,y1,z1) = ', fa
           WRITE(*,1000)'f(x2,y2,z2) = ', fb
           WRITE(*,1000)'Current Location (x3,y3,z3)= ', x3,y3,z3
           WRITE(*,1000)'Current value of abs(f) = ', DABS(f3)
           WRITE(*,1000)'Tolerance = ', TOL_F
           WRITE(*,*)'Maximum number of iterations = ', ITERMAX_INT
           WRITE(*,*)   'Please increase the intersection tolerance, '
           WRITE(*,*)   'or the maximum number of iterations, and try again.'
           WRITE(*,*)   'MFiX will exit now.'             
           CALL MFIX_EXIT(myPE) 
           x_intersection = UNDEFINED
           INTERSECT_FLAG = .FALSE.

        endif
      else
        xc = UNDEFINED
        yc = UNDEFINED
        zc = UNDEFINED
        INTERSECT_FLAG = .FALSE.
      endif

 1000 FORMAT(A,3(2X,G12.5)) 


      RETURN

      END SUBROUTINE INTERSECT_LINE

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: INTERSECT                                              C
!  Purpose: Intersects quadric with grid                               C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE INTERSECT(IJK,TYPE_OF_CELL,Xi,Yi,Zi)
    
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
      
      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,Q_ID,Q_ID2,N_int_x,N_int_y,N_int_z,N_USR
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb,xc,yc,zc,Fc
      DOUBLE PRECISION :: Xi,Yi,Zi,Xc_backup,Yc_backup,Zc_backup
      LOGICAL :: INTERSECT_FLAG,CLIP_FLAG

      include "function.inc"

      Xi = UNDEFINED
      Yi = UNDEFINED
      Zi = UNDEFINED

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

!======================================================================
!  Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================
      xa = X_NODE(7)
      ya = Y_NODE(7)
      za = Z_NODE(7)

      xb = X_NODE(8)
      yb = Y_NODE(8)
      zb = Z_NODE(8)


      N_int_x = 0
      INTERSECT_X(IJK) = .FALSE.
      Q_ID = 1
         CALL INTERSECT_LINE('QUADRIC',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
         IF ((INTERSECT_FLAG).AND.(xc/=Xi)) THEN
            N_int_x = N_int_x + 1
            INTERSECT_X(IJK) = .TRUE.
            xc_backup = Xi
            Xi = xc
         ENDIF

      IF(N_int_x /= 1) THEN
         Xi = UNDEFINED
         INTERSECT_X(IJK) = .FALSE.
      ENDIF


      DO Q_ID = 1, N_POLYGON
         CALL INTERSECT_LINE('POLYGON',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_X(IJK),xc,yc,zc)
         IF(INTERSECT_X(IJK)) Xi = xc
      ENDDO

      DO N_USR= 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_X(IJK),xc,yc,zc)
         IF(INTERSECT_X(IJK)) Xi = xc
      ENDDO

!      IF(USE_STL) THEN
!         CALL INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_X(IJK),xc,yc,zc)
!         IF(INTERSECT_X(IJK)) Xi = xc 
!      ENDIF

      IF(TYPE_OF_CELL=='U_MOMENTUM') THEN
         IF(SNAP(IJK)) THEN
            INTERSECT_X(IJK) = .TRUE.
            I = I_OF(IJK) 
            Xi = XG_E(I)
         ENDIF
      ENDIF



!======================================================================
!  Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================
      xa = X_NODE(6)
      ya = Y_NODE(6)
      za = Z_NODE(6)

      N_int_y = 0
      INTERSECT_Y(IJK) = .FALSE.
      Q_ID = 1
         CALL INTERSECT_LINE('QUADRIC',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
         IF ((INTERSECT_FLAG).AND.(yc/=Yi)) THEN
            N_int_y = N_int_y + 1
            INTERSECT_Y(IJK) = .TRUE.
            yc_backup = Yi
            Yi = yc            
         ENDIF

      IF(N_int_y /= 1) THEN
         Yi = UNDEFINED
         INTERSECT_Y(IJK) = .FALSE.
      ENDIF

      DO Q_ID = 1, N_POLYGON
         CALL INTERSECT_LINE('POLYGON',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Y(IJK),xc,yc,zc)
         IF(INTERSECT_Y(IJK)) Yi = yc
      ENDDO

      DO N_USR= 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Y(IJK),xc,yc,zc)
         IF(INTERSECT_Y(IJK)) Yi = yc
      ENDDO

!      IF(USE_STL) THEN
!         CALL INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_Y(IJK),xc,yc,zc)
!         IF(INTERSECT_Y(IJK)) Yi = yc 
!      ENDIF

      IF(TYPE_OF_CELL=='V_MOMENTUM') THEN
         IF(SNAP(IJK)) THEN
            INTERSECT_Y(IJK) = .TRUE.
            J = J_OF(IJK) 
            Yi = YG_N(J)
         ENDIF
      ENDIF

      IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================
         xa = X_NODE(4)
         ya = Y_NODE(4)
         za = Z_NODE(4)

         N_int_z = 0
         INTERSECT_Z(IJK) = .FALSE.
         Q_ID = 1
            CALL INTERSECT_LINE('QUADRIC',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_FLAG,xc,yc,zc)
            IF ((INTERSECT_FLAG).AND.(zc/=Zi)) THEN
               N_int_z = N_int_z + 1
               INTERSECT_Z(IJK) = .TRUE.
               zc_backup = Zi
               Zi = zc            
            ENDIF

            IF(N_int_z /= 1) THEN
               Zi = UNDEFINED
               INTERSECT_Z(IJK) = .FALSE.
            ENDIF

         DO Q_ID = 1, N_POLYGON
            CALL INTERSECT_LINE('POLYGON',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Z(IJK),xc,yc,zc)
            IF(INTERSECT_Z(IJK)) Zi = zc
         ENDDO

         DO N_USR= 1, N_USR_DEF
            CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Z(IJK),xc,yc,zc)
            IF(INTERSECT_Z(IJK)) Zi = zc
         ENDDO

   !      IF(USE_STL) THEN
   !         CALL INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_Z(IJK),xc,yc,zc)
   !         IF(INTERSECT_Z(IJK)) Zi = zc 
   !      ENDIF

         IF(TYPE_OF_CELL=='W_MOMENTUM') THEN
            IF(SNAP(IJK)) THEN
               INTERSECT_Z(IJK) = .TRUE.
               K = K_OF(IJK) 
               Zi = ZG_T(K)
            ENDIF
         ENDIF

      ENDIF


      IF(INTERSECT_X(IJK)) THEN
         POTENTIAL_CUT_CELL_AT(IJK) = .TRUE.
         POTENTIAL_CUT_CELL_AT(EAST_OF(IJK)) = .TRUE.
      ENDIF

      IF(INTERSECT_Y(IJK)) THEN
         POTENTIAL_CUT_CELL_AT(IJK) = .TRUE.
         POTENTIAL_CUT_CELL_AT(NORTH_OF(IJK)) = .TRUE.
      ENDIF


      IF(INTERSECT_Z(IJK)) THEN
         POTENTIAL_CUT_CELL_AT(IJK) = .TRUE.
         POTENTIAL_CUT_CELL_AT(TOP_OF(IJK)) = .TRUE.
      ENDIF


      RETURN
      
      END SUBROUTINE INTERSECT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLEAN_INTERSECT                                        C
!  Purpose: Remove Intersection flags in preparation of small cell     C
!           removal                                                    C 
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLEAN_INTERSECT(IJK,TYPE_OF_CELL,Xi,Yi,Zi)
    
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
      USE STL
      
      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,IM,JM,KM,IP,JP,KP
      INTEGER :: BCID
      INTEGER :: IMJK,IPJK,IJMK,IJPK,IJKM,IJKP,IMJPK,IMJKP,IPJMK,IJMKP,IPJKM,IJPKM
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb,xc,yc,zc
      DOUBLE PRECISION :: Xi,Yi,Zi
      DOUBLE PRECISION :: DFC,DFC_MAX,Fa,Fb,F4,F6,F7,F8
      LOGICAL :: CLIP_FLAG

      include "function.inc"

!======================================================================
!  Get coordinates of eight nodes
!======================================================================

      CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

      I = I_OF(IJK) 
      J = J_OF(IJK) 
      K = K_OF(IJK) 

      IM = I - 1 
      JM = J - 1 
      KM = K - 1

      IP = I + 1 
      JP = J + 1 
      KP = K + 1 

      IMJK = FUNIJK(IM,J,K)
      IPJK = FUNIJK(IP,J,K)
      IJMK = FUNIJK(I,JM,K)
      IJPK = FUNIJK(I,JP,K)
      IJKM = FUNIJK(I,J,KM)
      IJKP = FUNIJK(I,J,KP)

      IMJPK = FUNIJK(IM,JP,K)
      IMJKP = FUNIJK(IM,J,KP)

      IPJMK = FUNIJK(IP,JM,K)
      IJMKP = FUNIJK(I,JM,KP)

      IPJKM = FUNIJK(IP,J,KM)
      IJPKM = FUNIJK(I,JP,KM)


      IF(IMJK<1.OR.IMJK>DIMENSION_3) IMJK = IJK
      IF(IPJK<1.OR.IPJK>DIMENSION_3) IPJK = IJK
      IF(IJMK<1.OR.IJMK>DIMENSION_3) IJMK = IJK
      IF(IJPK<1.OR.IJPK>DIMENSION_3) IJPK = IJK
      IF(IJKM<1.OR.IJKM>DIMENSION_3) IJKM = IJK
      IF(IJKP<1.OR.IJKP>DIMENSION_3) IJKP = IJK

      IF(IMJPK<1.OR.IMJPK>DIMENSION_3) IMJPK = IJK
      IF(IMJKP<1.OR.IMJKP>DIMENSION_3) IMJKP = IJK

      IF(IPJMK<1.OR.IPJMK>DIMENSION_3) IPJMK = IJK
      IF(IJMKP<1.OR.IJMKP>DIMENSION_3) IJMKP = IJK

      IF(IPJKM<1.OR.IPJKM>DIMENSION_3) IPJKM = IJK
      IF(IJPKM<1.OR.IJPKM>DIMENSION_3) IJPKM = IJK

!======================================================================
!  Clean Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================

      xa = X_NODE(7)
      ya = Y_NODE(7)
      za = Z_NODE(7)

      xb = X_NODE(8)
      yb = Y_NODE(8)
      zb = Z_NODE(8)

      DFC_MAX = TOL_SNAP(1) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER

      IF(INTERSECT_X(IJK)) THEN

         DFC = DABS(Xi-xa) ! DISTANCE FROM CORNER (NODE 7)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 7'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            INTERSECT_X(IJK)  = .FALSE.
            IF(I>=IMIN1) THEN
               INTERSECT_X(IMJK) = .FALSE.
               INTERSECT_Y(IMJK)  = .FALSE.
               INTERSECT_Y(IMJPK) = .FALSE.
               IF(DO_K) INTERSECT_Z(IMJK)  = .FALSE.
               IF(DO_K) INTERSECT_Z(IMJKP) = .FALSE.

               SNAP(IMJK) = .TRUE.
            ENDIF
                
         ENDIF


         DFC = DABS(Xi-xb) ! DISTANCE FROM CORNER (NODE 8)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF


            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            IF(DO_K) INTERSECT_Z(IJK)  = .FALSE.
            SNAP(IJK) = .TRUE.

            IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
            IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
            IF(DO_K.AND.(K<=KMAX1)) INTERSECT_Z(IJKP) = .FALSE.

                
         ENDIF

         IF(USE_STL.OR.USE_MSH) THEN
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,7,F7,CLIP_FLAG,BCID)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
            IF(F7*F8>TOL_STL**2) INTERSECT_X(IJK)  = .FALSE.
         ENDIF

      ENDIF


      

!======================================================================
!  Clean Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================

      xa = X_NODE(6)
      ya = Y_NODE(6)
      za = Z_NODE(6)

      DFC_MAX = TOL_SNAP(2) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER


      IF(INTERSECT_Y(IJK)) THEN

         DFC = DABS(Yi-ya) ! DISTANCE FROM CORNER (NODE 6)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 6'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            INTERSECT_Y(IJK)  = .FALSE.
            IF(J>=JMIN1) THEN
               INTERSECT_X(IJMK)  = .FALSE.
               INTERSECT_X(IPJMK) = .FALSE.
               INTERSECT_Y(IJMK) = .FALSE.
               IF(DO_K) INTERSECT_Z(IJMK)  = .FALSE.
               IF(DO_K) INTERSECT_Z(IJMKP) = .FALSE.

               SNAP(IJMK) = .TRUE.
            ENDIF
                
         ENDIF


         DFC = DABS(Yi-yb) ! DISTANCE FROM CORNER (NODE 8)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF 

            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            IF(DO_K) INTERSECT_Z(IJK)  = .FALSE.
            SNAP(IJK) = .TRUE.

            IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
            IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
            IF(DO_K.AND.(K<=KMAX1)) INTERSECT_Z(IJKP) = .FALSE.

                
         ENDIF


         IF(USE_STL.OR.USE_MSH) THEN
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,6,F6,CLIP_FLAG,BCID)
            CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
            IF(F6*F8>TOL_STL**2) INTERSECT_Y(IJK)  = .FALSE.
         ENDIF

      ENDIF


      IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================

         xa = X_NODE(4)
         ya = Y_NODE(4)
         za = Z_NODE(4)

         DFC_MAX = TOL_SNAP(3) * DSQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)  ! MAXIMUM DISTANCE FROM CORNER

         IF(INTERSECT_Z(IJK)) THEN

            DFC = DABS(Zi-Za) ! DISTANCE FROM CORNER (NODE 4)

            IF(DFC < DFC_MAX) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 4'
                  WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
               ENDIF 

               INTERSECT_Z(IJK)  = .FALSE.

               IF(K>=KMIN1) THEN
                  INTERSECT_X(IJKM)  = .FALSE.
                  INTERSECT_X(IPJKM) = .FALSE.
                  INTERSECT_Y(IJKM)  = .FALSE.
                  INTERSECT_Y(IJPKM) = .FALSE.
                  INTERSECT_Z(IJKM) = .FALSE.
                   
                  SNAP(IJKM) = .TRUE.
               ENDIF

            ENDIF


            DFC = DABS(Zi-Zb) ! DISTANCE FROM CORNER (NODE 8)

            IF(DFC < DFC_MAX) THEN
               IF(PRINT_WARNINGS) THEN
                  WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 8'
                  WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
               ENDIF

               INTERSECT_X(IJK)  = .FALSE.
               INTERSECT_Y(IJK)  = .FALSE.
               INTERSECT_Z(IJK)  = .FALSE.
               SNAP(IJK) = .TRUE.

               IF(I<=IMAX1) INTERSECT_X(IPJK) = .FALSE.
               IF(J<=JMAX1) INTERSECT_Y(IJPK) = .FALSE.
               IF(K<=KMAX1) INTERSECT_Z(IJKP) = .FALSE.

                   
            ENDIF

            IF(USE_STL.OR.USE_MSH) THEN
               CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,4,F4,CLIP_FLAG,BCID)
               CALL EVAL_STL_FCT_AT(TYPE_OF_CELL,IJK,8,F8,CLIP_FLAG,BCID)
               IF(F4*F8>TOL_STL**2) INTERSECT_Z(IJK)  = .FALSE.
            ENDIF

         ENDIF

      ENDIF


      RETURN
      
      END SUBROUTINE CLEAN_INTERSECT



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: OPEN_CUT_CELL_FILES                                    C
!  Purpose: Open CUT CELL related file                                 C
!           and writes headers                                         C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE OPEN_CUT_CELL_FILES
    
      USE cutcell
      USE compar
      
      IMPLICIT NONE

      IF(MyPE == PE_IO)  THEN
         OPEN(UNIT = UNIT_CUT_CELL_LOG, FILE = 'CUT_CELL.LOG')
      ENDIF

      RETURN

      
      END SUBROUTINE OPEN_CUT_CELL_FILES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CLOSE_CUT_CELL_FILES                                   C
!  Purpose: Close CUT CELL related file                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CLOSE_CUT_CELL_FILES
    
      USE compar
      USE cutcell
      
      IMPLICIT NONE

      IF(MyPE == PE_IO) THEN
         CLOSE(UNIT_CUT_CELL_LOG)
      ENDIF

      RETURN

      
      END SUBROUTINE CLOSE_CUT_CELL_FILES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CAD_INTERSECT                                          C
!  Purpose: Intersects CAD (STL file or MSH file) geometry with grid   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CAD_INTERSECT(TYPE_OF_CELL,Xint,Yint,Zint)
    
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


      USE mpi_utility



      
      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,Q_ID,Q_ID2,N_int_x,N_int_y,N_int_z, COUNT
      INTEGER :: IM,IP,JM,JP,KM,KP,IMJK,IPJK,IJMK,IJPK,IJKM,IJKP
      INTEGER :: IJPKP,IPJKP,IPJPK
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb,xc,yc,zc,Fc
      DOUBLE PRECISION :: Xi,Yi,Zi,Xc_backup,Yc_backup,Zc_backup
      LOGICAL :: INTERSECT_FLAG,CLIP_FLAG,INSIDE_FACET

      DOUBLE PRECISION :: X1,X2,Y1,Y2,Z1,Z2

      DOUBLE PRECISION, DIMENSION(DIMENSION_3) :: Xint,Yint,Zint 

      INTEGER :: N,I1,I2,J1,J2,K1,K2

      DOUBLE PRECISION :: X_OFFSET, Y_OFFSET, Z_OFFSET

      DOUBLE PRECISION, DIMENSION(3) :: N4,N6,N7,N8
      DOUBLE PRECISION :: CURRENT_F,dotproduct

      INTEGER :: IJK2,CURRENT_I,CURRENT_J,CURRENT_K

      INTEGER :: N_UNDEFINED, N_PROP
      INTEGER, PARAMETER :: N_PROPMAX=1000
      LOGICAL:: F_FOUND


!      CHARACTER (LEN=3) :: CAD_PROPAGATE_ORDER

      include "function.inc"


      INTERSECT_X = .FALSE.
      INTERSECT_Y = .FALSE.
      INTERSECT_Z = .FALSE.

      Xint = UNDEFINED
      Yint = UNDEFINED
      Zint = UNDEFINED

      F_AT = UNDEFINED


      SELECT CASE (TYPE_OF_CELL)
         CASE('SCALAR')


            X_OFFSET = ZERO
            Y_OFFSET = ZERO
            Z_OFFSET = ZERO

         CASE('U_MOMENTUM')

            X_OFFSET = HALF
            Y_OFFSET = ZERO
            Z_OFFSET = ZERO

         CASE('V_MOMENTUM')

            X_OFFSET = ZERO
            Y_OFFSET = HALF
            Z_OFFSET = ZERO

         CASE('W_MOMENTUM')

            X_OFFSET = ZERO
            Y_OFFSET = ZERO
            Z_OFFSET = HALF


         CASE DEFAULT
            WRITE(*,*)'SUBROUTINE: GET_CELL_NODE_COORDINATES'
            WRITE(*,*)'UNKNOWN TYPE OF CELL:',TYPE_OF_CELL
            WRITE(*,*)'ACCEPTABLE TYPES ARE:' 
            WRITE(*,*)'SCALAR' 
            WRITE(*,*)'U_MOMENTUM' 
            WRITE(*,*)'V_MOMENTUM' 
            WRITE(*,*)'W_MOMENTUM' 
            CALL MFIX_EXIT(myPE)
      END SELECT


      DO N = 1,N_FACETS


         X1 = MINVAL(VERTEX(N,1:3,1))
         X2 = MAXVAL(VERTEX(N,1:3,1))
         Y1 = MINVAL(VERTEX(N,1:3,2))
         Y2 = MAXVAL(VERTEX(N,1:3,2))
         Z1 = MINVAL(VERTEX(N,1:3,3))
         Z2 = MAXVAL(VERTEX(N,1:3,3))


         I1 = IEND3
         I2 = ISTART3

         IF(X2>=ZERO.AND.X1<=XLENGTH) THEN
            DO I = ISTART3, IEND3
               IP = I+1
               IF(XG_E(I)+X_OFFSET*DX(IP)>=X1-TOL_STL) THEN
                  I1=I
                  EXIT
               ENDIF
            ENDDO

            DO I = IEND3, ISTART3,-1
               IP = I+1
               IF(XG_E(I)-DX(I)+X_OFFSET*DX(IP)<=X2+TOL_STL) THEN
                  I2=I
                  EXIT
               ENDIF
            ENDDO
         ENDIF


         J1 = JEND3
         J2 = JSTART3

         IF(Y2>=ZERO.AND.Y1<=YLENGTH) THEN
            DO J = JSTART3, JEND3
               JP = J+1
               IF(YG_N(J)+Y_OFFSET*DY(JP)>=Y1-TOL_STL) THEN
                  J1=J
                  EXIT
               ENDIF
            ENDDO

            DO J = JEND3, JSTART3,-1
               JP=J+1
               IF(YG_N(J)-DY(J)+Y_OFFSET*DY(JP)<=Y2+TOL_STL) THEN
                  J2=J
                  EXIT
               ENDIF
            ENDDO
         ENDIF

         K1 = KEND3
         K2 = KSTART3

         IF(Z2>=ZERO.AND.Z1<=ZLENGTH) THEN
            DO K = KSTART3, KEND3
               KP=K+1

               IF(ZG_T(K)+Z_OFFSET*DZ(KP)>=Z1-TOL_STL) THEN
                  K1=K
                  EXIT
               ENDIF
            ENDDO

            DO K = KEND3, KSTART3,-1
               KP = K+1
               IF(ZG_T(K)-DZ(K)+Z_OFFSET*DZ(KP)<=Z2+TOL_STL) THEN
                  K2=K
                  EXIT
               ENDIF
            ENDDO
         ENDIF




         DO K=K1,K2
            DO J=J1,J2
               DO I=I1,I2

                  IJK = FUNIJK(I,J,K)


                  IM = MAX0(I - 1 ,ISTART3)
                  JM = MAX0(J - 1 ,JSTART3)
                  KM = MAX0(K - 1 ,KSTART3)

                  IP = MIN0(I + 1 ,IEND3)
                  JP = MIN0(J + 1 ,JEND3)
                  KP = MIN0(K + 1 ,KEND3)


                  IMJK = FUNIJK(IM,J,K)
                  IPJK = FUNIJK(IP,J,K)
                  IJMK = FUNIJK(I,JM,K)
                  IJPK = FUNIJK(I,JP,K)
                  IJKM = FUNIJK(I,J,KM)
                  IJKP = FUNIJK(I,J,KP)

                  IJPKP = FUNIJK(I,JP,KP)
                  IPJKP = FUNIJK(IP,J,KP)
                  IPJPK = FUNIJK(IP,JP,K)
                  
                  
!======================================================================
!  Get coordinates of eight nodes
!======================================================================

                  CALL GET_CELL_NODE_COORDINATES(IJK,TYPE_OF_CELL)

!======================================================================
!  Intersection with Edge 7 (node 7-8, Face North-Top):
!======================================================================
                  xa = X_NODE(7)
                  ya = Y_NODE(7)
                  za = Z_NODE(7)

                  xb = X_NODE(8)
                  yb = Y_NODE(8)
                  zb = Z_NODE(8)

! Check if intersection occurs at corners

                  CALL IS_POINT_INSIDE_FACET(xa,ya,za,N,INSIDE_FACET)

                  IF(INSIDE_FACET) THEN   ! corner intersection at node 7

                     F_AT(IMJK) = ZERO

!                     BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR')  CALL ADD_FACET_AND_SET_BC_ID(IJK,N)

                     N8(1) = xb-xa
                     N8(2) = yb-ya
                     N8(3) = zb-za

                     dotproduct = DOT_PRODUCT(N8,NORM_FACE(N,:)) 

                     IF (DABS(dotproduct)>TOL_F) THEN
                        IF (F_AT(IJK)==UNDEFINED) F_AT(IJK) = -dotproduct
                     ENDIF

                  ENDIF


                  CALL IS_POINT_INSIDE_FACET(xb,yb,zb,N,INSIDE_FACET)

                  IF(INSIDE_FACET) THEN   ! corner intersection at node 8

                     F_AT(IJK) = ZERO

!                     BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,N)

                     N7(1) = xa-xb
                     N7(2) = ya-yb
                     N7(3) = za-zb

                     dotproduct = DOT_PRODUCT(N7,NORM_FACE(N,:)) 

                     IF (DABS(dotproduct)>TOL_F) THEN
                        IF (F_AT(IMJK)==UNDEFINED) F_AT(IMJK) = -dotproduct
                     ENDIF

                  ENDIF



! Check intersection within line 7-8, excluding corners


                  INTERSECT_FLAG = .FALSE.

                  IF(.NOT.INTERSECT_X(IJK)) CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,N,INTERSECT_FLAG,xc,yc,zc)


                  IF(INTERSECT_FLAG) THEN

                     IF(INTERSECT_X(IJK)) THEN

                        IF(DABS(Xint(IJK)-xc)>TOL_STL) THEN

                           INTERSECT_X(IJK) = .FALSE.        ! Ignore intersections when two intersections are detected on the same edge

                        ENDIF                  

                     ELSE

                        INTERSECT_X(IJK) = .TRUE.
                        Xint(IJK) = xc 

! Set values at corners if they are not zero

                        N7(1) = xa-xc
                        N7(2) = ya-yc
                        N7(3) = za-zc

                        IF(DABS(F_AT(IMJK))>TOL_F)   F_AT(IMJK) = -DOT_PRODUCT(N7,NORM_FACE(N,:))  

                        N8(1) = xb-xc
                        N8(2) = yb-yc
                        N8(3) = zb-zc

                        IF(DABS(F_AT(IJK))>TOL_F)   F_AT(IJK) = -DOT_PRODUCT(N8,NORM_FACE(N,:))  


!                        BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
!                        IF(JP<=J2) BC_ID(IJPK) = BC_ID_STL_FACE(N)
!                        IF(KP<=K2) BC_ID(IJKP) = BC_ID_STL_FACE(N) 
!                        IF(JP<=J2.AND.KP<=K2)BC_ID(IJPKP) = BC_ID_STL_FACE(N)  

                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,N)
                        IF(JP<=J2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJPK,N) 
                        IF(KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJKP,N) 
                        IF(JP<=J2.AND.KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJPKP,N)


                     ENDIF

                  ENDIF


                  IF(TYPE_OF_CELL=='U_MOMENTUM') THEN
                     IF(SNAP(IJK)) THEN
                        INTERSECT_X(IJK) = .TRUE.
                        Xn_int(IJK) = XG_E(I)
                     ENDIF
                  ENDIF




!======================================================================
!  Intersection with Edge 6 (node 6-8, Face East-Top):
!======================================================================
                  xa = X_NODE(6)
                  ya = Y_NODE(6)
                  za = Z_NODE(6)


! Check if intersection occurs at corners

                  CALL IS_POINT_INSIDE_FACET(xa,ya,za,N,INSIDE_FACET)

                  IF(INSIDE_FACET) THEN   ! corner intersection at node 6

                     F_AT(IJMK) = ZERO

!                     BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,N)

                     N8(1) = xb-xa
                     N8(2) = yb-ya
                     N8(3) = zb-za

                     dotproduct = DOT_PRODUCT(N8,NORM_FACE(N,:)) 

                     IF (DABS(dotproduct)>TOL_F) THEN
                        IF (F_AT(IJK)==UNDEFINED) F_AT(IJK) = -dotproduct
                     ENDIF

                  ENDIF


                  CALL IS_POINT_INSIDE_FACET(xb,yb,zb,N,INSIDE_FACET)

                  IF(INSIDE_FACET) THEN   ! corner intersection at node 8

                     F_AT(IJK) = ZERO

!                     BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
                     IF(TRIM(TYPE_OF_CELL).eq.'SCALAR')  CALL ADD_FACET_AND_SET_BC_ID(IJK,N)

                     N6(1) = xa-xb
                     N6(2) = ya-yb
                     N6(3) = za-zb

                     dotproduct = DOT_PRODUCT(N6,NORM_FACE(N,:)) 

                     IF (DABS(dotproduct)>TOL_F) THEN
                        IF (F_AT(IJMK)==UNDEFINED) F_AT(IJMK) = -dotproduct
                     ENDIF

                  ENDIF



! Check intersection within line 6-8, excluding corners



                  INTERSECT_FLAG = .FALSE.

                  IF(.NOT.INTERSECT_Y(IJK)) CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,N,INTERSECT_FLAG,xc,yc,zc)


                  IF(INTERSECT_FLAG) THEN

                     IF(INTERSECT_Y(IJK)) THEN

                        IF(DABS(Yint(IJK)-yc)>TOL_STL) THEN

                           INTERSECT_Y(IJK) = .FALSE. ! Ignore intersections when two intersections are detected on the same edge

                        ENDIF

                     ELSE


                        INTERSECT_Y(IJK) = .TRUE.
                        Yint(IJK) = yc 

! Set values at corners if they are not zero

                        N6(1) = xa-xc
                        N6(2) = ya-yc
                        N6(3) = za-zc

                        IF(DABS(F_AT(IJMK))>TOL_F)   F_AT(IJMK) = -DOT_PRODUCT(N6,NORM_FACE(N,:))  

                        N8(1) = xb-xc
                        N8(2) = yb-yc
                        N8(3) = zb-zc

                        IF(DABS(F_AT(IJK))>TOL_F)   F_AT(IJK) = -DOT_PRODUCT(N8,NORM_FACE(N,:))  


!                        BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
!                        IF(IP<=I2) BC_ID(IPJK) = BC_ID_STL_FACE(N) 
!                        IF(KP<=K2) BC_ID(IJKP) = BC_ID_STL_FACE(N)
!                        IF(IP<=I2.AND.KP<=K2)BC_ID(IPJKP) = BC_ID_STL_FACE(N) 

                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,N)
                        IF(IP<=I2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJK,N) 
                        IF(KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJKP,N) 
                        IF(IP<=I2.AND.KP<=K2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJKP,N)

                     ENDIF

                  ENDIF

                  IF(TYPE_OF_CELL=='V_MOMENTUM') THEN
                     IF(SNAP(IJK)) THEN
                        INTERSECT_Y(IJK) = .TRUE.
                        Ye_int(IJK) = YG_N(J)
                     ENDIF
                  ENDIF

                  IF(DO_K) THEN
!======================================================================
!  Intersection with Edge 11 (node 4-8, Face East-North):
!======================================================================
                     xa = X_NODE(4)
                     ya = Y_NODE(4)
                     za = Z_NODE(4)

! Check if intersection occurs at corners

                     CALL IS_POINT_INSIDE_FACET(xa,ya,za,N,INSIDE_FACET)

                     IF(INSIDE_FACET) THEN   ! corner intersection at node 4

                        F_AT(IJKM) = ZERO

!                        BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,N)

                        N8(1) = xb-xa
                        N8(2) = yb-ya
                        N8(3) = zb-za

                        dotproduct = DOT_PRODUCT(N8,NORM_FACE(N,:)) 

                        IF (DABS(dotproduct)>TOL_F) THEN
                           IF (F_AT(IJK)==UNDEFINED) F_AT(IJK) = -dotproduct
                        ENDIF

                     ENDIF


                     CALL IS_POINT_INSIDE_FACET(xb,yb,zb,N,INSIDE_FACET)

                     IF(INSIDE_FACET) THEN   ! corner intersection at node 8

                        F_AT(IJK) = ZERO

!                        BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
                        IF(TRIM(TYPE_OF_CELL).eq.'SCALAR')  CALL ADD_FACET_AND_SET_BC_ID(IJK,N)

                        N4(1) = xa-xb
                        N4(2) = ya-yb
                        N4(3) = za-zb

                        dotproduct = DOT_PRODUCT(N4,NORM_FACE(N,:)) 

                        IF (DABS(dotproduct)>TOL_F) THEN
                           IF (F_AT(IJKM)==UNDEFINED) F_AT(IJKM) = -dotproduct
                        ENDIF

                     ENDIF



! Check intersection within line 4-8, excluding corners


                     INTERSECT_FLAG = .FALSE. 

                     IF(.NOT.INTERSECT_Z(IJK)) CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,N,INTERSECT_FLAG,xc,yc,zc)

                     IF(INTERSECT_FLAG) THEN

                        IF(INTERSECT_Z(IJK)) THEN

                           IF(DABS(Zint(IJK)-zc)>TOL_STL) THEN

                              INTERSECT_Z(IJK) = .FALSE. ! Ignore intersections when two intersections are detected on the same edge

                           ENDIF

                        ELSE


                           INTERSECT_Z(IJK) = .TRUE.
                           Zint(IJK) = zc 


! Set values at corners if they are not zero

                           N4(1) = xa-xc
                           N4(2) = ya-yc
                           N4(3) = za-zc

                           IF(DABS(F_AT(IJKM))>TOL_F)   F_AT(IJKM) = -DOT_PRODUCT(N4,NORM_FACE(N,:))  

                           N8(1) = xb-xc
                           N8(2) = yb-yc
                           N8(3) = zb-zc

                           IF(DABS(F_AT(IJK))>TOL_F)   F_AT(IJK) = -DOT_PRODUCT(N8,NORM_FACE(N,:))  


!                           BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
!                           IF(IP<=I2)  BC_ID(IPJK) = BC_ID_STL_FACE(N)    
!                           IF(JP<=J2) BC_ID(IJPK) = BC_ID_STL_FACE(N)  
!                           IF(IP<=I2.AND.JP<=J2) BC_ID(IPJPK) = BC_ID_STL_FACE(N) 
                           
                           IF(TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJK,N)
                           IF(IP<=I2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJK,N) 
                           IF(JP<=J2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IJPK,N) 
                           IF(IP<=I2.AND.JP<=J2.AND.TRIM(TYPE_OF_CELL).eq.'SCALAR') CALL ADD_FACET_AND_SET_BC_ID(IPJPK,N)

                        ENDIF

                     ENDIF


                     IF(TYPE_OF_CELL=='W_MOMENTUM') THEN
                        IF(SNAP(IJK)) THEN
                           INTERSECT_Z(IJK) = .TRUE.
                           Zt_int(IJK) = ZG_T(K)
                        ENDIF
                     ENDIF

                  ENDIF


               ENDDO  ! I loop
            ENDDO  ! J loop
         ENDDO  ! K loop


      ENDDO  ! Loop over facets

      CURRENT_F = UNDEFINED



! Overwrite small values to set them to zero

      DO IJK = IJKSTART3, IJKEND3
         IF(DABS(F_AT(IJK))<TOL_STL) THEN
            F_AT(IJK)=ZERO
         ENDIF
      END DO

! Propagates node values to all interior cells
! in order defined by CAD_PROPAGATE_ORDER (outer loop)


      call send_recv(F_AT,2)


      SELECT CASE (CAD_PROPAGATE_ORDER)  

      CASE ('   ')


         N_UNDEFINED = UNDEFINED_I

         N_PROP = 0

         DO WHILE(N_UNDEFINED>0) 

            N_UNDEFINED = 0

            N_PROP = N_PROP + 1

            DO IJK = IJKSTART3, IJKEND3

               IF(INTERIOR_CELL_AT(IJK).AND.F_AT(IJK)==UNDEFINED) THEN

                  N_UNDEFINED = N_UNDEFINED + 1 

               ENDIF

            ENDDO


            DO IJK = IJKSTART3, IJKEND3

               IF(INTERIOR_CELL_AT(IJK).AND.F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN


                     IMJK = IM_OF(IJK)
                     IF(F_AT(IMJK)==UNDEFINED.AND.F_AT(IMJK)/=ZERO)  F_AT(IMJK)=F_AT(IJK)

                     IPJK = IP_OF(IJK)
                     IF(F_AT(IPJK)==UNDEFINED.AND.F_AT(IPJK)/=ZERO)  F_AT(IPJK)=F_AT(IJK)    

                     IJMK = JM_OF(IJK)
                     IF(F_AT(IJMK)==UNDEFINED.AND.F_AT(IJMK)/=ZERO)  F_AT(IJMK)=F_AT(IJK)

                     IJPK = JP_OF(IJK)
                     IF(F_AT(IJPK)==UNDEFINED.AND.F_AT(IJPK)/=ZERO)  F_AT(IJPK)=F_AT(IJK)

                     IJKM = KM_OF(IJK)
                     IF(F_AT(IJKM)==UNDEFINED.AND.F_AT(IJKM)/=ZERO)  F_AT(IJKM)=F_AT(IJK)

                     IJKP = KP_OF(IJK)
                     IF(F_AT(IJKP)==UNDEFINED.AND.F_AT(IJKP)/=ZERO)  F_AT(IJKP)=F_AT(IJK)


               ENDIF

            ENDDO


            call send_recv(F_AT,2)

            IF(N_PROP>N_PROPMAX.AND.N_UNDEFINED>0) THEN

               WRITE(*,*)'UNABLE TO PROPAGATE F_AT ARRAY FROM myPE=.', MyPE
               CALL MFIX_EXIT(myPE)

            ENDIF


         ENDDO ! WHILE Loop


         call send_recv(F_AT,2)





      CASE ('IJK')  

         DO I=ISTART3,IEND3
            DO K=KSTART3,KEND3
               F_FOUND = .FALSE.
               DO J=JSTART3,JEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO J=JSTART3,JEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO J=JSTART3,JEND3
            DO I=ISTART3,IEND3
               F_FOUND = .FALSE.
               DO K=KSTART3,KEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO K=KSTART3,KEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO K=KSTART3,KEND3
            DO J=JSTART3,JEND3
               F_FOUND = .FALSE.
               DO I=ISTART3,IEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO I=ISTART3,IEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

      CASE ('JKI')  

         DO J=JSTART3,JEND3
            DO I=ISTART3,IEND3
               F_FOUND = .FALSE.
               DO K=KSTART3,KEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO K=KSTART3,KEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO K=KSTART3,KEND3
            DO J=JSTART3,JEND3
               F_FOUND = .FALSE.
               DO I=ISTART3,IEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO I=ISTART3,IEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO I=ISTART3,IEND3
            DO K=KSTART3,KEND3
               F_FOUND = .FALSE.
               DO J=JSTART3,JEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO J=JSTART3,JEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)


      CASE ('KIJ')  



         DO K=KSTART3,KEND3
            DO J=JSTART3,JEND3
               F_FOUND = .FALSE.
               DO I=ISTART3,IEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO I=ISTART3,IEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO

         call send_recv(F_AT,2)

         DO I=ISTART3,IEND3
            DO K=KSTART3,KEND3
               F_FOUND = .FALSE.
               DO J=JSTART3,JEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO J=JSTART3,JEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO


         call send_recv(F_AT,2)

         DO J=JSTART3,JEND3
            DO I=ISTART3,IEND3
               F_FOUND = .FALSE.
               DO K=KSTART3,KEND3
                  IJK=FUNIJK(I,J,K)
                  IF(F_AT(IJK)/=UNDEFINED.AND.F_AT(IJK)/=ZERO) THEN
                     F_FOUND = .TRUE.
                     CURRENT_F = F_AT(IJK)
                     EXIT                
                  ENDIF         
               ENDDO
               IF(F_FOUND) THEN
                  DO K=KSTART3,KEND3
                     IJK=FUNIJK(I,J,K)
                     IF(F_AT(IJK)==UNDEFINED) THEN
                        F_AT(IJK)=CURRENT_F
                     ELSEIF(F_AT(IJK)/=ZERO) THEN
                        CURRENT_F = F_AT(IJK)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDDO


         call send_recv(F_AT,2)

         CASE DEFAULT
            IF(myPE == PE_IO) THEN
               WRITE(*,*)'CAD_INTERSECT.'
               WRITE(*,*)'UNKNOWN CAD_PROPAGATE_ORDER:',CAD_PROPAGATE_ORDER
               WRITE(*,*)'ACCEPTABLE VALUES:'
               WRITE(*,*)'IJK'
               WRITE(*,*)'JKI'
               WRITE(*,*)'KIJ'
            ENDIF
!            CALL MFIX_EXIT(myPE)

      END SELECT 



!      call send_recv(F_AT,2)

      RETURN
      
      END SUBROUTINE CAD_INTERSECT


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: ADD_FACET                                              C
!  Purpose: Add facet to list in IJK scalar cell                       C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 15-Oct-13  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE ADD_FACET_AND_SET_BC_ID(IJK,N)
    
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
      USE mpi_utility
      
      IMPLICIT NONE
      INTEGER :: IJK,N


      BC_ID(IJK) = BC_ID_STL_FACE(N)             ! Set tentative BC_ID
      
      IF(N_FACET_AT(IJK)<DIM_FACETS_PER_CELL) THEN

         N_FACET_AT(IJK) = N_FACET_AT(IJK) + 1
         LIST_FACET_AT(IJK,N_FACET_AT(IJK)) = N

      ELSE

         WRITE(*,*) ' FATAL ERROR: TOO MANY FACETS IN CELL: ', IJK
         CALL MFIX_EXIT(myPE)

      ENDIF

    
      
      END SUBROUTINE ADD_FACET_AND_SET_BC_ID


