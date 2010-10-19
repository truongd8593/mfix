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
      
      USE fldvar

      IMPLICIT NONE

      INTEGER :: IJK
      INTEGER :: I,J,K
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

!======================================================================
! Gather Data  and writes surface(s) defined by all cut cells 
!======================================================================

      CALL GATHER_DATA

      IF(WRITE_VTK_FILES) THEN
         CALL WRITE_CUT_SURFACE_VTK
         IF(NO_K) THEN
            WRITE_ANI_CUTCELL = .TRUE.
            CALL WRITE_VTK_FILE
         ELSE
            WRITE_ANI_CUTCELL = .FALSE.
         ENDIF
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
                  CALL GET_F_QUADRIC(x1,x2,x3,Q_ID,F_G(GROUP,0),CLIP_FLAG)
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

         CASE DEFAULT

            WRITE(*,*)'ERROR IN SUBROUTINE EVAL_F.'
            WRITE(*,*)'UNKNOWN METHOD:',METHOD
            WRITE(*,*)'ACCEPTABLE METHODS:'
            WRITE(*,*)'QUADRIC'
            WRITE(*,*)'POLYGON'
            WRITE(*,*)'USR_DEF'
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

      if (CLIP_FLAG) then
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
      INTEGER :: IJK,I,J,K,Q_ID,Q_ID2,N_int_x,N_int_y,N_int_z
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

      DO Q_ID = 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_X(IJK),xc,yc,zc)
         IF(INTERSECT_X(IJK)) Xi = xc
      ENDDO

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

      DO Q_ID = 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Y(IJK),xc,yc,zc)
         IF(INTERSECT_Y(IJK)) Yi = yc
      ENDDO

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

      DO Q_ID = 1, N_USR_DEF
         CALL INTERSECT_LINE('USR_DEF',xa,ya,za,xb,yb,zb,Q_ID,INTERSECT_Z(IJK),xc,yc,zc)
         IF(INTERSECT_Z(IJK)) Zi = zc
      ENDDO


         IF(TYPE_OF_CELL=='W_MOMENTUM') THEN
            IF(SNAP(IJK)) THEN
               INTERSECT_Z(IJK) = .TRUE.
               K = K_OF(IJK) 
               Zi = ZG_T(K)
            ENDIF
         ENDIF

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
      
      IMPLICIT NONE
      CHARACTER (LEN=*) :: TYPE_OF_CELL
      INTEGER :: IJK,I,J,K,IM,JM,KM,IP,JP,KP
      INTEGER :: IMJK,IPJK,IJMK,IJPK,IJKM,IJKP,IMJPK,IMJKP,IPJMK,IJMKP,IPJKM,IJPKM
      DOUBLE PRECISION :: xa,ya,za,xb,yb,zb,xc,yc,zc
      DOUBLE PRECISION :: Xi,Yi,Zi
      DOUBLE PRECISION :: DFC,DFC_MAX,Fa,Fb
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
            INTERSECT_X(IMJK) = .FALSE.
            INTERSECT_Y(IMJK)  = .FALSE.
            INTERSECT_Y(IMJPK) = .FALSE.
            INTERSECT_Z(IMJK)  = .FALSE.
            INTERSECT_Z(IMJKP) = .FALSE.

            SNAP(IMJK) = .TRUE.
                
         ENDIF


         DFC = DABS(Xi-xb) ! DISTANCE FROM CORNER (NODE 8)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING X-INTERSECTION ALONG EDGE 7 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF


            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_X(IPJK) = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            INTERSECT_Y(IJPK) = .FALSE.
            INTERSECT_Z(IJK)  = .FALSE.
            INTERSECT_Z(IJKP) = .FALSE.

            SNAP(IJK) = .TRUE.
                
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

            INTERSECT_X(IJMK)  = .FALSE.
            INTERSECT_X(IPJMK) = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            INTERSECT_Y(IJMK) = .FALSE.
            INTERSECT_Z(IJMK)  = .FALSE.
            INTERSECT_Z(IJMKP) = .FALSE.

            SNAP(IJMK) = .TRUE.
                
         ENDIF


         DFC = DABS(Yi-yb) ! DISTANCE FROM CORNER (NODE 8)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Y-INTERSECTION ALONG EDGE 6 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF 

            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_X(IPJK) = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            INTERSECT_Y(IJPK) = .FALSE.
            INTERSECT_Z(IJK)  = .FALSE.
            INTERSECT_Z(IJKP) = .FALSE.

            SNAP(IJK) = .TRUE.
                
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

      IF(INTERSECT_Y(IJK)) THEN

         DFC = DABS(Zi-Za) ! DISTANCE FROM CORNER (NODE 4)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 4'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF 

            INTERSECT_X(IJKM)  = .FALSE.
            INTERSECT_X(IPJKM) = .FALSE.
            INTERSECT_Y(IJKM)  = .FALSE.
            INTERSECT_Y(IJPKM) = .FALSE.
            INTERSECT_Z(IJK)  = .FALSE.
            INTERSECT_Z(IJKM) = .FALSE.
                
            SNAP(IJKM) = .TRUE.

         ENDIF


         DFC = DABS(Zi-Zb) ! DISTANCE FROM CORNER (NODE 8)

         IF(DFC < DFC_MAX) THEN
            IF(PRINT_WARNINGS) THEN
               WRITE(*,*)'MERGING Z-INTERSECTION ALONG EDGE 11 ONTO NODE 8'
               WRITE(*,*)'AT IJK,I,J,K=',IJK,I,J,K
            ENDIF

            INTERSECT_X(IJK)  = .FALSE.
            INTERSECT_X(IPJK) = .FALSE.
            INTERSECT_Y(IJK)  = .FALSE.
            INTERSECT_Y(IJPK) = .FALSE.
            INTERSECT_Z(IJK)  = .FALSE.
            INTERSECT_Z(IJKP) = .FALSE.

            SNAP(IJK) = .TRUE.
                
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



