!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: DEFINE_QUADRICS                                        C
!  Purpose: Defines all matrices used to evaluate quadrics             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE DEFINE_QUADRICS

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

      IMPLICIT NONE
      DOUBLE PRECISION:: x1,x2,x3
      DOUBLE PRECISION, DIMENSION(3,3) :: Rx,Ry,Rz,C_QUADRIC,R_QUADRIC
      DOUBLE PRECISION :: Two_H
      LOGICAL :: CLIP_FLAG,INTERSECT_FLAG

!======================================================================
! Quadric Normal form : lambda_x * x^2 + lambda_y * y^2 + lambda_z * z^2 + d = 0
!======================================================================

      DO QUADRIC_ID = 1 , N_QUADRIC

!       Build translation matrices
        CALL BUILD_1x3_MATRIX(t_x(QUADRIC_ID),t_y(QUADRIC_ID),t_z(QUADRIC_ID),T_QUADRIC(:,:,QUADRIC_ID))
!       Build characteristic matrices
        CALL BUILD_C_QUADRIC_MATRIX(lambda_x(QUADRIC_ID),lambda_y(QUADRIC_ID),lambda_z(QUADRIC_ID),C_QUADRIC)
!       Build Rotation matrices
        CALL BUILD_X_ROTATION_MATRIX(Theta_x(QUADRIC_ID), Rx)
        CALL BUILD_Y_ROTATION_MATRIX(Theta_y(QUADRIC_ID), Ry)
        CALL BUILD_Z_ROTATION_MATRIX(Theta_z(QUADRIC_ID), Rz)
        R_QUADRIC = MATMUL(Rz,MATMUL(Ry,Rx))
!       Build A-matrices
        A_QUADRIC(:,:,QUADRIC_ID) = MATMUL(TRANSPOSE(R_QUADRIC),MATMUL(C_QUADRIC,R_QUADRIC))

      END DO

      ! Activate Quadric
      QUADRIC_ID = N_QUADRIC


      RETURN


      END SUBROUTINE DEFINE_QUADRICS


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_F_QUADRIC                                          C
!  Purpose: Evaluate the function f(x,y,z) defining the quadric        C
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
        SUBROUTINE GET_F_QUADRIC(x1,x2,x3,Q_ID,f,CLIP_FLAG)

      USE parallel
      USE compar
      USE sendrecv
      USE quadric
      USE cutcell

      USE quadric


      IMPLICIT NONE

      DOUBLE PRECISION :: x1,x2,x3,xt,yt,zt,R1,R2
      DOUBLE PRECISION :: f,fq,fxe,fxw,fyn,fys,fzt,fzb,fclip
      DOUBLE PRECISION :: fxmin,fxmax,fymin,fymax,fzmin,fzmax
      DOUBLE PRECISION, DIMENSION(1,3) :: X_VECTOR,XMT
      DOUBLE PRECISION, DIMENSION(3,1) :: TXMT
      DOUBLE PRECISION, DIMENSION(1,1) :: TEMP_1x1
      INTEGER :: Q,Q_ID,BCID
      LOGICAL :: CLIP_X,CLIP_Y,CLIP_Z,CLIP_FLAG
      LOGICAL :: PIECE_X,PIECE_Y,PIECE_Z,PIECE_FLAG
      CHARACTER (LEN = 7) :: METHOD

!======================================================================

       PIECE_X = (piece_xmin(Q_ID) <= x1).AND.( x1 <= piece_xmax(Q_ID))
       PIECE_Y = (piece_ymin(Q_ID) <= x2).AND.( x2 <= piece_ymax(Q_ID))
       PIECE_Z = (piece_zmin(Q_ID) <= x3).AND.( x3 <= piece_zmax(Q_ID))

       PIECE_FLAG = (PIECE_X.AND.PIECE_Y.AND.PIECE_Z)

       IF(.NOT.PIECE_FLAG) THEN
          f = UNDEFINED
          RETURN
       ENDIF

       CLIP_X = (clip_xmin(Q_ID) <= x1).AND.( x1 <= clip_xmax(Q_ID))
       CLIP_Y = (clip_ymin(Q_ID) <= x2).AND.( x2 <= clip_ymax(Q_ID))
       CLIP_Z = (clip_zmin(Q_ID) <= x3).AND.( x3 <= clip_zmax(Q_ID))

       CLIP_FLAG = (CLIP_X.AND.CLIP_Y.AND.CLIP_Z)


         IF(TRIM(quadric_form(Q_ID))=='PLANE') THEN

            f = lambda_x(Q_ID)*x1 + lambda_y(Q_ID)*x2 +lambda_z(Q_ID)*x3 + dquadric(Q_ID)

         ELSEIF(TRIM(quadric_form(Q_ID))=='TORUS_INT') THEN

            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

            R1 = Torus_R1(Q_ID)
            R2 = Torus_R2(Q_ID)

            f = -(4*(xt**2+zt**2)*R1**2-(xt**2+yt**2+zt**2+R1**2-R2**2)**2)


         ELSEIF(TRIM(quadric_form(Q_ID))=='TORUS_EXT') THEN

            xt = x1-t_x(Q_ID)
            yt = x2-t_y(Q_ID)
            zt = x3-t_z(Q_ID)

            R1 = Torus_R1(Q_ID)
            R2 = Torus_R2(Q_ID)

            f = 4*(xt**2+zt**2)*R1**2-(xt**2+yt**2+zt**2+R1**2-R2**2)**2


         ELSE

            CALL BUILD_1x3_MATRIX(x1,x2,x3,X_VECTOR)

            XMT = X_VECTOR - T_QUADRIC(:,:,Q_ID)
            TXMT = TRANSPOSE(XMT)

            TEMP_1x1 = MATMUL(XMT,MATMUL(A_QUADRIC(:,:,Q_ID),TXMT))

            f = TEMP_1x1(1,1) + dquadric(Q_ID)

         ENDIF

! Each clipping limit is treated as a plane. For example, fxmin is
! the equation of the plane describing x=xmin, and a value of fxmin
! is compared with the current value of f to determine if the location
! is part of the computational domain. The comparison (min of max)
! follows the same logis as the 'AND' (max) , or 'OR' (min)
! logic when combining two quadrics.
! The clipping procedure is ignored when CLIP_FLAG is .FALSE.
! This will happen when we are in a 'PIECEWISE' group



         IF(FLUID_IN_CLIPPED_REGION(Q_ID)) THEN

            IF(clip_xmin(Q_ID)/=UNDEFINED) THEN
               fxmin = -(clip_xmin(Q_ID)-x1)
               f = dmin1(f,fxmin)
            ENDIF

            IF(clip_xmax(Q_ID)/=UNDEFINED) THEN
               fxmax = -(x1-clip_xmax(Q_ID))
               f = dmin1(f,fxmax)
            ENDIF

            IF(clip_ymin(Q_ID)/=UNDEFINED) THEN
               fymin = -(clip_ymin(Q_ID)-x2)
               f = dmin1(f,fymin)
            ENDIF

            IF(clip_ymax(Q_ID)/=UNDEFINED) THEN
               fymax = -(x2-clip_ymax(Q_ID))
               f = dmin1(f,fymax)
            ENDIF

            IF(clip_zmin(Q_ID)/=UNDEFINED) THEN
               fzmin = -(clip_zmin(Q_ID)-x3)
               f = dmin1(f,fzmin)
            ENDIF

            IF(clip_zmax(Q_ID)/=UNDEFINED) THEN
               fzmax = -(x3-clip_zmax(Q_ID))
               f = dmin1(f,fzmax)
            ENDIF

         ELSE

            IF(clip_xmin(Q_ID)/=UNDEFINED) THEN
               fxmin = clip_xmin(Q_ID)-x1
               f = dmax1(f,fxmin)
            ENDIF

            IF(clip_xmax(Q_ID)/=UNDEFINED) THEN
               fxmax = x1-clip_xmax(Q_ID)
               f = dmax1(f,fxmax)
            ENDIF

            IF(clip_ymin(Q_ID)/=UNDEFINED) THEN
               fymin = clip_ymin(Q_ID)-x2
               f = dmax1(f,fymin)
            ENDIF

            IF(clip_ymax(Q_ID)/=UNDEFINED) THEN
               fymax = x2-clip_ymax(Q_ID)
               f = dmax1(f,fymax)
            ENDIF

            IF(clip_zmin(Q_ID)/=UNDEFINED) THEN
               fzmin = clip_zmin(Q_ID)-x3
               f = dmax1(f,fzmin)
            ENDIF

            IF(clip_zmax(Q_ID)/=UNDEFINED) THEN
               fzmax = x3-clip_zmax(Q_ID)
               f = dmax1(f,fzmax)
            ENDIF

         ENDIF



      RETURN
      END SUBROUTINE GET_F_QUADRIC


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: REASSIGN_QUADRIC                                       C
!  Purpose: Reassign the quadric based on location                     C
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
      SUBROUTINE REASSSIGN_QUADRIC(x1,x2,x3,GROUP,Q_ID)

      USE compar
      USE quadric

      IMPLICIT NONE

      DOUBLE PRECISION x1,x2,x3
      INTEGER :: I,Q_ID,GROUP,GS,P
      LOGICAL :: PIECE_X,PIECE_Y,PIECE_Z,PIECE_FLAG
      CHARACTER(LEN=9) :: GR

      Q_ID = 0

      GS = GROUP_SIZE(GROUP)
      GR = TRIM(GROUP_RELATION(GROUP))

      IF( GR /= 'PIECEWISE') RETURN

      DO P = 1 , GS

         I = GROUP_Q(GROUP,P)

         PIECE_X = (piece_xmin(I) <= x1).AND.( x1 <= piece_xmax(I))
         PIECE_Y = (piece_ymin(I) <= x2).AND.( x2 <= piece_ymax(I))
         PIECE_Z = (piece_zmin(I) <= x3).AND.( x3 <= piece_zmax(I))

         PIECE_FLAG = (PIECE_X.AND.PIECE_Y.AND.PIECE_Z)

         IF (PIECE_FLAG) Q_ID = I

      ENDDO

      IF(Q_ID == 0 ) THEN
         WRITE(*,*)' No Quadric defined at current location x,y,z=', x1,x2,x3
         WRITE(*,*)' Please Check piecewise limits of quadric(s)'
         WRITE(*,*)' Mfix will exit now.'
         CALL MFIX_EXIT(myPE)
      ENDIF

      RETURN


      END SUBROUTINE REASSSIGN_QUADRIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_1x3_MATRIX                                       C
!  Purpose: Catesian Grid - Build a (1x3) matrix                       C
!           from 3 scalars                                             C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_1x3_MATRIX(scalar1,scalar2,scalar3,M1x3)

      IMPLICIT NONE

      DOUBLE PRECISION:: scalar1,scalar2,scalar3
      DOUBLE PRECISION, DIMENSION(1,3) :: M1x3
      M1x3(1,1) = scalar1
      M1x3(1,2) = scalar2
      M1x3(1,3) = scalar3

      RETURN

  END SUBROUTINE BUILD_1x3_MATRIX


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_C_QUADRIC_MATRIX                                 C
!  Purpose: Catesian Grid - Build a (3x3) diagonal matrix              C
!           whose diagonal elements are the characteristic values of   C
!           the quadric                                                C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_C_QUADRIC_MATRIX(lambda1,lambda2,lambda3,C_QUADRIC)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: lambda1,lambda2,lambda3
      DOUBLE PRECISION, DIMENSION(3,3) :: C_QUADRIC

!      Transpose is used because matrices are stored column-wise

      C_QUADRIC = TRANSPOSE(RESHAPE((/                                           &
                                          lambda1  ,   ZERO      ,  ZERO        ,&
                                          ZERO     ,   lambda2   ,  ZERO        ,&
                                          ZERO     ,   ZERO      ,  lambda3   /),&
                                                                                  (/3,3/)))
      RETURN

      END SUBROUTINE BUILD_C_QUADRIC_MATRIX



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_X_ROTATION_MATRIX                                C
!  Purpose: Catesian Grid - Build a (3x3) rotation matrix about x-axis,C
!           given the rotation angle in degrees                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_X_ROTATION_MATRIX(Theta, R)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: Theta
      DOUBLE PRECISION, DIMENSION(3,3) :: R

      theta = theta * (pi/180.0D0) ! Rotation angle about x-axis (radians)

!      Transpose is used because matrices are stored column-wise

      R = TRANSPOSE(RESHAPE((/                                                  &
                                  ONE   ,   ZERO         ,  ZERO               ,&
                                  ZERO  ,   dcos(theta)  ,  dsin(theta)        ,&
                                  ZERO  ,  -dsin(theta)  ,  dcos(theta)      /),&
                                                                                  (/3,3/)))

      RETURN

      END SUBROUTINE BUILD_X_ROTATION_MATRIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_Y_ROTATION_MATRIX                                C
!  Purpose: Catesian Grid - Build a (3x3) rotation matrix about y-axis,C
!           given the rotation angle in degrees                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_Y_ROTATION_MATRIX(Theta, R)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: Theta
      DOUBLE PRECISION, DIMENSION(3,3) :: R

      theta = theta * (pi/180.0D0) ! Rotation angle about x-axis (radians)

!      Transpose is used because matrices are stored column-wise

      R = TRANSPOSE(RESHAPE((/                                                &
                                  dcos(theta)  ,  ZERO  ,  dsin(theta)       ,&
                                  ZERO         ,  ONE   ,  ZERO              ,&
                                 -dsin(theta)  ,  ZERO  ,  dcos(theta)     /),&
                                                                                (/3,3/)))

      RETURN

      END SUBROUTINE BUILD_Y_ROTATION_MATRIX


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: BUILD_Z_ROTATION_MATRIX                                C
!  Purpose: Catesian Grid - Build a (3x3) rotation matrix about z-axis,C
!           given the rotation angle in degrees                        C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE BUILD_Z_ROTATION_MATRIX(Theta, R)

      USE param1, only: one, zero
      USE constant, only: pi

      IMPLICIT NONE

      DOUBLE PRECISION:: Theta
      DOUBLE PRECISION, DIMENSION(3,3) :: R

      theta = theta * (pi/180.0D0) ! Rotation angle about x-axis (radians)

!      Transpose is used because matrices are stored column-wise

      R = TRANSPOSE(RESHAPE((/                                                 &
                                   dcos(theta)  ,  dsin(theta)  ,  ZERO       ,&
                                  -dsin(theta)  ,  dcos(theta)  ,  ZERO       ,&
                                   ZERO         ,  ZERO         ,  ONE      /),&
                                                                                 (/3,3/)))

      RETURN

      END SUBROUTINE BUILD_Z_ROTATION_MATRIX

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CROSS_PRODUCT                                          C
!  Purpose: Performs the cross product between two vectors             C
!           C = A x B                                                  C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE CROSS_PRODUCT(A,B,C)

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(3) :: A,B,C

      C(1) = A(2) * B(3) - A(3) * B(2)
      C(2) = A(3) * B(1) - A(1) * B(3)
      C(3) = A(1) * B(2) - A(2) * B(1)

      RETURN

      END SUBROUTINE CROSS_PRODUCT


