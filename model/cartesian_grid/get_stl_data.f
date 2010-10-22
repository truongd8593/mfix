!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_STL_DATA                                           C
!  Purpose: reads face verticesa and normal vectors from an STL file   C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!  Reviewer:                                          Date: **-***-**  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_STL_DATA

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------

      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE run
      USE scalars
      USE funits 
      USE rxns
      USE compar             
      USE mpi_utility        
      USE progress_bar
      USE stl
      USE vtk
      USE quadric
      USE constant
      IMPLICIT NONE

      INTEGER :: POLY,V,N,NSKIP,IGNORED_FACETS
      LOGICAL :: PRESENT,KEEP_READING,IGNORE_CURRENT_FACET
      DOUBLE PRECISION ::v1x,v1y,v1z
      DOUBLE PRECISION ::v2x,v2y,v2z
      DOUBLE PRECISION ::v3x,v3y,v3z
      DOUBLE PRECISION ::x12,y12,z12,x13,y13,z13,dp,d12,d13
      DOUBLE PRECISION ::cos_angle,cos_small_angle
      DOUBLE PRECISION ::n1,n2,n3
      DOUBLE PRECISION ::ABSTRANS
      CHARACTER(LEN=32) ::TEST_CHAR,BUFF_CHAR





      WRITE(*,2000) 'READING geometry from geometry.stl...'

      INQUIRE(FILE='geometry.stl',EXIST=PRESENT)
      IF(.NOT.PRESENT) THEN
         IF(MyPE == PE_IO) THEN
            WRITE(*,"('(PE ',I3,'): input data file, ',A11,' is missing: run aborted')") &
            myPE,'geometry.stl'
         ENDIF
         CALL MFIX_EXIT(MYPE) 
      ENDIF
!     
!     
!     OPEN geometry.stl ASCII FILE
!     
      OPEN(UNIT=333, FILE='geometry.stl', STATUS='OLD', ERR=910) 

      WRITE(*,2000)'STL file opened. Starting reading data...'

      KEEP_READING = .TRUE.
      
      N_FACETS = 0
      IGNORED_FACETS = 0
      
      DO WHILE(KEEP_READING)
      
         READ(333,*,ERR=920,END=930) TEST_CHAR

         IF(TRIM(TEST_CHAR) == 'facet') THEN
 
            BACKSPACE(333)
            IGNORE_CURRENT_FACET = .FALSE.

            READ(333,*,ERR=920,END=930) BUFF_CHAR,BUFF_CHAR,N1,N2,N3  ! Read unit normal vector
            READ(333,*,ERR=920,END=930) 
            READ(333,*,ERR=920,END=930) BUFF_CHAR, V1x,V1y,V1z
            READ(333,*,ERR=920,END=930) BUFF_CHAR, V2x,V2y,V2z
            READ(333,*,ERR=920,END=930) BUFF_CHAR, V3x,V3y,V3z

            x12 = V2x - V1x
            y12 = V2y - V1y
            z12 = V2z - V1z

            x13 = V3x - V1x
            y13 = V3y - V1y
            z13 = V3z - V1z


            dp  = x12*x13 + y12*y13 + z12*z13
            d12 = dsqrt(x12**2+y12**2+z12**2)
            d13 = dsqrt(x13**2+y13**2+z13**2)

            IF((d12*d13)>TOL_STL) THEN
               cos_angle = dp/(d12*d13)
            ELSE
               cos_angle = ONE
            ENDIF

            cos_small_angle = dcos(STL_SMALL_ANGLE / 180.0 * PI)

            IF(DABS(cos_angle)>cos_small_angle) THEN
               IGNORE_CURRENT_FACET = .TRUE.    ! Ignore small facets
            ENDIF
            IF(DABS((N1**2+N2**2+N3**2)-ONE)>0.1) IGNORE_CURRENT_FACET = .TRUE.  ! Ignore facets with normal vector that is not normalized



            IF(IGNORE_CURRENT_FACET) THEN    
               IGNORED_FACETS = IGNORED_FACETS + 1
            ELSE                                                      ! Save vertex coordinates for valid facets
                                                                      ! and performs translation
               N_FACETS = N_FACETS + 1
     
               NORM_FACE(N_FACETS,1) = N1
               NORM_FACE(N_FACETS,2) = N2
               NORM_FACE(N_FACETS,3) = N3

               VERTEX(N_FACETS,1,1) = SCALE_STL*V1x + TX_STL  
               VERTEX(N_FACETS,1,2) = SCALE_STL*V1y + TY_STL  
               VERTEX(N_FACETS,1,3) = SCALE_STL*V1z + TZ_STL  

               VERTEX(N_FACETS,2,1) = SCALE_STL*V2x + TX_STL  
               VERTEX(N_FACETS,2,2) = SCALE_STL*V2y + TY_STL  
               VERTEX(N_FACETS,2,3) = SCALE_STL*V2z + TZ_STL  

               VERTEX(N_FACETS,3,1) = SCALE_STL*V3x + TX_STL  
               VERTEX(N_FACETS,3,2) = SCALE_STL*V3y + TY_STL  
               VERTEX(N_FACETS,3,3) = SCALE_STL*V3z + TZ_STL  

            ENDIF
 
         ELSEIF(TRIM(TEST_CHAR) == 'endsolid') THEN
 
            KEEP_READING = .FALSE.
 
         ENDIF
      
      ENDDO
      

      XMIN_STL = MINVAL(VERTEX(1:N_FACETS,:,1))
      XMAX_STL = MAXVAL(VERTEX(1:N_FACETS,:,1))
      YMIN_STL = MINVAL(VERTEX(1:N_FACETS,:,2))
      YMAX_STL = MAXVAL(VERTEX(1:N_FACETS,:,2))
      ZMIN_STL = MINVAL(VERTEX(1:N_FACETS,:,3))
      ZMAX_STL = MAXVAL(VERTEX(1:N_FACETS,:,3))


      WRITE(*,2000)'STL file successfully read.'
      WRITE(*,*)' Total number of facets read =',N_FACETS + IGNORED_FACETS
      WRITE(*,*)' Number of valid facets      =',N_FACETS
      WRITE(*,*)' Number of ignored facets    =',IGNORED_FACETS
      WRITE(*,*)' RANGE OF STL FILE:'
      IF(SCALE_STL/=ONE) THEN
         WRITE(*,5000)' AFTER SCALING BY A FACTOR OF ',SCALE_STL
      ENDIF
      ABSTRANS = dabs(TX_STL)+dabs(TY_STL)+dabs(TZ_STL)
      IF(ABSTRANS>TOL_STL) THEN
         WRITE(*,3000)' AFTER TRANSLATION OF (X,Y,Z)=',TX_STL,TY_STL,TZ_STL
      ENDIF
      WRITE(*,4000)'X-RANGE = ', XMIN_STL,XMAX_STL
      WRITE(*,4000)'Y-RANGE = ', YMIN_STL,YMAX_STL
      WRITE(*,4000)'Z-RANGE = ', ZMIN_STL,ZMAX_STL
      WRITE(*,4000)''


      XMIN_STL = XMIN_STL - 10.0*TOL_STL
      XMAX_STL = XMAX_STL + 10.0*TOL_STL
      YMIN_STL = YMIN_STL - 10.0*TOL_STL
      YMAX_STL = YMAX_STL + 10.0*TOL_STL
      ZMIN_STL = ZMIN_STL - 10.0*TOL_STL
      ZMAX_STL = ZMAX_STL + 10.0*TOL_STL


!      IF(XMIN_STL<ZERO) XMIN_STL=ZERO 
!      IF(XMAX_STL>XLENGTH) XMAX_STL=XLENGTH
!      IF(YMIN_STL<ZERO) YMIN_STL=ZERO 
!      IF(YMAX_STL>YLENGTH) YMAX_STL=YLENGTH
!      IF(ZMIN_STL<ZERO) ZMIN_STL=ZERO 
!      IF(ZMAX_STL>ZLENGTH) ZMAX_STL=ZLENGTH


      CLOSE(333)

      RETURN  

!======================================================================
!     HERE IF AN ERROR OCCURED OPENNING/READING THE FILE
!======================================================================
!     
 910  CONTINUE 
      WRITE (*, 1500) 
      CALL MFIX_EXIT(myPE) 
 920  CONTINUE 
      WRITE (*, 1600) 
      CALL MFIX_EXIT(myPE) 
 930  CONTINUE 
      WRITE (*, 1700) 
      CALL MFIX_EXIT(myPE) 
!     
 1500 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'Unable to open stl file',/1X,70('*')/) 
 1600 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'Error while reading stl file',/1X,70('*')/) 
 1700 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'End of file reached while reading stl file',/1X,70('*')/) 
 2000 FORMAT(1X,A) 
 2010 FORMAT(1X,A,I4,A)
 3000 FORMAT(1X,A,'(',F10.4,';',F10.4,';',F10.4,')')
 4000 FORMAT(1X,A,F10.4,' to ',F10.4)
 5000 FORMAT(1X,A,F10.4)
      END SUBROUTINE GET_STL_DATA

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: EVAL_STL_FCT                                           C
!  Purpose: Assigns a value to f_stl at any given point (X1,X2,X3)     C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!  Reviewer:                                          Date: **-***-**  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE EVAL_STL_FCT(X1,X2,X3,Q,f_stl,CLIP_FLAG,BCID)

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------

      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE run
      USE scalars
      USE funits 
      USE rxns
      USE compar             
      USE mpi_utility        
      USE progress_bar
      USE stl
      IMPLICIT NONE

      INTEGER :: POLY,V,N,NSKIP,N_FACES,CLOSEST_FACET,Q,BCID,C,COUNTER
      LOGICAL :: PRESENT,KEEP_READING,CLIP_FLAG,INSIDE_FACET
      LOGICAL :: X_OUT,Y_OUT,Z_OUT,POINT_IS_OUT,INTERSECT_FLAG
      DOUBLE PRECISION :: X1,X2,X3,XV1,XV2,XV3,XC1,XC2,XC3,DC,xc,yc,zc
      DOUBLE PRECISION :: RAY1,RAY2,RAY3
      DOUBLE PRECISION :: XC1_COPY,XC2_COPY,XC3_COPY
      DOUBLE PRECISION :: D(3),D_MIN,MINVAL_D,f_stl
      DOUBLE PRECISION :: DIST_INT(0:1000),D_I
      LOGICAL :: DUPLICATED_INTERSECTION
      DOUBLE PRECISION :: VFP1,VFP2,VFP3,NF1,NF2,NF3,DP,two
      INTEGER :: MINLOC_D(1),CLOSEST_VERTEX


      IF(N_FACETS < 1) RETURN

      BCID = -1

      D_MIN = UNDEFINED
      CLOSEST_FACET = 1

      X_OUT = (X1<XMIN_STL).OR.(X1>XMAX_STL)
      Y_OUT = (X2<YMIN_STL).OR.(X2>YMAX_STL)
      Z_OUT = (X3<ZMIN_STL).OR.(X3>ZMAX_STL)

      POINT_IS_OUT = X_OUT.OR.Y_OUT.OR.Z_OUT

      IF(POINT_IS_OUT) THEN
         f_stl = OUT_STL_VALUE
         CLIP_FLAG = .TRUE.
         BCID = -99
      ENDIF



      COUNTER = 0

      DIST_INT(COUNTER) = -1.0D0

      DO N = 1,N_FACETS

         CALL IS_POINT_INSIDE_FACET(X1,X2,X3,N,INSIDE_FACET)

         IF(INSIDE_FACET) THEN
            f_stl = ZERO
            CLIP_FLAG = .TRUE.
            BCID = STL_BC_ID
            RETURN
         ENDIF


         CALL GET_RAY_ORIGIN(X1,X2,X3,RAY1,RAY2,RAY3)

         CALL INTERSECT_LINE_WITH_FACET(X1,X2,X3,RAY1,RAY2,RAY3,N,INTERSECT_FLAG,xc,yc,zc)

         IF(INTERSECT_FLAG) THEN

            D_I = DSQRT((X1-xc)**2+(X2-yc)**2+(X3-zc)**2)

            DUPLICATED_INTERSECTION = .FALSE.
            DO C = 0,COUNTER
               IF(DABS(DIST_INT(C)-D_I)<TOL_STL) DUPLICATED_INTERSECTION = .TRUE.
            ENDDO

            IF(.NOT.DUPLICATED_INTERSECTION) THEN
               COUNTER = COUNTER + 1
               DIST_INT(COUNTER) = D_I
            ENDIF


         ENDIF

      ENDDO

      IF((COUNTER/2) * 2 == COUNTER) THEN    
         f_stl = OUT_STL_VALUE    ! even ==> outside
      ELSE                                   
         f_stl = -OUT_STL_VALUE   ! odd  ==> inside
      ENDIF


      CLIP_FLAG = .TRUE.
      BCID = STL_BC_ID

      RETURN  

!======================================================================
!     HERE IF AN ERROR OCCURED OPENNING/READING THE FILE
!======================================================================
!     
 910  CONTINUE 
      WRITE (*, 1500) 
      CALL MFIX_EXIT(myPE) 
 920  CONTINUE 
      WRITE (*, 1600) 
      CALL MFIX_EXIT(myPE) 
 930  CONTINUE 
      WRITE (*, 1700) 
      CALL MFIX_EXIT(myPE) 
!     
 1500 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'Unable to open stl file',/1X,70('*')/) 
 1600 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'Error while reading stl file',/1X,70('*')/) 
 1700 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'End of file reached while reading stl file',/1X,70('*')/) 
 2000 FORMAT(1X,A) 
 2010 FORMAT(1X,A,I4,A)

      END SUBROUTINE EVAL_STL_FCT




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: intersect_line_with_stl                                C
!  Purpose: Finds the intersection between a facet                     C
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
  SUBROUTINE INTERSECT_LINE_WITH_STL(xa,ya,za,xb,yb,zb,INTERSECT_FLAG,xc,yc,zc)
    
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
      USE STL
      
      IMPLICIT NONE
      DOUBLE PRECISION:: xa,ya,za,xb,yb,zb,xc,yc,zc
      DOUBLE PRECISION:: f1,f2
      INTEGER :: N,BCID
      LOGICAL :: CLIP_FLAG,CLIP_FLAG1,CLIP_FLAG2,INTERSECT_FLAG


!======================================================================
! This subroutine tries to find the intersection of a line AB with all  
! facets defined in the stl file
!
!
! 1) Verify that intersection is possible. Only facets and line AB that are
!
!======================================================================

      CALL EVAL_STL_FCT(xa,ya,za,N_FACETS,f1,CLIP_FLAG1,BCID)
      CALL EVAL_STL_FCT(xb,yb,zb,N_FACETS,f2,CLIP_FLAG2,BCID)

      CLIP_FLAG = (CLIP_FLAG1).AND.(CLIP_FLAG2)

      IF (CLIP_FLAG) THEN
         IF(DABS(f1)<TOL_STL) THEN  ! ignore intersection at corner
            xc = UNDEFINED
            yc = UNDEFINED
            zc = UNDEFINED
            INTERSECT_FLAG = .FALSE.
         ELSEIF(DABS(f2)<TOL_STL) THEN ! ignore intersection at corner
            xc = UNDEFINED
            yc = UNDEFINED
            zc = UNDEFINED
            INTERSECT_FLAG = .FALSE.
         ELSEIF(f1*f2 < ZERO) THEN      ! There is at least one intersection point
                                        ! Loop through all facets, find the first one that intersects the line and exit

            DO N = 1, N_FACETS
               CALL INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,N,INTERSECT_FLAG,xc,yc,zc)
               IF(INTERSECT_FLAG) EXIT
            ENDDO

            IF(.NOT.INTERSECT_FLAG) THEN
!               WRITE(*,*)'   Subroutine intersect_line_with_stl:'
!               WRITE(*,*)   'Unable to find the intersection of stl geometry'
!               WRITE(*,1000)'between (x1,y1,z1)= ', xa,ya,za
!               WRITE(*,1000)'   and  (x2,y2,z2)= ', xb,yb,zb
!               WRITE(*,1000)'f(x1,y1,z1) = ', f1
!               WRITE(*,1000)'f(x2,y2,z2) = ', f2
!               WRITE(*,1000)'Tolerance = ', TOL_STL
!               WRITE(*,*)   'MFiX will exit now.'             
!               CALL MFIX_EXIT(myPE) 
            xc = UNDEFINED
            yc = UNDEFINED
            zc = UNDEFINED
            INTERSECT_FLAG = .FALSE.



            ENDIF
         ELSE                          ! f1 and f2 have same sign => no intersection
            xc = UNDEFINED
            yc = UNDEFINED
            zc = UNDEFINED
            INTERSECT_FLAG = .FALSE.
         ENDIF
      ELSE                            ! ignore clipped region (not active)
          xc = UNDEFINED
          yc = UNDEFINED
          zc = UNDEFINED
          INTERSECT_FLAG = .FALSE.
      ENDIF

 1000 FORMAT(A,3(2X,G12.5)) 


      RETURN

      END SUBROUTINE INTERSECT_LINE_WITH_STL


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: intersect_line_with_facet                              C
!  Purpose: Finds the intersection between a facet                     C
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
  SUBROUTINE INTERSECT_LINE_WITH_FACET(xa,ya,za,xb,yb,zb,FACET,INTERSECT_FLAG,xc,yc,zc)
    
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
      USE STL
      
      IMPLICIT NONE
      DOUBLE PRECISION:: xa,ya,za,xb,yb,zb,xc,yc,zc
      INTEGER :: FACET
      DOUBLE PRECISION:: NFx,NFy,NFz,nabx,naby,nabz
      DOUBLE PRECISION :: dot_denom,dot_num,Inv_denom
      DOUBLE PRECISION :: VP1Ax,VP1Ay,VP1Az
      DOUBLE PRECISION :: Px,Py,Pz,V0x,V0y,V0z,V1x,V1y,V1z,V2x,V2y,V2z
      DOUBLE PRECISION :: dot00,dot01,dot02,dot11,dot12
      DOUBLE PRECISION :: t,u,v
      LOGICAL :: INSIDE_FACET,INTERSECT_FLAG


!======================================================================
! This subroutine tries to find the intersection of a line AB with one  
! of the facets defined in the stl file
!
!
! 1) Verify that intersection is possible. Only facets and line AB that are
!    not parallel are acceptable candidates
! 2) Find intersection point P of line AB with plane containing the facet. 
!    Valid intersection point must be between A and B to continue.
! 3) Verify that point P is inside triangular facet
!======================================================================
      
      INTERSECT_FLAG = .FALSE.

!======================================================================
!  Facet normal vector (normalized)
!======================================================================
      NFx = NORM_FACE(FACET,1)
      NFy = NORM_FACE(FACET,2)
      NFz = NORM_FACE(FACET,3)

!======================================================================
!  AB vector (NOT normalized)
!======================================================================
      nabx = xb - xa
      naby = yb - ya
      nabz = zb - za

      dot_denom = NFx*nabx + NFy*naby + NFz*nabz
!======================================================================
! 1) Verify that intersection is possible. Only facets and line AB that are
!    not parallel are acceptable candidates
!======================================================================
      IF(dabs(dot_denom)<TOL_STL) THEN

         INTERSECT_FLAG = .FALSE.              ! No intersection (facet nornal 
                                               ! and AB are perpendicular
         RETURN

      ELSE

!======================================================================
! 2) Find intersection point P of line AB with plane containing the facet
!    Line AB is parametrized with parameter t (from t=0 at A to t=1 at B)
!======================================================================
         VP1Ax = VERTEX(FACET,1,1) - xa
         VP1Ay = VERTEX(FACET,1,2) - ya
         VP1Az = VERTEX(FACET,1,3) - za

         dot_num = NFx*VP1Ax + NFy*VP1Ay + NFz*VP1Az

         t = dot_num / dot_denom

!======================================================================
! 3) Verify that point P is inside triangular facet
!    Facet is parametrized with (u,v), a point is inside triange if:
!    u   >= 0  , and
!    v   >= 0  , and
!    u+v <= 1
!======================================================================
         IF((t>=ZERO).AND.(t<=ONE)) THEN       ! Intersection between A and B
                                               ! Now test if intersection point is inside triangle

            Px = xa + t*nabx
            Py = ya + t*naby
            Pz = za + t*nabz

            CALL IS_POINT_INSIDE_FACET(Px,Py,Pz,FACET,INSIDE_FACET)

            IF(INSIDE_FACET) THEN
               INTERSECT_FLAG = .TRUE.              ! Valid intersection
               xc = Px                              ! point P is inside triangle)
               yc = Py
               zc = Pz
            ELSE
               INTERSECT_FLAG = .FALSE.             ! Invalid intersection
               RETURN                               ! point P is outside triangle)
            ENDIF

         ELSE

            INTERSECT_FLAG = .FALSE.           ! No intersection between A and B
            RETURN

         ENDIF

      ENDIF


 1000 FORMAT(A,3(2X,G12.5)) 


      RETURN

      END SUBROUTINE INTERSECT_LINE_WITH_FACET

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: IS_POINT_INSIDE_FACET                                  C
!  Purpose: Verifies that point P is inside facet                      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 21-Feb-08  C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: #                                                           C
!  Purpose: #                                                          C
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
  SUBROUTINE IS_POINT_INSIDE_FACET(Px,Py,Pz,FACET,INSIDE_FACET)
    
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
      USE STL
      
      IMPLICIT NONE
      INTEGER :: FACET,VV
      DOUBLE PRECISION :: NFx,NFy,NFz
      DOUBLE PRECISION :: Px,Py,Pz,V0x,V0y,V0z,V1x,V1y,V1z,V2x,V2y,V2z
      DOUBLE PRECISION :: dot00,dot01,dot02,dot11,dot12,dot_check
      DOUBLE PRECISION :: Inv_denom
      DOUBLE PRECISION :: t,u,v
      LOGICAL :: U_POSITIVE,V_POSITIVE,UPVL1,INSIDE_FACET


      DOUBLE PRECISION :: Vx,Vy,Vz
      DOUBLE PRECISION :: D(3),MINVAL_D

!======================================================================
!  If point P is very close to one of the Facet vertices, 
!  consider that P belongs to facet, and return.
!======================================================================
         DO VV = 1,3
            Vx = VERTEX(FACET,VV,1)
            Vy = VERTEX(FACET,VV,2)
            Vz = VERTEX(FACET,VV,3)
            D(VV) = DSQRT((Px - Vx)**2 + (Py - Vy)**2 + (Pz - Vz)**2 )
         ENDDO

         MINVAL_D = MINVAL(D)

         IF(MINVAL_D < TOL_STL) THEN
            INSIDE_FACET = .TRUE.
            RETURN
         ENDIF

!======================================================================
!  Facet normal vector (normalized)
!======================================================================
      NFx = NORM_FACE(FACET,1)
      NFy = NORM_FACE(FACET,2)
      NFz = NORM_FACE(FACET,3)

!======================================================================
! This subroutine verifies that point P is inside triangular facet
! Facet is parametrized with (u,v), a point is inside triange if:
! u   >= 0  , and
! v   >= 0  , and
! u+v <= 1
!======================================================================

      V0x = VERTEX(FACET,2,1) - VERTEX(FACET,1,1) 
      V0y = VERTEX(FACET,2,2) - VERTEX(FACET,1,2) 
      V0z = VERTEX(FACET,2,3) - VERTEX(FACET,1,3) 

      V1x = VERTEX(FACET,3,1) - VERTEX(FACET,1,1) 
      V1y = VERTEX(FACET,3,2) - VERTEX(FACET,1,2) 
      V1z = VERTEX(FACET,3,3) - VERTEX(FACET,1,3) 

      V2x = Px - VERTEX(FACET,1,1) 
      V2y = Py - VERTEX(FACET,1,2) 
      V2z = Pz - VERTEX(FACET,1,3) 

      dot_check = NFx*V2x + NFy*V2y + NFz*V2z

      IF(dabs(dot_check)>TOL_STL) THEN          ! reject points that do not 
         INSIDE_FACET = .FALSE.                 ! belong to plane containing facet
         RETURN
      ENDIF 

      dot00 = V0x*V0x + V0y*V0y + V0z*V0z
      dot01 = V0x*V1x + V0y*V1y + V0z*V1z
      dot02 = V0x*V2x + V0y*V2y + V0z*V2z
      dot11 = V1x*V1x + V1y*V1y + V1z*V1z
      dot12 = V1x*V2x + V1y*V2y + V1z*V2z
      
      Inv_denom = ONE / (dot00*dot11 - dot01*dot01)   

      u = (dot11*dot02 - dot01*dot12) * Inv_denom
      v = (dot00*dot12 - dot01*dot02) * Inv_denom

      U_POSITIVE = (u>=-TOL_STL)
      V_POSITIVE = (v>=-TOL_STL)
      UPVL1 = ((u+v)<=ONE+TOL_STL)

      INSIDE_FACET = (U_POSITIVE.AND.V_POSITIVE.AND.UPVL1) 


      RETURN

      END SUBROUTINE IS_POINT_INSIDE_FACET





!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GET_RAY_ORIGIN                                         C
!  Purpose: Defines the origin of the ray that is traced to            C
!           determie if a point is inside or outside stl geometry      C
!                                                                      C
!  Author: Jeff Dietiker                              Date: 30-JAN-09  C
!  Reviewer:                                          Date: **-***-**  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!                                                                      C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE GET_RAY_ORIGIN(X1,X2,X3,RAY1,RAY2,RAY3)

!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------

      USE param 
      USE param1 
      USE physprop
      USE fldvar
      USE run
      USE scalars
      USE funits 
      USE rxns
      USE compar             
      USE mpi_utility        
      USE progress_bar
      USE stl
      IMPLICIT NONE

      DOUBLE PRECISION :: X1,X2,X3
      DOUBLE PRECISION :: RAY1,RAY2,RAY3


         SELECT CASE (TRIM(STL_RAY_DIR))  
            CASE ('X-')
               RAY1 = -XLENGTH
               RAY2 = X2
               RAY3 = X3
            CASE ('Y-')  
               RAY1 = X1
               RAY2 = -YLENGTH
               RAY3 = X3
            CASE ('Z-') 
               RAY1 = X1
               RAY2 = X2
               RAY3 = -ZLENGTH

            CASE ('X+')
               RAY1 = 2.0D0*XLENGTH
               RAY2 = X2
               RAY3 = X3
            CASE ('Y+')  
               RAY1 = X1
               RAY2 = 2.0D0*YLENGTH
               RAY3 = X3
            CASE ('Z+') 
               RAY1 = X1
               RAY2 = X2
               RAY3 = 2.0D0*ZLENGTH


            CASE ('MIN')
               RAY1 = -XLENGTH
               RAY2 = -YLENGTH
               RAY3 = -ZLENGTH

            CASE ('MAX')  
               RAY1 = 2.0D0*XLENGTH
               RAY2 = 2.0D0*YLENGTH
               RAY3 = 2.0D0*ZLENGTH

            CASE DEFAULT
               WRITE(*,*)'INPUT ERROR: STL_RAY_DIR HAS INCORRECT VALUE: ',STL_RAY_DIR
               WRITE(*,*)'PLEASE CORRECT MFIX.DAT AND TRY AGAIN.'
               CALL MFIX_EXIT(MYPE)
         END SELECT

      RETURN  

!======================================================================
!     HERE IF AN ERROR OCCURED OPENNING/READING THE FILE
!======================================================================
!     
 910  CONTINUE 
      WRITE (*, 1500) 
      CALL MFIX_EXIT(myPE) 
 920  CONTINUE 
      WRITE (*, 1600) 
      CALL MFIX_EXIT(myPE) 
 930  CONTINUE 
      WRITE (*, 1700) 
      CALL MFIX_EXIT(myPE) 
!     
 1500 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'Unable to open stl file',/1X,70('*')/) 
 1600 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'Error while reading stl file',/1X,70('*')/) 
 1700 FORMAT(/1X,70('*')//' From: GET_STL_DATA',/' Message: ',&
      'End of file reached while reading stl file',/1X,70('*')/) 
 2000 FORMAT(1X,A) 
 2010 FORMAT(1X,A,I4,A)

      END SUBROUTINE GET_RAY_ORIGIN


