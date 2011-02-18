!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GENERATE_PARTICLE_CONFIG                               C
!
!  Purpose: Generate particle configuration based on input             
!
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 23-Oct-08  C
!  Comments: Added a new routine for clarity of functions              C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GENERATE_PARTICLE_CONFIG
      
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------       
      INTEGER L, K, M
      INTEGER CHECK_MPI
      INTEGER PART_COUNT
      DOUBLE PRECISION DIST, R_LM, DOML(DIMN)
!-----------------------------------------------       

      WRITE(*,'(3X,A)') &
         '---------- START GENERATE_PARTICLE_CONFIG ---------->'

      PART_COUNT = 0
      DO M = 1, MMAX
         DO L = 1, PART_MPHASE(M) 
            PART_COUNT = PART_COUNT + 1
            DES_RADIUS(PART_COUNT) = D_P0(M)*HALF
            RO_Sol(PART_COUNT) = RO_S(M)
         ENDDO
      ENDDO
      
      DOML(1) = DES_EPS_XSTART
      DOML(2) = DES_EPS_YSTART
      IF(DIMN.EQ.3) THEN
         DOML(3) = DES_EPS_ZSTART
      ENDIF

      IF (DES_EPS_XSTART > XLENGTH) THEN
         WRITE(UNIT_LOG,1001) 'X', 'X'
         WRITE(*,1003)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF (DES_EPS_YSTART > YLENGTH) THEN
         WRITE(UNIT_LOG,1001) 'Y', 'Y'
         WRITE(*,1003)
         CALL MFIX_EXIT(myPE)
      ENDIF
      IF (DIMN .EQ. 3 .AND. DES_EPS_ZSTART > ZLENGTH) THEN
         WRITE(UNIT_LOG,1001) 'Z', 'Z'
         WRITE(*,1003)
         CALL MFIX_EXIT(myPE)
      ENDIF

      
      CALL GENER_LATTICE_MOD(PARTICLES,doml(1:DIMN),&
         DES_POS_OLD(1:PARTICLES,1:DIMN),DES_RADIUS(1:PARTICLES))   

      OPEN(unit=24, file="particle_gener_conf.dat",&
         form="formatted")
      
      DO L = 1, PARTICLES
         OMEGA_OLD(L,:) = ZERO
         DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
         DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:)
         OMEGA_NEW(L,:) = OMEGA_OLD(L,:)
         WRITE(24,'(10(2X,ES12.5))') (DES_POS_OLD(L,K),K=1,DIMN),&
            DES_RADIUS(L),RO_Sol(L), (DES_VEL_OLD(L,K),K=1,DIMN) 
      ENDDO
      CLOSE(24)

      IF(MAXVAL(DES_POS_NEW(1:PARTICLES,2)).GT.&
      YLENGTH-2.d0*MAXVAL(DES_RADIUS(1:PARTICLES))) THEN 
         WRITE(UNIT_LOG,1002) MAXVAL(DES_POS_NEW(1:PARTICLES,2)), &
            YLENGTH-2.d0*MAXVAL(DES_RADIUS(1:PARTICLES))
         WRITE(*,1003)
         CALL MFIX_EXIT(myPE)  
      ENDIF

      WRITE(*,'(3X,A)') &
         '<---------- END GENERATE_PARTICLE_CONFIG ----------'
      
 1001 FORMAT(/1X,70('*')//' From: GENERATE_PARTICLE_CONFIG',/,&
         ' Message: DES_EPS_',A1,'START exceeds ',A1, 'LENGTH',/10X,&
         'Particles cannot be seeded outside the simulation ', &
         'domain',/1X,70('*')/)

 1002 FORMAT(/1X,70('*')//' From: GENERATE_PARTICLE_CONFIG',/,&
         ' Message: Positive overlap with walls in y-dir. Max. ',&
         'y-position of',/10X, 'particle (=', G12.5, &
         ') > YLENGTH-DMAX = ', G12.5,/10X, 'This may occur if',&
         'starting with close packing. Increase',/10X, 'the domain ',&
         'length in the y-dir or generate the particle',/10X,&
         'configuration in a bigger box and shrink it to fit',/10X,&
         'in the desired box size',/1X,70('*')/)

 1003 FORMAT(5X,'An error has occured see the *.LOG FILE for details')

      END SUBROUTINE GENERATE_PARTICLE_CONFIG
      


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE GENER_LATTICE_MOD(NBODY,DOMLIN,XC,DBDY)

      USE discretelement

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------  
! specified domain length for generating particle configuration
      double precision, INTENT(IN) ::  DOMLIN(DIMN)
      double precision :: DOML(DIMN)
! radius of particles/diameter (DES_RADIUS(PARTICLES))
      double precision, INTENT(INOUT) ::  DBDY(NBODY)
! number of particles  (PARTICLES)
      INTEGER, INTENT(IN) :: NBODY
! old position of particles (DES_POS_OLD(PARTICLES,DIMN)) 
      DOUBLE PRECISION, INTENT(OUT) :: xc(nbody,DIMN)
! maximum diameter
      DOUBLE PRECISION DMAX
      double precision :: FAC, XP, YP, ZP
      integer :: nprob1 , i, k, nx,  nz, np1, n
!-----------------------------------------------  

      WRITE(*,'(5X,A)') '---------- START GENER_LATTICE_MOD ---------->'

      DOML(:) = DOMLIN(:)
! convert radius to diameter      
      dbdy(1:NBODY) = 2.d0*dbdy(1:NBODY)
      dmax =  MAXVAL(dbdy(1:nbody))     

      FAC = 1.05
      N = 1
! slightly greater (factor1.05) than the max particle radius      
      YP = dmax*fac/2.d0
      
! loop over N-1 particles      
      DO WHILE (N.LT.NBODY) 

! grid the domain in x and z directions according to 1.05
! times the max particle diameter
         IF(DIMN.EQ.3) THEN 
            nz = floor((real(DOML(3))/(dmax*fac)))
         ELSE 
            nz = 1
         ENDIF
         nx = floor((real(doml(1))/(dmax*fac)))

! fill the domain starting from the bottom:         
         DO K = 1, nz
            ZP =  dmax*half*fac +  (K-1)*dmax*fac
            DO I = 1, nx
               XP = dmax*half*fac + (I-1)*dmax*fac
               XC(N,1) = XP
               XC(N,2) = YP
               IF(DIMN.EQ.3) XC(N,3) = ZP
               IF(N.EQ.NBODY) GOTO 200
               n = n+1
               
            ENDDO
         ENDDO

         YP = YP + dmax*fac
         
 200     CONTINUE 
      ENDDO

      WRITE(*,'(7X,A,I10)') &
         'Number of particles in gener_lattice_mod = ', N
      WRITE(*,'(5X,A)') '<---------- END GENER_LATTICE_MOD ----------'

! convert back to radius
      DBDY(1:NBODY) = 0.5d0*DBDY(1:NBODY)

      END SUBROUTINE GENER_LATTICE_MOD


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
! This routine is no longer called     
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE SET_INITIAL_VELOCITY

      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE toleranc
      USE geometry
      USE indices
      USE is
      USE tau_g
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE interpolation
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER IPJK, IJPK, IJKP, IMJK, IJMK, IJKM, IPJPK, IPJKP, IJPKP, &
              IPJPKP, IJK
      INTEGER I, J, K, II, JJ, KK, IJK_IMAX1, IJK_JMAX1, IJK_KMAX1
      INTEGER IW, IE, JS, JN, KB, KTP

! local variables to temporarily store gas velocities near periodic boundaries      
      DOUBLE PRECISION TEMP_U_G_X(JMAX2, KMAX2),&
         TEMP_V_G_X(JMAX2,KMAX2), TEMP_W_G_X(JMAX2, KMAX2)
      DOUBLE PRECISION TEMP_U_G_Y(IMAX2, KMAX2),&
         TEMP_V_G_Y(IMAX2,KMAX2), TEMP_W_G_Y(IMAX2, KMAX2)
      DOUBLE PRECISION TEMP_U_G_Z(IMAX2, JMAX2),&
          TEMP_V_G_Z(IMAX2,JMAX2), TEMP_W_G_Z(IMAX2, JMAX2)

      DOUBLE PRECISION AVG_FACTOR, AVG_FACTOR2

      INTEGER, DIMENSION(3):: PCELL
      INTEGER ONEW
      INTEGER NP, M
!-----------------------------------------------  
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'

      WRITE(*,'(7X,A)') &
         '---------- START SET_INITIAL_VELOCITY ---------->'
      WRITE(*,'(9X,A,/,A)') &
         'Setting the initial particle velocity equal ',&
         'to the local fluid velocity'
      
      IF (DES_PERIODIC_WALLS) THEN
         IF(DES_PERIODIC_WALLS_X) THEN 
            I = 1   ! ghost cell
            DO J = 1, JMAX2
               DO K = 1, KMAX2
                  IJK = funijk(I,J,K)
                  IJK_IMAX1 = funijk(IMAX1,J,K)   ! imax1=fluid cell
                  TEMP_U_G_X(J,K) = U_G(IJK) 
                  U_G(IJK) = U_G(IJK_IMAX1)
                  TEMP_V_G_X(J,K) = V_G(IJK) 
                  V_G(IJK) = V_G(IJK_IMAX1)
                  IF(DIMN.EQ.3) THEN
                     TEMP_W_G_X(J,K) = W_G(IJK) 
                     W_G(IJK) = W_G(IJK_IMAX1)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         IF(DES_PERIODIC_WALLS_Y) THEN
            J = 1   ! ghost cell
            DO I = 1, IMAX2
               DO K = 1, KMAX2
                  IJK = funijk(I,J,K)
                  IJK_JMAX1 = funijk(I,JMAX1,K)   ! jmax1=fluid cell
                  TEMP_U_G_Y(I,K) = U_G(IJK) 
                  U_G(IJK) = U_G(IJK_JMAX1)
                  TEMP_V_G_Y(I,K) = V_G(IJK) 
                  V_G(IJK) = V_G(IJK_JMAX1)
                  IF(DIMN.EQ.3) THEN 
                     TEMP_W_G_Y(I,K) = W_G(IJK) 
                     W_G(IJK) = W_G(IJK_JMAX1)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF 
         IF(DES_PERIODIC_WALLS_Z .AND. DIMN .EQ. 3) THEN
            K = 1   ! ghost cell
            DO J = 1, JMAX2
               DO I = 1, IMAX2
                  IJK = funijk(I,J,K)
                  IJK_KMAX1 = funijk(I,J,KMAX1)   ! kmax1=fluid cell
                  TEMP_U_G_Z(I,J) = U_G(IJK) 
                  U_G(IJK) = U_G(IJK_KMAX1)
                  TEMP_V_G_Z(I,J) = V_G(IJK) 
                  V_G(IJK) = V_G(IJK_KMAX1)
                  TEMP_W_G_Z(I,J) = W_G(IJK) 
                  W_G(IJK) = W_G(IJK_KMAX1)
               ENDDO
            ENDDO
         ENDIF 
      ENDIF   ! end if des_periodic_wallsF

! sets several quantities including interp_scheme, scheme, and order and
! allocates arrays necessary for interpolation
      CALL SET_INTERPOLATION_SCHEME(2)

! avg_factor =0.25  (in 3D) or =0.50 (in 2D)  
! avg_factor2=0.125 (in 3D) or =0.25 (in 2D)  
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)
      AVG_FACTOR2 = 0.125D0*(DIMN-2) + 0.25D0*(3-DIMN)
      
      DO NP = 1, PARTICLES

         I = PIJK(NP, 1)
         J = PIJK(NP, 2)
         K = PIJK(NP, 3)
         
         PCELL(1) = I-1
         PCELL(2) = J-1
         PCELL(3) = (3-DIMN)*1+(DIMN-2)*(K-1)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.         
         CALL SET_INTERPOLATION_STENCIL(PCELL, IW, IE, JS, JN, KB, &
            KTP, INTERP_SCHEME, DIMN, ORDERNEW = ONEW) 

         DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
            DO J = 1, ONEW
               DO I = 1, ONEW
                  II = IW+I-1
                  JJ = JS+J-1
                  KK = KB+K-1

                  GSTENCIL(I,J,K,1) = XE(II)
                  GSTENCIL(I,J,K,2) = YN(JJ)
! =K-1 (in 3D) or =1 (in 2D)
                  GSTENCIL(I,J,K,3) = ZT(KK)*(DIMN-2) + DZ(1)*(3-DIMN)

! adjust for periodic boundaries
                  IF (DES_PERIODIC_WALLS) THEN
                     IF (DES_PERIODIC_WALLS_X) THEN
                        IF(II.LT.1)     II = IMAX1+II-1
                        IF(II.GT.IMAX1) II = II-IMAX1+1
                     ENDIF
                     IF (DES_PERIODIC_WALLS_Y) THEN
                        IF(JJ.LT.1)     JJ = JMAX1+JJ-1
                        IF(JJ.GT.JMAX1) JJ = JJ-JMAX1+1
                     ENDIF
                     IF (DIMN.EQ.3 .AND. DES_PERIODIC_WALLS_Z) THEN
                        IF(KK.LT.1)     KK = KMAX1+KK-1
                        IF(KK.GT.KMAX1) KK = KK-KMAX1+1
                     ENDIF
                  ENDIF

                  IJK = FUNIJK(II,JJ,KK)
                  
                  IPJK = IP_OF (IJK)    
                  IJPK = JP_OF (IJK)
                  IJKP = KP_OF (IJK)
                  IJPKP = KP_OF(IJPK)
                  IPJKP = KP_OF(IPJK)
                  IPJPK = JP_OF(IPJK)
                  IPJPKP = KP_OF(IPJPK)                  
                  
                  vstencil(i,j,k,1) = AVG_FACTOR*( u_g(ijk) + u_g(ijpk) &
                     + (u_g(ijkp) + u_g(ijpkp))*(DIMN-2) )
                  
                  vstencil(i,j,k,2) = AVG_FACTOR*( v_g(ijk) + v_g(ipjk)&
                     + (v_g(ijkp) + v_g(ipjkp))*(DIMN-2) )
                  
                  IF(DIMN.EQ.3) THEN 
                     vstencil(i,j,k,3) = AVG_FACTOR*(w_g(ijk) +&
                        w_g(ijpk) + w_g(ipjk) + w_g(ipjpk) ) 
                  ELSE 
! doesn't matter what value is put here                          
                     vstencil(i,j,k,3) = 0.d0 
                  ENDIF                  
                  
               ENDDO
            ENDDO
         ENDDO
         
         IF(DIMN.EQ.2) THEN
            CALL interpolator( gstencil(1:onew,1:onew,1,1:DIMN), &
               vstencil(1:onew,1:onew,1,1:2), &
               DES_POS_NEW(np,1:2), vel_fp(np,1:2), &
               onew, interp_scheme, weightp)
         ELSE
            CALL interpolator( gstencil(1:onew,1:onew,1:onew,1:DIMN), &
               vstencil(1:onew,1:onew,1:onew,1:DIMN), &
               DES_POS_NEW(np, 1:3),vel_fp(np,1:3), &
               onew, interp_scheme, weightp )   
         ENDIF
         
         DES_VEL_NEW(NP,:) = VEL_FP(NP,:)
         DES_VEL_OLD(NP,:) = VEL_FP(NP,:)
         
      ENDDO
      
! reset the velocities which were previously adjusted for periodic boundaries
      IF (DES_PERIODIC_WALLS) THEN      
         IF(DES_PERIODIC_WALLS_X) THEN 
            I = 1
            DO J = 1, JMAX2
               DO K = 1, KMAX2
                  IJK = funijk(I,J,K)
                  U_G(IJK) = TEMP_U_G_X(J,K)
                  V_G(IJK) = TEMP_V_G_X(J,K)
                  IF(DIMN.EQ.3) W_G(IJK) = TEMP_W_G_X(J,K)
               ENDDO
            ENDDO
         ENDIF
         IF(DES_PERIODIC_WALLS_Y) THEN 
            J = 1
            DO K = 1, KMAX2
               DO I = 1, IMAX2
                  IJK = funijk(I,J,K)
                  U_G(IJK) = TEMP_U_G_Y(I,K)
                  V_G(IJK) = TEMP_V_G_Y(I,K)
                  IF(DIMN.EQ.3) W_G(IJK) = TEMP_W_G_Y(I,K)
               ENDDO
            ENDDO
         ENDIF 
         IF(DES_PERIODIC_WALLS_Z.AND.DIMN.EQ.3) THEN 
            K = 1
            DO J = 1, JMAX2
               DO I = 1, IMAX2
                  IJK = funijk(I,J,K)
                  U_G(IJK) = TEMP_U_G_Z(I,J)
                  V_G(IJK) = TEMP_V_G_Z(I,J)
                  W_G(IJK) = TEMP_W_G_Z(I,J)
               ENDDO
            ENDDO
         ENDIF 
      ENDIF   ! end if des_periodic_walls


      WRITE(*,'(7X,A)')&
         '<---------- END SET_INITIAL_VELOCITY ----------'   

      END SUBROUTINE SET_INITIAL_VELOCITY



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE writeic
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER:: NP   
!-----------------------------------------------

      OPEN(25,file='ic.dat',form='formatted')

      DO NP = 1,particles
         IF (DIMN .EQ. 2) THEN
            WRITE(25,32) DES_POS_NEW(NP,1:2), DES_VEL_NEW(NP,1:2)
         ELSE
            WRITE(25,32) DES_POS_NEW(NP,1:3), DES_VEL_NEW(NP,1:3)
         ENDIF
      ENDDO

      CLOSE(25)

 32   FORMAT(10(2XE17.8))

      END SUBROUTINE writeic



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C      
      SUBROUTINE init_particles_jn
      USE constant 
      USE geometry
      USE randomno
      USE discretelement
      IMPLICIT NONE 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I, J, K, L
      LOGICAL FILE_EXIST
      REAL*8 :: umf0(dimn), rsf(DIMN, DIMN)
! Quantities for LE BC      
! local variable for shear direction
      CHARACTER*4 SHEAR_DIR      
! shear rate
      DOUBLE PRECISION SHEAR_RATE
! distance between shear boundaries      
      DOUBLE PRECISION SHEAR_LENGTH
! temporary variable - velocity of particle based on position between
! shearing boundaries and shear rate
      DOUBLE PRECISION TMP_VEL
! indices for position and velocity
      INTEGER POSI, VELI
!-----------------------------------------------      

      WRITE(*,'(3X,A)') '---------- START INIT_PARTICLES_JN ---------->'
     
      IF (.NOT.DES_LE_BC) THEN
         WRITE(*,'(5X,A)') 'Initializing normal velocity distribution:'
         WRITE(*,'(5X,A,ES17.8,A,ES17.8)') 'mean = ', pvel_mean,&
            ' and standard deviation = ', PVEL_StDev
      ELSE
         SHEAR_DIR = TRIM(DES_LE_SHEAR_DIR)
         IF(SHEAR_DIR.EQ.'DUDY') THEN
            SHEAR_LENGTH = YLENGTH
            POSI = 2
            VELI = 1
         ELSEIF(SHEAR_DIR.EQ.'DWDY') THEN
            SHEAR_LENGTH = YLENGTH
            POSI = 2
            VELI = 3
         ELSEIF(SHEAR_DIR.EQ.'DVDX') THEN
            SHEAR_LENGTH = XLENGTH
            POSI = 1
            VELI = 2
         ELSEIF(SHEAR_DIR.EQ.'DWDX') THEN
            SHEAR_LENGTH = XLENGTH
            POSI = 1
            VELI = 3
         ELSEIF(SHEAR_DIR.EQ.'DUDZ') THEN
            SHEAR_LENGTH = ZLENGTH 
            POSI = 3
            VELI = 1
         ELSEIF(SHEAR_DIR.EQ.'DVDZ') THEN
            SHEAR_LENGTH = ZLENGTH 
            POSI = 3
            VELI = 2
         ENDIF  
         SHEAR_RATE =  (2.d0*DES_LE_REL_VEL)/SHEAR_LENGTH
         pvel_mean = 0.0d0
         PVEL_StDev = DABS(SHEAR_RATE)
         WRITE(*,'(5X,A,A)') 'Setting up velocity profile consistent',&
            'with shear'
         WRITE(*,'(5X,A,ES17.8,A,ES17.8)') 'mean = ', pvel_mean,&
            ' and standard deviation = ', PVEL_StDev
      ENDIF

      DO J=1,DIMN
         umf0(j)=pvel_mean
         DO I=1,DIMN
            IF(I.EQ.J) THEN
               rsf(I,J)=PVEL_StDev
            ELSE
               rsf(I,J)=0.0
            ENDIF
         ENDDO
         WRITE(*,'(5X,A,I5,2X,A)') 'FOR DIRECTION= ', J,&
            ' where (1=X,2=Y,3=Z):   '
         CALL nor_rno(DES_VEL_OLD(1:PARTICLES,J),umf0(J),rsf(J,J))
      ENDDO


! Adjust initial condition: change position and velocity according to
! shear direction and rate
      IF (DES_LE_BC) THEN
         DO L = 1, PARTICLES
            DES_VEL_OLD(L,VELI) = SHEAR_RATE*DES_POS_OLD(L,POSI)-&
               DES_LE_REL_VEL
!            TMP_VEL = SHEAR_RATE*DES_POS_OLD(L,POSI)-&
!               DES_LE_REL_VEL
!            DES_VEL_OLD(L,VELI) = DES_VEL_OLD(L,VELI) +  TMP_VEL
         ENDDO
      ENDIF
      
      DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)

! updating/writing initial particle configuration files      
      IF (GENER_PART_CONFIG) THEN
         INQUIRE(FILE='particle_gener_conf.dat',exist=FILE_EXIST)
         IF (FILE_EXIST) THEN
            OPEN(UNIT=24,FILE='particle_gener_conf.dat',&
                 STATUS='REPLACE')
            DO L = 1, PARTICLES
               WRITE(24,'(10(X,ES12.5))')&
                  (DES_POS_OLD(L,K),K=1,DIMN), DES_RADIUS(L),&
                  RO_Sol(L), (DES_VEL_OLD(L,K),K=1,DIMN) 
            ENDDO
            CLOSE(24)
         ENDIF
      ELSE
         OPEN(UNIT=24,FILE='particle_input2.dat',&
              STATUS='REPLACE')
         DO L = 1, PARTICLES
            WRITE(24,'(10(X,ES15.5))')&
               (DES_POS_OLD(L,K),K=1,DIMN), DES_RADIUS(L),&
               RO_Sol(L), (DES_VEL_OLD(L,K),K=1,DIMN) 
         ENDDO
         CLOSE(24)
      ENDIF

      WRITE(*,'(3X,A)') '<---------- END INIT_PARTICLES_JN ----------'


      END SUBROUTINE init_particles_jn
