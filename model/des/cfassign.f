!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CFASSIGN                                                C
!                                                                      C
!  Purpose: Assign the necessary values for DEM computation. For       C
!           example:                                                   C
!     - assigning DEM boundaries from the values entered for           C
!       MFIX input in mfix.datat                                       C
!     - assigning DEM gravity vector from MFIX input.                  C
!     - calculating DTSOLID based on particle properties: spring       C
!       coefficient, damping factor & mass                             C
!                                                                      C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 25-Mar-14  C
!  Comments: Breaking this subroutine into several subroutines for DEM C
!            and PIC models                                            C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE constant, only : GRAVITY_X, GRAVITY_Y, GRAVITY_Z
      USE discretelement 
      USE mfix_pic
      use error_manager
! Flag: DEM solids present.
      use run, only: DEM_SOLIDS
! Runtime flag specifying MPPIC solids
      use run, only: PIC_SOLIDS     

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: I, J, K, L
      INTEGER :: IJK, M  ! needed for calling bfx_s, etc
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------

      CALL INIT_ERR_MSG("CFASSIGN")

      
! The common assignments are done in this routine. 
! The model spcific assignmets are moved to the specific subroutines 

      GRAV(1) = GRAVITY_X
      GRAV(2) = GRAVITY_Y
      IF(DIMN.EQ.3) GRAV(3) = GRAVITY_Z
      WRITE(ERR_MSG, '(A, 2x, 3(ES15.7))') 'GRAVITY components for DES:', & 
      GRAVITY_X, GRAVITY_Y, GRAVITY_Z

      CALL FLUSH_ERR_MSG 
      
      IF(DEM_SOLIDS) CALL cfassign_dem
      IF(PIC_SOLIDS) CALL cfassign_pic

! Finalize the error manager.
      CALL FINL_ERR_MSG

      RETURN
      END SUBROUTINE CFASSIGN



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CFASSIGN_PIC                                            C
!                                                                      C
!  Purpose: PIC related cfassign source code moved here                C
!           example:                                                   C
!     - calculating DTSOLID based on particle response time            C
!                                                                      C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 25-Mar-14  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN_PIC

      
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE constant, only: Pi
      USE discretelement, only: des_mmax, dtsolid 
      USE discretelement, only: DES_RO_s, DES_D_P0
      USE physprop, only: MU_g0
      USE mfix_pic, only : dtpic_taup, des_tau_p
      use error_manager
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: M 

      CALL INIT_ERR_MSG("CFASSIGN_PIC")

      DO M = 1, DES_MMAX
         DES_TAU_P(M) = DES_RO_S(M)*(DES_D_P0(M)**2.d0)/(18.d0*MU_g0)
         WRITE(err_msg,'(/A,I2,A,G17.8)') &
         'TAU_P FOR ', M,'th SOLID PHASE= ', DES_TAU_P(M)
         CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)

      ENDDO

      DTSOLID = MINVAL(DES_TAU_P(1:DES_MMAX))
      DTPIC_TAUP = DTSOLID      !maximum dt for point-particles based on taup
      
      WRITE(ERR_MSG, 1000) DTSolid 
      CALL FLUSH_ERR_MSG(Header = .false.)
      
 1000 format('MPPIC: Point-particle ',&
      'approximation for particle-particle and ', /, & 
      'particle-wall collisions', /, &
      'DTSOLID based on particle time response Taup', /, & 
      'DTSOLID = ', 2x, E17.10)

! Finalize the error manager.
      CALL FINL_ERR_MSG

      END SUBROUTINE CFASSIGN_PIC 




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CFASSIGN_DEM                                            C
!                                                                      C
!  Purpose: DEM related cfassign source code moved here                C
!           example:                                                   C
!     - calculating DTSOLID based on particle properties: spring       C
!       coefficient, damping factor & mass                             C
!                                                                      C
!                                                                      C
!  Reviewer: Rahul Garg                               Date: 25-Mar-14  C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN_DEM 

      
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param1
      USE constant, only: Pi
      USE discretelement
      use error_manager
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: I, J, K, L, M 
      INTEGER :: COUNT_E
      DOUBLE PRECISION :: MASS_I, MASS_J, &
                          MASS_EFF, RED_MASS_EFF
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
! local variables for calculation of hertzian contact parameters
      DOUBLE PRECISION :: R_EFF, E_EFF, G_MOD_EFF

      CALL INIT_ERR_MSG("CFASSIGN_DEM")

      TCOLL = LARGE_NUMBER

      WRITE(ERR_MSG,'(A)') 'Setting collision model parameters for DEM'
      CALL FLUSH_ERR_MSG (Footer = .false.)
      
      IF (DES_COLL_MODEL_ENUM == HERTZIAN) THEN

         WRITE(ERR_MSG,'(A)') 'COLLISION MODEL: Hertzian'
         CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)

! particle-particle contact -------------------->
         DO I=1,DES_MMAX
            G_MOD(I) = 0.5d0*e_young(I)/(1.d0+v_poisson(I)) ! shear modulus 
            WRITE(err_msg,'(A,I5,X,A,X,2(ES15.7))') &
            'E_YOUNG AND V_POISSON FOR M = ', I, '=',&
            E_YOUNG(I), V_POISSON(I)
            CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)

         ENDDO

         COUNT_E = 0
         DO I=1,DES_MMAX
            DO J=I,DES_MMAX
! Arrange the coefficient of restitution matrix from en_input values
! use coef of rest to determine damping coefficient
               COUNT_E = COUNT_E + 1
               REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
               REAL_ET(I,J) = DES_ET_INPUT(COUNT_E)
               MASS_I = (PI/6.d0)*(DES_D_P0(I)**3)*DES_RO_S(I)
               MASS_J = (PI/6.d0)*(DES_D_P0(J)**3)*DES_RO_S(J)
               MASS_EFF = (MASS_I*MASS_J)/(MASS_I+MASS_J)
! In the Hertzian model introduce a factor of 2/7 to the effective mass
! for tangential direction to get a reduced mass.  See reference:
! Van der Hoef et al., Advances in Chemical Engineering, 2006, 31, 65-149
!   (see page 94-95)
               RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF
               R_EFF = 0.5d0*(DES_D_P0(I)*DES_D_P0(J)/&
                  (DES_D_P0(I)+DES_D_P0(J)))
               E_EFF = e_young(I)*e_young(J)/ &
               (e_young(I)*(1.d0-v_poisson(J)**2)+&
               e_young(J)*(1.d0-v_poisson(I)**2))
               G_MOD_EFF = G_MOD(I)*G_MOD(J)/&
               (G_MOD(I)*(2.d0-v_poisson(J))+&
               G_MOD(J)*(2.d0-v_poisson(I)))

               hert_kn(I,J)=(4.d0/3.d0)*E_EFF*SQRT(R_EFF)
               hert_kt(I,J)=(16.d0/3.d0)*G_MOD_EFF*SQRT(R_EFF)

               IF (REAL_EN(I,J) .NE. ZERO) THEN
                  DES_ETAN(I,J) = 2.d0*SQRT(hert_kn(I,J)*MASS_EFF)*&
                  ABS(LOG(REAL_EN(I,J)))
                  DES_ETAN(I,J) = DES_ETAN(I,J)/&
                  SQRT(PI*PI + (LOG(REAL_EN(I,J)))**2)
               ELSE
                  DES_ETAN(I,J) = 2.d0*SQRT(hert_kn(I,J)*MASS_EFF)
               ENDIF
               IF (REAL_ET(I,J) .NE. ZERO) THEN
                  DES_ETAT(I,J) = 2.d0*SQRT(hert_kt(I,J)*RED_MASS_EFF)*&
                  ABS(LOG(REAL_ET(I,J)))
                  DES_ETAT(I,J) = DES_ETAT(I,J)/&
                  SQRT(PI*PI + (LOG(REAL_ET(I,J)))**2)
               ELSE
                  DES_ETAT(I,J) = 2.d0*SQRT(hert_kt(I,J)*RED_MASS_EFF)
               ENDIF
               hert_kn(J,I) = hert_kn(I,J)
               hert_kt(J,I) = hert_kt(I,J)

               TCOLL_TMP = PI/SQRT(hert_kn(I,J)/MASS_EFF - ((DES_ETAN(I,J)/MASS_EFF)**2)/4.d0)
               TCOLL = MIN(TCOLL_TMP, TCOLL)

               
               WRITE(ERR_MSG,'(A,I5,X,I5,X,A,X,2(ES15.7))') &
               'KN AND KT FOR PAIR ',I, J, '=', &
               hert_kn(I,J), hert_kt(I,J)
               CALL FLUSH_ERR_MSG

            ENDDO
         ENDDO

! particle-wall contact -------------------->
         COUNT_E = 0
         DO I = 1, DES_MMAX
            COUNT_E = COUNT_E + 1
            REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
            REAL_ET_WALL(I) = DES_ET_WALL_INPUT(COUNT_E)
            MASS_I = (PI/6.d0)*(DES_D_P0(I)**3)*DES_RO_S(I)
            MASS_J = MASS_I
            MASS_EFF = MASS_I
            RED_MASS_EFF = (2.d0/7.d0)*MASS_I
            R_EFF = 0.5d0*DES_D_P0(I)
            E_EFF = e_young(I)*ew_young/ &
            (e_young(I)*(1.d0-vw_poisson**2)+&
            ew_young*(1.d0-v_poisson(I)**2))
            G_MOD_EFF = G_MOD(I)/(2.d0-v_poisson(I))

            hert_kwn(I) =(4.d0/3.d0)*E_EFF*SQRT(R_EFF)
            hert_kwt(I) = (16.d0/3.d0)*G_MOD_EFF*SQRT(R_EFF)

            DES_ETAN_WALL(I) = 2.d0*SQRT(hert_kwn(I)*MASS_EFF)*&
            ABS(LOG(REAL_EN_WALL(I)))
            DES_ETAN_WALL(I) = DES_ETAN_WALL(I)/&
            SQRT(PI*PI + (LOG(REAL_EN_WALL(I)))**2)
            DES_ETAT_WALL(I) = 2.d0*SQRT(hert_kwt(I)*RED_MASS_EFF)*&
            ABS(LOG(REAL_ET_WALL(I)))
            DES_ETAT_WALL(I) = DES_ETAT_WALL(I)/&
            SQRT(PI*PI + (LOG(REAL_ET_WALL(I)))**2)

            TCOLL_TMP = PI/SQRT(hert_kwn(i)/MASS_EFF - ((DES_ETAN_WALL(I)/MASS_EFF)**2)/4.d0) 
         ENDDO

      ELSE                      ! Linear spring-dashpot model

         WRITE(ERR_MSG,'(A)') &
         'COLLISION MODEL: Linear Spring-Dashpot (default)'
            
         CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)

! User's input for KT_FAC and KT_W_FAC will be used, otherwise these values are
! estimated using set factors.  See following references:
!   Schafer et al., J. Phys. I France, 1996, 6, 5-20 (see page 7&13), or
!   Van der Hoef et al., Advances in Chemical Engineering, 2006, 31, 65-149,
! (see page 94-95), or
!   Silbert et al., Physical Review E, 2001, 64, 051302 1-14
! (see page 051302-5)
         IF(KT_FAC == UNDEFINED) THEN
! in LSD model a factor of 2/7 makes period of tangential and normal
! oscillation equal for uniform spheres when en=1 (no dissipation)
            KT = (2.d0/7.d0)*KN
         ELSE
            KT = KT_FAC*KN
         ENDIF
         IF(KT_W_FAC == UNDEFINED) THEN
            KT_W = (2.d0/7.d0)*KN_W
         ELSE
            KT_W = KT_W_FAC*KN_W
         ENDIF
         WRITE(err_msg,'(A,ES17.10,2X,ES15.7)') &
         'KN AND KT = ', KN, KT
         CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)


! particle-particle contact -------------------->
         COUNT_E = 0
         DO I = 1, DES_MMAX
            DO J = I, DES_MMAX
! Arrange the coefficient of restitution matrix from en_input values
! use coef of rest to determine damping coefficient
               COUNT_E = COUNT_E + 1
               REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
               MASS_I = (PI/6.d0)*(DES_D_P0(I)**3.d0)*DES_RO_S(I)
               MASS_J = (PI/6.d0)*(DES_D_P0(J)**3.d0)*DES_RO_S(J)
               MASS_EFF = (MASS_I*MASS_J)/(MASS_I + MASS_J)

               IF (REAL_EN(I,J) .NE. ZERO) THEN
                  DES_ETAN(I,J) = 2.D0*SQRT(KN*MASS_EFF)*&
                  ABS(LOG(REAL_EN(I,J)))
                  DES_ETAN(I,J) = DES_ETAN(I,J)/&
                  SQRT(PI*PI + (LOG(REAL_EN(I,J)))**2)
               ELSE
                  DES_ETAN(I,J) = 2.D0*SQRT(KN*MASS_EFF)
               ENDIF

! User's input for DES_ETAT_FAC will be used, otherwise these values are
! estimated using set factors: See following reference:
!   Silbert et al., Physics of Fluids, 2003, 15, no. 1, 1-10 (see page 3)
               IF(DES_ETAT_FAC == UNDEFINED) THEN
                  DES_ETAT(I,J) = HALF*DES_ETAN(I,J)
               ELSE
                  DES_ETAT(I,J) = DES_ETAT_FAC*DES_ETAN(I,J)
               ENDIF

               TCOLL_TMP = PI/SQRT(KN/MASS_EFF - ((DES_ETAN(I,J)/MASS_EFF)**2)/4.d0)
               TCOLL = MIN(TCOLL_TMP, TCOLL)
            ENDDO
         ENDDO

! particle-wall contact -------------------->
         COUNT_E = 0
         DO I = 1, DES_MMAX
            COUNT_E = COUNT_E + 1
            REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
            MASS_I = (PI*(DES_D_P0(I)**3)*DES_RO_S(I))/6.d0
            MASS_J = MASS_I
            MASS_EFF = MASS_I
            IF (REAL_EN_WALL(I) .NE. ZERO) THEN
               DES_ETAN_WALL(I) = 2.d0*SQRT(KN_W*MASS_EFF)*&
               ABS(LOG(REAL_EN_WALL(I)))
               DES_ETAN_WALL(I) = DES_ETAN_WALL(I)/&
               SQRT(PI*PI + (LOG(REAL_EN_WALL(I)))**2)
            ELSE
               DES_ETAN_WALL(I) = 2.D0*SQRT(KN_W*MASS_EFF)
            ENDIF

            IF(DES_ETAT_W_FAC == UNDEFINED) THEN
               DES_ETAT_WALL(I) = HALF*DES_ETAN_WALL(I)
            ELSE
               DES_ETAT_WALL(I) = DES_ETAT_W_FAC*DES_ETAN_WALL(I)
            ENDIF


            TCOLL_TMP = PI/SQRT(KN_W/MASS_EFF - ((DES_ETAN_WALL(I)/MASS_EFF)**2.d0)/4.d0)
         ENDDO

      ENDIF    ! end if/else(des_coll_model_enum)
!-----------------------------------------------------------------<<<


! filling rest of arrays (symmetric matrix)
      DO I = 1, DES_MMAX
         DO J = I, DES_MMAX
            REAL_EN(J, I) = REAL_EN(I,J)
            REAL_ET(J, I) = REAL_ET(J,I)
            DES_ETAN(J,I) = DES_ETAN(I,J)
            DES_ETAT(J,I) = DES_ETAT(I,J)
         ENDDO
      ENDDO

! reporting information to logs
      DO I = 1, DES_MMAX
         DO J = I, DES_MMAX
            WRITE(err_msg, '(A,I10,2X,I10,A,2(ES15.7))') &
            'ETAN AND ETAT FOR PAIR ',&
            I, J, ' = ', DES_ETAN(I,J), DES_ETAT(I,J)
            CALL FLUSH_ERR_MSG (Header = .false., Footer = .false.)
         ENDDO
      ENDDO
      DTSOLID = TCOLL/50.d0

      WRITE(err_msg,'(A,E17.10,2X,E17.10)') &
      'MIN TCOLL AND DTSOLID = ',TCOLL, DTSOLID
      CALL FLUSH_ERR_MSG (Header = .false.)

! Finalize the error manager.
      CALL FINL_ERR_MSG
      END SUBROUTINE CFASSIGN_DEM 

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! subroutine: compute_volume_of_nodes                                      C
! Author: Rahul Garg                                                       C
! Dare: Dec 20, 2011                                                       C
! Purpose:  Due to the staggered grid, the interpolaiton of mean fields    C
! is always first done at the nodes (of the scalar cell) for any quantity. C
! For a quantity like drag force or ep_s, one needs to assign a geometric  C
! volume to a node. In the past this was done on-the-fly in drag_fgs.f.    C
! VCELL was the variable that was used and it used to be incorrecty set to C
! the volume of the scalar cell that the particle belongs to.              C
! This will be ok for uniform grid and will not hold for non-uniform grids.C
! In addition, at the edges (in 2-D, the node volume is only half of       C
! the standard cell volume. And for corner (in 2-D), it is only one fourth.C
! This was also done on-the-fly earlier                                    C
! But again the volume of the cell was used, which was not correct         C
! also not extendable to cut-cell. So this routine computes the geoemetric C
! volume of the nodes.                                                     C
!                                                                          C
! des_vol_node_ratio was defined in order to develop an estimation         C
! algorithm for computing mean fields when Cartesian grid is on.           C
! Subsequently, another algorithm was developed that negated the need for  C
! des_vol_node_ratio.                                                      C
! However, still keeping this variable for now (RG: July 27, 2012)         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE compute_volume_of_nodes

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE cutcell
      implicit none
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! ijk index of fluid grid and corresponding i, j, k indices
      integer :: ijk, iraw, jraw, kraw
! i, j, k and (i+1, j+1, k+1) indices corrected for any
! cyclic ghost boundaries on fluid grid
      INTEGER :: I, J, K, ip, jp, kp
      integer :: ipjk, ijpk, ipjpk, ijkp, ipjkp, ijpkp, ipjpkp
! volume of indicated grid
      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp
      double precision :: vol_node_count, vol_node_actual_count
! weighting factor
      double precision :: avg_factor
! not used?
      double precision :: vol_node_uncorr
! used for reporting information
      integer :: FLUID_IND, CUTCELL_IND, WALL_IND
      character*100 :: filename
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------

      avg_factor = 0.125d0*(DIMN-2) + 0.25d0*(3-DIMN)

! compute the volume at the grid nodes
! grid nodes start from istart2 to iend1
      IF(DIMN.EQ.2) vol_node_count = 4.
      IF(DIMN.EQ.3) vol_node_count = 8.

      DO ijk = ijkstart3,ijkend3
         des_vol_node(ijk) = zero
         iraw = i_of(ijk)
         jraw = j_of(ijk)
         kraw = k_of(ijk)



! start at 1 (ghost cell) and go to last fluid cell. why start at a
! ghost cell and not a fluid cell?
! See below

! Since we are interested in nodes making up the interpolation stencil,
! their numbering goes from 1 to iend1.
! Think of a case with IMAX = 3. Here the nodes where the interpolation will be
! done will run from 1 (=istart2) to 4 (=iend1).
         IF(iraw.LT.istart2 .OR. iraw.GT.iend1) CYCLE
         IF(jraw.LT.jstart2 .OR. jraw.GT.jend1) CYCLE
         IF(kraw.LT.kstart2 .OR. kraw.GT.kend1) CYCLE

! this will force indices of ghost cells on cyclic borders back to
! the corresponding fluid cells. since we are using i, j, k indices and
! not just a composite ijk index we need these to be shifted to account
! for periodicity
         I = imap_c(iraw)
         J = jmap_c(jraw)
         K = kmap_c(kraw)
         IP = imap_c(iraw+1)
         JP = jmap_c(jraw+1)

! using a function like ip_of(ijk) should work the same as getting funijk
! of the shifted i, j, k indices.  however, small differences will
! occur on the 'edges/corners'. so retaining the latter method at this
! time. see j. galvin for discussion
         ipjk = funijk(IP,J,K)
         ijpk = funijk(I,JP,K)
         ipjpk = funijk(IP,JP,K)

! the existing variable vol(ijk) is not used here for cut-cell reasons
! see r. garg for discussion
         vol_ijk   = dx(I) *dy(J) *dz(K)
         vol_ipjk  = dx(IP)*dy(J) *dz(K)
         vol_ijpk  = dx(I) *dy(JP)*dz(K)
         vol_ipjpk = dx(IP)*dy(JP)*dz(K)

         vol_node_uncorr = avg_factor*(vol_ijk + vol_ipjk + vol_ijpk + vol_ipjpk)
         vol_node_actual_count = vol_node_count

         IF(.NOT.FLUID_AT(ijk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijk  = zero
         ENDIF

         IF(.NOT.FLUID_AT(ipjk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjk  = zero
         ENDIF

         IF(.NOT.FLUID_AT(ijpk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijpk  = zero
         ENDIF

         IF(.NOT.FLUID_AT(ipjpk)) THEN
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjpk = zero
         ENDIF

! this will have non-zero values for non-fluid cells at the
! west/south/bottom borders but not for east/north/top borders?
         des_vol_node(ijk) = avg_factor*(vol_ijk + vol_ipjk + &
            vol_ijpk + vol_ipjpk)

         IF (DIMN.EQ.3) THEN
            KP     = kmap_c(kraw+1)
            ijkp   = funijk(I, J, KP)
            ipjkp  = funijk(IP,J, KP)
            ijpkp  = funijk(I, JP,KP)
            ipjpkp = funijk(IP,JP,KP)

            vol_ijkp   = dx(I) *dy(J) *dz(KP)
            vol_ipjkp  = dx(IP)*dy(J) *dz(KP)
            vol_ijpkp  = dx(I) *dy(JP)*dz(KP)
            vol_ipjpkp = dx(IP)*dy(JP)*dz(KP)

            vol_node_uncorr = avg_factor*(vol_node_uncorr + vol_ijkp + &
               vol_ipjkp + vol_ijpkp + vol_ipjpkp)

            IF(.NOT.FLUID_AT(ijkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ijkp   = zero
            ENDIF

            IF(.NOT.FLUID_AT(ipjkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ipjkp  = zero
            ENDIF

            IF(.NOT.FLUID_AT(ijpkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ijpkp  = zero
            ENDIF

            IF(.NOT.FLUID_AT(ipjpkp)) THEN
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ipjpkp = zero
            ENDIF

            des_vol_node(ijk) = des_vol_node(ijk) + avg_factor*&
               (vol_ijkp + vol_ipjpkp + vol_ijpkp + vol_ipjkp)

         ENDIF

         IF(vol_node_actual_count.GT.zero) THEN
            DES_VOL_NODE_RATIO(ijk)  = vol_node_count/vol_node_actual_count
         ELSE
            DES_VOL_NODE_RATIO(ijk) = zero
         ENDIF

      ENDDO   ! end do ijk=ijkstart3,ijkend3

      RETURN
      RETURN
      END SUBROUTINE compute_volume_of_nodes

