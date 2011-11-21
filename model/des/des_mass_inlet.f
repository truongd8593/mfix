!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_MASS_INLET(BCV_I)                                    !
!                                                                      !
!  Purpose:  This routine fills in the necessary information for new   !
!  particles entereing the system.                                     !
!                                                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 13-Jul-09  !
!  Revision: Modified the loop structure to account for parallel vers  !
!            First compute the position of particles and if inside the !
!            processor then add the particles                          !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_MASS_INLET(BCV_I)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop
      use desgrid 
      use mpi_utility 
 
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER IJK             ! Necessary for function.inc
      INTEGER I, J, K         ! Indices 
      INTEGER IP, LL          ! Loop indices
      INTEGER LS              ! Loop counter start value
      INTEGER TMPI            ! Loop counter
      INTEGER NPG             ! Number of particles in cell
      INTEGER NP              ! Particle id. number
      INTEGER BCV_I, BCV      ! Boundary Condition ID
      INTEGER M               ! Mass phase of new particle

! this variable temporarily stores the number ids of all particles 
! in the ijk location of interest
      INTEGER, DIMENSION(:), ALLOCATABLE :: HOLDER
! a random number between 0 and 1
      DOUBLE PRECISION RAND 
! random integer number between 1 and NUMFRAC_LIMIT
      INTEGER RAND_I 
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS      
! size of mesh for grid based search      
      DOUBLE PRECISION SIZEDX, SIZEDY, SIZEDZ       
! Pradeep altering the structure 
      double precision lpar_rad 
      double precision, dimension(dimn) :: lpar_pos
      logical :: lflag 
      integer :: lglobal_id
      logical :: ltouching
      double precision,dimension(dimn):: lrange_min,lrange_max,lrand3
!-----------------------------------------------
      INCLUDE 'function.inc'
      LS = 1
      BCV = DES_BC_MI_ID(BCV_I)

      DO IP = 1, PI_COUNT(BCV_I) !Loop over the particles being injected

! Pradeep - altering the loop structure 
! First the particle is assigned its position and type ; Based on the position
! it is added to the corresponding processor
         if (mype.eq.0) then 
            if(des_bc_poly(bcv_i)) then
               call random_number(rand)
               rand_i = ceiling(dble(numfrac_limit)*rand)
               m = des_bc_poly_layout(bcv_i,rand_i)
            else
               m = des_bc_poly_layout(bcv_i,1)
            endif
            lpar_rad = (d_p0(m) * half)
            call des_place_new_particle(bcv_i,lpar_rad,lpar_pos)
         end if 
         call bcast(lpar_rad)
         call bcast(lpar_pos)

         imax_global_id = imax_global_id + 1
         lglobal_id = imax_global_id 

         lflag = .false. 
         select case (des_mi_class(bcv_i))
         case ('XW','XE', 'YZw','YZe')
            if (dimn .eq. 2) then   
               if (    lpar_pos(2) .ge. yn(jstart1-1) &
                 .and. lpar_pos(2) .lt. yn(jend1)) lflag = .true.   
            else 
               if (    lpar_pos(2) .ge. yn(jstart1-1) &
                 .and. lpar_pos(2) .lt. yn(jend1)     &
                 .and. lpar_pos(3) .ge. zt(kstart1-1) &
                 .and. lpar_pos(3) .lt. zt(kend1)) lflag = .true.  
            end if 
         case ('YN','YS', 'XZn','XZs')
            if (dimn .eq. 2) then   
               if (    lpar_pos(1) .ge. xe(istart1-1) &
                 .and. lpar_pos(1) .lt. xe(iend1)) lflag = .true.   
            else 
               if (    lpar_pos(1) .ge. xe(istart1-1) &
                 .and. lpar_pos(1) .lt. xe(iend1)     &
                 .and. lpar_pos(3) .ge. zt(kstart1-1) &
                 .and. lpar_pos(3) .lt. zt(kend1)) lflag = .true.  
            end if 
         case ('ZT','ZB', 'XYt','XYb')
            if (    lpar_pos(1) .ge. xe(istart1-1) &
              .and. lpar_pos(1) .lt. xe(iend1)     &
              .and. lpar_pos(2) .ge. yn(jstart1-1) &
              .and. lpar_pos(2) .lt. yn(jend1)) lflag = .true.  
         end select 

         if (lflag) then 

            if(particle_plcmnt(bcv_i) == 'RAND')then
! In case of random particle position, check if the inserted particle 
! overlabs any existing particles; if so then adjust the position  
! set range where particles can be inserted  
               select case (des_mi_class(bcv_i))
               case ('XW','XE', 'YZw','YZe')
                  lrange_min(1) = lpar_pos(1);lrange_max(1) = lpar_pos(1) 
                  lrange_min(2) = max(DES_BC_Y_s(BCV),yn(jstart1-1)) + lpar_rad
                  lrange_max(2) = min(DES_BC_Y_n(BCV),yn(jend1)) - lpar_rad
                  if (dimn .eq. 3) then   
                     lrange_min(3) = max(DES_BC_Z_b(BCV),zt(kstart1-1)) + lpar_rad
                     lrange_max(3) = min(DES_BC_Z_t(BCV),zt(kend1)) - lpar_rad
                  end if 
               case ('YN','YS', 'XZn','XZs')
                  lrange_min(1) = max(DES_BC_X_w(BCV),xe(istart1-1)) + lpar_rad
                  lrange_max(1) = min(DES_BC_X_e(BCV),xe(iend1)) - lpar_rad
                  lrange_min(2) = lpar_pos(2); lrange_max(2)=lpar_pos(2)
                  if (dimn .eq. 3) then   
                     lrange_min(3) = max(DES_BC_Z_b(BCV),zt(kstart1-1)) + lpar_rad
                     lrange_max(3) = min(DES_BC_Z_t(BCV),zt(kend1)) - lpar_rad
                  end if 
               case ('ZT','ZB', 'XYt','XYb')
                  lrange_min(1) = max(DES_BC_X_w(BCV),xe(istart1-1)) + lpar_rad
                  lrange_max(1) = min(DES_BC_X_e(BCV),xe(iend1)) - lpar_rad
                  lrange_min(2) = max(DES_BC_Y_s(BCV),yn(jstart1-1)) + lpar_rad
                  lrange_max(2) = min(DES_BC_Y_n(BCV),yn(jend1)) - lpar_rad
                  lrange_min(3) = lpar_pos(3); lrange_max(3)=lpar_pos(3)
               end select 
! if particle does not belong to this range change the position
! this ensures that the particle will be placed at a radius distance 
! from the processor interface, i,e no contact with particles 
! inserted in neighbour proc
               ltouching = .false.
               do li = 1,dimn
                  if (lpar_pos(li).lt.lrange_min(li) .or. &
                      lpar_pos(li).gt.lrange_max(li)) ltouching =.true. 
               end do 
               do while (.true.)
                  if (ltouching) then 
                     call random_number(lrand3(1))
                     call random_number(lrand3(2))
                     if(dimn.eq.3) call random_number(lrand3(3))
                     lpar_pos= lrange_min+(lrange_max-lrange_min)*lrand3
                  end if 
                  call des_new_particle_test(bcv_i,lpar_rad,lpar_pos,ltouching)
                  if(.not.ltouching) exit
              end do 
            end if 

! Check to see if MAX_PIP has been exceeded, if so, STOP
            IF(PIP .GE. MAX_PIP) THEN
               WRITE(UNIT_LOG, 1000)
                WRITE(*,1000)
                CALL MFIX_EXIT(myPE)
            ENDIF 

! Find the first free space in the particle existance array
            DO NP = LS, MAX_PIP
               IF(.NOT.PEA(NP,1)) THEN
                  LS = NP
                  EXIT
               ENDIF
            ENDDO
! Set the flag in the particle existance array
            PEA(NP,1) = .TRUE. 

! Set the flag that the particle is new.  This allows it to be ignored
! by various subroutines
            PEA(NP,2) = .TRUE.
! Set outlet and ghost particle flase 
            PEA(NP,3) = .FALSE.
            PEA(NP,4) = .FALSE.

! Increment the particle in system value by one
            PIP = PIP + 1
! Set the initial velocity values
            DES_VEL_OLD(NP,1) = DES_BC_U_s(BCV)
            DES_VEL_OLD(NP,2) = DES_BC_V_s(BCV)
            IF(DIMN == 3) DES_VEL_OLD(NP,3) = DES_BC_W_s(BCV)
            DES_VEL_NEW(NP,:) = DES_VEL_OLD(NP,:)
! Set the initial angular velocity values
            OMEGA_OLD(NP,:) = 0
            OMEGA_NEW(NP,:) = 0

! Set the particle radius value
            DES_RADIUS(NP) = (D_P0(M) * HALF)

! Set the particle density value
            RO_Sol(NP) = RO_S(M)

! Set the particle mass phase
            PIJK(NP,5) = M

! Calculate the new particle's Volume, Mass, OMOI
            PVOL(NP) = (4.0d0/3.0d0) * PI * DES_RADIUS(NP)**3
            PMASS(NP) = PVOL(NP) * RO_Sol(NP)
            OMOI(NP) = 5.d0 / (2.d0 * PMASS(NP) * DES_RADIUS(NP)**2) 

! Set the initial position values based on mass inlet class
            DES_POS_NEW(NP,:) = lpar_pos
            DES_POS_OLD(NP,:) = lpar_pos 
            iglobal_id(np) = lglobal_id 

            XPOS = DES_POS_NEW(NP,1)
            YPOS = DES_POS_NEW(NP,2)
            IF (DIMN == 3) THEN
               ZPOS = DES_POS_NEW(NP,3)
            ENDIF               

! Determine the i, j, k indices of the cell containing the new 
! particle(s) by checking the cells near the mass inlet using 
! GS_ARRAY. Note the indices will place the particle in a ghost 
! cell. 
! ------------------------------------------------------------
            DO I = GS_ARRAY(BCV_I,1), GS_ARRAY(BCV_I,2)
               IF(XPOS < XE(1))THEN 
                  PIJK(NP,1) = 1
                  EXIT
               ELSEIF(XPOS >= XE(IMAX1))THEN 
                  PIJK(NP,1) = IMAX2
                  EXIT
               ELSEIF((XPOS >= XE(I-1)) .AND. (XPOS < XE(I))) THEN
                  PIJK(NP,1) = I
                  EXIT
               ENDIF
            ENDDO

            DO J = GS_ARRAY(BCV_I,3), GS_ARRAY(BCV_I,4)
               IF(YPOS < YN(1))THEN
                  PIJK(NP,2) = 1
                  EXIT
               ELSEIF(YPOS >= YN(JMAX1))THEN 
                  PIJK(NP,2) = JMAX2
                  EXIT
               ELSEIF((YPOS >= YN(J-1)) .AND. (YPOS < YN(J))) THEN
                  PIJK(NP,2) = J
                  EXIT
               ENDIF
            ENDDO

            IF(DIMN == 2) THEN
               PIJK(NP,3)  = 1
            ELSE
               DO K = GS_ARRAY(BCV_I,5), GS_ARRAY(BCV_I,6)
                  IF(ZPOS < ZT(1))THEN
                     PIJK(NP,3) = 1
                     EXIT
                  ELSEIF(ZPOS >= ZT(KMAX1))THEN 
                     PIJK(NP,3) = KMAX2
                     EXIT               
                  ELSEIF((ZPOS >= ZT(K-1)) .AND. (ZPOS < ZT(K))) THEN 
                     PIJK(NP,3) = K
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

! Update the PIC array for new particles so that any subsequent
! particles that are to be injected will be checked to prevent 
! overlap with previously injected particles 
! pradeep - PIC is now one dimensional for parallel processing 
            I = PIJK(NP,1)
            J = PIJK(NP,2)
            K = PIJK(NP,3)
            IJK = FUNIJK(I,J,K)
            PIJK(NP,4) = IJK

            IF (ASSOCIATED(PIC(IJK)%P)) THEN
               NPG = SIZE(PIC(IJK)%P)
               ALLOCATE( HOLDER (NPG) )

! store the particle no. id of all particles at the ijk location            
               TMPI = 1
               DO LL = 1, NPG
                  HOLDER(TMPI) = PIC(IJK)%P(LL)
                  TMPI = TMPI + 1
               ENDDO
               DEALLOCATE(PIC(IJK)%P)

! essentially increasing the no. of particles at the ijk location by 1
               ALLOCATE(PIC(IJK)%P(TMPI))
               DO LL = 1, TMPI - 1
                  PIC(IJK)%P(LL) = HOLDER(LL)
               ENDDO
! storing the new particle no. id in the list
               PIC(IJK)%P(TMPI) = NP
               DEALLOCATE(HOLDER)
            ELSE
! no other particles were at this ijk location
               TMPI = 1
               ALLOCATE(PIC(IJK)%P(TMPI))
               PIC(IJK)%P(TMPI) = NP
            ENDIF 
            PINC(IJK) = TMPI
         end if 
! All the computation moved to desgri mod. The following code is not 
! required 
!! If using des_neighbor_search option 4 (cell/grid based search) then
!! determine the i,j,k indices of the cell containing the new particle
!! based on the mesh for the grid based search. If cell is outside the 
!! domain then either set the index to 1 or add 2 to the index to 
!! account for ghost cells.  
!! Note that this section is probably unnecessary since the routine
!! particles_in_cell will pickup the same information before neighbor
!! search is ever called and it is only in the neighbor search routine 
!! that this information should be needed
!! ------------------------------------------------------------         
!         IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN
!            SIZEDX = XLENGTH/DESGRIDSEARCH_IMAX
!            SIZEDY = YLENGTH/DESGRIDSEARCH_JMAX
!            IF (DIMN .EQ. 2) THEN
!               SIZEDZ = ONE
!            ELSE
!               SIZEDZ = ZLENGTH/DESGRIDSEARCH_KMAX
!            ENDIF
!
!            IF (XPOS < 0) THEN
!               I = 1
!            ELSEIF (XPOS >= XLENGTH) THEN
!               I = DESGS_IMAX2
!            ELSE         
!               I = INT(XPOS/SIZEDX)+2
!            ENDIF
!
!            IF (YPOS < 0 ) THEN
!               J = 1
!            ELSEIF (YPOS >= YLENGTH) THEN
!               J = DESGS_JMAX2
!            ELSE
!               J = INT(YPOS/SIZEDY)+2
!            ENDIF
!
!            IF (DIMN .EQ. 2) THEN
!               K = 1
!            ELSE
!               IF (ZPOS < 0 ) THEN
!                  K = 1
!               ELSEIF (ZPOS >= ZLENGTH) THEN
!                  K = DESGS_KMAX2
!               ELSE
!                  K = INT(ZPOS/SIZEDZ)+2
!               ENDIF
!            ENDIF
!   
!            DESGRIDSEARCH_PIJK(NP,1) = I
!            DESGRIDSEARCH_PIJK(NP,2) = J
!            DESGRIDSEARCH_PIJK(NP,3) = K
!   
!! Update the DESGRIDSEARCH_PIC array for new particles. Note that 
!! is not needed to prevent subsequent particles that are injected
!! from overlapping with previous particles as the variable PIC is 
!! used for that.
!            IF (ASSOCIATED(DESGRIDSEARCH_PIC(I,J,K)%P)) THEN
!               NPG = SIZE(DESGRIDSEARCH_PIC(I,J,K)%P)
!               ALLOCATE( HOLDER (NPG) )
!
!! store the particle no. id of all particles at the ijk location            
!               TMPI = 1
!               DO LL = 1, NPG
!                  HOLDER(TMPI) = DESGRIDSEARCH_PIC(I,J,K)%P(LL)
!                  TMPI = TMPI + 1
!               ENDDO
!               DEALLOCATE(DESGRIDSEARCH_PIC(I,J,K)%P)
!   
!! essentially increasing the no. of particles at the ijk location by 1
!               ALLOCATE(DESGRIDSEARCH_PIC(I,J,K)%P(TMPI))
!               DO LL = 1, TMPI - 1
!                  DESGRIDSEARCH_PIC(I,J,K)%P(LL) = HOLDER(LL)
!               ENDDO
!! storing the new particle no. id in the list
!               DESGRIDSEARCH_PIC(I,J,K)%P(TMPI) = NP
!               DEALLOCATE(HOLDER)
!            ELSE
!! no other particles were at this ijk location
!               TMPI = 1
!               ALLOCATE(DESGRIDSEARCH_PIC(I,J,K)%P(TMPI))
!               DESGRIDSEARCH_PIC(I,J,K)%P(TMPI) = NP
!            ENDIF 
!         ENDIF   ! end if des_neighbor_search = 4  

      ENDDO   ! end loop over the no. of injected particles (NP)

 1000 FORMAT(/1X,70('*')//,' From: DES_MASS_INLET -',/&
         ' Message: Maximum number of particles in the system MAX_PIS',&
         /10X,' has been exceeded. Increase the value in mfix.dat',/&
         1X,70('*')/)
 
      RETURN
      END SUBROUTINE DES_MASS_INLET



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: DES_PLACE_NEW_PARTICLE                                 !
!                                                                      !
!  Purpose:  This routine uses the classification information to place !
!  a new particle in the proper location.                              !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_PLACE_NEW_PARTICLE(BCV_I,lpar_rad,lpar_pos)

      USE compar
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE param1
      USE physprop

      IMPLICIT NONE
! Dummy Variables 
      double precision lpar_rad
      double precision, dimension(dimn) :: lpar_pos
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! the associated bc no.      
      INTEGER BCV_I, BCV  
! a random number between 0-1
      DOUBLE PRECISION RAND1, RAND2
! 'lower' x,y,z location of current window
      DOUBLE PRECISION MI_WIN
! a random index number for I_OF_MI or J_OF_MI
      INTEGER MI_ORD
!      
      LOGICAL TOUCHING 
!-----------------------------------------------

      TOUCHING = .TRUE.
      BCV = DES_BC_MI_ID(BCV_I)

      SELECT CASE(DES_MI_CLASS(BCV_I))   

! 2D domain with verticle mass inlet: (west face)
! ----------------------------------------
         CASE ('XW')
! Offset Y axis by 1 particle diameter (max diameter) into ghost cell
            lpar_pos(1) = -(DES_BC_OFFSET(BCV_I))
            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  CALL RANDOM_NUMBER(RAND1)
                  lpar_pos(2) = DES_BC_Y_s(BCV) + lpar_rad +&
                     RAND1*(DES_BC_Y_n(BCV)-DES_BC_Y_s(BCV) -&
                     lpar_rad*2.0d0)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window 
                  DES_BC_Y_s(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 2D domain with verticle mass inlet: (east face)
! ----------------------------------------
         CASE ('XE')
! Offset Y axis by 1 particle diameter into ghost cell
            lpar_pos(1) = XLENGTH + (DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(2) = DES_BC_Y_s(BCV) + lpar_rad +&
                     RAND1*(DES_BC_Y_n(BCV)-DES_BC_Y_s(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 2D domain with horizontal mass inlet: (south face)
! ----------------------------------------
         CASE ('YS')
! Offset X axis by 1 particle diameter into ghost cell
            lpar_pos(2) = -(DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(1) = DES_BC_X_w(BCV) + lpar_rad +&
                     RAND1*(DES_BC_X_e(BCV) - DES_BC_X_w(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!                ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 2D domain with horizontal mass inlet: (north face)
! ----------------------------------------
         CASE ('YN')
! Offset X axis by 1 particle diameter into ghost cell
            lpar_pos(2) = YLENGTH + (DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(1) = DES_BC_X_w(BCV) + lpar_rad +&
                     RAND1*(DES_BC_X_e(BCV) - DES_BC_X_w(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) + &   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF


! 3D domain with mass inlet on XZ plane: (south face)
! ----------------------------------------
         CASE ('XZs')
! Offset XZ plane by 1 particle diameter into ghost cell
            lpar_pos(2) = -(DES_BC_OFFSET(BCV_I))
            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(1) = DES_BC_X_w(BCV)+ lpar_rad +&
                     RAND1*(DES_BC_X_e(BCV) - DES_BC_X_w(BCV) -&
                     2.0d0*lpar_rad)
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  lpar_pos(3) = DES_BC_Z_b(BCV) + lpar_rad +&
                     RAND2*(DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
! Determine x-coordinate from ordered grid position I
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on XZ plane: (north face)
! ----------------------------------------
         CASE ('XZn')
! Offset XZ plane by 1 particle diameter into ghost cell
            lpar_pos(2) = YLENGTH + (DES_BC_OFFSET(BCV_I))
            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(1) = DES_BC_X_w(BCV) + lpar_rad +&
                     RAND1*(DES_BC_X_e(BCV) - DES_BC_X_w(BCV) -&
                     2.0d0*lpar_rad)
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  lpar_pos(3) = DES_BC_Z_b(BCV) + lpar_rad +&
                     RAND2*(DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine x-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)                  
               lpar_pos(3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on XY plane: (bottom face)
! ----------------------------------------
         CASE ('XYb')
! Offset XY plane by 1 particle diameter into ghost cell
            lpar_pos(3) = -(DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(1) = DES_BC_X_w(BCV) + lpar_rad +&
                     RAND1*(DES_BC_X_e(BCV) - DES_BC_X_w(BCV) -&
                     2.0d0*lpar_rad)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  lpar_pos(2) = DES_BC_Y_s(BCV) + lpar_rad +&
                     RAND2*(DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine x-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine y-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(2) = RAND2 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
                 IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                    MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
                 ELSE
                    MI_FACTOR(BCV_I) = 1
                 ENDIF
              ENDIF

! 3D domain with mass inlet on XY plane: (top face)
! ----------------------------------------
         CASE ('XYt')
! Offset XY plane by 1 particle diameter into ghost cell
            lpar_pos(3) = ZLENGTH + (DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_X_w(BCV) and DES_BC_X_e(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(1) = DES_BC_X_w(BCV) + lpar_rad +&
                     RAND1*(DES_BC_X_e(BCV) - DES_BC_X_w(BCV) -&
                     2.0d0*lpar_rad)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV)
                  lpar_pos(2) = DES_BC_Y_s(BCV) + lpar_rad +&
                     RAND2*(DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine x-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(1) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_X_w(BCV)   !shift relative to BC location
! Determine y-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(2) = RAND2 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on YZ plane: (west face)
! ----------------------------------------
         CASE ('YZw')
! Offset YZ plane by 1 particle diameter into ghost cell
            lpar_pos(1) = -(DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
! Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(2) = DES_BC_Y_s(BCV) + lpar_rad +&
                     RAND1*(DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV) -&
                     2.0d0*lpar_rad)
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  lpar_pos(3) = DES_BC_Z_b(BCV) + lpar_rad +&
                     RAND2*(DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine y-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) + &   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
               IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                  MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
               ELSE
                  MI_FACTOR(BCV_I) = 1
               ENDIF
            ENDIF

! 3D domain with mass inlet on YZ plane: (east face)
! ----------------------------------------
         CASE ('YZe')
! Offset YZ plane by 1 particle diamter into ghost cell
            lpar_pos(1) = XLENGTH + (DES_BC_OFFSET(BCV_I))

            IF(PARTICLE_PLCMNT(BCV_I) == 'RAND')THEN
!               DO WHILE (TOUCHING)
!  Obtain a random number from (0,1]
                  CALL RANDOM_NUMBER(RAND1)
                  CALL RANDOM_NUMBER(RAND2)
! Place randomly between DES_BC_Y_s(BCV) and DES_BC_Y_n(BCV) and ensure
! that the particle center will be placed a radius away from either edge
! of the inlet
                  lpar_pos(2) = DES_BC_Y_s(BCV) + lpar_rad +&
                     RAND1*(DES_BC_Y_n(BCV) - DES_BC_Y_s(BCV) -&
                     2.0d0*lpar_rad)
! Place randomly between DES_BC_Z_b(BCV) and DES_BC_Z_t(BCV)
                  lpar_pos(3) = DES_BC_Z_b(BCV) + lpar_rad +&
                     RAND2*(DES_BC_Z_t(BCV) - DES_BC_Z_b(BCV) -&
                     2.0d0*lpar_rad)
! Test that no new particles are touching
!                  CALL DES_NEW_PARTICLE_TEST(BCV_I,lpar_rad,lpar_pos,TOUCHING)
!               ENDDO
            ENDIF
            IF(PARTICLE_PLCMNT(BCV_I) == 'ORDR')THEN
! Obtain a random number from (0,1]
               CALL RANDOM_NUMBER(RAND1)
               CALL RANDOM_NUMBER(RAND2)
               MI_ORD = MI_ORDER(BCV_I)%VALUE(MI_FACTOR(BCV_I))
! Determine y-coordinate from ordered grid position I
               MI_WIN = dble(I_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(2) = RAND1 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Y_s(BCV)   !shift relative to BC location
! Determine z-coordinate from ordered grid position J
               MI_WIN = dble(J_OF_MI(BCV_I)%VALUE(MI_ORD))*&
                  MI_WINDOW(BCV_I)
               lpar_pos(3) = RAND2 *&
                  (MI_WINDOW(BCV_I) - (lpar_rad*2.0d0)) +&   !play room
                  MI_WIN + (lpar_rad*2.0d0)*HALF +&   !shift loc. into window
                  DES_BC_Z_b(BCV)   !shift relative to BC location
! Increment MI_FACTOR by one or loop back to one
                 IF(MI_FACTOR(BCV_I) < SIZE(MI_ORDER(BCV_I)%VALUE))THEN
                    MI_FACTOR(BCV_I) = MI_FACTOR(BCV_I) + 1
                 ELSE
                    MI_FACTOR(BCV_I) = 1
                 ENDIF
              ENDIF

         CASE DEFAULT
            PRINT*,'INVALID DES MASS INLET CLASSIFICATION'
            CALL MFIX_EXIT(myPE)
      END SELECT

      RETURN
      END SUBROUTINE DES_PLACE_NEW_PARTICLE


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name:  DES_NEW_PARTICLE_TEST                                 !
!                                                                      !
!  Purpose:  This routine checks if a new particle placed using the    !
!  random inlet was placed in contact with an existing particle.  If   !
!  so a flag is set indicating contact, and the new particle is        !
!  repositioned within the inlet domain.                               !
!                                                                      !
!  Author: J.Musser                                   Date: 14-Aug-09  !
!                                                                      !
!  Purpose: This routine has to be modified for parallel version
!           the parameter now accepts the lpar_rad and lpar_pos and tests
!           if it touches any particles 
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE DES_NEW_PARTICLE_TEST(BCV_I,ppar_rad,ppar_pos,TOUCHING)

      USE compar
      USE constant
      USE des_bc
      USE discretelement
      USE funits
      USE geometry
      USE indices
      USE param1
      USE physprop

      IMPLICIT NONE
! Dummy variables
! index of boundary condition 
      INTEGER BCV_I
      double precision, dimension(dimn) :: ppar_pos 
      double precision :: ppar_rad
      LOGICAL TOUCHING
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! particle number id of a potential overlapping/contacting particle
      INTEGER NP2      
! total number of particles in current ijk cell and loop counter
      INTEGER NPG, LL
! i, j, k indices along boundary used for loop counters
      INTEGER I, J, K, IJK
! for parallel processing 
      integer listart,liend,ljstart,ljend,lkstart,lkend 

      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
      DOUBLE PRECISION  DISTVEC(DIMN), DIST, R_LM
!-----------------------------------------------

      INCLUDE 'function.inc'

      TOUCHING = .FALSE.

! For parallel processing the arrays has to be limited 
      select case (des_mi_class(bcv_i))
      case ('XW','XE', 'YZw','YZe')
         listart = gs_array(bcv_i,1)
         liend = gs_array(bcv_i,2)
         ljstart = max(gs_array(bcv_i,3),jstart)
         ljend = min(gs_array(bcv_i,4),jend)
         lkstart = max(gs_array(bcv_i,5),jstart)
         lkend = min(gs_array(bcv_i,6),jend)
      case ('YN','YS', 'XZn','XZs')
         listart = max(gs_array(bcv_i,1),istart)
         liend = min(gs_array(bcv_i,2),iend)
         ljstart = gs_array(bcv_i,3)
         ljend = gs_array(bcv_i,4)
         lkstart = max(gs_array(bcv_i,5),jstart)
         lkend = min(gs_array(bcv_i,6),jend)
      case ('ZT','ZB', 'XYt','XYb')
         listart = max(gs_array(bcv_i,1),istart)
         liend = min(gs_array(bcv_i,2),iend)
         ljstart = max(gs_array(bcv_i,3),jstart)
         ljend = min(gs_array(bcv_i,4),jend)
         lkstart = gs_array(bcv_i,5)
         lkend = gs_array(bcv_i,6)
      end select 
    
      DO k = lkstart,lkend 
      DO j = ljstart,ljend 
      DO i = listart,liend 
!      DO K = GS_ARRAY(BCV_I,5), GS_ARRAY(BCV_I,6)
!         DO J = GS_ARRAY(BCV_I,3), GS_ARRAY(BCV_I,4)
!           DO I =  GS_ARRAY(BCV_I,1), GS_ARRAY(BCV_I,2)
             IJK = FUNIJK(I,J,K)
             IF(ASSOCIATED(PIC(IJK)%P)) THEN
               NPG =  SIZE(PIC(IJK)%P)                     
               DO LL = 1, NPG
                  NP2 = PIC(IJK)%P(LL)
                  DISTVEC(:) = ppar_pos(:) - DES_POS_NEW(NP2,:)
                  DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))
                  R_LM = ppar_rad + DES_RADIUS(NP2)
                  IF(DIST .LE. R_LM) TOUCHING = .TRUE.
               ENDDO
             ENDIF
           ENDDO
         ENDDO
       ENDDO

      RETURN
      END SUBROUTINE DES_NEW_PARTICLE_TEST
