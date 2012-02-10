!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
! Module name: CFASSIGN                                                C
!
! Purpose:
! Assign the necessary values for DEM  computation. For example:
! - assigning DEM boundaries from the values entered for MFIX input
!   in mfix.datat
! - assigning DEM gravity vector from MFIX input. 
! - calculating DTSOLID based on particle properties: spring 
!   coefficient, damping factor & mass 
!      
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CFASSIGN

!-----------------------------------------------
! Modules 
!-----------------------------------------------
      USE param1
      USE physprop
      USE geometry
      USE constant
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: I, J, K, L
      INTEGER :: IJK, M  ! needed for calling bfx_s, etc      
      INTEGER :: COUNT_E
      DOUBLE PRECISION :: MASS_I, MASS_J, &
                          MASS_EFF, RED_MASS_EFF
      DOUBLE PRECISION :: TCOLL, TCOLL_TMP
! local variables for calculation of hertzian contact parameters
      DOUBLE PRECISION :: R_EFF, E_EFF, G_MOD_EFF     
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'b_force1.inc'
      INCLUDE 'b_force2.inc'
!-----------------------------------------------

      IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(3X,A)') &
         '---------- START CFASSIGN ---------->'

! compute the volume of nodes - more description would be nice
      CALL compute_volume_of_nodes

      TCOLL = LARGE_NUMBER

! Set misc quantities defining the system
!----------------------------------------------------------------->>>
! Set boundary edges 
! In some instances wx1,ex2, etc have been used and in others
! xlength,zero, etc are used. the code should be modified for
! consistency throughout      
      WX1 = ZERO 
      EX2 = XLENGTH 
      BY1 = ZERO
      TY2 = YLENGTH 
      SZ1 = ZERO 
      NZ2 = ZLENGTH

! the DEM variable grav(:) will not accomodate a body force that varies
! in space or on phases unlike the implementation in the continuum 
! simulations
      GRAV(1) = BFX_s(1,1)
      GRAV(2) = BFY_s(1,1)
      IF(DIMN.EQ.3) GRAV(3) = BFZ_s(1,1)     

! Note : the quantities xe, zt cannot be readily replaced with the
! similar appearing variables x_e, z_t in main mfix code as they 
! are not the same.  also the variable y_n does not exist in main
! mfix code. 
! Each loop starts at 2 and goes to max+2 (i.e., imin1=2, imax2=imax+2) 
! However, the indices range to include ghost cells (0-imax2) to avoid 
! multiple if statements in particles_in_cell
      XE(IMIN2-1) = ZERO-DX(IMIN2)
      DO I = IMIN2, IMAX2
         XE(I) = XE(I-1) + DX(I)
      ENDDO
      YN(JMIN2-1) = ZERO-DY(JMIN2)
      DO J  = JMIN2, JMAX2
         YN(J) = YN(J-1) + DY(J)
      ENDDO
      IF(DIMN.EQ.3) THEN
         ZT(KMIN2-1) = ZERO-DZ(KMIN2)
         DO K = KMIN2, KMAX2
            ZT(K) = ZT(K-1) + DZ(K)
         ENDDO
      ENDIF

! used in octree and quadtree neighbor search algorithms      
      RADIUS_EQ = MAX_RADIUS*1.05d0
      INQC = INIT_QUAD_COUNT
!-----------------------------------------------------------------<<<


! Calculate collision parameters
!----------------------------------------------------------------->>>
      IF(DMP_LOG.AND..NOT.MPPIC) WRITE(UNIT_LOG,'(/2X,A,A)') &
      'SOFT-SPRING FOR PARTICLE-PARTICLE AND PARTICLE-WALL COLLISIONS'
      
      IF (TRIM(DES_COLL_MODEL) == 'HERTZIAN') THEN 

         IF(DMP_LOG.AND..NOT.MPPIC) &
            WRITE(UNIT_LOG,'(5X,A)') 'COLLISION MODEL: Hertzian'

! particle-particle contact -------------------->
         DO I=1,MMAX
            G_MOD(I) = 0.5d0*e_young(I)/(1.d0+v_poisson(I)) ! shear modulus 
            if(dmp_log)write(unit_log,'(5X,A,I5,X,A,X,2(ES15.7))') &
               'E_YOUNG AND V_POISSON FOR M = ', I, '=',&
               E_YOUNG(I), V_POISSON(I) 
         ENDDO
            
         COUNT_E = 0
         DO I=1,MMAX
            DO J=I,MMAX
! Arrange the coefficient of restitution matrix from en_input values
! use coef of rest to determine damping coefficient 
               COUNT_E = COUNT_E + 1
               REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
               REAL_ET(I,J) = DES_ET_INPUT(COUNT_E)            
               MASS_I = (PI/6.d0)*(D_P0(I)**3)*RO_S(I)
               MASS_J = (PI/6.d0)*(D_P0(J)**3)*RO_S(J)
               MASS_EFF = (MASS_I*MASS_J)/(MASS_I+MASS_J)
! In the Hertzian model introduce a factor of 2/7 to the effective mass 
! for tangential direction to get a reduced mass.  See reference: 
! Van der Hoef et al., Advances in Chemical Engineering, 2006, 31, 65-149
!   (see page 94-95)                
               RED_MASS_EFF = (2.d0/7.d0)*MASS_EFF               
               R_EFF = 0.5d0*(D_P0(I)*D_P0(J)/(D_P0(I)+D_P0(J)))
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
                  
               IF(DMP_LOG) &
                  WRITE(UNIT_LOG,'(5X,A,I5,X,I5,X,A,X,2(ES15.7))') &
                  'KN AND KT FOR PAIR ',I, J, '=', &
                  hert_kn(I,J), hert_kt(I,J)
            ENDDO
         ENDDO

! particle-wall contact -------------------->          
         COUNT_E = 0
         DO I = 1, MMAX
            COUNT_E = COUNT_E + 1  
            REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
            REAL_ET_WALL(I) = DES_ET_WALL_INPUT(COUNT_E)
            MASS_I = (PI/6.d0)*(D_P0(I)**3)*RO_S(I)
            MASS_J = MASS_I
            MASS_EFF = MASS_I
            RED_MASS_EFF = (2.d0/7.d0)*MASS_I
            R_EFF = 0.5d0*D_P0(I)
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

         IF(DMP_LOG) WRITE(UNIT_LOG,'(5X,A)') &
            'COLLISION MODEL: Linear Spring-Dashpot (default)'

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
         IF(DMP_LOG) WRITE(UNIT_LOG,'(5X,A,ES17.10,2X,ES15.7)') &
            'KN AND KT = ', KN, KT

! particle-particle contact -------------------->
         COUNT_E = 0
         DO I = 1, MMAX
            DO J = I, MMAX
! Arrange the coefficient of restitution matrix from en_input values
! use coef of rest to determine damping coefficient 
               COUNT_E = COUNT_E + 1
               REAL_EN(I,J) = DES_EN_INPUT(COUNT_E)
               MASS_I = (PI/6.d0)*(D_P0(I)**3.d0)*RO_S(I)
               MASS_J = (PI/6.d0)*(D_P0(J)**3.d0)*RO_S(J)
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
         DO I = 1, MMAX
            COUNT_E = COUNT_E + 1  
            REAL_EN_WALL(I) = DES_EN_WALL_INPUT(COUNT_E)
            MASS_I = (PI*(D_P0(I)**3)*RO_S(I))/6.d0
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

      ENDIF    ! end if/else(des_coll_model)
!-----------------------------------------------------------------<<<


! filling rest of arrays (symmetric matrix)
      DO I = 1, MMAX
         DO J = I, MMAX
            REAL_EN(J, I) = REAL_EN(I,J)
            REAL_ET(J, I) = REAL_ET(J,I)
            DES_ETAN(J,I) = DES_ETAN(I,J)
            DES_ETAT(J,I) = DES_ETAT(I,J)
         ENDDO
      ENDDO
         
! reporting information to logs         
      DO I = 1, MMAX
         DO J = I, MMAX
         IF(DMP_LOG) WRITE(UNIT_LOG,'(5X,A,I10,2X,I10,A,2(ES15.7))') &
            'ETAN AND ETAT FOR PAIR ',&
            I, J, ' = ', DES_ETAN(I,J), DES_ETAT(I,J)
         ENDDO
      ENDDO
      DTSOLID = TCOLL/50.d0

      IF(DMP_LOG.AND..NOT.MPPIC) &
         WRITE(UNIT_LOG,'(5X,A,E17.10,2X,E17.10)') &
         'MIN TCOLL AND DTSOLID = ',TCOLL, DTSOLID
        
         
      IF(MPPIC) THEN 
         DO M = 1, MMAX 
            DES_TAU_P(M) = RO_S(M)*(D_P0(M)**2.d0)/(18.d0*MU_g0)
            IF(DMP_LOG) WRITE(UNIT_LOG,'(/2X,A,I2,A,G17.8)') &
               'TAU_P FOR ', M,'th SOLID PHASE= ', DES_TAU_P(M)
         ENDDO

         DTSOLID = DTSOLID_FACTOR*MINVAL(DES_TAU_P(1:MMAX))
         DTPIC_TAUP = DTSOLID   !maximum dt for point-particles based on taup

         IF(DMP_LOG) THEN 
            WRITE(UNIT_LOG,'(/2X,A,A,A)') 'MPPIC: POINT PARTICLE ',&
               'APPROXIMATION FOR PARTICLE-PARTICLE AND ',&
               'PARTICLE-WALL COLLISIONS'
            WRITE(UNIT_LOG,'(2X,A)') &
               'DTSOLID BASED ON PARTICLE TIME RESPONSE TAUP'
            WRITE(UNIT_LOG,'(2X,A,E17.10)') 'DTSOLID = ', DTSOLID 
         ENDIF
      ENDIF
      
      IF(DMP_LOG.AND.DEBUG_DES) &
         WRITE(UNIT_LOG,'(3X,A)') '<---------- END CFASSIGN ----------'

      RETURN
      END SUBROUTINE CFASSIGN


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! subroutine: compute_volume_of_nodes
! Purpose:  please provide more apt description
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C      

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
      INTEGER :: I, J, K, IJK
      integer :: ipjk, ijpk, ipjpk, ijkp, ipjkp, ijpkp, ipjpkp
      double precision :: vol_ijk, vol_ipjk, vol_ijpk, vol_ipjpk
      double precision :: vol_ijkp, vol_ipjkp, vol_ijpkp, vol_ipjpkp
      double precision :: avg_factor, zcor, vol_node_uncorr
      integer :: FLUID_IND, CUTCELL_IND, WALL_IND
      double precision :: vol_node_count, vol_node_actual_count
      character*100 :: filename
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------            
      
      avg_factor = 0.125d0*(dimn-2) + 0.25d0*(3-dimn)

! compute the volume at the grid nodes 
! grid nodes start from istart2 to iend1
      if(dimn.eq.2) vol_node_count = 4.
      if(dimn.eq.3) vol_node_count = 8.
      do ijk = ijkstart3,ijkend3
         des_vol_node(ijk) = zero
         i = i_of(ijk)
         j = j_of(ijk)
         k = k_of(ijk)
         if(i.lt.istart2 .or. i.gt.iend1) cycle 
         if(j.lt.jstart2 .or. j.gt.jend1) cycle 
         if(k.lt.kstart2 .or. k.gt.kend1) cycle 

         ipjk = funijk(i+1,j,k)
         ijpk = funijk(i,j+1,k)
         ipjpk = funijk(i+1,j+1,k)
         
         vol_ijk   = dx(i)*dy(j)*dz(k)
         vol_ipjk  = dx(i+1)*dy(j)*dz(k)
         vol_ijpk  = dx(i)*dy(j+1)*dz(k)
         vol_ipjpk = dx(i+1)*dy(j+1)*dz(k)
         
         vol_node_uncorr = avg_factor*(vol_ijk + vol_ipjk + vol_ijpk + vol_ipjpk)

         vol_node_actual_count = vol_node_count 

         if(.not.fluid_at(ijk))   then 
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijk  = zero  
         endif

         if(.not.fluid_at(ipjk))    then 
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ipjk  = zero 
         endif

         if(.not.fluid_at(ijpk))    then 
            vol_node_actual_count = vol_node_actual_count - 1
            vol_ijpk  = zero 
         endif

         if(.not.fluid_at(ipjpk))   then 
            vol_node_actual_count = vol_node_actual_count - 1 
            vol_ipjpk = zero 
         endif
         
         des_vol_node(ijk) = avg_factor*(vol_ijk + vol_ipjk + vol_ijpk + vol_ipjpk)

         if (dimn.eq.3) then  
            ijkp = funijk(i,j,k+1)
            ipjkp = funijk(i+1,j,k+1)
            ijpkp = funijk(i,j+1,k+1)
            ipjpkp = funijk(i+1,j+1,k+1)
            
            vol_ijkp   = dx(i)*dy(j)*dz(k+1)
            vol_ipjkp  = dx(i+1)*dy(j)*dz(k+1)
            vol_ijpkp  = dx(i)*dy(j+1)*dz(k+1)
            vol_ipjpkp = dx(i+1)*dy(j+1)*dz(k+1)

            vol_node_uncorr = avg_factor*(vol_node_uncorr + vol_ijkp + vol_ipjkp + vol_ijpkp + vol_ipjpkp)

            if(.not.fluid_at(ijkp))   then 
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ijkp   = zero 
            endif

            if(.not.fluid_at(ipjkp))   then 
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ipjkp  = zero 
            endif

            if(.not.fluid_at(ijpkp))   then 
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ijpkp  = zero 
            endif

            if(.not.fluid_at(ipjpkp))   then 
               vol_node_actual_count = vol_node_actual_count - 1
               vol_ipjpkp = zero 
            endif
                        
            des_vol_node(ijk) = des_vol_node(ijk) + avg_factor*(vol_ijkp + &
            & vol_ipjpkp + vol_ijpkp + vol_ipjkp)
            
         end if 
         if(vol_node_actual_count.gt.zero) then 
            DES_VOL_NODE_RATIO(ijk)  = vol_node_count/vol_node_actual_count
         ELSE
            DES_VOL_NODE_RATIO(ijk) = zero 
         endif
      end do
      
      WRITE(filename,'(A,"_",I5.5,".dat")') 'VOL_NODE_',myPE
      OPEN(1000,file=TRIM(filename),form ='formatted',status='unknown')
      write(1000,'(9(A,3X),A)') 'VARIABLES=',&
         '"I"', '"J"', '"K"', '"FLUID"', '"CUTCELL"', '"WALL"', &
         '"VOL_CELL"', '"VOL_NODE"', '"VOL_NODE_RATIO"'
      write(1000,'(A,2X,3(A,I5,2X))') 'ZONE F=POINT, ',&
         'I=', IEND2-ISTART2+1, ', J=', JEND2-JSTART2+1, &
         ', K=', KEND2-KSTART2+1
      
      DO K=KSTART2, KEND2
         DO J=JSTART2, JEND2
            DO I=ISTART2, IEND2
               IJK  = FUNIJK(I,J,K)
               FLUID_IND = 0 
               CUTCELL_IND = 0 
               WALL_IND = 0
               IF(FLUID_AT(IJK)) FLUID_IND = 1
               IF(CUT_CELL_AT(IJK)) CUTCELL_IND = 1
               if(cartesian_grid)then
                  IF(SCALAR_NODE_ATWALL(IJK)) WALL_IND = 1 
               endif
               if(dimn.eq.2) zcor = zt(k)
               if(dimn.eq.3) zcor = zt(k-1) + dz(k)
               write(1000,'(3(2X,I10),2((2X,I4)),3(2X,G17.8))')&
                  I, J, K, fluid_ind, cutcell_ind, wall_ind, &
                  (vol(ijk)/(DX(I)*DY(J)*DZ(K)))*100., &
                  des_vol_node(ijk), des_vol_node_ratio(ijk) 
            enddo
         enddo
      enddo
      close(1000, status='keep')

      RETURN      
      END SUBROUTINE compute_volume_of_nodes
    
