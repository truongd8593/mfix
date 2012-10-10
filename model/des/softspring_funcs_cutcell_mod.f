module softspring_funcs_cutcell
  
  CONTAINS
    
      SUBROUTINE CFSLIDEWALL2(TANGNT,PARTICLE_SLIDE)
      
      USE discretelement
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: TANGNT
      
      logical PARTICLE_SLIDE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      INTEGER :: K
      DOUBLE PRECISION FTMD, FNMD

!-----------------------------------------------      
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
!----------------------------------------------- 

      PARTICLE_SLIDE = .FALSE.
      FTMD = SQRT(DES_DOTPRDCT(FTAN, FTAN))
      FNMD = SQRT(DES_DOTPRDCT(FNORM,FNORM))

      IF (FTMD.GT.(MEW_W*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FTAN(:) =  MEW_W * FNMD * FTAN(:)/FTMD
         ELSE
            FTAN(:) = -MEW_W * FNMD * TANGNT(:)
         ENDIF
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDEWALL2.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))') &
         'FTMD, mu_w*FNMD = ', FTMD, MEW_W*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDEWALL2.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDEWALL2



      SUBROUTINE CFFCTOWALL2(L, NORM, DIST_LI)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  L
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMN) ::  NORM

      
! distance between centers of particles
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION :: CROSSP(DIMN)

! distance from the contact point to the particle centers 
      DOUBLE PRECISION :: DIST_CL, DIST_CI      

      FC(L,:) = FC(L,:) + FNORM(:) + FTAN(:) 

      
! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI - DES_RADIUS(L)
! remember that the dist_CL for the particle-wall case is actually the distance 
! from the center of L particle to wall and not the substituted particle as was 
! done previously

      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FTAN)
         TOW(L,:) = TOW(L,:) + DIST_CL*CROSSP(:)
      ELSE 
         CROSSP(1) = NORM(1)*FTAN(2) - NORM(2)*FTAN(1)
         TOW(L,1) =  TOW(L,1) + DIST_CL*CROSSP(1)
      ENDIF 

      RETURN
      END SUBROUTINE CFFCTOWALL2

      SUBROUTINE CFRELVEL2(L, II, VRN, VRT, TANGNT, NORM, DIST_LI, &
                          WALLCONTACT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L, II

! marker for particle/wall contact (1=p/w)
      INTEGER WALLCONTACT

      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL, DIST_CI      

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!----------------------------------------------- 
      
      VRELTRANS(:) = (DES_VEL_NEW(L,:) - DES_VEL_NEW(II,:))

! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2       
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
      (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL
      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL + &
         OMEGA_NEW(II,:)*DIST_CI
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DIST_CL + &
         OMEGA_NEW(II,1)*DIST_CI
         OMEGA_SUM(2) = ZERO
      ENDIF

      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! slip velocity of the contact point 
! Equation (8) in Tsuji et al. 1992      
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
! the magnitude of the tangential vector      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction  
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)
      
      
      IF(DEBUG_DES) THEN 
         WRITE(*,*) 'IN CFRELVEL2 ---------------------------------'
         
         WRITE(*,'(3(2x,g17.8))') 'VEL LL = ', DES_VEL_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'VEL I = ', DES_VEL_NEW(II,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA LL = ',OMEGA_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA I = ', OMEGA_NEW(II,:)
         WRITE(*,'(3(2x,g17.8))') 'NORMAL = ', NORM(:)
         WRITE(*,'(3(2x,g17.8))') 'TANGENT = ', TANGNT(:)
         
         WRITE(*,*) 'DIST_CL, DIST_CI = ', DIST_CL, DIST_CI
         
         WRITE(*,'(3(2x,g17.8))') 'VRN, VRT = ', VRN, VRT
         WRITE(*,*) 'OUT OF CFRELVEL2 ---------------------------------'
      ENDIF
      RETURN
      END SUBROUTINE CFRELVEL2
      
      SUBROUTINE CFSLIDE2(TANGNT,PARTICLE_SLIDE)    
      USE discretelement
      USE param1
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(DIMN), INTENT(IN) :: TANGNT
!-----------------------------------------------
! Local Variables
!-----------------------------------------------

      INTEGER :: K
      DOUBLE PRECISION FTMD, FNMD
      LOGICAL PARTICLE_SLIDE

!-----------------------------------------------      
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT
      
!----------------------------------------------- 


      FTMD = SQRT(DES_DOTPRDCT(FTAN, FTAN))
      FNMD = SQRT(DES_DOTPRDCT(FNORM,FNORM))

      IF (FTMD.GT.(MEW*FNMD)) THEN
         PARTICLE_SLIDE = .TRUE.
         IF(DES_DOTPRDCT(TANGNT,TANGNT).EQ.0) THEN
            FTAN(:) =  MEW * FNMD * FTAN(:)/FTMD
         ELSE
            FTAN(:) = -MEW * FNMD * TANGNT(:)
         ENDIF
      ENDIF

      IF(DEBUG_DES .AND. PARTICLE_SLIDE) THEN
         WRITE(*,'(7X,A)') &
            'FROM CFSLIDE2.F ---------->'
         WRITE(*,'(9X,A)') 'PARTICLE_SLIDE = T'
         WRITE(*,'(9X,A,2(ES15.7,X))')&
         'FTMD, mu*FNMD = ', FTMD, MEW*FNMD
         WRITE(*,'(7X,A)') '<----------END CFSLIDE2.F'
      ENDIF

      RETURN
      END SUBROUTINE CFSLIDE2

      SUBROUTINE CFFCTOW2(L, II,  NORM, DIST_LI)
     
      USE param1 
      USE discretelement
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  L, II
      DOUBLE PRECISION, INTENT(IN), DIMENSION(DIMN) ::  NORM(DIMN)

      
! distance between particles
      DOUBLE PRECISION, INTENT(IN) :: DIST_LI

!-----------------------------------------------
! Local variables
!-----------------------------------------------

      DOUBLE PRECISION :: CROSSP(DIMN)

! distance from the contact point to the particle centers 
      DOUBLE PRECISION :: DIST_CL, DIST_CI      

!-----------------------------------------------

      FC(L,:) = FC(L,:) + FNORM(:) + FTAN(:) 
!      FC(II,:) = FC(II,:) - FNORM(:) - FTAN(:)


! calculate the distance from the particle center to the contact point,
! which is taken as the radical line
! dist_ci+dist_cl=dist_li; dist_ci^2+a^2=ri^2;  dist_cl^2+a^2=rl^2       
      DIST_CL = (DIST_LI**2 + DES_RADIUS(L)**2 - DES_RADIUS(II)**2)/&
         (2.d0*DIST_LI)
      DIST_CI = DIST_LI - DIST_CL

      IF(DIMN.EQ.3) THEN 
         CALL DES_CROSSPRDCT(CROSSP, NORM, FTAN)
         TOW(L,:)  = TOW(L,:)  + DIST_CL*CROSSP(:)
!         TOW(II,:) = TOW(II,:) + DIST_CI*CROSSP(:)
! remember torque is R cross FT, which, compared to I particle, are
! both negative for the J particle.  Therefore, the toqrue, unlike tangential
! and normal contact forces, will be in the same direction for both the 
! particles making the pair 
      ELSE 
         CROSSP(1) = NORM(1)*FTAN(2) - NORM(2)*FTAN(1)
         TOW(L,1)  = TOW(L,1)  + DIST_CL*CROSSP(1)
!        TOW(II,1) = TOW(II,1) + DIST_CI*CROSSP(1)
      ENDIF 


      RETURN
      END SUBROUTINE CFFCTOW2

      SUBROUTINE CFRELVEL_WALL2(L, II, VRN, VRT, TANGNT, NORM, DIST_LI, &
                          WALLCONTACT)
      
      USE discretelement
      USE param1
      IMPLICIT NONE     

!-----------------------------------------------
! Local variables
!----------------------------------------------- 
      INTEGER L, II

! marker for particle/wall contact (1=p/w)
      INTEGER WALLCONTACT

      DOUBLE PRECISION TANGNT(DIMN), NORM(DIMN)
      DOUBLE PRECISION TANMOD, VRN, VRT
      DOUBLE PRECISION VRELTRANS(DIMN)
      DOUBLE PRECISION VSLIP(DIMN), &
                       V_ROT(DIMN), OMEGA_SUM(DIMN)

! distance between particles
      DOUBLE PRECISION DIST_LI
! distance from the contact point to the particle centers 
      DOUBLE PRECISION DIST_CL, DIST_CI      

!----------------------------------------------- 
! Functions
!----------------------------------------------- 
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 
!----------------------------------------------- 

! translational relative velocity 
      VRELTRANS(:) = DES_VEL_NEW(L,:)

! rotational contribution  : v_rot
! calculate the distance from the particle center to the wall
      DIST_CL = DIST_LI         !- DES_RADIUS(L)
      IF(DIMN.EQ.3) THEN
         OMEGA_SUM(:) = OMEGA_NEW(L,:)*DIST_CL
      ELSE
         OMEGA_SUM(1) = OMEGA_NEW(L,1)*DIST_CL
         OMEGA_SUM(2) = ZERO
      ENDIF
      
      CALL DES_CROSSPRDCT(V_ROT, OMEGA_SUM, NORM)

! total relative velocity 
      VRELTRANS(:) =  VRELTRANS(:) + V_ROT(:)

! normal component of relative velocity (scalar)
      VRN = DES_DOTPRDCT(VRELTRANS,NORM)

! slip velocity of the contact point 
! Equation (8) in Tsuji et al. 1992      
      VSLIP(:) =  VRELTRANS(:) - VRN*NORM(:)
      
! the magnitude of the tangential vector      
      TANMOD = SQRT(DES_DOTPRDCT(VSLIP,VSLIP))     
      IF(TANMOD.NE.ZERO) THEN
! the unit vector in the tangential direction  
         TANGNT(:) = VSLIP(:)/TANMOD
      ELSE
         TANGNT(:) = ZERO
      ENDIF

! tangential component of relative surface velocity (scalar)
      VRT  = DES_DOTPRDCT(VRELTRANS,TANGNT)
      
      
      IF(DEBUG_DES) THEN 
         WRITE(*,*) 'IN CFRELVEL_WALL2------------------------------'
         
         WRITE(*,'(3(2x,g17.8))') 'VEL LL = ', DES_VEL_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'VEL I = ', DES_VEL_NEW(II,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA LL = ',OMEGA_NEW(L,:)
         WRITE(*,'(3(2x,g17.8))') 'OMEGA I = ', OMEGA_NEW(II,:)
         WRITE(*,'(3(2x,g17.8))') 'NORMAL = ', NORM(:)
         WRITE(*,'(3(2x,g17.8))') 'TANGENT = ', TANGNT(:)
         
         WRITE(*,*) 'DIST_CL, DIST_CI = ', DIST_CL, DIST_CI
         
         WRITE(*,'(3(2x,g17.8))') 'VRN, VRT = ', VRN, VRT
         WRITE(*,*) 'OUT OF CFRELVEL_WALL2 ---------------------------------'
      ENDIF
      RETURN
      END SUBROUTINE CFRELVEL_WALL2

 end module softspring_funcs_cutcell
    
