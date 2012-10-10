!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: CHECK_DATA_06                                           C
!  Purpose: check the initial conditions input section                 C
!     - check geometry of any specified IC region                      C
!     - check specification of physical quantities                     C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 27-JAN-92  C
!                                                                      C
!  Revision Number: 1                                                  C
!  Purpose: Check the specification of physical quantities             C
!  Author: M. Syamlal                                 Date: 24-JUL-92  C
!  Reviewer: W. Rogers                                Date: 11-DEC-92  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE CHECK_DATA_06 

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1 
      USE geometry
      USE ic
      USE fldvar
      USE physprop
      USE run
      USE indices
      USE funits 
      USE scalars
      USE compar      
      USE mpi_utility      
      USE sendrecv    
      USE rxns
      USE discretelement      
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------

!    ICV,I,J,K  - loop/variable indices
!    I_w,I_e    -  cell indices in X direction from IC_X_e,IC_X_w
!    J_s,J_n    -  cell indices in Y direction from IC_Y_s,IC_Y_n
!    K_b,K_t    -  cell indices in Z direction from IC_Z_b,IC_Z_t
!    IC2 - Last two digits of IC

      INTEGER :: ICV
      INTEGER :: IC2
      INTEGER :: I_w, I_e, J_s, J_n, K_b, K_t
      INTEGER :: I, J, K, IJK
      INTEGER :: M, N, IER
      DOUBLE PRECISION SUM, SUM_EP

!-----------------------------------------------
! External functions
!-----------------------------------------------
      LOGICAL, EXTERNAL :: COMPARE 
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
!-----------------------------------------------


! Read thermochemical database
      If(.not.database_read .and. species_name(1) /= UNDEFINED_C) &
        call read_database(IER)


! Initialize the icbc_flag array.  If not a NEW run then do not
! check the initial conditions.
       DO K = Kstart3, Kend3       
         DO J = Jstart3, Jend3 
           DO I = Istart3, Iend3 
               IJK = FUNIJK(I,J,K)
               IF (RUN_TYPE == 'NEW') THEN 
                  ICBC_FLAG(IJK) = '   ' 
               ELSE 
                  ICBC_FLAG(IJK) = '.--' 
               ENDIF 

! If at domain boundaries then set default values (wall or, if
! specified, cyclic)               
               IF (DO_K) THEN 
                  IF (K==KMIN3 .OR. K==KMIN2 .OR. &
                      K==KMAX2 .OR. K==KMAX3) THEN 
                     IF (CYCLIC_Z_PD) THEN 
                        ICBC_FLAG(IJK) = 'C--' 
                     ELSEIF (CYCLIC_Z) THEN 
                        ICBC_FLAG(IJK) = 'c--' 
                     ELSE 
                        ICBC_FLAG(IJK) = 'W--' 
                     ENDIF 
                  ENDIF 
               ENDIF 

               IF (DO_J) THEN 
                  IF (J==JMIN3 .OR. J==JMIN2 .OR. &
                      J==JMAX2 .OR. J==JMAX3) THEN 
                     IF (CYCLIC_Y_PD) THEN 
                        ICBC_FLAG(IJK) = 'C--' 
                     ELSEIF (CYCLIC_Y) THEN 
                        ICBC_FLAG(IJK) = 'c--' 
                     ELSE 
                        ICBC_FLAG(IJK) = 'W--' 
                     ENDIF 
                  ENDIF 
               ENDIF 

               IF (DO_I) THEN 
                  IF (I==IMIN3 .OR. I==IMIN2 .OR. &
                      I==IMAX2 .OR. I==IMAX3) THEN 
                     IF (CYCLIC_X_PD) THEN 
                        ICBC_FLAG(IJK) = 'C--' 
                     ELSEIF (CYCLIC_X) THEN 
                        ICBC_FLAG(IJK) = 'c--' 
                     ELSE 
                        ICBC_FLAG(IJK) = 'W--' 
                     ENDIF 
                  ENDIF 
                  IF (I==1 .AND. CYLINDRICAL .AND. XMIN==ZERO) &
                     ICBC_FLAG(IJK) = 'S--' 
               ENDIF 

! corner cells are wall cells
               IF ((I==IMIN3 .OR. I==IMIN2 .OR. I==IMAX2 .OR. I==IMAX3) .AND. &
                   (J==JMIN3 .OR. J==JMIN2 .OR. J==JMAX2 .OR. J==JMIN3) .AND. &
                   (K==KMIN3 .OR. K==KMIN2 .OR. K==KMAX2 .OR. K==KMAX3)) THEN 
                  IF (ICBC_FLAG(IJK) /= 'S--') ICBC_FLAG(IJK) = 'W--' 
          
               ENDIF 
       
            ENDDO   ! end do loop (i=istart3, iend3)
         ENDDO   ! end do loop (j=jstart3, jend3)
      ENDDO   ! end do loop (k=kstart3, kend3)

      

! Check geometry of any specified IC region      
      DO ICV = 1, DIMENSION_IC 
         IC_DEFINED(ICV) = .FALSE. 
         IF (IC_X_W(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_X_E(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Y_S(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Y_N(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Z_B(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_Z_T(ICV) /= UNDEFINED) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_I_W(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_I_E(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_J_S(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_J_N(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_K_B(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 
         IF (IC_K_T(ICV) /= UNDEFINED_I) IC_DEFINED(ICV) = .TRUE. 

         IF (IC_DEFINED(ICV)) THEN 

            IF (IC_X_W(ICV)==UNDEFINED .AND. IC_I_W(ICV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  IC_X_W(ICV) = ZERO 
               ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) &
                     'IC_X_w and IC_I_w ', ICV 
                     call mfix_exit(myPE) 
               ENDIF 
            ENDIF 
            IF (IC_X_E(ICV)==UNDEFINED .AND. IC_I_E(ICV)==UNDEFINED_I) THEN 
               IF (NO_I) THEN 
                  IC_X_E(ICV) = XLENGTH 
               ELSE 
                 IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_X_e and IC_I_e ', ICV 
                 call mfix_exit(myPE) 
               ENDIF 
            ENDIF 
            IF (IC_Y_S(ICV)==UNDEFINED .AND. IC_J_S(ICV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  IC_Y_S(ICV) = ZERO 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_Y_s and IC_J_s ', ICV 
                  call mfix_exit(myPE) 
               ENDIF 
            ENDIF 
            IF (IC_Y_N(ICV)==UNDEFINED .AND. IC_J_N(ICV)==UNDEFINED_I) THEN 
               IF (NO_J) THEN 
                  IC_Y_N(ICV) = YLENGTH 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_Y_n and IC_J_n ', ICV 
                  call mfix_exit(myPE) 
               ENDIF 
            ENDIF 
            IF (IC_Z_B(ICV)==UNDEFINED .AND. IC_K_B(ICV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  IC_Z_B(ICV) = ZERO 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_Z_b and IC_K_b ', ICV 
                  call mfix_exit(myPE) 
               ENDIF 
            ENDIF 
            IF (IC_Z_T(ICV)==UNDEFINED .AND. IC_K_T(ICV)==UNDEFINED_I) THEN 
               IF (NO_K) THEN 
                  IC_Z_T(ICV) = ZLENGTH 
               ELSE 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_Z_t and IC_K_t ', ICV 
                  call mfix_exit(myPE) 
               ENDIF 
            ENDIF 

         ENDIF   ! end if (ic_defined(icv))
      ENDDO   ! end loop over (icv = 1,dimension_ic)


      DO ICV = 1, DIMENSION_IC 
         IF (IC_X_W(ICV)/=UNDEFINED .AND. IC_X_E(ICV)/=UNDEFINED) THEN 
            IF (NO_I) THEN 
               I_W = 1 
               I_E = 1 
            ELSE 
               CALL CALC_CELL (XMIN, IC_X_W(ICV), DX, IMAX, I_W) 
               I_W = I_W + 1 
               CALL CALC_CELL (XMIN, IC_X_E(ICV), DX, IMAX, I_E) 
            ENDIF 
            IF (IC_I_W(ICV)/=UNDEFINED_I .OR. IC_I_E(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_I_W(ICV), I_W, ICV, 'IC - west') 
               CALL LOCATION_CHECK (IC_I_E(ICV), I_E, ICV, 'IC - east') 
            ELSE 
               IC_I_W(ICV) = I_W 
               IC_I_E(ICV) = I_E 
            ENDIF 
         ENDIF 

         IF (IC_Y_S(ICV)/=UNDEFINED .AND. IC_Y_N(ICV)/=UNDEFINED) THEN 
            IF (NO_J) THEN 
               J_S = 1 
               J_N = 1 
            ELSE 
               CALL CALC_CELL (ZERO, IC_Y_S(ICV), DY, JMAX, J_S) 
               J_S = J_S + 1 
               CALL CALC_CELL (ZERO, IC_Y_N(ICV), DY, JMAX, J_N) 
            ENDIF 
            IF (IC_J_S(ICV)/=UNDEFINED_I .OR. IC_J_N(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_J_S(ICV), J_S, ICV, 'IC - south') 
               CALL LOCATION_CHECK (IC_J_N(ICV), J_N, ICV, 'IC - north') 
            ELSE 
               IC_J_S(ICV) = J_S 
               IC_J_N(ICV) = J_N 
            ENDIF 
         ENDIF 

         IF (IC_Z_B(ICV)/=UNDEFINED .AND. IC_Z_T(ICV)/=UNDEFINED) THEN 
            IF (NO_K) THEN 
               K_B = 1 
               K_T = 1 
            ELSE 
               CALL CALC_CELL (ZERO, IC_Z_B(ICV), DZ, KMAX, K_B) 
               K_B = K_B + 1 
               CALL CALC_CELL (ZERO, IC_Z_T(ICV), DZ, KMAX, K_T) 
            ENDIF 
            IF (IC_K_B(ICV)/=UNDEFINED_I .OR. IC_K_T(ICV)/=UNDEFINED_I) THEN 
               CALL LOCATION_CHECK (IC_K_B(ICV), K_B, ICV, 'IC - bottom') 
               CALL LOCATION_CHECK (IC_K_T(ICV), K_T, ICV, 'IC - top') 
            ELSE 
               IC_K_B(ICV) = K_B 
               IC_K_T(ICV) = K_T 
            ENDIF 
         ENDIF 

      ENDDO   ! end loop over (icv=1,dimension_ic)

      

      DO ICV = 1, DIMENSION_IC 

         IF (IC_DEFINED(ICV)) THEN 
!----------------------------------------------------------------->>>
! For restart runs IC is defined only if IC_TYPE='PATCH'
            IF (RUN_TYPE/='NEW' .AND. IC_TYPE(ICV)/='PATCH') THEN 
               IC_DEFINED(ICV) = .FALSE. 
               CYCLE  
            ENDIF 

            IF (IC_I_W(ICV) > IC_I_E(ICV)) GO TO 900 
            IF (IC_J_S(ICV) > IC_J_N(ICV)) GO TO 900 
            IF (IC_K_B(ICV) > IC_K_T(ICV)) GO TO 900 
            IF (IC_I_W(ICV)<IMIN1 .OR. IC_I_W(ICV)>IMAX1) GO TO 900 
            IF (IC_I_E(ICV)<IMIN1 .OR. IC_I_E(ICV)>IMAX1) GO TO 900 
            IF (IC_J_S(ICV)<JMIN1 .OR. IC_J_S(ICV)>JMAX1) GO TO 900 
            IF (IC_J_N(ICV)<JMIN1 .OR. IC_J_N(ICV)>JMAX1) GO TO 900 
            IF (IC_K_B(ICV)<KMIN1 .OR. IC_K_B(ICV)>KMAX1) GO TO 900 
            IF (IC_K_T(ICV)<KMIN1 .OR. IC_K_T(ICV)>KMAX1) GO TO 900 


! If IC_TYPE is 'PATCH' then do not need to check whether all the
! variables are specified (otherwise check)
            IF (IC_TYPE(ICV) /= 'PATCH') THEN 

! Check the specification of physical quantities 
! GAS PHASE Quantities
! -------------------------------------------->>>
! check that gas phase velocity components are initialized
               IF (IC_U_G(ICV) == UNDEFINED) THEN 
                  IF (NO_I) THEN 
                     IC_U_G(ICV) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_U_g', ICV 
                     call mfix_exit(myPE) 
                  ENDIF 
               ENDIF 
               IF (IC_V_G(ICV) == UNDEFINED) THEN 
                  IF (NO_J) THEN 
                     IC_V_G(ICV) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_V_g', ICV 
                     call mfix_exit(myPE) 
                  ENDIF 
               ENDIF 
               IF (IC_W_G(ICV) == UNDEFINED) THEN 
                  IF (NO_K) THEN 
                     IC_W_G(ICV) = ZERO 
                  ELSE 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_W_g', ICV 
                     call mfix_exit(myPE) 
                  ENDIF 
               ENDIF 

! check that gas phase void fraction is initialized
               IF (IC_EP_G(ICV) == UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_EP_g', ICV 
                  call mfix_exit(myPE) 
               ENDIF 

! check that if the gas phase pressure is initialized and the gas is
! compressible that the gas phase pressure is not zero or negative
               IF (IC_P_G(ICV) /= UNDEFINED) THEN 
                  IF (RO_G0==UNDEFINED .AND. IC_P_G(ICV)<=ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1010) ICV, IC_P_G(ICV) 
                     call mfix_exit(myPE) 
                  ENDIF 
               ENDIF 

               IF ((ENERGY_EQ .OR. RO_G0==UNDEFINED .OR. MU_G0==UNDEFINED)&
                   .AND. IC_T_G(ICV)==UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_T_g', ICV 
                     call mfix_exit(myPE) 
               ENDIF 

               IF (ENERGY_EQ) THEN 
                  IF (IC_GAMA_RG(ICV) < ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1001) 'IC_GAMA_Rg', ICV 
                     call mfix_exit(myPE) 
                  ELSEIF (IC_GAMA_RG(ICV) > ZERO) THEN 
                     IF (IC_T_RG(ICV) == UNDEFINED) THEN 
                       IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_T_Rg', ICV 
                       call mfix_exit(myPE) 
                     ENDIF 
                  ENDIF 
               ENDIF 

! check sum of the gas species mass fractions
! sum together those gas phase species mass fractions that are defined
               SUM = ZERO
               DO N = 1, NMAX(0) 
                  IF (IC_X_G(ICV,N) /= UNDEFINED) &
                     SUM = SUM + IC_X_G(ICV,N) 
               ENDDO 
! if the gas phase has more than one species, then set any undefined
! species mass fractions to zero. Warn the user if the sum of species 
! mass fractions is not one and indicate which species is undefined
! (mfix will exit in next check).   
               DO N = 1, NMAX(0) 
                  IF (IC_X_G(ICV,N) == UNDEFINED) THEN 
                     IF (.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)&
                        WRITE (UNIT_LOG, 1050) ICV, N 
                     IC_X_G(ICV,N) = ZERO 
                  ENDIF 
               ENDDO 
! if sum of gas phase species mass fraction not 1....
               IF (.NOT.COMPARE(ONE,SUM)) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1055) ICV 
                  IF (SPECIES_EQ(0) .OR. RO_G0==UNDEFINED .AND. &
                      MW_AVG==UNDEFINED) THEN
                     call mfix_exit(myPE)  
                  ENDIF
               ENDIF 

               DO N = 1, NScalar
                  IF (IC_Scalar(ICV,N) == UNDEFINED) THEN 
                     IC_Scalar(ICV,N) = ZERO 
                  ENDIF 
               ENDDO  
               
               IF(K_Epsilon) THEN
                  IF (IC_K_Turb_G(ICV) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_K_Turb_G', ICV 
                     call mfix_exit(myPE) 
                  ENDIF
                  IF (IC_E_Turb_G(ICV) == UNDEFINED) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1000) 'IC_E_Turb_G', ICV 
                     call mfix_exit(myPE) 
                  ENDIF
               ENDIF
! GAS PHASE Quantities               
! --------------------------------------------<<<


! SOLIDS PHASE Quantities               
! -------------------------------------------->>>
               IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID.OR.MPPIC) THEN

! at this point ic_ep_g must be defined
                  SUM_EP = IC_EP_G(ICV) 
                  DO M = 1, SMAX 
                     IF (IC_ROP_S(ICV,M) == UNDEFINED) THEN 
                        IF (IC_EP_G(ICV) == ONE) THEN 
                           IC_ROP_S(ICV,M) = ZERO 
                        ELSEIF (SMAX == 1 .AND. &
                                .NOT.DES_CONTINUUM_HYBRID) THEN 
                           IC_ROP_S(ICV,M) = (ONE - IC_EP_G(ICV))*RO_S(M)
                        ELSE
! bulk density must be explicitly defined for hybrid model and cannot be
! defined from 1-ic_ep_g
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'IC_ROP_s', ICV, M 
                           call mfix_exit(myPE) 
                        ENDIF
                     ENDIF
! at this point ic_rop_s must be defined
! sum of void fraction and solids volume fractions                  
                     SUM_EP = SUM_EP + IC_ROP_S(ICV,M)/RO_S(M) 

! check sum of the solids phase m species mass fractions                  
                     SUM = ZERO 
                     DO N = 1, NMAX(M)
! sum together those solids phase M species mass fractions that are
! defined
                        IF(IC_X_S(ICV,M,N)/=UNDEFINED) &
                           SUM=SUM+IC_X_S(ICV,M,N) 
                     ENDDO
! if no solids phase M present and no species are defined for phase M 
! then set mass fraction of solids phase M species 1 to one                
                     IF (IC_ROP_S(ICV,M)==ZERO .AND. SUM==ZERO) THEN 
                        IC_X_S(ICV,M,1) = ONE 
                        SUM = ONE 
                     ENDIF 
! set any undefined species mass fractions to zero. Warn the user if 
! the sum of species mass fractions is not one and indicate which 
! species is undefined (mfix will exit in next check). 
                     DO N = 1, NMAX(M) 
                        IF (IC_X_S(ICV,M,N) == UNDEFINED) THEN 
                           IF(.NOT.COMPARE(ONE,SUM) .AND. DMP_LOG)&
                              WRITE (UNIT_LOG, 1110)ICV,M,N 
                           IC_X_S(ICV,M,N) = ZERO 
                        ENDIF 
                     ENDDO
! if sum of solids phase M species mass fraction not 1...
                     IF (.NOT.COMPARE(ONE,SUM)) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1120) ICV, M 
                        IF (SPECIES_EQ(M)) call mfix_exit(myPE)  
                     ENDIF 

! check that solids phase m velocity components are initialized
                     IF (IC_U_S(ICV,M) == UNDEFINED) THEN 
                        IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_I) THEN 
                           IC_U_S(ICV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'IC_U_s', ICV, M 
                           call mfix_exit(myPE) 
                        ENDIF 
                     ENDIF 
                     IF (IC_V_S(ICV,M) == UNDEFINED) THEN 
                        IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_J) THEN 
                           IC_V_S(ICV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'IC_V_s', ICV, M 
                           call mfix_exit(myPE) 
                        ENDIF 
                     ENDIF 
                     IF (IC_W_S(ICV,M) == UNDEFINED) THEN 
                        IF (IC_ROP_S(ICV,M)==ZERO .OR. NO_K) THEN 
                           IC_W_S(ICV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'IC_W_s', ICV, M 
                           call mfix_exit(myPE) 
                        ENDIF 
                     ENDIF 

                     IF (ENERGY_EQ .AND. IC_T_S(ICV,M)==UNDEFINED) THEN 
                        IF (IC_ROP_S(ICV,M) == ZERO) THEN 
                           IC_T_S(ICV,M) = IC_T_G(ICV) 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'IC_T_s', ICV, M 
                           call mfix_exit(myPE) 
                        ENDIF 
                     ENDIF 

                     IF (GRANULAR_ENERGY .AND. &
                         IC_THETA_M(ICV,M)==UNDEFINED) THEN 
                        IF (IC_ROP_S(ICV,M) == ZERO) THEN 
                           IC_THETA_M(ICV,M) = ZERO 
                        ELSE 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                              'IC_Theta_m', ICV, M 
                           call mfix_exit(myPE) 
                        ENDIF 
                     ENDIF 
   
                     IF (ENERGY_EQ) THEN 
                        IF (IC_GAMA_RS(ICV,M) < ZERO) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1101) &
                              'IC_GAMA_Rs', ICV, M 
                           call mfix_exit(myPE) 
                        ELSEIF (IC_GAMA_RS(ICV,M) > ZERO) THEN 
                           IF (IC_T_RS(ICV,M) == UNDEFINED) THEN 
                              IF(DMP_LOG)WRITE (UNIT_LOG, 1100) &
                                 'IC_T_Rs', ICV, M 
                              call mfix_exit(myPE) 
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDDO   ! end loop over (m=1,smax)

! check sum of gas void fraction and all solids volume fractions                  
                  IF (.NOT.DES_CONTINUUM_HYBRID) THEN
                     IF (.NOT.COMPARE(ONE,SUM_EP)) THEN 
                           IF(DMP_LOG)WRITE (UNIT_LOG, 1125) ICV 
                           call mfix_exit(myPE)
                     ENDIF 
                  ELSE
! sum_ep not necessarily one at this point since discrete particles
! present in hybrid model
                     IF (SUM_EP>ONE .OR. SUM_EP<ZERO) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1130) ICV 
                        call mfix_exit(myPE)
                     ENDIF
                  ENDIF
               ELSE   ! discrete_element and .not.des_continuum_hybrid
                  SUM_EP = IC_EP_G(ICV)                        
                  IF (SUM_EP>ONE .OR. SUM_EP<ZERO) THEN 
                     IF(DMP_LOG)WRITE (UNIT_LOG, 1130) ICV 
                     call mfix_exit(myPE)
                  ENDIF
               ENDIF   !  end if/else (.not.discrete_element)
! SOLIDS PHASE Quantities
! --------------------------------------------<<<

            ENDIF      ! end if (ic_type(icv) /= 'patch')


!  Set ICBC flag
            DO K = IC_K_B(ICV), IC_K_T(ICV) 
               DO J = IC_J_S(ICV), IC_J_N(ICV) 
                  DO I = IC_I_W(ICV), IC_I_E(ICV) 
                     IF (.NOT.IS_ON_myPE_plus2layers(I,J,K)) CYCLE
                     IJK = FUNIJK(I,J,K) 
                     ICBC_FLAG(IJK)(1:1) = '.' 
                     IC2 = MOD(ICV,100) 
                     WRITE (ICBC_FLAG(IJK)(2:3), 1150) IC2 
                  ENDDO 
               ENDDO 
            ENDDO 

!-----------------------------------------------------------------<<<
         ELSE      ! else branch if(ic_defined(icv))
!----------------------------------------------------------------->>>

! Check whether physical quantities are specified for undefined
! initial conditions and if so flag error

! GAS PHASE quantities
! -------------------------------------------->>>                
            IF (IC_U_G(ICV) /= UNDEFINED) THEN 
                IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_U_g', ICV 
                call mfix_exit(myPE) 
            ENDIF 
            IF (IC_V_G(ICV) /= UNDEFINED) THEN 
                IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_V_g', ICV 
                call mfix_exit(myPE) 
            ENDIF 
            IF (IC_W_G(ICV) /= UNDEFINED) THEN 
                IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_W_g', ICV 
                call mfix_exit(myPE) 
            ENDIF 
            IF (IC_EP_G(ICV) /= UNDEFINED) THEN 
                IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_EP_g', ICV 
                call mfix_exit(myPE) 
            ENDIF 
            IF (IC_T_G(ICV) /= UNDEFINED) THEN 
                IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_T_g', ICV 
                call mfix_exit(myPE) 
            ENDIF 
            IF (IC_T_RG(ICV) /= UNDEFINED) THEN 
                IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_T_Rg', ICV 
                call mfix_exit(myPE) 
            ENDIF 
            DO N = 1, DIMENSION_N_G 
               IF (IC_X_G(ICV,N) /= UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_X_g', ICV 
                  call mfix_exit(myPE) 
               ENDIF 
            ENDDO 
    
            DO N = 1, NScalar 
               IF (IC_Scalar(ICV,N) /= UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_Scalar', ICV 
                  CALL MFIX_EXIT(myPE)
               ENDIF 
            ENDDO 
    
            IF( K_Epsilon ) THEN
               IF (IC_K_Turb_G(ICV) /= UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_K_Turb_G', ICV 
                  CALL MFIX_EXIT(myPE)
               ENDIF 
    
               IF (IC_E_Turb_G(ICV) /= UNDEFINED) THEN 
                  IF(DMP_LOG)WRITE (UNIT_LOG, 1200) 'IC_E_Turb_G', ICV 
                  CALL MFIX_EXIT(myPE)
               ENDIF 
            ENDIF
! GAS PHASE quantities
! --------------------------------------------<<<

! SOLIDS PHASE quantities
! -------------------------------------------->>>
            IF (.NOT.DISCRETE_ELEMENT .OR. DES_CONTINUUM_HYBRID .OR. MPPIC) THEN

               DO M = 1, DIMENSION_M 
                  IF (IC_ROP_S(ICV,M) /= UNDEFINED) THEN 
                      IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_ROP_s', ICV, M 
                      call mfix_exit(myPE) 
                  ENDIF 
                  DO N = 1, DIMENSION_N_S 
                     IF (IC_X_S(ICV,M,N) /= UNDEFINED) THEN 
                         IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_X_s', ICV, M
                         call mfix_exit(myPE) 
                     ENDIF 
                  ENDDO 
                  IF (IC_U_S(ICV,M) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_U_s', ICV, M 
                        call mfix_exit(myPE) 
                  ENDIF 
                  IF (IC_V_S(ICV,M) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_V_s', ICV, M 
                        call mfix_exit(myPE) 
                  ENDIF 
                  IF (IC_W_S(ICV,M) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_W_s', ICV, M 
                        call mfix_exit(myPE) 
                  ENDIF 
                  IF (IC_T_S(ICV,M) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_T_s', ICV, M 
                        call mfix_exit(myPE) 
                  ENDIF 
                  IF (IC_T_RS(ICV,M) /= UNDEFINED) THEN 
                        IF(DMP_LOG)WRITE (UNIT_LOG, 1300) 'IC_T_Rs', ICV, M
                        call mfix_exit(myPE) 
                  ENDIF 
               ENDDO  ! end loop over (m=1,dimension_m)
            ENDIF   ! end if (.not.discrete_element.or.des_continuum_hybrid)
! SOLIDS PHASE Quantities
! --------------------------------------------<<<

         ENDIF     ! end if/else (ic_defined(icv) 
!-----------------------------------------------------------------<<<
      ENDDO        ! end loop over (icv=1,dimension_IC)

      call send_recv(icbc_flag,2)


      RETURN  

! if error in indices go here:
  900 CONTINUE 

      CALL ERROR_ROUTINE ('check_data_06', 'Invalid IC region specified', 0, 2) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC number = ', ICV
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC_I_w(ICV) = ', IC_I_W(ICV) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC_I_e(ICV) = ', IC_I_E(ICV) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC_J_s(ICV) = ', IC_J_S(ICV) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC_J_n(ICV) = ', IC_J_N(ICV) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC_K_b(ICV) = ', IC_K_B(ICV) 
         IF(DMP_LOG)WRITE (UNIT_LOG, *) ' IC_K_t(ICV) = ', IC_K_T(ICV) 
      CALL ERROR_ROUTINE (' ', ' ', 1, 3) 

 1000 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,&
         ') not specified',/1X,70('*')/) 
 1001 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,&
         ') value is unphysical',/1X,70('*')/) 
 1010 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC_P_g( ',I2,&
         ') = ',G12.5,/,&
         ' Pressure should be greater than 0.0 for compressible flow',/1X,70(&
         '*')/) 
 1050 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC_X_g(',I2,',',I2&
         ,') not specified',/1X,70('*')/) 
 1055 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - Sum of gas mass fractions is NOT equal to one',/1X,70('*')/) 
 1100 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,',',I1,&
         ') not specified',/1X,70('*')/) 
 1101 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,',',I1,&
         ') unphysical',/1X,70('*')/) 
 1110 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC_X_s(',I2,',',I2&
         ,',',I2,') not specified',/1X,70('*')/) 
 1120 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - Sum of solids-',I1,' mass fractions is NOT equal to one',/1X,70(&
         '*')/) 
 1125 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - Sum of volume fractions is NOT equal to one',/1X,70('*')/) 

 1130 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: IC number:',I2,&
         ' - void fraction is unphysical (>1 or <0)',/1X,70('*')/) 

 1150 FORMAT(I2.2) 
 1200 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,&
         ') specified',' for an undefined IC region',/1X,70('*')/) 
 1300 FORMAT(/1X,70('*')//' From: CHECK_DATA_06',/' Message: ',A,'(',I2,',',I1,&
         ') specified',' for an undefined IC region',/1X,70('*')/)

      END SUBROUTINE CHECK_DATA_06 


