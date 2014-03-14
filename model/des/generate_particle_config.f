!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 23-Oct-08  C
!  Comments: Added a new routine for clarity of functions              C
!                                                                      C
!  Revision: Modified subroutine for parallel processing               C
!  Authour : Pradeep G                                Date:28-Feb-11   C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      
      SUBROUTINE GENERATE_PARTICLE_CONFIG

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE mfix_pic, only: MPPIC 
      IMPLICIT NONE

        
      IF(MPPIC) THEN 
         CALL GENERATE_PARTICLE_CONFIG_MPPIC
      ELSE
         CALL GENERATE_PARTICLE_CONFIG_DEM 
      ENDIF
      
      END SUBROUTINE GENERATE_PARTICLE_CONFIG

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  SUBROUTINE: GENERATE_PARTICLE_CONFIG                                C
!                                                                      C
!  Purpose: Generate particle configuration based on maximum particle  C
!           radius and filling from top to bottom within specified     C
!           bounds                                                     C
!                                                                      C
!                                                                      C
!  Authors: Rahul Garg                                Date: 08-10-2013 C
!  Now seeding the particles based on the IC regions 
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      
      SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM

      
! Global Variables:
! Runtime Flag: Generate initial particle configuation.
      USE discretelement, only: GENER_PART_CONFIG
! particle radius and density 
      USE discretelement, only: DES_RADIUS, RO_Sol
! particle position new and old  
      USE discretelement, only: des_pos_new, des_pos_old
! particle velocity new and old  
      USE discretelement, only: des_vel_new, des_vel_old
! Simulation dimension (2D/3D)
      USE discretelement, only: DIMN
! X , Y, and Z position of cell faces of computational fluid grid
      USE discretelement, only: XE, YN, ZT
! Number of particles in the system (current)
      USE discretelement, only:  PIP
! Number of DEM solids phases.
      USE discretelement, only: DES_MMAX
! DEM solid phase diameters and densities.
      USE discretelement, only: DES_D_p0, DES_RO_s
! Computed volume of IC region for seeding
      USE discretelement, only: VOL_IC_REGION
! Number of particles seeded, per phase in each IC region
      USE discretelement, only: PART_MPHASE_BYIC
! Number of particles to read from input file.
      USE discretelement, only: PARTICLES

      USE discretelement, only: PEA
! Constant: 3.14159...
      USE constant, only: PI
! Flag indicating that the IC region is defined.
      USE ic, only: IC_DEFINED
! IC Region bulk density (RO_s * EP_s)
      USE ic, only: IC_ROP_s
! IC Region gas volume fraction.
      USE ic, only: IC_EP_G
! min and max physical co-ordinates of IC regions in each direction 
      USE ic, only: IC_X_w, IC_X_e, IC_Y_s, IC_Y_n, IC_Z_b, IC_Z_t
! initally specified velocity field and granular temperature 
      USE ic, only: IC_U_s, IC_V_s, IC_W_s, IC_Theta_M
! Flag to extend the lattice distribution in a given IC to available area
      Use ic, only: IC_DES_FIT_TO_REGION
! Parameter for detecting unspecified values, zero, and one
      USE param1, only: UNDEFINED, UNDEFINED_I, ZERO, ONE, Half
! Parameter for small and large numbers
      USE param1, only: SMALL_NUMBER, LARGE_NUMBER
! Maximum number of initial conditions
      USE param, only: DIMENSION_IC
! first and last index of the physical cells in regular MFIX grid 
      USe compar, only: istart1, iend1, jstart1, jend1, kstart1, kend1
! to access random number generator subroutines 
      USE randomno
      use error_manager

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER :: M, ICV, I,J, K, NP, idim, IC_COUNT
      INTEGER :: lproc_parcount, pcount_byic_byphase(dimension_ic, DES_MMAX) 
      INTEGER :: seed_x, seed_y, seed_z 
      INTEGER :: TOTAL_PARTS_IC, last_counter, TMP_PART_COUNT_INTHIS_IC
      integer, dimension(:), allocatable :: map_to_proc
      double precision lmax_dia,lfac,xp,yp,zp, parts_temp 
      double precision :: XSTART_IC, YSTART_IC, ZSTART_IC, adj_dia, ep_sm 
      double precision :: XEND_IC, YEND_IC, ZEND_IC
      double precision :: xinit, yinit, zinit, ymax 
      !max point where the particle center cud be placed 
      double precision ::  max_ic_pt(3)
      !factor to re-fit the configuration to the IC region size 
      double precision :: refit_fac(3) 
      !particles min and max positions
      double precision :: PART_CEN_MIN(DIMN), PART_CEN_MAX(DIMN)

      !Particle mean velocity and standard deviation
      double precision :: vel_mean(3), vel_sig(3)
       
      double precision, dimension(:,:), allocatable :: pvel_temp
      
!-----------------------------------------------

      CALL INIT_ERR_MSG("GENERATE_PARTICLE_CONFIG_DEM")


! initializing particle count
      lproc_parcount = 0 

      PCOUNT_BYIC_BYPHASE(:,:)  = 0 
! setting particle seed spacing grid to be slightly greater than
! the maximum particle diameter. seed at ~particle radii
      lfac = 1.05d0      

      write(err_msg, 2015)
      CALL FLUSH_ERR_MSG(FOOTER = .false.)

 2015 FORMAT('IC region wise report on particle initialization')
      

      IC_LOOP : DO ICV = 1, DIMENSION_IC 

         IF (IC_DEFINED(ICV)) THEN 

            PART_CEN_MIN(:) = LARGE_NUMBER
            PART_CEN_MAX(:) = SMALL_NUMBER

            TOTAL_PARTS_IC = 0 

            LMAX_DIA = SMALL_NUMBER 
            DO M = 1, DES_MMAX
               
! setting a local value of maximum diameter for each IC region      
               if(PART_MPHASE_BYIC(ICV,M).gt.0) then
                  total_parts_ic = total_parts_ic + PART_MPHASE_BYIC(ICV,M)
                  LMAX_DIA = MAX(LMAX_DIA, DES_D_P0(M))
               endif
            ENDDO

            IF(total_parts_ic.eq.0) cycle IC_LOOP 

            WRITE(ERR_MSG,2016) ICV 

 2016       FORMAT(/1X,70('-')/, 5x, &
            'IC number         = ', I4)

            CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)
            
            !write(*,*) 'total_parts_ic', ICV, total_parts_ic,PART_MPHASE_BYIC(ICV,1:DES_MMAX)
            if(.not.allocated(map_to_proc)) allocate(map_to_proc(total_parts_ic))

            if(Allocated(map_to_proc)) then 
               deallocate(map_to_proc)
               allocate(map_to_proc(total_parts_ic))
            endif
            
            last_counter = 0
            DO M = 1, DES_MMAX 
               if(PART_MPHASE_BYIC(ICV,M).gt.0) then
                  map_to_proc(last_counter+1:last_counter+PART_MPHASE_BYIC(ICV,M)) = M 
                  last_counter = last_counter+PART_MPHASE_BYIC(ICV,M)
               endif
            ENDDO
                             
            XSTART_IC = IC_X_W(ICV) 
            YSTART_IC = IC_Y_S(ICV) 
            ZSTART_IC = IC_Z_B(ICV) 
            XEND_IC   = IC_X_E(ICV)
            YEND_IC   = IC_Y_N(ICV)
            ZEND_IC   = IC_Z_T(ICV)
            
            
            MAX_IC_PT(1) = XEND_IC - 0.5d0*LMAX_DIA*LFAC
            MAX_IC_PT(2) = YEND_IC - 0.5d0*LMAX_DIA*LFAC
            MAX_IC_PT(3) = ZEND_IC - 0.5d0*LMAX_DIA*LFAC
            XINIT = XSTART_IC
            YINIT = YSTART_IC
            ZINIT = ZSTART_IC

            
            DO M = 1, DES_MMAX 
               if(PART_MPHASE_BYIC(ICV,M).eq.0) cycle 

               SEED_X = 1
               SEED_Y = 1
               SEED_Z = 1
               
               ADJ_DIA = LFAC*DES_D_P0(M)               
               
               SEED_X = FLOOR((XEND_IC - XINIT)/ADJ_DIA)
               SEED_Y = FLOOR((YEND_IC - YINIT)/ADJ_DIA)         
               SEED_Z = FLOOR((ZEND_IC - ZINIT)/ADJ_DIA)
               !write(*,*) 'adj_dia = ', adj_dia, lfac, lmax_dia
               !write(*,*) 'seedx  = ', seed_x, seed_y, seed_z
               if(dimn.eq.2) seed_z = 1 

               outerloop :  DO  J = 1, SEED_Y
                  DO  K = 1, SEED_Z 
                     DO  I = 1, SEED_X 
                        XP = XINIT + I*ADJ_DIA - DES_D_P0(M)*0.5D0
                        YP = YINIT + J*ADJ_DIA - DES_D_P0(M)*0.5D0
                        ZP = ZINIT + K*ADJ_DIA - DES_D_P0(M)*0.5D0
                        
                        TMP_PART_COUNT_INTHIS_IC = (PCOUNT_BYIC_BYPHASE(ICV,M)) + 1
                        
                        
                        IF(TMP_PART_COUNT_INTHIS_IC.GT.PART_MPHASE_BYIC(ICV,M)) EXIT outerloop 
                     
!Associate this new particle to the solid phase id based on the map_to_proc defined 
!earlier 

                        
                        PCOUNT_BYIC_BYPHASE(ICV,M)  = PCOUNT_BYIC_BYPHASE(ICV,M) + 1
!Set it to -1 here. It will be set to non '-1' value (equal to particle id)
!on the processor it belongs to

                        MAP_TO_PROC(SUM(PCOUNT_BYIC_BYPHASE(ICV,:))) = -1 
!computational FLUID grid (xe,yn, zt)
                        IF ( XP.GE.XE(ISTART1-1) .AND. XP.LT.XE(IEND1) .AND. &
                             YP.GE.YN(JSTART1-1) .AND. YP.LT.YN(JEND1) .AND. &
                             ZP.GE.ZT(KSTART1-1) .AND. ZP.LT.ZT(KEND1)) THEN
                           LPROC_PARCOUNT = LPROC_PARCOUNT + 1
                           
                           MAP_TO_PROC(SUM(PCOUNT_BYIC_BYPHASE(ICV,:))) = LPROC_PARCOUNT
                           
                           PEA(LPROC_PARCOUNT,1) = .TRUE.
                           DES_RADIUS(LPROC_PARCOUNT) = DES_D_P0(M)*HALF 
                           RO_SOL(LPROC_PARCOUNT) = DES_RO_S(M)
                           DES_POS_NEW(LPROC_PARCOUNT,1) = XP 
                           DES_POS_NEW(LPROC_PARCOUNT,2) = YP 
                           IF(DIMN.EQ.3) DES_POS_NEW(LPROC_PARCOUNT,3) = ZP 
                           
                        ENDIF
                        PART_CEN_MIN(1)  = MIN(XP , PART_CEN_MIN(1))
                        PART_CEN_MIN(2)  = MIN(YP , PART_CEN_MIN(2))
                        
                        PART_CEN_MAX(1)  = MAX(XP , PART_CEN_MAX(1))
                        PART_CEN_MAX(2)  = MAX(YP , PART_CEN_MAX(2))
                        IF(DIMN.EQ.3) THEN 
                           PART_CEN_MIN(3)  = MIN(ZP, PART_CEN_MIN(3))
                           PART_CEN_MAX(3)  = MAX(ZP, PART_CEN_MAX(3))
                        ENDIF

                        YMAX = YP + DES_D_P0(M)*0.5D0
                        
                     end DO
                  end DO
               end DO outerloop
            
               YINIT = YMAX 
            end DO

            IF(IC_DES_FIT_TO_REGION(ICV)) THEN 
               refit_fac = 1.d0 

               DO IDIM = 1, DIMN 
                  IF((PART_CEN_MAX(IDIM)-PART_CEN_MIN(IDIM).GT.LMAX_DIA).AND. &
                       (MAX_IC_PT(IDIM) - PART_CEN_MAX(IDIM).GT.LMAX_DIA)) THEN 
                     
                     REFIT_FAC(IDIM)  = MAX_IC_PT(IDIM)/PART_CEN_MAX(IDIM)
                     !write(*,*) ' REFI, IDIM =', IDIM, REFIT_FAC(IDIM)
                  END IF
               END DO
               
               write(err_msg, 2020) ICV, refit_fac(1:dimn)
               CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)


               DO NP = 1, SUM(PCOUNT_BYIC_BYPHASE(ICV,:))
                  IF (MAP_TO_PROC(NP).NE.-1) THEN 
                     DES_POS_NEW(MAP_TO_PROC(NP), 1:DIMN) = DES_POS_NEW(MAP_TO_PROC(NP), 1:DIMN)*REFIT_FAC(1:DIMN)
                  END IF
               END DO
            END IF !DES_IC_FIT_TO_REGION

            
            last_counter = 0 
            
            DO M = 1, DES_MMAX
               if(PCOUNT_BYIC_BYPHASE(ICV,M).eq.0) cycle 

               VEL_MEAN = ZERO
               VEL_SIG = ZERO 
               VEL_MEAN(1) = IC_U_S(ICV, M)
               VEL_MEAN(2) = IC_V_S(ICV, M)
               IF(DIMN.eq.3) VEL_MEAN(3) = IC_W_S(ICV, M)
               !granular temp is defined as (Variance uprime + Variance vprime + Variance wprime)/3 
               !assuming equal energy in each direction 
               !Variance uprime  = IC_Theta
               !Stdev (or sigma) = sqrt(Variance)
               VEL_SIG(:) = ZERO

               IF(IC_Theta_M(ICV,M).ne.UNDEFINED) & 
               & VEL_SIG(:) = sqrt(IC_Theta_M(ICV, M))
               
               
               write(ERR_MSG,2022) M,  &
               vel_mean(:), IC_theta_m(ICV, M), vel_sig(:)
               CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)

               
               allocate(pvel_temp( PCOUNT_BYIC_BYPHASE(ICV,M) , DIMN))
               do IDIM = 1, DIMN 
                  
                  pvel_temp(:, idim) = vel_mean(idim)

                  IF(vel_sig(idim).gt.zero) then 
                     CALL nor_rno(pvel_temp(1:PCOUNT_BYIC_BYPHASE(ICV,M),IDIM), &
                     vel_mean(idim),vel_sig(idim))
                  ENDIF
               ENDDO
               
               
               DO NP = 1 , PCOUNT_BYIC_BYPHASE(ICV,M)
                  IF (MAP_TO_PROC(NP+last_counter).eq.-1) cycle 

                  DES_VEL_NEW(MAP_TO_PROC(NP+last_counter), 1:DIMN) = pvel_temp(np,1:dimn)
                  DES_VEL_OLD(MAP_TO_PROC(NP+last_counter), 1:DIMN) = pvel_temp(np,1:dimn)
               enddo
               
               last_counter = last_counter + PCOUNT_BYIC_BYPHASE(ICV,M)
               
               deallocate(pvel_temp)

               EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
               WRITE(ERR_MSG,2017) EP_SM, PART_MPHASE_BYIC(ICV, M), & 
               PCOUNT_BYIC_BYPHASE(ICV,M)  

               CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)
               
               IF(PCOUNT_BYIC_BYPHASE(ICV,M).ne.PART_MPHASE_BYIC(ICV, M)) then 
                  WRITE(ERR_MSG, 2018) M
                  CALL FLUSH_ERR_MSG(HEADER = .false., Footer = .false.)
               ENDIF

            ENDDO
         END IF
         
      END DO IC_LOOP
      
! setting pip to particle count
      pip = lproc_parcount 

      !WRITE(*,*) 'PIP on proc ', mype , '  = ', pip 
!      DO ICV = 1, DIMENSION_IC 
!         
!         IF (IC_DEFINED(ICV).AND.IC_EP_G(ICV).lt.ONE) THEN 
            
!            DO M = 1, DES_MMAX
!               EP_SM = IC_ROP_S(ICV,M)/DES_RO_s(M)
!               IF(DMP_LOG)        WRITE(UNIT_LOG,2017) M,  &
!                    EP_SM, PART_MPHASE_BYIC(ICV, M), PCOUNT_BYIC_BYPHASE(ICV,M)  
!               IF(MYPE.EQ.PE_IO)  WRITE(*,2017) M,  &
!                    EP_SM, PART_MPHASE_BYIC(ICV, M), PCOUNT_BYIC_BYPHASE(ICV,M)
               
!               IF(PCOUNT_BYIC_BYPHASE(ICV,M).ne.PART_MPHASE_BYIC(ICV, M)) then 
!                  IF(DMP_LOG)       write(unit_log, 2018) M
!                  if(mype.eq.pe_io) write(*, 2018) M 
!               ENDIF
!            ENDDO
!            
!         ENDIF
         
!      ENDDO
            
      
 2017 FORMAT(5x, &
      'Solids vol fraction for M   =  ', (G15.8,2X), /5x, & 
      '# of particles implied from IC for phase M  = ',  I10, /5x, &
      '# of particles actually seeded for phase M  = ',  I10)
      
 2018 FORMAT(1X,70('.'),/,5x, &
      '####  Warning for phase Index, M  ###', I5,2X, /5x, & 
      'Difference in mass of solids initialized and desired')
      
 2019 FORMAT( 1X,70('*')/)
      
 2020 Format(/5x, 'Refitting to box for IC', I4, /5x,   &
      'Refitting factors (1-DIMN): ', 3(2x,g17.8))
      

 2022 FORMAT(1X,70('.'),/5x, & 
      'PHASE INDEX, M                              =  ', I5,2X, /5x, & 
      'INITIALIZING SOLIDS VELOCITY FIELD', /5x, & 
      'Mean velocity direction wise                =  ', 3(G15.8,2X), /5x, & 
      'Use specified initial granular temperature  =  ', (G15.8,2X), /5x, & 
      'Velocity standard deviation direction wise  =  ', 3(G15.8,2X))
      
      CALL FINL_ERR_MSG
          
      RETURN
    END SUBROUTINE GENERATE_PARTICLE_CONFIG_DEM



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Subroutine: GENERATE_PARTICLE_CONFIG_MMPPIC                         C
!  Purpose: generates particle position distribution for MPPIC         C
!  Author: Rahul Garg                                 Date: 3-May-2011 C      
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      USE param 
      USE param1
      USE geometry
      USE funits
      USE compar      
      USE discretelement
      USE run
      USE constant
      USE physprop
      USE fldvar 
      USE indices 
      USE randomno
      USE mfix_pic 
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------       
      INTEGER :: L, I, J, K, M, IDIM , IJK
      INTEGER :: PART_COUNT, CNP_CELL_COUNT, IPCOUNT
      DOUBLE PRECISION :: DOML(DIMN), CORD_START(DIMN) 
      DOUBLE PRECISION :: REAL_PARTS(DIM_M), STAT_WT
      INTEGER :: LPROC_PARCOUNT
      DOUBLE PRECISION :: VOLIJK
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: RANDPOS
!-----------------------------------------------  
! Include statement functions
!-----------------------------------------------  
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------  
      
      IF(dmp_log.and.debug_des) WRITE(unit_log,'(3X,A)') &
      '---------- START GENERATE_PARTICLE_CONFIG MPPIC ---------->'
      PART_COUNT = 0
      
      DO K = KSTART1, KEND1 
         DO J = JSTART1, JEND1
            DO I = ISTART1, IEND1 
               
               IJK  = FUNIJK(I,J,K)
               CORD_START(1) = XE(I-1)
               CORD_START(2) = YN(J-1)
               IF(DIMN.EQ.3) CORD_START(3) = ZT(K-1)
               DOML(1) = DX(I)
               DOML(2) = DY(J)
               IF(DIMN.EQ.3) DOML(3) = DZ(K)
               DO M = 1, DES_MMAX 
                  CNP_CELL_COUNT = CNP_ARRAY(IJK, M)
                  IF(CNP_CELL_COUNT.EQ.0) CYCLE 
                  ALLOCATE(RANDPOS(CNP_CELL_COUNT, DIMN))
                  RANDPOS = ZERO 
                  DO IDIM = 1, DIMN
                     CALL UNI_RNO(RANDPOS(1:CNP_CELL_COUNT, IDIM))
                  ENDDO

                  VOLIJK = VOL(IJK)
                 
! this is based ep_s that would be set based on IC_EPS from
! two-fluid IC specification. 
! this should be replaced with a true dem quantity                 
                  REAL_PARTS(M) = 6.d0*EP_S(IJK,M)*VOLIJK/&
                     (PI*(DES_D_p0(M)**3.d0))
                  
                  IF(MPPIC_CONSTANTNPC) THEN 
! calculate the statistical weight for CP's belonging to this solids
! phase
                     STAT_WT = REAL_PARTS(M)/REAL(CNP_CELL_COUNT)
                  ELSEIF(MPPIC_CONSTANTWT) THEN
! although the weight was specified in the input file, the integer 
! number of CP's requires recalculating statistical weight. This 
! slightly different statistical weight will ensure that the initial 
! volume fraction is as inputted. If the input statwt_pic is used,
! then the initial volume fraction might be slightly less than the 
! input initial volume fraction

                     STAT_WT = REAL_PARTS(M)/REAL(CNP_CELL_COUNT)
                  ENDIF
                  
                  DO IPCOUNT = 1, CNP_CELL_COUNT
                        DES_POS_OLD(PART_COUNT + IPCOUNT, :) = CORD_START(:) + RANDPOS(IPCOUNT, :)*DOML(:)
                     
                        DES_POS_NEW(PART_COUNT + IPCOUNT, :) =  DES_POS_OLD(PART_COUNT + IPCOUNT, :)
                        
                        DES_VEL_NEW(PART_COUNT + IPCOUNT, :) =  ZERO
                        DES_VEL_OLD(PART_COUNT + IPCOUNT, :) =  ZERO
                        DES_RADIUS(PART_COUNT + IPCOUNT) = DES_D_p0(M)*HALF
                        RO_Sol(PART_COUNT + IPCOUNT) = DES_RO_S(M)
                     
                        DES_STAT_WT(PART_COUNT + IPCOUNT) = STAT_WT

                        MARK_PART(PART_COUNT + IPCOUNT) = 1
                        
                        PIJK(PART_COUNT + IPCOUNT,1) = I
                        PIJK(PART_COUNT + IPCOUNT,2) = J
                        PIJK(PART_COUNT + IPCOUNT,3) = K
                        PIJK(PART_COUNT + IPCOUNT,4) = IJK
                        PIJK(PART_COUNT + IPCOUNT,5) = M
                        
                        IF(DES_POS_NEW(PART_COUNT + IPCOUNT,2).LE.YLENGTH/2.d0) MARK_PART(PART_COUNT + IPCOUNT) = 0
                        !MARK_PART(PART_COUNT + IPCOUNT) = myPE 

                        PEA(PART_COUNT + IPCOUNT,1) = .true.
                  ENDDO
                  PART_COUNT = PART_COUNT + CNP_CELL_COUNT 
                  DEALLOCATE(RANDPOS)
               ENDDO
! set the cnp_array to zero. It will be used for handling inflow later                
               CNP_ARRAY(IJK,:) = 0 


            end DO
         end DO
      end DO

! setting pip to part_count 
      PIP = PART_COUNT 

      IF(DMP_LOG) WRITE(UNIT_LOG,'(2X, A,/,10X, A, i5, /,2X, A, i10,/)') &
      'In GENERATE_PARTICLE_CONFIG MPPIC ',  & 
      'FROM pe =', mype, & 
      'NUMBER OF PARCELS INITIATED ON THIS PROC = ', PIP
      
      WRITE(*,'(2X, A,/,10X, A, i5, /,2X, A, i10,/)') &
      'In GENERATE_PARTICLE_CONFIG MPPIC ',  & 
      'FROM pe =', mype, & 
      'NUMBER OF PARCELS INITIATED ON THIS PROC = ', PIP
      
      IF(PART_COUNT.NE.SUM(CNP_PIC(1:DES_MMAX))) THEN 
         IF(DMP_LOG) THEN 
            WRITE(UNIT_LOG,*) 'ERROR IN GENERATE_PARTICLE_CONFIG_MPPIC'
            WRITE(UNIT_LOG,*) &
          & 'NUMBER OF PARTICLES INITIALIZED (', PART_COUNT, '), NOT &
          & EQUAL TO EARLIER CALCULATED TOTAL NUMBER OF PARTICLES (', & 
          & SUM(CNP_PIC(1:DES_MMAX)), ' )'

            WRITE(UNIT_LOG,*) 'TERMINAL ERROR: STOPPING'
         ENDIF
         CALL mfix_exit(myPE) 
         
      END IF
      
      IF(DMP_LOG.and.debug_des) WRITE(UNIT_LOG,'(3X,A)') &
      '---------- END GENERATE_PARTICLE_CONFIG MPPIC ---------->'
      !CALL mfix_exit(mypE)
      END SUBROUTINE GENERATE_PARTICLE_CONFIG_MPPIC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
! Subroutine : writeic
! Purpose    : write the initial position and velocity of the particles 
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE writeic
!-----------------------------------------------
! Modules
!-----------------------------------------------
      use compar 
      use mpi_utility 
      USE discretelement
      use desmpi
      use cdist 
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER:: NP,lunit,li,lparcnt,lcurpar
      integer :: lproc,llocalcnt,lglocnt,lgathercnts(0:numpes-1)
      character(30) :: lfilename
      double precision,dimension(:,:),allocatable :: ltemparray
!-----------------------------------------------

! based on distributed or single IO open the input file
      WRITE(*,*) 'WRITING IC.dat'
      lunit = 25 
      if (bdist_io) then 
         write(lfilename,'("ic_",I4.4,".dat")')mype
         open(unit=lunit, file=lfilename, form="formatted")
      else 
         if(mype.eq.pe_io) then  
            lfilename= "ic.dat"
            open(unit=lunit, file=lfilename, form="formatted")
         endif 
      endif  

! Write the information      
!----------------------------------------------------------------->>>      

      if (bdist_io) then 
         lparcnt = 1 
         do lcurpar = 1,pip
            if(lparcnt.gt.pip) exit
            if(.not.pea(lcurpar,1)) cycle 
            lparcnt = lparcnt+1
            if(pea(lcurpar,4)) cycle 
            write(lunit,'(10(2x,es12.5))') (des_pos_new(lcurpar,li),li=1,dimn),&
            des_radius(lcurpar),ro_sol(lcurpar), (des_vel_new(lcurpar,li),li=1,dimn)
         enddo

!-----------------------------------------------------------------<<<
      else   ! else branch if not bdist_IO
!----------------------------------------------------------------->>>

! set parameters required for gathering info at PEIO and write as single file   
         lglocnt = 10
         llocalcnt = pip - ighost_cnt 
         call global_sum(llocalcnt,lglocnt) 
         allocate (ltemparray(lglocnt,2*dimn+2)) 
         allocate (dprocbuf(llocalcnt),drootbuf(lglocnt))
         igath_sendcnt = llocalcnt 
         lgathercnts = 0
         lgathercnts(mype) = llocalcnt
         call global_sum(lgathercnts,igathercnts)
         idispls(0) = 0 
         do lproc = 1,numpes-1 
            idispls(lproc) = idispls(lproc-1) + igathercnts(lproc-1)  
         end do 
         do li = 1,dimn 
            call des_gather(des_pos_new(:,li))
            if(mype.eq.pe_io) ltemparray(1:lglocnt,li) = drootbuf(1:lglocnt)
         end do 
         do li = 1,dimn 
            call des_gather(des_vel_new(:,li))
            if(mype.eq.pe_io)ltemparray(1:lglocnt,dimn+li) = drootbuf(1:lglocnt)
         end do 
         call des_gather(des_radius)
         if(mype.eq.pe_io)ltemparray(1:lglocnt,2*dimn+1) = drootbuf(1:lglocnt)
         call des_gather(ro_sol)
         if(mype.eq.pe_io)ltemparray(1:lglocnt,2*dimn+2) = drootbuf(1:lglocnt)
         if (mype.eq.pe_io) then  
            do lcurpar = 1,lglocnt
               write(lunit,'(10(2x,es12.5))') (ltemparray(lcurpar,li),li=1,2*dimn+2)
            enddo
         end if 
         deallocate(ltemparray)
         deallocate(dprocbuf,drootbuf)

      endif   ! end if/else bdist_io
!-----------------------------------------------------------------<<<

      if(bdist_io .or. mype.eq.pe_io) close(lunit)

      RETURN
      END SUBROUTINE writeic



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!      
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE init_particles_jn
      USE randomno
      USE discretelement
      USE constant 
      USE compar 
      USE geometry
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
! number of total nodes
      INTEGER NNODES      
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
!-----------------------------------------------      


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
         ENDDO
      ENDIF
      
      DES_VEL_NEW(:,:) = DES_VEL_OLD(:,:)

! The writing of files below needs to be updated for MPI case 
! updating/writing initial particle configuration files      
      NNODES = NODESI*NODESJ*NODESK
      IF(NNODES.GT.1)   RETURN

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

      RETURN
      END SUBROUTINE init_particles_jn
