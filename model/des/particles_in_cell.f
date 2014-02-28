!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Subroutine: PARTICLES_IN_CELL                                       !
!                                                                      !
!  Purpose:                                                            !
!     - For each particle find the computational fluid cell            !
!       containing the particle center.                                !
!     - Calculate the bulk density in each computational fluid         !
!       cell.                                                          !
!     - Calculate the volume average solids velocity in each           !
!       computational fluid cell.                                      !
!     - For parallel processing indices are altered                    !
!                                                                      !
!                                                                      !
!                                                                      !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE PARTICLES_IN_CELL

!-----------------------------------------------
! Modules
!-----------------------------------------------      
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid 
      use desmpi
      USE cutcell 
      USE mfix_pic
      USE des_rxns
      USE run
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase no.    
      INTEGER M
! ijk indices      
      INTEGER I, J, K, IJK, IPROC
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos
! particle count in ijk fluid cell 
      INTEGER, DIMENSION(DIMENSION_3):: particle_count
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

! count for number of particles that were found in the ghost cell and
! hence removed       
      INTEGER :: PIP_DEL_COUNT  
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1)
! number of i,j,k nodes (parallelization)
      INTEGER :: NNODES
! IER for error reporting
      INTEGER IER 
      INTEGER epg_min_loc(1)
      double precision :: epg_min2
! Difference between a particles diameter (density) and the diameter
! (density) of a phase specified in the data file.
      DOUBLE PRECISION dDp, dRho
      
!!$      double precision omp_start, omp_end
!!$      double precision omp_get_wtime	      
!-----------------------------------------------
! Include statement functions      
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------


! following quantities are reset every call to particles_in_cell
      PIP_DEL_COUNT = 0 
      PINC(:) = 0
      IF(FIRST_PASS) THEN
         IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(3X,A)') &
            '---------- START FIRST PASS PARTICLES_IN_CELL ---------->'
      ENDIF

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
      CALL DES_PAR_EXCHANGE


! Assigning PIJK(L,1), PIJK(L,2) and PIJK(L,3) the i, j, k indices 
! of particle L (locating particle on fluid grid). Also determine
! composite ijk index. If first_pass, also assigning PIJK(L,5) the
! solids phase index of particle.
! ---------------------------------------------------------------->>>
!Handan Liu commented here on June 6 2013 and revised as below:
! 1)Directly adding the directives in this loop will result in the datarace;
! 2)Using Reduction-clause for PINC to avoid race condition for OpenMP
!   at 'PINC(IJK)=PINC(IJK)+1';
! 3)Using Reduction-clause for PIP_DEL_COUNT for OpenMP
!   at 'PIP_DEL_COUNT = PIP_DEL_COUNT + 1' for MPPIC. 
! ----------------------------------------------------------------<<<<
!!$omp parallel do default(shared)                               &
!!$omp private(l,m,xpos,ypos,zpos,i,j,k,ijk) schedule (guided,50)  
!!$      omp_start=omp_get_wtime()
!$omp parallel default(shared)            &
!$omp private(l,m,xpos,ypos,zpos,i,j,k,ijk,dDp,dRho)
!$omp do reduction(+:PINC) reduction(+:PIP_DEL_COUNT) schedule (guided,50)  

      DO L = 1, MAX_PIP
! skipping particles that do not exist
         IF(.NOT.PEA(L,1)) CYCLE
! skipping ghost particles
         IF(PEA(L,4)) CYCLE
        
! Assigning local aliases for particle position
         XPOS = DES_POS_NEW(L,1)
         YPOS = DES_POS_NEW(L,2)
         IF (DIMN .EQ. 3) THEN
            ZPOS = DES_POS_NEW(L,3)
         ENDIF


         IF(FIRST_PASS) THEN 
! Determining the solids phase of each particle by matching the diameter
! and density to those specified in the data file. Since reacting flows
! may change a particle's density (variable_density) or diameter 
! (shrinking_particle), each particle's phase is stored in the restart
! when solving the species equations.
            IF((RUN_TYPE == 'NEW') .OR. .NOT.ANY_DES_SPECIES_EQ) THEN
               M_LP: DO M = 1, DES_MMAX
                  dDp  = ABS(2.0d0*DES_RADIUS(L)-DES_D_P0(M))
                  dRho = ABS( RO_Sol(L)-DES_RO_S(M))
                  IF( dDp < SMALL_NUMBER .AND. dRho < SMALL_NUMBER) THEN
                     PIJK(L,5) = M
                     EXIT M_LP
                  ENDIF
               ENDDO M_LP

     
               IF(PIJK(L,5).EQ.0) THEN
                  IF(DMP_LOG) THEN
                     write(unit_log,'(5X,A,A,I10)') &
                        'Problem determining the solids ',&
                        'association for particle: ',L
                     write(unit_log,'(7X,A,(ES15.9))') &
                        'Particle position = ', DES_POS_NEW(L,:)
                     write(unit_log,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                        'Particle diameter = ', 2.0*DES_RADIUS(L),&
                        'and D_P0(1:MMAX)= ', DES_D_P0(1:DES_MMAX)
                     write(unit_log,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                        'Particle density = ', Ro_Sol(L), &
                        'and RO_s(1:MMAX) = ', DES_RO_S(1:DES_MMAX)
                  ENDIF
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF   ! end if(first_pass)


! Using a brute force technique to determine the particle locations in
! the Eulerian fluid grid                 
! In case of parallel processing the cell index is from istart1 to iend1
! XE will have values from istart2-1 to iend2             
! Also particles may be found in ghost cells in the event of particles
! particles entering/exiting the domain (this might occur for restart) 
            DO I = ISTART2,IEND2
               IF( XPOS >= XE(I-1) .and. XPOS < XE(I)) THEN
                  PIJK(L,1) = I
                  EXIT
               ENDIF
            ENDDO
            DO J = jstart2,jend2
               IF(YPOS >= YN(J-1) .and. YPOS < YN(J)) THEN
                  PIJK(L,2) = J
                  EXIT
               ENDIF
            ENDDO
            IF(DIMN.EQ.2) THEN
               K=1
               PIJK(L,3) = 1
            ELSE
               DO K = KSTART2, KEND2
                  IF(ZPOS >= ZT(K-1) .and. ZPOS < ZT(K)) THEN 
                     PIJK(L,3) = K
                     EXIT
                 ENDIF
               ENDDO
            ENDIF

         ELSE    ! if not first_pass

! Using an incremental approach to determine the new location of the 
! particles.
! Note a particle cannot move by more than one cell distance or more
! than its radius in single time step as these will be captured in 
! cfnewvalues. So, first check if the particle is in its previous 
! location, and if not, add +/- 1 based on location. 

! assigning local aliases for particle i, j, k fluid grid indices
            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
    
            IF(XPOS >= XE(i-1) .and. XPOS < XE(i)) THEN 
               PIJK(L,1) = i
            ELSEIF(XPOS >= XE(i)) THEN 
               PIJK(L,1) = i+1
            ELSE 
               PIJK(L,1) = i-1 
            ENDIF 

            IF(YPOS >= YN(j-1) .and. YPOS < YN(j)) THEN 
               PIJK(L,2) = j
            ELSEIF(YPOS >= YN(j))THEN 
               PIJK(L,2) = j+1
            ELSE
               PIJK(L,2) = j-1
            ENDIF 

            IF(DIMN.EQ.2) THEN
               PIJK(L,3) = 1
            ELSE
               IF(ZPOS >= ZT(k-1) .and. ZPOS < ZT(k)) THEN
                  PIJK(L,3) = k
               ELSEIF(ZPOS >= ZT(k)) THEN 
                  PIJK(L,3) = k+1
               ELSE
                  PIJK(L,3) = k-1
               ENDIF 
            ENDIF 
         ENDIF   ! end of (if/else first_pass)


! Conducting checks on particle movement
! ---------------------------------------------------------------->>>

! if particle is not new or exiting then check its location
! if particle is not new or exiting it shouldn't be in a ghost cell
         IF (.NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
               
            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
! this ijk is still an old value as it has not been updated      
            IJK=PIJK(L,4)            

            IF(.NOT.MPPIC) THEN
               IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
                  IF(DMP_LOG) WRITE(UNIT_LOG,1007) L,'I',I,'X',&
                     XPOS,'X',DES_VEL_NEW(L,1)
                  CALL WRITE_DES_DATA
                  CALL MFIX_EXIT(MYPE)
               ENDIF
               IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
                  IF(DMP_LOG) WRITE(UNIT_LOG,1007) L,'J',J,'Y',&
                     YPOS,'Y',DES_VEL_NEW(L,2)
                  CALL WRITE_DES_DATA
                  CALL MFIX_EXIT(MYPE)
               ENDIF   
               IF ((DIMN.EQ.3) .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
                  IF(DMP_LOG) WRITE(UNIT_LOG,1007) L,'K',K,'Z',&
                     ZPOS,'Z',DES_VEL_NEW(L,3)
                  CALL WRITE_DES_DATA
                  CALL MFIX_EXIT(MYPE)
               ENDIF

            ELSE   ! mppic
! in MPPIC a particle can lie on the surface of the wall as only the 
! centers are tracked.
               IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
                  IF(I.EQ.IEND1+1 .AND. &
                    (XPOS >= XE(IEND1-1) .AND. XPOS <= XE(IEND1)) )THEN
! This could happen if the cell adjacent to the ghost cell is a cut-cell
! and due to some tolerance issues, the particle is not detected outside
! the system. This will be a rare occurence, but it may occur and there
! is no point in stalling the simulation here.
! Implementing an easy fix for now: delete this particle.
! To add more stuff later
! 1. re distribute particle's weight among other particles in the domain
!    so that mass is conserved
! 2. rather than deactivating the particle, reflect the particle
!    inside the domain using the ghost cell bc's instead of cut-face bc 
                    
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'I',I,'X',&
                        XPOS,DES_POS_OLD(L,1),'X',DES_VEL_NEW(L,1)
                     PIJK(L,1) = IEND1
                  ELSE
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'I',I,'X',&
                        XPOS,DES_POS_OLD(L,1),'X',DES_VEL_NEW(L,1),&
                        DES_VEL_OLD(L,1), CUT_CELL_AT(IJK), &
                        FLUID_AT(IJK)
                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE 
                  ENDIF
               ENDIF
               IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
                  IF(J.EQ.JEND1+1.AND.&
                    (YPOS >= YN(JEND1-1) .AND. YPOS <= YN(JEND1)) ) THEN
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'J',J,'Y',&
                        YPOS,DES_POS_OLD(L,2),'Y',DES_VEL_NEW(L,2)
                     PIJK(L,2) = JEND1
                  ELSE
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'J',J,'Y',&
                        YPOS,DES_POS_OLD(L,2),'Y',DES_VEL_NEW(L,2),&
                        DES_VEL_OLD(L,2), CUT_CELL_AT(IJK),&
                        FLUID_AT(IJK)
                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE
                  ENDIF
               ENDIF   
               IF ((DIMN.EQ.3) .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
                  IF(K.EQ.KEND1+1.AND.&
                    (ZPOS >= ZT(KEND1-1) .AND. ZPOS <= ZT(KEND1)) ) THEN
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'K',K,'Z',&
                        ZPOS, DES_POS_OLD(L,3),'Z',DES_VEL_NEW(L,3)
                     PIJK(L,3) = KEND1
                  ELSE
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'K',K,'Z',&
                        ZPOS,DES_POS_OLD(L,2),'Z',DES_VEL_NEW(L,3),&
                        DES_VEL_OLD(L,3), CUT_CELL_AT(IJK), &
                        FLUID_AT(IJK)
                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE 
                  ENDIF
               ENDIF
            ENDIF  ! end if/else(.NOT.MPPIC)

         ENDIF   ! end if(.not.pea(l,2) and .not.pea(l.3)) then
! End checks on particle movement         
! ----------------------------------------------------------------<<<

! Assigning PIJK(L,4) now that particles have been located on the fluid
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = FUNIJK(I,J,K)
         PIJK(L,4) = IJK



		 
         PINC(IJK) = PINC(IJK) + 1
		 
      ENDDO   ! end loop over L = 1,particles
!!$omp end parallel do 
!$omp end parallel 

! ----------------------------------------------------------------<<<  


      IF (FIRST_PASS) THEN
! Calling exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
! unclear why this needs to be called again if first_pass - update?
         CALL DES_PAR_EXCHANGE

         IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(3X,A)') &
            '<---------- END FIRST PASS PARTICLES_IN_CELL ----------'
      ENDIF



! Assigning the variable PIC(IJK)%p(:). For each computational fluid
! cell compare the number of current particles in the cell to what was
! in the cell previously. If different reallocate. Store the particle
! ids
! ---------------------------------------------------------------->>>
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)           &
!$omp private(ijk,npic) !schedule (guided,50)     
      DO IJK = IJKSTART3, IJKEND3

! checking all cells (including ghost cells); updating entering/exiting
! particle regions
         NPIC =  PINC(IJK)
         IF (ASSOCIATED(PIC(IJK)%p)) THEN
            IF (NPIC.NE.SIZE(PIC(IJK)%p)) THEN
               DEALLOCATE(PIC(IJK)%p)
               IF (NPIC.GT.0) ALLOCATE(PIC(IJK)%p(NPIC))
            ENDIF
         ELSE
            IF (NPIC.GT.0) ALLOCATE(PIC(IJK)%p(NPIC))
         ENDIF
      ENDDO
!$omp end parallel do 

      particle_count(:) = 1
      PC = 1
      DO L = 1, MAX_PIP
! exiting loop if reached max number of particles in processor
         IF(PC.GT.PIP) exit
! skipping indices with no particles (non-existent particles)         
         IF(.NOT.PEA(L,1)) CYCLE
! incrementing particle account when particle exists         
         pc = pc+1
! skipping ghost particles
         IF(PEA(L,4)) CYCLE
         IJK = PIJK(L,4)
         pos = particle_count(IJK)
         pic(IJK)%p(pos) = L
         particle_count(IJK) = particle_count(IJK) + 1
      ENDDO
! ----------------------------------------------------------------<<<


! calculating mean fields using either interpolation methods or a simple
! cell average technique. the void fraction is also calculated during
! this step
      IF(DES_INTERP_MEAN_FIELDS) THEN 
         CALL COMP_MEAN_FIELDS_INTERP
      ELSE
         CALL COMP_MEAN_FIELDS_ZERO_ORDER
      ENDIF

      IF(MPPIC) THEN 
         PIP = PIP - PIP_DEL_COUNT
         LPIP_DEL_COUNT_ALL(:) = 0
         LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT 
         CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL) 
         IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN 
            IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10)') &
               'TOTAL NUMBER OF PARTICLES OUSIDE DOMAIN IN PIC = ', &
               SUM(LPIP_DEL_COUNT_ALL(:))
            WRITE(UNIT_LOG,'(/,2x,A,2x,i10)') &
               'TOTAL NUMBER OF PARTICLES OUTSIDE DOMAIN IN PIC= ',&
               SUM(LPIP_DEL_COUNT_ALL(:))
            DO IPROC = 0, NUMPES-1 
               WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') &
               'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,&
               ' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
            ENDDO            
         ENDIF

         NNODES = NODESI*NODESJ*NODESK
         IF (NNODES.EQ.1 ) THEN
            EPG_MIN2 = MINVAL(EP_G(:))
            epg_min_loc = MINLOC(EP_G(:))
            IJK = epg_min_loc(1)
            I = I_OF(IJK)
            J = J_OF(IJK)
            K = K_OF(IJK)
            WRITE(*,1014) epg_min2, I_OF(IJK), j_of(ijk), k_of(ijk), &
            & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
            & PINC(IJK), cut_cell_at(ijk), fluid_at(ijk)
         ENDIF
      ENDIF
         

      FIRST_PASS = .FALSE.

 1000 FORMAT(3X,'---------- FROM PARTICLES_IN_CELL ---------->')
 1001 FORMAT(3X,'<---------- END PARTICLES_IN_CELL ----------') 

 1006 FORMAT(/1X,70('*')//&
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/&
         1X,70('*')/)

 1007 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/,&         
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/& 
         1X,70('*')/)

 1010     FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/,&         
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,/& 
         '-velocity new and old: ',2(ES17.9,4x),/& 
         ' CUT_CELL and FLUID AT IJK_OLD ?', 2(L2,2x),/& 
         ' Marking this particle as inactive',/&          
          1X,70('*')/)

 1011 FORMAT(/1X,70('*')//,&
         ' From: PARTICLES_IN_CELL: Particle recovered',&
         ' from ghost cell -',/,&         
         ' Message: Particle ',I8,' had moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,'-velocity: ',&
          ES17.9,/1X,70('*')/)

 1014 FORMAT( /, &
      &      10x,'EPGMIN NORMAL = ', 2x,g17.8, / & 
      &      10x,'EPG_MIN_LOCATION, I, J, K = ', 3(2x,i5),/, &
      &      10x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8),/ & 
      &      10x,'No of paricles in cell = ',I10, & 
      &      10x,'CUT CELL, FLUID AT IJK ?    ', 2(2x, L2)) !,/& 
!      &      1X,70('*')/)
     
      RETURN
      END SUBROUTINE PARTICLES_IN_CELL



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid 
      use desmpi
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! particle index
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase index
      INTEGER M, CM 
! ijk indices      
      INTEGER I, J, K, IJK
! Variable to distribute particle volume
      DOUBLE PRECISION ::  WTP
! 1 over volume of fluid cell      
      DOUBLE PRECISION :: OVOL
! total volume of mth phase solids in cell and 1 over that value      
      DOUBLE PRECISION SOLVOLINC(DIMENSION_3,DES_MMAX), OSOLVOL
! solids volume fraction of mth solids phase
      DOUBLE PRECISION EP_SM
! total solids volume fraction of continuum solids phases
      DOUBLE PRECISION SUM_EPS
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

      SOLVOLINC(:,:) = ZERO
      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO

! Add particle values (volume, velocity) to ongoing summations      
!!$      omp_start1=omp_get_wtime()	  
!!$omp single private(l,wtp,i,j,k,ijk,m) 
      PC = 1
      DO L = 1, MAX_PIP
! exiting loop if reached max number of particles in processor
         IF(PC.GT.PIP) EXIT
! skipping indices with no particles (non-existent particles) 
         IF(.NOT.PEA(L,1)) CYCLE
! incrementing particle account when particle exists 
         PC = PC + 1
! skipping ghost particles
         IF(PEA(L,4)) CYCLE
         
! assigning local aliases for particle i, j, k fluid grid indices
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = FUNIJK(I,J,K)
         M = PIJK(L,5)

         WTP = ONE
         IF(MPPIC) WTP = DES_STAT_WT(L) 
! adding particle volume to ongoing summation of solids volume
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) + PVOL(L)*WTP
! adding particle velocity to ongoing summation of solids velocity
         DES_U_S(IJK,M) = DES_U_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,1)*WTP
         DES_V_S(IJK,M) = DES_V_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,2)*WTP
         IF(DIMN.EQ.3) DES_W_S(IJK,M) = DES_W_S(IJK,M) + & 
            PVOL(L)*DES_VEL_NEW(L,3)*WTP
      ENDDO      ! end loop over L = 1,particles


! Calculate the cell average solids velocity, the bulk density,
! and the void fraction. 
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(ijk,i,j,k,cm,m,sum_eps,ep_sm,                    &
!$omp         osolvol,ovol)                                    
      DO IJK = ijkstart3, ijkend3
         IF(.NOT.FLUID_AT(IJK)) CYCLE 

         IF (.NOT.DES_CONTINUUM_HYBRID) THEN
            EP_G(IJK) = ONE   
         ELSE
! summing together total continuum solids volume
            SUM_EPS = ZERO
            DO CM = 1,SMAX
               SUM_EPS = SUM_EPS + EP_S(IJK,CM) 
            ENDDO
            EP_G(IJK) = ONE - SUM_EPS
         ENDIF  ! end if/else (.not.des_continuum_hybrid)


! calculating the cell average solids velocity for each solids phase         
          DO M = 1, DES_MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OSOLVOL = ONE/SOLVOLINC(IJK,M)   
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OSOLVOL
               IF(DIMN.EQ.3) THEN
                  DES_W_s(IJK,M) = DES_W_s(IJK,M)*OSOLVOL
               ENDIF
            ENDIF

! calculating the bulk density of solids phase m based on the total
! number of particles having their center in the cell
            IF(VOL(IJK).GT.0) THEN 
               OVOL = ONE/(VOL(IJK))
               DES_ROP_S(IJK,M) = DES_RO_S(M)*SOLVOLINC(IJK,M)*OVOL
            ENDIF
            
! calculating void fraction in fluid cell based on value of bulk density
! calculated above
            IF(DES_ROP_S(IJK,M) >= ZERO) THEN
! calculating solids volume fraction based on bulk density                    
               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
               IF(.NOT.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_SM
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)

! ep_g does not matter if granular flow simulation (i.e. no fluid)               
               IF(EP_G(IJK)<ZERO .AND.DES_CONTINUUM_COUPLED) THEN 
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1000)
                     WRITE(UNIT_LOG,1004) IJK, I_OF(IJK), J_OF(IJK), &
                        EP_SM, PINC(IJK)
                     WRITE(UNIT_LOG,1001)
                  ENDIF
                  
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDDO   ! end loop over M=1,DES_MMAX

      ENDDO     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do 

            
 1000 FORMAT(3X,'---------- FROM COMPUTE_MEAN_FIELDS_ZERO_ORDER ',&
             '---------->')
 1001 FORMAT(3X,'<--------- END COMPUTE_MEAN_FIELDS_ZERO_ORDER ',&
             '----------') 

 1004 FORMAT(5X,'WARNING: EP_G < 0 at IJK=', I10,' I=', I10, &
         ' J=', I10,/5X,'EP_S=', ES15.9, ' & PINC (number of ',&
         'particles in cell)= ',I10)


      END SUBROUTINE COMP_MEAN_FIELDS_ZERO_ORDER
   


!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE COMP_MEAN_FIELDS_INTERP

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE constant
      USE physprop
      USE fldvar
      USE run
      USE geometry
      USE indices
      USE bc
      USE compar
      USE sendrecv
      USE discretelement
      USE drag
      USE interpolation
      use desmpi 
      USE cutcell 
      USE mfix_pic
      USE mpi_utility
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: I, J, K, IJK, &
                 IPJK, IJPK, IJKP, IMJK, IJMK, IJKM,&
                 IPJPK, IPJKP, IJPKP, IPJPKP, &
                 II, JJ, KK              
      INTEGER :: I1, I2, J1, J2, K1, K2
      INTEGER :: IDIM, IJK2
      INTEGER :: ICUR, JCUR, KCUR, CUR_IJK
! indices used for interpolation stencil (unclear why IE, JN, KTP are
! needed)
      INTEGER :: IW, IE, JS, JN, KB, KTP
! i,j,k indices of the fluid cell the particle resides in minus 1 
! (e.g., shifted 1 in west, south, bottom direction)
      INTEGER, DIMENSION(3):: PCELL
! order of interpolation set in the call to set_interpolation_scheme
! unless it is re/set later through the call to 
! set_interpolation_stencil
      INTEGER :: ONEW
! constant whose value depends on dimension of system 
! avg_factor=0.250 (in 3D) or =0.50 (in 2D) 
      DOUBLE PRECISION :: AVG_FACTOR
! index of solid phase that particle NP belongs to      
      INTEGER :: M
! particle number index, used for looping      
      INTEGER :: NP, NINDX
! index to track accounted for particles 
      INTEGER :: PC 
! Statistical weight of the particle. Equal to one for DEM 
      DOUBLE PRECISION :: WTP
! one over the solids volume fraction 
      DOUBLE PRECISION :: OEPS

      DOUBLE PRECISION :: VOL_SURR
      DOUBLE PRECISION :: MASS_SOL1, MASS_SOL2 
! sum of mass_sol1 and mass_sol2 across all processors     
      DOUBLE PRECISION :: MASS_SOL1_ALL, MASS_SOL2_ALL
      
      DOUBLE PRECISION :: TEMP1(MAX_PIP)

! for error messages      
      INTEGER :: IER

      DOUBLE PRECISION :: JUNK_VAL(3)

      INTEGER :: COUNT_NODES_OUTSIDE, COUNT_NODES_INSIDE, &
                 COUNT_NODES_INSIDE_MAX, COUNT_TEMP
      double precision :: RESID_ROPS(DES_MMAX), &
                          RESID_VEL(DIMN, DES_MMAX)
      double precision :: NORM_FACTOR
!Handan Liu added on Jan 17 2013     
	  DOUBLE PRECISION, DIMENSION(2,2,2,3) :: gst_tmp,vst_tmp
	  DOUBLE PRECISION, DIMENSION(2,2,2) :: weight_ft
	  DOUBLE PRECISION :: desposnew(dimn)      
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
!-----------------------------------------------

! initializing            
      MASS_SOL1 = ZERO
      MASS_SOL2 = ZERO 
      MASS_SOL1_ALL = ZERO 
      MASS_SOL2_ALL = ZERO 
! avg_factor=0.25 (in 3D) or =0.5 (in 2D)
      AVG_FACTOR = 0.25D0*(DIMN-2) + 0.5D0*(3-DIMN)

! cartesian_grid related quantities
      IF(DIMN.EQ.2) COUNT_NODES_INSIDE_MAX = 4
      IF(DIMN.EQ.3) COUNT_NODES_INSIDE_MAX = 8

      DO IJK = IJKSTART3,IJKEND3
         DES_VEL_NODE(IJK, :, :) = ZERO
         DES_ROPS_NODE(IJK, :) = ZERO
         IF(FLUID_AT(IJK)) THEN
            DES_ROP_S(IJK,:) = zero 
            DES_U_S(IJK, :) = ZERO
            DES_V_S(IJK, :) = ZERO
            IF(DIMN.EQ.3) DES_W_S(IJK, :) = ZERO
         ENDIF
      ENDDO

      
! sets several quantities including interp_scheme, scheme, and 
! order and allocates arrays necessary for interpolation      
      CALL SET_INTERPOLATION_SCHEME(2)

!Handan Liu added on Jan 17 2013; again on June 2013	  
!$omp	parallel default(shared)								&
!$omp	private(IJK,I,J,K,PCELL,COUNT_NODES_INSIDE,II,JJ,KK,IW,	&
!$omp		IE,JS,JN,KB,KTP,ONEW,CUR_IJK,IPJK,IJPK,IPJPK,IJKP,	&
!$omp		IJPKP,IPJKP,IPJPKP,gst_tmp,vst_tmp,nindx,np,wtp,m,	&
!$omp		JUNK_VAL,desposnew,weight_ft,icur,jcur,kcur,		&	
!$omp		I1, I2, J1, J2, K1, K2, IDIM,IJK2,NORM_FACTOR,		&
!$omp		RESID_ROPS,RESID_VEL,COUNT_NODES_OUTSIDE)
!$omp do reduction(+:MASS_SOL1) reduction(+:DES_ROPS_NODE,DES_VEL_NODE)
      !IJKLOOP: DO IJK = IJKSTART3,IJKEND3	! Removed by Handan Liu
      DO IJK = IJKSTART3,IJKEND3

! Cycle this cell if not in the fluid domain or if it contains no
! particle/parcel
         IF(.NOT.FLUID_AT(IJK) .OR. PINC(IJK).EQ.0) CYCLE 

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         
         IF(.NOT.IS_ON_myPE_owns(I,J,K)) THEN 
! PINC array count reflects only the acutal particles and does not 
! include the ghost particles. Therefore, PINC can only be non-zero
! for scalar cells that belong to this processor. So this is just
! a sanity check that will ensure an error is flagged if the logic
! is broken elsewhere in the code.             
            IF(DMP_LOG) THEN 
               WRITE(UNIT_LOG,*) &
                  'Critical Error in compute_mean_fields_interp:', &
                  'This cell does not belong to this proc ', myPE 
               WRITE(UNIT_LOG,'(A,10(2x,i5))') &
                  'PINC in I, J, K, IJP, NP = ', & 
                  I_OF(IJK), J_OF(IJK), K_OF(IJK), PINC(IJK), &
                  PINC(IM_OF(IJK)), PINC(IP_OF(IJK))
            ENDIF
            WRITE(*,*) 'Critical Error in compute_mean_fields_interp:', &
               'This cell does not belong to this proc ', myPE, & 
               'Exiting the code, check the log file for details'
            CALL MFIX_EXIT(myPE)
         ENDIF

         PCELL(1) = I-1
         PCELL(2) = J-1
         PCELL(3) = (3-DIMN)*1+(DIMN-2)*(K-1) ! =K-1 (IN 3D) OR =1 (IN 2D)

! setup the stencil based on the order of interpolation and factoring in
! whether the system has any periodic boundaries. sets onew to order.
         CALL SET_INTERPOLATION_STENCIL(PCELL,IW,IE,JS,JN,KB,&
            KTP,INTERP_SCHEME,DIMN,ORDERNEW = ONEW) 

         COUNT_NODES_OUTSIDE = 0 		
! Computing/setting the geometric stencil 
         DO K = 1,(3-DIMN)*1+(DIMN-2)*ONEW
            DO J = 1,ONEW
               DO I = 1,ONEW
                  II = IW + I-1
                  JJ = JS + J-1
                  KK = KB + K-1
                  CUR_IJK = funijk(IMAP_C(II),JMAP_C(JJ),KMAP_C(KK))
                  IPJK    = funijk(IMAP_C(II+1),JMAP_C(JJ),KMAP_C(KK))
                  IJPK    = funijk(IMAP_C(II),JMAP_C(JJ+1),KMAP_C(KK))
                  IPJPK   = funijk(IMAP_C(II+1),JMAP_C(JJ+1),KMAP_C(KK))
                  IF(DIMN.EQ.3) THEN 
                     IJKP    = funijk(IMAP_C(II),JMAP_C(JJ),KMAP_C(KK+1))
                     IJPKP   = funijk(IMAP_C(II),JMAP_C(JJ+1),KMAP_C(KK+1))
                     IPJKP   = funijk(IMAP_C(II+1),JMAP_C(JJ),KMAP_C(KK+1))
                     IPJPKP  = funijk(IMAP_C(II+1),JMAP_C(JJ+1),KMAP_C(KK+1))
                  ENDIF

                  GST_TMP(I,J,K,1) = XE(II)
                  GST_TMP(I,J,K,2) = YN(JJ)
                  GST_TMP(I,J,K,3) = ZT(KK)*(DIMN-2) + DZ(1)*(3-DIMN)
                  VST_TMP(I,J,K,:) = ZERO
!===================================================================>>> Handan Liu	

                  IF(CARTESIAN_GRID) THEN 
                     IF(SCALAR_NODE_ATWALL(CUR_IJK)) COUNT_NODES_OUTSIDE = &
                     & COUNT_NODES_OUTSIDE + 1 
                  ENDIF
               ENDDO
            ENDDO
         ENDDO


! Calculate des_rops_node so des_rop_s, and in turn, ep_g can be updated 
!----------------------------------------------------------------->>>            
           
! looping through particles in the cell  
         DO NINDX = 1,PINC(IJK)
            NP = PIC(IJK)%P(NINDX)

! conducting some sanity checks             
            IF(PEA(NP, 4)) THEN 
               IF(DMP_LOG) THEN 
                  WRITE(UNIT_LOG,*) 'Encountered a ghost ',&
                     'particle in compute_mean_fields_interp'
                  WRITE(UNIT_LOG,'(A,10(2x,i5))') &
                     'PINC in I, J, K, IJP, NP = ', I_OF(IJK), &
                     J_OF(IJK), K_OF(IJK), PINC(IJK), PINC(IM_OF(IJK)),&
                     PINC(IP_OF(IJK))
               ENDIF
               WRITE(*,*) 'Encountered a ghost particle ',&
                  'in compute_mean_fields_interp'
               WRITE(*,'(A,10(2x,i5))') &
                  'PINC in I, J, K, IJP, NP = ', I_OF(IJK), J_OF(IJK), &
                  K_OF(IJK), PINC(IJK), PINC(IM_OF(IJK)), &
                  PINC(IP_OF(IJK))
               CALL MFIX_EXIT(myPE)
            ENDIF

            desposnew(:) = des_pos_new(np,:)
            call DRAG_INTERPOLATION(dimn,gst_tmp,vst_tmp,desposnew,JUNK_VAL,weight_ft)
!===================================================================>>> Handan Liu 

            M = PIJK(NP,5)
            WTP = ONE
            IF(MPPIC) WTP = DES_STAT_WT(NP)
            MASS_SOL1 = MASS_SOL1 + PMASS(NP)*WTP
            
            DO K = 1, (3-DIMN)*1+(DIMN-2)*ONEW
               DO J = 1, ONEW
                  DO I = 1, ONEW
! shift loop index to new variables for manipulation
                     II = IW + I-1
                     JJ = JS + J-1
                     KK = KB + K-1

! The interpolation is done using node. so one should use consistent
! numbering system in the current version imap_c is used instead of 
! ip_of or im_of                    
                     ICUR = IMAP_C(II)
                     JCUR = JMAP_C(JJ)
                     KCUR = KMAP_C(KK)
                     CUR_IJK = funijk(ICUR, JCUR, KCUR) 

                     !TEMP1 = WEIGHTP(I,J,K)*DES_RO_S(M)*PVOL(NP)*WTP
! Changed TEMP1 as an array TEMP1(NP) to ensure different TEMP1 
! for each particle in an ijk cell <June 18 2013>
                     TEMP1(NP) = WEIGHT_FT(I,J,K)*DES_RO_S(M)*PVOL(NP)*WTP					 
                     DES_ROPS_NODE(CUR_IJK,M) = DES_ROPS_NODE(CUR_IJK,M) + TEMP1(NP) 
                     DES_VEL_NODE(CUR_IJK, 1:DIMN,M) = &
                        DES_VEL_NODE(CUR_IJK, 1:DIMN,M) + TEMP1(NP)*DES_VEL_NEW(NP, 1:DIMN)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO   ! end do (nindx=1,pinc(ijk))
!-----------------------------------------------------------------<<<

         
! Only for cutcell cases may count_nodes_inside become less than its
! original set value. In such an event, the contribution of scalar nodes
! that do not reside in the domain is added to a residual array. This
! array is then redistribited equally to the nodes that are in the fluid
! domain. These steps are done to conserve mass.
!----------------------------------------------------------------->>>
         IF (CARTESIAN_GRID) THEN

! only for cartesian_grid will count_nodes_outside be modified from zero         
            COUNT_NODES_INSIDE = &
               COUNT_NODES_INSIDE_MAX - COUNT_NODES_OUTSIDE

            IF(COUNT_NODES_INSIDE.LT.COUNT_NODES_INSIDE_MAX) THEN 

! initializing                    
               RESID_ROPS(1:DES_MMAX) = ZERO 
               RESID_VEL(1:DIMN, 1:DES_MMAX) = ZERO

! Convention used to number node numbers
! i=1, j=2           i=2, j=2
!   _____________________
!   |                   |
!   |  I = 2, J = 2     |
!   |___________________|
! i=1, j=1           i=2, j=1
! setting indices based on convention
               I = I_OF(IJK)
               J = J_OF(IJK)
               K = K_OF(IJK)
               I1 = I-1
               I2 = I
               J1 = J-1
               J2 = J
! K1 = K if DIMN = 2, amd K1 = K-1 DIMN = 3
               K1 = (3-DIMN)*K+(DIMN-2)*(K-1)
               K2 = K
! first calculate the residual des_rops_node and des_vel_node that was
! computed on nodes that do not belong to the domain

               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IJK2 = funijk(II, JJ, KK)

                        IF(SCALAR_NODE_ATWALL(IJK2)) THEN 
                           RESID_ROPS(1:DES_MMAX) = &
                              RESID_ROPS(1:DES_MMAX) +&
                              DES_ROPS_NODE(IJK2,1:DES_MMAX)
                           DES_ROPS_NODE(IJK2,1:DES_MMAX) = ZERO 
                           DO IDIM = 1, DIMN
                              RESID_VEL(IDIM, 1:DES_MMAX) = &
                                 RESID_VEL(IDIM, 1:DES_MMAX) + & 
                                 DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX)
                              DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX) = ZERO
                           ENDDO
                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO
                  
! now add this residual equally to the remaining nodes
               NORM_FACTOR = ONE/REAL(COUNT_NODES_INSIDE)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IJK2 = funijk(II, JJ, KK)

                        IF(.NOT.SCALAR_NODE_ATWALL(IJK2)) THEN 
                           DES_ROPS_NODE(IJK2,1:DES_MMAX) = &
                              DES_ROPS_NODE(IJK2,1:DES_MMAX) + &
                              RESID_ROPS(1:DES_MMAX)*NORM_FACTOR
                           DO IDIM = 1, DIMN
                              DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX) = &
                                 DES_VEL_NODE(IJK2,IDIM, 1:DES_MMAX) + &
                                 RESID_VEL(IDIM, 1:DES_MMAX)*NORM_FACTOR
                           ENDDO
                        ENDIF

                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
         ENDIF   ! end if (cartesian_grid)
!-----------------------------------------------------------------<<<
      !ENDDO IJKLOOP   ! end do ijkloop (ijk=ijkstart3,ijkend3)
      ENDDO
!$omp end parallel   
 
! At the interface des_rops_node has to be added since particles 
! across the processors will contribute to the same scalar node. 
! sendrecv will be called and the node values will be added 
! at the junction. des_rops_node is altered by the routine when
! periodic boundaries are invoked
      CALL DES_ADDNODEVALUES_MEAN_FIELDS
      

! Now go from node to scalar center. Same convention as sketched
! earlier
!----------------------------------------------------------------->>>
! Explanation by RG: 08/17/2012 
! the approach used here is to make it general enough for cutcells to be
! included as well. The new changes do not alter earlier calculations
! but make the technique general as to include cartesian grid (cut-cell)
! simulations.
! Previously, the volume of the node (by array des_vol_node) was used to
! first scale the nodal values. Subsequently, these nodal values were
! equally added to compute the cell centered values for the scalar cell.

! Consider an internal node next to an edge node (a node adjacent to a
! boundary). In 2D, the volume of an edge node will be half that of an
! internal node. And an edge node will contribute double compared to
! an internal node to the value of the scalar cell they share. These
! calculations were previously accomplished via the variable volume of
! node.  Now this is accomplished by the ratio vol(ijk2)/vol_sur, where
! vol(ijk2) is the volume of the scalar cell in consideration and
! vol_sur is the sum of all the scalar cell volumes that have this node
! as the common node.

! looping over all fluid cells    
!Handan Liu added here on Feb. 28 2013
!$omp 	parallel do default(shared)		&
!$omp   private(K,J,I,IJK,I1,I2,J1,J2,K1,K2,	&
!$omp			II,JJ,KK,IJK2,M,VOL_SURR)   collapse (3)
      DO K = KSTART2, KEND1
         DO J = JSTART2, JEND1
            DO I = ISTART2, IEND1
               IF (DEAD_CELL_AT(I,J,K)) CYCLE  ! skip dead cells
               IJK = funijk(I,J,K)
               I1 = I
               I2 = I+1
               J1 = J
               J2 = J+1
               K1 = K
               K2 = (3-DIMN)*K+(DIMN-2)*(K+1)

               VOL_SURR = ZERO 

! looping over stencil points (node values)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells
                        IJK2 = funijk(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK))
                        IF(FLUID_AT(IJK2)) VOL_SURR = VOL_SURR+VOL(IJK2)
                     ENDDO
                  ENDDO
               ENDDO

! looping over stencil points (NODE VALUES)
               DO KK = K1, K2
                  DO JJ = J1, J2
                     DO II = I1, I2
                        IF (DEAD_CELL_AT(II,JJ,KK)) CYCLE  ! skip dead cells

                        IJK2 = funijk(IMAP_C(II), JMAP_C(JJ), KMAP_C(KK))
                        IF(FLUID_AT(IJK2).and.(IS_ON_myPE_wobnd(II, JJ, KK))) THEN 
! Since the data in the ghost cells is spurious anyway and overwritten during 
! subsequent send receives, do not compute any value here as this will 
! mess up the total mass value that is computed below to ensure mass conservation
! between Lagrangian and continuum representations 
                           DO M = 1, DES_MMAX
                              DES_ROP_S(IJK2, M) = DES_ROP_S(IJK2, M) + &
                                 DES_ROPS_NODE(IJK,M)*VOL(IJK2)/VOL_SURR
                              DES_U_S(IJK2, M) = DES_U_S(IJK2, M) + & 
                                 DES_VEL_NODE(IJK, 1, M)*VOL(IJK2)/VOL_SURR
                              DES_V_S(IJK2, M) = DES_V_S(IJK2, M) + & 
                                 DES_VEL_NODE(IJK, 2, M)*VOL(IJK2)/VOL_SURR
                              IF(DIMN.eq.3) DES_W_S(IJK2, M) = DES_W_S(IJK2, M) + & 
                                 DES_VEL_NODE(IJK, 3, M)*VOL(IJK2)/VOL_SURR
                           ENDDO
                        ENDIF
                     ENDDO  ! end do (ii=i1,i2)
                  ENDDO  ! end do (jj=j1,j2)
               ENDDO  ! end do (kk=k1,k2)

            ENDDO   ! end do (i=istart2,iend1)
         ENDDO   ! end do (j=jstart2,jend1)
      ENDDO   ! end do (k=kstart2,kend1)
		  
!-----------------------------------------------------------------<<<      

!Handan Liu added here on Feb. 28 2013
!$omp 	parallel do default(shared)		&
!$omp   private(K,J,I,IJK,M) reduction(+:MASS_SOL2)      
      DO IJK = IJKSTART3, IJKEND3
         IF(.not.FLUID_AT(IJK)) cycle
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
       
         DO M = 1, DES_MMAX
            IF(DES_ROP_S(IJK, M).GT.ZERO) THEN 
               DES_U_S(IJK, M) = DES_U_S(IJK,M)/DES_ROP_S(IJK, M)
               DES_V_S(IJK, M) = DES_V_S(IJK,M)/DES_ROP_S(IJK, M)
               IF(DIMN.eq.3) DES_W_S(IJK, M) = DES_W_S(IJK,M)/DES_ROP_S(IJK, M)
               
! Finally divide by scalar cell volume to obtain \eps * \rho_s
               DES_ROP_S(IJK, M) = DES_ROP_S(IJK, M)/VOL(IJK)

! Note that for cut-cell, it is important to check both fluid_at and
! is_on_mype_wobnd.  Fluid_at is not enough as it is shared between procs
! and fluid_at(ijkend3) might be true when in fact it does not belong to
! that proc 
               IF(IS_ON_myPE_wobnd(I, J, K)) MASS_SOL2 = MASS_SOL2 + & 
               DES_ROP_S(IJK, M)*VOL(IJK)
               
               !WRITE(*,*) 'ROP_S = ', DES_ROP_S(IJK, M), DES_ROP_S(IJK, M)/DES_RO_S(M)

            ENDIF
         ENDDO   ! end loop over M=1,DES_MMAX
      ENDDO                     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do

! RG: Aug, 20, 2012. 
! ------------------------------------------------------------------- 
! Decoupling the computation of ep_g and rop_g from the computation
! of des_rop_s. In the current implementation, des_rop_s (or rop_s) 
! and epg and also rop_g are first calculated and then sent 
! received. However, ep_g (and also rop_g) are simply functions of 
! des_rop_s. Therefore, once des_rop_s is sent received, ep_g and 
! rop_g can be correctly computed even for the ghost cells; thus, not
! requiring additional send_recv calls for ep_g and rop_g.       

! DES_ROP_S is not exchanged across the interface. 
! ROP_S and EP_G are done at the end of the DEM time step in 
! des_time_march.f. 
! -------------------------------------------------------------------

      IF (.NOT.MPPIC) THEN
! According to the current DEM implementation, the mean fields are
! communicated (i.e., sent received) following the inner DEM 
! time steps. Since the number of these inner time steps could be large 
! in DEM and the cell centered DES mean fields (such as des_rop_s, 
! des_u_s, etc.) in the ghost cells (across processor boundaries) 
! are not used during these substeps, the sent receive of these mean
! fields at the end of DEM substeps does not result in any error. 

! It should be noted above that the ep_g values in the ghost cells 
! are not correct, and if any future development requires the correct
! value of ep_g in ghost cells during DEM substeps, a send_recv
! of des_rop_s should be done before calling the following routine. 
! This will impose a substantial penalty on computational time since 
! communication overheads will go up. 
         CALL COMP_EPG_ROP_G_FROM_ROP_S


      ELSE
      
! share the des_rop_s array across the processors          
         CALL SEND_RECV(DES_ROP_S,2)
! compute the arrays epg and rop_g. This negates the need for further 
! communication of these arrays.          
         CALL COMP_EPG_ROP_G_FROM_ROP_S

! Now calculate Eulerian mean velocity fields like U_S, V_S, and W_S. 
         CALL SEND_RECV(DES_U_S,2)
         CALL SEND_RECV(DES_V_S,2) 
         IF(DIMN.EQ.3) CALL SEND_RECV(DES_W_S,2) 

! The Eulerian velocity field is used to set up the stencil to interpolate
! mean solid velocity at the parcel's location. DES_U_S could have also been
! used, but that also would have require the communication at this stage.
! The final interpolated value does not change if the stencil is formed by
! first obtaining face centered Eulerian velocities (U_S, etc.) 
! and then computing the node velocities from them or directly computing 
! the node velocities from cell centered average velocity field (DES_U_S, 
! etc.). We are using the first approach as it is more natural to set 
! BC's on solid velocity field in the face centered represenation (U_S,
! etc.)

         IF(.NOT.CARTESIAN_GRID) THEN 
            CALL MPPIC_COMP_EULERIAN_VELS_NON_CG
         ELSE
            CALL MPPIC_COMP_EULERIAN_VELS_CG
         ENDIF
      ENDIF   ! end if (.not.mppic)
     
! turn on the below statements to check if the mass is conserved 
! between discrete and continuum representations. Should be turned to
! false for any production runs.      
      IF(DES_REPORT_MASS_INTERP) THEN 
         CALL GLOBAL_SUM(MASS_SOL1, MASS_SOL1_ALL)
         CALL GLOBAL_SUM(MASS_SOL2, MASS_SOL2_ALL)
         if(myPE.eq.pe_IO) THEN
            WRITE(*,'(10x,A,4(2x,g17.8))') & 
                 'SOLIDS MASS DISCRETE AND CONTINUUM =  ', & 
                 MASS_SOL1_ALL, MASS_SOL2_ALL
         ENDIF
      ENDIF
      END SUBROUTINE COMP_MEAN_FIELDS_INTERP



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!      

      SUBROUTINE COMP_EPG_ROP_G_FROM_ROP_S

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param
      USE param1
      USE parallel
      USE constant
      USE geometry
      USE funits
      USE indices
      USE compar
      USE physprop
      USE fldvar
      USE discretelement
      USE cutcell 
      USE mfix_pic
      implicit none 
!-----------------------------------------------
! Local variables
!-----------------------------------------------
! general i, j, k indices
      INTEGER :: IJK
! solids phase indices
      INTEGER :: CM, M
! solids volume fraction of mth solids phase
      DOUBLE PRECISION EP_SM
! total solids volume fraction of continuum solids phases
      DOUBLE PRECISION SUM_EPS
!-----------------------------------------------
! Include statement functions
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'
!-----------------------------------------------

!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(ijk,cm,m,sum_eps,ep_sm)    
      DO IJK = IJKSTART3, IJKEND3
         IF(.NOT.FLUID_AT(IJK)) CYCLE
         
         IF (.NOT.DES_CONTINUUM_HYBRID) THEN
            EP_G(IJK) = ONE   
         ELSE
! summing together total continuum solids volume
            SUM_EPS = ZERO
            DO CM = 1,SMAX
               SUM_EPS = SUM_EPS + EP_S(IJK,CM) 
            ENDDO
            EP_G(IJK) = ONE - SUM_EPS
         ENDIF  ! end if/else (.not.des_continuum_hybrid)
         
         DO M = 1, DES_MMAX
! calculating void fraction in fluid cell based on value of bulk density
            IF(DES_ROP_S(IJK, M).GT.ZERO) THEN 
               EP_SM = DES_ROP_S(IJK,M)/DES_RO_S(M)
               IF(.NOT.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_SM
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)

! ep_g does not matter if granular flow simulation (i.e. no fluid)
              IF(EP_G(IJK)<ZERO .AND.DES_CONTINUUM_COUPLED) THEN 
                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1004) EP_G(IJK), IJK, I_OF(IJK), & 
                        J_OF(IJK), K_OF(IJK), EP_SM, PINC(IJK), &
                        CUT_CELL_AT(IJK)
                  ENDIF
                  IF(mype.eq.pe_IO) WRITE(*,1004) EP_G(IJK), IJK, &
                     I_OF(IJK), J_OF(IJK), K_OF(IJK), EP_SM, &
                     PINC(IJK), CUT_CELL_AT(IJK)
                  
                  IF(CARTESIAN_GRID) THEN 
                     CALL WRITE_DES_DATA
                     CALL WRITE_VTK_FILE
                  ENDIF
                     
                  CALL MFIX_EXIT(myPE)
               ENDIF
            ENDIF
         ENDDO                  ! end loop over M=1,DES_MMAX
      ENDDO                     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do 
 1004 FORMAT(/5X, 'Message from comp_epg_rop_g_from_rop_s', & 
      /,5X,'Warning, EP_G = ', g17.8, 2x, 'LT Zero at IJK', I20, & 
      /,5X,'I,J,K = ', I10, 2X, I10, 2x, I10, & 
      /,5X,'EP_S=', ES15.9, &
      /,5X,'PINC (number of particles in cell)= ',I10, & 
      /,5X,'Cut cell ? ', L2,/)


      END SUBROUTINE COMP_EPG_ROP_G_FROM_ROP_S
      
