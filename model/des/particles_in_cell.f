!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  PARTICLES_IN_CELL                                     C
!
!  Purpose: DES - Finding the fluid computational cell in which      
!           a particle lies, to calculte void fraction and also       
!           the volume averaged solids velocity of the cell            
!
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C 
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Removed the separate volume definitions and added pic     C
!            array formulation and bed height calculation.             C
!                                                                      C
!            For parallel processing indices are altered and changes   C
!            to variables related to desgridsearch are made; steps     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

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
! 1 over volume of fluid cell      
      DOUBLE PRECISION :: OVOL
! total volume of mth phase solids in cell and 1 over that value      
      DOUBLE PRECISION SOLVOLINC(DIMENSION_3,MMAX), OSOLVOL
! variables that count/store the number of particles in ijk cell
      INTEGER:: npic, pos
! particle count in ijk fluid cell 
      INTEGER, DIMENSION(DIMENSION_3):: particle_count
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
! Variable to distribute particle volume
      DOUBLE PRECISION ::  WTP
! count for number of particles that were found in the ghost cell and
! hence removed       
      INTEGER :: PIP_DEL_COUNT, IMJK, ijmk, ijkm
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1)
! number of i,j,k nodes (parallelization)
      INTEGER :: NNODES
! IER for error reporting
      INTEGER IER 
!$      double precision omp_start1, omp_end1
!$      double precision omp_get_wtime	      
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
      SOLVOLINC(:,:) = ZERO
      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO

      IF(FIRST_PASS) THEN
         IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(3X,A)') &
            '---------- START FIRST PASS PARTICLES_IN_CELL ---------->'
      ENDIF

! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
      CALL DES_PAR_EXCHANGE


! Assigning PIJK(L,1), PIJK(L,2) and PIJK(L,3) the i, j, k indices 
! of particle L (locating particle on fluid grid). Also determine
! composite ijk index.  If first_pass, also assigning PIJK(L,5) the
! solids phase index of particle.
! ---------------------------------------------------------------->>>

!$omp parallel do default(shared)                               &
!$omp private(l,m,xpos,ypos,zpos,i,j,k,ijk) schedule (guided,50)  
      DO L = 1, MAX_PIP
! skip particles that do not exist
         IF(.NOT.PEA(L,1)) CYCLE
! skip ghost particles
         IF(PEA(L,4)) CYCLE
        
         IF(FIRST_PASS) THEN 

! Determining which solids phase a particle belongs based on matching each
! particle diameter and density to a solids phase particle diameter and
! density                 
            DO M = 1, MMAX
               IF(ABS(2.0d0*DES_RADIUS(L)-D_P0(M)).LT.SMALL_NUMBER.AND. &
               ABS( RO_Sol(L)-RO_S(M)).LT.SMALL_NUMBER) THEN
               PIJK(L,5) = M 
               EXIT
               ENDIF
            ENDDO
            IF(PIJK(L,5).EQ.0) THEN
               IF(DMP_LOG) THEN
                  write(unit_log,'(5X,A,A,I10)') &
                  'Problem determining the solids ',&
                  'association for particle: ',L
                  write(unit_log,'(7X,A,(ES15.9))') &
                  'Particle position = ', DES_POS_NEW(L,:)
                  write(unit_log,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                  'Particle diameter = ', 2.0*DES_RADIUS(L),&
                  'and D_P0(1:MMAX)= ', D_P0(1:MMAX)
                  write(unit_log,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                  'Particle density = ', Ro_Sol(L), &
                  'and RO_s(1:MMAX) = ', RO_S(1:MMAX)
               ENDIF
            ENDIF

! Brute force technique to determine the particle locations in the
! Eulerian fluid grid
            XPOS = DES_POS_NEW(L,1)
            YPOS = DES_POS_NEW(L,2)
            IF (DIMN .EQ. 3) THEN
               ZPOS = DES_POS_NEW(L,3)
            ENDIF
! In case of parallel processing the cell index is from istart1 to iend1
! XE will have values from istart2-1 to iend2             
! Also particles may be found in ghost cells in the event of particles
! particles entering/exiting the domain (this might occur for restart) 
            do I = istart2, iend2
               if( XPOS >= XE(I-1) .and. XPOS < XE(I)) then
                  PIJK(L,1) = i
                  exit 
               endif
            enddo
            do J = jstart2,jend2
               if(YPOS >= YN(J-1) .and. YPOS < YN(J)) then
                  PIJK(L,2) = J
                  exit
               endif
            enddo
            if(DIMN.EQ.2) then
               K=1
               PIJK(L,3) = 1
            else
               do K = kstart2,kend2
                  if(ZPOS >= ZT(K-1) .and. ZPOS < ZT(K)) then 
                     PIJK(L,3) = k
                     exit
                  endif
               enddo
            endif

         ELSE    ! if not first_pass

! Incremental approach to determine the new location of the particles  
! Note a particle cannot move by more than one cell distance or more than 
! its radius in single time step as these will be captured in cfnewvalues.
! So, first check if the particle is in its previous location, and if
! not, add +/- 1 based on location. 

! The loop is modified b.c using I-2 index for particle entering which has I=1 
! might result in segmentation error
            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
                        
            XPOS = DES_POS_NEW(L,1) 
            YPOS = DES_POS_NEW(L,2)
            IF (DIMN .eq. 3) THEN
               ZPOS = DES_POS_NEW(L,3)
            ENDIF
 
            if(XPOS >= XE(i-1) .and. XPOS < XE(i)) then 
               PIJK(L,1) = i
            elseif(XPOS >= XE(i)) then 
               PIJK(L,1) = i+1
            else 
               PIJK(L,1) = i-1 
            end if 

            if(YPOS >= YN(j-1) .and. YPOS < YN(j)) then 
               PIJK(L,2) = j
            elseif(YPOS >= YN(j))then 
               PIJK(L,2) = j+1
            else
               PIJK(L,2) = j-1
            end if 

            if(dimn.eq.2) then
               PIJK(L,3) = 1
            else
               if(ZPOS >= ZT(k-1) .and. ZPOS < ZT(k)) then
                  PIJK(L,3) = k
               elseif(ZPOS >= ZT(k)) then 
                  PIJK(L,3) = k+1
               else
                  PIJK(L,3) = k-1
               end if 
            endif 


         ENDIF   ! end of (if/else first_pass)

! Conduct checks on particle movement
! ---------------------------------------------------------------->>>

! if particle is not new or exiting then check its location
! if particle is not new or exiting it shouldn't be in a ghost cell
         IF (.NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
               
            I = PIJK(L,1)
            J = PIJK(L,2)
            K = PIJK(L,3)
! this ijk is still an old value as it has not been updated      
            IJK=PIJK(L,4)            

            IF (I.GT.IEND1 .OR. I.LT.ISTART1) THEN
               IF(MPPIC) THEN
! in MPPIC a particle can lie on the surface of the wall as only the 
! centers are tracked.
                  IF(I.EQ.IEND1+1.AND. &
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
                     WRITE(*,*) 'PIC: PIJK = ', PIJK(L,:)
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'I',I,'X',&
                        XPOS,DES_POS_OLD(L,1),'X',DES_VEL_NEW(L,1),&
                        DES_VEL_OLD(L,1), MPPIC, CARTESIAN_GRID, &
                        CUT_CELL_AT(IJK), FLUID_AT(IJK)
                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE 
                  ENDIF
               ELSE   ! not mppic
                  IF(DMP_LOG) WRITE(UNIT_LOG,1007) L,'I',I,'X',&
                     XPOS,'X',DES_VEL_NEW(L,1)
                  CALL MFIX_EXIT(MYPE)
               ENDIF   ! end if/else mppic
            ENDIF

            IF(J.GT.JEND1 .OR. J.LT.JSTART1) THEN
               IF(MPPIC) THEN
                  IF(J.EQ.JEND1+1.AND.&
                    (YPOS >= YN(JEND1-1) .AND. YPOS <= YN(JEND1)) ) THEN
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'J',J,'Y',&
                        YPOS,DES_POS_OLD(L,2),'Y',DES_VEL_NEW(L,2)
                     PIJK(L,2) = JEND1
                  ELSE
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'J',J,'Y',&
                        YPOS,DES_POS_OLD(L,2),'Y',DES_VEL_NEW(L,2),&
                        DES_VEL_OLD(L,2), MPPIC, CARTESIAN_GRID,&
                        CUT_CELL_AT(IJK), FLUID_AT(IJK)
                        PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE
                  ENDIF
               ELSE   ! not mppic
                  IF(DMP_LOG) WRITE(UNIT_LOG,1007) L,'J',J,'Y',&
                     YPOS,'Y',DES_VEL_NEW(L,2)
                  CALL MFIX_EXIT(MYPE)
               ENDIF   ! end if/else mppic
            ENDIF 
   
            IF ((DIMN.EQ.3) .AND. (K.GT.KEND1 .OR. K.LT.KSTART1)) THEN
               IF(MPPIC) THEN                 
                  IF(K.EQ.KEND1+1.AND.&
                    (ZPOS >= ZT(KEND1-1) .AND. ZPOS <= ZT(KEND1)) ) THEN
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'K',K,'Z',&
                        ZPOS, DES_POS_OLD(L,3),'Z',DES_VEL_NEW(L,3)
                     PIJK(L,3) = KEND1
                  ELSE
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'K',K,'Z',&
                     ZPOS,DES_POS_OLD(L,2),'Z',DES_VEL_NEW(L,3),&
                     DES_VEL_OLD(L,3), MPPIC, CARTESIAN_GRID, &
                     CUT_CELL_AT(IJK), FLUID_AT(IJK)
                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE 
                  ENDIF
               ELSE
                  IF(DMP_LOG) WRITE(UNIT_LOG,1007) L,'K',K,'Z',&
                     ZPOS,'Z',DES_VEL_NEW(L,3)
                  CALL MFIX_EXIT(MYPE)
               ENDIF
            ENDIF

         ENDIF   ! end if(.not.pea(l,2) and .not.pea(l.3)) then
! End checks on particle movement         
! ----------------------------------------------------------------<<<         

      ENDDO   ! end loop over l = 1,particles
!$omp end parallel do 
! ----------------------------------------------------------------<<<  


! Assigning PIJK(L,4) now that particles have been located on the fluid
! grid and summing accumulation terms for calculating a cell average 
! solids phase velocity
! ---------------------------------------------------------------->>>
      PC = 1
      DO L = 1, MAX_PIP
         IF(PC.GT.PIP) EXIT
         IF(.NOT.PEA(L,1)) CYCLE
         PC = PC + 1
         IF(PEA(L,4)) CYCLE
         
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
         IJK = FUNIJK(I,J,K)
         PIJK(L,4) = IJK        

         PINC(IJK) = PINC(IJK) + 1
         M = PIJK(L,5)

         WTP = ONE
         IF(MPPIC) WTP = DES_STAT_WT(L) 
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) +  PVOL(L)*WTP
         DES_U_S(IJK,M) = DES_U_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,1)*WTP
         DES_V_S(IJK,M) = DES_V_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,2)*WTP
         IF(DIMN.EQ.3) DES_W_S(IJK,M) = DES_W_S(IJK,M) + &
            PVOL(L)*DES_VEL_NEW(L,3)*WTP

      ENDDO   ! end loop over L = 1,particles
! ----------------------------------------------------------------<<<


! Call exchange particles - this will exchange particle crossing
! boundaries as well as updates ghost particles information
! unclear why this needs to be called again if first_pass - update?      
      IF (FIRST_PASS) CALL DES_PAR_EXCHANGE

      IF (FIRST_PASS) THEN
         IF(DMP_LOG.AND.DEBUG_DES) WRITE(UNIT_LOG,'(3X,A)') &
            '---------- END FIRST PASS PARTICLES_IN_CELL ---------->'
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
      ENDIF


! Assign/allocate the variable PIC(IJK)%p(:). For each cell compare 
! the number of current particles in the cell to what was in the
! cell previously.  If different reallocate.  Store the particle ids
! ---------------------------------------------------------------->>>
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)           &
!$omp private(ijk,npic) !schedule (guided,50)     
      DO IJK = IJKSTART3, IJKEND3

! check all cells (including ghost cells); update entering/exiting 
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
! exit loop if reached max number of particles in processor
         IF(PC.GT.PIP) exit
! skip indices with no particles (non-existent particles)         
         IF(.NOT.PEA(L,1)) CYCLE
! incrementing particle account when particle exists         
         pc = pc+1
! skip ghost particles
         IF(PEA(L,4)) CYCLE
         IJK = PIJK(L,4)
         pos = particle_count(IJK)
         pic(IJK)%p(pos) = L
         particle_count(IJK) = particle_count(IJK) + 1
      ENDDO
! ----------------------------------------------------------------<<<


! if using drag interpolation then the bulk density and corresponding
! void fraction should be updated based on the interpolated particle
! position as is used in the drag calculation.  hence the bulk density
! should be calculated/updated here using interpolation routines

      
! Calculate the cell average solids velocity, the bulk density (if not
! des_interp_on and not first_pass), and the void fraction. 
! ---------------------------------------------------------------->>>
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(ijk,i,j,k,m,osolvol,ovol) !schedule (guided,50)     
      DO IJK = IJKSTART3,IJKEND3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)

         IF(.NOT.FLUID_AT(IJK).OR..NOT.IS_ON_myPE_owns(I,J,K)) CYCLE 
! Note that for cut-cell, it is important to check both fluid_at and
! is_on_mype_owns.  Fluid_at is not enough as it is shared between procs
! and fluid_at(ijkend3) might be true when in fact it does not belong to
! that proc 

         ROP_SO(IJK,:) = ZERO 
         EP_G(IJK) = ONE
         
         DO M = 1, MMAX

! calculating average solids velocity in a fluid cell         
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OSOLVOL = ONE/SOLVOLINC(IJK,M)   
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OSOLVOL
               IF(DIMN.EQ.3) THEN
                  DES_W_s(IJK,M) = DES_W_s(IJK,M)*OSOLVOL
               ENDIF
            ENDIF

! calculating bulk density in a fluid cell            
! this bulk density is based simply on locating a particle center in a
! corresponding fluid cell.  So this calculation creates discrepancies
! if interpolated drag calculations are used.  As the interpolated drag
! will be using oudated void fraction information. In that case the 
! bulk density (and void fraction) should be updated based on the 
! interpolation scheme
            IF(VOL(IJK).GT.0) THEN 
               OVOL = ONE/(VOL(IJK))
               IF(FIRST_PASS .OR. &
                 ((.NOT.FIRST_PASS).AND.(.NOT.DES_INTERP_ON))) THEN
                  ROP_S(IJK,M)  = RO_S(M)*SOLVOLINC(IJK,M)*OVOL
               ENDIF
            ENDIF

! calculating void fraction in fluid cell based on solids bulk density       
            IF(ROP_S(IJK,M) > ZERO) THEN
               ROP_SO(IJK,M)  = ROP_S(IJK,M) 
               IF(.not.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)

! ep_g does not matter if granular flow simulation (i.e. no fluid)               
               IF(EP_G(IJK).LT.ZERO .AND.DES_CONTINUUM_COUPLED.AND.&
                 (.NOT.MPPIC)) THEN 

                  IF(DMP_LOG) THEN
                     WRITE(UNIT_LOG,1000)
                     WRITE(UNIT_LOG,1004) IJK, I_OF(IJK), J_OF(IJK), &
                        EP_S(IJK,M), PINC(IJK)
                     WRITE(UNIT_LOG,1001)
                     WRITE(UNIT_LOG,*) 'Cut cell? ', cut_cell_at(IJK)
                  ENDIF
                  
                  CALL MFIX_EXIT(myPE)
               ENDIF

            ENDIF
         ENDDO   ! end loop over M=1,MMAX
      ENDDO     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do 
! ----------------------------------------------------------------<<<

      FIRST_PASS = .FALSE.

      IF(MPPIC) THEN 
         CALL MPPIC_COMPUTE_MEAN_FIELDS2
         NNODES = NODESI*NODESJ*NODESK
         IF (NNODES.EQ.1 .AND. DMP_LOG) THEN
! writing out useful information on minimum epg
! Currently, it is only valid for serial runs.                 
            WRITE(UNIT_LOG,1014) I_OF(IJK), j_of(ijk), k_of(ijk), &
               xe(I)-0.5*dx(i), yn(J)-0.5*DY(J), zt(K)-0.5*DZ(K), & 
               cut_cell_at(ijk),fluid_at(ijk)
         ENDIF 
      ENDIF 
      
 1000 FORMAT(3X,'---------- FROM PARTICLES_IN_CELL ---------->')
 1001 FORMAT(3X,'<---------- END PARTICLES_IN_CELL ----------') 

 1004 FORMAT(5X,'WARNING: EP_G < 0 at IJK=', I10,' I=', I10, &
         ' J=', I10,/5X,'EP_S=', ES15.9, ' & PINC (number of ',&
         'particles in cell)= ',I10)

 1007 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/,&         
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/& 
         1X,70('*')/)


 1010     FORMAT(/1X,70('*')//,&
         ' From: PARTICLES_IN_CELL -',/,&         
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,/& 
         '-velocity new and old: ',ES17.9,/& 
         ' MPPIC and Cartesian Grid ?', 2(L2,2x),/& 
         ' CUT_CELL and FLUID AT IJK ?', 2(L2,2x),/& 
         ' Marking this particle as inactive',/&          
          1X,70('*')/)

 1011     FORMAT(/1X,70('*')//,&
         ' From: PARTICLES_IN_CELL: Particle recovered',&
         ' from ghost cell -',/,&         
         ' Message: Particle ',I8,' had moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,'-velocity: ',&
          ES17.9,/1X,70('*')/)

 1014 FORMAT(10x,'EPG_MIN_LOCATION, I, J, K = ', 3(2X,I5),/,10X, &
         'XMID, YMID, ZMID FOR CELL = ', 3(2X,G17.8),/,10X, & 
         'CUT CELL, FLUID AT IJK ?    ', 2(2X, L2), /)

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL


