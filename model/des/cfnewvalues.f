!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,       
!           position, angular velocity etc                            
!
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C 
!
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

      USE param
      USE param1
      USE parallel
      USE matrix
      USE scales
      USE constant
      USE compar
      USE sendrecv
      USE physprop
      USE fldvar
      USE visc_g
      USE rxns
      USE run
      USE geometry
      USE indices
      USE drag
      USE discretelement
      USE des_bc

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, I
      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST

! index to track accounted for particles  
      INTEGER PC 

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

!-----------------------------------------------
! Functions 
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------      
      PC = 1
      DO L = 1, MAX_PIS
         IF(PC.GT. PIS) EXIT
         IF(.NOT.PEA(L,1)) CYCLE


         DES_LOC_DEBUG = .FALSE.

! If a particle is classified as new, then forces are ignored. 
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.PEA(L,2))THEN 
            FC(L, :) = FC(L,:)/PMASS(L) + GRAV(:)
         ELSE 
            FC(L,:) = ZERO
            TOW(L,:) = ZERO         
         ENDIF


! Advance particle position, velocity
         IF (TRIM(DES_INTG_METHOD) .EQ. 'EULER') THEN 
! first-order method              
            DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + FC(L,:)*DTSOLID
            DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
               DES_VEL_NEW(L,:)*DTSOLID 
! following is equivalent to x=xold + vold*dt + 1/2acc*dt^2
!         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + 0.5d0*&
!             (DES_VEL_NEW(L,:)+DES_VEL_OLD(L,:))*DTSOLID 
            OMEGA_NEW(L,:)   = OMEGA_OLD(L,:) + TOW(L,:)*OMOI(L)*DTSOLID

         ELSEIF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') THEN
! T.Li:  second-order Adams-Bashforth scheme
            DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + 0.5d0*&
               ( 3.d0*DES_VEL_OLD(L,:)-DES_VEL_OOLD(L,:) )*DTSOLID
            DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + 0.5d0*&
               ( 3.d0*FC(L,:)-DES_ACC_OLD(L,:) )*DTSOLID
            OMEGA_NEW(L,:)   =  OMEGA_OLD(L,:) + 0.5d0*&
               ( 3.d0*TOW(L,:)*OMOI(L)-ROT_ACC_OLD(L,:) )*DTSOLID
            DES_ACC_OLD(L,:) = FC(L,:)
            ROT_ACC_OLD(L,:) = TOW(L,:)*OMOI(L)
         ENDIF


! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so, 
! make sure that neighbor is called in des_time_march
         IF(.NOT.DO_NSEARCH) THEN
            D(:) = DES_POS_NEW(L,:) - PPOS(L,:)
            DIST = SQRT(DES_DOTPRDCT(D,D))
            NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*&
               DES_RADIUS(L)
            IF(DIST.GE.NEIGHBOR_SEARCH_DIST) DO_NSEARCH = .TRUE.
         ENDIF


! Check if the particle has moved a distance greater than or equal to 
! its radius during one solids time step. if so, call stop
         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         DIST = SQRT(DES_DOTPRDCT(D,D))
         IF(DIST.GE.DES_RADIUS(L)) THEN
            WRITE(*,1002) L, DIST, DES_RADIUS(L)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'old particle pos = ', DES_POS_OLD(L,:)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'new particle pos = ', DES_POS_NEW(L,:)
            WRITE(*,'(5X,A,3(ES17.9))')&
               'new particle vel = ', DES_VEL_NEW(L,:) 
            WRITE(*,1003)
            STOP
         ENDIF


! Lees & Edwards BC or Periodic BC.  Note current implementation of the
! LE BC implies periodic treatment of remaining boundaries but is is
! handled in the LE BC routine
         IF(DES_LE_BC) THEN
            CALL DES_LEES_EDWARDS_BC(L)
         ELSEIF(DES_PERIODIC_WALLS) THEN
            CALL DES_PERIODIC_BC(L)
         ENDIF


! Warning message for particles moving into ghost cells:
! Note that if this occurs then the particle_in_cell 
! subroutine will call a stop
         IF((DES_POS_NEW(L,1) < ZERO .OR. DES_POS_NEW(L,1) > XLENGTH) .AND.&
         .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! A new or exiting particle may exist in ghost cells
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000) 
            ENDIF         
            WRITE(*,'(5X,A,I10)') &
               'X position outside domain for particle ', L
            WRITE(*,'(7X,A,3(ES17.9))')&
               'particle pos = ', DES_POS_NEW(L,:)
            WRITE(*,'(7X,A,3(ES17.9))')&
               'particle vel = ', DES_VEL_NEW(L,:)
         ENDIF 

         IF((DES_POS_NEW(L,2) < ZERO .OR. DES_POS_NEW(L,2) > YLENGTH) .AND.&
         .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! A new or exiting particle may exist in ghost cells
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000) 
            ENDIF         
            WRITE(*,'(5X,A,I10)') &
               'Y position outside domain for particle=: ', L
            WRITE(*,'(7X,A,3(ES17.9))')&
               'particle pos = ', DES_POS_NEW(L,:)
            WRITE(*,'(7X,A,3(ES17.9))')&
               'particle vel = ', DES_VEL_NEW(L,:)
         ENDIF 

         IF (DIMN .EQ. 3) THEN
            IF((DES_POS_NEW(L,3) < ZERO .OR. &
                DES_POS_NEW(L,3) > ZLENGTH) .AND.&
            .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
! A new or exiting particle may exist in ghost cells
               IF (.NOT.DES_LOC_DEBUG) THEN
                  DES_LOC_DEBUG = .TRUE.
                  WRITE(*,1000) 
               ENDIF         
               WRITE(*,'(5X,A,I10)') &
                  'Z position outside domain for particle ', L
               WRITE(*,'(7X,A,3(ES17.9))')&
                  'particle pos = ', DES_POS_NEW(L,:)
               WRITE(*,'(7X,A,3(ES17.9))')&
                  'particle vel = ', DES_VEL_NEW(L,:)
            ENDIF
         ENDIF 


! Reset total contact force and torque      
         FC(L,:) = ZERO
         TOW(L,:) = ZERO


         IF (DES_LOC_DEBUG) WRITE(*,1001)

         PC = PC + 1
      ENDDO

 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')  

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)         

      RETURN
      END SUBROUTINE CFNEWVALUES



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: DES_PERIODIC_BC
!
!  Purpose: At this point the particles have already been advanced in
!     time according to the force balance.  Now, update the particle 
!     position/velocity according to simple periodic boundary conditions.
!     Only particles that have crossed a periodic BC will have their
!     positon modified accordingly. 
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE DES_PERIODIC_BC(L)

!-----------------------------------------------
! Modules
!-----------------------------------------------
     
      USE geometry
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! given particle ID number
      INTEGER, INTENT (IN) :: L

! local variables for system dimensions
      DOUBLE PRECISION LX, LXE, LXW, LY, LYN, LYS, LZ, LZT, LZB

! local variables for x, y, z position of the particle
      DOUBLE PRECISION XPOS, YPOS, ZPOS 

! local variables for i, j, k indices associated with particle position
      INTEGER I_INDEX, J_INDEX, K_INDEX, I, J, K      
!-----------------------------------------------
! Functions 
!-----------------------------------------------

!-----------------------------------------------

! assign temporary local variables for quick reference
      LXE = EX2
      LXW = WX1
      LX = LXE - LXW
      LYN = TY2
      LYS = BY1
      LY = LYN - LYS
      LZT = NZ2
      LZB = SZ1
      LZ = LZT - LZB

! initialize local variables 
      XPOS = DES_POS_NEW(L,1)
      YPOS = DES_POS_NEW(L,2)
      I_INDEX = PIJK(L,1)
      J_INDEX = PIJK(L,2)
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(L,3)
         K_INDEX = PIJK(L,3)
      ENDIF


      IF(DES_PERIODIC_WALLS_X) THEN
         IF(XPOS.GE.LXE) THEN
            XPOS = XPOS - LX
            I_INDEX = IMIN1
         ELSEIF(XPOS.LT.LXW) THEN
            XPOS = XPOS + LX
            I_INDEX = IMAX1
         ENDIF
      ENDIF
      IF(DES_PERIODIC_WALLS_Y) THEN
         IF(YPOS.GE.LYN) THEN
            YPOS = YPOS - LY
            J_INDEX = JMIN1
         ELSEIF(YPOS.LT.LYS) THEN
            YPOS = YPOS + LY
            J_INDEX = JMAX1
         ENDIF
      ENDIF
      IF(DIMN.EQ.3 .AND. DES_PERIODIC_WALLS_Z) THEN
         IF(ZPOS.GE.LZT) THEN
            ZPOS = ZPOS - LZ
            K_INDEX = KMIN1
         ELSEIF(ZPOS.LT.LZB) THEN
            ZPOS = ZPOS + LZ
            K_INDEX = KMAX1
         ENDIF
      ENDIF

! set particle position and index according to periodicity      
      DES_POS_NEW(L,1) = XPOS
      DES_POS_NEW(L,2) = YPOS
      PIJK(L,1) = I_INDEX
      PIJK(L,2) = J_INDEX
      IF (DIMN.EQ.3) THEN
         DES_POS_NEW(L,3) = ZPOS
         PIJK(L,3) = K_INDEX
      ENDIF

      RETURN

      END SUBROUTINE DES_PERIODIC_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: DES_LEES_EDWARDS_BC
!
!  Purpose: At this point the particles have already been advanced in
!     time according to the force balance.  Now, update the particle 
!     position/velocity according to Lees & Edwards boundary conditions.
!     Only particles that have crossed a LE BC will have their positon/
!     velocity modified accordingly. 
!
!
!  Author: Janine Galvin                              Date: 
!  Reviewer:                                          Date: 
!
!  Comments: For further details recommend Computer Simulation of
!     Liquids, by M.P. Allen and D. J. Tildesley (section 8.2, 
!     pages 246-247). 
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE DES_LEES_EDWARDS_BC(L)

!-----------------------------------------------
! Modules
!-----------------------------------------------
     
      USE geometry
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! given particle ID number
      INTEGER, INTENT (IN) :: L

! x, y, z position of the particle
      DOUBLE PRECISION XPOS, YPOS, ZPOS

! x, y, z velocity of the particle      
      DOUBLE PRECISION XVEL, YVEL, ZVEL

! local variables for system dimensions
      DOUBLE PRECISION LX, LY, LZ, LXE, LXW, LYN, LYS, LZT, LZB
      
! i, j, k indices associated with particle position
      INTEGER I_INDEX, J_INDEX, K_INDEX, I, J, K

! the value of the i, j, k index minus two (i.e., I=I-2) or 1 (I=1)
! whichever is greater.  this will prevent errors in referencing the
! particle position by indices outside the array bounds.  For example,
! it is possible that a particle is in a cell with i = 2 to move to a
! cell with i = 4, wherein XE(I-2) may be used causing the code to fail.
      INTEGER IMINUS2, JMINUS2, KMINUS2

! local variable for shear direction
      CHARACTER*4 SHEAR_DIR      

! local variable for relative velocity of shear
      DOUBLE PRECISION REL_VEL

! quantity to indicate direction of particle movement resulting from
! crossing the LE BC
!   0 = no move (i.e. particle did not cross LE BC)      
!  -1 = west/south/bottom move
!   1 = east/north/top move
      INTEGER MOVE_DIR

! the distance between the periodic boundaries corresponding to the
! direction the shear is acting. for du/dy shear this corresponds to the
! x domain length
      DOUBLE PRECISION DOMAIN_SIZE

! determined by first calculating the distance the LE boundary (cell) 
! that was originally aligned with the center cell traveled in a given
! time step.  then integer multiples of the domain size are subtracted
! from this quantity until a distance less than the domain size remains
      DOUBLE PRECISION OFFSET_DISTANCE

! logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

!-----------------------------------------------
! Functions 
!-----------------------------------------------

!-----------------------------------------------      

! Given that a particle's movement is currently restrained to be less 
! than its own radius during a single time step, the following
! manipulations only need to be conducted once in a given step.  If a
! particle could travel more than the domain length in a single time
! step then the 'IF' statements should be replaced with 'WHILE' loops.

! For example, in the case of DUDY shear the looping would need to be 
! performed until the new y position is smaller than ylength and/or 
! the distance between the periodic boundaries corresponding to the
! direction the shear is acting. For du/dy shear this corresponds to the
! x domain length greater than 0.  In such a case, looping would also be
! needed for adjusting the x position and additional measures would be 
! needed for locating the particle's i and j indices. Such calculations
! are not accounted for in the current code.

! If the particle does not cross the LE boundary nor the periodic
! boundary its indices will be correctly determined using the existing
! method found in particles_in_cell.f

! assign temporary local variables for quick reference
      LXE = EX2
      LXW = WX1
      LX = LXE - LXW
      LYN = TY2
      LYS = BY1
      LY = LYN - LYS
      LZT = NZ2
      LZB = SZ1
      LZ = LZT - LZB

      REL_VEL = DES_LE_REL_VEL
      SHEAR_DIR = TRIM(DES_LE_SHEAR_DIR)

! assign temporary local variables for manipulation/use
      XPOS = DES_POS_NEW(L,1)
      YPOS = DES_POS_NEW(L,2)
      XVEL = DES_VEL_NEW(L,1)
      YVEL = DES_VEL_NEW(L,2)
      I_INDEX = PIJK(L,1)
      J_INDEX = PIJK(L,2)      
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(L,3)
         ZVEL = DES_VEL_NEW(L,3)
         K_INDEX = PIJK(L,3)
      ENDIF 


! initialize local quantities
      OFFSET_DISTANCE = 0
      IF (REL_VEL > 0 ) THEN
         MOVE_DIR = 1   ! particle is moved east, north or up
      ELSEIF (REL_VEL < 0 ) THEN
         MOVE_DIR = -1  ! particle is moved west, south, or down
      ELSE
         MOVE_DIR = 0
      ENDIF


      IF (DIMN .EQ. 2) THEN

! 2D shear : du/dy
! ----------------------------------------               
         IF(SHEAR_DIR.EQ.'DUDY') THEN
            DOMAIN_SIZE = LX
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF
                        
            IF (YPOS >= LYN) THEN
! particle crossed north Lees & Edwards boundary                    
               YPOS = YPOS - LY 
               J_INDEX = JMIN1
               XPOS = XPOS - OFFSET_DISTANCE
               XVEL = XVEL - REL_VEL
               MOVE_DIR = -1*MOVE_DIR
            ELSEIF (YPOS < LYS) THEN
! particle crossed south Lees & Edwards boundary            
               YPOS= YPOS +  LY
               J_INDEX = JMAX1
               XPOS = XPOS + OFFSET_DISTANCE
               XVEL = XVEL + REL_VEL
            ELSE
! particle did not cross the Lees & Edwards boundary                    
               MOVE_DIR = 0
            ENDIF

            IF (XPOS >= LXE) THEN
! particle crossed east periodic boundary 
               XPOS = XPOS - DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN 
! particle did not cross the LE boundary so the value of PIJK can
! readily be assigned as a particle is constrained to move less than 
! its own radius during a single time step
                  I_INDEX = IMIN1
               ELSE
! particle crossed the LE boundary and depending on the shear rate and
! particle position it is possible for the particle to have been moved
! more than its own radius and outside the domain.  the particle was
! moved east (move_dir = 1). 
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
! in this case the particle must have been on the east side of the
! domain for it to have crossed the east periodic boundary and the move
! will shift it to the west side of the domain. so to minimize the
! search start at the west most domain index and increment up toward the
! index corresponding to the old x position
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF

            ELSEIF (XPOS < LXW) THEN
! particle crossed west periodic boundary 
               XPOS = XPOS + DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN
! particle did not cross the LE boundary
                  I_INDEX = IMAX1
               ELSE
! particle crossed the LE boundary. the particle was moved left 
! (move_dir=-1).  
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
! in this case the particle must have been on the west side of the
! domain for it to have crossed the west periodic boundary and the move
! will shift it to the east side of the domain
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF

            ELSEIF (MOVE_DIR .NE. 0) THEN
! particle did not cross the periodic boundary but did cross the LE
! boundary (i.e. particle was moved but not across domain). identify the
! index by first searching if particle is in the 2 outer most cells then
! searching the the cells neighboring the old position

               I = I_INDEX   ! for shorthand/quick reference
               IMINUS2 = I - 2
               IF (IMINUS2 < 1) IMINUS2 = 1

               IF(XPOS >= XE(1) .AND. XPOS < XE(IMIN1)) THEN
                  I_INDEX = IMIN1
               ELSEIF (XPOS >= XE(IMAX1-1) .AND. XPOS < XE(IMAX1)) THEN
                  I_INDEX = IMAX1
               ELSEIF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN 
                  I_INDEX = I
               ELSEIF(XPOS >= XE(I) .AND. XPOS < XE(I+1)) THEN 
                  I_INDEX = I+1
               ELSEIF(XPOS >= XE(IMINUS2) .AND. XPOS < XE(I-1)) THEN
                     I_INDEX = I-1
               ELSEIF (MOVE_DIR .EQ. -1) THEN
! particle was moved west
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
! in this case the particle must have been on the east side of the domain 
! and the move will shift it to the west of the domain
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ELSEIF (MOVE_DIR .EQ. 1) THEN
! particle was moved east
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
! in this case the particle must have been on the LHS of the domain and
! the move will shift it to the RHS of the domain
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF  ! end final block to identify index

            ENDIF   ! the particle did not cross a LE or periodic boundary 


! 2D shear : dv/dx
! see 2D DUDY shear section for more details on code
! ---------------------------------------- 
         ELSEIF(SHEAR_DIR.EQ.'DVDX') THEN
            DOMAIN_SIZE = LY            
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF
          
            IF (XPOS >= LXE) THEN
! particle crossed east Lees & Edwards boundary
               XPOS = XPOS - LX 
               I_INDEX = IMIN1
               YPOS = YPOS - OFFSET_DISTANCE
               YVEL = YVEL - REL_VEL
               MOVE_DIR = -1*MOVE_DIR
            ELSEIF (XPOS < LXW) THEN
! particle crossed west Lees & Edwards boundary
               XPOS = XPOS + LX 
               I_INDEX = IMAX1
               YPOS = YPOS + OFFSET_DISTANCE
               YVEL = YVEL + REL_VEL
            ENDIF

            IF (YPOS >= LYN) THEN
! particle crossed north periodic boundary
               YPOS = YPOS - DOMAIN_SIZE
               IF (MOVE_DIR .EQ. 0) THEN
! particle did not cross the LE boundary
                  J_INDEX = JMIN1
               ELSE
! the particle crossed the LE boundary
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMIN1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMIN1,YPOS,&
                        YN,'y','j')
                  ENDIF
               ENDIF

            ELSEIF (YPOS < LYS) THEN
! particle crossed south periodic boundary
               YPOS = YPOS + DOMAIN_SIZE
               IF (MOVE_DIR .EQ. 0) THEN
! particle did not cross the LE boundary
                  J_INDEX = JMAX1
               ELSE
! particle crossed the LE boundary
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMAX1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMAX1,YPOS,&
                        YN,'y','j')
                  ENDIF
               ENDIF

            ELSEIF (MOVE_DIR .NE. 0) THEN
! particle did not cross the periodic boundary but did cross the LE
! boundary
               J = J_INDEX   ! for shorthand/quick reference
               JMINUS2 = J - 2
               IF (JMINUS2 < 1) JMINUS2 = 1

               IF(YPOS >= YN(1) .AND. YPOS < YN(JMIN1)) THEN
                  J_INDEX = JMIN1
               ELSEIF(YPOS >= YN(JMAX1-1) .AND. YPOS < YN(JMAX1)) THEN
                  J_INDEX = JMAX1
               ELSEIF(YPOS >= YN(J-1) .AND. YPOS < YN(J)) THEN
                  J_INDEX = J
               ELSEIF(YPOS >= YN(J) .AND. YPOS < YN(J+1)) THEN
                  J_INDEX = J+1
               ELSEIF(YPOS >= YN(JMINUS2) .AND. YPOS < YN(J-1)) THEN
                  J_INDEX = J-1
               ELSEIF (MOVE_DIR .EQ. -1) THEN
! particle was moved south
                  IF(DABS(OFFSET_DISTANCE) >= 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMIN1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMIN1,YPOS,&
                        YN,'y','j')
                  ENDIF
               ELSEIF(MOVE_DIR .EQ. 1) THEN
! the particle was moved north
                  IF (DABS(OFFSET_DISTANCE) >= 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMAX1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMAX1,YPOS,&
                        YN,'y','j')
                  ENDIF                 
               ENDIF   ! end final block to identify index

            ENDIF   ! the particle did not cross a LE or periodic boundary

         ENDIF   ! endif des_le_shear_dir == dudy or dvdx
      ENDIF  ! if dimn == 2


      IF (DIMN .EQ. 3) THEN
! 3D shear : du/dy
! ---------------------------------------- 
         IF(TRIM(DES_LE_SHEAR_DIR).EQ.'DUDY') THEN
            DOMAIN_SIZE = LX
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
               XPOS = XPOS - OFFSET_DISTANCE
               XVEL = XVEL - REL_VEL
               MOVE_DIR = -1*MOVE_DIR
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
               XPOS = XPOS + OFFSET_DISTANCE
               XVEL = XVEL + REL_VEL
            ELSE
               MOVE_DIR = 0
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN 
                  I_INDEX = IMIN1
               ELSE
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN
                  I_INDEX = IMAX1
               ELSE
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF
            ELSEIF (MOVE_DIR .NE. 0) THEN
               I = I_INDEX   ! for shorthand/quick reference
               IMINUS2 = I - 2
               IF (IMINUS2 < 1) IMINUS2 = 1

               IF(XPOS >= XE(1) .AND. XPOS < XE(IMIN1)) THEN
                  I_INDEX = IMIN1
               ELSEIF (XPOS >= XE(IMAX1-1) .AND. XPOS < XE(IMAX1)) THEN
                  I_INDEX = IMAX1
               ELSEIF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN 
                  I_INDEX = I
               ELSEIF(XPOS >= XE(I) .AND. XPOS < XE(I+1)) THEN 
                  I_INDEX = I+1
               ELSEIF(XPOS >= XE(IMINUS2) .AND. XPOS < XE(I-1)) THEN
                     I_INDEX = I-1
               ELSEIF (MOVE_DIR .EQ. -1) THEN
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ELSEIF (MOVE_DIR .EQ. 1) THEN
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF  ! end final block to identify i index

            ENDIF   ! the particle did not cross a LE or periodic boundary 

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
            ENDIF

! 3D shear : du/dz
! ----------------------------------------             
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DUDZ') THEN
            DOMAIN_SIZE = LX
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
               XPOS = XPOS - OFFSET_DISTANCE
               XVEL = XVEL - REL_VEL
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
               XPOS = XPOS + OFFSET_DISTANCE
               XVEL = XVEL + REL_VEL
            ENDIF

! work out code to find i index 
            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
            ENDIF
            
            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
            ENDIF

! 3D shear : dv/dx
! ---------------------------------------- 
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DVDX') THEN
            DOMAIN_SIZE = LY
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
               I_INDEX = IMIN1
               YPOS = YPOS - OFFSET_DISTANCE
               YVEL = YVEL - REL_VEL
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
               I_INDEX = IMAX1
               YPOS = YPOS + OFFSET_DISTANCE
               YVEL = YVEL + REL_VEL
            ENDIF

! work out code to find j index            
            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
            ENDIF

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
            ENDIF

! 3D shear : dv/dz
! ----------------------------------------             
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DVDZ') THEN
            DOMAIN_SIZE = LY
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
               YPOS = YPOS - OFFSET_DISTANCE
               YVEL = YVEL - REL_VEL
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
               YPOS = YPOS + OFFSET_DISTANCE
               YVEL = YVEL + REL_VEL
            ENDIF

! work out code to find j index                   
            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
               I_INDEX = IMIN1
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
               I_INDEX = IMAX1
            ENDIF


! 3D shear : dw/dx
! ---------------------------------------- 
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DWDX') THEN
            DOMAIN_SIZE = LZ
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - (EX2 - LXW)
               I_INDEX = IMIN1
               ZPOS = ZPOS - OFFSET_DISTANCE
               ZVEL = ZVEL - REL_VEL
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + (EX2 - LXW)
               I_INDEX = IMAX1
               ZPOS = ZPOS + OFFSET_DISTANCE
               ZVEL = ZVEL + REL_VEL
            ENDIF

! work out code to find k index            
            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
            ENDIF

            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
            ENDIF


! 3D shear : dw/dy
! ----------------------------------------             
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DWDY') THEN
            DOMAIN_SIZE = LZ
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
               ZPOS = ZPOS - OFFSET_DISTANCE
               ZVEL = ZVEL - REL_VEL
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
               ZPOS = ZPOS + OFFSET_DISTANCE
               ZVEL = ZVEL + REL_VEL
            ENDIF

! work out code to find k index
            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
               I_INDEX = IMIN1
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
               I_INDEX = IMAX1
            ENDIF


         ENDIF   ! endif des_le_shear_dir == dudy, dudz, dvdx, dvdz, dwdx or dwdz

      ENDIF   ! endif dimn == 3

! set particle position and index according to periodicity/lees & edwards BC
      DES_POS_NEW(L,1) = XPOS  
      DES_POS_NEW(L,2) = YPOS
      DES_VEL_NEW(L,1) = XVEL
      DES_VEL_NEW(L,2) = YVEL
      PIJK(L,1) = I_INDEX
      PIJK(L,2) = J_INDEX
      IF (DIMN .EQ. 3) THEN
         DES_POS_NEW(L,3) = ZPOS  
         DES_VEL_NEW(L,3) = ZVEL
         PIJK(L,3) = K_INDEX
      ENDIF

      RETURN

 1002 FORMAT(/1X,70('*')//,' From: DES_DES_LEES_EDWARDS_BC',/,&
         ' Message: Calculation for the ', A, ' index associated with ',&
         'the particles new',/10X, A, '-position= ', ES12.5,&
         ' where the old ',A, '-position= ',/10X, ES12.5 ' and old ',&
         A, ' index= ', I5, /1X,70('*')/)



      END SUBROUTINE DES_LEES_EDWARDS_BC
      
