!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: GRID_BASED_NEIGHBOR_SEARCH                             C
!  Purpose: Cell linked search 
!                                                                      C
!  Author: Rahul Garg                                 Date: 01-Aug-07  C
!  Reviewer: Sreekanth Pannala                        Date: 24-OCT-08  C
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      
      SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH
      
      USE discretelement

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! loop counter, particle no id. (current target particle)
      INTEGER LL

! index to track accounted for particles
      INTEGER PC

!-----------------------------------------------
! Functions
!-----------------------------------------------         

!-----------------------------------------------  


! Loop only over particles in the system
!   cycle if there is no particle in the LL seat

      PC = 1
      DO LL = 1, MAX_PIS
         IF(PC .GT. PIS) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE

         IF(DES_LE_BC) THEN
            CALL GRIDBASED_LEBC_NEIGHBOR_CELL_SEARCH(LL)
         ELSE
            CALL GRIDBASED_STANDARD_NEIGHBOR_CELL_SEARCH(LL)
         ENDIF

         PC = PC + 1

      ENDDO  ! Particles in system loop

 
      RETURN
      END SUBROUTINE GRID_BASED_NEIGHBOR_SEARCH




!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
!  Module name: GRIDBASED_STANDARD_NEIGHBOR_CELL_SEARCH
!  Purpose: Conduct a search of each of the surrounding grid cells for 
!     potential neighbors, if the indices of the current target cell
!     fall beyond a periodic boundary then adjust the target particle
!     position and shift the indices of the cell to be searched for
!     neighbors accordingly.  
!
!     Note if a boundary is periodic then the boundary cannot have a 
!     DEM inlet/outlet (that is, DEM_MI/MO cannot also be true) or the
!     code will have already exited.      
!      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC      
      
      SUBROUTINE GRIDBASED_STANDARD_NEIGHBOR_CELL_SEARCH(LL)
      
      USE discretelement
      USE des_bc

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! given loop counter, particle no id. (current target particle)
      INTEGER LL

! local variables for system dimensions
      DOUBLE PRECISION LX, LY, LZ

! loop indices for x, y, z coordinates      
      INTEGER EXTRA_CELLS
      INTEGER I_INDEX, J_INDEX, K_INDEX
      INTEGER II, IP, IM, JJ, JP, JM, KK, KP, KM

! min and max fluid cell index for grid based search mesh
!  min fluid cell index no. (i.e. 2)      
      INTEGER DESGS_IMIN1, DESGS_JMIN1, DESGS_KMIN1
!  max fluid cell index no. (max+1)
      INTEGER DESGS_IMAX1, DESGS_JMAX1, DESGS_KMAX1
!  min cell index no. (i.e. 1) (ghost cell)
      INTEGER DESGS_IMIN2,DESGS_JMIN2, DESGS_KMIN2

! local variables for x, y, z position of the particle
      DOUBLE PRECISION XPOS, YPOS, ZPOS 

! temporary variable to store particle LL position (target particle)
! which may be modified from actual particle position due to BC
      DOUBLE PRECISION POS_LL(DIMN) 

!-----------------------------------------------
! Functions
!-----------------------------------------------         

!-----------------------------------------------  
! initialize 
      POS_LL(:) = DES_POS_NEW(LL,:)

! assign temporary local variables for quick reference
      LX = EX2 - WX1 !=XLENGTH =XE(IMAX1) - XE(1)
      LY = TY2 - BY1 !=YLENGTH =YN(JMAX1) - YN(1)
      LZ = NZ2 - SZ1 !=ZLENGTH =ZT(KMAX1) - ZT(1)

      DESGS_IMAX1 = DESGRIDSEARCH_IMAX+1
      DESGS_IMIN1 = 2
      DESGS_IMIN2 = 1
      DESGS_JMAX1 = DESGRIDSEARCH_JMAX+1
      DESGS_JMIN1 = 2
      DESGS_JMIN2 = 1
      IF (DIMN .EQ. 2) THEN
         DESGS_KMAX1 = DESGRIDSEARCH_KMAX
         DESGS_KMIN1 = DESGRIDSEARCH_KMAX
         DESGS_KMIN2 = DESGRIDSEARCH_KMAX
      ELSE
         DESGS_KMAX1 = DESGRIDSEARCH_KMAX+1
         DESGS_KMIN1 = 2
         DESGS_KMIN2 = 1
      ENDIF

      XPOS = DES_POS_NEW(LL,1)
      YPOS = DES_POS_NEW(LL,2) 
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(LL,3)
      ENDIF      

! Increment indices by +/- indicated quantity (extra_cells). So in this
! case search only the immediate surrounding cells
      EXTRA_CELLS = 1
      I_INDEX  = DESGRIDSEARCH_PIJK(LL,1)
      IP = I_INDEX + EXTRA_CELLS
      IM = I_INDEX - EXTRA_CELLS
      J_INDEX  = DESGRIDSEARCH_PIJK(LL,2)
      JP = J_INDEX + EXTRA_CELLS
      JM = J_INDEX - EXTRA_CELLS
      K_INDEX  = DESGRIDSEARCH_PIJK(LL,3)      
      KP = K_INDEX 
      KM = K_INDEX
      IF(DIMN.EQ.3)THEN 
         KP = K_INDEX + EXTRA_CELLS
         KM = K_INDEX - EXTRA_CELLS
      ENDIF

! Loop over indices for potential neighbors
      KLOOP : DO KK = KM, KP
         JLOOP : DO JJ = JM, JP
            ILOOP : DO II = IM, IP

! Check i index
! --------------------               
               IF (.NOT. DES_MI_X .AND. .NOT. DES_MO_X) THEN
                  IF(II.GT.DESGS_IMAX1) THEN 
! If particle crossed over a periodic boundary, correct the index for
! the neighbor cells and shift the target position accordingly
                     IF(DES_PERIODIC_WALLS_X) THEN
! Generally I becomes DESGS_IMIN1 (the min fluid cell no.) but rather
! than simply assigning the index as such, this calculation allows
! greater flexibility in the number of elements/cells to be searched
                        I_INDEX = II - DESGRIDSEARCH_IMAX  
                        POS_LL(1) = XPOS - LX
                     ELSE
! No need to search since the index would be set such that 
! an earlier/upcoming search would be repeated
                        CYCLE ILOOP
                     ENDIF
                  ELSEIF(II.LT.DESGS_IMIN1) THEN 
                     IF(DES_PERIODIC_WALLS_X) THEN 
! Generally this becomes DESGS_IMAX1 the max fluid cell no. 
                        I_INDEX = II + DESGRIDSEARCH_IMAX 
                        POS_LL(1) = XPOS + LX
                     ELSE
                        CYCLE ILOOP
                     ENDIF
                  ELSE
                     I_INDEX = II
                     POS_LL(1) = XPOS
                  ENDIF
               ELSE
! If DEM inlet or outlet, then a particle may exist in a ghost cell.
                  IF(II.GT.DESGS_IMAX2 .OR. II.LT.DESGS_IMIN2) &
                     CYCLE ILOOP   
               ENDIF

! Check j index
! --------------------               
               IF(.NOT.DES_MI_Y .AND. .NOT.DES_MO_Y) THEN
                  IF(JJ.GT.DESGS_JMAX1) THEN
                     IF(DES_PERIODIC_WALLS_Y) THEN 
                        J_INDEX = JJ - DESGRIDSEARCH_JMAX  ! DESGS_JMIN1
                        POS_LL(2) = YPOS - LY
                     ELSE
                        CYCLE JLOOP
                     ENDIF
                   ELSEIF(JJ.LT.DESGS_JMIN1) THEN
                     IF(DES_PERIODIC_WALLS_Y) THEN 
                        J_INDEX = JJ + DESGRIDSEARCH_JMAX   ! DESGS_JMAX1    
                        POS_LL(2) = YPOS + LY
                     ELSE
                        CYCLE JLOOP
                     ENDIF
                  ELSE
                     J_INDEX = JJ
                     POS_LL(2) = YPOS
                  ENDIF
               ELSE
                  IF(JJ.GT.DESGS_JMAX2 .OR. JJ.LT.DESGS_JMIN2) &
                     CYCLE JLOOP
               ENDIF

! Check k index
! --------------------               
               IF(DIMN.EQ.3) THEN 
                  IF(.NOT.DES_MI_Z .AND. .NOT.DES_MO_Z)THEN
                     IF(KK.GT.DESGS_KMAX1) THEN 
                        IF(DES_PERIODIC_WALLS_Z) THEN
                           K_INDEX = KK - DESGRIDSEARCH_KMAX   !DESGS_KMIN1
                           POS_LL(3) = ZPOS - LZ
                        ELSE
                           CYCLE KLOOP
                        ENDIF
                     ELSEIF(KK.LT.DESGS_KMIN1) THEN 
                        IF(DES_PERIODIC_WALLS_Z) THEN
                           K_INDEX = KK + DESGRIDSEARCH_KMAX   !DESGS_KMAX1
                           POS_LL(3) = ZPOS + LZ
                        ELSE
                           CYCLE KLOOP
                        ENDIF
                     ELSE
                        K_INDEX = KK
                        POS_LL(3) = ZPOS
                     ENDIF
                  ELSE
                     IF(KK.GT.DESGS_KMAX2 .OR. KK.LT.DESGS_KMIN2) &
                        CYCLE KLOOP
                  ENDIF
               ELSE
                  K_INDEX = KK
               ENDIF

! Search indicated IJK cell for neighbors with particle LL (the
! current target particle).  LL's position has already been shifted 
! according to the boundary rules (e.g., wrap around for periodic)
               CALL GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL(LL, &
                  POS_LL,I_INDEX,J_INDEX,K_INDEX)
               

            ENDDO ILOOP   ! II cell loop
         ENDDO JLOOP   ! JJ cell loop
      ENDDO KLOOP   ! KK cell loop


      RETURN
      END SUBROUTINE GRIDBASED_STANDARD_NEIGHBOR_CELL_SEARCH



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!
!  Module name: GRIDBASED_LEBC_NEIGHBOR_CELL_SEARCH
!  Purpose: Conduct a search of each of the surrounding grid cells for 
!     potential neighbors, if the indices of the current target cell
!     fall beyond a periodic boundary then adjust the target particle
!     position and shift the indices of the cell to be searched for
!     neighbors accordingly.  
!
!     Note if a boundary is periodic then the boundary cannot have a 
!     DEM inlet/outlet (that is, DEM_MI/MO cannot also be true) or the
!     code will have already exited.      
!      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC      
      
      SUBROUTINE GRIDBASED_LEBC_NEIGHBOR_CELL_SEARCH(LL)

      USE compar       
      USE discretelement

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! given loop counter, particle no id. (current target particle)
      INTEGER LL

! loop indices for x, y, z coordinates      
      INTEGER EXTRA_CELLS
      INTEGER I_INDEX, J_INDEX, K_INDEX
      INTEGER II, IP, IM, JJ, JP, JM, KK, KP, KM

! min and max fluid cell index for grid based search mesh
!  min fluid cell index no. (i.e. 2)      
      INTEGER DESGS_IMIN1, DESGS_JMIN1, DESGS_KMIN1
!  max fluid cell index no. (max+1)
      INTEGER DESGS_IMAX1, DESGS_JMAX1, DESGS_KMAX1

! temporary variable to store particle LL position (target particle)
! which may be modified from actual particle position due to BC
      DOUBLE PRECISION XPOS, YPOS, ZPOS, XPOS_TMP, YPOS_TMP, &
                       ZPOS_TMP, POS_LL(DIMN) 

! dimension/size of grid cells      
      DOUBLE PRECISION SIZEDX, SIZEDY, SIZEDZ

! local variables for system dimensions
      DOUBLE PRECISION LX, LY, LZ, LXE, LXW, LYN, LYS, LZT, LZB

! local variable for shear direction
      CHARACTER*4 SHEAR_DIR      

! local variable for relative velocity of shear
      DOUBLE PRECISION REL_VEL

! the distance between the periodic boundaries corresponding to the
! direction the shear is acting. for du/dy shear this corresponds to the
! x domain length
      DOUBLE PRECISION DOMAIN_SIZE

! determined by first calcauling the distance the LE boundary (cell) 
! that was originally aligned with the center cell traveled in a given
! time step.  then integer multiples of the domain_size are subtracted
! from this quantity until a distance less than the domain size remains
      DOUBLE PRECISION OFFSET_DISTANCE    

!-----------------------------------------------
! Functions
!-----------------------------------------------         

!-----------------------------------------------  
! initialize
      POS_LL(:) = DES_POS_NEW(LL,:)


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

      DESGS_IMAX1 = DESGRIDSEARCH_IMAX+1
      DESGS_IMIN1 = 2
      DESGS_JMAX1 = DESGRIDSEARCH_JMAX+1
      DESGS_JMIN1 = 2
      IF (DIMN .EQ. 2) THEN
         DESGS_KMAX1 = DESGRIDSEARCH_KMAX
         DESGS_KMIN1 = DESGRIDSEARCH_KMAX
      ELSE
         DESGS_KMAX1 = DESGRIDSEARCH_KMAX+1
         DESGS_KMIN1 = 2
      ENDIF

      REL_VEL = DES_LE_REL_VEL
      SHEAR_DIR = TRIM(DES_LE_SHEAR_DIR)

      XPOS = DES_POS_NEW(LL,1)
      YPOS = DES_POS_NEW(LL,2) 
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(LL,3)
      ENDIF      

! Increment indices by +/-1 (search only immediate surrounding cells)
! The number of 'extra_cells' cannot be in excess of the maximum number
! of search cells (imax, jmax, kmax) 
      EXTRA_CELLS = 1
      I_INDEX = DESGRIDSEARCH_PIJK(LL,1)
      J_INDEX = DESGRIDSEARCH_PIJK(LL,2)
      K_INDEX = DESGRIDSEARCH_PIJK(LL,3) 

! Size of the mesh
      SIZEDX = LX/DESGRIDSEARCH_IMAX
      SIZEDY = LY/DESGRIDSEARCH_JMAX
      IF (DIMN .EQ. 2) THEN
         SIZEDZ = ONE
      ELSE
         SIZEDZ = LZ/DESGRIDSEARCH_KMAX
      ENDIF

      IF (DIMN .EQ. 2) THEN

! 2D shear : du/dy
! ----------------------------------------
         IF(TRIM(SHEAR_DIR).EQ.'DUDY') THEN
            DOMAIN_SIZE = LX            
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE*&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

! Determine j indices to loop over during neighbor search
            JP = J_INDEX + EXTRA_CELLS
            JM = J_INDEX - EXTRA_CELLS

! Loop over j indices                
            JLOOP1 : DO JJ = JM, JP

               IF(JJ.GT.DESGS_JMAX1) THEN
                  J_INDEX = JJ - DESGRIDSEARCH_JMAX
                  POS_LL(2) = YPOS - LY
                  XPOS_TMP = XPOS - OFFSET_DISTANCE
               ELSEIF(JJ.LT.DESGS_JMIN1) THEN
                  J_INDEX = JJ + DESGRIDSEARCH_JMAX
                  POS_LL(2) = YPOS + LY
                  XPOS_TMP = XPOS + OFFSET_DISTANCE
               ELSE
                  J_INDEX = JJ
                  POS_LL(2) = YPOS
                  XPOS_TMP = XPOS
               ENDIF

! Shift particle x-position to within (can be above or below) x-domain
               XPOS_TMP = XPOS_TMP - FLOOR(XPOS_TMP/LX)*LX

! Find i index associated with this particle x-position 
               IF (XPOS_TMP < LXW .OR. XPOS_TMP >= LXE) THEN
                  WRITE(*,1001) 'X', LL
                  CALL MFIX_EXIT(myPE)                  
               ELSE         
                  I_INDEX = INT(XPOS_TMP/SIZEDX)+2
               ENDIF

! Calculate i indices to loop over during neighbor search
               IP=I_INDEX+EXTRA_CELLS
               IM=I_INDEX-EXTRA_CELLS    

! Loop over i indices                
               ILOOP1 : DO II = IM, IP

! Store current loop index to new variable for manipulation
                  IF(II.GT.DESGS_IMAX1) THEN
                     I_INDEX = II - DESGRIDSEARCH_IMAX
                     POS_LL(1) = XPOS_TMP - LX
                  ELSEIF(II.LT.DESGS_IMIN1) THEN
                     I_INDEX = II + DESGRIDSEARCH_IMAX
                     POS_LL(1) = XPOS_TMP + LX
                  ELSE
                     I_INDEX = II
                     POS_LL(1) = XPOS_TMP
                  ENDIF

! Error checking : the ijk cell being searched for neighbors must be a
! 'fluid cell' (not a boundary or ghost cell)

                  IF (I_INDEX .GT. DESGS_IMAX1 .OR. &
                      I_INDEX .LT. DESGS_IMIN1) THEN
                     WRITE(*,1002) 'I', I_INDEX, DESGS_IMIN1, &
                        DESGS_IMAX1, LL
                     CALL MFIX_EXIT(myPE)         
                  ENDIF
                  IF (J_INDEX .GT. DESGS_JMAX1 .OR. &
                      J_INDEX .LT. DESGS_JMIN1) THEN
                     WRITE(*,1002) 'J', J_INDEX, DESGS_JMIN1, &
                        DESGS_JMAX1, LL
                     CALL MFIX_EXIT(myPE)         
                  ENDIF

! Search indicated IJK cell for neighbors with particle LL (the
! current target particle).  LL's position has already been shifted 
! according to the boundary rules (e.g., wrap around for periodic)                  
                  CALL GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL(LL, &
                     POS_LL,I_INDEX,J_INDEX,K_INDEX)

               ENDDO ILOOP1
            ENDDO JLOOP1

! 2D shear : dv/dx
! ----------------------------------------            
         ELSEIF(TRIM(SHEAR_DIR).EQ.'DVDX') THEN
            DOMAIN_SIZE = LY            
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE*&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

! Determine i indices to loop over during neighbor search
            IP = I_INDEX + EXTRA_CELLS
            IM = I_INDEX - EXTRA_CELLS

! Loop over I indices                
            ILOOP2 : DO II = IM, IP

               IF(II.GT.DESGS_IMAX1) THEN
                  I_INDEX = II - DESGRIDSEARCH_IMAX
                  POS_LL(1) = XPOS - LX
                  YPOS_TMP = YPOS - OFFSET_DISTANCE
               ELSEIF(II.LT.DESGS_IMIN1) THEN
                  I_INDEX = II + DESGRIDSEARCH_IMAX
                  POS_LL(1) = XPOS + LX
                  YPOS_TMP = YPOS + OFFSET_DISTANCE
               ELSE
                  I_INDEX = II
                  POS_LL(1) = XPOS
                  YPOS_TMP = YPOS
               ENDIF

! Shift particle y-position to within (can be above or below) y-domain
               YPOS_TMP = YPOS_TMP - FLOOR(YPOS_TMP/LY)*LY

! Find j index associated with this particle x-position 
               IF (YPOS_TMP < LYS .OR. YPOS_TMP >= LYN) THEN
                  WRITE(*,1001) 'Y', LL
                  CALL MFIX_EXIT(myPE)                  
               ELSE         
                  J_INDEX = INT(YPOS_TMP/SIZEDY)+2
               ENDIF

! Calculate j indices to loop over during neighbor search
               JP=J_INDEX+EXTRA_CELLS
               JM=J_INDEX-EXTRA_CELLS    

! Loop over j indices                
               JLOOP2 : DO JJ = JM, JP

! Store current loop index to new variable for manipulation
                  IF(JJ.GT.DESGS_JMAX1) THEN
                     J_INDEX = JJ - DESGRIDSEARCH_JMAX
                     POS_LL(2) = YPOS_TMP - LY
                  ELSEIF(JJ.LT.DESGS_JMIN1) THEN
                     J_INDEX = JJ + DESGRIDSEARCH_JMAX
                     POS_LL(2) = YPOS_TMP + LY
                  ELSE
                     J_INDEX = JJ
                     POS_LL(2) = YPOS_TMP
                  ENDIF

! Error checking : the ijk cell being searched for neighbors must be a
! fluid cell
                  IF (I_INDEX .GT. DESGS_IMAX1 .OR. &
                      I_INDEX .LT. DESGS_IMIN1) THEN
                     WRITE(*,1002) 'I', I_INDEX, DESGS_IMIN1, &
                        DESGS_IMAX1, LL
                     CALL MFIX_EXIT(myPE)         
                  ENDIF
                  IF (J_INDEX .GT. DESGS_JMAX1 .OR. &
                      J_INDEX .LT. DESGS_JMIN1) THEN
                     WRITE(*,1002) 'J', J_INDEX, DESGS_JMIN1, &
                        DESGS_JMAX1, LL
                     CALL MFIX_EXIT(myPE)         
                  ENDIF

! Search indicated IJK cell for neighbors with particle LL (the
! current target particle).  LL's position has already been shifted 
! according to the boundary rules (e.g., wrap around for periodic)                  
                  CALL GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL(LL, &
                     POS_LL,I_INDEX,J_INDEX,K_INDEX)

               ENDDO JLOOP2
            ENDDO ILOOP2
         
         ENDIF   ! endif des_le_shear_dir == dudy or dvdx
      ENDIF   ! if (des_dimn == 2)


      IF (DIMN.EQ.3) THEN  

! 3D shear : du/dy
! ----------------------------------------
         IF(TRIM(SHEAR_DIR).EQ.'DUDY') THEN
            DOMAIN_SIZE = LX            
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE*&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

! Determine k indices to loop over during neighbor search
            KP = K_INDEX + EXTRA_CELLS
            KM = K_INDEX - EXTRA_CELLS

! Loop over k indices
            KLOOP3 : DO KK = KM, KP
               IF(KK.GT.DESGS_KMAX1) THEN
                  K_INDEX = KK - DESGRIDSEARCH_KMAX
                  POS_LL(3) = ZPOS - LZ
               ELSEIF(KK.LT.DESGS_KMIN1) THEN
                  K_INDEX = KK + DESGRIDSEARCH_KMAX
                  POS_LL(3) = ZPOS + LZ
               ELSE
                  K_INDEX = KK
                  POS_LL(3) = ZPOS
               ENDIF

! Determine j indices to loop over during neighbor search
               JP = J_INDEX + EXTRA_CELLS
               JM = J_INDEX - EXTRA_CELLS

! Loop over j indices                
               JLOOP3 : DO JJ = JM, JP

                  IF(JJ.GT.DESGS_JMAX1) THEN
                     J_INDEX = JJ - DESGRIDSEARCH_JMAX
                     POS_LL(2) = YPOS - LY
                     XPOS_TMP = XPOS - OFFSET_DISTANCE
                  ELSEIF(JJ.LT.DESGS_JMIN1) THEN
                     J_INDEX = JJ + DESGRIDSEARCH_JMAX
                     POS_LL(2) = YPOS + LY
                     XPOS_TMP = XPOS + OFFSET_DISTANCE
                  ELSE
                     J_INDEX = JJ
                     POS_LL(2) = YPOS
                     XPOS_TMP = XPOS
                  ENDIF

! Shift particle x-position to within (can be above or below) x-domain
                  XPOS_TMP = XPOS_TMP - FLOOR(XPOS_TMP/LX)*LX

! Find i index associated with this particle x-position 
                  IF (XPOS_TMP < LXW .OR. XPOS_TMP >= LXE) THEN
                     WRITE(*,1001) 'X', LL
                     CALL MFIX_EXIT(myPE)                  
                  ELSE         
                     I_INDEX = INT(XPOS_TMP/SIZEDX)+2
                  ENDIF

! Calculate i indices to loop over during neighbor search
                  IP=I_INDEX+EXTRA_CELLS
                  IM=I_INDEX-EXTRA_CELLS    

! Loop over i indices                
                  ILOOP3 : DO II = IM, IP

! Store current loop index to new variable for manipulation
                     IF(II.GT.DESGS_IMAX1) THEN
                        I_INDEX = II - DESGRIDSEARCH_IMAX
                        POS_LL(1) = XPOS_TMP - LX
                     ELSEIF(II.LT.DESGS_IMIN1) THEN
                        I_INDEX = II + DESGRIDSEARCH_IMAX
                        POS_LL(1) = XPOS_TMP + LX
                     ELSE
                        I_INDEX = II
                        POS_LL(1) = XPOS_TMP
                     ENDIF

! Error checking : the ijk cell being searched for neighbors must be a
! 'fluid cell' (not a boundary or ghost cell)

                     IF (I_INDEX .GT. DESGS_IMAX1 .OR. &
                         I_INDEX .LT. DESGS_IMIN1) THEN
                        WRITE(*,1002) 'I', I_INDEX, DESGS_IMIN1, &
                           DESGS_IMAX1, LL
                        CALL MFIX_EXIT(myPE)         
                     ENDIF
                     IF (J_INDEX .GT. DESGS_JMAX1 .OR. &
                         J_INDEX .LT. DESGS_JMIN1) THEN
                        WRITE(*,1002) 'J', J_INDEX, DESGS_JMIN1, &
                           DESGS_JMAX1, LL
                        CALL MFIX_EXIT(myPE)         
                     ENDIF
                     IF (K_INDEX .GT. DESGS_KMAX1 .OR. &
                         K_INDEX .LT. DESGS_KMIN1) THEN
                        WRITE(*,1002) 'K', K_INDEX, DESGS_KMIN1, &
                           DESGS_KMAX1, LL
                        CALL MFIX_EXIT(myPE)         
                     ENDIF

! Search indicated IJK cell for neighbors with particle LL (the
! current target particle).  LL's position has already been shifted 
! according to the boundary rules (e.g., wrap around for periodic)                  
                     CALL GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL(LL, &
                        POS_LL,I_INDEX,J_INDEX,K_INDEX)

                  ENDDO ILOOP3
               ENDDO JLOOP3
            ENDDO KLOOP3

         ENDIF   ! endif des_le_shear_dir == dudy 
      ENDIF  ! if des_dimn == 3

              
      RETURN

 1001 FORMAT(/1X,70('*')//,' From: GRID_BASED_NEIGHBOR_SEARCH ->',/,6X,&
         'GRIDBASED_LEBC_NEIGHBOR_CELL_SEARCH - ',/,&
         ' Message: ',A,' position is out of bounds for particle ',&
         I10,/,10X,' and should not be at this point',/,1X,70('*')/)
      
 1002 FORMAT(/1X,70('*')//,' From: GRID_BASED_NEIGHBOR_SEARCH ->',/,6X,&
         'GRIDBASED_LEBC_NEIGHBOR_CELL_SEARCH - ',/,&
         ' Message: the search index ',A,'=',I6, 'is out of the',&
         ' acceptable',/10X, 'range (',I6,'-',I6,') for particle ',&
         I10,/,1X,70('*')/)     

      END SUBROUTINE GRIDBASED_LEBC_NEIGHBOR_CELL_SEARCH



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                    
!  Module name: GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL
!  Purpose: Check given i, j, k cell for neighbors contacts with a given
!    particle at the given position.  If cell contains potential
!    neighbor contacts for given particle, store the information accordingly
!  
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
      
      SUBROUTINE GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL(LL, POS_LL,&
         I, J, K)

      USE compar         
      USE discretelement

      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------
! given particle number id for current target particle
      INTEGER LL   

! given position of target particle LL (adjusted for BC)
      DOUBLE PRECISION POS_LL(DIMN) 

! given loop index for x, y, z coordinate
      INTEGER I, J, K

! number of particles in target ijk cell      
      INTEGER NPG

! loop counters
      INTEGER NP, NEIGH_L

! particle number id for potential neighbor particle in cell ijk      
      INTEGER PNO  

! current number of neighbours of particle      
      INTEGER NLIM 

! contact variables
      DOUBLE PRECISION DISTVEC(DIMN), DIST, R_LM 

! if true particle LL and PNO have already been classified as neighbors      
      LOGICAL ALREADY_NEIGHBOURS

!-----------------------------------------------
! Functions
!-----------------------------------------------         
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------        


! If cell IJK contains particles, store the amount in NPG
      IF(ASSOCIATED(DESGRIDSEARCH_PIC(I,J,K)%P))THEN
         NPG = SIZE(DESGRIDSEARCH_PIC(I,J,K)%P)
      ELSE
         NPG = 0
      ENDIF

! Loop over all particles in IJK cell to determine if any are neighbors
! to particle LL (current target particle)
      DO NP = 1,NPG
         PNO = DESGRIDSEARCH_PIC(I,J,K)%P(NP)
         IF(.NOT.PEA(PNO,1)) CYCLE 

         IF(PNO.GT.LL)THEN 
! Calculate the distance between particles
            R_LM = DES_RADIUS(LL) + DES_RADIUS(PNO)
            R_LM = FACTOR_RLM*R_LM
            DISTVEC(:) = DES_POS_NEW(PNO,:) - POS_LL(:)
            DIST = SQRT(DES_DOTPRDCT(DISTVEC,DISTVEC))

! Check if LL and PNO are already neighbors
            ALREADY_NEIGHBOURS = .FALSE.
            DO NEIGH_L = 2, NEIGHBOURS(LL,1)+1
               IF(PNO .EQ. NEIGHBOURS(LL,NEIGH_L)) ALREADY_NEIGHBOURS=.TRUE.
            ENDDO

! Check the contact condition       
            IF( (DIST .LE. R_LM) .AND. &
                (.NOT.ALREADY_NEIGHBOURS)) THEN 

! Check that the array is not full and store PNO                           
               NEIGHBOURS(LL,1) = NEIGHBOURS(LL,1) + 1
               NLIM  = NEIGHBOURS(LL,1) + 1
               IF(NLIM.GT.MN) THEN 
                  WRITE(*,1000) NLIM, MN, LL
                  CALL MFIX_EXIT(myPE)
               ELSE
                  NEIGHBOURS(LL,NLIM) = PNO
               ENDIF

! Check that the array is not full and store LL
               NEIGHBOURS(PNO,1) = NEIGHBOURS(PNO,1) + 1
               NLIM  = NEIGHBOURS(PNO,1) + 1
               IF(NLIM.GT.MN) THEN 
                  WRITE(*,1000) NLIM, MN, PNO
                  CALL MFIX_EXIT(myPE)
               ELSE
                  NEIGHBOURS(PNO,NLIM) = LL
               ENDIF
            ENDIF  !contact condition

         ENDIF  !PNO.GT.LL

      ENDDO  !NP
 
      RETURN

 1000 FORMAT(/1X,70('*')//,' From: GRID_BASED_NEIGHBOR_SEARCH ->',/,10X,&
         'GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL - ',/,&
         ' Message: Neighbors(=', I4,') > MN(=', I4,&
         ') for particle: ', I10,/10X,&
         'Either reduce the factor R_LM or increase ',&
         'MN in mfix.dat',/1X,70('*')/)               

      END SUBROUTINE GRIDBASED_FIND_LL_NEIGHBORS_IN_NEIGHBOR_CELL
      



