!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: UPDATE_SEARCH_GRIDS                                    C
!  Purpose: This module updates the search grids used to identify      C
!           cohesion particle interactions                             C
!                                                                      C
!   Author: Mike Weber                              Date: 9/30/04      C
!   Reviewer:                                       Date:              C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE UPDATE_SEARCH_GRIDS

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE discretelement
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER i,j, NP
      INTEGER i_ind,j_ind,k_ind,np_ijk
! new position      
      INTEGER new(3)
      LOGICAL ALREADY_EXISTS
!-----------------------------------------------

      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**START UPDATE SEARCH GRIDS'
      END IF

! START LOOP OVER ALL PARTICLES
      DO NP=1,PARTICLES

         DO j=1,DIMN
            new(j)=INT(DES_POS_NEW(NP,j)/SEARCH_GRID_SIZE(j))+1
         ENDDO 
         IF(DIMN.eq.2)THEN
           new(3)=1
         ENDIF
        
! Did particle change grids?
         IF((new(1).ne.PART_GRID(NP,1)).OR.&
            (new(2).ne.PART_GRID(NP,2)).OR.&
            (new(3).ne.PART_GRID(NP,3)))THEN

            IF(COHESION_DEBUG.gt.2)THEN
               PRINT *,'***Particle=',NP
               PRINT *,'*****Old NP j=',PART_GRID(NP,1), PART_GRID(NP,2)
               PRINT *,'*****New NP j=',new(1), new(j)
            ENDIF

! Remove particle from old grid 
! ---------------------------------------------------------------->>>

! identify indices associated with particle old position 
            i_ind = PART_GRID(NP,1)
            j_ind = PART_GRID(NP,2)
            k_ind = PART_GRID(NP,3)
! number of particles at this ijk location on grid 
            np_ijk = PART_IN_GRID(i_ind,j_ind,k_ind,1)

! Check if particle is the last particle in the list at old grid
            IF(PART_GRID(NP,4).LT.(np_ijk+1)) THEN

! number of particles in grid+1?
               J=part_grid(NP,4)
               DO WHILE(j.le.(np_ijk+1))

! Shift particles to fill in vacant slot in grid
                  PART_IN_GRID(i_ind,j_ind,k_ind,j)=&
                     PART_IN_GRID(i_ind,j_ind,k_ind,j+1)

! Update grid slot recorded with each particle
                  IF(PART_IN_GRID(i_ind,j_ind,k_ind,j+1).GT.0)THEN
                     PART_GRID(PART_IN_GRID(i_ind,j_ind,k_ind,j+1),4)=j
                  ENDIF

                  j=j+1
               ENDDO  
            ELSE
               PART_IN_GRID(i_ind,j_ind,k_ind,PART_GRID(NP,4))=-1
            ENDIF

! Update count of particles in grid
            PART_IN_GRID(i_ind,j_ind,k_ind,1)=&
               PART_IN_GRID(i_ind,j_ind,k_ind,1)-1

! Remove particle from old grid 
! ----------------------------------------------------------------<<<

! Add particle to new grid 
! ---------------------------------------------------------------->>>
            DO j=1,DIMN
               PART_GRID(NP,j)=new(j)
            ENDDO

! identify indices associated with particle old position         
            i_ind = new(1)
            j_ind = new(2)
            k_ind = new(3)
! number of particles at this ijk location on grid 
            np_ijk = PART_IN_GRID(i_ind,j_ind,k_ind,1)

! place particle in last position in new grid   
            PART_IN_GRID(i_ind,j_ind,k_ind,np_ijk+2)=NP
   
! record particle's position in new grid with it's own info array
            PART_GRID(NP,4)=np_ijk+2

! Update count of particles in new grid
            PART_IN_GRID(i_ind,j_ind,k_ind,1)=np_ijk+1

! End add particle to new grid            
! ----------------------------------------------------------------<<<            


         ENDIF  ! end if particle moved out of current grid?  

      ENDDO  ! end do loop over all particles


      IF(COHESION_DEBUG.gt.0)THEN
         PRINT *,'**END UPDATE SEARCH GRIDS'
      ENDIF

      RETURN
      END SUBROUTINE UPDATE_SEARCH_GRIDS
