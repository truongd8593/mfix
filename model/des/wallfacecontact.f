!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WALLFACECONTACT(WALL, L, WALLCONTACTI)                 C
!                                                                      C
!  Author: Tingwen Li                                 Date: 15-Jan-08  C
!>  Purposes: Check if particle contacts with the face of wall cell, if 
!>            it does, generate the wall particle                       
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!  Any questions, please contact litw@interchange.ubc.ca or dem@lists.mfix.org

      SUBROUTINE WALLFACECONTACT(WALL, L, WALLCONTACTI)

      Use discretelement
      USE param
      USE param1
      USE parallel
      USE fldvar
      USE run
      USE geometry
      USE matrix
      USE indices
      USE physprop
      USE drag
      USE constant
      USE compar
      IMPLICIT NONE
      
      INTEGER L,i,j,k,WALL,WALLCONTACTI
      double precision DES_R	
      
!     
!---------------------------------------------------------------------
!     Checking if a particle is in contact with any of the walls
!---------------------------------------------------------------------

      DES_R = DES_RADIUS(L)
      IF(WALLREFLECT) THEN
         DES_R = ZERO
      END IF

      DES_WALL_VEL(wall,1) = ZERO
      DES_WALL_VEL(wall,2) = ZERO
      IF(DIMN .EQ. 3) DES_WALL_VEL(wall,3) = ZERO      

      WALLCONTACTI = 0
      
      i=pijk(l,1)
      j=pijk(l,2)
!     Debugging
!     if(des_pos_new(l,1)-xe(i-1).le. 0 .or. des_pos_new(l,2)-yn(j-1).le.0  &
!     .or. xe(i)-des_pos_new(l,1).le.0 .or. yn(j)-des_pos_new(l,2).le.0)then
!     write(*,*)'particles exit domain'
!     write(*,109)wall,xe(i),yn(j),DES_POS_NEW(L,1),DES_POS_NEW(L,2)
!     stop
!     end if
!     if(l.eq.26)then
!     write(*,*)'particle 26'
!     write(*,109)wall,xe(i),yn(j),DES_POS_NEW(L,1),DES_POS_NEW(L,2)
!     end if
      
      IF(WALL .EQ. 1) THEN
         i=pijk(L,1) 
         IF((DES_POS_NEW(L,1)-xe(i-1)) .LE. DES_R) THEN
            WALLCONTACTI = 1
            DES_WALL_POS(wall,1) = xe(i-1) - DES_R	
            DES_WALL_POS(wall,2) = DES_POS_NEW(L,2)
            IF(DIMN .EQ. 3) DES_WALL_POS(wall,3) = DES_POS_NEW(L,3)
            WALL_NORMAL(wall,1) = -ONE
            WALL_NORMAL(wall,2) = ZERO
            IF(DIMN .EQ. 3) WALL_NORMAL(wall,3) = ZERO                        
         END IF
      ELSE IF(WALL .EQ. 2) THEN
         i=pijk(L,1)
         IF((xe(i)-DES_POS_NEW(L,1)) .LE. DES_R) THEN
            WALLCONTACTI = 1
            DES_WALL_POS(wall,1) = xe(i) + DES_R
            DES_WALL_POS(wall,2) = DES_POS_NEW(L,2)
            IF(DIMN .EQ. 3) DES_WALL_POS(wall,3) = DES_POS_NEW(L,3)
            WALL_NORMAL(wall,1) = ONE
            WALL_NORMAL(wall,2) = ZERO
            IF(DIMN .EQ. 3) WALL_NORMAL(wall,3) = ZERO            
         END IF

      ELSE IF(WALL .EQ. 3) THEN
         j=pijk(L,2)
         IF((DES_POS_NEW(L,2)-yn(j-1)) .LE. DES_R) THEN
            WALLCONTACTI = 1
            DES_WALL_POS(wall,1) = DES_POS_NEW(L,1)
            DES_WALL_POS(wall,2) = yn(j-1) - DES_R 
            IF(DIMN .EQ. 3) DES_WALL_POS(wall,3) = DES_POS_NEW(L,3)
            WALL_NORMAL(wall,1) = ZERO
            WALL_NORMAL(wall,2) = -ONE
            IF(DIMN .EQ. 3) WALL_NORMAL(wall,3) = ZERO            
         END IF

      ELSE IF(WALL .EQ. 4) THEN
         j=pijk(L,2)
         IF((yn(j)-DES_POS_NEW(L,2)) .LE. DES_R) THEN
            WALLCONTACTI = 1
            DES_WALL_POS(wall,1) = DES_POS_NEW(L,1)
            DES_WALL_POS(wall,2) = yn(j) + DES_R
            IF(DIMN .EQ. 3) DES_WALL_POS(wall,3) = DES_POS_NEW(L,3)
            WALL_NORMAL(wall,1) = ZERO
            WALL_NORMAL(wall,2) = ONE
            IF(DIMN .EQ. 3) WALL_NORMAL(wall,3) = ZERO            
         END IF

      ELSE IF(WALL .EQ. 5) THEN
         k=pijk(L,3)
         IF((DES_POS_NEW(L,3)-zt(k-1)) .LE. DES_R) THEN
            WALLCONTACTI = 1
            DES_WALL_POS(wall,1) = DES_POS_NEW(L,1)
            DES_WALL_POS(wall,2) = DES_POS_NEW(L,2)
            IF(DIMN .EQ. 3) DES_WALL_POS(wall,3) = zt(k-1) - DES_R
            WALL_NORMAL(wall,1) = ZERO
            WALL_NORMAL(wall,2) = ZERO
            IF(DIMN .EQ. 3) WALL_NORMAL(wall,3) = -ONE            
         END IF

      ELSE IF(WALL .EQ. 6) THEN
         k=pijk(L,3)
         IF((zt(k)-DES_POS_NEW(L,3)) .LE. DES_R) THEN
            WALLCONTACTI = 1
            DES_WALL_POS(wall,1) = DES_POS_NEW(L,1)
            DES_WALL_POS(wall,2) = DES_POS_NEW(L,2)
            IF(DIMN .EQ. 3) DES_WALL_POS(wall,3) = zt(k) + DES_R
            WALL_NORMAL(wall,1) = ZERO
            WALL_NORMAL(wall,2) = ZERO
            IF(DIMN .EQ. 3) WALL_NORMAL(wall,3) = ONE            
         END IF
      END IF

      RETURN
      END SUBROUTINE WALLFACECONTACT


