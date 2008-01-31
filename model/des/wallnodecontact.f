!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WALLNODECONTACT(CIJK, L, WALLCONTACTI)                 C
!                                                                      C
!  Author: Tingwen Li                                 Date: 15-Jan-08  C
!  Purposes: Check if particle contacts with the node of wall cell, if C 
!            it does, generate the wall particle                       C  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!  Any questions, please contact litw@interchange.ubc.ca or dem@lists.mfix.org

      SUBROUTINE WALLNODECONTACT(CIJK, L, WALLCONTACTI)

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
      
      INTEGER L, i,j,k, cijk, WALLCONTACTI
      double precision dispnd, DES_R	
      
!     
!---------------------------------------------------------------------
!     Checking if a particle is in contact with any of the walls
!---------------------------------------------------------------------
!     
      DES_R = DES_RADIUS(L)
      IF(WALLREFLECT) THEN
         DES_R = ZERO
      END IF

      DES_WALL_VEL(1,1) = ZERO
      DES_WALL_VEL(1,2) = ZERO
      IF(DIMN .EQ. 3) DES_WALL_VEL(1,3) = ZERO      

      WALLCONTACTI = 0
      
      i=pijk(L,1)+c_near_w(cijk,1)
      j=pijk(L,2)+c_near_w(cijk,2)
      k=pijk(L,3)+c_near_w(cijk,3)
      DISPND=sqrt((DES_POS_NEW(L,1)-xe(i))**2+(DES_POS_NEW(L,2)-yn(j))**2+(DES_POS_NEW(L,3)-zt(k))**2 ) 
      if(dispnd .le. DES_R)then
         wallcontacti=1
         DES_WALL_POS(1,1) = xe(i) + DES_R/dispnd*(xe(i)-DES_POS_NEW(L,1))
         DES_WALL_POS(1,2) = yn(j) + DES_R/dispnd*(yn(j)-DES_POS_NEW(L,2))						
         DES_WALL_POS(1,3) = zt(k) + DES_R/dispnd*(zt(k)-DES_POS_NEW(L,3))						
         WALL_NORMAL(1,1) = (xe(i)-DES_POS_NEW(L,1))/dispnd
         WALL_NORMAL(1,2) = (yn(j)-DES_POS_NEW(L,2))/dispnd
         WALL_NORMAL(1,3) = (zt(k)-DES_POS_NEW(L,3))/dispnd
      end if

      RETURN
      END SUBROUTINE WALLNODECONTACT
