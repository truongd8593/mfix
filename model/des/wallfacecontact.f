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
      


      RETURN
      END SUBROUTINE WALLFACECONTACT


