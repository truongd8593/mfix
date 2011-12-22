!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WALLNODECONTACT(CIJK, L, WALLCONTACTI)                 C
!                                                                      C
!  Author: Tingwen Li                                 Date: 15-Jan-08  C
!>  Purposes: Check if particle contacts with the node of wall cell, if 
!>            it does, generate the wall particle                       
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
      

      RETURN
      END SUBROUTINE WALLNODECONTACT
