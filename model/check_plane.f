!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CHECK_PLANE(X_CONSTANT,Y_CONSTANT,Z_CONSTANT,BC,NAME)  C
!  Purpose: make sure the flow boundary condition or internal surface  C
!           is a plane                                                 C
!                                                                      C
!  Author: P. Nicoletti                               Date: 10-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: None                                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: N                                                  C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CHECK_PLANE(X_CONSTANT, Y_CONSTANT, Z_CONSTANT, BC, NAME) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE funits 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!     surface indicators
      LOGICAL  X_CONSTANT,Y_CONSTANT,Z_CONSTANT
!
!     boundary condition or internal surface index
      INTEGER  BC
!
!     BC or IS
      CHARACTER*2 NAME
!
! local variables
!
!     number of directions that are not constant (must equal 2)
      INTEGER N
!-----------------------------------------------
!
!     boundary condition or internal surface index
!
!     BC or IS
!
! local variables
!
!     number of directions that are not constant (must equal 2)
!
      N = 3 
      IF (X_CONSTANT) N = N - 1 
      IF (Y_CONSTANT) N = N - 1 
      IF (Z_CONSTANT) N = N - 1 
!
      IF (N /= 2) THEN 
         WRITE (UNIT_LOG, 1000) NAME, BC 
         STOP  
      ENDIF 
!
      RETURN  
 1000 FORMAT(/70('*')//' From: CHECK_PLANE',/'Message: ',A,' No ',I3,&
         ' is not a plane',/70('*')/) 
      END SUBROUTINE CHECK_PLANE 
