!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_CELL (RMIN, REACTOR_LOC,D_DIR,N_DIR,CELL_LOC)     C
!  Purpose: calculate the cell index for a reactor location            C
!                                                                      C
!  Author: P. Nicoletti                               Date: 02-DEC-91  C
!  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:  None                                         C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, CELL_START, CELL_END                           C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE CALC_CELL(RMIN, REACTOR_LOC, D_DIR, N_DIR, CELL_LOC) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! passed arguments
!
!    RMIN       - starting value of the axis -- XMIN need not be zero,
!                   YMIN and ZMIN are assumed to be zero
!
!    REACTOR_LOC - location along one of the axis for which the cell
!                  index is to be found
!    D_DIR       - cell lengths (DX,DY,DZ)
!    N_DIR       - number of cells in this direction (IMAX,JMAX,KMAX)
!    CELL_LOC    - cell index corresponding to REACTOR_LOC
!
! local variables
!    LC         -  loop counter
!    CELL_START -  start coordinate for cell
!    CELL_END   -  end   coordinate for cell
      INTEGER N_DIR, CELL_LOC 
      DOUBLE PRECISION RMIN, REACTOR_LOC 

      DOUBLE PRECISION, DIMENSION(0:(N_DIR+3)) :: D_DIR 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LC 
      DOUBLE PRECISION :: CELL_START, CELL_END 
!-----------------------------------------------
!
      CELL_LOC = -1 
      CELL_START = RMIN 
      DO LC = 2, N_DIR + 1 
         CELL_END = CELL_START + D_DIR(LC) 
         IF (REACTOR_LOC <= CELL_START + HALF*D_DIR(LC)) THEN 
            CELL_LOC = LC - 1 
            RETURN  
         ELSE IF (REACTOR_LOC <= CELL_END + HALF*D_DIR(LC+1)) THEN 
            CELL_LOC = LC 
            RETURN  
         ENDIF 
         CELL_START = CELL_END 
      END DO 
      RETURN  
      END SUBROUTINE CALC_CELL 
!
!
! calculate reactor location from cell index
!
      SUBROUTINE CALC_LOC(RMIN, D_DIR, CELL_LOC, REACTOR_LOC) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
! passed arguments
!    RMIN       - starting value of the axis -- XMIN need not be zero,
!                   YMIN and ZMIN are assumed to be zero
!
!    REACTOR_LOC - location along one of the axis for which the cell
!                  index is to be found
!    D_DIR       - cell lengths (DX,DY,DZ)
!    CELL_LOC    - cell index corresponding to REACTOR_LOC
!
! local variables
!    LC         -  loop counter
!
      INTEGER CELL_LOC 
      DOUBLE PRECISION RMIN, REACTOR_LOC 

      DOUBLE PRECISION, DIMENSION(0:CELL_LOC) :: D_DIR
!  
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LC 
!-----------------------------------------------
!
      REACTOR_LOC = RMIN 
      LC = 2 
      IF (CELL_LOC - 1 > 0) THEN 
         REACTOR_LOC = REACTOR_LOC + SUM(D_DIR(2:CELL_LOC)) 
         LC = CELL_LOC + 1 
      ENDIF 
      RETURN  
      END SUBROUTINE CALC_LOC 
      
!// Comments on the modifications for DMP version implementation      
!// 100 EFD Replaced inheritance based dimensioning of arrays, i.e. D_DIR(*) 
