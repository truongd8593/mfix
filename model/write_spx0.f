!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_SPX0(L, unit_add)                                C
!  Purpose: write out the initial restart records (REAL)               C
!                                                                      C
!  Author: P. Nicoletti                               Date: 13-DEC-91  C
!  Reviewer: P. Nicoletti, W. Rogers, M. Syamlal      Date: 24-JAN-92  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced: RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR  C
!                        ID_MINUTE, ID_SECOND                          C
!  Variables modified: None                                            C
!                                                                      C
!  Local variables: LC, VERSION                                        C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!
      SUBROUTINE WRITE_SPX0(L, unit_add) 
!...Translated by Pacific-Sierra Research VAST-90 2.06G5  12:17:31  12/09/98  
!...Switches: -xf
!
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE param 
      USE param1 
      USE run
      USE funits 
      USE compar           !// 
      USE mpi_utility      !//
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER L 

!              offset for use in post_mfix
      INTEGER  unit_add
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!                file version ID
      CHARACTER :: VERSION*512 
      INTEGER  uspx   ! UNIT_SPX + offset from post_mfix
!-----------------------------------------------
!
        uspx = UNIT_SPX + unit_add

!
      if (myPE.ne.PE_IO) return    !// 
!
      VERSION = 'SPx = 02.00' 
      WRITE (VERSION(3:3), 1000) L 
      WRITE (uspx + L, REC=1) VERSION 
      WRITE (uspx + L, REC=2) RUN_NAME, ID_MONTH, ID_DAY, ID_YEAR, ID_HOUR&
         , ID_MINUTE, ID_SECOND 
!
!  The first field contains the pointer to the next record.
!  The second field contains the number of records written each time step
!  (The 4 and -1 will be overwritten in WRITE_SPX1)
!
      WRITE (uspx + L, REC=3) 4, -1 
      if(unit_add == 0) CALL FLUSH (uspx + L) 
 1000 FORMAT(I1) 
      RETURN  
      END SUBROUTINE WRITE_SPX0 
