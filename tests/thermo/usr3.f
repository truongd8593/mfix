!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is     C
!           user-definable.  The user may insert code in this routine   
!           or call appropriate user defined subroutines.  
!           This routine is not called from an IJK loop, hence  
!           all indices are undefined.                                 C
!                                                                      C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Revision Number:                                                    C
!  Purpose:                                                            C
!  Author:                                            Date: dd-mmm-yy  C
!  Reviewer:                                          Date: dd-mmm-yy  C
!                                                                      C
!  Literature/Document References:                                     C
!                                                                      C
!  Variables referenced:                                               C
!  Variables modified:                                                 C
!                                                                      C
!  Local variables:                                                    C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE USR3 

      use energy
      use rxns
      use usr

      IMPLICIT NONE

! Error index
      INTEGER :: IER
! fluid cell index and loop counter
      INTEGER :: IJK, LC
! total heat of reaction.
      DOUBLE PRECISION :: tHORg, tHORs(1:3)
! File unit
      INTEGER, parameter :: lUnit = 668

!      INCLUDE 'function.inc'

! Calculate user defined reaction rates.
      CALL RRATES0(IER)

      tHORg = 0.0d0
      tHORs = 0.0d0

      OPEN(lUnit,FILE='POST_Thermo.dat',access='APPEND')
      WRITE(lUnit,1100) 'Reaction','HORg', 'HORs1', 'HORs2', 'HORs3'

      DO LC=1, NO_OF_RXNS
         IJK = LC + 49
         WRITE(lUnit,1200) Reaction(LC)%Name(1:18), &
            HOR_g(IJK), HOR_s(IJK,1:3)
         tHORg = tHORg + HOR_g(IJK)
         tHORs = tHORs + HOR_s(IJK,1:3)
      ENDDO

      WRITE(lUnit,1300)'Total', tHORg, tHORs(1:3)
      CLOSE(lUnit)

      RETURN

 1100 FORMAT(3X,A,17X,A,3(11X,A))
 1200 FORMAT(3X,A,4(3X,G12.3))
 1300 FORMAT(90('-'),/3X,A,16X,G12.3,3(3X,G12.3))
      END SUBROUTINE USR3 
