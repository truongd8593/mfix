CvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
C                                                                      C
C  Module name: SHIFT_DXYZ                                             C
C  Purpose:  shift the data in the dx,dy,dz arrays from 1:IMAX to      C
C            IMIN1:IMAX1,  1:JMAX to JMIN1:JMAX1 ,                     C
C            1:KMAX to KMIN1:KMAX1                                     C
C                                                                      C
C  Author: P. Nicoletti                               Date: 03-DEC-91  C
C  Reviewer: M.SYAMLAL, W.ROGERS, P.NICOLETTI         Date: 24-JAN-92  C
C                                                                      C
C  Revision Number:                                                    C
C  Purpose:                                                            C
C  Author:                                            Date: dd-mmm-yy  C
C  Reviewer:                                          Date: dd-mmm-yy  C
C                                                                      C
C  Literature/Document References:                                     C
C                                                                      C
C  Variables referenced: IMAX, IMAX1, IMAX2, JMAX, JMAX1, JMAX2, KMAX  C
C                        KMAX1 , KMAX2, IMIN1, JMIN1, KMIN1, NO_I,     C
C                        NO_J, NO_K                                    C
C  Variables modified:  DX, DY, DZ                                     C
C                                                                      C
C  Local variables: LC                                                 C
C                                                                      C
C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
C
      SUBROUTINE SHIFT_DXYZ
C
      IMPLICIT NONE
      INCLUDE  'param.inc'
      INCLUDE 'param1.inc'
      INCLUDE  'geometry.inc'
C
C              loop counter
      INTEGER  LC
C
      IF( DO_I) THEN
        DX(IMAX2) = DX(IMAX)
        DO 100 LC = IMAX1,IMIN1,-1
          DX(LC) = DX(LC-1)
100     CONTINUE
        DX(1)     = DX(IMIN1)
      ENDIF
C
      IF(DO_J) THEN
        DY(JMAX2) = DY(JMAX)
        DO 200 LC = JMAX1,JMIN1,-1
          DY(LC) = DY(LC-1)
200     CONTINUE
        DY(1)     = DY(JMIN1)
      ENDIF
C
      IF(DO_K) THEN
        DZ(KMAX2) = DZ(KMAX)
        DO 300 LC = KMAX1,KMIN1,-1
          DZ(LC) = DZ(LC-1)
300     CONTINUE
        DZ(1)     = DZ(KMIN1)
      ENDIF
C
      RETURN
      END
