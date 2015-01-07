!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: USR3                                                   C
!  Purpose: This routine is called after the time loop ends and is
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
!
      SUBROUTINE USR3

        USE param
        USE param1
        USE parallel
        USE constant
        USE run
        USE toleranc
        USE geometry
        USE indices
        USE compar
        USE mpi_utility
        USE sendrecv
        USE fldvar
        USE usr
        USE bc
        USE mms

      IMPLICIT NONE
!-----------------------------------------------
!
!  Include files defining common blocks here
!
!
!  Define local variables here
!
!
!  Include files defining statement functions here
!
!
!  Insert user-defined code here
!

!        Call WRITE_TECPLOT_DATA

        Call CALCULATE_DE_NORMS

        Deallocate(MMS_Ep_G)
        Deallocate(MMS_P_G)
        Deallocate(MMS_U_G)
        Deallocate(MMS_V_G)
        Deallocate(MMS_W_G)
        Deallocate(MMS_T_G)

        Deallocate(MMS_ROP_S)
        Deallocate(MMS_U_S)
        Deallocate(MMS_V_S)
        Deallocate(MMS_W_S)        
        Deallocate(MMS_T_S)        
        Deallocate(MMS_Theta_m)

        Deallocate(MMS_ROP_G_SRC)
        Deallocate(MMS_U_G_SRC)
        Deallocate(MMS_V_G_SRC)
        Deallocate(MMS_W_G_SRC)
        Deallocate(MMS_T_G_SRC)

        Deallocate(MMS_ROP_S_SRC)
        Deallocate(MMS_U_S_SRC)
        Deallocate(MMS_V_S_SRC)
        Deallocate(MMS_W_S_SRC)        
        Deallocate(MMS_T_S_SRC)        
        Deallocate(MMS_Theta_m_SRC)  

        Deallocate(P_G_Sh)

        Deallocate(xtr)
        Deallocate(ytr)
        Deallocate(ztr)

      RETURN
      END SUBROUTINE USR3

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: WRITE_TECPLOT_DATA                                     C
!  Purpose: Write Tecplot data for visualization		                   C
!                                                                      C
!  Author: Aniruddha Choudhary                        Date: 17-Oct-11  C
!  email: anirudd@vt.edu					                                     C
!  Reviewer:                                          Date:            C
!                                                                      C
!  Revision Number #                                  Date: ##-###-##  C
!  Author: anirudd@vt.edu                                              C
!  Purpose: Modified to write tecplot output for parallel runs         C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
!-----------------------------------------------

      SUBROUTINE WRITE_TECPLOT_DATA
        USE param
        USE param1
        USE parallel
        USE constant
        USE run
        USE toleranc
        USE geometry
        USE indices
        USE compar
        Use mpi_utility
        USE sendrecv
        USE fldvar
        USE mms
        USE bc
        USE physprop
      IMPLICIT NONE

        INTEGER :: I,J,K,IJK

        ! Temporary arrays to gather variables !
        Double precision, allocatable :: arr_Pg(:)
        Double precision, allocatable :: arr_MMSPg(:)
        Double precision, allocatable :: arr_PgSh(:)
        Double precision, allocatable :: arr_Epg(:)
        Double precision, allocatable :: arr_MMSEpg(:)
        Double precision, allocatable :: arr_Ug(:)
        Double precision, allocatable :: arr_MMSUg(:)
        Double precision, allocatable :: arr_Vg(:)
        Double precision, allocatable :: arr_MMSVg(:)
        Double precision, allocatable :: arr_Wg(:)
        Double precision, allocatable :: arr_MMSWg(:)
        Double precision, allocatable :: arr_ROPs(:)
        Double precision, allocatable :: arr_MMSROPs(:)
        Double precision, allocatable :: arr_Us(:)
        Double precision, allocatable :: arr_MMSUs(:)
        Double precision, allocatable :: arr_Vs(:)
        Double precision, allocatable :: arr_MMSVs(:)
        Double precision, allocatable :: arr_Ws(:)
        Double precision, allocatable :: arr_MMSWs(:)
        Double precision, allocatable :: arr_Tg(:)
        Double precision, allocatable :: arr_MMSTg(:)
        Double precision, allocatable :: arr_Ts(:)
        Double precision, allocatable :: arr_MMSTs(:)
        Double precision, allocatable :: arr_Ths(:)
        Double precision, allocatable :: arr_MMSThs(:)
        Double precision, allocatable :: arr_xtr(:)
        Double precision, allocatable :: arr_ytr(:)
        Double precision, allocatable :: arr_ztr(:)

        INTEGER :: IMJK, IJMK, IJKM

        DOUBLE PRECISION :: xt, yt, zt
                           ! temporary variables holding x, y, z locations

        !! X, Y coordinate variables !!
        DOUBLE PRECISION, DIMENSION(1:IMAX+1,1:JMAX+1,1:KMAX+1) :: &
                          X_temp, Y_temp, Z_temp

        ! Temporary 3D arrays !
        DOUBLE PRECISION, DIMENSION(1:IMAX+2,1:JMAX+2,1:KMAX+2) :: &
          PgSh_temp, Epg_temp, Ug_temp, Vg_temp, Wg_temp, &
          ROPS_temp, Us_temp, Vs_temp, Ws_temp, &
          Tg_temp, Ts_temp, Ths_temp, &
          MMSPg_temp, MMSEpg_temp, MMSUg_temp, MMSVg_temp, MMSWg_temp, &
          MMSROPS_temp, MMSUs_temp, MMSVs_temp, MMSWs_temp, &
          MMSTg_temp, MMSTs_temp, MMSThs_temp, &
          TempArr_x, TempArr_y, TempArr_z ! Temporary staggered mesh

        ! =.TRUE. plots vector data (U_G,U_S,etc.) at their computation
        ! location. Useful while debugging because this way staggered
        ! variables are actual values and not interpolated/avergaged
        ! values.
        LOGICAL         :: VECTOR_DATA_FLAG = .TRUE.

        include "function.inc"

        allocate(arr_Pg(ijkmax3))
        allocate(arr_MMSPg(ijkmax3))
        allocate(arr_PgSh(ijkmax3))

        allocate(arr_Epg(ijkmax3))
        allocate(arr_MMSEpg(ijkmax3))

        allocate(arr_Ug(ijkmax3))
        allocate(arr_MMSUg(ijkmax3))

        allocate(arr_Vg(ijkmax3))
        allocate(arr_MMSVg(ijkmax3))

        allocate(arr_Wg(ijkmax3))
        allocate(arr_MMSWg(ijkmax3))

        allocate(arr_ROPs(ijkmax3))
        allocate(arr_MMSROPs(ijkmax3))

        allocate(arr_Us(ijkmax3))
        allocate(arr_MMSUs(ijkmax3))

        allocate(arr_Vs(ijkmax3))
        allocate(arr_MMSVs(ijkmax3))

        allocate(arr_Ws(ijkmax3))
        allocate(arr_MMSWs(ijkmax3))

        allocate(arr_Tg(ijkmax3))
        allocate(arr_MMSTg(ijkmax3))

        allocate(arr_Ts(ijkmax3))
        allocate(arr_MMSTs(ijkmax3))

        allocate(arr_Ths(ijkmax3))
        allocate(arr_MMSThs(ijkmax3))

        allocate(arr_xtr(ijkmax3))
        allocate(arr_ytr(ijkmax3))
        allocate(arr_ztr(ijkmax3))        

        call gather(P_G,arr_Pg,root)
        call gather(MMS_P_G,arr_MMSPg,root)
        call gather(P_G_Sh,arr_PgSh,root)

        call gather(EP_g,arr_Epg,root)
        call gather(MMS_EP_G,arr_MMSEpg,root)

        call gather(U_G,arr_Ug,root)
        call gather(MMS_U_G,arr_MMSUg,root)

        call gather(V_G,arr_Vg,root)
        call gather(MMS_V_G,arr_MMSVg,root)

        call gather(W_G,arr_Wg,root)
        call gather(MMS_W_G,arr_MMSWg,root)

        call gather(ROP_s(:,1),arr_ROPs,root)
        call gather(MMS_ROP_s,arr_MMSROPs,root)

        call gather(U_s(:,1),arr_Us,root)
        call gather(MMS_U_s,arr_MMSUs,root)

        call gather(V_s(:,1),arr_Vs,root)
        call gather(MMS_V_s,arr_MMSVs,root)

        call gather(W_s(:,1),arr_Ws,root)
        call gather(MMS_W_s,arr_MMSWs,root)

        call gather(T_g,arr_Tg,root)
        call gather(MMS_T_g,arr_MMSTg,root)

        call gather(T_s(:,1),arr_Ts,root)
        call gather(MMS_T_s,arr_MMSTs,root)

        call gather(Theta_m(:,1),arr_Ths,root)
        call gather(MMS_Theta_m,arr_MMSThs,root)

        call gather(xtr,arr_xtr,root)
        call gather(ytr,arr_ytr,root)
        call gather(ztr,arr_ztr,root)

        if(myPE == PE_IO) then

          !! Shift pressure to reference value
          !! This will not be correct when cases are compressible.
          !IF(RO_G0 .NE. UNDEFINED) THEN
            K=1
            !DO K = 1, KMAX+2            
            DO J = 1, JMAX+2
            DO I = 1, IMAX+2
              IJK = FUNIJK_GL(I,J,K)
              arr_PgSh(IJK) = arr_Pg(IJK) + arr_MMSPg(IJK_Sh) - arr_Pg(IJK_Sh)
            END DO
            END DO
            !END DO
          !ELSE
          !  DO K = 1, KMAX+2
          !  DO J = 1, JMAX+2
          !  DO I = 1, IMAX+2
          !    IJK = FUNIJK_GL(I,J,K)
          !    arr_PgSh(IJK) = arr_Pg(IJK)
          !  END DO
          !  END DO
          !  END DO            
          !ENDIF

          !! Probably needs a better way to find empty units for file
          !output.

          OPEN(UNIT = 777, FILE="solution_cellc.dat", Status='unknown')
          Write(777,"(27a)") 'variables = "x""y""z"&
                &"Pg""Ug""Vg""Wg""Us""Vs""Ws"&
                &"Tg""Ts""epg""rops""Ths"&                
                &"MMSPg""MMSUg""MMSVg""MMSWg"&
                &"MMSUs""MMSVs""MMSWs"&
                &"MMSTg""MMSTs""MMSEpg""MMSRops""MMSThs"'
          Write(777,*) 'zone T="',0,'" '
          !Write(777,*) 'I=',IMAX,' J=',JMAX,' K=',KMAX
          Write(777,*) 'I=',IMAX,' J=',JMAX

!          Write(777,*) 'DATAPACKING=BLOCK'
!          WRITE(777,*) "VARLOCATION=([4,5,6,7,8,9,10,11,12,13,14,15,16,&
!          &17,18,19,20,21,22,23,24,25,26,27]=CELLCENTERED)"

          !We calculate and output cell-centered data at the cell-centers locations. 
          !Else, when cell-centered data is outputted with VARLOCATION=CELLCENTERED
          !option, tecplot does inaccurate interpolation near the boundary cells.
          !So current version is only useful for plotting all variables
          !in a single file but it is not the traditional tecplot
          !cell-centered data file. Some modification needed.

          ! Create temporary 3D arrays (mesh nodes)
          K=1
          !DO K = 1, KMAX+1
          DO J = 1, JMAX+1
          DO I = 1, IMAX+1
             IJK = FUNIJK_GL(I,J,K)
             X_temp(I,J,K) = arr_xtr(IJK)
             Y_temp(I,J,K) = arr_ytr(IJK)
             Z_temp(I,J,K) = arr_ztr(IJK)
          END DO
          END DO
          !END DO

          ! Create temporary 3D arrays (solution variables)
          K = 1
          !DO K = 2, KMAX+1
          DO J = 2, JMAX+1
          DO I = 2, IMAX+1

            ijk = funijk_gl(i,j,k)
            imjk = funijk_gl(i-1,j,k)
            ijmk = funijk_gl(i,j-1,k)
            !ijkm = funijk_gl(i,j,k-1)
            ijkm = funijk_gl(i,j,k)
            !imjk = im_of(ijk)
            !ijmk = jm_of(ijk)
            !ijkm = km_of(ijk)
   
            PgSh_temp(i,j,k) = arr_PgSh(IJK)
            MMSPg_temp(i,j,k) = arr_MMSPg(IJK)

            Epg_temp(i,j,k) = arr_Epg(IJK)
            MMSEpg_temp(i,j,k) = arr_MMSEpg(IJK)

            Ug_temp(i,j,k) = HALF*(arr_Ug(ijk)+arr_Ug(imjk))
            MMSUg_temp(i,j,k) = HALF*(arr_MMSUg(ijk)+arr_MMSUg(imjk))
            
            Vg_temp(i,j,k) = HALF*(arr_Vg(ijk)+arr_Vg(ijmk))
            MMSVg_temp(i,j,k) = HALF*(arr_MMSVg(ijk)+arr_MMSVg(ijmk))

            Wg_temp(i,j,k) = HALF*(arr_Wg(ijk)+arr_Wg(ijkm))
            MMSWg_temp(i,j,k) = HALF*(arr_MMSWg(ijk)+arr_MMSWg(ijkm))

            Rops_temp(i,j,k) = arr_ROPs(IJK)       
            MMSRops_temp(i,j,k) = arr_MMSROPs(IJK)

            Us_temp(i,j,k) = HALF*(arr_Us(ijk)+arr_Us(imjk))
            MMSUs_temp(i,j,k) = HALF*(arr_MMSUs(ijk)+arr_MMSUs(imjk))
            
            Vs_temp(i,j,k) = HALF*(arr_Vs(ijk)+arr_Vs(ijmk))
            MMSVs_temp(i,j,k) = HALF*(arr_MMSVs(ijk)+arr_MMSVs(ijmk))

            Ws_temp(i,j,k) = HALF*(arr_Ws(ijk)+arr_Ws(ijkm))
            MMSWs_temp(i,j,k) = HALF*(arr_MMSWs(ijk)+arr_MMSWs(ijkm))

            Tg_temp(i,j,k) = arr_Tg(IJK)
            MMSTg_temp(i,j,k) = arr_MMSTg(IJK)
            
            Ts_temp(i,j,k) = arr_Ts(IJK)
            MMSTs_temp(i,j,k) = arr_MMSTs(IJK)

            Ths_temp(i,j,k) = arr_Ths(IJK)
            MMSThs_temp(i,j,k) = arr_MMSThs(IJK)       

          END DO
          END DO
          !END DO

          ! write cell-centered data to file
          K = 1
          !DO K = 2, KMAX+1
          DO J = 2, JMAX+1
          DO I = 2, IMAX+1
            xt = HALF*(X_temp(i,j,k)+X_temp(i-1,j,k))
            yt = HALF*(Y_temp(i,j,k)+Y_temp(i,j-1,k))
            !zt = HALF*(Z_temp(i,j,k)+Z_temp(i,j,k-1))
            zt = HALF*(Z_temp(i,j,k)+Z_temp(i,j,k))
            Write(777,*) xt, yt, zt, &
              PgSh_temp(i,j,k), &
              Ug_temp(i,j,k), Vg_temp(i,j,k), Wg_temp(i,j,k),&
              Us_temp(i,j,k), Vs_temp(i,j,k), Ws_temp(i,j,k),&
              Tg_temp(i,j,k), Ts_temp(i,j,k), &
              Epg_temp(i,j,k), Rops_temp(i,j,k),&  
              Ths_temp(i,j,k), &
              MMSPg_temp(i,j,k), &
              MMSUg_temp(i,j,k), MMSVg_temp(i,j,k), MMSWg_temp(i,j,k),&
              MMSUs_temp(i,j,k), MMSVs_temp(i,j,k), MMSWs_temp(i,j,k),&
              MMSTg_temp(i,j,k), MMSTs_temp(i,j,k), &
              MMSEpg_temp(i,j,k), MMSRops_temp(i,j,k),&
              MMSThs_temp(i,j,k)
          END DO
          END DO
          !END DO

          CLOSE(777)


          IF(VECTOR_DATA_FLAG) THEN
          !! Output vector data (where calculated) !!
          !! UNDER CONSTRUCTION !!

          OPEN(UNIT = 782, FILE="solution_U.dat", Status='unknown')
          Write(782,"(9a)") 'variables = "x""y""z""Ug""Us""MMSUg""MMSUs"&
                &"UgErr""UsErr"'
          Write(782,*) 'zone T="',0,'" '
          !Write(782,*) 'I=',IMAX+2,' J=',JMAX+2,' K=',KMAX+2
          Write(782,*) 'I=',IMAX+2,' J=',JMAX+2

          K = 1
          !DO K = 1, KMAX+2
          DO J = 1, JMAX+2
          DO I = 1, IMAX+2

            IJK = FUNIJK_GL(I,J,K)

            TempArr_x(I,J,K) = arr_xtr(IJK)
            TempArr_y(I,J,K) = arr_ytr(IJK) - DY(J)*HALF
            TempArr_z(I,J,K) = arr_ztr(IJK) - DZ(K)*HALF

            Ug_temp(I,J,K) = arr_Ug(IJK)
            MMSUg_temp(I,J,K) = arr_MMSUg(IJK)

            Us_temp(I,J,K) = arr_Us(IJK)
            MMSUs_temp(I,J,K) = arr_MMSUs(IJK)

            Write(782,*) TempArr_x(I,J,K), TempArr_y(I,J,K), &
              TempArr_z(I,J,K), Ug_temp(I,J,K), Us_temp(I,J,K), &
              MMSUg_temp(I,J,K), MMSUs_temp(I,J,K), &
              Ug_temp(I,J,K)-MMSUg_temp(I,J,K), &
              Us_temp(I,J,K)-MMSUs_temp(I,J,K)

          END DO
          END DO
          !END DO

          CLOSE(782)

          OPEN(UNIT = 783, FILE="solution_V.dat", Status='unknown')
          Write(783,"(9a)") 'variables = "x""y""z""Vg""Vs""MMSVg""MMSVs"&
                &"VgErr""VsErr"'
          Write(783,*) 'zone T="',0,'" '
          !Write(783,*) 'I=',IMAX+2,' J=',JMAX+2,' K=',KMAX+2
          Write(783,*) 'I=',IMAX+2,' J=',JMAX+2

          K = 1
          !DO K = 1, KMAX+2
          DO J = 1, JMAX+2
          DO I = 1, IMAX+2

            IJK = FUNIJK_GL(I,J,K)

            TempArr_x(I,J,K) = arr_xtr(IJK) - DX(I)*HALF
            TempArr_y(I,J,K) = arr_ytr(IJK)
            TempArr_z(I,J,K) = arr_ztr(IJK) - DZ(K)*HALF

            Vg_temp(I,J,K) = arr_Vg(IJK)
            MMSVg_temp(I,J,K) = arr_MMSVg(IJK)

            Vs_temp(I,J,K) = arr_Vs(IJK)
            MMSVs_temp(I,J,K) = arr_MMSVs(IJK)

            Write(783,*) TempArr_x(I,J,K), TempArr_y(I,J,K), &
              TempArr_z(I,J,K), Vg_temp(I,J,K), Vs_temp(I,J,K), &
              MMSVg_temp(I,J,K), MMSVs_temp(I,J,K), &
              Vg_temp(I,J,K)-MMSVg_temp(I,J,K), &
              Vs_temp(I,J,K)-MMSVs_temp(I,J,K)

          END DO
          END DO
          !END DO

          CLOSE(783)

          OPEN(UNIT = 784, FILE="solution_W.dat", Status='unknown')
          Write(784,"(9a)") 'variables = "x""y""z""Wg""Ws""MMSWg""MMSWs"&
                &"WgErr""WsErr"'
          Write(784,*) 'zone T="',0,'" '
          !Write(784,*) 'I=',IMAX+2,' J=',JMAX+2,' K=',KMAX+2
          Write(784,*) 'I=',IMAX+2,' J=',JMAX+2

          K = 1
          !DO K = 1, KMAX+2
          DO J = 1, JMAX+2
          DO I = 1, IMAX+2

            IJK = FUNIJK_GL(I,J,K)

            TempArr_x(I,J,K) = arr_xtr(IJK) - DX(I)*HALF
            TempArr_y(I,J,K) = arr_ytr(IJK) - DY(J)*HALF
            TempArr_z(I,J,K) = arr_ztr(IJK)

            Wg_temp(I,J,K) = arr_Wg(IJK)
            MMSWg_temp(I,J,K) = arr_MMSWg(IJK)

            Ws_temp(I,J,K) = arr_Ws(IJK)
            MMSWs_temp(I,J,K) = arr_MMSWs(IJK)

            Write(784,*) TempArr_x(I,J,K), TempArr_y(I,J,K), &
              TempArr_z(I,J,K), Wg_temp(I,J,K), Ws_temp(I,J,K), &
              MMSWg_temp(I,J,K), MMSWs_temp(I,J,K), &
              Wg_temp(I,J,K)-MMSWg_temp(I,J,K), &
              Ws_temp(I,J,K)-MMSWs_temp(I,J,K)

          END DO
          END DO
          !END DO

          CLOSE(784)
         
          ENDIF ! end of IF(VECTOR_DATA_FLAG)

        Endif

        deallocate(arr_Pg)
        deallocate(arr_MMSPg)
        deallocate(arr_PgSh)

        deallocate(arr_Epg)
        deallocate(arr_MMSEpg)

        deallocate(arr_Ug)
        deallocate(arr_MMSUg)

        deallocate(arr_Vg)
        deallocate(arr_MMSVg)

        deallocate(arr_Wg)
        deallocate(arr_MMSWg)

        deallocate(arr_ROPs)
        deallocate(arr_MMSROPs)

        deallocate(arr_Us)
        deallocate(arr_MMSUs)

        deallocate(arr_Vs)
        deallocate(arr_MMSVs)

        deallocate(arr_Ws)
        deallocate(arr_MMSWs)

        deallocate(arr_Tg)
        deallocate(arr_MMSTg)

        deallocate(arr_Ts)
        deallocate(arr_MMSTs)

        deallocate(arr_Ths)
        deallocate(arr_MMSThs)

        deallocate(arr_xtr)
        deallocate(arr_ytr)
        deallocate(arr_ztr)

      END SUBROUTINE WRITE_TECPLOT_DATA
!-----------------------------------------------


      SUBROUTINE CALCULATE_DE_NORMS
        USE param
        USE param1
        USE parallel
        USE constant
        USE run
        USE toleranc
        USE geometry
        USE indices
        USE compar
        USE mpi_utility
        USE sendrecv
        USE fldvar
        USE mms
        USE bc
      IMPLICIT NONE

        INTEGER :: I,J,K,IJK

        ! Temporary arrays to gather variables !
        Double precision, allocatable :: arr_Pg(:)
        Double precision, allocatable :: arr_MMSPg(:)
        Double precision, allocatable :: arr_PgSh(:)
        Double precision, allocatable :: arr_Epg(:)
        Double precision, allocatable :: arr_MMSEpg(:)
        Double precision, allocatable :: arr_Ug(:)
        Double precision, allocatable :: arr_MMSUg(:)
        Double precision, allocatable :: arr_Vg(:)
        Double precision, allocatable :: arr_MMSVg(:)
        Double precision, allocatable :: arr_Wg(:)
        Double precision, allocatable :: arr_MMSWg(:)
        Double precision, allocatable :: arr_ROPs(:)
        Double precision, allocatable :: arr_MMSROPs(:)
        Double precision, allocatable :: arr_Us(:)
        Double precision, allocatable :: arr_MMSUs(:)
        Double precision, allocatable :: arr_Vs(:)
        Double precision, allocatable :: arr_MMSVs(:)
        Double precision, allocatable :: arr_Ws(:)
        Double precision, allocatable :: arr_MMSWs(:)
        Double precision, allocatable :: arr_Tg(:)
        Double precision, allocatable :: arr_MMSTg(:)
        Double precision, allocatable :: arr_Ts(:)
        Double precision, allocatable :: arr_MMSTs(:)
        Double precision, allocatable :: arr_Ths(:)
        Double precision, allocatable :: arr_MMSThs(:)

        DOUBLE PRECISION, DIMENSION(1:12) :: L1DE, L2DE, LinfDE
        ! 1 = P_g, 2 = U_g, 3 = V_g, 4 = W_g, 5 = U_s, 6 = V_s, 7 = W_s !
        ! 8 = T_g, 9 = T_s, 10 = EP_g, 11 = ROP_s, 12 = Th_s

        include "function.inc"

        allocate(arr_Pg(ijkmax3))
        allocate(arr_MMSPg(ijkmax3))
        allocate(arr_PgSh(ijkmax3))

        allocate(arr_Epg(ijkmax3))
        allocate(arr_MMSEpg(ijkmax3))

        allocate(arr_Ug(ijkmax3))
        allocate(arr_MMSUg(ijkmax3))

        allocate(arr_Vg(ijkmax3))
        allocate(arr_MMSVg(ijkmax3))

        allocate(arr_Wg(ijkmax3))
        allocate(arr_MMSWg(ijkmax3))

        allocate(arr_ROPs(ijkmax3))
        allocate(arr_MMSROPs(ijkmax3))

        allocate(arr_Us(ijkmax3))
        allocate(arr_MMSUs(ijkmax3))

        allocate(arr_Vs(ijkmax3))
        allocate(arr_MMSVs(ijkmax3))

        allocate(arr_Ws(ijkmax3))
        allocate(arr_MMSWs(ijkmax3))

        allocate(arr_Tg(ijkmax3))
        allocate(arr_MMSTg(ijkmax3))

        allocate(arr_Ts(ijkmax3))
        allocate(arr_MMSTs(ijkmax3))

        allocate(arr_Ths(ijkmax3))
        allocate(arr_MMSThs(ijkmax3))  

        call gather(P_G,arr_Pg,root)
        call gather(MMS_P_G,arr_MMSPg,root)
        call gather(P_G_Sh,arr_PgSh,root)

        call gather(EP_g,arr_Epg,root)
        call gather(MMS_EP_G,arr_MMSEpg,root)

        call gather(U_G,arr_Ug,root)
        call gather(MMS_U_G,arr_MMSUg,root)

        call gather(V_G,arr_Vg,root)
        call gather(MMS_V_G,arr_MMSVg,root)

        call gather(W_G,arr_Wg,root)
        call gather(MMS_W_G,arr_MMSWg,root)

        call gather(ROP_s(:,1),arr_ROPs,root)
        call gather(MMS_ROP_s,arr_MMSROPs,root)

        call gather(U_s(:,1),arr_Us,root)
        call gather(MMS_U_s,arr_MMSUs,root)

        call gather(V_s(:,1),arr_Vs,root)
        call gather(MMS_V_s,arr_MMSVs,root)

        call gather(W_s(:,1),arr_Ws,root)
        call gather(MMS_W_s,arr_MMSWs,root)

        call gather(T_g,arr_Tg,root)
        call gather(MMS_T_g,arr_MMSTg,root)

        call gather(T_s(:,1),arr_Ts,root)
        call gather(MMS_T_s,arr_MMSTs,root)

        call gather(Theta_m(:,1),arr_Ths,root)
        call gather(MMS_Theta_m,arr_MMSThs,root)

        If(myPE == PE_IO) then

          !! Shift pressure to reference value
          !! This will not be correct when cases are compressible.
          !IF(RO_G0 .NE. UNDEFINED) THEN
            K = 1
            !DO K = 1, KMAX+2
            DO J = 1, JMAX+2
            DO I = 1, IMAX+2
              IJK = FUNIJK_GL(I,J,K)
              arr_PgSh(IJK) = arr_Pg(IJK) + arr_MMSPg(IJK_Sh) - arr_Pg(IJK_Sh)
            END DO
            END DO
            !END DO
          !ELSE
          !  DO K = 1, KMAX+2
          !  DO J = 1, JMAX+2
          !  DO I = 1, IMAX+2
          !    IJK = FUNIJK_GL(I,J,K)
          !    arr_PgSh(IJK) = arr_Pg(IJK)
          !  END DO
          !  END DO
          !  END DO            
          !ENDIF

          !! Calculate and output DE Norms !!
          ! For P_G !
          L1DE(1) = ZERO
          L2DE(1) = ZERO
          LinfDE(1) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
            IJK = FUNIJK_GL(I,J,K)
            L1DE(1) = L1DE(1) + abs(arr_PgSh(IJK) - arr_MMSPg(IJK))
            L2DE(1) = L2DE(1) + abs(arr_PgSh(IJK) - arr_MMSPg(IJK))**2
            LinfDE(1) = MAX(LinfDE(1), abs(arr_PgSh(IJK) - arr_MMSPg(IJK)))
          !END DO
          END DO
          END DO
          L1DE(1) = L1DE(1)/REAL((IMAX)*(JMAX)*(KMAX))
          L2DE(1) = SQRT(L2DE(1)/REAL((IMAX)*(JMAX)*(KMAX)))

          ! For U_g !
          L1DE(2) = ZERO
          L2DE(2) = ZERO
          LinfDE(2) = ZERO
          DO I = 1, IMAX+1
          DO J = 2, JMAX+1
          !DO K = 2, KMAX+1
          K = 1
             IJK = funijk_gl(i,j,k)
             L1DE(2) = L1DE(2) + abs(arr_Ug(IJK) - arr_MMSUg(IJK))
             L2DE(2) = L2DE(2) + abs(arr_Ug(IJK) - arr_MMSUg(IJK))**2
             LinfDE(2) = MAX(LinfDE(2), abs(arr_Ug(IJK) - arr_MMSUg(IJK)))
          !END DO
          END DO
          END DO
          L1DE(2) = L1DE(2)/REAL((IMAX-1)*(JMAX)*(KMAX))
          L2DE(2) = SQRT(L2DE(2)/REAL((IMAX-1)*(JMAX)*(KMAX)))

          ! For V_g !
          L1DE(3) = ZERO
          L2DE(3) = ZERO
          LinfDE(3) = ZERO
          DO I = 2, IMAX+1
          DO J = 1, JMAX+1
          !DO K = 2, KMAX+1
          K = 1
             IJK = funijk_gl(i,j,k)
             L1DE(3) = L1DE(3) + abs(arr_Vg(IJK) - arr_MMSVg(IJK))
             L2DE(3) = L2DE(3) + abs(arr_Vg(IJK) - arr_MMSVg(IJK))**2
             LinfDE(3) = MAX(LinfDE(3), abs(arr_Vg(IJK) - arr_MMSVg(IJK)))
          !END DO
          END DO
          END DO
          L1DE(3) = L1DE(3)/REAL((IMAX)*(JMAX-1)*(KMAX))
          L2DE(3) = SQRT(L2DE(3)/REAL((IMAX)*(JMAX-1)*(KMAX)))

          ! For W_g !
          L1DE(4) = ZERO
          L2DE(4) = ZERO
          LinfDE(4) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          !DO K = 2, KMAX
          K = 1
             IJK = funijk_gl(i,j,k)
             L1DE(4) = L1DE(4) + abs(arr_Wg(IJK) - arr_MMSWg(IJK))
             L2DE(4) = L2DE(4) + abs(arr_Wg(IJK) - arr_MMSWg(IJK))**2
             LinfDE(4) = MAX(LinfDE(4), abs(arr_Wg(IJK) - arr_MMSWg(IJK)))
          !END DO
          END DO
          END DO
          !L1DE(4) = L1DE(4)/REAL((IMAX)*(JMAX)*(KMAX-1))
          !L2DE(4) = SQRT(L2DE(4)/REAL((IMAX)*(JMAX)*(KMAX-1)))
          L1DE(4) = L1DE(4)/REAL((IMAX)*(JMAX))
          L2DE(4) = SQRT(L2DE(4)/REAL((IMAX)*(JMAX)))

          ! For U_s !
          L1DE(5) = ZERO
          L2DE(5) = ZERO
          LinfDE(5) = ZERO
          DO I = 1, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
             IJK = funijk_gl(i,j,k)
             L1DE(5) = L1DE(5) + abs(arr_Us(IJK) - arr_MMSUs(IJK))
             L2DE(5) = L2DE(5) + abs(arr_Us(IJK) - arr_MMSUs(IJK))**2
             LinfDE(5) = MAX(LinfDE(5), abs(arr_Us(IJK) - arr_MMSUs(IJK)))
          !END DO
          END DO
          END DO
          L1DE(5) = L1DE(5)/REAL((IMAX-1)*(JMAX)*(KMAX))
          L2DE(5) = SQRT(L2DE(5)/REAL((IMAX-1)*(JMAX)*(KMAX)))

          ! For V_s !
          L1DE(6) = ZERO
          L2DE(6) = ZERO
          LinfDE(6) = ZERO
          DO I = 2, IMAX+1
          DO J = 1, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
             IJK = funijk_gl(i,j,k)
             L1DE(6) = L1DE(6) + abs(arr_Vs(IJK) - arr_MMSVs(IJK))
             L2DE(6) = L2DE(6) + abs(arr_Vs(IJK) - arr_MMSVs(IJK))**2
             LinfDE(6) = MAX(LinfDE(6), abs(arr_Vs(IJK) - arr_MMSVs(IJK)))
          !END DO
          END DO
          END DO
          L1DE(6) = L1DE(6)/REAL((IMAX)*(JMAX-1)*(KMAX))
          L2DE(6) = SQRT(L2DE(6)/REAL((IMAX)*(JMAX-1)*(KMAX)))

          ! For W_s !
          L1DE(7) = ZERO
          L2DE(7) = ZERO
          LinfDE(7) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX
             IJK = funijk_gl(i,j,k)
             L1DE(7) = L1DE(7) + abs(arr_Ws(IJK) - arr_MMSWs(IJK))
             L2DE(7) = L2DE(7) + abs(arr_Ws(IJK) - arr_MMSWs(IJK))**2
             LinfDE(7) = MAX(LinfDE(7), abs(arr_Ws(IJK) - arr_MMSWs(IJK)))
          !END DO
          END DO
          END DO
          !L1DE(7) = L1DE(7)/REAL((IMAX)*(JMAX)*(KMAX-1))
          !L2DE(7) = SQRT(L2DE(7)/REAL((IMAX)*(JMAX)*(KMAX-1)))
          L1DE(7) = L1DE(7)/REAL((IMAX)*(JMAX))
          L2DE(7) = SQRT(L2DE(7)/REAL((IMAX)*(JMAX)))

          ! For T_G !
          L1DE(8) = ZERO
          L2DE(8) = ZERO
          LinfDE(8) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
            IJK = FUNIJK_GL(I,J,K)
            L1DE(8) = L1DE(8) + abs(arr_Tg(IJK) - arr_MMSTg(IJK))
            L2DE(8) = L2DE(8) + abs(arr_Tg(IJK) - arr_MMSTg(IJK))**2
            LinfDE(8) = MAX(LinfDE(8), abs(arr_Tg(IJK) - arr_MMSTg(IJK)))
          !END DO
          END DO
          END DO
          L1DE(8) = L1DE(8)/REAL((IMAX)*(JMAX)*(KMAX))
          L2DE(8) = SQRT(L2DE(8)/REAL((IMAX)*(JMAX)*(KMAX)))

          ! For T_S !
          L1DE(9) = ZERO
          L2DE(9) = ZERO
          LinfDE(9) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
            IJK = FUNIJK_GL(I,J,K)
            L1DE(9) = L1DE(9) + abs(arr_Ts(IJK) - arr_MMSTs(IJK))
            L2DE(9) = L2DE(9) + abs(arr_Ts(IJK) - arr_MMSTs(IJK))**2
            LinfDE(9) = MAX(LinfDE(9), abs(arr_Ts(IJK) - arr_MMSTs(IJK)))
          !END DO
          END DO
          END DO
          L1DE(9) = L1DE(9)/REAL((IMAX)*(JMAX)*(KMAX))
          L2DE(9) = SQRT(L2DE(9)/REAL((IMAX)*(JMAX)*(KMAX)))

          ! For EP_G !
          L1DE(10) = ZERO
          L2DE(10) = ZERO
          LinfDE(10) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
            IJK = FUNIJK_GL(I,J,K)
            L1DE(10) = L1DE(10) + abs(arr_Epg(IJK) - arr_MMSEpg(IJK))
            L2DE(10) = L2DE(10) + abs(arr_Epg(IJK) - arr_MMSEpg(IJK))**2
            LinfDE(10) = MAX(LinfDE(10), abs(arr_Epg(IJK) - arr_MMSEpg(IJK)))
          !END DO
          END DO
          END DO
          L1DE(10) = L1DE(10)/REAL((IMAX)*(JMAX)*(KMAX))
          L2DE(10) = SQRT(L2DE(10)/REAL((IMAX)*(JMAX)*(KMAX)))

          ! For ROP_S !
          L1DE(11) = ZERO
          L2DE(11) = ZERO
          LinfDE(11) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
            IJK = FUNIJK_GL(I,J,K)
            L1DE(11) = L1DE(11) + abs(arr_ROPs(IJK) - arr_MMSROPs(IJK))
            L2DE(11) = L2DE(11) + abs(arr_ROPs(IJK) - arr_MMSROPs(IJK))**2
            LinfDE(11) = MAX(LinfDE(11), abs(arr_ROPs(IJK) - arr_MMSROPs(IJK)))
          !END DO
          END DO
          END DO
          L1DE(11) = L1DE(11)/REAL((IMAX)*(JMAX)*(KMAX))
          L2DE(11) = SQRT(L2DE(11)/REAL((IMAX)*(JMAX)*(KMAX)))

          ! For TH_S !
          L1DE(12) = ZERO
          L2DE(12) = ZERO
          LinfDE(12) = ZERO
          DO I = 2, IMAX+1
          DO J = 2, JMAX+1
          K = 1
          !DO K = 2, KMAX+1
            IJK = FUNIJK_GL(I,J,K)
            L1DE(12) = L1DE(12) + abs(arr_Ths(IJK) - arr_MMSThs(IJK))
            L2DE(12) = L2DE(12) + abs(arr_Ths(IJK) - arr_MMSThs(IJK))**2
            LinfDE(12) = MAX(LinfDE(12), abs(arr_Ths(IJK) - arr_MMSThs(IJK)))
          !END DO
          END DO
          END DO
          L1DE(12) = L1DE(12)/REAL((IMAX)*(JMAX)*(KMAX))
          L2DE(12) = SQRT(L2DE(12)/REAL((IMAX)*(JMAX)*(KMAX)))

          OPEN(UNIT = 779, FILE="de_norms.dat", Status='unknown')
          Write(779,*) "IMAX= ",IMAX, " JMAX=", JMAX, " KMAX=", KMAX
          !Write(779,*) "1st line: L1 Norms, 2nd line: L2 Norms,&
          !              & 3rd line: Linf Norms"
          !Write(779,*) "Columns: P_g : U_g : V_g : W_g : U_s : V_s : W_s &
          !  & : T_g : T_s : Ep_g : ROP_s : Th_S "
          Write(779,200) L1DE(1), L1DE(2), L1DE(3) !, L1DE(4) , L1DE(5), &
            !&L1DE(6), L1DE(7), L1DE(8), L1DE(9), L1DE(10), L1DE(11), &
            !&L1DE(12)
          Write(779,200) L2DE(1), L2DE(2), L2DE(3) !, L2DE(4) , L2DE(5), &
            !&L2DE(6), L2DE(7), L2DE(8), L2DE(9), L2DE(10), L2DE(11), &
            !&L2DE(12)
          Write(779,200) LinfDE(1), LinfDE(2), LinfDE(3) !, LinfDE(4) , &
            !&LinfDE(5), LinfDE(6), LinfDE(7), LinfDE(8), LinfDE(9), &
            !&LinfDE(10), LinfDE(11), LinfDE(12)
          CLOSE(779)

        200 Format (12Es24.16)

        Endif

        deallocate(arr_Pg)
        deallocate(arr_MMSPg)
        deallocate(arr_PgSh)

        deallocate(arr_Epg)
        deallocate(arr_MMSEpg)

        deallocate(arr_Ug)
        deallocate(arr_MMSUg)

        deallocate(arr_Vg)
        deallocate(arr_MMSVg)

        deallocate(arr_Wg)
        deallocate(arr_MMSWg)

        deallocate(arr_ROPs)
        deallocate(arr_MMSROPs)

        deallocate(arr_Us)
        deallocate(arr_MMSUs)

        deallocate(arr_Vs)
        deallocate(arr_MMSVs)

        deallocate(arr_Ws)
        deallocate(arr_MMSWs)

        deallocate(arr_Tg)
        deallocate(arr_MMSTg)

        deallocate(arr_Ts)
        deallocate(arr_MMSTs)

        deallocate(arr_Ths)
        deallocate(arr_MMSThs)

      END SUBROUTINE CALCULATE_DE_NORMS
!-----------------------------------------------
