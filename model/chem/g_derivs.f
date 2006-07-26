!                           DISCLAIMER
!
!   This file was generated on 10/26/04 by the version of
!   ADIFOR compiled on June, 1998.
!
!   ADIFOR was prepared as an account of work sponsored by an
!   agency of the United States Government, Rice University, and
!   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
!   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
!   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
!   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
!   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
!   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
!   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
!
      subroutine g_derivs(g_p_, y, g_y, ldg_y, yt, g_yt, ldg_yt, n_sh)
!
!
        implicit none
        double precision x, y(12), yt(12)
!
        double precision n_sh
        double precision ro_s
        double precision mw_g, mw_g1, mw_g2, mw_g3
        double precision mw_g4, mw_g5, mw_s1, mw_max
        double precision d_p
        double precision p_g
        double precision diff_sih2, diff_sih4
        double precision r_sih2, r_sih4
        double precision p_h2, p_sih2
        double precision p_sih4
        double precision raf, rab, rbf, rbb, rc, rd
        double precision rg(5), rs, rs1
        double precision c_pg, c_ps
        double precision cpsih4, cpsih2, cph2, cpsi2h6, cpn2
        double precision cpsi, cpal2o3
!
        integer g_pmax_
        parameter (g_pmax_ = 12)
        integer g_i_, g_p_, ldg_y, ldg_yt
        double precision d22_b, d1_w, d2_w, d19_b, d18_b, d17_b, d22_v, &
     &d21_v, d2_v, d3_v
        double precision d4_v, d5_v, d6_v, d7_v, d8_v, d9_v, d10_v, d11_&
     &v, d12_v, d17_v
        double precision d14_v, d18_v, d2_b, d3_b, d4_b, d5_b, d6_b, d7_&
     &b, d8_b, d9_b
        double precision d10_b, d11_b, d12_b, d13_b, d14_b, d1_p, d2_p, &
     &d15_v, d16_v, d15_b
        double precision d16_b, g_mw_g(g_pmax_), g_y(ldg_y, 12), g_p_g(g&
     &_pmax_), g_diff_sih2(g_pmax_), g_p_sih4(g_pmax_), g_p_h2(g_pmax_),&
     & g_p_sih2(g_pmax_), g_cph2(g_pmax_), g_cpn2(g_pmax_)
        double precision g_c_pg(g_pmax_), g_c_ps(g_pmax_), g_d1_w(g_pmax&
     &_), g_raf(g_pmax_), g_rab(g_pmax_), g_rbf(g_pmax_), g_rbb(g_pmax_)&
     &, g_d2_w(g_pmax_), g_rc(g_pmax_), g_rd(g_pmax_)
        double precision g_rg(g_pmax_, 5), g_rs(g_pmax_), g_yt(ldg_yt, 1&
     &2)
        integer g_ehfid
        save g_d1_w, g_raf, g_rab, g_rbf, g_rbb, g_d2_w, g_rc, g_rd, g_r&
     &g, g_rs
        save g_mw_g, g_p_g, g_diff_sih2, g_p_sih4, g_p_h2, g_p_sih2, g_c&
     &ph2, g_cpn2, g_c_pg, g_c_ps
        intrinsic dble
        data g_ehfid /0/

!        call ehsfid(g_ehfid, 'derivs','g_derivs.f')

        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        ro_s = 3.9d0
        d_p = 0.0082
        mw_g1 = 32.13
        mw_g2 = 30.11
        mw_g3 = 2.02
        mw_g4 = 62.24
        mw_g5 = 28.01
        mw_s1 = 28.09

        mw_max = 1.0d0 / 500.0d0
        d4_b = 1.0d0 / mw_g5
        d7_b = 1.0d0 / mw_g4
        d10_b = 1.0d0 / mw_g3
        d13_b = 1.0d0 / mw_g2
        d14_b = 1.0d0 / mw_g1
        do g_i_ = 1, g_p_
          g_mw_g(g_i_) = d4_b * g_y(g_i_, 7) + d7_b * g_y(g_i_, 6) + d10&
     &_b * g_y(g_i_, 5) + d13_b * g_y(g_i_, 4) + d14_b * g_y(g_i_, 3)
        enddo
        mw_g = y(3) / mw_g1 + y(4) / mw_g2 + y(5) / mw_g3 + y(6) / mw_g4&
     & + y(7) / mw_g5
!--------
        d2_v = max (mw_g, mw_max)

        if (mw_g .gt.  mw_max) then
           d1_p = 1.0d0
           d2_p = 0.0d0
        else if (mw_g .lt.  mw_max) then
           d1_p = 0.0d0
           d2_p = 1.0d0
!        else
!           call ehbfDO (7,mw_g, mw_max, d2_v, d1_p, d2_p,
!     +g_ehfid,
!     +104)
!           d2_p = 1.0d0 -  d1_p
        endif
        d3_v = 1.0d0 / d2_v
        d3_b = (-d3_v) / d2_v * d1_p
        do g_i_ = 1, g_p_
          g_mw_g(g_i_) = d3_b * g_mw_g(g_i_)
        enddo
        mw_g = d3_v
!--------
        d3_v = 8314.56d4 * y(2)
        d6_v = y(1) * d3_v / mw_g
        d2_b = 1.0d0 / mw_g
        d3_b = (-d6_v) / mw_g
        d4_b = d2_b * d3_v
        d6_b = d2_b * y(1) * 8314.56d4
        do g_i_ = 1, g_p_
          g_p_g(g_i_) = d3_b * g_mw_g(g_i_) + d6_b * g_y(g_i_, 2) + d4_b&
     & * g_y(g_i_, 1)
        enddo
        p_g = d6_v
!--------

        r_sih2 = (8314.7295 / 30.11)
        r_sih4 = (8314.7295 / 32.12)
        d2_v = dble(1.5)

!        if ( y(2) .ne. 0.0d0 ) then
           d3_v = y(2) ** ( d2_v - 2.0d0)
           d3_v =  d3_v * y(2)
           d1_p =  d2_v *  d3_v
           d3_v =  d3_v * y(2)
!        else
!          (y(2) = 0)
!           d3_v = y(2) **  d2_v

!           if (  d2_v .lt. 1.0d0 ) then
!              call ehbfDO (10,y(2), d2_v, d3_v, d1_p, 0.0d0,
!     +g_ehfid,
!     +143)
!           else if (  d2_v .lt. 2.0d0 ) then
!              d1_p = 0.0d0
!              call ehbfDO (10,y(2), d2_v, d3_v, d1_p, 0.0d0,
!     +g_ehfid,
!     +148)
!           else
!              d1_p = 0.0d0
!           endif
!        endif
        d4_b = 0.0d0
        d2_b = dble(3.45e-5)
        d3_b = d2_b * d1_p
        do g_i_ = 1, g_p_
          g_diff_sih2(g_i_) = d3_b * g_y(g_i_, 2)
        enddo
        diff_sih2 = dble(3.45e-5) * d3_v
!--------
        diff_sih4 = 3.45e-5 * y(2) ** 1.5

        d2_v = p_g * 0.1d-3
        d4_v = d2_v * y(3)
        d2_b = 1.0d0 / mw_g1
        d3_b = d2_b * mw_g
        d4_b = d2_b * d4_v
        d6_b = d3_b * d2_v
        d7_b = d3_b * y(3) * 0.1d-3
        do g_i_ = 1, g_p_
          g_p_sih4(g_i_) = d4_b * g_mw_g(g_i_) + d6_b * g_y(g_i_, 3) + d&
     &7_b * g_p_g(g_i_)
        enddo
        p_sih4 = d4_v * mw_g / mw_g1
!--------
        d2_v = p_g * 0.1d-3
        d4_v = d2_v * y(5)
        d2_b = 1.0d0 / mw_g3
        d3_b = d2_b * mw_g
        d4_b = d2_b * d4_v
        d6_b = d3_b * d2_v
        d7_b = d3_b * y(5) * 0.1d-3
        do g_i_ = 1, g_p_
          g_p_h2(g_i_) = d4_b * g_mw_g(g_i_) + d6_b * g_y(g_i_, 5) + d7_&
     &b * g_p_g(g_i_)
        enddo
        p_h2 = d4_v * mw_g / mw_g3
!--------
        d2_v = p_g * 0.1d-03
        d4_v = d2_v * y(4)
        d2_b = 1.0d0 / mw_g2
        d3_b = d2_b * mw_g
        d4_b = d2_b * d4_v
        d6_b = d3_b * d2_v
        d7_b = d3_b * y(4) * 0.1d-03
        do g_i_ = 1, g_p_
          g_p_sih2(g_i_) = d4_b * g_mw_g(g_i_) + d6_b * g_y(g_i_, 4) + d&
     &7_b * g_p_g(g_i_)
        enddo
        p_sih2 = d4_v * mw_g / mw_g2
!--------

        cpsih4 = 14.13 / 32.12
        cpsih2 = 9.413 / 30.11
        d4_b = 1.0d0 / dble(2.02) * dble(0.00081)
        do g_i_ = 1, g_p_
          g_cph2(g_i_) = d4_b * g_y(g_i_, 2)
        enddo
        cph2 = (dble(6.62) + dble(0.00081) * y(2)) / dble(2.02)
!--------
        cpsi2h6 = 43.68 / 62.22
        d4_b = 1.0d0 / dble(28.01) * dble(0.001)
        do g_i_ = 1, g_p_
          g_cpn2(g_i_) = d4_b * g_y(g_i_, 2)
        enddo
        cpn2 = (dble(6.5) + dble(0.001) * y(2)) / dble(28.01)
!--------
        cpsi = 0.1442
        cpal2o3 = 0.1442
        do g_i_ = 1, g_p_
          g_c_pg(g_i_) = y(7) * g_cpn2(g_i_) + cpn2 * g_y(g_i_, 7) + cps&
     &i2h6 * g_y(g_i_, 6) + y(5) * g_cph2(g_i_) + cph2 * g_y(g_i_, 5) + &
     &cpsih2 * g_y(g_i_, 4) + cpsih4 * g_y(g_i_, 3)
        enddo
        c_pg = y(3) * cpsih4 + y(4) * cpsih2 + y(5) * cph2 + y(6) * cpsi&
     &2h6 + y(7) * cpn2
!--------

        if (y(8) .gt. 0.0) then
          do g_i_ = 1, g_p_
            g_c_ps(g_i_) = cpal2o3 * g_y(g_i_, 11) + cpsi * g_y(g_i_, 10&
     &)
          enddo
          c_ps = y(10) * cpsi + y(11) * cpal2o3
!--------
        else
          do g_i_ = 1, g_p_
            g_c_ps(g_i_) = 0.0d0
          enddo
          c_ps = cpal2o3
!--------
        endif

        d2_v = dble(-25600.) / y(2)
        d2_b = (-d2_v) / y(2)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_y(g_i_, 2)
        enddo
        d1_w = d2_v
        d3_v = (1.0d0 - y(8)) * 1.08d13
        d5_v = exp(d1_w)
        d1_p =  d5_v
        d6_v = d3_v * d5_v
        d10_v = y(1) * y(3) / mw_g1
        d4_b = d6_v * (1.0d0 / mw_g1)
        d5_b = d4_b * y(3)
        d6_b = d4_b * y(1)
        d9_b = d10_v * d3_v * d1_p
        d11_b = -(d10_v * d5_v * 1.08d13)
        do g_i_ = 1, g_p_
          g_raf(g_i_) = d6_b * g_y(g_i_, 3) + d5_b * g_y(g_i_, 1) + d9_b&
     & * g_d1_w(g_i_) + d11_b * g_y(g_i_, 8)
        enddo
        raf = d6_v * d10_v
!--------
        d2_v = dble(-2520.) / y(2)
        d2_b = (-d2_v) / y(2)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_y(g_i_, 2)
        enddo
        d1_w = d2_v
        d3_v = (1.0d0 - y(8)) * 1.26d12
        d5_v = exp(d1_w)
        d1_p =  d5_v
        d6_v = d3_v * d5_v
        d10_v = y(1) * y(4) / mw_g2
        d11_v = d6_v * d10_v
        d14_v = y(1) * y(5) / mw_g3
        d4_b = d11_v * (1.0d0 / mw_g3)
        d6_b = d4_b * y(1)
        d7_b = d14_v * d10_v
        d9_b = d14_v * d6_v * (1.0d0 / mw_g2)
        d5_b = d4_b * y(5) + d9_b * y(4)
        d10_b = d9_b * y(1)
        d13_b = d7_b * d3_v * d1_p
        d15_b = -(d7_b * d5_v * 1.26d12)
        do g_i_ = 1, g_p_
          g_rab(g_i_) = d6_b * g_y(g_i_, 5) + d10_b * g_y(g_i_, 4) + d5_&
     &b * g_y(g_i_, 1) + d13_b * g_d1_w(g_i_) + d15_b * g_y(g_i_, 8)
        enddo
        rab = d11_v * d14_v
!--------

        d2_v = dble(-2300.) / y(2)
        d2_b = (-d2_v) / y(2)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_y(g_i_, 2)
        enddo
        d1_w = d2_v
        d3_v = (1.0d0 - y(8)) * 1.0d14
        d5_v = exp(d1_w)
        d1_p =  d5_v
        d6_v = d3_v * d5_v
        d10_v = y(1) * y(3) / mw_g1
        d11_v = d6_v * d10_v
        d14_v = y(1) * y(4) / mw_g2
        d4_b = d11_v * (1.0d0 / mw_g2)
        d6_b = d4_b * y(1)
        d7_b = d14_v * d10_v
        d9_b = d14_v * d6_v * (1.0d0 / mw_g1)
        d5_b = d4_b * y(4) + d9_b * y(3)
        d10_b = d9_b * y(1)
        d13_b = d7_b * d3_v * d1_p
        d15_b = -(d7_b * d5_v * 1.0d14)
        do g_i_ = 1, g_p_
          g_rbf(g_i_) = d6_b * g_y(g_i_, 4) + d10_b * g_y(g_i_, 3) + d5_&
     &b * g_y(g_i_, 1) + d13_b * g_d1_w(g_i_) + d15_b * g_y(g_i_, 8)
        enddo
        rbf = d11_v * d14_v
!--------
        d2_v = dble(-26200.) / y(2)
        d2_b = (-d2_v) / y(2)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_y(g_i_, 2)
        enddo
        d1_w = d2_v
        d3_v = (1.0d0 - y(8)) * 5.62d15
        d5_v = exp(d1_w)
        d1_p =  d5_v
        d6_v = d3_v * d5_v
        d10_v = y(1) * y(6) / mw_g4
        d4_b = d6_v * (1.0d0 / mw_g4)
        d5_b = d4_b * y(6)
        d6_b = d4_b * y(1)
        d9_b = d10_v * d3_v * d1_p
        d11_b = -(d10_v * d5_v * 5.62d15)
        do g_i_ = 1, g_p_
          g_rbb(g_i_) = d6_b * g_y(g_i_, 6) + d5_b * g_y(g_i_, 1) + d9_b&
     & * g_d1_w(g_i_) + d11_b * g_y(g_i_, 8)
        enddo
        rbb = d6_v * d10_v
!--------
        d2_v = dble(-23016.0) / y(2)
        d2_b = (-d2_v) / y(2)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_y(g_i_, 2)
        enddo
        d1_w = d2_v
        d2_v = dble(3954.0) / y(2)
        d2_b = (-d2_v) / y(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d2_b * g_y(g_i_, 2)
        enddo
        d2_w = d2_v
        d4_v = 6.d0 * y(8) / d_p * 2.15d10
        d6_v = exp(d1_w)
        d2_p =  d6_v
        d7_v = d4_v * d6_v
        d11_v = y(1) * y(3) / mw_g1
        d17_v = exp(d2_w)
        d1_p =  d17_v
        d18_v = 7.6d-3 * d17_v
        d21_v = 1.0d0 + dble(0.034) * p_h2 + d18_v * p_sih4
        d22_v = d7_v * d11_v / d21_v
        d2_b = 1.0d0 / d21_v
        d3_b = (-d22_v) / d21_v
        d7_b = d3_b * d18_v
        d9_b = d3_b * p_sih4 * 7.6d-3 * d1_p
        d11_b = d3_b * dble(0.034)
        d12_b = d2_b * d11_v
        d14_b = d2_b * d7_v * (1.0d0 / mw_g1)
        d15_b = d14_b * y(3)
        d16_b = d14_b * y(1)
        d19_b = d12_b * d4_v * d2_p
        d22_b = d12_b * d6_v * 2.15d10 * (1.0d0 / d_p) * 6.d0
        do g_i_ = 1, g_p_
          g_rc(g_i_) = d7_b * g_p_sih4(g_i_) + d9_b * g_d2_w(g_i_) + d11&
     &_b * g_p_h2(g_i_) + d16_b * g_y(g_i_, 3) + d15_b * g_y(g_i_, 1) + &
     &d19_b * g_d1_w(g_i_) + d22_b * g_y(g_i_, 8)
        enddo
        rc = d22_v
!--------
        d2_v = diff_sih2 * n_sh
        d5_v = d2_v * p_sih2 * dble(6.0)
        d10_v = d_p * mw_g2 * (d_p * r_sih2 * y(2))
        d11_v = d5_v * y(8) / d10_v
        d2_b = 1.0d0 / d10_v
        d5_b = (-d11_v) / d10_v * (d_p * mw_g2) * (d_p * r_sih2)
        d7_b = d2_b * d5_v
        d8_b = d2_b * y(8) * dble(6.0)
        d10_b = d8_b * d2_v
        d11_b = d8_b * p_sih2 * n_sh
        do g_i_ = 1, g_p_
          g_rd(g_i_) = d5_b * g_y(g_i_, 2) + d7_b * g_y(g_i_, 8) + d10_b&
     & * g_p_sih2(g_i_) + d11_b * g_diff_sih2(g_i_)
        enddo
        rd = d11_v
!--------

        do g_i_ = 1, g_p_
          g_rg(g_i_, 1) = (-mw_g1) * g_rc(g_i_) + mw_g1 * g_rbb(g_i_) + &
     &(-mw_g1) * g_rbf(g_i_) + mw_g1 * g_rab(g_i_) + (-mw_g1) * g_raf(g_&
     &i_)
        enddo
        rg(1) = (-raf + rab - rbf + rbb - rc) * mw_g1
!--------
        do g_i_ = 1, g_p_
          g_rg(g_i_, 2) = (-mw_g2) * g_rd(g_i_) + mw_g2 * g_rbb(g_i_) + &
     &(-mw_g2) * g_rbf(g_i_) + (-mw_g2) * g_rab(g_i_) + mw_g2 * g_raf(g_&
     &i_)
        enddo
        rg(2) = (raf - rab - rbf + rbb - rd) * mw_g2
!--------
        d7_b = mw_g3 * dble(2.0)
        do g_i_ = 1, g_p_
          g_rg(g_i_, 3) = mw_g3 * g_rd(g_i_) + d7_b * g_rc(g_i_) + (-mw_&
     &g3) * g_rab(g_i_) + mw_g3 * g_raf(g_i_)
        enddo
        rg(3) = (raf - rab + dble(2.0) * rc + rd) * mw_g3
!--------
        do g_i_ = 1, g_p_
          g_rg(g_i_, 4) = (-mw_g4) * g_rbb(g_i_) + mw_g4 * g_rbf(g_i_)
        enddo
        rg(4) = (rbf - rbb) * mw_g4
!--------
        do g_i_ = 1, g_p_
          g_rg(g_i_, 5) = 0.0d0
        enddo
        rg(5) = 0.0d0
!--------
        do g_i_ = 1, g_p_
          g_rs(g_i_) = mw_s1 * g_rd(g_i_) + mw_s1 * g_rc(g_i_)
        enddo
        rs = (rc + rd) * mw_s1
!--------
        rs1 = 0.0d0

        d2_v = 1.0d0 - y(8)
        d3_v = 1.0d0 / d2_v
        d12_v = y(1) / ro_s
        d15_v = rg(1) + rg(2) + rg(3) + rg(4) + d12_v * rs
        d7_b = d3_v * d12_v
        d8_b = d3_v * rs * (1.0d0 / ro_s)
        d16_b = -(d15_v * ((-d3_v) / d2_v))
        do g_i_ = 1, g_p_
          g_yt(g_i_, 1) = d7_b * g_rs(g_i_) + d8_b * g_y(g_i_, 1) + d3_v&
     & * g_rg(g_i_, 4) + d3_v * g_rg(g_i_, 3) + d3_v * g_rg(g_i_, 2) + d&
     &3_v * g_rg(g_i_, 1) + d16_b * g_y(g_i_, 8)
        enddo
        yt(1) = d3_v * d15_v
!--------
        d12_v = (-(dble(57060.0) * (raf - rab) - dble(54350.0) * (rbf - &
     &rbb))) / y(1)
        d14_v = 1.0d0 - y(8)
        d15_v = d12_v / d14_v
        d17_v = d15_v / c_pg
        d2_b = 1.0d0 / c_pg
        d3_b = (-d17_v) / c_pg
        d4_b = d2_b * (1.0d0 / d14_v)
        d6_b = -(d2_b * ((-d15_v) / d14_v))
        d8_b = d4_b * ((-d12_v) / y(1))
        d9_b = -(d4_b * (1.0d0 / y(1)))
        d12_b = (-d9_b) * dble(54350.0)
        d15_b = d9_b * dble(57060.0)
        do g_i_ = 1, g_p_
          g_yt(g_i_, 2) = d3_b * g_c_pg(g_i_) + d6_b * g_y(g_i_, 8) + d8&
     &_b * g_y(g_i_, 1) + (-d12_b) * g_rbb(g_i_) + d12_b * g_rbf(g_i_) +&
     & (-d15_b) * g_rab(g_i_) + d15_b * g_raf(g_i_)
        enddo
        yt(2) = d17_v
!--------
        d3_v = 1.0d0 - y(8)
        d5_v = d3_v * y(1)
        d6_v = rg(1) / d5_v
        d8_v = 1.0d0 - y(8)
        d9_v = d8_v * y(1)
        d10_v = y(3) / d9_v
        d16_v = rg(1) + rg(2) + rg(3) + rg(4)
        d12_b = (-d16_v) * (1.0d0 / d9_v)
        d13_b = (-d16_v) * ((-d10_v) / d9_v)
        d10_b = -d10_v + 1.0d0 / d5_v
        d17_b = (-d6_v) / d5_v
        d15_b = d13_b * d8_v + d17_b * d3_v
        d16_b = -(d13_b * y(1)) + (-(d17_b * y(1)))
        do g_i_ = 1, g_p_
          g_yt(g_i_, 3) = (-d10_v) * g_rg(g_i_, 4) + (-d10_v) * g_rg(g_i&
     &_, 3) + (-d10_v) * g_rg(g_i_, 2) + d12_b * g_y(g_i_, 3) + d15_b * &
     &g_y(g_i_, 1) + d16_b * g_y(g_i_, 8) + d10_b * g_rg(g_i_, 1)
        enddo
        yt(3) = d6_v - d10_v * d16_v
!--------
        d3_v = 1.0d0 - y(8)
        d5_v = d3_v * y(1)
        d6_v = rg(2) / d5_v
        d8_v = 1.0d0 - y(8)
        d9_v = d8_v * y(1)
        d10_v = y(4) / d9_v
        d16_v = rg(1) + rg(2) + rg(3) + rg(4)
        d12_b = (-d16_v) * (1.0d0 / d9_v)
        d13_b = (-d16_v) * ((-d10_v) / d9_v)
        d11_b = -d10_v + 1.0d0 / d5_v
        d17_b = (-d6_v) / d5_v
        d15_b = d13_b * d8_v + d17_b * d3_v
        d16_b = -(d13_b * y(1)) + (-(d17_b * y(1)))
        do g_i_ = 1, g_p_
          g_yt(g_i_, 4) = (-d10_v) * g_rg(g_i_, 4) + (-d10_v) * g_rg(g_i&
     &_, 3) + (-d10_v) * g_rg(g_i_, 1) + d12_b * g_y(g_i_, 4) + d15_b * &
     &g_y(g_i_, 1) + d16_b * g_y(g_i_, 8) + d11_b * g_rg(g_i_, 2)
        enddo
        yt(4) = d6_v - d10_v * d16_v
!--------
        d3_v = 1.0d0 - y(8)
        d5_v = d3_v * y(1)
        d6_v = rg(3) / d5_v
        d8_v = 1.0d0 - y(8)
        d9_v = d8_v * y(1)
        d10_v = y(5) / d9_v
        d16_v = rg(1) + rg(2) + rg(3) + rg(4)
        d12_b = (-d16_v) * (1.0d0 / d9_v)
        d13_b = (-d16_v) * ((-d10_v) / d9_v)
        d9_b = -d10_v + 1.0d0 / d5_v
        d17_b = (-d6_v) / d5_v
        d15_b = d13_b * d8_v + d17_b * d3_v
        d16_b = -(d13_b * y(1)) + (-(d17_b * y(1)))
        do g_i_ = 1, g_p_
          g_yt(g_i_, 5) = (-d10_v) * g_rg(g_i_, 4) + (-d10_v) * g_rg(g_i&
     &_, 2) + (-d10_v) * g_rg(g_i_, 1) + d12_b * g_y(g_i_, 5) + d15_b * &
     &g_y(g_i_, 1) + d16_b * g_y(g_i_, 8) + d9_b * g_rg(g_i_, 3)
        enddo
        yt(5) = d6_v - d10_v * d16_v
!--------
        d3_v = 1.0d0 - y(8)
        d5_v = d3_v * y(1)
        d6_v = rg(4) / d5_v
        d8_v = 1.0d0 - y(8)
        d9_v = d8_v * y(1)
        d10_v = y(6) / d9_v
        d16_v = rg(1) + rg(2) + rg(3) + rg(4)
        d12_b = (-d16_v) * (1.0d0 / d9_v)
        d13_b = (-d16_v) * ((-d10_v) / d9_v)
        d7_b = -d10_v + 1.0d0 / d5_v
        d17_b = (-d6_v) / d5_v
        d15_b = d13_b * d8_v + d17_b * d3_v
        d16_b = -(d13_b * y(1)) + (-(d17_b * y(1)))
        do g_i_ = 1, g_p_
          g_yt(g_i_, 6) = (-d10_v) * g_rg(g_i_, 3) + (-d10_v) * g_rg(g_i&
     &_, 2) + (-d10_v) * g_rg(g_i_, 1) + d12_b * g_y(g_i_, 6) + d15_b * &
     &g_y(g_i_, 1) + d16_b * g_y(g_i_, 8) + d7_b * g_rg(g_i_, 4)
        enddo
        yt(6) = d6_v - d10_v * d16_v
!--------
        d3_v = 1.0d0 - y(8)
        d5_v = d3_v * y(1)
        d6_v = rg(5) / d5_v
        d8_v = 1.0d0 - y(8)
        d9_v = d8_v * y(1)
        d10_v = y(7) / d9_v
        d17_v = rg(1) + rg(2) + rg(3) + rg(4)
        d12_b = (-d17_v) * (1.0d0 / d9_v)
        d13_b = (-d17_v) * ((-d10_v) / d9_v)
        d17_b = 1.0d0 / d5_v
        d18_b = (-d6_v) / d5_v
        d15_b = d13_b * d8_v + d18_b * d3_v
        d16_b = -(d13_b * y(1)) + (-(d18_b * y(1)))
        do g_i_ = 1, g_p_
          g_yt(g_i_, 7) = (-d10_v) * g_rg(g_i_, 4) + (-d10_v) * g_rg(g_i&
     &_, 3) + (-d10_v) * g_rg(g_i_, 2) + (-d10_v) * g_rg(g_i_, 1) + d12_&
     &b * g_y(g_i_, 7) + d15_b * g_y(g_i_, 1) + d16_b * g_y(g_i_, 8) + d&
     &17_b * g_rg(g_i_, 5)
        enddo
        yt(7) = d6_v - d10_v * d17_v
!--------
        if (y(8) .ge. 1.0d-8) then
          d2_b = 1.0d0 / ro_s
          do g_i_ = 1, g_p_
            g_yt(g_i_, 8) = d2_b * g_rs(g_i_)
          enddo
          yt(8) = rs / ro_s
!--------
          d9_v = (-(dble(-8210.0) * rc - dble(65360.0) * rd)) / ro_s / c_ps
          d11_v = d9_v / y(8)
          d2_b = 1.0d0 / y(8)
          d3_b = (-d11_v) / y(8)
          d5_b = d2_b * ((-d9_v) / c_ps)
          d7_b = -(d2_b * (1.0d0 / c_ps) * (1.0d0 / ro_s))
          d10_b = (-d7_b) * dble(65360.0)
          d11_b = d7_b * dble(-8210.0)
          do g_i_ = 1, g_p_
            g_yt(g_i_, 9) = d3_b * g_y(g_i_, 8) + d5_b * g_c_ps(g_i_) + &
     &d10_b * g_rd(g_i_) + d11_b * g_rc(g_i_)
          enddo
          yt(9) = d11_v
!--------
          d3_v = y(8) * ro_s
          d4_v = rs / d3_v
          d6_v = y(8) * ro_s
          d7_v = y(10) / d6_v
          d6_b = (-rs) * (1.0d0 / d6_v)
          d5_b = -d7_v + 1.0d0 / d3_v
          d8_b = (-rs) * ((-d7_v) / d6_v) * ro_s + (-d4_v) / d3_v * ro_s
          do g_i_ = 1, g_p_
            g_yt(g_i_, 10) = d6_b * g_y(g_i_, 10) + d8_b * g_y(g_i_, 8) &
     &+ d5_b * g_rs(g_i_)
          enddo
          yt(10) = d4_v - d7_v * rs
!--------
          d2_v = y(8) * ro_s
          d3_v = rs1 / d2_v
          d5_v = y(8) * ro_s
          d6_v = y(11) / d5_v
          d6_b = (-rs) * (1.0d0 / d5_v)
          d8_b = (-rs) * ((-d6_v) / d5_v) * ro_s + (-d3_v) / d2_v * ro_s
          do g_i_ = 1, g_p_
            g_yt(g_i_, 11) = (-d6_v) * g_rs(g_i_) + d6_b * g_y(g_i_, 11)&
     & + d8_b * g_y(g_i_, 8)
          enddo
          yt(11) = d3_v - d6_v * rs
!--------
        else
          do g_i_ = 1, g_p_
            g_yt(g_i_, 8) = 0.0d0
          enddo
          yt(8) = 0.0d0
!--------
          do g_i_ = 1, g_p_
            g_yt(g_i_, 9) = 0.0d0
          enddo
          yt(9) = 0.0d0
!--------
          do g_i_ = 1, g_p_
            g_yt(g_i_, 10) = 0.0d0
          enddo
          yt(10) = 0.0d0
!--------
          do g_i_ = 1, g_p_
            g_yt(g_i_, 11) = 0.0d0
          enddo
          yt(11) = 0.0d0
!--------
        endif
        do g_i_ = 1, g_p_
          g_yt(g_i_, 12) = 0.0d0
        enddo
        yt(12) = 0.0d0
!--------

        return
      end
