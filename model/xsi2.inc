!-----------------------------------------------
!     Special inline for xsi
!     should be similar to
!
!      xsi(v,dw) = merge( v >= 0, dwf, one-dwf)
!
!      slight difference when v is exactly zero but
!      may not be significant
!-----------------------------------------------
      xsi(XXXv,XXXdwf) = (sign(one,(-XXXv))+one)/(one+one) + &
                         sign(one,XXXv)*XXXdwf
