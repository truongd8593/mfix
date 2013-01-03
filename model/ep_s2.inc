!QX RO_S changed form 1d to 2
!      EP_s(IJK, xxM) = ROP_s(IJK, xxM)/RO_SV(xxM)
!      if(RO_SV(IJK,xxM) .gt. 0.d0) then
      EP_s(IJK, xxM) = ROP_s(IJK, xxM)/RO_SV(IJK,xxM)
!      endif
!end
