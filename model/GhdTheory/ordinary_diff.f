!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  subroutine name: ordinary_diff(s,mi,sigmai,alpha,phii,T,phi,ni,n,rhoi,rho,
!             m,mu,sigma,chi,zeta0,theta,Ti,DT,nu,Dij,I_ilj,dTl_dnj,
!             dzeta0_dnj,dchi0il_dnj)
!
!  author:  C. Hrenya, Feb 2009
!
!  Purpose: find ordinary diffusion coefficient according to GHD polydisperse KT
!
!  Literature/References:  
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!                                                         
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine ordinary_diff(s,mi,sigmai,alpha,phii,T,phi,ni,n,rhoi, &
          rho,m,mu,sigma,chi,zeta0,theta,Ti,DT,nu,Dij,I_ilj,dTl_dnj, &
          dzeta0_dnj,dchi0il_dnj)
      Implicit NONE

      integer s 

      double precision mi(s),sigmai(s),alpha(s,s),phii(s),T,phi, &
                      ni(s),n,rhoi(s),rho,m,mu(s,s),sigma(s,s), &
                      chi(s,s),zeta0, &
                      theta(s),Ti(s),DT(s),nu(s,s),Dij(s,s), &
                      I_ilj(s,s,s),dTl_dnj(s,s),dzeta0_dnj(s), &
                      dchi0il_dnj(s,s,s)

      integer i,j,k,l,kk
      double precision kronecker(s,s),M1,M2,M3,dphi_dnl(s), &
                      dM1_dnl(s),dM2_dnl(s),dM3_dnl(s), &
                      dni_dnl(s,s), &
                      dmuioverT_dnl(s,s),perturbation,delta_ni(s), &
                      delta_phii(s),delta_rhoi(s),chi_ij, &
                      theta_pos(s),Ti_pos(s), &
                      zeta0_pos,group1(s,s),group2(s,s), &
                      p_pos,niTi_pos(s),theta_neg(s), &
                      Ti_neg(s),zeta0_neg,p_neg,niTi_neg(s), &
                      dp_dnj(s),dniTi_dnj(s,s),sum1(s,s), &
                      Amat(s,s),bmat(s,s), &
                      Amat0(s,s),bmat0(s)
      double precision indx(s),d, twoDeltaNi

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve

      double precision pi
      parameter (pi=3.14159265458979323846d0)

      do i=1,s
         do j=1,s
            if (i.eq.j) then
               kronecker(i,j) = 1.d0
            else
               kronecker(i,j) = 0.d0
            endif
         enddo
      enddo

!calculate Mn's - p 6.1 CMH notes
      M1=0.d0
      M2=0.d0
      M3=0.d0
      do i=1,s
         M1 = M1 + ni(i)*sigmai(i) 
         M2 = M2 + ni(i)*sigmai(i)**2 
         M3 = M3 + ni(i)*sigmai(i)**3 
      enddo
      M1 = M1/n 
      M2 = M2/n 
      M3 = M3/n 

!calculate partial derivative of phi and Mn's wrt n_l - p 6.2 CMH notes
      do l=1,s
         dphi_dnl(l) = pi/6.d0*sigmai(l)**3 
         dM1_dnl(l)  = (sigmai(l)-M1)/n 
         dM2_dnl(l)  = (sigmai(l)**2-M2)/n 
         dM3_dnl(l)  = (sigmai(l)**3-M3)/n 
         do i=1,s
            dni_dnl(i,l)  = kronecker(i,l) 
         enddo
      enddo

!calculate partial derivative of mu_i/T wrt n_l - p 6.4 CMH notes;
!this Matlab code was generated by Mathematica, where the symbolic derivative was evaluated
! note that 1/ni term will cancel with -1/ni defined in I_ilj
      do i=1,s
         do l=1,s
            dmuioverT_dnl(i,l) = &
             (3.d0*sigmai(i)*phi*dM2_dnl(l))/(M3*(1.d0 - phi)) -  &
       (3.d0*sigmai(i)*M2*phi*dM3_dnl(l))/(M3**2*(1.d0 - phi)) +  &
       dphi_dnl(l)/(1.d0 - phi) + (3.d0*sigmai(i)*M2* &
       dphi_dnl(l))/(M3*(1.d0 - phi)) +  &
       (3.d0*sigmai(i)*M2*phi*dphi_dnl(l))/(M3*(1.d0 - phi)**2)+  &
       3.d0*sigmai(i)**2*((phi*dM1_dnl(l))/(M3*(1.d0 - phi)) +  &
          (2.d0*DLog(1.d0 - phi)*M2*dM2_dnl(l))/M3**2 +  &
          (2.d0*M2*phi*dM2_dnl(l))/(M3**2*(1.d0 - phi)**2) -  &
          (2.d0*DLog(1.d0 - phi)*M2**2*dM3_dnl(l))/M3**3 -  &
          (2.d0*M2**2*phi*dM3_dnl(l))/(M3**3*(1.d0 - phi)**2) -  &
          (M1*phi*dM3_dnl(l))/(M3**2*(1.d0 - phi)) +  &
          (M2**2*dphi_dnl(l))/(M3**2*(1.d0 - phi)**2) -  &
          (M2**2*dphi_dnl(l))/(M3**2*(1.d0 - phi)) +  &
          (M1*dphi_dnl(l))/(M3*(1.d0 - phi)) +  &
          (2.d0*M2**2*phi*dphi_dnl(l))/(M3**2*(1.d0 - phi)**3) +  &
          (M1*phi*dphi_dnl(l))/(M3*(1.d0 - phi)**2)) +  &
       sigmai(i)**3*((3.d0*M2*phi**2*dM1_dnl(l))/ &
          (M3**2*(1.d0 - phi)**2) -  &
          (6.d0*DLog(1.d0 - phi)*M2**2*dM2_dnl(l))/M3**3 +  &
          (3.d0*M1*phi**2*dM2_dnl(l))/(M3**2*(1.d0 - phi)**2) -  &
          (3.d0*M2**2*phi*(2.d0 - 5.d0*phi + phi**2)* &
          dM2_dnl(l))/(M3**3*(1.d0 - phi)**3) +  &
          (6.d0*DLog(1.d0 - phi)*M2**3*dM3_dnl(l))/M3**4 -  &
          (phi*dM3_dnl(l))/(M3**2*(1.d0 - phi)) -  &
          (6.d0*M1*M2*phi**2*dM3_dnl(l))/ &
          (M3**3*(1.d0 - phi)**2) +  &
          (3.d0*M2**3*phi*(2.d0 - 5.d0*phi + phi**2)* &
          dM3_dnl(l))/(M3**4*(1.d0 - phi)**3) +  &
          (2.d0*M2**3*dphi_dnl(l))/(M3**3*(1.d0 - phi)) +  &
           dphi_dnl(l)/(M3*(1.d0 - phi)) +  &
          (6.d0*M1*M2*phi*dphi_dnl(l))/(M3**2*(1.d0 - phi)**2) +  &
          (phi*dphi_dnl(l))/(M3*(1.d0 - phi)**2) +  &
          (6.d0*M1*M2*phi**2*dphi_dnl(l))/ &
          (M3**2*(1.d0 - phi)**3) -  &
          (M2**3*(2.d0 - 5.d0*phi + phi**2)*dphi_dnl(l))/ &
          (M3**3*(1.d0 - phi)**3) -  &
          (3.d0*M2**3*phi*(2.d0 - 5.d0*phi +  &
          phi**2)*dphi_dnl(l))/(M3**3*(1.d0 - phi)**4) -  &
          (M2**3*phi*(-5.d0*dphi_dnl(l) +  &
          2.d0*phi*dphi_dnl(l)))/(M3**3*(1.d0 - phi)**3))
         enddo
      enddo

!calculate partial derivative of pair correlation function wrt n_j - p 6.5 CMH notes;
!this Matlab code was generated by Mathematica, where the symbolic derivative was evaluated
      do i=1,s
         do l=1,s
            do j=1,s
              dchi0il_dnj(i,l,j) = &
             (1.5d0*sigmai(i)*sigmai(l)*phi*dM2_dnl(j))/ &
       (sigma(i,l)*M3*(1.d0 - phi)**2) +  &
       (1.d0*sigmai(i)**2*sigmai(l)**2*M2*phi**2*dM2_dnl(j))/ &
        (sigma(i,l)**2*M3**2*(1.d0 - phi)**3) -  &
       (1.5d0*sigmai(i)*sigmai(l)*M2*phi*dM3_dnl(j))/ &
       (sigma(i,l)*M3**2*(1.d0 - phi)**2) -  &
       (1.d0*sigmai(i)**2*sigmai(l)**2*M2**2*phi**2* &
       dM3_dnl(j))/(sigma(i,l)**2*M3**3*(1.d0 - phi)**3) +  &
       dphi_dnl(j)/(1.d0 - phi)**2 + (1.5d0*sigmai(i)* &
       sigmai(l)*M2*dphi_dnl(j))/ &
        (sigma(i,l)*M3*(1.d0 - phi)**2) + (1.d0*sigmai(i)**2* &
       sigmai(l)**2*M2**2*phi*dphi_dnl(j))/ &
        (sigma(i,l)**2*M3**2*(1.d0 - phi)**3) +  &
       (3.d0*sigmai(i)*sigmai(l)*M2*phi*dphi_dnl(j))/ &
       (sigma(i,l)*M3*(1.d0 - phi)**3) + (1.5d0*sigmai(i)**2* &
       sigmai(l)**2*M2**2*phi**2*dphi_dnl(j))/ &
        (sigma(i,l)**2*M3**2*(1.d0 - phi)**4)
            enddo
         enddo
      enddo

!determination of I_ilj - p 6.1 of CMH notes
!      do i=1,s
!         do l=1,s
!            do j=1,s
!               I_ilj(i,l,j) = (1.5d0/pi*(kronecker(j,l)*ni(l)* &
!                 dmuioverT_dnl(i,l)  &
!                 -ni(l)/ni(i)*kronecker(j,l)*kronecker(i,j)  &
!                 -4.d0*pi/3.d0*kronecker(j,l)*ni(l)*chi(i,j)* &
!                 sigma(i,j)**3)-ni(l)**2*sigma(i,l)**3* &
!                 dchi0il_dnj(i,l,j)) 
!            enddo
!         enddo
!      enddo
!
!CMH correction (10/9/09) to previous error in notes - this I_ilj for binary mix only
      I_ilj(1,2,1) = ( 3.d0/2.d0/pi*(dmuioverT_dnl(1,1)-4.d0*pi/3.d0*chi(1,1)*sigma(1,1)**3) &
          -ni(1)*sigma(1,1)**3*dchi0il_dnj(1,1,1) &
          -ni(2)*sigma(1,2)**3*dchi0il_dnj(1,2,1) ) &
          /(chi(1,2)*sigma(1,2)**3)
     
      I_ilj(1,2,2) = ( 3.d0/2.d0/pi*(dmuioverT_dnl(1,2)- 4.d0*pi/3.d0*chi(1,2)*sigma(1,2)**3) &
          -ni(1)*chi(1,1)*sigma(1,1)**3/chi(1,1)*dchi0il_dnj(1,1,2) &
          -ni(2)*chi(1,2)*sigma(1,2)**3/chi(1,2)*dchi0il_dnj(1,2,2) ) &
          /(chi(1,2)*sigma(1,2)**3)
     
      I_ilj(2,1,1) = ( 3.d0/2.d0/pi*(dmuioverT_dnl(2,1)- 4.d0*pi/3.d0*chi(2,1)*sigma(2,1)**3) &
         -ni(2)*chi(2,2)*sigma(2,2)**3/chi(2,2)*dchi0il_dnj(2,2,1) & 
         -ni(1)*chi(2,1)*sigma(2,1)**3/chi(2,1)*dchi0il_dnj(2,1,1) ) &
          /(chi(2,1)*sigma(2,1)**3)
     
      I_ilj(2,1,2) = ( 3.d0/2.d0/pi*(dmuioverT_dnl(2,2)-4.d0*pi/3.d0*chi(2,2)*sigma(2,2)**3) &
          -ni(2)*chi(2,2)*sigma(2,2)**3/chi(2,2)*dchi0il_dnj(2,2,2) &
          -ni(1)*chi(2,1)*sigma(2,1)**3/chi(2,1)*dchi0il_dnj(2,1,2) ) &
          /(chi(2,1)*sigma(2,1)**3)
      
      I_ilj(1,1,1) = 0.d0
      I_ilj(1,1,2) = 0.d0
      I_ilj(2,2,1) = 0.d0
      I_ilj(2,2,2) = 0.d0

!Numerical evaulation of several partial derivatives - p 6.6 of CMH notes.
!Each of the partial derivatives is with respect to nj, where T and other
!ni are constant, so need to 
!   1) perturb nj a small positive amount (and correspondingly n,phii,phi)
!   2) evaulate dependent quantity of interest
!   3) perturb nj a small negative amount (and correspondingly n,phii,phi)
!   4) evaluate dependent quantity of interest
!   5) calculate partial derivative
!   6) return nj & n to original values

      perturbation = 0.000001d0    ! small perturb of current value in each direction
      do j=1,s
         !calculate perturbation amounts
	    if(ni(j)/=0d0) then
              delta_ni(j)  = perturbation*ni(j) 
              delta_phii(j)= perturbation*ni(j)*pi/6.d0*sigmai(j)**3 
              delta_rhoi(j)   = perturbation*ni(j)*mi(j) !;
	      twoDeltaNi        = 2d0*delta_ni(j)
	    else
	      delta_ni(j)   = perturbation
              delta_phii(j) = perturbation*pi/6.d0*sigmai(j)**3 
              delta_rhoi(j) = perturbation*mi(j)
	      twoDeltaNi      = delta_ni(j)
	    endif
         !perturb current nj values in + direction
            ni(j)   = ni(j) + delta_ni(j) 
            n       = n     + delta_ni(j) 
            phii(j) = phii(j)+ delta_phii(j) 
            phi     = phi    + delta_phii(j) 
            rhoi(j) = rhoi(j) + delta_rhoi(j) !;
            do i=1,s
               do k=1,s
                  call chi_ij_GHD(s,i,k,sigmai,phi,ni,chi_ij)
                  chi(i,k) = chi_ij
               enddo
            enddo
         !evaluate dependent quantites at increased value of nj et al.
            call cooling_rate(s,mi,ni,n,m,T,chi,sigmai,alpha,rhoi, &
                             theta_pos)
            do i=1,s
               Ti_pos(i) = mi(i)*T/(m*theta_pos(i))
            enddo
            do l=1,s                                   
               group1(1,l)  = chi(1,l)*ni(l)*mu(l,1)* &
                             sigma(1,l)**2*(1d0+alpha(1,l))
               group2(1,l)  = mu(l,1)/2.d0*(1d0+alpha(1,l))
            enddo
            zeta0_pos = 0d0 !;
            do k=1,s
               zeta0_pos = zeta0_pos + group1(1,k)* &
                       dsqrt((theta_pos(1)+theta_pos(k)) &
                       /(theta_pos(1)*theta_pos(k))) &
                       * ((1.d0-group2(1,k)*(theta_pos(1)+ &
                       theta_pos(k))/theta_pos(k)))
            enddo
            zeta0_pos = 8.d0/3.d0*dsqrt(2.d0*pi*T/m)*zeta0_pos            
            call pressure (s,alpha,ni,n,mu,sigma,chi,T,Ti_pos, &
                          p_pos)
            do i=1,s
               niTi_pos(i) = ni(i)*Ti_pos(i) 
            enddo
         !perturb current nj values in - direction
	    if(ni(j) /= delta_ni(j)) then
              ni(j)   = ni(j) - 2.d0*delta_ni(j)      
              n       = n     - 2.d0*delta_ni(j) 
              phii(j) = phii(j) - 2.d0*delta_phii(j) 
              phi     = phi     - 2.d0*delta_phii(j) 
              rhoi(j) = rhoi(j) - 2.d0*delta_rhoi(j) !;
	    else ! bring ni back to zero and do a one-sided (1st order) derivative
              ni(j)   = ni(j) - delta_ni(j)      
              n       = n     - delta_ni(j) 
              phii(j) = phii(j) - delta_phii(j) 
              phi     = phi     - delta_phii(j) 
              rhoi(j) = rhoi(j) - delta_rhoi(j)
	    endif
            do i=1,s
               do k=1,s
                  call chi_ij_GHD(s,i,k,sigmai,phi,ni,chi_ij)
                  chi(i,k) = chi_ij
               enddo
            enddo
         !evaluate dependent quantites at decreased value of nj et al.
            call cooling_rate(s,mi,ni,n,m,T,chi,sigmai,alpha,rhoi, &
                             theta_neg)
            do i=1,s
               Ti_neg(i) = mi(i)*T/(m*theta_neg(i))
            enddo
            do l=1,s                                   
               group1(1,l)  = chi(1,l)*ni(l)*mu(l,1)* &
                             sigma(1,l)**2*(1d0+alpha(1,l))
               group2(1,l)  = mu(l,1)/2.d0*(1d0+alpha(1,l))
            enddo
            zeta0_neg = 0d0 !;
            do k=1,s
               zeta0_neg = zeta0_neg + group1(1,k)* &
                       dsqrt((theta_neg(1)+theta_neg(k)) &
                       /(theta_neg(1)*theta_neg(k))) &
                       * ((1.d0-group2(1,k)*(theta_neg(1)+ &
                       theta_neg(k))/theta_neg(k)))
            enddo
            zeta0_neg = 8.d0/3.d0*dsqrt(2.d0*pi*T/m)*zeta0_neg            
            call pressure (s,alpha,ni,n,mu,sigma,chi,T,Ti_neg, &
                          p_neg)
            do i=1,s
               niTi_neg(i) = ni(i)*Ti_neg(i) 
            enddo
         !find partial derivatives
            dzeta0_dnj(j) = (zeta0_pos - zeta0_neg)/ &
                           ( twoDeltaNi) 
            dp_dnj(j) = (p_pos - p_neg)/( twoDeltaNi) 
            do l=1,s
               dTl_dnj(l,j) = (Ti_pos(l) - Ti_neg(l))/ &
                             ( twoDeltaNi)  
            enddo
            do i=1,s
               dniTi_dnj(i,j) = (niTi_pos(i) - niTi_neg(i))/ &
                               ( twoDeltaNi)  
            enddo
         !reset nj et al. to original values
	    if(ni(j) /= 0d0) then
              ni(j)   = ni(j) + delta_ni(j)      
              n       = n     + delta_ni(j) 
              phii(j) = phii(j) + delta_phii(j) 
              phi     = phi     + delta_phii(j) 
              rhoi(j) = rhoi(j) + delta_rhoi(j) !;
	    endif
            do i=1,s
               do k=1,s
                  call chi_ij_GHD(s,i,k,sigmai,phi,ni,chi_ij)
                  chi(i,k) = chi_ij
               enddo
            enddo
      enddo

      do i=1,s
         do j=1,s
            sum1(i,j) = 0.d0    !calculate summation used in b vector - p 6 CMH notes
         enddo
      enddo
      do i=1,s
         do j=1,s
            do l=1,s
               sum1(i,j) = sum1(i,j) + chi(i,l)*sigma(i,l)**3* &
                 mu(l,i)*(1.d0+alpha(i,l))  &
                 *((Ti(i)/mi(i)+Ti(l)/mi(l))*(kronecker(j,l)  &
                 +0.5d0*(ni(l)/chi(i,l)* &
                 dchi0il_dnj(i,l,j)+I_ilj(i,l,j)))  &
                 +ni(l)/mi(l)*dTl_dnj(l,j)) 
            enddo
         enddo
      enddo

      do i=1,s
         do j=1,s
            Amat(i,j) = (nu(i,j)-0.5d0*zeta0*kronecker(i,j))*  &        !A matrix for solution of Dij (p 7 CMH notes)
                        mi(j)/mi(i)   
            bmat(i,j) = rho**2/mi(i)/mi(j)*dzeta0_dnj(j)*DT(i)  &       !b matrix for solution of Dij (p 7 CMH notes)
              + rho/mi(i)/mi(j)*dniTi_dnj(i,j)  &
              - ni(i)/mi(j)*dp_dnj(j)  &
              + 2.d0*pi/3.d0*rho*ni(i)/mi(j)*sum1(i,j) 
         enddo
      enddo

! this extra kk loop and addition of Amat0 and bmat0 is necessary 
! since x & b in Ax=b are s by s matrices rather than vectors of length s,
! whereas LUBSKB is specific to x & b vectors of length s
      do kk=1,s 
         do i=1,s
            do j=1,s
                Amat0(i,j) = Amat(i,j)
            enddo
         enddo
         do i=1,s
            bmat0(i) = bmat(i,kk)
         enddo

         CALL LUDCMP(Amat0, s, NP, indx, d) ! solve system of s linear equations using
         CALL LUBKSB(Amat0, s, NP, indx, bmat0) ! LU decomposition

         do i=1,s
            Dij(i,kk) = bmat0(i)
         enddo
      enddo
        
      if(s==2) then ! for binary only
        Dij(1,1) = -mi(2)/mi(1)*Dij(2,1)
        Dij(1,2) = -mi(2)/mi(1)*Dij(2,2)
      endif

      return
      end subroutine ordinary_diff

