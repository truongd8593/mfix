!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                      
!  subroutine name: thermal_mobility(s,mi,alpha,ni,mu,sigma,chi,zeta0,
!                                theta,Ti,DF,gammaij,omega,Lij)
!
!  author:  C. Hrenya, Feb 2009
!
!  Purpose: find thermal mobility coefficient according to GHD polydisperse KT
!
!  Literature/References:  
!     C. Hrenya handwritten notes & Garzo, Hrenya, Dufty papers (PRE, 2007)
!                                                         
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine thermal_mobility(s,mi,m,sigmai,v0,alpha,ni,n,T,phii,phi,rho,rhoi,mu,sigma,chi,zeta0, &
                                 theta,Ti,DF,gammaij,omega,dTl_dnj,dzeta0_dnj,domegaij_dnj,domegaii_dni,Lij)
      Implicit NONE

      integer s 

      double precision mi(s),alpha(s,s),ni(s),mu(s,s),sigma(s,s), &
                      chi(s,s),zeta0,theta(s),Ti(s),DF(s,s), &
                      gammaij(s,s),omega(s,s),phii(s),rhoi(s), &
                      n,phi,rho,chi_ij,m, T,theta_pos(s), &
                      Ti_pos(s),group1(s,s),group2(s,s),zeta0_pos, &
                      p_pos, DT_pos(s),nu_pos(s,s),DF_pos(s,s), theta_neg(s), &
                      Ti_neg(s), zeta0_neg, p_neg,DT_neg(s), nu_neg(s,s), DF_neg(s,s), &
                      dTl_dnj(s,s),dzeta0_dnj(s),sigmai(s), v0, &
                      dDf_dnk(s,s,s),domegaij_dnj(s,s),domegaii_dni(s),Lij(s,s)

      integer i,j,k,p,kk,counter,l
      double precision kronecker(s,s),sum1(s,s),l_bar(s,s),lkj(s,s), &
                      Lkin(s,s),Lcol(s,s),Amat(s,s),bmat(s,s), &
                      Amat0(s,s),bmat0(s)
      double precision indx(s),d

!variables used to calculate dDF_dnk
      double precision perturbation,delta_ni(s),delta_phii(s),delta_rhoi(s)

      integer NP
      parameter (NP=15)     !max no. of linear equations to solve

      double precision pi
      parameter (pi=3.14159265458979323846d0)
      
      counter = 0
      
      do i=1,s
          if(ni(i) .eq. 0.d0) then
	     counter=counter+1
	  endif
         do j=1,s
            if (i.eq.j) then
               kronecker(i,j) = 1.d0
            else
               kronecker(i,j) = 0.d0
            endif
         enddo
      enddo


!Numerical evaulation of several partial derivatives - p 6.6 of CMH notes.
!Each of the partial derivatives is with respect to nj, where T and other
!ni are constant, so need to 
!   1) perturb nj a small positive amount (and correspondingly n,phii,phi)
!   2) evaulate dependent quantity of interest
!   3) perturb nj a small negative amount (and correspondingly n,phii,phi)
!   4) evaluate dependent quantity of interest
!   5) calculate partial derivative
!   6) return nj & n to original values
! This is done to evaluate a limit of thermal mobility as ni->0, but nj > 0
   if(counter > 0) then
     perturbation = 0.001d0    !perturb 0.1% of current value in each direction
      do j=1,s
         !conditional in the event one of the particle number densities is 0
         do l=1,s
         if(ni(l) .eq. 0.d0)then
           ni(l) = 1D-9
           phii(l) = pi*ni(l)*sigmai(l)**3/6.d0
           rhoi(l) = ni(l)*mi(l)
         endif
         enddo
         !calculate perturbation amounts
            delta_ni(j)  = perturbation*ni(j) 
            delta_phii(j)= perturbation*ni(j)*pi/6.d0*sigmai(j)**3 
            delta_rhoi(j)   = perturbation*ni(j)*mi(j) !;
         !perturb current nj values in + direction
            ni(j)   = ni(j) + delta_ni(j) 
            n       = n     + delta_ni(j) 
            phii(j) = phii(j)+ delta_phii(j) 
            phi     = phi    + delta_phii(j) 
            rhoi(j) = rhoi(j) + delta_rhoi(j) 
            rho = rho + delta_rhoi(j)!;
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
            
	    call thermal_diffusivity(s,alpha,ni,mi,rho,v0,mu,sigma,chi,zeta0_pos,&
                                     theta_pos,Ti_pos,p_pos,DT_pos,nu_pos)

            call mass_mobility(s,mi,ni,rho,zeta0_pos,theta_pos,nu_pos,DF_pos)

         !perturb current nj values in - direction
            ni(j)   = ni(j) - 2.d0*delta_ni(j)      
            n       = n     - 2.d0*delta_ni(j) 
            phii(j) = phii(j) - 2.d0*delta_phii(j) 
            phi     = phi     - 2.d0*delta_phii(j) 
            rhoi(j) = rhoi(j) - 2.d0*delta_rhoi(j) !;
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
	    
            call thermal_diffusivity(s,alpha,ni,mi,rho,v0,mu,sigma,chi,zeta0_neg,&
                                     theta_neg,Ti_neg,p_neg,DT_neg,nu_neg)

            call mass_mobility(s,mi,ni,rho,zeta0_neg,theta_neg,nu_neg,DF_neg)
         do k=1,s
            do l=1,s
            dDF_dnk(k,l,j) = (DF_pos(k,l) - DF_neg(k,l))/(2.d0*delta_ni(j))
            enddo
         enddo

         !reset nj et al. to original values
            ni(j)   = ni(j) + delta_ni(j)      
            n       = n     + delta_ni(j) 
            phii(j) = phii(j) + delta_phii(j) 
            phi     = phi     + delta_phii(j) 
            rhoi(j) = rhoi(j) + delta_rhoi(j) !;

            if(ni(j) .eq. 1D-9)then
 		ni(j) = 0.d0
		phii(j) = 0.d0
                rhoi(j) = 0.d0
            endif

            do i=1,s
               do k=1,s
                  call chi_ij_GHD(s,i,k,sigmai,phi,ni,chi_ij)
                  chi(i,k) = chi_ij
               enddo
            enddo
      enddo
   endif
!calculate summation used in l_bar (b vector) - p. 19 of CMH notes
      do i=1,s
         do j=1,s
            sum1(i,j) = 0.d0
         enddo
      enddo


      do i=1,s
         do j=1,s
            do k=1,s
               if(counter .eq. 0)then
                  sum1(i,j) = sum1(i,j) + (omega(i,k)-zeta0* &
                   kronecker(i,k))/(ni(k)*Ti(k))*DF(k,j)
               else
                 if(i .ne. k)then
		  sum1(i,j) = sum1(i,j) + &
                  (kronecker(i,k)*Ti(i)**2*DF(k,j)*(omega(i,k)-zeta0*kronecker(i,k))+ &
                   2*Ti(i)*ni(i)*DF(k,j)*dTl_dnj(i,k)*(omega(i,k)-kronecker(i,k))+ &
                   ni(i)*Ti(i)**2*dDF_dnk(k,j,k)*(omega(i,k)-zeta0*kronecker(i,k)) + &
                   ni(i)*Ti(i)**2*DF(k,j)*(domegaij_dnj(i,k)-dzeta0_dnj(k)*kronecker(i,k)))/(Ti(k)+ni(k)*dTl_dnj(k,k))
                 else
		  sum1(i,j) = sum1(i,j) + &
                  (kronecker(i,k)*Ti(i)**2*DF(k,j)*(omega(i,k)-zeta0*kronecker(i,k))+ &
                   2*Ti(i)*ni(i)*DF(k,j)*dTl_dnj(i,k)*(omega(i,k)-kronecker(i,k))+ &
                   ni(i)*Ti(i)**2*dDF_dnk(k,j,k)*(omega(i,k)-zeta0*kronecker(i,k)) + &
                   ni(i)*Ti(i)**2*DF(k,j)*(domegaii_dni(i)-dzeta0_dnj(k)*kronecker(i,k)))/(Ti(k)+ni(k)*dTl_dnj(k,k))
                 endif
               endif
            enddo
         enddo
      enddo

!calculate l_bar (p 19 of CMH notes)
      do i=1,s     
         do j=1,s

            if (counter .eq. 0)then
               l_bar(i,j) = -2.5d0*ni(i)*Ti(i)**2/mi(i)*sum1(i,j)
            else
               l_bar(i,j) = -2.5d0/mi(i)*sum1(i,j)
            endif

            Amat(i,j) = gammaij(i,j) - 0.5d0*zeta0*kronecker(i,j)     !A matrix for solution of lkj (p 19 CMH notes)
            bmat(i,j) = l_bar(i,j)                                    !B matrix for solution of lkj (p 19 CMH notes)
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
            lkj(i,kk) = bmat0(i)
         enddo
      enddo

!kinetic contribution to thermal mobility (p 19 CMH notes)
      do i = 1,s
         do j=1,s
            Lkin(i,j) = lkj(i,j) + 2.5d0*Ti(i)/mi(i)*DF(i,j) 
         enddo
      enddo

!collisional contribution to thermal mobility (p 20 CMH notes)
      do i=1,s
         do j=1,s
            Lcol(i,j) = 0.d0
         enddo
      enddo
      do i=1,s
         do j=1,s
            do p=1,s
               Lcol(i,j) = Lcol(i,j) + (1.d0+alpha(i,p))/8.d0*mi(p)* &
                 mu(i,p)*sigma(i,p)**3*chi(i,p)*(4.d0*pi/5.d0* &
                 (1.d0-alpha(i,p))*(mu(i,p)-mu(p,i))*ni(i)  &
                 *(2.d0/mi(p)*Lkin(p,j)+5.d0*Ti(i)/mi(i)/mi(p)* &
                 DF(p,j))+48.d0*pi/15.d0*ni(i)*(2.d0*mu(p,i)/mi(p)* &
                 Lkin(p,j)-5.d0*(2.d0*mu(i,p)-mu(p,i))*Ti(i)/mi(i)/ &
                 mi(p)*DF(p,j))) 
            enddo

            Lij(i,j) = Lkin(i,j) + Lcol(i,j)          !thermal mobility (p 19 CMH notes)
         enddo
      enddo

      return
      end subroutine thermal_mobility

