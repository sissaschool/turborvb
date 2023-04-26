!TL off
         subroutine makefun_pbc(iopt,iocc,indt,i0,indtmin,indtm,typec,indpar&
     &,indorb,indshell,nelskip,z,dd,r,rmu,distp &
     &,iflagnorm_fake,c,rmucos,rmusin,sinphase,cosphase)                            
         use constants 
         use Cell,      only:cellscale,rphase 
         implicit none 
         integer, intent(in) :: iopt,iocc(*),nelskip,indt       &
     &,iflagnorm_fake,typec 
         double precision, intent(out) ::  z(nelskip,i0:*)              
                                                                        
         double precision, intent(in) :: rmu(3,0:indtm),r(0:indtm)&
     &,rmucos(3,0:indtm),rmusin(3,0:indtm),sinphase(3,0:indtm)&
     &,cosphase(0:indtm)                                      
                                                                        
         integer i,k,indpar,indorbp,indorb,indshell,indshellp,ic,indtmin&
     &,indparp,indtm,i0                                   
                                                                        
         real*8 distp(0:indtm,20),fun,fun0                            &
     &,fun2,rp0,rp1,rp2,dd1,dd2,c,dsqrt,c0,c1,dd(*),hess3d(3,5)         &
     &,radhess4f(3),der4f(2),hess4f(3)         
                                                                        
         real(8) hess(3),fun3
                                                                        
         integer npower
                                                                        
!        indorb are the number of orbitals occupied before calling      
!        this subroutine                                                
!                                                                       
!        indpar is the number of variational parameters used            
!        before calling this subroutine                                 
!                                                                       
!        indshell is the index of the last  occupied orbital            
!        in the shell, characterized by occupation number iocc(indshell+
!                                                                       
!        z(i,indt+4)   contains the laplacian of the orbital i          
!        z(i,indt+mu) contains the gradient of the orbital i (mu=1,2,3) 
!       fun0 = f(r)   ! the radial part without the angular component
!       fun =1/r  d/dr f(r)                                             
!       fun3 =(1/r  d/dr ) ( 1/r d/dr ) f(r)                            
                                                                        
       
      
                                                                        
      select case(iopt) 
                     ! 1s single Z NO CUSP!                             
      case(1) 
                                                                        
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            dd1=dd(indpar+1) 
!           if(iflagnorm.gt.2) then 
!           c=dd1*dsqrt(dd1)/dsqrt(pi) 
            c=dd1*dsqrt(dd1)*0.56418958354775628695d0
!           endif 

                                                                        
            indorbp=indorb+1 
               do k=indtmin,indtm 
                  distp(k,1)=c*dexp(-dd1*r(k)) 
               enddo 
                                                                        
            do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*cosphase(i) 
            enddo 
                                                                        
            if(typec.ne.1) then 
                                                                        
               fun=-dd1*distp(0,1)/r(0) 
               fun3=(-dd1-1.d0/r(0))*fun/r(0) 
                                                                        
               do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &          
                & -distp(0,1)*sinphase(i,0)*rphase(i)
               enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2       &
      &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                                   
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*distp(0,1)*rphase(:)**2                  
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        


                                                                        
            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                                                                        
                    ! 2s gaussian periodic equal to 16                  
      case(161,16) 
! R(r)=exp(-z*r**2) single zeta                                         
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
!            if(iflagnorm.gt.2) then 
!           we use here the normalization of the                        
!           non periodic gaussian if dd1>0                              
!           unnormalized for dd1< 0 (hopefully small)                   
                  c=.71270547035499016035d0*dd1**0.75d0 
!                 if(c.gt.0) then 
!                 c=c**0.75d0 
!                 else 
!                 c=1.d0 
!                 endif 
!            endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
            do i=i0,indtm 
            z(indorbp,i)=distp(i,1)*cosphase(i) 
            enddo 
                                                                        
            if(typec.ne.1) then 
            fun=-2.d0*dd1*distp(0,1) 
            fun3=4.d0*dd1**2*distp(0,1) 
                                                                        
           do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -distp(0,1)*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2    &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)  
     
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
    &*sinphase(:,0)*rphase(:)-cosphase(0)*distp(0,1)*rphase(:)**2     
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
            indorb=indorbp 
         endif 
                                                                        
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
      


       case(20) ! single exp p orbital periodic                    
         dd1=dd(indpar+1) 
!        if(iflagnorm.gt.2) then
!         if(dd1.gt.0) then
          c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0
!         c=dd1**2.5d0*0.5641895835477562d0
!         else
!         c=1.d0
!         endif
!        endif
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                z(indorbp,i)=distp(i,1)*rmu(ic,i)*rmucos(ic,i)*cosphase(i)
              enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            rp1=dd1*(1.d0+dd1*r(0)) 
            fun0=distp(0,1) 
            fun =-dd1*distp(0,1)/r(0) 
            fun3=rp1*distp(0,1)/r(0)**3 
            indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
           do i=1,3
z(indorbp,indt+i)=rmu(i,0)*fun*rmucos(i,0) &
   &*rmu(ic,0)*rmucos(ic,0)*cosphase(0) & ! 2) 
   & -fun0*sinphase(i,0)*rphase(i)  &  
   &*rmu(ic,0)*rmucos(ic,0)  ! 1) 
            enddo 

z(indorbp,indt+ic)=z(indorbp,indt+ic) &
     &+fun0*(2.d0*rmucos(ic,0)**2-1.d0)*cosphase(0) ! 3)                    
        hess(:)=(rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)))*cosphase(0)                              
        hess(ic)=hess(ic)                                               &
     &+(fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-fun0*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic))*cosphase(0)                                                 
        hess(:)=hess(:)-rmu(:,0)*rmucos(:,0)  &    !from 1) 
     &*rphase(:)*fun*sinphase(:,0)                  &
     &*rmucos(ic,0)*rmu(ic,0)    

         hess(:)=hess(:)+rphase(:)*(     & ! from 2)      
     &-fun*rmu(:,0)*rmucos(:,0)*sinphase(:,0)    &
     &-fun0*cosphase(0)    & ! from 1)   
     &*rphase(:))*rmucos(ic,0)*rmu(ic,0)

         hess(ic)=hess(ic)-2.d0*fun0*(2.d0*rmucos(ic,0)**2-1.d0) &
     &*rphase(ic)*sinphase(ic,0)        ! from 3) and 1) 

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
               enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp            



          case(40) ! derivative of 20 with respect zeta 
! R(r)=c*(c1-r)*exp(-z*r)                                     
                                                                        
         dd1=dd(indpar+1) 
!        if(iflagnorm.gt.2) then
!         if(dd1.gt.0) then
          c=dd1**2.5d0*0.5641895835477562d0
!         c=dsqrt((2.d0*dd1)**5/8.d0/pi)/2.d0
!         else
!         c=1.d0
!         endif
!        endif

                                                                        
            if(dd1.gt.0.d0) then 
            c1=2.5d0/dd1 
            else 
            c1=0.d0 
            endif 
                                                                        
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
         indorbp=indorb 
                                                                        
                                                                        
!                                                                       
    do ic=1,3 
       if(iocc(indshell+ic).eq.1) then 
          indorbp=indorbp+1 
          do i=i0,indtm 
             z(indorbp,i)=rmu(ic,i)*distp(i,1)*(c1-r(i))*rmucos(ic,i)*cosphase(i)           
          enddo 
       endif 
    enddo 
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
                                                                        
            rp0=r(0) 
            rp1=dd1*rp0 
            rp2=rp1*rp1 
                                                                        
            fun0=distp(0,1)*(c1-rp0) 
            fun=distp(0,1)*(rp1-1.d0-c1*dd1)/r(0) 
            fun3=distp(0,1)*(1.d0+rp1-rp2+c1*dd1*(1.d0+rp1))/r(0)**3                           
                                                                        
           indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 

           do i=1,3
z(indorbp,indt+i)=rmu(i,0)*fun*rmucos(i,0) &
   &*rmu(ic,0)*rmucos(ic,0)*cosphase(0) & ! 2) 
   & -fun0*sinphase(i,0)*rphase(i)  &  
   &*rmu(ic,0)*rmucos(ic,0)  ! 1) 
            enddo 

z(indorbp,indt+ic)=z(indorbp,indt+ic) &
     &+fun0*(2.d0*rmucos(ic,0)**2-1.d0)*cosphase(0) ! 3)                    
        hess(:)=(rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)))*cosphase(0)                              
        hess(ic)=hess(ic)                                               &
     &+(fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-fun0*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic))*cosphase(0)                                                 
        hess(:)=hess(:)-rmu(:,0)*rmucos(:,0)  &    !from 1) 
     &*rphase(:)*fun*sinphase(:,0)                  &
     &*rmucos(ic,0)*rmu(ic,0)    

         hess(:)=hess(:)+rphase(:)*(     & ! from 2)      
     &-fun*rmu(:,0)*rmucos(:,0)*sinphase(:,0)    &
     &-fun0*cosphase(0)    & ! from 1)   
     &*rphase(:))*rmucos(ic,0)*rmu(ic,0)

         hess(ic)=hess(ic)-2.d0*fun0*(2.d0*rmucos(ic,0)**2-1.d0) &
     &*rphase(ic)*sinphase(ic,0)        ! from 3) and 1) 

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
               enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp            


                ! 2s gaussian periodic                                  
      case(17) 
! R(r)=r^2*exp(-z*r**2) single zeta                                     
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
!           if(iflagnorm.gt.2) then 
!           we use here the normalization of the                        
!           non periodic gaussian if dd1>0                              
!           unnormalized for dd1< 0 (hopefully small)                   
!                 if(dd1.gt.0) then 
!        c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0) 
         c=0.73607904464954686606d0*dd1**1.75d0
!                 else 
!                 c=1.d0 
!                 endif 
!           endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
            do i=i0,indtm 
            z(indorbp,i)=distp(i,1)*r(i)**2*cosphase(i) 
            enddo 
                                                                        
            if(typec.ne.1) then 
            rp1=dd1*r(0)**2 
                                                                        
            fun0=distp(0,1)*r(0)**2
            fun=2.d0*(1.d0-rp1)*distp(0,1) 
            fun3=4.d0*dd1*(rp1-2.d0)*distp(0,1) 
                                                                        
           do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2 

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
            indorb=indorbp 
         endif 
                                                                        
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                 !   derivative of 17                                   
       case(46) 
                                                                        
! R(r)=exp(-z*r**2)*(-c*r^4  +c1 r^2)                                   
                                                                        
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
!           if(iflagnorm.gt.2) then 
!                 if(dd1.gt.0) then 
!        c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0) 
         c=.73607904464954686606d0*dd1**1.75d0
!        c=4.d0*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)/dsqrt(15.d0) 
!                 else 
!                 c=1.d0 
!                 endif 
!           endif 
                                                                        
            if(dd1.gt.0.d0) then 
            c1=1.75d0/dd1 
            else 
            c1=0.d0 
            endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
                                                                        
  do i=i0,indtm 
    z(indorbp,i)=distp(i,1)*r(i)**2*(c1-r(i)**2)*cosphase(i) 
   enddo 
            if(typec.ne.1) then 
                                                                        
            rp1=r(0)**2 
            rp2=dd1*rp1 
                                                                        
           fun0=distp(0,1)*r(0)**2*(c1-r(0)**2)
           fun=2.d0*distp(0,1)*(c1*(1-rp2)+rp1*rp2-2.d0*rp1) 
           fun3=4.d0*distp(0,1)*(4*rp2-2.d0-rp2**2+c1*dd1*(rp2-2.d0)) 
                                                                        
           do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
! ********** DIAGONAL PART HESSIAN ******************                   
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(rmucos(:,0)**2                                          &
     &-PI*rmu(:,0)*rmusin(:,0)/cellscale(1:3)))*cosphase(0)                    
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2 

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                 ! r^3 x Gaussian                                       
       case(60) 
! R(r)=r^3*exp(-z*r**2) single zeta                                     
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
                                                                        
                                                                        
!          if(iflagnorm.gt.2) then 
!                 if(dd1.gt.0) then 
!     c=2.d0/pi**(3.d0/4.d0)*(2.d0*dd1)**(9.d0/4.d0)*dsqrt(2.d0/105.d0) 
      c=dd1**2.25d0*.55642345640820284397d0
!                 else 
!                 c=1.d0 
!                 endif 
!           endif 
                                                                        
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
            do i=i0,indtm 
            z(indorbp,i)=distp(i,1)*r(i)**3*cosphase(i) 
            enddo 
                                                                        
            if(typec.ne.1) then 
            rp1=dd1*r(0)**2 
                                                                        

            fun0=distp(0,1)*r(0)**3     
            fun=r(0)*(3.d0-2.d0*rp1)*distp(0,1) 
            fun3=(3.d0-12.d0*rp1+4.d0*rp1**2)*distp(0,1)/r(0) 
                                                                        
      do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &          
               & -fun0*sinphase(i,0)*rphase(i)
      enddo                                                              
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2 

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
            indorb=indorbp 
         endif 
                                                                        
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                 !   derivative of 60                                   
       case(61) 
                                                                        
! R(r)=exp(-z*r**2)*(-c *r^5  +c1 r^3)                                  
                                                                        
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
!          if(iflagnorm.gt.2) then 
!                 if(dd1.gt.0) then 
     c=dd1**2.25d0*.55642345640820284397d0
!
!     c=2.d0/pi**(3.d0/4.d0)*(2.d0*dd1)**(9.d0/4.d0)*dsqrt(2.d0/105.d0) 
!                 else 
!                 c=1.d0 
!                 endif 
!           endif 
                                                                        
!           if(dd1.gt.0.d0) then 
            c1=2.25d0/dd1 
!           else 
!           c1=0.d0 
!           endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
                                                                        
    do i=i0,indtm 
      z(indorbp,i)=distp(i,1)*r(i)**3*(c1-r(i)**2)*cosphase(i) 
    enddo 
            if(typec.ne.1) then 
                                                                        
                                                                        
                                                                        
                                                                        
            rp1=r(0)**2 
            rp2=dd1*rp1 
                                                                        
           fun0=distp(0,1)*r(0)**3*(c1-r(0)**2)
           fun=r(0)*distp(0,1)*(c1*(3.d0-2.d0*rp2)                  &
     &+2.d0*rp1*rp2-5.d0*rp1)                                           
           fun3=distp(0,1)*(-15.d0*rp1+20.d0*rp1*rp2-4.d0*rp2**2*rp1    &
     &+c1*(3.d0-12.d0*rp2+4.d0*rp2**2))/r(0)                        
                                                                        
              do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &          
                & -fun0*sinphase(i,0)*rphase(i)
              enddo                                                     
! ********** DIAGONAL PART HESSIAN ******************                   
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(rmucos(:,0)**2                                          &
     &-PI*rmu(:,0)*rmusin(:,0)/cellscale(1:3)))*cosphase(0)                    
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2 

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                                                                        
                 ! single gaussian p orbital periodic                   
       case(36) 
                                                                        
        dd1=dd(indpar+1) 
! normalization for non periodic one                                    
!        if(iflagnorm.gt.2) then
!         if(dd1.gt.0) then
!         c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0
!         else
!         c=dsqrt(2.d0)*pi*(-0.75d0)
!         endif
          c=dd1**1.25d0*1.42541094070998d0
!        endif                                                                  
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
z(indorbp,i)=distp(i,1)*rmu(ic,i)*rmucos(ic,i)*cosphase(i)
              enddo 
            endif 
         enddo 
                                                                        
    if(typec.ne.1) then 
            fun =-2.d0*dd1*distp(0,1) 
            fun3=4.d0*dd1**2*distp(0,1) 
                                                                        
               indorbp=indorb 
                                                                        
    do ic=1,3 
       if(iocc(indshell+ic).eq.1) then 
          indorbp=indorbp+1 
            do i=1,3
z(indorbp,indt+i)=rmu(i,0)*fun*rmucos(i,0) &
   &*rmu(ic,0)*rmucos(ic,0)*cosphase(0) & ! 2) 
   & -distp(0,1)*sinphase(i,0)*rphase(i)  &  
   &*rmu(ic,0)*rmucos(ic,0)  ! 1) 
            enddo 

z(indorbp,indt+ic)=z(indorbp,indt+ic) &
     &+distp(0,1)*(2.d0*rmucos(ic,0)**2-1.d0)*cosphase(0) ! 3)                    
        hess(:)=(rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)))*cosphase(0)                              
        hess(ic)=hess(ic)                                               &
     &+(fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-distp(0,1)*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic))*cosphase(0)                                                 
        hess(:)=hess(:)-rmu(:,0)*rmucos(:,0)  &    !from 1) 
     &*rphase(:)*fun*sinphase(:,0)                  &
     &*rmucos(ic,0)*rmu(ic,0)    

         hess(:)=hess(:)+rphase(:)*(     & ! from 2)      
     &-fun*rmu(:,0)*rmucos(:,0)*sinphase(:,0)    &
     &-distp(0,1)*cosphase(0)    & ! from 1)   
     &*rphase(:))*rmucos(ic,0)*rmu(ic,0)

         hess(ic)=hess(ic)-2.d0*distp(0,1)*(2.d0*rmucos(ic,0)**2-1.d0) &
     &*rphase(ic)*sinphase(ic,0)        ! from 3) and 1)

!      z_xyz(indorbp,:)=hess(:) 
       z(indorbp,indt+4)=sum(hess(:)) 
       endif 
     enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                                                                        
                 ! single gaussian p orbital periodic                   
       case(62) 
                                                                        
         dd1=dd(indpar+1) 
!        if(iflagnorm.gt.2) then 
! normalization for non periodic one                                    
!        if(dd1.gt.0) then 
       c=2.d0/pi**0.75d0*(2.d0*dd1)**1.75d0/dsqrt(5.d0) 
!      c=dd1**1.75d0*1.2749263037197753d0
!        else 
!        c=1.d0 
!        endif 
!        endif 
                                                                        
                                                                        
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
      z(indorbp,i)=distp(i,1)*rmu(ic,i)*rmucos(ic,i)*r(i)*cosphase(i) 
              enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            rp1=2.d0*dd1*r(0)**2 
            fun0=distp(0,1)*r(0) 
            fun =(1.d0-rp1)*distp(0,1)/r(0) 
            fun3=(rp1**2-1.d0-2.d0*rp1)*distp(0,1)/r(0)**3 
                                                                        
               indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                                                                        
           do i=1,3
z(indorbp,indt+i)=rmu(i,0)*fun*rmucos(i,0) &
   &*rmu(ic,0)*rmucos(ic,0)*cosphase(0) & ! 2) 
   & -fun0*sinphase(i,0)*rphase(i)  &  
   &*rmu(ic,0)*rmucos(ic,0)  ! 1) 
            enddo 

z(indorbp,indt+ic)=z(indorbp,indt+ic) &
     &+fun0*(2.d0*rmucos(ic,0)**2-1.d0)*cosphase(0) ! 3)                    
        hess(:)=(rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)))*cosphase(0)                              
        hess(ic)=hess(ic)                                               &
     &+(fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-fun0*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic))*cosphase(0)                                                 
        hess(:)=hess(:)-rmu(:,0)*rmucos(:,0)  &    !from 1) 
     &*rphase(:)*fun*sinphase(:,0)                  &
     &*rmucos(ic,0)*rmu(ic,0)    

         hess(:)=hess(:)+rphase(:)*(     & ! from 2)      
     &-fun*rmu(:,0)*rmucos(:,0)*sinphase(:,0)    &
     &-fun0*cosphase(0)    & ! from 1)   
     &*rphase(:))*rmucos(ic,0)*rmu(ic,0)

         hess(ic)=hess(ic)-2.d0*fun0*(2.d0*rmucos(ic,0)**2-1.d0) &
     &*rphase(ic)*sinphase(ic,0)        ! from 3) and 1)

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
            enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                    ! derivative of 36 with respect zeta                
           case(63) 
! R(r)=c*x*r*exp(-z*r**2)*(c1-r**2)                                     
                                                                        
          
         dd1=dd(indpar+1) 
                                                                        
!         if(iflagnorm.gt.2) then 
! normalization for non periodic one                                    
!         if(dd1.gt.0) then 
    c=dd1**1.75d0*1.2749263037197753d0
!       c=2.d0/pi**0.75d0*(2.d0*dd1)**1.75d0/dsqrt(5.d0) 
!         else 
!         c=1.d0 
!         endif 
!         endif 
                                                                        
                                                                        
                                                                        
            if(dd1.gt.0.d0) then 
            c1=1.75d0/dd1 
            else 
            c1=0.d0 
            endif 
                                                                        
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
                                                                        
                                                                        
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                   z(indorbp,i)=rmu(ic,i)*distp(i,1)*               &
     &     r(i)*(c1-r(i)**2)*rmucos(ic,i)*cosphase(i)            
               enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
                                                                        
            rp0=r(0)**2 
            rp1=2*dd1*rp0 
            rp2=rp1*rp0 
                                                                        
                                                                        
                                                                        
            fun=distp(0,1)*(c1*(1.d0-rp1)+rp2-3.d0*rp0)/r(0) 
            fun0=r(0)*distp(0,1)*(c1-rp0) 
        fun3=-distp(0,1)*(c1*(1.d0+2.d0*rp1-2.d0*rp2*dd1)               &
     &+3.d0*rp0-6.d0*rp2+rp1*rp2)/r(0)**3                           
                                                                        
           indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                    
          do i=1,3
z(indorbp,indt+i)=rmu(i,0)*fun*rmucos(i,0) &
   &*rmu(ic,0)*rmucos(ic,0)*cosphase(0) & ! 2) 
   & -fun0*sinphase(i,0)*rphase(i)  &  
   &*rmu(ic,0)*rmucos(ic,0)  ! 1) 
            enddo 

        z(indorbp,indt+ic)=z(indorbp,indt+ic) &
     &+fun0*(2.d0*rmucos(ic,0)**2-1.d0)*cosphase(0) ! 3)                    
        hess(:)=(rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)))*cosphase(0)                              
        hess(ic)=hess(ic)                                               &
     &+(fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-fun0*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic))*cosphase(0)                                                 
        hess(:)=hess(:)-rmu(:,0)*rmucos(:,0)  &    !from 1) 
     &*rphase(:)*fun*sinphase(:,0)                  &
     &*rmucos(ic,0)*rmu(ic,0)    

         hess(:)=hess(:)+rphase(:)*(     & ! from 2)      
     &-fun*rmu(:,0)*rmucos(:,0)*sinphase(:,0)    &
     &-fun0*cosphase(0)    & ! from 1)   
     &*rphase(:))*rmucos(ic,0)*rmu(ic,0)

         hess(ic)=hess(ic)-2.d0*fun0*(2.d0*rmucos(ic,0)**2-1.d0) &
     &*rphase(ic)*sinphase(ic,0)        ! from 3) and 1)

!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
               enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                    ! derivative of 36 with respect zeta                
           case(44) 
! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)                                      
                                                                        
          
         dd1=dd(indpar+1) 
!       if(iflagnorm.gt.2) then 
!         if(dd1.gt.0) then
!         c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0 
          c=dd1**1.25d0*1.42541094070998d0
!         else
!         c=dsqrt(2.d0)*pi*(-0.75d0) 
!         endif 
!        endif 

         if(dd1.gt.0) then
         c1=1.25d0/dd1
         else
         c1=0.d0
         endif
                                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                   z(indorbp,i)=rmu(ic,i)*distp(i,1)*               &
     &             (c1-r(i)**2)*rmucos(ic,i)*cosphase(i)            
               enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
        
            fun=2.d0*distp(0,1)*(dd1*r(0)**2-1.d0-c1*dd1) 
            fun0=distp(0,1)*(c1-r(0)**2) 
            fun3=4.d0*dd1*distp(0,1)*(2.d0-r(0)**2*dd1+c1*dd1)

                                                                        
           indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 

           do i=1,3
z(indorbp,indt+i)=rmu(i,0)*fun*rmucos(i,0) &
   &*rmu(ic,0)*rmucos(ic,0)*cosphase(0) & ! 2) 
   & -fun0*sinphase(i,0)*rphase(i)  &  
   &*rmu(ic,0)*rmucos(ic,0)  ! 1) 
            enddo 

z(indorbp,indt+ic)=z(indorbp,indt+ic) &
     &+fun0*(2.d0*rmucos(ic,0)**2-1.d0)*cosphase(0) ! 3)                    
        hess(:)=(rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)))*cosphase(0)                              
        hess(ic)=hess(ic)                                               &
     &+(fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-fun0*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic))*cosphase(0)                                                 
        hess(:)=hess(:)-rmu(:,0)*rmucos(:,0)  &    !from 1) 
     &*rphase(:)*fun*sinphase(:,0)                  &
     &*rmucos(ic,0)*rmu(ic,0)    

         hess(:)=hess(:)+rphase(:)*(     & ! from 2)      
     &-fun*rmu(:,0)*rmucos(:,0)*sinphase(:,0)    &
     &-fun0*cosphase(0)    & ! from 1)   
     &*rphase(:))*rmucos(ic,0)*rmu(ic,0)

         hess(ic)=hess(ic)-2.d0*fun0*(2.d0*rmucos(ic,0)**2-1.d0) &
     &*rphase(ic)*sinphase(ic,0)        ! from 3) and 1)
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
               enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                                                                        
                   ! derivative of 16 with respect to z                 
         case(19) 
! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)                                      
                                                                        
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
!           if(iflagnorm.gt.2) then 
!             c=2*dd1*INV_PI 
!             if(c.gt.0) then 
!             c=c**0.75d0 
!             else 
!             c=1.d0 
!             endif 
             c=0.71270547035499016d0*dd1**0.75d0
!           endif 
                                                                        
!           if(dd1.gt.0.d0) then 
            c1=0.75d0/dd1 
!           else 
!           c1=0.d0 
!           endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
                                                                        
            do i=i0,indtm 
             z(indorbp,i)=distp(i,1)*(c1-r(i)**2)*cosphase(i) 
            enddo 
            if(typec.ne.1) then 
                                                                        
           fun0=distp(0,1)*(c1-r(0)**2) 
           fun=2.d0*distp(0,1)*(dd1*r(0)**2-c1*dd1-1.d0) 
           fun3=-4.d0*dd1*distp(0,1)*(dd1*r(0)**2-c1*dd1-2.d0) 
                                                                        
           do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                       ! 2s single  Z WITH CUSP zero periodic           
          case(34) 
                                    ! normalized                        
! exp(-dd1*r) + dd1*r*exp(-dd1*r)                                       
                                                                        
          
         indshellp=indshell+1 
                                                                        
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
!           if(iflagnorm.gt.2) then 
!            if(dd1.gt.0) then 
!           c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+      &
!    &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)                           
!            else 
!            c1=1.d0 
!            endif 
             c=dd1*dsqrt(dd1)*.2132436186229231d0
!           endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)) 
            enddo 
                                                                        
            do i=i0,indtm 
            z(indorbp,i)=distp(i,1)*(1.d0+r(i)*dd1)*cosphase(i) 
            enddo 
                                                                        
            if(typec.ne.1) then 
                                                                        
             fun0=distp(0,1)*(1.d0+r(0)*dd1) 
             fun=-dd1**2*distp(0,1) 
             fun3=dd1**3*distp(0,1)/r(0) 


           do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                                                                        
                    ! 3s -derivative of 34 with respect to dd1          
      case(38) 
! R(r)=r**2*exp(-z1*r)                                                  
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
!           if(iflagnorm.gt.2) then 
!           if(dd1.gt.0.d0) then 
!           c=1.d0/dsqrt(1.d0/4.d0/dd1**3+12.d0*dd1/(2.d0*dd1)**4+      &
!    &3.d0*dd1**2/4.d0/dd1**5)/dsqrt(4.d0*pi)                           
!           else 
            c=dd1*dsqrt(dd1)*.2132436186229231d0
!           c=1.d0 
!           endif 
!           endif 
                                                                        
            c0=-c*dd1 
                                                                        
                                                                        
!           if(dd1.gt.0) then 
            c1=1.5d0*c/dd1 
!           else 
!           c1=0.d0 
!           endif 
                                                                        
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=dexp(-dd1*r(k)) 
            enddo 
                                                                        
            do i=i0,indtm 
            z(indorbp,i)=(c0*r(i)**2+c1*(1.d0+dd1*r(i)))        &
     &*distp(i,1)*cosphase(i)                                                       
            enddo 
                                                                        
            c1=c1*dd1**2 
                                                                        
            if(typec.ne.1) then 
                                                                        
                                                                        
     fun0=(c0*r(0)**2+c1/dd1**2*(1.d0+dd1*r(0)))*distp(0,1)

               fun=(c0*(2.d0-dd1*r(0))-c1)*distp(0,1) 
               fun3=-dd1*(fun+c0*distp(0,1))/r(0) 
                                                                        
           do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2      &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        
            endif 
            indorb=indorbp 
         endif 
         indpar=indpar+1 
         indshell=indshellp 
                       ! 2s single  Z WITH CUSP zero periodic           
          case(12) 
                                    ! normalized                        
! r**3*exp(-dd1*r)                                                      
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
                                                                        
!           if(iflagnorm.gt.2) then 
!           if(dd1.gt.0) then 
!           c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0 
            c=dd1**4.5d0*.03178848180059307346d0
!           else 
!           c=1.d0 
!           endif 
!           endif 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=c*dexp(-dd1*r(k)) 
            enddo 
                                                                        
            do i=i0,indtm 
            z(indorbp,i)=distp(i,1)*r(i)**3*cosphase(i) 
            enddo 
                                                                        
                                                                        
                                                                        
            if(typec.ne.1) then 
                                                                        
                                                                        
             fun0=distp(0,1)*r(0)**3 
             rp1=dd1*r(0) 
             fun=-distp(0,1)*r(0)*(rp1-3.d0) 
             fun3=distp(0,1)*(3.d0-5.d0*rp1+rp1**2)/r(0) 
                                                                        
                                                                        
            do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
             endif 
                                                                        
                                                                        
            indorb=indorbp 
                                                                        
         endif 
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                                                                        
                       ! derivative of 12 with respect to dd1           
          case(50) 
! c0 *r**4*exp(-dd1*r)+c1*r^3*exp(-dd1*r)                               
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
                                                                        
!           if(iflagnorm.gt.2) then 
!           if(dd1.gt.0.d0) then 
!           c=dsqrt((2*dd1)**9/40320.d0/pi)/2.d0 
            c=dd1**4.5d0*.03178848180059307346d0
!           else 
!           c=1.d0 
!           endif 
!           endif 
                                                                        
            c0=-c 
!           if(dd1.gt.0) then 
            c1=4.5d0*c/dd1 
!           else 
!           c1=0.d0 
!           endif 
            do k=indtmin,indtm 
            distp(k,1)=dexp(-dd1*r(k)) 
            enddo 
                                                                        
            do i=i0,indtm 
     z(indorbp,i)=(c0*r(i)**4+c1*r(i)**3)*distp(i,1)*cosphase(i) 
            enddo 
                                                                        
                                                                        
                                                                        
            if(typec.ne.1) then 
                                                                        
                                                                        
             rp1=dd1*r(0) 
             rp2=rp1*rp1 
                                                                        
           fun0=(c0*r(0)**4+c1*r(0)**3)*distp(0,1)
           fun=-distp(0,1)*r(0)                                   &
     &*(c0*r(0)*(rp1-4.d0)+c1*(rp1-3.d0))                           
                                                                        
             fun3=distp(0,1)*(c0*(8.d0-7.d0*rp1+rp2)                    &
     & +c1*(3.d0-5.d0*rp1+rp2)/r(0))                                

            do i=1,3 
z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)*cosphase(0) &
                  & -fun0*sinphase(i,0)*rphase(i)
            enddo
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=(fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))*cosphase(0)                               
hess(:)=hess(:)-2.d0*fun*rmu(:,0)*rmucos(:,0) &
              &*sinphase(:,0)*rphase(:)             &
              &-cosphase(0)*fun0*rphase(:)**2                                                        
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        
                                                                        
                                                                        
            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
         indpar=indpar+1 
         indshell=indshellp 
                                                                        
                                                                        
                                                                        
       case(37) 
! d orbitals                                                            
! R(r)= exp(-alpha r^2)                                                 
! each gaussian term is normalized                                      
                                                                        
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd1=dd(indparp) 
                                                                        
!        if(iflagnorm.gt.2) then 
! overall normalization                                                 
!        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0) 
         c=dd1**1.75d0*1.64592278064948967213d0
!        endif 
                                                                        
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=(3.d0*(rmu(3,i)*rmucos(3,i))**2                &
     &-r(i)**2)*cost1d                                              

                                                                        
      distp(i,3)=((rmu(1,i)*rmucos(1,i))**2                     &
     &-(rmu(2,i)*rmucos(2,i))**2)*cost2d                        
                                                 ! lz=+/-2              
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d                           
                                               ! lz=+/-2                
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d                           

                                               ! lz=+/-1                
          enddo 
                                                                        
                                                                        
         do ic=1,5 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k) 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            dd1=dd(indparp) 
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=4.d0*dd1**2*distp(0,1) 
                                                                        
       hess3d=0.d0 
                                                                        
       hess3d(:,1)=-2.d0*rmu(:,0)*rmucos(:,0)*cost1d 
                                                                        
       hess3d(3,1)=hess3d(3,1)                                          &
     &+6.d0*rmu(3,0)*rmucos(3,0)                           &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)*cost1d                   
                                                                        
       hess3d(1,2)=2.d0*rmu(1,0)*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)*cost2d                   
                                                                        
       hess3d(2,2)=-2.d0*rmu(2,0)*rmucos(2,0)                   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)*cost2d                   
                                                                        
                                                        ! lz=+/-2       
       hess3d(2,3)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-2       
       hess3d(1,3)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(2,4)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,4)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(1,5)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,5)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
               indorbp=indorb 
               do ic=1,5 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
           do i=1,3 
      z(indorbp,indt+i)=distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)   & ! 1)
     &*cosphase(0)-distp(0,1+ic)*fun0*sinphase(i,0)*rphase(i)  ! 2)
                                                                        
          hess(i)=((fun2*(rmu(i,0)*rmucos(i,0))**2               &
     &+fun*(rmucos(i,0)**2-rmusin(i,0)**2))*distp(0,1+ic)       &
     &+hess3d(i,ic)*fun*rmu(i,0)*rmucos(i,0))*cosphase(0)             

                  !enddo for i                                          
   hess(i)=hess(i)-distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)  & ! from 1
   &*sinphase(i,0)*rphase(i)  &
   &-distp(0,1+ic)*fun0*cosphase(0)*rphase(i)**2   & ! from 2
   &-distp(0,1+ic)*fun*sinphase(i,0)*rphase(i)*rmu(i,0)*rmucos(i,0)   & ! from 2
   &-fun0*sinphase(i,0)*rphase(i)*hess3d(i,ic)    ! from 2
           enddo 
                       if(ic.eq.1) then 
                    do i=1,3 

       z(indorbp,indt+i)=z(indorbp,indt+i)                              &
     &-2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d*cosphase(0) ! 3 
                                                                        
   hess(i)=hess(i)+2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d &
     &*sinphase(i,0)*rphase(i) ! from  3 

            hess(i)=hess(i)+(-2.d0*fun0*cost1d                                 &
     &*(rmucos(i,0)**2-rmusin(i,0)**2)                          &
     &-2.d0*fun*rmu(i,0)**2*rmucos(i,0)**2*cost1d)*cosphase(0)   

                    enddo 

       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0) ! 4                  

       hess(3)=hess(3)                              &
     &-(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) ! from 4                  
                                                                        
       hess(3)=hess(3)+(6.d0*cost1d*fun0                            &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)**2                       &
     &-6.d0*rmu(3,0)*rmucos(3,0)*cost1d                    &
     &*fun0*4.d0*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)        &
     &+6.d0*rmu(3,0)**2*rmucos(3,0)**2*cost1d              &
     &*fun*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      
                                                                        
                       elseif(ic.eq.2) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

       hess(1)=hess(1)-(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

      hess(2)=hess(2)-(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       hess(1)=hess(1)+(2.d0*rmu(1,0)**2                             &
     &*fun*cost2d*rmucos(1,0)**2                                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &+2.d0*fun0*cost2d*(rmucos(1,0)**2-rmusin(1,0)**2)         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &-2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*4.d0*PI/cellscale(1)*rmusin(1,0)*rmucos(1,0))*cosphase(0)             
                                                                        
       hess(2)=hess(2)+(-2.d0*rmu(2,0)**2                             &
     &*fun*cost2d*rmucos(2,0)**2                                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &-2.d0*fun0*cost2d*(rmucos(2,0)**2-rmusin(2,0)**2)         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &+2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*4.d0*PI/cellscale(2)*rmusin(2,0)*rmucos(2,0))*cosphase(0)             
                                                                        
                       elseif(ic.eq.3) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

        hess(1)=hess(1)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)           &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
        hess(2)=hess(2)+(-rmu(1,0)*cost3d*rmucos(1,0)             &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(1)=hess(1)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
                       elseif(ic.eq.4) then 
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)        &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                          

        hess(3)=hess(3)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)       &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3)                          
        hess(2)=hess(2)+(-rmu(3,0)*cost3d*rmucos(3,0)           &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                                                        
                       elseif(ic.eq.5) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)              

        hess(1)=hess(1)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      

        hess(3)=hess(3)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)      &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) 
                                                                        
       hess(1)=hess(1)+(-rmu(3,0)*cost3d*rmucos(3,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(1,0)*cost3d*rmucos(1,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                !endif for ic                           
                       endif 
                                                                        
!     z_xyz(indorbp,:)=hess(:) 
      z(indorbp,indt+4)=sum(hess(:)) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+5 
         indorb=indorbp 
                                                                        
                                                                        
       case(45) 
! derivative of the orbital 37 d orbitals                               
! d R(r) /d z = d exp(-alpha r^2) / dz = c*exp(-z r^2)*(7/4/z-r^2)      
! each gaussian term is normalized                                      
                                                                        
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0) 
          c=dd1**1.75d0*1.64592278064948967213d0
!         endif 
                                                                        
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=(3.d0*(rmu(3,i)*rmucos(3,i))**2                &
     &-r(i)**2)*cost1d                                              
                           ! lz=0                                       
                                                                        
      distp(i,3)=((rmu(1,i)*rmucos(1,i))**2                     &
     &-(rmu(2,i)*rmucos(2,i))**2)*cost2d                        
                                                 ! lz=+/-2              
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d                           
                                               ! lz=+/-2                
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
          enddo 
                                                                        
                                                                        
         do ic=1,5 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*(7.d0/4.d0/dd1-r(k)**2)*      &
     &        distp(k,1+ic)*cosphase(k)                        
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            dd1=dd(indparp) 
            fun0=distp(0,1)*(7.d0/4.d0/dd1-r(0)**2) 
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-11.d0/2.d0) 
            fun2=distp(0,1)*4.d0*dd1                                    &
     &-(2.d0*dd1*r(0)**2-11.d0/2.d0)*2.d0*dd1*distp(0,1)            
                                                                        
       hess3d=0.d0 
                                                                        
       hess3d(:,1)=-2.d0*rmu(:,0)*rmucos(:,0)*cost1d 
                                                                        
       hess3d(3,1)=hess3d(3,1)                                          &
     &+6.d0*rmu(3,0)*rmucos(3,0)                           &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)*cost1d                   
                                                                        
       hess3d(1,2)=2.d0*rmu(1,0)*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)*cost2d                   
                                                                        
       hess3d(2,2)=-2.d0*rmu(2,0)*rmucos(2,0)                   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)*cost2d                   
                                                                        
                                                        ! lz=+/-2       
       hess3d(2,3)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-2       
       hess3d(1,3)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(2,4)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,4)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(1,5)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,5)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
               indorbp=indorb 
               do ic=1,5 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                  

                   do i=1,3 
      z(indorbp,indt+i)=distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)   & ! 1)
     &*cosphase(0)-distp(0,1+ic)*fun0*sinphase(i,0)*rphase(i)  ! 2)
                                                                        
          hess(i)=((fun2*(rmu(i,0)*rmucos(i,0))**2               &
     &+fun*(rmucos(i,0)**2-rmusin(i,0)**2))*distp(0,1+ic)       &
     &+hess3d(i,ic)*fun*rmu(i,0)*rmucos(i,0))*cosphase(0)             

                  !enddo for i                                          
   hess(i)=hess(i)-distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)  & ! from 1
   &*sinphase(i,0)*rphase(i)  &
   &-distp(0,1+ic)*fun0*cosphase(0)*rphase(i)**2   & ! from 2
   &-distp(0,1+ic)*fun*sinphase(i,0)*rphase(i)*rmu(i,0)*rmucos(i,0)   & ! from 2
   &-fun0*sinphase(i,0)*rphase(i)*hess3d(i,ic)    ! from 2
           enddo 
                       if(ic.eq.1) then 
                    do i=1,3 

       z(indorbp,indt+i)=z(indorbp,indt+i)                              &
     &-2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d*cosphase(0) ! 3 
                                                                        
   hess(i)=hess(i)+2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d &
     &*sinphase(i,0)*rphase(i) ! from  3 

            hess(i)=hess(i)+(-2.d0*fun0*cost1d                                 &
     &*(rmucos(i,0)**2-rmusin(i,0)**2)                          &
     &-2.d0*fun*rmu(i,0)**2*rmucos(i,0)**2*cost1d)*cosphase(0)   

                    enddo 

       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0) ! 4                  

       hess(3)=hess(3)                              &
     &-(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) ! from 4                  
                                                                        
       hess(3)=hess(3)+(6.d0*cost1d*fun0                            &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)**2                       &
     &-6.d0*rmu(3,0)*rmucos(3,0)*cost1d                    &
     &*fun0*4.d0*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)        &
     &+6.d0*rmu(3,0)**2*rmucos(3,0)**2*cost1d              &
     &*fun*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      
                                                                        
                       elseif(ic.eq.2) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

       hess(1)=hess(1)-(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

      hess(2)=hess(2)-(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       hess(1)=hess(1)+(2.d0*rmu(1,0)**2                             &
     &*fun*cost2d*rmucos(1,0)**2                                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &+2.d0*fun0*cost2d*(rmucos(1,0)**2-rmusin(1,0)**2)         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &-2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*4.d0*PI/cellscale(1)*rmusin(1,0)*rmucos(1,0))*cosphase(0)             
                                                                        
       hess(2)=hess(2)+(-2.d0*rmu(2,0)**2                             &
     &*fun*cost2d*rmucos(2,0)**2                                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &-2.d0*fun0*cost2d*(rmucos(2,0)**2-rmusin(2,0)**2)         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &+2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*4.d0*PI/cellscale(2)*rmusin(2,0)*rmucos(2,0))*cosphase(0)             
                                                                        
                       elseif(ic.eq.3) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

        hess(1)=hess(1)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)           &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
        hess(2)=hess(2)+(-rmu(1,0)*cost3d*rmucos(1,0)             &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(1)=hess(1)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
                       elseif(ic.eq.4) then 
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)        &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                          

        hess(3)=hess(3)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)       &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3)                          
        hess(2)=hess(2)+(-rmu(3,0)*cost3d*rmucos(3,0)           &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                                                        
                       elseif(ic.eq.5) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)              

        hess(1)=hess(1)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      

        hess(3)=hess(3)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)      &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) 
                                                                        
       hess(1)=hess(1)+(-rmu(3,0)*cost3d*rmucos(3,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(1,0)*cost3d*rmucos(1,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                !endif for ic    

                       endif 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
      z(indorbp,indt+4)=sum(hess(:)) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+5 
         indorb=indorbp 
                                                                        
      case(68) 
! d single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 

         dd(indparp)=abs(dd(indparp))                                                                         
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization 
!         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0) 
         c=dd1**1.75d0*1.64592278064948967213d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=(3.d0*rmu(3,i)**2-rp0)*cost1d  ! lz=0
                                                                        
      distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d  ! lz=+/-2      
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d   ! lz=+/-2
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d  ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d  ! lz=+/-1    
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,5
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,5
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then

                   der4f(2)=-2.d0*cost1d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(1,0)**2 + rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost1d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                   der4f(2)=4.d0*cost1d*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*4.d0*(rmucos(3,0)**2-rmusin(3,0)**2)

                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then

                   der4f(2)=2.d0*cost2d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(rmucos(1,0)**2 - rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost2d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*rmusin(1,0)

                        elseif(i.eq.2) then


                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*rmusin(2,0)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        elseif(i.eq.2) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*rmusin(2,0)

                        else


                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(2,0)*rmucos(2,0)*rmucos(3,0)*rmusin(3,0)

                        endif                        


                     else


                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)*rmusin(1,0)

                        elseif(i.eq.2) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        else

                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(1,0)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)

                        endif                 

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+5
         indorb=indorbp 


      case(69) 
! d single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 
! orbital with the minimal power of rmucos
! orbital derivative of 68

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 

         dd(indparp)=abs(dd(indparp))                                                                         
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization 
!         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0) 
          c=dd1**1.75d0*1.64592278064948967213d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=(3.d0*rmu(3,i)**2-rp0)*cost1d  ! lz=0
                                                                        
      distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d  ! lz=+/-2      
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d   ! lz=+/-2
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d  ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d  ! lz=+/-1    
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,5
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*(7.d0/4.d0/dd1-r(k)**2)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1)*(7.d0/4.d0/dd1-r(0)**2)  
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-11.d0/2.d0)  
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +17.d0*dd1*r(0)**2-11.d0/2.d0)  

            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,5
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then

                   der4f(2)=-2.d0*cost1d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(1,0)**2 + rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost1d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                   der4f(2)=4.d0*cost1d*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*4.d0*(rmucos(3,0)**2-rmusin(3,0)**2)

                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then

                   der4f(2)=2.d0*cost2d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(rmucos(1,0)**2 - rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost2d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*rmusin(1,0)

                        elseif(i.eq.2) then


                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*rmusin(2,0)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        elseif(i.eq.2) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*rmusin(2,0)

                        else


                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(2,0)*rmucos(2,0)*rmucos(3,0)*rmusin(3,0)

                        endif                        


                     else


                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)*rmusin(1,0)

                        elseif(i.eq.2) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        else

                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(1,0)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)

                        endif                 

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(7.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(7.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+5
         indorb=indorbp 
                                                                        
                                                                        

       case(30) 
! d orbitals without cusp and one parmater                                                             
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd1=dd(indparp) 
                                                                        
!        if(iflagnorm.gt.2) then 
! overall normalization                                                 
!        c=1.d0/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0)
         c=dd1**3.5d0*0.26596152026762178d0
!        endif 
                                                                        
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=(3.d0*(rmu(3,i)*rmucos(3,i))**2                &
     &-r(i)**2)*cost1d                                              
                           ! lz=0                                       
                                                                        
      distp(i,3)=((rmu(1,i)*rmucos(1,i))**2                     &
     &-(rmu(2,i)*rmucos(2,i))**2)*cost2d                        
                                                 ! lz=+/-2              
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d                           
                                               ! lz=+/-2                
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
          enddo 
                                                                        
                                                                        
         do ic=1,5 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k) 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            dd1=dd(indparp) 
            fun0=distp(0,1) 
            fun=-dd1*distp(0,1)/r(0) 
            fun2=(dd1**2*distp(0,1)+dd1*distp(0,1)/r(0))/r(0)**2
                                                                        
       hess3d=0.d0 
                                                                        
       hess3d(:,1)=-2.d0*rmu(:,0)*rmucos(:,0)*cost1d 
                                                                        
       hess3d(3,1)=hess3d(3,1)                                          &
     &+6.d0*rmu(3,0)*rmucos(3,0)                           &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)*cost1d                   
                                                                        
       hess3d(1,2)=2.d0*rmu(1,0)*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)*cost2d                   
                                                                        
       hess3d(2,2)=-2.d0*rmu(2,0)*rmucos(2,0)                   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)*cost2d                   
                                                                        
                                                        ! lz=+/-2       
       hess3d(2,3)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-2       
       hess3d(1,3)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(2,4)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,4)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(1,5)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,5)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
               indorbp=indorb 
               do ic=1,5 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
           do i=1,3 
      z(indorbp,indt+i)=distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)   & ! 1)
     &*cosphase(0)-distp(0,1+ic)*fun0*sinphase(i,0)*rphase(i)  ! 2)
                                                                        
          hess(i)=((fun2*(rmu(i,0)*rmucos(i,0))**2               &
     &+fun*(rmucos(i,0)**2-rmusin(i,0)**2))*distp(0,1+ic)       &
     &+hess3d(i,ic)*fun*rmu(i,0)*rmucos(i,0))*cosphase(0)             

                  !enddo for i                                          
   hess(i)=hess(i)-distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)  & ! from 1
   &*sinphase(i,0)*rphase(i)  &
   &-distp(0,1+ic)*fun0*cosphase(0)*rphase(i)**2   & ! from 2
   &-distp(0,1+ic)*fun*sinphase(i,0)*rphase(i)*rmu(i,0)*rmucos(i,0)   & ! from 2
   &-fun0*sinphase(i,0)*rphase(i)*hess3d(i,ic)    ! from 2
           enddo 
                       if(ic.eq.1) then 
                    do i=1,3 

       z(indorbp,indt+i)=z(indorbp,indt+i)                              &
     &-2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d*cosphase(0) ! 3 
                                                                        
   hess(i)=hess(i)+2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d &
     &*sinphase(i,0)*rphase(i) ! from  3 

            hess(i)=hess(i)+(-2.d0*fun0*cost1d                                 &
     &*(rmucos(i,0)**2-rmusin(i,0)**2)                          &
     &-2.d0*fun*rmu(i,0)**2*rmucos(i,0)**2*cost1d)*cosphase(0)   

                    enddo 

       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0) ! 4                  

       hess(3)=hess(3)                              &
     &-(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) ! from 4                  
                                                                        
       hess(3)=hess(3)+(6.d0*cost1d*fun0                            &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)**2                       &
     &-6.d0*rmu(3,0)*rmucos(3,0)*cost1d                    &
     &*fun0*4.d0*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)        &
     &+6.d0*rmu(3,0)**2*rmucos(3,0)**2*cost1d              &
     &*fun*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      
                                                                        
                       elseif(ic.eq.2) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

       hess(1)=hess(1)-(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

      hess(2)=hess(2)-(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       hess(1)=hess(1)+(2.d0*rmu(1,0)**2                             &
     &*fun*cost2d*rmucos(1,0)**2                                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &+2.d0*fun0*cost2d*(rmucos(1,0)**2-rmusin(1,0)**2)         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &-2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*4.d0*PI/cellscale(1)*rmusin(1,0)*rmucos(1,0))*cosphase(0)             
                                                                        
       hess(2)=hess(2)+(-2.d0*rmu(2,0)**2                             &
     &*fun*cost2d*rmucos(2,0)**2                                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &-2.d0*fun0*cost2d*(rmucos(2,0)**2-rmusin(2,0)**2)         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &+2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*4.d0*PI/cellscale(2)*rmusin(2,0)*rmucos(2,0))*cosphase(0)             
                                                                        
                       elseif(ic.eq.3) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

        hess(1)=hess(1)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)           &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
        hess(2)=hess(2)+(-rmu(1,0)*cost3d*rmucos(1,0)             &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(1)=hess(1)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
                       elseif(ic.eq.4) then 
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)        &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                          

        hess(3)=hess(3)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)       &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3)                          
        hess(2)=hess(2)+(-rmu(3,0)*cost3d*rmucos(3,0)           &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                                                        
                       elseif(ic.eq.5) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)              

        hess(1)=hess(1)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      

        hess(3)=hess(3)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)      &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) 
                                                                        
       hess(1)=hess(1)+(-rmu(3,0)*cost3d*rmucos(3,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(1,0)*cost3d*rmucos(1,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                !endif for ic                           
                       endif 
                                                                        
!     z_xyz(indorbp,:)=hess(:) 
      z(indorbp,indt+4)=sum(hess(:)) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+5 
         indorb=indorbp  



       case(42) 
! derivative of the orbital 30 d orbitals                               
! d R(r) /d z = d c*exp(-z*r) / dz = c*exp(-z r)*(-r+ d c / d z)      
! each gaussian term is normalized                                      
                                                                        
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!         c=1.d0/(2.d0**3*3.d0)*dsqrt(1.d0/pi)*(2.d0*dd1)**(7.d0/2.d0)
         c=dd1**3.5d0*0.26596152026762178d0
!         endif 
         c0=-c 
         c1=3.5d0*c/dd1                                                                  
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k))
         distp(k,7)=distp(k,1)*(c0*r(k)+c1) 
         enddo 
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=(3.d0*(rmu(3,i)*rmucos(3,i))**2                &
     &-r(i)**2)*cost1d                                              
                           ! lz=0                                       
                                                                        
      distp(i,3)=((rmu(1,i)*rmucos(1,i))**2                     &
     &-(rmu(2,i)*rmucos(2,i))**2)*cost2d                        
                                                 ! lz=+/-2              
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d                           
                                               ! lz=+/-2                
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d                           
                                               ! lz=+/-1                
          enddo 
                                                                        
                                                                        
         do ic=1,5 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,7)*distp(k,1+ic)*cosphase(k)                        
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            dd1=dd(indparp) 
            fun0=distp(0,7)
            fun=(-dd1*distp(0,1)*(c0*r(0)+c1)+distp(0,1)*c0)/r(0)
            fun2=(-fun+(dd1**2*distp(0,1)*(c0*r(0)+c1) &
     & -dd1*distp(0,1)*c0-dd1*distp(0,1)*c0))/r(0)**2
            
       hess3d=0.d0 
                                                                        
       hess3d(:,1)=-2.d0*rmu(:,0)*rmucos(:,0)*cost1d 
                                                                        
       hess3d(3,1)=hess3d(3,1)                                          &
     &+6.d0*rmu(3,0)*rmucos(3,0)                           &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)*cost1d                   
                                                                        
       hess3d(1,2)=2.d0*rmu(1,0)*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)*cost2d                   
                                                                        
       hess3d(2,2)=-2.d0*rmu(2,0)*rmucos(2,0)                   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)*cost2d                   
                                                                        
                                                        ! lz=+/-2       
       hess3d(2,3)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-2       
       hess3d(1,3)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(2,4)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,4)=rmu(2,0)*rmucos(2,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(1,5)=rmu(3,0)*rmucos(3,0)*cost3d                  &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          
                                                                        
                                                        ! lz=+/-1       
       hess3d(3,5)=rmu(1,0)*rmucos(1,0)*cost3d                  &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)                          
                                                                        
               indorbp=indorb 
               do ic=1,5 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                  

                   do i=1,3 
      z(indorbp,indt+i)=distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)   & ! 1)
     &*cosphase(0)-distp(0,1+ic)*fun0*sinphase(i,0)*rphase(i)  ! 2)
                                                                        
          hess(i)=((fun2*(rmu(i,0)*rmucos(i,0))**2               &
     &+fun*(rmucos(i,0)**2-rmusin(i,0)**2))*distp(0,1+ic)       &
     &+hess3d(i,ic)*fun*rmu(i,0)*rmucos(i,0))*cosphase(0)             

                  !enddo for i                                          
   hess(i)=hess(i)-distp(0,1+ic)*fun*rmu(i,0)*rmucos(i,0)  & ! from 1
   &*sinphase(i,0)*rphase(i)  &
   &-distp(0,1+ic)*fun0*cosphase(0)*rphase(i)**2   & ! from 2
   &-distp(0,1+ic)*fun*sinphase(i,0)*rphase(i)*rmu(i,0)*rmucos(i,0)   & ! from 2
   &-fun0*sinphase(i,0)*rphase(i)*hess3d(i,ic)    ! from 2
           enddo 
                       if(ic.eq.1) then 
                    do i=1,3 

       z(indorbp,indt+i)=z(indorbp,indt+i)                              &
     &-2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d*cosphase(0) ! 3 
                                                                        
   hess(i)=hess(i)+2.d0*rmu(i,0)*rmucos(i,0)*fun0*cost1d &
     &*sinphase(i,0)*rphase(i) ! from  3 

            hess(i)=hess(i)+(-2.d0*fun0*cost1d                                 &
     &*(rmucos(i,0)**2-rmusin(i,0)**2)                          &
     &-2.d0*fun*rmu(i,0)**2*rmucos(i,0)**2*cost1d)*cosphase(0)   

                    enddo 

       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0) ! 4                  

       hess(3)=hess(3)                              &
     &-(6.d0*rmu(3,0)*rmucos(3,0)*cost1d                   &
     &*fun0*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) ! from 4                  
                                                                        
       hess(3)=hess(3)+(6.d0*cost1d*fun0                            &
     &*(rmucos(3,0)**2-rmusin(3,0)**2)**2                       &
     &-6.d0*rmu(3,0)*rmucos(3,0)*cost1d                    &
     &*fun0*4.d0*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)        &
     &+6.d0*rmu(3,0)**2*rmucos(3,0)**2*cost1d              &
     &*fun*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      
                                                                        
                       elseif(ic.eq.2) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

       hess(1)=hess(1)-(2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

      hess(2)=hess(2)-(-2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)   &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       hess(1)=hess(1)+(2.d0*rmu(1,0)**2                             &
     &*fun*cost2d*rmucos(1,0)**2                                    &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &+2.d0*fun0*cost2d*(rmucos(1,0)**2-rmusin(1,0)**2)         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2)                          &
     &-2.d0*rmu(1,0)*fun0*cost2d*rmucos(1,0)                    &
     &*4.d0*PI/cellscale(1)*rmusin(1,0)*rmucos(1,0))*cosphase(0)             
                                                                        
       hess(2)=hess(2)+(-2.d0*rmu(2,0)**2                             &
     &*fun*cost2d*rmucos(2,0)**2                                    &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &-2.d0*fun0*cost2d*(rmucos(2,0)**2-rmusin(2,0)**2)         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2)                          &
     &+2.d0*rmu(2,0)*fun0*cost2d*rmucos(2,0)                    &
     &*4.d0*PI/cellscale(2)*rmusin(2,0)*rmucos(2,0))*cosphase(0)             
                                                                        
                       elseif(ic.eq.3) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)                          

        hess(1)=hess(1)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)           &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
        hess(2)=hess(2)+(-rmu(1,0)*cost3d*rmucos(1,0)             &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(1)=hess(1)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
                       elseif(ic.eq.4) then 
       z(indorbp,indt+2)=z(indorbp,indt+2)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*cosphase(0)                          

       hess(2)=hess(2)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)        &
     &*(rmucos(2,0)**2-rmusin(2,0)**2))*sinphase(2,0)*rphase(2) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(2,0)*fun0*cost3d*rmucos(2,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                          

        hess(3)=hess(3)-(rmu(2,0)*fun0*cost3d*rmucos(2,0)       &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3)                          
        hess(2)=hess(2)+(-rmu(3,0)*cost3d*rmucos(3,0)           &
     &*(4.d0*PI/cellscale(2)*fun0*rmusin(2,0)*rmucos(2,0)       &
     &+fun*(rmusin(2,0)**2-rmucos(2,0)**2)                      &
     &*rmu(2,0)*rmucos(2,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(2,0)*cost3d*rmucos(2,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                                                        
                       elseif(ic.eq.5) then 
       z(indorbp,indt+1)=z(indorbp,indt+1)                              &
     &+(rmu(3,0)*fun0*cost3d*rmucos(3,0)                         &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*cosphase(0)              

        hess(1)=hess(1)-(rmu(3,0)*fun0*cost3d*rmucos(3,0)     &
     &*(rmucos(1,0)**2-rmusin(1,0)**2))*sinphase(1,0)*rphase(1) 
                                                                        
       z(indorbp,indt+3)=z(indorbp,indt+3)                              &
     &+(rmu(1,0)*fun0*cost3d*rmucos(1,0)                         &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*cosphase(0)                      

        hess(3)=hess(3)-(rmu(1,0)*fun0*cost3d*rmucos(1,0)      &
     &*(rmucos(3,0)**2-rmusin(3,0)**2))*sinphase(3,0)*rphase(3) 
                                                                        
       hess(1)=hess(1)+(-rmu(3,0)*cost3d*rmucos(3,0)              &
     &*(4.d0*PI/cellscale(1)*fun0*rmusin(1,0)*rmucos(1,0)       &
     &+fun*(rmusin(1,0)**2-rmucos(1,0)**2)                      &
     &*rmu(1,0)*rmucos(1,0)))*cosphase(0)                                    
                                                                        
       hess(3)=hess(3)+(-rmu(1,0)*cost3d*rmucos(1,0)              &
     &*(4.d0*PI/cellscale(3)*fun0*rmusin(3,0)*rmucos(3,0)       &
     &+fun*(rmusin(3,0)**2-rmucos(3,0)**2)                      &
     &*rmu(3,0)*rmucos(3,0)))*cosphase(0)                                    
                                                                        
                                !endif for ic    

                       endif 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
      z(indorbp,indt+4)=sum(hess(:)) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+5 
         indorb=indorbp    



      case(48) 
! f single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R
! Y_i are the stretched coordinates for the angular part A

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(Y_i(x_i)) \phi(x_i)
! R depends on X_i, A depends on Y_i, the phase \phi depends on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! where  d^2/dx^2_k F(X_i(x_i)) = (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

! X_i = rmu
! Y_i = rmu * rmucos
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin
! d/dx_k Y_i  =  2 * rmucos**2  - 1 
! d^2/dx^2_k Y_i = - 4*Pi/L * rmusin * rmucos

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!         c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0) 
          c=dd1**2.25d0*1.47215808929909374563d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

      rp0=(rmu(1,i)*rmucos(1,i))**2+(rmu(2,i)*rmucos(2,i))**2 &
     &   +(rmu(3,i)*rmucos(3,i))**2

!     proposal  SSS
!     rp0=r(k)**2


      distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
     &  *(5.d0*(rmu(3,i)*rmucos(3,i))**2-3.d0*rp0)

!     proposal SSS
!     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
!    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)

                                                          ! lz=0        
      distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
     &  *(5.d0*(rmu(3,i)*rmucos(3,i))**2-rp0)

!     proposal SSS
!     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
!    &  *(5.d0*rmu(3,i)**2-rp0)

                                                      ! lz=+/-1         
      distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
     &  *(5.d0*(rmu(3,i)*rmucos(3,i))**2-rp0)                    

!     proposal SSS
!     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
!    &  *(5.d0*rmu(3,i)**2-rp0)

                                                      ! lz=+/-1         
      distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
     &  *((rmu(1,i)*rmucos(1,i))**2-(rmu(2,i)*rmucos(2,i))**2)                     

!     proposal SSS
!     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
!    &  *(rmu(1,i)**2-rmu(2,i)**2)

                                                    ! lz=+/-2           
      distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)           &
     &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             
!     proposal SSS --> the same this is the most singular orbital with 3 rmucos.



                                            ! lz=+/-2                   
      distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)                &
     &  *((rmu(1,i)*rmucos(1,i))**2-3.d0*(rmu(2,i)*rmucos(2,i))**2)                
                                                          ! lz=+/-3     
!     proposal SSS
!     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)                &
!    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 


      distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)                &
     &  *(3.d0*(rmu(1,i)*rmucos(1,i))**2-(rmu(2,i)*rmucos(2,i))**2)                
                                                          ! lz=+/-3     
!     proposal SSS
!     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)                &
!    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  


          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=(rmu(1,0)*rmucos(1,0))**2+(rmu(2,0)*rmucos(2,0))**2 &
     &         +(rmu(3,0)*rmucos(3,0))**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              
! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     der4f(2)=0.d0
                     hess4f(2)=0.d0
                     hess4f(3)=0.d0

! contribution of the angular part
                     if(ic.eq.1) then 
! rp1 =   d/dx_k angular
! rp2 =   d^2/dx^2_k angular
                        rp1=-6.d0*cost1f*rmu(i,0)*rmucos(i,0)*rmu(3,0)*rmucos(3,0)
                        rp2=-6.d0*cost1f*rmu(3,0)*rmucos(3,0)
! der4f(2)  = der4f(2)  +  fun0 * d/dx_k angular
! hess4f(2) = hess4f(2) +  fun0 * d^2/dx^2_k angular
! hess4f(3) = hess4f(3) +  2 * fun * x_k * d/dx_k angular
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        if(i.eq.3) then
                           rp1=cost1f*(15.d0*(rmu(i,0)*rmucos(i,0))**2-3.d0*rp0)
                           rp2=cost1f*18.d0*rmu(i,0)*rmucos(i,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.2) then 
                        rp1=-2.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(1,0)*rmucos(1,0)
                        rp2=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        if(i.eq.1) then
                           rp1=cost2f*(5.d0*(rmu(3,0)*rmucos(3,0))**2-rp0)
                           rp2=-cost2f*4.d0*rmu(i,0)*rmucos(i,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.3) then 
                           rp1=10.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(1,0)*rmucos(1,0)
                           rp2=10.d0*cost2f*rmu(1,0)*rmucos(1,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.3) then
                        rp1=-2.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(2,0)*rmucos(2,0)
                        rp2=-2.d0*cost2f*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        if(i.eq.2) then
                           rp1=cost2f*(5.d0*(rmu(3,0)*rmucos(3,0))**2-rp0)
                           rp2=-cost2f*4.d0*rmu(i,0)*rmucos(i,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.3) then
                           rp1=10.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(2,0)*rmucos(2,0)
                           rp2=10.d0*cost2f*rmu(2,0)*rmucos(2,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.4) then 
                        if(i.eq.1) then
                           rp1=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           rp2=2.d0*cost3f*rmu(3,0)*rmucos(3,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           rp2=-2.d0*cost3f*rmu(3,0)*rmucos(3,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                           rp1=cost3f*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.5) then 
                        if(i.eq.1) then
                           rp1=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                           rp1=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.6) then 
                        if(i.eq.1) then
                           rp1=3.d0*cost4f*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                           rp2=6.d0*cost4f*rmu(1,0)*rmucos(1,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           rp2=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     else 
                        if(i.eq.1) then
                           rp1=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) 
                           rp2=6.d0*cost4f*rmu(2,0)*rmucos(2,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=3.d0*cost4f*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                           rp2=-6.d0*cost4f*rmu(2,0)*rmucos(2,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                        !endif for ic                           
                     endif
                       !enddo for i                                 

! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 

                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 


      case(49) 
! f single gaussian orbital derivative of 48
! R(r)= c*exp(-z r^2)*(9/4/z-r^2)  

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R
! Y_i are the stretched coordinates for the angular part A

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(Y_i(x_i)) \phi(x_i)
! R depends on X_i, A depends on Y_i, the phase \phi depends on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! where  d^2/dx^2_k F(X_i(x_i)) = (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

! X_i = rmu
! Y_i = rmu * rmucos
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin
! d/dx_k Y_i  =  2 * rmucos**2  - 1 
! d^2/dx^2_k Y_i = - 4*Pi/L * rmusin * rmucos

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
          c=dd1**2.25d0*1.47215808929909374563d0
!         c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0) 
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

      rp0=(rmu(1,i)*rmucos(1,i))**2+(rmu(2,i)*rmucos(2,i))**2 &
     &   +(rmu(3,i)*rmucos(3,i))**2

      distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
     &  *(5.d0*(rmu(3,i)*rmucos(3,i))**2-3.d0*rp0)
                                                          ! lz=0        
      distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
     &  *(5.d0*(rmu(3,i)*rmucos(3,i))**2-rp0)
                                                      ! lz=+/-1         
      distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
     &  *(5.d0*(rmu(3,i)*rmucos(3,i))**2-rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
     &  *((rmu(1,i)*rmucos(1,i))**2-(rmu(2,i)*rmucos(2,i))**2)                     
                                                    ! lz=+/-2           
      distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)           &
     &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             
                                            ! lz=+/-2                   
      distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)                &
     &  *((rmu(1,i)*rmucos(1,i))**2-3.d0*(rmu(2,i)*rmucos(2,i))**2)                
                                                          ! lz=+/-3     
      distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)                &
     &  *(3.d0*(rmu(1,i)*rmucos(1,i))**2-(rmu(2,i)*rmucos(2,i))**2)                
                                                          ! lz=+/-3     
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*(9.d0/4.d0/dd1-r(k)**2)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''

            fun0=distp(0,1)*(9.d0/4.d0/dd1-r(0)**2) 
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-13.d0/2.d0) 
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +19.d0*dd1*r(0)**2-13.d0/2.d0)  
            rp0=(rmu(1,0)*rmucos(1,0))**2+(rmu(2,0)*rmucos(2,0))**2 &
     &         +(rmu(3,0)*rmucos(3,0))**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              
! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     der4f(2)=0.d0
                     hess4f(2)=0.d0
                     hess4f(3)=0.d0

! contribution of the angular part
                     if(ic.eq.1) then 
! rp1 =   d/dx_k angular
! rp2 =   d^2/dx^2_k angular
                        rp1=-6.d0*cost1f*rmu(i,0)*rmucos(i,0)*rmu(3,0)*rmucos(3,0)
                        rp2=-6.d0*cost1f*rmu(3,0)*rmucos(3,0)
! der4f(2)  = der4f(2)  +  fun0 * d/dx_k angular
! hess4f(2) = hess4f(2) +  fun0 * d^2/dx^2_k angular
! hess4f(3) = hess4f(3) +  2 * fun * x_k * d/dx_k angular
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        if(i.eq.3) then
                           rp1=cost1f*(15.d0*(rmu(i,0)*rmucos(i,0))**2-3.d0*rp0)
                           rp2=cost1f*18.d0*rmu(i,0)*rmucos(i,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.2) then 
                        rp1=-2.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(1,0)*rmucos(1,0)
                        rp2=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        if(i.eq.1) then
                           rp1=cost2f*(5.d0*(rmu(3,0)*rmucos(3,0))**2-rp0)
                           rp2=-cost2f*4.d0*rmu(i,0)*rmucos(i,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.3) then 
                           rp1=10.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(1,0)*rmucos(1,0)
                           rp2=10.d0*cost2f*rmu(1,0)*rmucos(1,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.3) then
                        rp1=-2.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(2,0)*rmucos(2,0)
                        rp2=-2.d0*cost2f*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        if(i.eq.2) then
                           rp1=cost2f*(5.d0*(rmu(3,0)*rmucos(3,0))**2-rp0)
                           rp2=-cost2f*4.d0*rmu(i,0)*rmucos(i,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.3) then
                           rp1=10.d0*cost2f*rmu(i,0)*rmucos(i,0)*rmu(2,0)*rmucos(2,0)
                           rp2=10.d0*cost2f*rmu(2,0)*rmucos(2,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.4) then 
                        if(i.eq.1) then
                           rp1=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           rp2=2.d0*cost3f*rmu(3,0)*rmucos(3,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           rp2=-2.d0*cost3f*rmu(3,0)*rmucos(3,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                           rp1=cost3f*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.5) then 
                        if(i.eq.1) then
                           rp1=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                           rp1=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           rp2=0.d0
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.6) then 
                        if(i.eq.1) then
                           rp1=3.d0*cost4f*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                           rp2=6.d0*cost4f*rmu(1,0)*rmucos(1,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           rp2=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     else 
                        if(i.eq.1) then
                           rp1=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) 
                           rp2=6.d0*cost4f*rmu(2,0)*rmucos(2,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                           rp1=3.d0*cost4f*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                           rp2=-6.d0*cost4f*rmu(2,0)*rmucos(2,0)
                           der4f(2)=der4f(2)+fun0*rp1
                           hess4f(2)=hess4f(2)+fun0*rp2
                           hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                        !endif for ic                           
                     endif
                       !enddo for i                                 

! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 







      case(58) 
! f single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!         c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0) 
          c=dd1**2.25d0*1.47215808929909374563d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

     rp0=r(i)**2

     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)
                                                          ! lz=0        
     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
    &  *(rmu(1,i)**2-rmu(2,i)**2)

     distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)          &
    &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             

     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)               &
    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 
     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)               &
    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                           der4f(2)=-6.d0*cost1f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(1,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost1f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(2,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.3) then
                           der4f(2)=4.d0*cost1f*rmu(3,0)**2*rmucos(3,0)**2  &
                                +cost1f*(2.d0*rmucos(3,0)**2-1.d0)*(5.d0*rmu(3,0)**2-3.d0*rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=cost1f*fun0*(12.d0*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(3,0)**2-1.d0) &
                                -4.d0*Pi**2/cellscale(3)**2*rmu(3,0)*rmucos(3,0)*(5.d0*rmu(3,0)**2-3.d0*rp0))
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)**2*rmucos(1,0)**2 &
                                    +cost2f*(2.d0*rmucos(1,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-6.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
          &-4.d0*PI**2/cellscale(1)**2*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
        der4f(2)=-2.d0*cost2f*rmu(2,0)**2*rmucos(2,0)**2 &
                                   & +cost2f*(2.d0*rmucos(2,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
        hess4f(2)=-6.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
                 &-4.d0*PI**2/cellscale(2)**2*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(1,0)**2-1.d0)            
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=cost3f*(2.d0*rmucos(3,0)**2-1.d0)*(rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-4.d0*fun0*cost3f*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)*(rmu(1,0)**2-rmu(2,0)**2)
                        endif

                     elseif(ic.eq.5) then 

                        if(i.eq.1) then
   der4f(2)=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        elseif(i.eq.2) then
   der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        else
   der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        endif

                     elseif(ic.eq.6) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost4f*rmu(1,0)**2*rmucos(1,0)**2  &
                                   +cost4f*(2.d0*rmucos(1,0)**2-1.d0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
                &-4.d0*fun0*cost4f*PI**2/cellscale(1)**2*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                     else 

                        if(i.eq.1) then
                           der4f(2)=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
     der4f(2)=-2.d0*cost4f*rmu(2,0)**2*rmucos(2,0)**2  &
          &+cost4f*(2.d0*rmucos(2,0)**2-1.d0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
          &-4.d0*fun0*cost4f*PI**2/cellscale(2)**2*rmu(2,0)*rmucos(2,0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                        !endif for ic                           
                     endif
                       !enddo for i                                 


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 



      case(59) 
! f single gaussian orbital derivative of 58                                             
! radial R(r)= c*exp(-z r^2)*(9/4/z-r^2)                                                  
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
          c=dd1**2.25d0*1.47215808929909374563d0
!         c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0) 
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

     rp0=r(i)**2

     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)
                                                          ! lz=0        
     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
    &  *(rmu(1,i)**2-rmu(2,i)**2)

     distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)          &
    &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             

     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)               &
    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 
     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)               &
    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*(9.d0/4.d0/dd1-r(k)**2)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
            fun0=distp(0,1)*(9.d0/4.d0/dd1-r(0)**2) 
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-13.d0/2.d0) 
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +19.d0*dd1*r(0)**2-13.d0/2.d0)  
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                           der4f(2)=-6.d0*cost1f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(1,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost1f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(2,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.3) then
   der4f(2)=4.d0*cost1f*rmu(3,0)**2*rmucos(3,0)**2  &
          &+cost1f*(2.d0*rmucos(3,0)**2-1.d0)*(5.d0*rmu(3,0)**2-3.d0*rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=cost1f*fun0*(12.d0*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(3,0)**2-1.d0) &
         &-4.d0*Pi**2/cellscale(3)**2*rmu(3,0)*rmucos(3,0)*(5.d0*rmu(3,0)**2-3.d0*rp0))
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
   der4f(2)=-2.d0*cost2f*rmu(1,0)**2*rmucos(1,0)**2 &
               &+cost2f*(2.d0*rmucos(1,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
              &-4.d0*PI**2/cellscale(1)**2*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
    der4f(2)=-2.d0*cost2f*rmu(2,0)**2*rmucos(2,0)**2 &
            &+cost2f*(2.d0*rmucos(2,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
            &-4.d0*PI**2/cellscale(2)**2*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(1,0)**2-1.d0)            
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=cost3f*(2.d0*rmucos(3,0)**2-1.d0)*(rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-4.d0*fun0*cost3f*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0) &
             &*(rmu(1,0)**2-rmu(2,0)**2)
                        endif

                     elseif(ic.eq.5) then 

                        if(i.eq.1) then
     der4f(2)=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)&
           &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0) &
           &*rmusin(i,0)*rmucos(i,0)
                        elseif(i.eq.2) then
     der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
           &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
         &*rmusin(i,0)*rmucos(i,0)
                        else
      der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
              &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
         &*rmusin(i,0)*rmucos(i,0)
                        endif

                     elseif(ic.eq.6) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost4f*rmu(1,0)**2*rmucos(1,0)**2  &
                                   +cost4f*(2.d0*rmucos(1,0)**2-1.d0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
        &-4.d0*fun0*cost4f*PI**2/cellscale(1)**2*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                     else 

                        if(i.eq.1) then
                           der4f(2)=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
    der4f(2)=-2.d0*cost4f*rmu(2,0)**2*rmucos(2,0)**2  &
      &+cost4f*(2.d0*rmucos(2,0)**2-1.d0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
       &-4.d0*fun0*cost4f*PI**2/cellscale(2)**2*rmu(2,0)*rmucos(2,0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                        !endif for ic                           
                     endif
                       !enddo for i                                 


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 

      case(70) 
! f single exp  orbital                                             
! radial R(r)= exp(-alpha r)                                                 
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
          c=dd1**4.5d0*0.084104417400672d0
!         c=1.d0/dsqrt(10.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(9.d0/2.d0)/3.d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

     rp0=r(i)**2

     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)
                                                          ! lz=0        
     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
    &  *(rmu(1,i)**2-rmu(2,i)**2)

     distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)          &
    &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             

     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)               &
    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 
     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)               &
    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-dd1*distp(0,1)/r(0) 
            fun2=distp(0,1)*dd1**2
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                           der4f(2)=-6.d0*cost1f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(1,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost1f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(2,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.3) then
                           der4f(2)=4.d0*cost1f*rmu(3,0)**2*rmucos(3,0)**2  &
                                +cost1f*(2.d0*rmucos(3,0)**2-1.d0)*(5.d0*rmu(3,0)**2-3.d0*rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=cost1f*fun0*(12.d0*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(3,0)**2-1.d0) &
                                -4.d0*Pi**2/cellscale(3)**2*rmu(3,0)*rmucos(3,0)*(5.d0*rmu(3,0)**2-3.d0*rp0))
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)**2*rmucos(1,0)**2 &
                                    +cost2f*(2.d0*rmucos(1,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-6.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
          &-4.d0*PI**2/cellscale(1)**2*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
        der4f(2)=-2.d0*cost2f*rmu(2,0)**2*rmucos(2,0)**2 &
                                   & +cost2f*(2.d0*rmucos(2,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
        hess4f(2)=-6.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
                 &-4.d0*PI**2/cellscale(2)**2*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(1,0)**2-1.d0)            
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=cost3f*(2.d0*rmucos(3,0)**2-1.d0)*(rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-4.d0*fun0*cost3f*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)*(rmu(1,0)**2-rmu(2,0)**2)
                        endif

                     elseif(ic.eq.5) then 

                        if(i.eq.1) then
   der4f(2)=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        elseif(i.eq.2) then
   der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        else
   der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        endif

                     elseif(ic.eq.6) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost4f*rmu(1,0)**2*rmucos(1,0)**2  &
                                   +cost4f*(2.d0*rmucos(1,0)**2-1.d0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
                &-4.d0*fun0*cost4f*PI**2/cellscale(1)**2*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                     else 

                        if(i.eq.1) then
                           der4f(2)=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
     der4f(2)=-2.d0*cost4f*rmu(2,0)**2*rmucos(2,0)**2  &
          &+cost4f*(2.d0*rmucos(2,0)**2-1.d0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
          &-4.d0*fun0*cost4f*PI**2/cellscale(2)**2*rmu(2,0)*rmucos(2,0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                        !endif for ic                           
                     endif
                       !enddo for i                                 


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 



      case(71) 
! f single gaussian orbital derivative of 70                                             
! radial R(r)= c*exp(-z r)*(9/2/z-r)                                                  
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!         c=1.d0/dsqrt(10.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(9.d0/2.d0)/3.d0
          c=dd1**4.5d0*0.084104417400672d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

     rp0=r(i)**2

     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)
                                                          ! lz=0        
     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
    &  *(rmu(1,i)**2-rmu(2,i)**2)

     distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)          &
    &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             

     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)               &
    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 
     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)               &
    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
           z(indorbp,k)=distp(k,1)*(4.5d0/dd1-r(k))*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
            fun0=distp(0,1)*(4.5d0/dd1-r(0)) 
            fun=distp(0,1)*(-5.5d0/r(0)+dd1) 
            fun2=distp(0,1)*(6.5d0*dd1-r(0)*dd1**2)  
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                           der4f(2)=-6.d0*cost1f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(1,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost1f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(2,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.3) then
   der4f(2)=4.d0*cost1f*rmu(3,0)**2*rmucos(3,0)**2  &
          &+cost1f*(2.d0*rmucos(3,0)**2-1.d0)*(5.d0*rmu(3,0)**2-3.d0*rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=cost1f*fun0*(12.d0*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(3,0)**2-1.d0) &
         &-4.d0*Pi**2/cellscale(3)**2*rmu(3,0)*rmucos(3,0)*(5.d0*rmu(3,0)**2-3.d0*rp0))
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
   der4f(2)=-2.d0*cost2f*rmu(1,0)**2*rmucos(1,0)**2 &
               &+cost2f*(2.d0*rmucos(1,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
              &-4.d0*PI**2/cellscale(1)**2*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
    der4f(2)=-2.d0*cost2f*rmu(2,0)**2*rmucos(2,0)**2 &
            &+cost2f*(2.d0*rmucos(2,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
            &-4.d0*PI**2/cellscale(2)**2*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(1,0)**2-1.d0)            
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=cost3f*(2.d0*rmucos(3,0)**2-1.d0)*(rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-4.d0*fun0*cost3f*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0) &
             &*(rmu(1,0)**2-rmu(2,0)**2)
                        endif

                     elseif(ic.eq.5) then 

                        if(i.eq.1) then
     der4f(2)=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)&
           &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0) &
           &*rmusin(i,0)*rmucos(i,0)
                        elseif(i.eq.2) then
     der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
           &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
         &*rmusin(i,0)*rmucos(i,0)
                        else
      der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
              &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
         &*rmusin(i,0)*rmucos(i,0)
                        endif

                     elseif(ic.eq.6) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost4f*rmu(1,0)**2*rmucos(1,0)**2  &
                                   +cost4f*(2.d0*rmucos(1,0)**2-1.d0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
        &-4.d0*fun0*cost4f*PI**2/cellscale(1)**2*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                     else 

                        if(i.eq.1) then
                           der4f(2)=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
    der4f(2)=-2.d0*cost4f*rmu(2,0)**2*rmucos(2,0)**2  &
      &+cost4f*(2.d0*rmucos(2,0)**2-1.d0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
       &-4.d0*fun0*cost4f*PI**2/cellscale(2)**2*rmu(2,0)*rmucos(2,0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                        !endif for ic                           
                     endif
                       !enddo for i                                 


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(9.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 



      case(51) 
! g single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R
! Y_i are the stretched coordinates for the angular part A

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(Y_i(x_i)) \phi(x_i)
! R depends on X_i, A depends on Y_i, the phase \phi depends on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! where  d^2/dx^2_k F(X_i(x_i)) = (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

! X_i = rmu
! Y_i = rmu * rmucos
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin
! d/dx_k Y_i  =  2 * rmucos**2  - 1 
! d^2/dx^2_k Y_i = - 4*Pi/L * rmusin * rmucos

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!          if(iflagnorm.gt.2) then
! overall normalization                                                 
!          c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)  
           c=dd1**2.75d0*1.11284691281640568826d0
!          endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

      rp0=(rmu(1,i)*rmucos(1,i))**2+(rmu(2,i)*rmucos(2,i))**2 &
     &   +(rmu(3,i)*rmucos(3,i))**2


      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4*rmucos(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*rmucos(3,i)**2*rp0+3.d0*rp0**2)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2*rmucos(1,i)**2-rmu(2,i)**2*rmucos(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-rp0)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                  &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-rp0)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(rmu(1,i)**2*rmucos(1,i)**2-3.0*rmu(2,i)**2*rmucos(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(3.d0*rmu(1,i)**2*rmucos(1,i)**2-rmu(2,i)**2*rmucos(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4*rmucos(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmucos(1,i)**2*rmu(2,i)**2*rmucos(2,i)**2+rmu(2,i)**4*rmucos(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                 &
     &   *(rmu(1,i)**2*rmucos(1,i)**2-rmu(2,i)**2*rmucos(2,i)**2)      
                                                      ! lz=+/-4

          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,9
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=(rmu(1,0)*rmucos(1,0))**2+(rmu(2,0)*rmucos(2,0))**2 &
     &         +(rmu(3,0)*rmucos(3,0))**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,9
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              
! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     der4f(2)=0.d0
                     hess4f(2)=0.d0
                     hess4f(3)=0.d0

! contribution of the angular part
                     if(ic.eq.1) then 
! rp1 =   d/dx_k angular
! rp2 =   d^2/dx^2_k angular
                        if(i.eq.1) then
                        rp1=cost1g*(-60.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)**2*rmucos(3,0)**2  &
     &                     +12.d0*rmu(1,0)*rmucos(1,0)*rp0)
                        rp2=cost1g*12.d0*(3.d0*(rmu(1,0)*rmucos(1,0))**2   &
     &                     +(rmu(2,0)*rmucos(2,0))**2-4.d0*rmu(3,0)**2*rmucos(3,0)**2)
! der4f(2)  = der4f(2)  +  fun0 * d/dx_k angular
! hess4f(2) = hess4f(2) +  fun0 * d^2/dx^2_k angular
! hess4f(3) = hess4f(3) +  2 * fun * x_k * d/dx_k angular
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then 
                        rp1=cost1g*(-60.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)**2*rmucos(3,0)**2  &
     &                      +12.d0*rmu(2,0)*rmucos(2,0)*rp0)
                        rp2=cost1g*12.d0*((rmu(1,0)*rmucos(1,0))**2+3.d0*(rmu(2,0)*rmucos(2,0))**2  &
     &                      -4.d0*rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost1g*(80.d0*rmu(3,0)**3*rmucos(3,0)**3-48.d0*rmu(3,0)*rmucos(3,0)*rp0)
                        rp2=-cost1g*48.d0*((rmu(1,0)*rmucos(1,0))**2+(rmu(2,0)*rmucos(2,0))**2  &
     &                      -2.d0*rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.2) then
                        if(i.eq.1) then
                        rp1=cost2g*(-9.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(3,0)*rmucos(3,0)   &
     &                      -3.d0*(rmu(2,0)*rmucos(2,0))**2*rmu(3,0)*rmucos(3,0)   &
     &                      +4.d0*rmu(3,0)**3*rmucos(3,0)**3)
                        rp2=-18.d0*cost2g*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost2g*(-6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)) 
                        rp2=-6.d0*cost2g*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost2g*(-3.d0*rmu(1,0)*rmucos(1,0)*((rmu(1,0)*rmucos(1,0))**2   &
     &                      +(rmu(2,0)*rmucos(2,0))**2 &
     &                      -4.d0*rmu(3,0)**2*rmucos(3,0)**2))   
                        rp2=24.d0*cost2g*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.3) then
                        if(i.eq.1) then
                        rp1=cost2g*(-6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0))  
                        rp2=-6.d0*cost2g*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost2g*(-3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(3,0)*rmucos(3,0)  &
     &                      -9.d0*(rmu(2,0)*rmucos(2,0))**2*rmu(3,0)*rmucos(3,0)  &
     &                      +4.d0*rmu(3,0)**3*rmucos(3,0)**3)
                        rp2=-18.d0*cost2g*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost2g*(-3.d0*rmu(2,0)*rmucos(2,0)*((rmu(1,0)*rmucos(1,0))**2  &
     &                      +(rmu(2,0)*rmucos(2,0))**2-4.d0*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=24.d0*cost2g*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.4) then
                        if(i.eq.1) then
                        rp1=cost3g*(-4.d0*((rmu(1,0)*rmucos(1,0))**3  &
     &                      -3.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=-12.d0*cost3g*((rmu(1,0)*rmucos(1,0))**2-rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost3g*(4.d0*((rmu(2,0)*rmucos(2,0))**3 &
     &                      -3.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=12.d0*cost3g*((rmu(2,0)*rmucos(2,0))**2-rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost3g*(12.d0*((rmu(1,0)*rmucos(1,0))**2 &
     &                      -(rmu(2,0)*rmucos(2,0))**2)*rmu(3,0)*rmucos(3,0))
                        rp2=12.d0*cost3g*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.5) then
                        if(i.eq.1) then
                        rp1=cost3g*(-2.d0*rmu(2,0)*rmucos(2,0)*(3.d0*(rmu(1,0)*rmucos(1,0))**2  &
     &                       +(rmu(2,0)*rmucos(2,0))**2-6.d0*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=-12.d0*cost3g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)    
                        elseif(i.eq.2) then
                        rp1=cost3g*(-2.d0*rmu(1,0)*rmucos(1,0)*((rmu(1,0)*rmucos(1,0))**2  &
     &                       +3.d0*(rmu(2,0)*rmucos(2,0))**2-6.d0*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=-12.d0*cost3g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)    
                        else
                        rp1=cost3g*24.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        rp2=cost3g*24.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)    
                        endif
                     elseif(ic.eq.6) then
                        if(i.eq.1) then
                        rp1=cost4g*3.d0*((rmu(1,0)*rmucos(1,0))**2  &
     &                     -(rmu(2,0)*rmucos(2,0))**2)*rmu(3,0)*rmucos(3,0)
                        rp2=cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)                            
                        elseif(i.eq.2) then
                        rp1=-cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        rp2=-cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost4g*((rmu(1,0)*rmucos(1,0))**3  &
     &                      -3.d0*rmu(1,0)*rmucos(1,0)*(rmu(2,0)*rmucos(2,0))**2)
                        rp2=0.d0
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.7) then
                        if(i.eq.1) then
                        rp1=cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        rp2=cost4g*6.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost4g*3.d0*((rmu(1,0)*rmucos(1,0))**2 &
     &                      -(rmu(2,0)*rmucos(2,0))**2)*rmu(3,0)*rmucos(3,0)
                        rp2=-cost4g*6.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost4g*(3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(2,0)*rmucos(2,0)  &
     &                       -(rmu(2,0)*rmucos(2,0))**3)
                        rp2=0.d0
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.8) then
                        if(i.eq.1) then
                        rp1=cost5g*4.d0*((rmu(1,0)*rmucos(1,0))**3  &
     &                        -3.d0*rmu(1,0)*rmucos(1,0)*(rmu(2,0)*rmucos(2,0))**2)
                        rp2=cost5g*12.d0*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost5g*4.d0*((rmu(2,0)*rmucos(2,0))**3   &
     &                        -3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(2,0)*rmucos(2,0))
                        rp2=cost5g*12.d0*((rmu(2,0)*rmucos(2,0))**2-(rmu(1,0)*rmucos(1,0))**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     else
                        if(i.eq.1) then
                        rp1=cost5g*4.d0*(3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(2,0)*rmucos(2,0) &
     &                         -(rmu(2,0)*rmucos(2,0))**3)
                        rp2=24.d0*cost5g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost5g*4.d0*((rmu(1,0)*rmucos(1,0))**3 &
     &                         -3.d0*rmu(1,0)*rmucos(1,0)*(rmu(2,0)*rmucos(2,0))**2)   
                        rp2=-24.d0*cost5g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     endif

! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 

                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+9 
         indorb=indorbp 



      case(52) 
! g single gaussian orbital 
! derivative of 51                                            
! normalized 

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R
! Y_i are the stretched coordinates for the angular part A

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(Y_i(x_i)) \phi(x_i)
! R depends on X_i, A depends on Y_i, the phase \phi depends on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! where  d^2/dx^2_k F(X_i(x_i)) = (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

! X_i = rmu
! Y_i = rmu * rmucos
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin
! d/dx_k Y_i  =  2 * rmucos**2  - 1 
! d^2/dx^2_k Y_i = - 4*Pi/L * rmusin * rmucos

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!          if(iflagnorm.gt.2) then
! overall normalization                                                 
           c=dd1**2.75d0*1.11284691281640568826d0
!          c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)  
!          endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

      rp0=(rmu(1,i)*rmucos(1,i))**2+(rmu(2,i)*rmucos(2,i))**2 &
     &   +(rmu(3,i)*rmucos(3,i))**2


      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4*rmucos(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*rmucos(3,i)**2*rp0+3.d0*rp0**2)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2*rmucos(1,i)**2-rmu(2,i)**2*rmucos(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-rp0)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                  &
     &           *(7.d0*rmu(3,i)**2*rmucos(3,i)**2-rp0)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(rmu(1,i)**2*rmucos(1,i)**2-3.0*rmu(2,i)**2*rmucos(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)                       &
     &           *(3.d0*rmu(1,i)**2*rmucos(1,i)**2-rmu(2,i)**2*rmucos(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4*rmucos(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmucos(1,i)**2*rmu(2,i)**2*rmucos(2,i)**2+rmu(2,i)**4*rmucos(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                 &
     &   *(rmu(1,i)**2*rmucos(1,i)**2-rmu(2,i)**2*rmucos(2,i)**2)      
                                                      ! lz=+/-4

          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,9
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*(11.d0/4.d0/dd1-r(k)**2)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
            fun0=distp(0,1)*(11.d0/4.d0/dd1-r(0)**2) 
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-15.d0/2.d0) 
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +21.d0*dd1*r(0)**2-15.d0/2.d0)  

            rp0=(rmu(1,0)*rmucos(1,0))**2+(rmu(2,0)*rmucos(2,0))**2 &
     &         +(rmu(3,0)*rmucos(3,0))**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,9
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              
! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     der4f(2)=0.d0
                     hess4f(2)=0.d0
                     hess4f(3)=0.d0

! contribution of the angular part
                     if(ic.eq.1) then 
! rp1 =   d/dx_k angular
! rp2 =   d^2/dx^2_k angular
                        if(i.eq.1) then
                        rp1=cost1g*(-60.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)**2*rmucos(3,0)**2  &
     &                     +12.d0*rmu(1,0)*rmucos(1,0)*rp0)
                        rp2=cost1g*12.d0*(3.d0*(rmu(1,0)*rmucos(1,0))**2   &
     &                     +(rmu(2,0)*rmucos(2,0))**2-4.d0*rmu(3,0)**2*rmucos(3,0)**2)
! der4f(2)  = der4f(2)  +  fun0 * d/dx_k angular
! hess4f(2) = hess4f(2) +  fun0 * d^2/dx^2_k angular
! hess4f(3) = hess4f(3) +  2 * fun * x_k * d/dx_k angular
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then 
                        rp1=cost1g*(-60.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)**2*rmucos(3,0)**2  &
     &                      +12.d0*rmu(2,0)*rmucos(2,0)*rp0)
                        rp2=cost1g*12.d0*((rmu(1,0)*rmucos(1,0))**2+3.d0*(rmu(2,0)*rmucos(2,0))**2  &
     &                      -4.d0*rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost1g*(80.d0*rmu(3,0)**3*rmucos(3,0)**3-48.d0*rmu(3,0)*rmucos(3,0)*rp0)
                        rp2=-cost1g*48.d0*((rmu(1,0)*rmucos(1,0))**2+(rmu(2,0)*rmucos(2,0))**2  &
     &                      -2.d0*rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.2) then
                        if(i.eq.1) then
                        rp1=cost2g*(-9.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(3,0)*rmucos(3,0)   &
     &                      -3.d0*(rmu(2,0)*rmucos(2,0))**2*rmu(3,0)*rmucos(3,0)   &
     &                      +4.d0*rmu(3,0)**3*rmucos(3,0)**3)
                        rp2=-18.d0*cost2g*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost2g*(-6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)) 
                        rp2=-6.d0*cost2g*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost2g*(-3.d0*rmu(1,0)*rmucos(1,0)*((rmu(1,0)*rmucos(1,0))**2   &
     &                      +(rmu(2,0)*rmucos(2,0))**2 &
     &                      -4.d0*rmu(3,0)**2*rmucos(3,0)**2))   
                        rp2=24.d0*cost2g*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.3) then
                        if(i.eq.1) then
                        rp1=cost2g*(-6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0))  
                        rp2=-6.d0*cost2g*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost2g*(-3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(3,0)*rmucos(3,0)  &
     &                      -9.d0*(rmu(2,0)*rmucos(2,0))**2*rmu(3,0)*rmucos(3,0)  &
     &                      +4.d0*rmu(3,0)**3*rmucos(3,0)**3)
                        rp2=-18.d0*cost2g*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost2g*(-3.d0*rmu(2,0)*rmucos(2,0)*((rmu(1,0)*rmucos(1,0))**2  &
     &                      +(rmu(2,0)*rmucos(2,0))**2-4.d0*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=24.d0*cost2g*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.4) then
                        if(i.eq.1) then
                        rp1=cost3g*(-4.d0*((rmu(1,0)*rmucos(1,0))**3  &
     &                      -3.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=-12.d0*cost3g*((rmu(1,0)*rmucos(1,0))**2-rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost3g*(4.d0*((rmu(2,0)*rmucos(2,0))**3 &
     &                      -3.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=12.d0*cost3g*((rmu(2,0)*rmucos(2,0))**2-rmu(3,0)**2*rmucos(3,0)**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost3g*(12.d0*((rmu(1,0)*rmucos(1,0))**2 &
     &                      -(rmu(2,0)*rmucos(2,0))**2)*rmu(3,0)*rmucos(3,0))
                        rp2=12.d0*cost3g*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.5) then
                        if(i.eq.1) then
                        rp1=cost3g*(-2.d0*rmu(2,0)*rmucos(2,0)*(3.d0*(rmu(1,0)*rmucos(1,0))**2  &
     &                       +(rmu(2,0)*rmucos(2,0))**2-6.d0*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=-12.d0*cost3g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)    
                        elseif(i.eq.2) then
                        rp1=cost3g*(-2.d0*rmu(1,0)*rmucos(1,0)*((rmu(1,0)*rmucos(1,0))**2  &
     &                       +3.d0*(rmu(2,0)*rmucos(2,0))**2-6.d0*rmu(3,0)**2*rmucos(3,0)**2))
                        rp2=-12.d0*cost3g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)    
                        else
                        rp1=cost3g*24.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        rp2=cost3g*24.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)    
                        endif
                     elseif(ic.eq.6) then
                        if(i.eq.1) then
                        rp1=cost4g*3.d0*((rmu(1,0)*rmucos(1,0))**2  &
     &                     -(rmu(2,0)*rmucos(2,0))**2)*rmu(3,0)*rmucos(3,0)
                        rp2=cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)                            
                        elseif(i.eq.2) then
                        rp1=-cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        rp2=-cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost4g*((rmu(1,0)*rmucos(1,0))**3  &
     &                      -3.d0*rmu(1,0)*rmucos(1,0)*(rmu(2,0)*rmucos(2,0))**2)
                        rp2=0.d0
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.7) then
                        if(i.eq.1) then
                        rp1=cost4g*6.d0*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        rp2=cost4g*6.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost4g*3.d0*((rmu(1,0)*rmucos(1,0))**2 &
     &                      -(rmu(2,0)*rmucos(2,0))**2)*rmu(3,0)*rmucos(3,0)
                        rp2=-cost4g*6.d0*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        else
                        rp1=cost4g*(3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(2,0)*rmucos(2,0)  &
     &                       -(rmu(2,0)*rmucos(2,0))**3)
                        rp2=0.d0
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     elseif(ic.eq.8) then
                        if(i.eq.1) then
                        rp1=cost5g*4.d0*((rmu(1,0)*rmucos(1,0))**3  &
     &                        -3.d0*rmu(1,0)*rmucos(1,0)*(rmu(2,0)*rmucos(2,0))**2)
                        rp2=cost5g*12.d0*((rmu(1,0)*rmucos(1,0))**2-(rmu(2,0)*rmucos(2,0))**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost5g*4.d0*((rmu(2,0)*rmucos(2,0))**3   &
     &                        -3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(2,0)*rmucos(2,0))
                        rp2=cost5g*12.d0*((rmu(2,0)*rmucos(2,0))**2-(rmu(1,0)*rmucos(1,0))**2)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     else
                        if(i.eq.1) then
                        rp1=cost5g*4.d0*(3.d0*(rmu(1,0)*rmucos(1,0))**2*rmu(2,0)*rmucos(2,0) &
     &                         -(rmu(2,0)*rmucos(2,0))**3)
                        rp2=24.d0*cost5g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        elseif(i.eq.2) then
                        rp1=cost5g*4.d0*((rmu(1,0)*rmucos(1,0))**3 &
     &                         -3.d0*rmu(1,0)*rmucos(1,0)*(rmu(2,0)*rmucos(2,0))**2)   
                        rp2=-24.d0*cost5g*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                        der4f(2)=der4f(2)+fun0*rp1
                        hess4f(2)=hess4f(2)+fun0*rp2
                        hess4f(3)=hess4f(3)+2.d0*fun*rp1*rmu(i,0)
                        endif
                     endif

! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dY^2_k A)  ((d/dx_k Y_i) (x_i))^2
!                      + R  (d/dY_k A)  (d^2/dx^2_k Y_i) (x_i)
!                      + 2 (d/dX_k R)  (d/dY_k A)   (d/dx_k X_i) (x_i)  (d/dx_k Y_i) (x_i) 

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 

                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+9 
         indorb=indorbp 



      case(53) 
! g single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization 
!         c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)                                                  
          c=dd1**2.75d0*1.11284691281640568826d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*rp0+3.d0*rp0**2)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2-rp0)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)  &
     &           *(7.d0*rmu(3,i)**2-rp0)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(rmu(1,i)**2-3.0*rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)   &
     &   *(rmu(1,i)**2-rmu(2,i)**2)      
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,9
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,9 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                   der4f(2)=cost1g*12*rmu(1,0)*rmucos(1,0)*(-5*rmu(3,0)**2 + rp0) 
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(1,0)**2*rmucos(1,0)**2 -   &
     &        (rmucos(1,0)**2 - rmusin(1,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.2) then
                   der4f(2)=cost1g*12*rmu(2,0)*rmucos(2,0)*(-5*rmu(3,0)**2 + rp0)  
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(2,0)**2*rmucos(2,0)**2 - &
     &        (rmucos(2,0)**2 - rmusin(2,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.3) then
                   der4f(2)=cost1g*16*rmu(3,0)*rmucos(3,0)*(5*rmu(3,0)**2 - 3*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*16*(rmu(3,0)**2*(9*rmucos(3,0)**2 - &
     &         5*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*(-6*rmu(1,0)**2*rmucos(1,0)**2 + &
     &           (rmucos(1,0)**2 - rmusin(1,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*(-2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (9*rmu(1,0)*(rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &         2*Pi/cellscale(1)*rmusin(1,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.2) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &            rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*6*rmu(1,0)*rmu(3,0)*   &
     &            rmucos(1,0)*rmucos(3,0)*(-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then 
                   der4f(2)=cost2g*rmu(1,0)*rmucos(1,0)* &
     &           (rmu(3,0)**2*(15*rmucos(3,0)**2 - 7*rmusin(3,0)**2) + &
     &           3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*4*rmu(1,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) +  & 
     &           3*Pi/cellscale(3)*rmusin(3,0)*rp0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &             rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*6*rmu(2,0)*rmu(3,0)*    &
     &             rmucos(2,0)*rmucos(3,0)*(-rmucos(1,0)**2 + rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*    &
     &         (-6*rmu(2,0)**2*rmucos(2,0)**2 + (rmucos(2,0)**2 - & 
     &          rmusin(2,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*  &
     &         (9*rmu(2,0)*(rmucos(2,0)**2 - rmusin(2,0)**2) +   &
     &      2*Pi/cellscale(2)*rmusin(2,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.3) then
                   der4f(2)=cost2g*rmu(2,0)*rmucos(2,0)*(rmu(3,0)**2*(15*rmucos(3,0)**2 -   &
     &      7*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(4*rmu(2,0)*rmucos(2,0)*rmucos(3,0)* &
     &     (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) + 3*Pi/cellscale(3)*rmusin(3,0)*rp0))
                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                   der4f(2)=-cost3g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(3,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*fun0*cost3g*(rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &             3*rmu(3,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost3g*4*rmu(2,0)*(rmu(2,0)**2 - 3*rmu(3,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=4*fun0*cost3g*(rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) + &
     &             3*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)* &
     &             (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                        


                     elseif(ic.eq.5) then 


                        if(i.eq.1) then
                   der4f(2)=-2*cost3g*rmu(2,0)*rmucos(2,0)*((rmu(2,0)**2 - 6*rmu(3,0)**2)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2) +   &
     &             rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*(2*(rmu(2,0)**2 - & 
     &              6*rmu(3,0)**2)*rmusin(1,0)*pi/cellscale(1) +  & 
     &         rmu(1,0)*(-3*rmucos(1,0)**2 + 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-2*cost3g*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &         rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) +  &
     &         6*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*  &
     &         (2*(rmu(1,0)**2 - 6*rmu(3,0)**2)*rmusin(2,0)*pi/cellscale(2) +  & 
     &         rmu(2,0)*(-3*rmucos(2,0)**2 + 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*24*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*24*rmu(1,0)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                 


                     elseif(ic.eq.6) then 


                        if(i.eq.1) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)* &
     &           (rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) + &
     &           3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (6*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &          (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-6*cost4g*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=6*cost4g*fun0*rmu(1,0)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then
                   der4f(2)=cost4g*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(-4*rmu(1,0)*(rmu(1,0)**2 - &
     &          3*rmu(2,0)**2)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 



                     elseif(ic.eq.7) then


                        if(i.eq.1) then
                   der4f(2)=cost4g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*6*rmu(2,0)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)*(3*rmu(1,0)**2* &
     &            (rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &            rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (6*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + & 
     &            rmu(2,0)*(3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=-cost4g*(rmu(2,0)*(-3*rmu(1,0)**2 +  &
     &            rmu(2,0)**2)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*  &
     &            rmucos(2,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 

                     elseif(ic.eq.8) then


                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*4*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &            rmusin(1,0)**2) + 3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*cost5g*fun0*(3*rmu(1,0)**2*(rmucos(2,0)**2 - &
     &            rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                     else 

                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(2,0)*rmucos(2,0)*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &             rmusin(1,0)**2) + rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(8*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &       (2*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &       (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - &
     &             rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(-8*rmu(1,0)*rmucos(1,0)*rmucos(2,0)* &
     &        (2*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + rmu(2,0)*  &
     &       (3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+9
         indorb=indorbp 


      case(54) 
! g single gaussian orbital                                             
! derivative of 53
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization 
!         c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)  
          c=dd1**2.75d0*1.11284691281640568826d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*rp0+3.d0*rp0**2)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2-rp0)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)  &
     &           *(7.d0*rmu(3,i)**2-rp0)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(rmu(1,i)**2-3.0*rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)   &
     &   *(rmu(1,i)**2-rmu(2,i)**2)      
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,9
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*(11.d0/4.d0/dd1-r(k)**2)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                 
! fun2 = radial''
            fun0=distp(0,1)*(11.d0/4.d0/dd1-r(0)**2) 
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-15.d0/2.d0) 
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +21.d0*dd1*r(0)**2-15.d0/2.d0)  
                                                       
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,9 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                   der4f(2)=cost1g*12*rmu(1,0)*rmucos(1,0)*(-5*rmu(3,0)**2 + rp0) 
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(1,0)**2*rmucos(1,0)**2 -   &
     &        (rmucos(1,0)**2 - rmusin(1,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.2) then
                   der4f(2)=cost1g*12*rmu(2,0)*rmucos(2,0)*(-5*rmu(3,0)**2 + rp0)  
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(2,0)**2*rmucos(2,0)**2 - &
     &        (rmucos(2,0)**2 - rmusin(2,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.3) then
                   der4f(2)=cost1g*16*rmu(3,0)*rmucos(3,0)*(5*rmu(3,0)**2 - 3*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*16*(rmu(3,0)**2*(9*rmucos(3,0)**2 - &
     &         5*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*(-6*rmu(1,0)**2*rmucos(1,0)**2 + &
     &           (rmucos(1,0)**2 - rmusin(1,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*(-2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (9*rmu(1,0)*(rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &         2*Pi/cellscale(1)*rmusin(1,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.2) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &            rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*6*rmu(1,0)*rmu(3,0)*   &
     &            rmucos(1,0)*rmucos(3,0)*(-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then 
                   der4f(2)=cost2g*rmu(1,0)*rmucos(1,0)* &
     &           (rmu(3,0)**2*(15*rmucos(3,0)**2 - 7*rmusin(3,0)**2) + &
     &           3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*4*rmu(1,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) +  & 
     &           3*Pi/cellscale(3)*rmusin(3,0)*rp0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &             rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*6*rmu(2,0)*rmu(3,0)*    &
     &             rmucos(2,0)*rmucos(3,0)*(-rmucos(1,0)**2 + rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*    &
     &         (-6*rmu(2,0)**2*rmucos(2,0)**2 + (rmucos(2,0)**2 - & 
     &          rmusin(2,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*  &
     &         (9*rmu(2,0)*(rmucos(2,0)**2 - rmusin(2,0)**2) +   &
     &      2*Pi/cellscale(2)*rmusin(2,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.3) then
                   der4f(2)=cost2g*rmu(2,0)*rmucos(2,0)*(rmu(3,0)**2*(15*rmucos(3,0)**2 -   &
     &      7*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(4*rmu(2,0)*rmucos(2,0)*rmucos(3,0)* &
     &     (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) + 3*Pi/cellscale(3)*rmusin(3,0)*rp0))
                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                   der4f(2)=-cost3g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(3,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*fun0*cost3g*(rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &             3*rmu(3,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost3g*4*rmu(2,0)*(rmu(2,0)**2 - 3*rmu(3,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=4*fun0*cost3g*(rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) + &
     &             3*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)* &
     &             (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                        


                     elseif(ic.eq.5) then 


                        if(i.eq.1) then
                   der4f(2)=-2*cost3g*rmu(2,0)*rmucos(2,0)*((rmu(2,0)**2 - 6*rmu(3,0)**2)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2) +   &
     &             rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*(2*(rmu(2,0)**2 - & 
     &              6*rmu(3,0)**2)*rmusin(1,0)*pi/cellscale(1) +  & 
     &         rmu(1,0)*(-3*rmucos(1,0)**2 + 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-2*cost3g*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &         rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) +  &
     &         6*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*  &
     &         (2*(rmu(1,0)**2 - 6*rmu(3,0)**2)*rmusin(2,0)*pi/cellscale(2) +  & 
     &         rmu(2,0)*(-3*rmucos(2,0)**2 + 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*24*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*24*rmu(1,0)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                 


                     elseif(ic.eq.6) then 


                        if(i.eq.1) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)* &
     &           (rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) + &
     &           3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (6*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &          (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-6*cost4g*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=6*cost4g*fun0*rmu(1,0)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then
                   der4f(2)=cost4g*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(-4*rmu(1,0)*(rmu(1,0)**2 - &
     &          3*rmu(2,0)**2)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 



                     elseif(ic.eq.7) then


                        if(i.eq.1) then
                   der4f(2)=cost4g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*6*rmu(2,0)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)*(3*rmu(1,0)**2* &
     &            (rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &            rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (6*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + & 
     &            rmu(2,0)*(3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=-cost4g*(rmu(2,0)*(-3*rmu(1,0)**2 +  &
     &            rmu(2,0)**2)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*  &
     &            rmucos(2,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 

                     elseif(ic.eq.8) then


                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*4*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &            rmusin(1,0)**2) + 3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*cost5g*fun0*(3*rmu(1,0)**2*(rmucos(2,0)**2 - &
     &            rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                     else 

                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(2,0)*rmucos(2,0)*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &             rmusin(1,0)**2) + rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(8*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &       (2*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &       (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - &
     &             rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(-8*rmu(1,0)*rmucos(1,0)*rmucos(2,0)* &
     &        (2*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + rmu(2,0)*  &
     &       (3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+9
         indorb=indorbp 

      case(55) 
! g single Slater orbital                                             
! radial R(r)= exp(-alpha r)                                                 
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization 
!         c=1.d0/dsqrt(7.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(11.d0/2.d0)/3.d0/5.d0
          c=dd1**5.5d0*.020104801169736915d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*rp0+3.d0*rp0**2)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2-rp0)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)  &
     &           *(7.d0*rmu(3,i)**2-rp0)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(rmu(1,i)**2-3.0*rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)   &
     &   *(rmu(1,i)**2-rmu(2,i)**2)      
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,9
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1)
            fun=-dd1*distp(0,1)/r(0)
            fun2=distp(0,1)*dd1**2

            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,9 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                   der4f(2)=cost1g*12*rmu(1,0)*rmucos(1,0)*(-5*rmu(3,0)**2 + rp0) 
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(1,0)**2*rmucos(1,0)**2 -   &
     &        (rmucos(1,0)**2 - rmusin(1,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.2) then
                   der4f(2)=cost1g*12*rmu(2,0)*rmucos(2,0)*(-5*rmu(3,0)**2 + rp0)  
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(2,0)**2*rmucos(2,0)**2 - &
     &        (rmucos(2,0)**2 - rmusin(2,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.3) then
                   der4f(2)=cost1g*16*rmu(3,0)*rmucos(3,0)*(5*rmu(3,0)**2 - 3*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*16*(rmu(3,0)**2*(9*rmucos(3,0)**2 - &
     &         5*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*(-6*rmu(1,0)**2*rmucos(1,0)**2 + &
     &           (rmucos(1,0)**2 - rmusin(1,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*(-2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (9*rmu(1,0)*(rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &         2*Pi/cellscale(1)*rmusin(1,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.2) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &            rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*6*rmu(1,0)*rmu(3,0)*   &
     &            rmucos(1,0)*rmucos(3,0)*(-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then 
                   der4f(2)=cost2g*rmu(1,0)*rmucos(1,0)* &
     &           (rmu(3,0)**2*(15*rmucos(3,0)**2 - 7*rmusin(3,0)**2) + &
     &           3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*4*rmu(1,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) +  & 
     &           3*Pi/cellscale(3)*rmusin(3,0)*rp0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &             rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*6*rmu(2,0)*rmu(3,0)*    &
     &             rmucos(2,0)*rmucos(3,0)*(-rmucos(1,0)**2 + rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*    &
     &         (-6*rmu(2,0)**2*rmucos(2,0)**2 + (rmucos(2,0)**2 - & 
     &          rmusin(2,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*  &
     &         (9*rmu(2,0)*(rmucos(2,0)**2 - rmusin(2,0)**2) +   &
     &      2*Pi/cellscale(2)*rmusin(2,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.3) then
                   der4f(2)=cost2g*rmu(2,0)*rmucos(2,0)*(rmu(3,0)**2*(15*rmucos(3,0)**2 -   &
     &      7*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(4*rmu(2,0)*rmucos(2,0)*rmucos(3,0)* &
     &     (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) + 3*Pi/cellscale(3)*rmusin(3,0)*rp0))
                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                   der4f(2)=-cost3g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(3,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*fun0*cost3g*(rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &             3*rmu(3,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost3g*4*rmu(2,0)*(rmu(2,0)**2 - 3*rmu(3,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=4*fun0*cost3g*(rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) + &
     &             3*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)* &
     &             (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                        


                     elseif(ic.eq.5) then 


                        if(i.eq.1) then
                   der4f(2)=-2*cost3g*rmu(2,0)*rmucos(2,0)*((rmu(2,0)**2 - 6*rmu(3,0)**2)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2) +   &
     &             rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*(2*(rmu(2,0)**2 - & 
     &              6*rmu(3,0)**2)*rmusin(1,0)*pi/cellscale(1) +  & 
     &         rmu(1,0)*(-3*rmucos(1,0)**2 + 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-2*cost3g*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &         rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) +  &
     &         6*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*  &
     &         (2*(rmu(1,0)**2 - 6*rmu(3,0)**2)*rmusin(2,0)*pi/cellscale(2) +  & 
     &         rmu(2,0)*(-3*rmucos(2,0)**2 + 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*24*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*24*rmu(1,0)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                 


                     elseif(ic.eq.6) then 


                        if(i.eq.1) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)* &
     &           (rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) + &
     &           3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (6*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &          (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-6*cost4g*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=6*cost4g*fun0*rmu(1,0)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then
                   der4f(2)=cost4g*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(-4*rmu(1,0)*(rmu(1,0)**2 - &
     &          3*rmu(2,0)**2)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 



                     elseif(ic.eq.7) then


                        if(i.eq.1) then
                   der4f(2)=cost4g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*6*rmu(2,0)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)*(3*rmu(1,0)**2* &
     &            (rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &            rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (6*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + & 
     &            rmu(2,0)*(3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=-cost4g*(rmu(2,0)*(-3*rmu(1,0)**2 +  &
     &            rmu(2,0)**2)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*  &
     &            rmucos(2,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 

                     elseif(ic.eq.8) then


                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*4*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &            rmusin(1,0)**2) + 3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*cost5g*fun0*(3*rmu(1,0)**2*(rmucos(2,0)**2 - &
     &            rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                     else 

                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(2,0)*rmucos(2,0)*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &             rmusin(1,0)**2) + rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(8*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &       (2*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &       (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - &
     &             rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(-8*rmu(1,0)*rmucos(1,0)*rmucos(2,0)* &
     &        (2*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + rmu(2,0)*  &
     &       (3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+9
         indorb=indorbp 


      case(56) 
! g single Slater orbital                                             
! derivative of 55
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization 
          c=dd1**5.5d0*.020104801169736915d0
!         c=1.d0/dsqrt(7.d0)*(2.d0/pi)**(1.d0/2.d0)*dd1**(11.d0/2.d0)/3.d0/5.d0
!         endif 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=c*dexp(-dd1*r(k)) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*rp0+3.d0*rp0**2)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(7.d0*rmu(3,i)**2-3.d0*rp0)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2-rp0)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)  &
     &           *(7.d0*rmu(3,i)**2-rp0)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmucos(1,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(rmu(1,i)**2-3.0*rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmucos(2,i)*rmu(3,i)*rmucos(3,i)  &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)   &
     &   *(rmu(1,i)**2-rmu(2,i)**2)      
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F * phase

         do ic=1,9
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
         z(indorbp,k)=distp(k,1)*(5.5d0/dd1-r(k))*distp(k,1+ic)*cosphase(k)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                 
! fun2 = radial''
            fun0=distp(0,1)*(5.5d0/dd1-r(0)) 
            fun=distp(0,1)*(dd1-6.5d0/r(0)) 
            fun2=distp(0,1)*(7.5d0*dd1-dd1**2*r(0))  
                                                       
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,9 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                   der4f(2)=cost1g*12*rmu(1,0)*rmucos(1,0)*(-5*rmu(3,0)**2 + rp0) 
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(1,0)**2*rmucos(1,0)**2 -   &
     &        (rmucos(1,0)**2 - rmusin(1,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.2) then
                   der4f(2)=cost1g*12*rmu(2,0)*rmucos(2,0)*(-5*rmu(3,0)**2 + rp0)  
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*12*(2*rmu(2,0)**2*rmucos(2,0)**2 - &
     &        (rmucos(2,0)**2 - rmusin(2,0)**2)*(5*rmu(3,0)**2 - rp0))
                        elseif(i.eq.3) then
                   der4f(2)=cost1g*16*rmu(3,0)*rmucos(3,0)*(5*rmu(3,0)**2 - 3*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost1g*16*(rmu(3,0)**2*(9*rmucos(3,0)**2 - &
     &         5*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*(-6*rmu(1,0)**2*rmucos(1,0)**2 + &
     &           (rmucos(1,0)**2 - rmusin(1,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*(-2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (9*rmu(1,0)*(rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &         2*Pi/cellscale(1)*rmusin(1,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.2) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &            rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*6*rmu(1,0)*rmu(3,0)*   &
     &            rmucos(1,0)*rmucos(3,0)*(-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then 
                   der4f(2)=cost2g*rmu(1,0)*rmucos(1,0)* &
     &           (rmu(3,0)**2*(15*rmucos(3,0)**2 - 7*rmusin(3,0)**2) + &
     &           3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2g*fun0*4*rmu(1,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) +  & 
     &           3*Pi/cellscale(3)*rmusin(3,0)*rp0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                   der4f(2)=-cost2g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*  &
     &             rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*6*rmu(2,0)*rmu(3,0)*    &
     &             rmucos(2,0)*rmucos(3,0)*(-rmucos(1,0)**2 + rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost2g*rmu(3,0)*rmucos(3,0)*    &
     &         (-6*rmu(2,0)**2*rmucos(2,0)**2 + (rmucos(2,0)**2 - & 
     &          rmusin(2,0)**2)*(7*rmu(3,0)**2 - 3*rp0))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*  &
     &         (9*rmu(2,0)*(rmucos(2,0)**2 - rmusin(2,0)**2) +   &
     &      2*Pi/cellscale(2)*rmusin(2,0)*(7*rmu(3,0)**2 - 3*rp0)))
                        elseif(i.eq.3) then
                   der4f(2)=cost2g*rmu(2,0)*rmucos(2,0)*(rmu(3,0)**2*(15*rmucos(3,0)**2 -   &
     &      7*rmusin(3,0)**2) + 3*(-rmucos(3,0)**2 + rmusin(3,0)**2)*rp0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost2g*(4*rmu(2,0)*rmucos(2,0)*rmucos(3,0)* &
     &     (rmu(3,0)*(6*rmucos(3,0)**2 - 13*rmusin(3,0)**2) + 3*Pi/cellscale(3)*rmusin(3,0)*rp0))
                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                   der4f(2)=-cost3g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(3,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*fun0*cost3g*(rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) +  &
     &             3*rmu(3,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost3g*4*rmu(2,0)*(rmu(2,0)**2 - 3*rmu(3,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=4*fun0*cost3g*(rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) + &
     &             3*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*12*(rmu(1,0)**2 - rmu(2,0)**2)* &
     &             (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                        


                     elseif(ic.eq.5) then 


                        if(i.eq.1) then
                   der4f(2)=-2*cost3g*rmu(2,0)*rmucos(2,0)*((rmu(2,0)**2 - 6*rmu(3,0)**2)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2) +   &
     &             rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*(2*(rmu(2,0)**2 - & 
     &              6*rmu(3,0)**2)*rmusin(1,0)*pi/cellscale(1) +  & 
     &         rmu(1,0)*(-3*rmucos(1,0)**2 + 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-2*cost3g*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &         rmu(2,0)**2*(3*rmucos(2,0)**2 - rmusin(2,0)**2) +  &
     &         6*rmu(3,0)**2*(-rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost3g*fun0*(4*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*  &
     &         (2*(rmu(1,0)**2 - 6*rmu(3,0)**2)*rmusin(2,0)*pi/cellscale(2) +  & 
     &         rmu(2,0)*(-3*rmucos(2,0)**2 + 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=cost3g*24*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost3g*24*rmu(1,0)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                        endif                 


                     elseif(ic.eq.6) then 


                        if(i.eq.1) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)* &
     &           (rmu(1,0)**2*(3*rmucos(1,0)**2 - rmusin(1,0)**2) + &
     &           3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(2*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &          (6*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &          (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=-6*cost4g*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=6*cost4g*fun0*rmu(1,0)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)* &
     &           (-rmucos(2,0)**2 + rmusin(2,0)**2)
                        elseif(i.eq.3) then
                   der4f(2)=cost4g*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)* &
     &          (rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(-4*rmu(1,0)*(rmu(1,0)**2 - &
     &          3*rmu(2,0)**2)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 



                     elseif(ic.eq.7) then


                        if(i.eq.1) then
                   der4f(2)=cost4g*6*rmu(1,0)*rmu(2,0)*rmu(3,0)*rmucos(1,0)*rmucos(2,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*6*rmu(2,0)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (rmucos(1,0)**2 - rmusin(1,0)**2)
                        elseif(i.eq.2) then
                   der4f(2)=cost4g*rmu(3,0)*rmucos(3,0)*(3*rmu(1,0)**2* &
     &            (rmucos(2,0)**2 - rmusin(2,0)**2) +  & 
     &            rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost4g*fun0*(-2*rmu(3,0)*rmucos(2,0)*rmucos(3,0)* &
     &            (6*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + & 
     &            rmu(2,0)*(3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        elseif(i.eq.3) then
                   der4f(2)=-cost4g*(rmu(2,0)*(-3*rmu(1,0)**2 +  &
     &            rmu(2,0)**2)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=fun0*cost4g*(4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*  &
     &            rmucos(2,0)*rmucos(3,0)*rmusin(3,0)*pi/cellscale(3))
                        endif                 

                     elseif(ic.eq.8) then


                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(1,0)*(rmu(1,0)**2 - 3*rmu(2,0)**2)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*4*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &            rmusin(1,0)**2) + 3*rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(2,0)*(-3*rmu(1,0)**2 + rmu(2,0)**2)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-4*cost5g*fun0*(3*rmu(1,0)**2*(rmucos(2,0)**2 - &
     &            rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                     else 

                       if(i.eq.1) then
                   der4f(2)=cost5g*4*rmu(2,0)*rmucos(2,0)*(rmu(1,0)**2*(3*rmucos(1,0)**2 - &
     &             rmusin(1,0)**2) + rmu(2,0)**2*(-rmucos(1,0)**2 + rmusin(1,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(8*rmu(2,0)*rmucos(1,0)*rmucos(2,0)* &
     &       (2*rmu(2,0)**2*rmusin(1,0)*pi/cellscale(1) + rmu(1,0)* &
     &       (3*rmucos(1,0)**2 - 5*rmusin(1,0)**2)))
                        elseif(i.eq.2) then
                   der4f(2)=cost5g*4*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2*(rmucos(2,0)**2 - &
     &             rmusin(2,0)**2) + rmu(2,0)**2*(-3*rmucos(2,0)**2 + rmusin(2,0)**2))
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost5g*fun0*(-8*rmu(1,0)*rmucos(1,0)*rmucos(2,0)* &
     &        (2*rmu(1,0)**2*rmusin(2,0)*pi/cellscale(2) + rmu(2,0)*  &
     &       (3*rmucos(2,0)**2 - 5*rmusin(2,0)**2)))
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


! rest of laplacian (with the phase)
! d^2/dx^2_k \Psi(X_i) = d^2/dx^2_k F(X_i) \phi(x_i) + 2 d/dx_k F(X_i)  d/dx_k \phi(x_i) 
!                        + F(X_i) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) already computed and stored in hess4f(k)
! d/dx_k F(X_i) already computed and stored in z(indorbp,indt+k)
! F(X_i) = distp(0,1)*distp(0,1+ic)
                      hess4f(1)=hess4f(1)*cosphase(0) &
                           -2.d0*z(indorbp,indt+i)*sinphase(i,0)*rphase(i) &
                           -distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*cosphase(0)*rphase(i)**2

! gradient (with the phase)
! compute d/dx_k \Psi(X_i)
!  (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
                      z(indorbp,indt+i)=z(indorbp,indt+i)*cosphase(0) &
                           -distp(0,1)*(11.d0/4.d0/dd1-r(0)**2)*distp(0,1+ic)*sinphase(i,0)*rphase(i)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+9
         indorb=indorbp 





! ******************* END GAUSSIAN BASIS ************************       
! ** ** ** ** ** ** **  JASTROW ORBITALS ** ** ** ** ** ** ** ** *      
       case(100) 
!     2s single gaussian                                                
!     exp(-dd2*r^2)                                                     
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1) 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
               fun=-dd2*distp(0,1)*2.d0 
               fun3=dd2**2*4.d0*distp(0,1) 
                                                                        
                  do i=1,3 
                     z(indorbp,indt+i)=fun*rmu(i,0)                 &
     & *rmucos(i,0)                                                 
                                                                        
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
!       hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 
!    1+fun*(rmucos(:,0)**2-PI*rmu(:,0)*rmusin(:,0))         
                                                                        
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
                                                                        
                                                                        
       case(145) 
!     2s without cusp condition  !derivative 100                        
!     -(r^2*exp(-dd2*r^2))                                              
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=-distp(i,1)*r(i)**2 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
               fun0=dd2*r(0)**2 
               fun=-2.d0*distp(0,1)*(1.d0-fun0) 
                fun3=+dd2*4.d0*distp(0,1)*(1.d0-fun0)                   &
     &+4.d0*distp(0,1)*dd2                                              
                                                                        
                  do i=1,3 
                  z(indorbp,indt+i)=fun*rmu(i,0)                    &
     & *rmucos(i,0)                                                 
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
!       hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 
!    1+fun*(rmucos(:,0)**2-PI*rmu(:,0)*rmusin(:,0))         
                                                                        
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(rmucos(:,0)**2                                          &
     &-PI*rmu(:,0)*rmusin(:,0)/cellscale(1:3))                    
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
       case(10) 
!     2s with cusp condition                                            
!     ( r^2*exp(-dd2*r))  ! with no cusp condition                      
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*r(i)**2 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            fun=(2.d0-dd2*r(0))*distp(0,1) 
            fun3=-dd2*distp(0,1)/r(0)*(3.d0-dd2*r(0)) 
                                                                        
                  do i=1,3 
                  z(indorbp,indt+i)=fun*rmu(i,0)                    &
     & *rmucos(i,0)                                                 
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
                                                                        

        case(131) 
!     2s without cusp condition                                         
!     dd1*(r^2*exp(-dd2*r^2))             
          
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)**2)
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*r(i)**2 
               enddo 
            endif 
                                                                        
         if(typec.ne.1) then 
               fun =2.d0*distp(0,1)*(1.d0-dd2*r(0)**2)
               fun3=-dd2*distp(0,1)*4.d0*(2.d0-dd2*r(0)**2)
                                                                        
           do i=1,3 
             z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)                 
           enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 



         case(149) 
! derivative of 131 with respect z_1                                    
! - r^4 exp(-z_1 r^2)                                                                         
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
            dd1=dd(indpar+1) 
                                                                        
            do k=indtmin,indtm 
            distp(k,1)=dexp(-dd1*r(k)**2) 
            enddo 
                                                                        
                                                                        
            do i=i0,indtm 
             z(indorbp,i)=-r(i)**4*distp(i,1)
            enddo 
            if(typec.ne.1) then 
               fun =-2.d0*distp(0,1)*r(0)**2*(2.d0-dd1*r(0)**2)

               fun3=-distp(0,1)*4.d0*(2.d0-2.d0*dd1*r(0)**2  &
           &   -(2.d0-dd1*r(0)**2)*r(0)**2*dd1)
                                                                        
            do i=1,3 
               z(indorbp,indt+i)=fun*rmu(i,0)*rmucos(i,0)                                                 
            enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:))                                                                 

            endif 
                                                                        
            indorb=indorbp 
                                                                        
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 


        case(137) 
!     2s with cusp condition                                            
!     dd1*(exp(-dd2*r)*(1+dd2*r))                                       
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
!           if(iflagnorm.gt.2) then                                     
!           c=1.d0/dsqrt(1/4.d0/dd2**3+12*dd2/(2.d0*dd2)**4+            
!    &3*dd2**2/4/dd2**5)/dsqrt(4.0*pi)                                  
!           endif                                                       
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*(1.d0+dd2*r(i)) 
               enddo 
            endif 
                                                                        
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
               fun=-dd2**2*distp(0,1) 
               fun3=dd2**3*distp(0,1)/r(0) 
                                                                        
                  do i=1,3 
                  z(indorbp,indt+i)=fun*rmu(i,0)                    &
     & *rmucos(i,0)                                                 
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
                                                                        
       case(138) 
!     2s with cusp condition                                            
!     ( -dd2*r^2*exp(-dd2*r))  ! with no cusp condition der of 137      
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=-dd2*dexp(-dd2*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*r(i)**2 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
               fun=(2.d0-dd2*r(0))*distp(0,1) 

               fun3=-dd2*distp(0,1)/r(0)                            &
     &-dd2*(2.d0-dd2*r(0))*distp(0,1)/r(0)                      
                                                                        
                 do i=1,3 
                  z(indorbp,indt+i)=fun*rmu(i,0)                    &
     & *rmucos(i,0)                                                 
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
                                                                        
                  ! single gaussian p orbital periodic                  
       case(103) 
                                                                        
         dd1=dd(indpar+1) 
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
         z(indorbp,i)=distp(i,1)*rmu(ic,i)*rmucos(ic,i) 
              enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
            fun =-2.d0*dd1*distp(0,1) 
            fun3=4.d0*dd1**2*distp(0,1) 
                                                                        
               indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
            z(indorbp,indt+1)=rmu(1,0)*fun*rmucos(1,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+2)=rmu(2,0)*fun*rmucos(2,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+3)=rmu(3,0)*fun*rmucos(3,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
                                                                        
       z(indorbp,indt+ic)=z(indorbp,indt+ic)                            &
     &   +distp(0,1)*(2.d0*rmucos(ic,0)**2-1.d0)                    
                                                                        
        hess(:)=rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))                              
                                                                        
        hess(ic)=hess(ic)                                               &
     &+fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-distp(0,1)*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic)                                                 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
            enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                ! derivative of 103 unnormalized                        
      case(146) 
! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)                                      
                                                                        
          
         dd1=dd(indpar+1) 
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                   z(indorbp,i)=-rmu(ic,i)*distp(i,1)*              &
     &             r(i)**2*rmucos(ic,i)                         
               enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-2.d0) 
            fun0=-distp(0,1)*r(0)**2 
            fun3=dd1*distp(0,1)*(-4.d0*r(0)**2*dd1+8.d0) 
                                                                        
           indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
         z(indorbp,indt+1)=rmu(ic,0)*rmu(1,0)*                  &
     &                       fun*rmucos(1,0)*rmucos(ic,0)       
         z(indorbp,indt+2)=rmu(ic,0)*rmu(2,0)*                  &
     &                       fun*rmucos(2,0)*rmucos(ic,0)       
         z(indorbp,indt+3)=rmu(ic,0)*rmu(3,0)*                  &
     &                       fun*rmucos(3,0)*rmucos(ic,0)       
                                                                        
       z(indorbp,indt+ic)=z(indorbp,indt+ic)                            &
     &   +fun0*(2.d0*rmucos(ic,0)**2-1.d0)                          
                                                                        
        hess(:)=rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(rmucos(:,0)**2-rmusin(:,0)**2))                     
                                                                        
        hess(ic)=hess(ic)                                               &
     &+2.d0*(2.d0*rmucos(ic,0)**2-1.d0)                             &
     &*rmu(ic,0)*rmucos(ic,0)*fun                               &
     &-4.d0*fun0*rmucos(ic,0)*rmusin(ic,0)*PI/cellscale(ic)     
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
               enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                  !  2p single exponential  r e^{-z r^2}                
       case(150) 
                                                                        
         dd1=dd(indpar+1) 
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
         z(indorbp,i)=distp(i,1)*rmu(ic,i)*rmucos(ic,i)*r(i) 
              enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            fun0=distp(0,1)*r(0) 
            fun =(1.d0/r(0)-2.d0*dd1*r(0))*distp(0,1) 
            rp1=2.d0*dd1*r(0)**2 
            fun3=(rp1**2-1.d0-2.d0*rp1)*distp(0,1)/r(0)**3 
                                                                        
                                                                        
               indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
            z(indorbp,indt+1)=rmu(1,0)*fun*rmucos(1,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+2)=rmu(2,0)*fun*rmucos(2,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+3)=rmu(3,0)*fun*rmucos(3,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
                                                                        
       z(indorbp,indt+ic)=z(indorbp,indt+ic)                            &
     &   +fun0*(2.d0*rmucos(ic,0)**2-1.d0)                          
                                                                        
                                                                        
       hess(:)=rmu(ic,0)*rmucos(ic,0)*                          &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(rmucos(:,0)**2-rmusin(:,0)**2))                     
                                                                        
        hess(ic)=hess(ic)                                               &
     &+2.d0*(2.d0*rmucos(ic,0)**2-1.d0)                             &
     &*rmu(ic,0)*rmucos(ic,0)*fun                               &
     &-4.d0*fun0*rmucos(ic,0)*rmusin(ic,0)*PI/cellscale(ic)     
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
            enddo 
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                ! derivative of 150 unnormalized                        
      case(151) 
! R(r)=exp(-z*r**2)*(-r**3)                                             
                                                                        
          
         dd1=dd(indpar+1) 
                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                   z(indorbp,i)=-rmu(ic,i)*distp(i,1)*              &
     &             r(i)**3*rmucos(ic,i)                         
               enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
                                                                        
            rp1=2.d0*dd1*r(0)**2 
            fun=distp(0,1)*(rp1-3.d0)*r(0) 
            fun0=-distp(0,1)*r(0)**3 
            fun3=distp(0,1)/r(0)*(-3.d0+6.d0*rp1-rp1**2) 
                                                                        
           indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
         z(indorbp,indt+1)=rmu(ic,0)*rmu(1,0)*                  &
     &                       fun*rmucos(1,0)*rmucos(ic,0)       
         z(indorbp,indt+2)=rmu(ic,0)*rmu(2,0)*                  &
     &                       fun*rmucos(2,0)*rmucos(ic,0)       
         z(indorbp,indt+3)=rmu(ic,0)*rmu(3,0)*                  &
     &                       fun*rmucos(3,0)*rmucos(ic,0)       
                                                                        
       z(indorbp,indt+ic)=z(indorbp,indt+ic)                            &
     &   +fun0*(2.d0*rmucos(ic,0)**2-1.d0)                          
                                                                        
        hess(:)=rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(rmucos(:,0)**2-rmusin(:,0)**2))                     
                                                                        
        hess(ic)=hess(ic)                                               &
     &+2.d0*(2.d0*rmucos(ic,0)**2-1.d0)                             &
     &*rmu(ic,0)*rmucos(ic,0)*fun                               &
     &-4.d0*fun0*rmucos(ic,0)*rmusin(ic,0)*PI/cellscale(ic)     
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
               enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
       case(152) 
!     2s single gaussian                                                
!     r^3 exp(-dd2*r^2)                                                 
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*r(i)**3 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
               rp1=2.d0*dd2*r(0)**2 
               fun=r(0)*(3.d0-rp1)*distp(0,1) 
               fun3=(3.d0-6.d0*rp1+rp1**2)*distp(0,1)/r(0) 
                                                                        
                  do i=1,3 
                     z(indorbp,indt+i)=fun*rmu(i,0)                 &
     & *rmucos(i,0)                                                 
                                                                        
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
                                                                        
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                
         endif  !endif for indt                                      
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
                                                                        
                                                                        
       case(153) 
!     2s without cusp condition  !derivative 152                        
!     -(r^5*exp(-dd2*r^2))                                              
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=dexp(-dd2*r(k)*r(k)) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=-distp(i,1)*r(i)**5 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
               rp1=2.d0*dd2*r(0)**2 
               fun=r(0)**3*distp(0,1)*(rp1-5.d0) 
               fun3=-r(0)*(15.d0-10.d0*rp1+rp1**2)*distp(0,1) 
                                                                        
                                                                        
                 do i=1,3 
                  z(indorbp,indt+i)=fun*rmu(i,0)                    &
     & *rmucos(i,0)                                                 
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
                                                                        
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 
                                                                        
       case(39) 
!     2s with cusp condition                                            
!     ( r^3*exp(-dd2*r))  !  der of   10                                
          
                                                                        
         dd2=dd(indpar+1) 
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
               do k=indtmin,indtm 
               distp(k,1)=-dexp(-dd2*r(k))*r(k) 
               enddo 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=distp(i,1)*r(i)**2 
               enddo 
            endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
               fun=(3.d0-dd2*r(0))*distp(0,1) 
               fun3=-dd2*distp(0,1)/r(0)                            &
     &+(3.d0-dd2*r(0))*(-dd2*distp(0,1)+distp(0,1)/r(0))        &
     &/r(0)                                                         
                                                                        
                 do i=1,3 
                  z(indorbp,indt+i)=fun*rmu(i,0)                    &
     & *rmucos(i,0)                                                 
                  enddo 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
                                                                        
                                                                        
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshellp 
         indorb=indorbp 


       case(147)
! d single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) 
! R depends on X_i, A depends on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i)
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 

         dd(indparp)=abs(dd(indparp))                                                                         
         dd1=dd(indparp) 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=(3.d0*rmu(3,i)**2-rp0)*cost1d  ! lz=0
                                                                        
      distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d  ! lz=+/-2      
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d   ! lz=+/-2
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d  ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d  ! lz=+/-1    
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular
! F = radial * angular 
! compute \Psi = F

         do ic=1,5
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,5
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then

                   der4f(2)=-2.d0*cost1d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(1,0)**2 + rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost1d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                   der4f(2)=4.d0*cost1d*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*4.d0*(rmucos(3,0)**2-rmusin(3,0)**2)

                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then

                   der4f(2)=2.d0*cost2d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(rmucos(1,0)**2 - rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost2d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*rmusin(1,0)

                        elseif(i.eq.2) then


                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*rmusin(2,0)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        elseif(i.eq.2) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*rmusin(2,0)

                        else


                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(2,0)*rmucos(2,0)*rmucos(3,0)*rmusin(3,0)

                        endif                        


                     else


                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)*rmusin(1,0)

                        elseif(i.eq.2) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        else

                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(1,0)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)

                        endif                 

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)
                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+5
         indorb=indorbp 


       case(148)
! d single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! normalized 
! orbital with the minimal power of rmucos
! orbital derivative of 147

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i)
! R depends on X_i, A depends on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i)))
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 

         dd(indparp)=abs(dd(indparp))                                                                         
         dd1=dd(indparp) 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos with the minimal rmucos strategy
          do i=indtmin,indtm 

      rp0=r(i)**2

      distp(i,2)=(3.d0*rmu(3,i)**2-rp0)*cost1d  ! lz=0
                                                                        
      distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d  ! lz=+/-2      
                                                                        
      distp(i,4)=rmu(1,i)*rmu(2,i)                              &
     &*rmucos(1,i)*rmucos(2,i)*cost3d   ! lz=+/-2
                                                                        
      distp(i,5)=rmu(2,i)*rmu(3,i)                              &
     &*rmucos(2,i)*rmucos(3,i)*cost3d  ! lz=+/-1                
                                                                        
      distp(i,6)=rmu(1,i)*rmu(3,i)                              &
     &*rmucos(1,i)*rmucos(3,i)*cost3d  ! lz=+/-1    
                                                      ! lz=+/-4
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular
! F = radial * angular 
! compute \Psi = F

         do ic=1,5
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=-r(k)**2*distp(k,1)*distp(k,1+ic)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=-r(0)**2*distp(0,1)
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-2.d0)  
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +10.d0*dd1*r(0)**2-2.d0)  

            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,5
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part ! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then

                   der4f(2)=-2.d0*cost1d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(1,0)**2 + rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost1d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                   der4f(2)=4.d0*cost1d*rmu(3,0)*rmucos(3,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost1d*fun0*4.d0*(rmucos(3,0)**2-rmusin(3,0)**2)

                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then

                   der4f(2)=2.d0*cost2d*rmu(1,0)*rmucos(1,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(rmucos(1,0)**2 - rmusin(1,0)**2)

                        elseif(i.eq.2) then

                   der4f(2)=-2.d0*cost2d*rmu(2,0)*rmucos(2,0)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=cost2d*fun0*2.d0*(-rmucos(2,0)**2 + rmusin(2,0)**2)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(2,0)*rmucos(1,0)*rmucos(2,0)*rmusin(1,0)

                        elseif(i.eq.2) then


                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(1,0)*rmucos(1,0)*rmucos(2,0)*rmusin(2,0)

                        else

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        endif                        

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        elseif(i.eq.2) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(2,0)**2 - rmusin(2,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(2)*rmu(3,0)*rmucos(2,0)*rmucos(3,0)*rmusin(2,0)

                        else


                   der4f(2)=cost3d*rmu(2,0)*rmucos(2,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(2,0)*rmucos(2,0)*rmucos(3,0)*rmusin(3,0)

                        endif                        


                     else


                        if(i.eq.1) then

                   der4f(2)=cost3d*rmu(3,0)*rmucos(3,0)*(rmucos(1,0)**2 - rmusin(1,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(1)*rmu(3,0)*rmucos(1,0)*rmucos(3,0)*rmusin(1,0)

                        elseif(i.eq.2) then

                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0

                        else

                   der4f(2)=cost3d*rmu(1,0)*rmucos(1,0)*(rmucos(3,0)**2 - rmusin(3,0)**2)
                   hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                   der4f(2)=fun0*der4f(2)
                   hess4f(2)=-cost3d*fun0*4.d0*Pi/cellscale(3)*rmu(1,0)*rmucos(1,0)*rmucos(3,0)*rmusin(3,0)

                        endif                 

                        !endif for ic                           
                     endif


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+5
         indorb=indorbp 


       case(154)
! f single gaussian orbital                                             
! radial R(r)= exp(-alpha r^2)                                                 
! unnormalized and NO PHASE (Jastrow orbital)
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

     rp0=r(i)**2

     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)
                                                          ! lz=0        
     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
    &  *(rmu(1,i)**2-rmu(2,i)**2)

     distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)          &
    &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             

     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)               &
    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 
     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)               &
    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
                                                                        
            fun0=distp(0,1) 
            fun=-2.d0*dd1*distp(0,1) 
            fun2=fun*(1.d0-2.d0*dd1*r(0)**2)
            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                           der4f(2)=-6.d0*cost1f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(1,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost1f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(2,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.3) then
                           der4f(2)=4.d0*cost1f*rmu(3,0)**2*rmucos(3,0)**2  &
                                +cost1f*(2.d0*rmucos(3,0)**2-1.d0)*(5.d0*rmu(3,0)**2-3.d0*rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=cost1f*fun0*(12.d0*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(3,0)**2-1.d0) &
                                -4.d0*Pi**2/cellscale(3)**2*rmu(3,0)*rmucos(3,0)*(5.d0*rmu(3,0)**2-3.d0*rp0))
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)**2*rmucos(1,0)**2 &
                                    +cost2f*(2.d0*rmucos(1,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-6.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
          &-4.d0*PI**2/cellscale(1)**2*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
        der4f(2)=-2.d0*cost2f*rmu(2,0)**2*rmucos(2,0)**2 &
                                   & +cost2f*(2.d0*rmucos(2,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
        hess4f(2)=-6.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
                 &-4.d0*PI**2/cellscale(2)**2*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(1,0)**2-1.d0)            
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=cost3f*(2.d0*rmucos(3,0)**2-1.d0)*(rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-4.d0*fun0*cost3f*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0)*(rmu(1,0)**2-rmu(2,0)**2)
                        endif

                     elseif(ic.eq.5) then 

                        if(i.eq.1) then
   der4f(2)=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        elseif(i.eq.2) then
   der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        else
   der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
                &*rmusin(i,0)*rmucos(i,0)
                        endif

                     elseif(ic.eq.6) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost4f*rmu(1,0)**2*rmucos(1,0)**2  &
                                   +cost4f*(2.d0*rmucos(1,0)**2-1.d0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
                &-4.d0*fun0*cost4f*PI**2/cellscale(1)**2*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                     else 

                        if(i.eq.1) then
                           der4f(2)=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
     der4f(2)=-2.d0*cost4f*rmu(2,0)**2*rmucos(2,0)**2  &
          &+cost4f*(2.d0*rmucos(2,0)**2-1.d0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
          &-4.d0*fun0*cost4f*PI**2/cellscale(2)**2*rmu(2,0)*rmucos(2,0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                        !endif for ic                           
                     endif
                       !enddo for i                                 


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)

!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 


       case(155)
! f single gaussian orbital derivative of 154                                             
! radial R(r)= -r^2*exp(-z r^2)                                                  
! unnormalized and NO PHASE (Jastrow orbital)
! orbital with the minimal power of rmucos

! metric set by boundary conditions
! general representation  X_i(x_i) where x_i are the plain cartesian coordinates
! X_i are the stretched coordinates for the radial part R

! general form of the orbital \Psi(x_i) = R(X_i(x_i)) A(x_i) \phi(x_i)
! R depends on X_i, A and the phase \phi depend on x_i

! gradient
! d/dx_k \Psi(x_i) = (d/dx_k F(X_i(x_i))) \phi(x_i) + F(X_i(x_i)) d/dx_k \phi(x_i) 
! where  d/dx_k F(X_i(x_i)) = (d/dX_k R) A (d/dx_k X_i) (x_i) + R  (d/dY_k R) A (d/dx_k Y_i) (x_i)
! (d/dX_k R) (d/dY_k A)  are calculated as in makefun but evaluated at stretched coordinates X_k Y_k

!laplacian 
! d^2/dx^2_k \Psi(x_i) = (d^2/dx^2_k F(X_i(x_i))) \phi(x_i) 
!                      + 2.d0 (d/dx_k F(X_i(x_i)))  d/dx_k \phi(x_i)
!                      + F(X_i(x_i)) d^2/dx^2_k \phi(x_i)
! d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

! X_i = rmu
! d/dx_k X_i  =  rmucos
! d^2/dx^2_k X_i = - Pi/L * rmusin

          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
! radial part                                                                        
         do k=indtmin,indtm 
         distp(k,1)=dexp(-dd1*r(k)**2) 
         enddo 
                                                                        
! angular part                                                                        
! rmu replaced by rmu*rmucos
          do i=indtmin,indtm 

     rp0=r(i)**2

     distp(i,2)=cost1f*rmu(3,i)*rmucos(3,i)               &
    &  *(5.d0*rmu(3,i)**2-3.d0*rp0)
                                                          ! lz=0        
     distp(i,3)=cost2f*rmu(1,i)*rmucos(1,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,4)=cost2f*rmu(2,i)*rmucos(2,i)               &
    &  *(5.d0*rmu(3,i)**2-rp0)

     distp(i,5)=cost3f*rmu(3,i)*rmucos(3,i)               &
    &  *(rmu(1,i)**2-rmu(2,i)**2)

     distp(i,6)=cost3f*2.d0*rmu(3,i)*rmucos(3,i)          &
    &  *rmu(1,i)*rmucos(1,i)*rmu(2,i)*rmucos(2,i)                             

     distp(i,7)=cost4f*rmu(1,i)*rmucos(1,i)               &
    &  *(rmu(1,i)**2-3.d0*rmu(2,i)**2)
 
     distp(i,8)=cost4f*rmu(2,i)*rmucos(2,i)               &
    &  *(3.d0*rmu(1,i)**2-rmu(2,i)**2)
  
          enddo 
                                                                        
! definition of the orbital \Psi
! radial * angular * phase
! F = radial * angular 
! compute \Psi = F

         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=-r(k)**2*distp(k,1)*distp(k,1+ic)
              enddo 
            endif 
         enddo 
                                                                        
! calculate (d/dx_k F) and store them in z(indorbp,indt+k)  
! calculate (d^2/dx^2_k F) and store them in hess4f(k,indorbp)
! d^2/dx^2_k R(x) = (fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2)) 

                                                                     
         if(typec.ne.1) then 

! fun0 = radial
! fun = radial'/r
! fun2 = radial''
            fun0=-r(0)**2*distp(0,1)
            fun=distp(0,1)*(2.d0*dd1*r(0)**2-2.d0)  
            fun2=distp(0,1)*(-4.d0*dd1**2*r(0)**4                   &
     &      +10.d0*dd1*r(0)**2-2.d0)  

            rp0=r(0)**2

! compute:  fun2 * x^2_k/r^2 + fun * (1 - x^2_k/r^2) and store it in radhess4f
            do i=1,3
               radhess4f(i)=fun2*rmu(i,0)**2/r(0)**2+fun*(1.d0-rmu(i,0)**2/r(0)**2)
            enddo
                                                                        
            indorbp=indorb 
            do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                  indorbp=indorbp+1 

                  z(indorbp,indt+4)=0.d0  

                  do i=1,3 

! derivatives of the radial part
! der4f(1) = (d/dx_k R) A 
                     der4f(1)=distp(0,1+ic)*rmu(i,0)*fun                                              

! hess4f(1) = (d^2/dx^2_k R) A
                     hess4f(1)=distp(0,1+ic)*radhess4f(i)

! contribution of the angular part
! der4f(2) = R (d/dx_k A)
! hess4f(2) = R (d^2/dx^2_k A)
! hess4f(3) = (d/dx_k R) (d/dx_k A)
                     if(ic.eq.1) then

                        if(i.eq.1) then
                           der4f(2)=-6.d0*cost1f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(1,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost1f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost1f*(2.d0*rmucos(2,0)**2-1.d0)*rmu(3,0)*rmucos(3,0)
                        elseif(i.eq.3) then
   der4f(2)=4.d0*cost1f*rmu(3,0)**2*rmucos(3,0)**2  &
          &+cost1f*(2.d0*rmucos(3,0)**2-1.d0)*(5.d0*rmu(3,0)**2-3.d0*rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
   hess4f(2)=cost1f*fun0*(12.d0*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(3,0)**2-1.d0) &
         &-4.d0*Pi**2/cellscale(3)**2*rmu(3,0)*rmucos(3,0)*(5.d0*rmu(3,0)**2-3.d0*rp0))
                        endif                        

                     elseif(ic.eq.2) then
                        
                        if(i.eq.1) then
   der4f(2)=-2.d0*cost2f*rmu(1,0)**2*rmucos(1,0)**2 &
               &+cost2f*(2.d0*rmucos(1,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
              &-4.d0*PI**2/cellscale(1)**2*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.3) then

                        if(i.eq.1) then
                           der4f(2)=-2.d0*cost2f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
    der4f(2)=-2.d0*cost2f*rmu(2,0)**2*rmucos(2,0)**2 &
            &+cost2f*(2.d0*rmucos(2,0)**2-1.d0)*(5.d0*rmu(3,0)**2-rp0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-6.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
            &-4.d0*PI**2/cellscale(2)**2*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(5.d0*rmu(3,0)**2-rp0)
                        elseif(i.eq.3) then 
                           der4f(2)=8.d0*cost2f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=8.d0*fun0*cost2f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(3,0)**2-1.d0)
                        endif

                     elseif(ic.eq.4) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(1,0)**2-1.d0)            
                        elseif(i.eq.2) then
                           der4f(2)=-2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-2.d0*fun0*cost3f*rmu(3,0)*rmucos(3,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=cost3f*(2.d0*rmucos(3,0)**2-1.d0)*(rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-4.d0*fun0*cost3f*PI/cellscale(3)*rmusin(3,0)*rmucos(3,0) &
             &*(rmu(1,0)**2-rmu(2,0)**2)
                        endif

                     elseif(ic.eq.5) then 

                        if(i.eq.1) then
     der4f(2)=2.d0*cost3f*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0)&
           &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(2,0)*rmucos(2,0)*rmu(3,0)*rmucos(3,0) &
           &*rmusin(i,0)*rmucos(i,0)
                        elseif(i.eq.2) then
     der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
           &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(3,0)*rmucos(3,0) &
         &*rmusin(i,0)*rmucos(i,0)
                        else
      der4f(2)=2.d0*cost3f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
              &*(2.d0*rmucos(i,0)**2-1.d0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
    hess4f(2)=-8.d0*cost3f*fun0*PI/cellscale(i)*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0) &
         &*rmusin(i,0)*rmucos(i,0)
                        endif

                     elseif(ic.eq.6) then 

                        if(i.eq.1) then
                           der4f(2)=2.d0*cost4f*rmu(1,0)**2*rmucos(1,0)**2  &
                                   +cost4f*(2.d0*rmucos(1,0)**2-1.d0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(1,0)**2-1.d0) &
        &-4.d0*fun0*cost4f*PI**2/cellscale(1)**2*rmu(1,0)*rmucos(1,0)*(rmu(1,0)**2-3.d0*rmu(2,0)**2)
                        elseif(i.eq.2) then
                           der4f(2)=-6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=-6.d0*fun0*cost4f*rmu(1,0)*rmucos(1,0)*(2.d0*rmucos(2,0)**2-1.d0)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif

                     else 

                        if(i.eq.1) then
                           der4f(2)=6.d0*cost4f*rmu(1,0)*rmucos(1,0)*rmu(2,0)*rmucos(2,0)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
                           hess4f(2)=6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(1,0)**2-1.d0)
                        elseif(i.eq.2) then
    der4f(2)=-2.d0*cost4f*rmu(2,0)**2*rmucos(2,0)**2  &
      &+cost4f*(2.d0*rmucos(2,0)**2-1.d0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                           hess4f(3)=2.d0*fun*rmu(i,0)*der4f(2)
                           der4f(2)=fun0*der4f(2)
     hess4f(2)=-6.d0*fun0*cost4f*rmu(2,0)*rmucos(2,0)*(2.d0*rmucos(2,0)**2-1.d0) &
       &-4.d0*fun0*cost4f*PI**2/cellscale(2)**2*rmu(2,0)*rmucos(2,0)*(3.d0*rmu(1,0)**2-rmu(2,0)**2)
                        else
                           der4f(2)=0.d0
                           hess4f(2)=0.d0
                           hess4f(3)=0.d0
                        endif


                        !endif for ic                           
                     endif
                       !enddo for i                                 


! gradient (without phase)
! compute (d/dX_k F) (X_i(x_i)) (d/dx_k X_i) (x_i)
!                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0) &
!     &                                 +der4f(2)*(2.d0*rmucos(i,0)**2-1.d0)
                      z(indorbp,indt+i)=der4f(1)*rmucos(i,0)+der4f(2)

! part of the laplacian (without phase)
! compute d^2/dx^2_k F(X_i) =  (d^2/dX^2_k R) A ((d/dx_k X_i) (x_i))^2 
!                      + (d/dX_k R) (X_i(x_i)) A (d^2/dx^2_k X_i) (x_i)
!                      + R (d^2/dx^2_k A)
!                      + 2 (d/dX_k R)  (d/dx_k A) (d/dx_k X_i) (x_i)

!                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
!     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
!     &                            +hess4f(2)*(2.d0*rmucos(i,0)**2-1.d0)**2  &
!     &                            -der4f(2)*4.d0*PI/cellscale(i)*rmusin(i,0)*rmucos(i,0) &
!     &                            +hess4f(3)*rmucos(i,0)*(2.d0*rmucos(i,0)**2-1.d0)

                      hess4f(1)=hess4f(1)*rmucos(i,0)**2  &
     &                            -der4f(1)*PI/cellscale(i)*rmusin(i,0) &
     &                            +hess4f(2)+hess4f(3)*rmucos(i,0)


!                     z_xyz(indorbp,i)=hess4f(1) 
                      z(indorbp,indt+4)=z(indorbp,indt+4)+hess4f(1) 

                   enddo
                    

                endif
                       ! enddo fot ic                                   
             enddo
         endif 

         indpar=indpar+1 
         indshell=indshell+7 
         indorb=indorbp 


                                                                        
       case(199) 
!     LA COSTANTE                                                       
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=0.d0 
               enddo 
            endif 
                                                                        
         if(typec.ne.1)  z(indorbp,indt+1:indt+4)=0.d0
                                                                        
         indshell=indshellp 
         indorb=indorbp 
                                                                        
                                                                        
        case(1000:1099) 
!     s gaussian  r**(2*npower)*exp(-alpha*r**2)                        
         npower=iopt-1000 
          
         indorbp=indorb+1 
         indshellp=indshell+1 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd2=dd(indpar+1) 
         do k=indtmin,indtm 
            distp(k,1)=r(k)**(2*npower)*dexp(-dd2*r(k)**2) 
         enddo 
                                                                        
         if(iocc(indshellp).eq.1) then 
            do i=i0,indtm 
               z(indorbp,i)=distp(i,1) 
            enddo 
         endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            rp1=r(0)**2 
            fun=-dd2*distp(0,1)*2.d0+npower*distp(0,1)*2.d0/rp1 
            fun3=-dd2*fun*2.d0+npower*fun*2.d0/rp1                      &
     &-npower*distp(0,1)*4.d0/rp1/rp1                                   
                                                                        
            if(iocc(indshellp).eq.1) then 
              z(indorbp,indt+1)=fun*rmu(1,0)                        &
     & *rmucos(1,0)                                                 
                  z(indorbp,indt+2)=fun*rmu(2,0)                    &
     & *rmucos(2,0)                                                 
                  z(indorbp,indt+3)=fun*rmu(3,0)                    &
     & *rmucos(3,0)                                                 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+1 
         indorb=indorbp 
                                                                        
                                                                        
        case(2000:2099) 
!     s gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative of 1000
                                                                        
         npower=iopt+1-2000 
          
         indorbp=indorb+1 
         indshellp=indshell+1 
                                                                        
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd2=dd(indpar+1) 
         do k=indtmin,indtm 
            distp(k,1)=-r(k)**(2*npower)*dexp(-dd2*r(k)**2) 
         enddo 
                                                                        
         if(iocc(indshellp).eq.1) then 
            do i=i0,indtm 
               z(indorbp,i)=distp(i,1) 
            enddo 
         endif 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            rp1=r(0)**2 
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1 
                                                                        
            fun3=fun*(npower-dd2*rp1)*2.d0/rp1                          &
     &-2.d0*dd2*distp(0,1)*2.d0/rp1                                     &
     &-(npower-dd2*rp1)*distp(0,1)*4.d0/rp1**2                          
                                                                        
        if(iocc(indshellp).eq.1) then 
              z(indorbp,indt+1)=fun*rmu(1,0)                        &
     & *rmucos(1,0)                                                 
                  z(indorbp,indt+2)=fun*rmu(2,0)                    &
     & *rmucos(2,0)                                                 
                  z(indorbp,indt+3)=fun*rmu(3,0)                    &
     & *rmucos(3,0)                                                 
                                                                        
!********** DIAGONAL PART HESSIAN ******************                    
        hess(:)=fun3*rmu(:,0)**2*rmucos(:,0)**2                 &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2)                               
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
            endif 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+1 
         indorb=indorbp 
                                                                        
                                                                        
         case(1100:1199) 
!     p gaussian  r**(2*npower)*exp(-alpha*r**2)                        
                                                                        
         npower=iopt-1100 
          
         indorbp=indorb 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd2=dd(indpar+1) 
         do k=indtmin,indtm 
            distp(k,1)=r(k)**(2*npower)*dexp(-dd2*r(k)**2) 
         enddo 
                                                                        
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                  z(indorbp,i)=rmu(ic,i)*distp(i,1)*rmucos(ic,i) 
               enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            rp1=r(0)**2 
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1 
                                                                        
            fun3=fun*(npower-dd2*rp1)*2.d0/rp1                          &
     &-4.d0*(npower-dd2*rp1)*distp(0,1)/rp1/rp1                         &
     &-dd2*2.d0*distp(0,1)*2.d0/rp1                                     
                                                                        
            indorbp=indorb 
                                                                        
                                                                        
                do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
            z(indorbp,indt+1)=rmu(1,0)*fun*rmucos(1,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+2)=rmu(2,0)*fun*rmucos(2,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+3)=rmu(3,0)*fun*rmucos(3,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
                                                                        
       z(indorbp,indt+ic)=z(indorbp,indt+ic)                            &
     &   +distp(0,1)*(2.d0*rmucos(ic,0)**2-1.d0)                    
                                                                        
        hess(:)=rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))                              
                                                                        
        hess(ic)=hess(ic)                                               &
     &+fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-distp(0,1)*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic)                                                 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
            enddo 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                                                                        
          case(2100:2199) 
!     p gaussian  -r**(2*(npower+1))*exp(-alpha*r**2) derivative 1100   
                                                                        
         npower=iopt+1-2100 
          
         indorbp=indorb 
                                                                        
         dd(indpar+1)=abs(dd(indpar+1)) 
         dd2=dd(indpar+1) 
         do k=indtmin,indtm 
            distp(k,1)=-r(k)**(2*npower)*dexp(-dd2*r(k)**2) 
         enddo 
                                                                        
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                  z(indorbp,i)=rmu(ic,i)*distp(i,1) 
               enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
                                                                        
            rp1=r(0)**2 
                                                                        
            fun=(npower-dd2*rp1)*distp(0,1)*2.d0/rp1 
                                                                        
            fun3=(npower-dd2*rp1)*fun*2.d0/rp1                          &
     &-2.d0*dd2*distp(0,1)*2.d0/rp1                                     &
     &-2.d0*(npower-dd2*rp1)*distp(0,1)*2.d0/rp1/rp1                    
                                                                        
            indorbp=indorb 
                                                                        
                do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
            z(indorbp,indt+1)=rmu(1,0)*fun*rmucos(1,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+2)=rmu(2,0)*fun*rmucos(2,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
            z(indorbp,indt+3)=rmu(3,0)*fun*rmucos(3,0)          &
     &   *rmu(ic,0)*rmucos(ic,0)                                
                                                                        
       z(indorbp,indt+ic)=z(indorbp,indt+ic)                            &
     &   +distp(0,1)*(2.d0*rmucos(ic,0)**2-1.d0)                    
                                                                        
        hess(:)=rmu(ic,0)*rmucos(ic,0)*                         &
     &(fun3*rmu(:,0)**2*rmucos(:,0)**2                          &
     &+fun*(1.d0-2.d0*rmusin(:,0)**2))                              
                                                                        
        hess(ic)=hess(ic)                                               &
     &+fun*(1.d0-2.d0*rmusin(ic,0)**2)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &+fun*(2.d0*rmucos(ic,0)**2-1.d0)                              &
     &*rmu(ic,0)*rmucos(ic,0)                                   &
     &-distp(0,1)*4.d0*rmusin(ic,0)*rmucos(ic,0)                &
     &*PI/cellscale(ic)                                                 
                                                                        
!                z_xyz(indorbp,:)=hess(:) 
                 z(indorbp,indt+4)=sum(hess(:)) 
               endif 
            enddo 
                !endif for indt                                         
         endif 
                                                                        
         indpar=indpar+1 
         indshell=indshell+3 
         indorb=indorbp 
                                                                        
                                                                        
                                                                        
                                                                        
       case(200) 
! derivative of 199 LA COSTANTE                                         
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=1.d0 
               enddo 
            endif 
                                                                        
         if(typec.ne.1) then 
!       z_xyz(indorbp,:)=0.d0 
        z(indorbp,indt+1:indt+4)=0.d0
                !endif for indt                                         
         endif 
                                                                        
         indshell=indshellp 
         indorb=indorbp 
                                                                        
      case default 
      write(6,*) 'WARNING makefun_pbc: orbital',iopt,'not found' 
                                                                        
      iflagerr=1
                                                                        
      end select 
! ** ** ** ** ** ** ** END OF JASTROW ORBITALS ** ** ** ** ** ** ** ** *
      return 
      END                                           