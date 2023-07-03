!TL off
        subroutine makefun_bump(iopt,iocc,indt,i0,indtmin,indtm,typec,indpar  &
     &,indorb,indshell,nelskip,z,dd,r,rmu,distp,iflagnorm_fake,c) 
     
         use allio, only: cutoff_p ! user-defined cutoff for bump orbitals (see file bump_orbitals_new.f90)
         use constants 
         
         implicit none 
         integer iopt,indt,i,k,nelskip,indpar,iocc(*),indorbp           &
     &,indorb,indshell,indshellp,ic,indtmin,i0              &
     &,iflagnorm_fake,indparp,indtm,npower,typec
         real(8) z(nelskip,i0:*),dd(*),rmu(3,0:indtm),r(0:indtm), &
         distp(0:indtm,20),peff,fun,fun0,fun2,rp1,rp2,rp3,dd1,dd2,          &
         dd3,dd4,dd5,c,dsqrt,funp,fun2p,fun0z,funz,fun2z,                   &
         peff2,arg,c0,c1,cost,cut,cost_fd
     
!     Derivative with respect to Z not implemented
!                                                                       
!        indorb are the number of orbitals occupied before calling      
!        this subroutine                                                
!                                                                       
!        indpar is the number of variational parameters used            
!        before calling this subroutine                                 
!                                                                       
!        indshell is the index of the last  occupied orbital            
!        in the shell, characterized by occupation number iocc(indshell+
!                                                                       
!        z(i,indt+4)  contains the laplacian of the orbital i          
!        z(i,indt+mu) contains the gradient of the orbital i (mu=1,2,3) 
!        In the following given a radial part of the orbital f(r)       
!        fun=1/r  d f(r)/d r                                            
!        fun2= d^2 f(r)/dr^2                                            
                                                                        
                                        
      select case(iopt)
      
      case(16,161) 
! R(r)=exp(-z*r**2) single zeta                                         
                                                                        
                                                                        
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
                                                                        
            dd1=dd(indpar+1) 
                                                                        
!            if(iflagnorm.gt.2) then 
!             c=(2.d0*dd1/pi)**(3.d0/4.d0)*ratios
             c=dd1**0.75d0*ratios
!            endif 
                                                                        
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 
                                                                        
            do i=i0,indtm 
              z(indorbp,i)=distp(i,1) 
            enddo 
                                                                        
            if(typec.ne.1) then 
            if(distp(0,1).ne.0.d0) then
!              the first derivative /r                                  

            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2
                                                                        
                  do i=1,3 
                     z(indorbp,indt+i)=fun*rmu(i,0) 
                  enddo 
                                                                        
                  z(indorbp,indt+4)=2.d0*fun+fun2 
            else
            z(indorbp,indt+1:indt+4)=0.d0
            endif 
                                                                        
            endif 
                                                                        
!           indorb=indorbp 
                                                                        
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshellp 
                                                                        
      case(19) ! derivative of bump gaussian
! R(r)=c*exp(-z*r**2)*(3/4/z-r**2)                                      
          
         indshellp=indshell+1 
                                                                        
         if(iocc(indshellp).eq.1) then 
                                                                        
            indorbp=indorb+1 
                                                                        
            dd1=dd(indpar+1) 
                                                                        
!         if(iflagnorm.gt.2) then 
!         ratios--> ratios*(2/pi)^3/4
!         c=(2.d0*dd1/pi)**(3.d0/4.d0)*ratios
          c=ratios*dd1**0.75d0
!         endif 
                                                                       
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 





                                                                        
            do i=i0,indtm 
            if(distp(i,1).ne.0.d0) then
            rp1=dd1*r(i)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            cost= (rp3+2.d0*rp2)/rp3
            z(indorbp,i)=distp(i,1)*(3.d0/4.d0/dd1-r(i)**2*cost) 
            else
            z(indorbp,i)=0.d0
            endif 
            enddo 
                                                                        
           if(typec.ne.1) then 
!              the first derivative /r                                  
            if(distp(0,1).ne.0.d0) then 
!              the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2



!          fun0z=-distp(0,1)*rp1/dd1*(rp3+2.d0*rp2)/rp3

           funz=2.d0*distp(0,1)*(-81.d0+97d0*rp1-19.d0*rp1**2+rp1**3+&
         &2.d0*rp2*(9.d0-18.d0*rp1+2.d0*rp1**2)+4.d0*rp1*rp2**2)/rp3**2

        fun2z=-2.d0*distp(0,1)*(-729.d0+3798d0*rp1-2594d0*rp1**2+610d0*rp1**3&
       &-59.d0*rp1**4+2.d0*rp1**5+2.d0*rp2*(81.d0-783d0*rp1+600.d0*rp1**2&
       &-112.d0*rp1**3+6.d0*rp1**4)+4.d0*rp1*rp2**2*(45.d0-53.d0*rp1+6.d0*&
       &rp1**2)+16.d0*rp1**2*rp2**3)/rp3**3


!          fun0=0.75d0*fun0+fun0z
           fun=0.75d0*fun/dd1+funz
           fun2=0.75d0*fun2/dd1+fun2z

                                                                        
                  do i=1,3 
                     z(indorbp,indt+i)=fun*rmu(i,0) 
                  enddo 
                                                                        
                  z(indorbp,indt+4)=2.d0*fun+fun2 
            else
            z(indorbp,indt+1:indt+4)=0.d0
            endif 
           endif 
                                                                        
!           indorb=indorbp 
                                                                        
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshellp 
                                                                        
      case(36) 
                                                                        
          
                                                                        
         dd1=dd(indpar+1) 
                                                                        
                                                                        
!        if(iflagnorm.gt.2) then 
!         c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0*ratiop
!        ratiop=ratiop*dsqrt(2.d0)*pi**(-0.75d0)*2^1.25
         c=dd1**1.25d0*ratiop
!        endif 
                                                                        




          do k=indtmin,indtm 
          cost=dd1*r(k)**2
           if(cost.ge.9.d0) then
            distp(k,1)=0.d0
           else
            distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
            if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
           endif
          enddo 



                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
                  z(indorbp,i)=rmu(ic,i)*distp(i,1) 
              enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
          if(distp(0,1).ne.0.d0) then
            fun0=distp(0,1) 
!              the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2




               indorbp=indorb 
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                     do i=1,3 
                        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*   &
     &                       fun                                        
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0 
                     enddo 
             z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2) 
                  endif 
               enddo 
            else
               indorbp=indorb 
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                    z(indorbp,indt+1:indt+4)=0.d0
                 endif 
               enddo 
            endif 
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+3 
!        indorb=indorbp 

       case(44)   ! derivative of 36 
! R(r)=x*exp(-z*r**2)*(5/4/z-r**2)                                      
                                                                        
          
         dd1=dd(indpar+1) 
!        if(iflagnorm.gt.2) then 
!        c=dsqrt(2.d0)*pi**(-0.75d0)*(2.d0*dd1)**1.25d0*ratiop
         c=dd1**1.25d0*ratiop
!        endif 
                                                                        

            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 
                                                                        
         indorbp=indorb 
!                                                                       
         do ic=1,3 
            if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               do i=i0,indtm 
               if(distp(i,1).ne.0.d0) then
               rp1=dd1*r(i)**2
               rp3=-9.d0+rp1
               rp2=dlog(-rp3/9.d0)
               cost= (rp3+2.d0*rp2)/rp3
                   z(indorbp,i)=rmu(ic,i)*distp(i,1)*               &
     &             (1.25d0/dd1-r(i)**2*cost)                   
               else
                   z(indorbp,i)=0.d0
               endif 
               enddo 
            endif 
         enddo 
                                                                        
         if(typec.ne.1) then 
            if(distp(0,1).ne.0.d0) then 
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2


           fun0z=-distp(0,1)*rp1/dd1*(rp3+2.d0*rp2)/rp3

           funz=2.d0*distp(0,1)*(-81.d0+97d0*rp1-19.d0*rp1**2+rp1**3+&
         &2.d0*rp2*(9.d0-18.d0*rp1+2.d0*rp1**2)+4.d0*rp1*rp2**2)/rp3**2

        fun2z=-2.d0*distp(0,1)*(-729.d0+3798d0*rp1-2594d0*rp1**2+610d0*rp1**3&
       &-59.d0*rp1**4+2.d0*rp1**5+2.d0*rp2*(81.d0-783d0*rp1+600.d0*rp1**2&
       &-112.d0*rp1**3+6.d0*rp1**4)+4.d0*rp1*rp2**2*(45.d0-53.d0*rp1+6.d0*&
       &rp1**2)+16.d0*rp1**2*rp2**3)/rp3**3


           fun0=1.25d0*fun0/dd1+fun0z
           fun=1.25d0*fun/dd1+funz
           fun2=1.25d0*fun2/dd1+fun2z


                                                                        
               indorbp=indorb 
                                                                        
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                     do i=1,3 
                        z(indorbp,indt+i)=rmu(ic,0)*rmu(i,0)*   &
     &                       fun                                        
                if(i.eq.ic) z(indorbp,indt+i)=z(indorbp,indt+i)+fun0 
                     enddo 
             z(indorbp,indt+4)=rmu(ic,0)*(4.d0*fun+fun2) 
                  endif 
               enddo 
            else
              indorbp=indorb 
               do ic=1,3 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                    z(indorbp,indt+1:indt+4)=0.d0
                 endif 
               enddo 
            endif 
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+3 
!        indorb=indorbp 
                                                                        
                                                                        
                                                                        
       case(37) 
! d orbitals                                                            
! R(r)= exp(-alpha r^2)                                                 
! each gaussian term is normalized                                      
                                                                        
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!         c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)*ratiod
!         ratiod->ratiod*4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0
          c=dd1**1.75d0*ratiod
!         endif 
                                                                        

            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 


                                                                        
                                                                        
          do i=indtmin,indtm 
                                                           ! lz=0       
      distp(i,2)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d 
                                                          ! lz=+/-2     
      distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d 
                                                   ! lz=+/-2            
      distp(i,4)=rmu(1,i)*rmu(2,i)*cost3d 
                                                   ! lz=+/-1            
      distp(i,5)=rmu(2,i)*rmu(3,i)*cost3d 
                                                   ! lz=+/-1            
      distp(i,6)=rmu(1,i)*rmu(3,i)*cost3d 
          enddo 
                                                                        
                                                                        
         do ic=1,5 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic) 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            if(distp(0,1).ne.0.d0) then
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2


               indorbp=indorb 
               do ic=1,5 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                     do i=1,3 
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)            &
     &                *fun                                              
                       if(ic.eq.1) then 
                          if(i.ne.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                            2.d0*rmu(i,0)*fun0*cost1d         
                          else 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            4.d0*rmu(i,0)*fun0*cost1d         
                          endif 
                       elseif(ic.eq.2) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            2.d0*rmu(i,0)*fun0*cost2d         
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                            2.d0*rmu(i,0)*fun0*cost2d         
                          endif 
                       elseif(ic.eq.3) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(2,0)*fun0*cost3d              
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(1,0)*fun0*cost3d              
                          endif 
                       elseif(ic.eq.4) then 
                          if(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(3,0)*fun0*cost3d              
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(2,0)*fun0*cost3d              
                          endif 
                       elseif(ic.eq.5) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(3,0)*fun0*cost3d              
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(1,0)*fun0*cost3d              
                                !endif for i                            
                          endif 
                                !endif for ic                           
                       endif 
                           !enddo for i                                 
                    enddo 
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
         else

               indorbp=indorb 
               do ic=1,5 
                 if(iocc(indshell+ic).eq.1) then 
                    indorbp=indorbp+1 
                    z(indorbp,indt+1:indt+4)=0.d0
                 endif 
               enddo

         endif 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+5 
!        indorb=indorbp 
                                                                        
                ! derivative of 37 with respect to z                    
       case(45) 
! d orbitals                                                            
! R(r)= c*exp(-z r^2)*(7/4/z-r^2)                                       
                                                                        
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
         c=dd1**1.75d0*ratiod
!        c=4.d0/dsqrt(3.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(7.d0/4.d0)*ratiod
!         endif 
                                                                        
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 
                                                                        
                                                                        
                                                                        
          do i=indtmin,indtm 
                                                           ! lz=0       
      distp(i,2)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d 
                                                          ! lz=+/-2     
      distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d 
                                                   ! lz=+/-2            
      distp(i,4)=rmu(1,i)*rmu(2,i)*cost3d 
                                                   ! lz=+/-1            
      distp(i,5)=rmu(2,i)*rmu(3,i)*cost3d 
                                                   ! lz=+/-1            
      distp(i,6)=rmu(1,i)*rmu(3,i)*cost3d 
          enddo 
                                                                        
                                                                        
         do ic=1,5 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
               if(distp(k,1).ne.0.d0) then
               rp1=dd1*r(k)**2
               rp3=-9.d0+rp1
               rp2=dlog(-rp3/9.d0)
               cost= (rp3+2.d0*rp2)/rp3
         z(indorbp,k)=distp(k,1)*(7.d0/4.d0/dd1-r(k)**2*cost)*&
     &         distp(k,1+ic)                                             
               else
               z(indorbp,k)=0.d0
               endif 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
            if(distp(0,1).ne.0.d0) then
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2

          fun0z=-distp(0,1)*rp1/dd1*(rp3+2.d0*rp2)/rp3

           funz=2.d0*distp(0,1)*(-81.d0+97d0*rp1-19.d0*rp1**2+rp1**3+&
         &2.d0*rp2*(9.d0-18.d0*rp1+2.d0*rp1**2)+4.d0*rp1*rp2**2)/rp3**2

        fun2z=-2.d0*distp(0,1)*(-729.d0+3798d0*rp1-2594d0*rp1**2+610d0*rp1**3&
       &-59.d0*rp1**4+2.d0*rp1**5+2.d0*rp2*(81.d0-783d0*rp1+600.d0*rp1**2&
       &-112.d0*rp1**3+6.d0*rp1**4)+4.d0*rp1*rp2**2*(45.d0-53.d0*rp1+6.d0*&
       &rp1**2)+16.d0*rp1**2*rp2**3)/rp3**3


           fun0=1.75d0*fun0/dd1+fun0z
           fun=1.75d0*fun/dd1+funz
           fun2=1.75d0*fun2/dd1+fun2z
                                                                        
               indorbp=indorb 
               do ic=1,5 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                     do i=1,3 
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)            &
     &                *fun                                              
                       if(ic.eq.1) then 
                          if(i.ne.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                            2.d0*rmu(i,0)*fun0*cost1d         
                          else 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            4.d0*rmu(i,0)*fun0*cost1d         
                          endif 
                       elseif(ic.eq.2) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            2.d0*rmu(i,0)*fun0*cost2d         
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                            2.d0*rmu(i,0)*fun0*cost2d         
                          endif 
                       elseif(ic.eq.3) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(2,0)*fun0*cost3d              
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(1,0)*fun0*cost3d              
                          endif 
                       elseif(ic.eq.4) then 
                          if(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(3,0)*fun0*cost3d              
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(2,0)*fun0*cost3d              
                          endif 
                       elseif(ic.eq.5) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(3,0)*fun0*cost3d              
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                            rmu(1,0)*fun0*cost3d              
                                !endif for i                            
                          endif 
                                !endif for ic                           
                       endif 
                           !enddo for i                                 
                    enddo 
      z(indorbp,indt+4)=distp(0,1+ic)*(6.d0*fun+fun2) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
           else
            indorbp=indorb 
            do ic=1,5 
             if(iocc(indshell+ic).eq.1) then 
               indorbp=indorbp+1 
               z(indorbp,indt+1:indt+4)=0.d0
             endif 
            enddo
           endif 
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+5 
!        indorb=indorbp 
                                                                        

                                                                        
      case(48) 
! f single gaussian orbital                                             
! R(r)= exp(-alpha r^2)                                                 
! normalized                                                            
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
!        c=8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(9.d0/4.d0)*ratiof
!        ratiof-->ratiof*8.d0/dsqrt(15.d0)*(2.d0/pi)**(3.d0/4.d0)
         c=dd1**2.25d0*ratiof
!         endif 
                                                                        
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=cost1f*rmu(3,i)                                    &
     &           *(5.d0*rmu(3,i)**2-3.d0*r(i)**2)               
                                                          ! lz=0        
      distp(i,3)=cost2f*rmu(1,i)                                    &
     &           *(5.d0*rmu(3,i)**2-r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2f*rmu(2,i)                                    &
     &           *(5.d0*rmu(3,i)**2-r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3f*rmu(3,i)                                    &
     &           *(rmu(1,i)**2-rmu(2,i)**2)                     
                                                    ! lz=+/-2           
      distp(i,6)=cost3f*2.d0*rmu(3,i)                               &
     &           *rmu(1,i)*rmu(2,i)                             
                                            ! lz=+/-2                   
      distp(i,7)=cost4f*rmu(1,i)                                    &
     &           *(rmu(1,i)**2-3.d0*rmu(2,i)**2)                
                                                          ! lz=+/-3     
      distp(i,8)=cost4f*rmu(2,i)                                    &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                
                                                          ! lz=+/-3     
          enddo 
                                                                        
                                                                        
         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic) 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            if(distp(0,1).ne.0.d0) then
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2


               indorbp=indorb 
               do ic=1,7 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                     do i=1,3 
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)            &
     &                *fun                                              
                       if(ic.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 6.d0*cost1f*fun0*rmu(i,0)*rmu(3,0)       
                          if(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &       cost1f*fun0*(15.d0*rmu(i,0)**2-3.d0*r(0)**2)       
                          endif 
                       elseif(ic.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 2.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)       
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &             cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)       
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                10.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)       
                          endif 
                       elseif(ic.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 2.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)       
                          if(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &             cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)       
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                10.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)       
                          endif 
                       elseif(ic.eq.4) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)       
                          else 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &              cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)       
                          endif 
                       elseif(ic.eq.5) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)       
                          else 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)       
                          endif 
                       elseif(ic.eq.6) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &         3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)       
                          endif 
                       else 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &         3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)       
                          endif 
                                !endif for ic                           
                       endif 
                           !enddo for i                                 
                    enddo 
      z(indorbp,indt+4)=distp(0,1+ic)*(8.d0*fun+fun2) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 

             else
             indorbp=indorb 
              do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                 indorbp=indorbp+1 
                 z(indorbp,indt+1:indt+4)=0.d0
               endif 
              enddo
             endif 

                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+7 
!        indorb=indorbp 
                                                                        
                                                                        
                ! derivative of 48 with respect to z                    
       case(49) 
! f orbitals                                                            
! R(r)= c*exp(-z r^2)*(9/4/z-r^2)                                       
                                                                        
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
         dd1=dd(indparp) 
                                                                        
!         if(iflagnorm.gt.2) then 
! overall normalization                                                 
         c=dd1**2.25d0*ratiof
!         endif 
                                                                        
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=cost1f*rmu(3,i)                                    &
     &           *(5.d0*rmu(3,i)**2-3.d0*r(i)**2)               
                                                          ! lz=0        
      distp(i,3)=cost2f*rmu(1,i)                                    &
     &           *(5.d0*rmu(3,i)**2-r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2f*rmu(2,i)                                    &
     &           *(5.d0*rmu(3,i)**2-r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3f*rmu(3,i)                                    &
     &           *(rmu(1,i)**2-rmu(2,i)**2)                     
                                                    ! lz=+/-2           
      distp(i,6)=cost3f*2.d0*rmu(3,i)                               &
     &           *rmu(1,i)*rmu(2,i)                             
                                            ! lz=+/-2                   
      distp(i,7)=cost4f*rmu(1,i)                                    &
     &           *(rmu(1,i)**2-3.d0*rmu(2,i)**2)                
                                                          ! lz=+/-3     
      distp(i,8)=cost4f*rmu(2,i)                                    &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                
                                                          ! lz=+/-3     
          enddo 
                                                                        
                                                                        
         do ic=1,7 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
               if(distp(k,1).ne.0.d0) then
               rp1=dd1*r(k)**2
               rp3=-9.d0+rp1
               rp2=dlog(-rp3/9.d0)
               cost= (rp3+2.d0*rp2)/rp3
               z(indorbp,k)=distp(k,1)*(9.d0/4.d0/dd1-r(k)**2*cost)*&
      &        distp(k,1+ic)
               else
               z(indorbp,k)=0.d0
               endif
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            if(distp(0,1).ne.0.d0) then
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2


          fun0z=-distp(0,1)*rp1/dd1*(rp3+2.d0*rp2)/rp3

           funz=2.d0*distp(0,1)*(-81.d0+97d0*rp1-19.d0*rp1**2+rp1**3+&
         &2.d0*rp2*(9.d0-18.d0*rp1+2.d0*rp1**2)+4.d0*rp1*rp2**2)/rp3**2

        fun2z=-2.d0*distp(0,1)*(-729.d0+3798d0*rp1-2594d0*rp1**2+610d0*rp1**3&
       &-59.d0*rp1**4+2.d0*rp1**5+2.d0*rp2*(81.d0-783d0*rp1+600.d0*rp1**2&
       &-112.d0*rp1**3+6.d0*rp1**4)+4.d0*rp1*rp2**2*(45.d0-53.d0*rp1+6.d0*&
       &rp1**2)+16.d0*rp1**2*rp2**3)/rp3**3


           fun0=2.25d0*fun0/dd1+fun0z
           fun=2.25d0*fun/dd1+funz
           fun2=2.25d0*fun2/dd1+fun2z



                                                                        
               indorbp=indorb 
               do ic=1,7 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                     do i=1,3 
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)            &
     &                *fun                                              
                       if(ic.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 6.d0*cost1f*fun0*rmu(i,0)*rmu(3,0)       
                          if(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &       cost1f*fun0*(15.d0*rmu(i,0)**2-3.d0*r(0)**2)       
                          endif 
                       elseif(ic.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 2.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)       
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &             cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)       
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                10.d0*cost2f*fun0*rmu(i,0)*rmu(1,0)       
                          endif 
                       elseif(ic.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 2.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)       
                          if(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &             cost2f*fun0*(5.d0*rmu(3,0)**2-r(0)**2)       
                          elseif(i.eq.3) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                10.d0*cost2f*fun0*rmu(i,0)*rmu(2,0)       
                          endif 
                       elseif(ic.eq.4) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)       
                          else 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &              cost3f*fun0*(rmu(1,0)**2-rmu(2,0)**2)       
                          endif 
                       elseif(ic.eq.5) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(2,0)*rmu(3,0)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(1,0)*rmu(3,0)       
                          else 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 2.d0*cost3f*fun0*rmu(1,0)*rmu(2,0)       
                          endif 
                       elseif(ic.eq.6) then 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &         3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)-       &
     &                 6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)       
                          endif 
                       else 
                          if(i.eq.1) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &                 6.d0*cost4f*fun0*rmu(1,0)*rmu(2,0)       
                          elseif(i.eq.2) then 
                             z(indorbp,indt+i)=z(indorbp,indt+i)+       &
     &         3.d0*cost4f*fun0*(rmu(1,0)**2-rmu(2,0)**2)       
                          endif 
                                !endif for ic                           
                       endif 
                           !enddo for i                                 
                    enddo 
      z(indorbp,indt+4)=distp(0,1+ic)*(8.d0*fun+fun2) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
            else
              indorbp=indorb 
              do ic=1,7 
               if(iocc(indshell+ic).eq.1) then 
                 indorbp=indorbp+1 
                 z(indorbp,indt+1:indt+4)=0.d0
               endif 
              enddo
            endif 
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+7 
!        indorb=indorbp 
                 


      case(51) 
! g single gaussian orbital                                             
! R(r)= exp(-alpha r^2)                                                 
! normalized                                                            
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!      if(iflagnorm.gt.2) then 
! overall normalization                                                 
!      c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)*ratiog
!      ratiog-->ratiog*16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)
       c=dd1**2.75d0*ratiog
!      endif 
                                                                        
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 
                                                                        
                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*r(i)**2+3.d0*r(i)**4)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmu(3,i)                       &
     &           *(7.d0*rmu(3,i)**2-3.d0*r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmu(3,i)                       &
     &           *(7.d0*rmu(3,i)**2-3.d0*r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2-r(i)**2)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmu(2,i)                  &
     &           *(7.d0*rmu(3,i)**2-r(i)**2)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmu(3,i)                       &
     &           *(rmu(1,i)**2-3.0*rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmu(3,i)                       &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmu(2,i)                 &
     &   *(rmu(1,i)**2-rmu(2,i)**2)      
                                                      ! lz=+/-4
          enddo 
                                                                        
                                                                        
         do ic=1,9 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic) 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 


            if(distp(0,1).ne.0.d0) then
                                                                        
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2


                                                                        
                                                                        
               indorbp=indorb 
               do ic=1,9 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                     do i=1,3 
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)*fun
                       if(ic.eq.1) then 
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost1g*fun0*(-60.d0*rmu(1,0)*rmu(3,0)**2+12.d0*rmu(1,0)*r(0)**2)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost1g*fun0*(-60.d0*rmu(2,0)*rmu(3,0)**2+12.d0*rmu(2,0)*r(0)**2)
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost1g*fun0*(80.d0*rmu(3,0)**3-48.d0*rmu(3,0)*r(0)**2)
                       endif
                       elseif(ic.eq.2) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-9.d0*rmu(1,0)**2*rmu(3,0)-3.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))                  
                       else 
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-3.d0*rmu(1,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))   
                       endif
                       elseif(ic.eq.3) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))                  
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-3.d0*rmu(1,0)**2*rmu(3,0)-9.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
                       else 
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-3.d0*rmu(2,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))
                       endif
                       elseif(ic.eq.4) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(-4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(3,0)**2))
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(4.d0*(rmu(2,0)**3-3.d0*rmu(2,0)*rmu(3,0)**2))
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(12.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0))
                       endif
                       elseif(ic.eq.5) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(-2.d0*rmu(2,0)*(3.d0*rmu(1,0)**2+rmu(2,0)**2-6.d0*rmu(3,0)**2))
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(-2.d0*rmu(1,0)*(rmu(1,0)**2+3.d0*rmu(2,0)**2-6.d0*rmu(3,0)**2))
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*24.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)                         
                       endif
                       elseif(ic.eq.6) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      -cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
                       endif
                       elseif(ic.eq.7) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
                       endif                          
                       elseif(ic.eq.8) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(rmu(2,0)**3-3.d0*rmu(1,0)**2*rmu(2,0))
                       endif
                       elseif(ic.eq.9) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
                       endif
                       endif
                           !enddo for i                                 
                    enddo 
      z(indorbp,indt+4)=distp(0,1+ic)*(10.d0*fun+fun2) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
          else
              indorbp=indorb 
              do ic=1,9 
               if(iocc(indshell+ic).eq.1) then 
                 indorbp=indorbp+1 
                 z(indorbp,indt+1:indt+4)=0.d0
               endif 
              enddo
          endif 
                                                                        
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+9
!        indorb=indorbp 
                                                                        
                                                                        
      case(52) 
! g single gaussian orbital                                             
! derivative of 51
! R(r)= exp(-alpha r^2)                                                 
! normalized                                                            
                                                                        
          
         indorbp=indorb 
         indparp=indpar+1 
                                                                        
         dd1=dd(indparp) 
                                                                        
!      if(iflagnorm.gt.2) then 
! overall normalization                                                 
      c=dd1**2.75d0*ratiog
!      c=16.d0/dsqrt(105.d0)*(2.d0/pi)**(3.d0/4.d0)*dd1**(11.d0/4.d0)*ratiog
!      endif 
                                                                        
            do k=indtmin,indtm 
            cost=dd1*r(k)**2
             if(cost.ge.9.d0) then
              distp(k,1)=0.d0
             else
              distp(k,1)=c*dexp(-cost-dlog(1.d0-cost/9.d0)**2) 
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
             endif
            enddo 

                                                                        
                                                                        
          do i=indtmin,indtm 
      distp(i,2)=cost1g*(35.d0*rmu(3,i)**4                          &
     &          -30.d0*rmu(3,i)**2*r(i)**2+3.d0*r(i)**4)
                                                      ! lz=0        
      distp(i,3)=cost2g*rmu(1,i)*rmu(3,i)                       &
     &           *(7.d0*rmu(3,i)**2-3.d0*r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,4)=cost2g*rmu(2,i)*rmu(3,i)                       &
     &           *(7.d0*rmu(3,i)**2-3.d0*r(i)**2)                    
                                                      ! lz=+/-1         
      distp(i,5)=cost3g*(rmu(1,i)**2-rmu(2,i)**2)               &
     &           *(7.d0*rmu(3,i)**2-r(i)**2)                    
                                                      ! lz=+/-2           
      distp(i,6)=cost3g*2.d0*rmu(1,i)*rmu(2,i)                  &
     &           *(7.d0*rmu(3,i)**2-r(i)**2)          
                                                      ! lz=+/-2                   
      distp(i,7)=cost4g*rmu(1,i)*rmu(3,i)                       &
     &           *(rmu(1,i)**2-3.0*rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,8)=cost4g*rmu(2,i)*rmu(3,i)                       &
     &           *(3.d0*rmu(1,i)**2-rmu(2,i)**2)                             
                                                      ! lz=+/-3
      distp(i,9)=cost5g*(rmu(1,i)**4                                &
     &   -6.d0*rmu(1,i)**2*rmu(2,i)**2+rmu(2,i)**4)      
                                                      ! lz=+/-4
      distp(i,10)=cost5g*4.d0*rmu(1,i)*rmu(2,i)                 &
     &   *(rmu(1,i)**2-rmu(2,i)**2)      
                                                      ! lz=+/-4
          enddo 
                                                                        
             
         do ic=1,9 
            if(iocc(indshell+ic).eq.1) then 
              indorbp=indorbp+1 
              do k=i0,indtm 
              if(distp(k,1).ne.0.d0) then
              rp1=dd1*r(k)**2
              rp3=-9.d0+rp1
              rp2=dlog(-rp3/9.d0)
              cost= (rp3+2.d0*rp2)/rp3
              z(indorbp,k)=distp(k,1)*(11.d0/4.d0/dd1-r(k)**2*cost)*&
     &        distp(k,1+ic)                                             
              else
              z(indorbp,k)=0.d0
              endif 
              enddo 
            endif 
         enddo 
                                                                        
                                                                        
         if(typec.ne.1) then 
                                                                        
            if(distp(0,1).ne.0.d0) then
            fun0=distp(0,1)
!           the first derivative /r                                  
            rp1=dd1*r(0)**2
            rp3=-9.d0+rp1
            rp2=dlog(-rp3/9.d0)
            fun=-2.d0*dd1*distp(0,1)*(rp1-9.d0+2.d0*rp2)/rp3 
                                                                        
!              the second derivative                                    
           fun2=2.d0*dd1*distp(0,1)*( -81.d0+176.d0*rp1-37.d0*rp1**2+2.d0*rp1**3&
         &+2.d0*rp2*(9.d0-35.d0*rp1+4.d0*rp1**2)+8.d0*rp1*rp2**2)/rp3**2

          fun0z=-distp(0,1)*rp1/dd1*(rp3+2.d0*rp2)/rp3

           funz=2.d0*distp(0,1)*(-81.d0+97d0*rp1-19.d0*rp1**2+rp1**3+&
         &2.d0*rp2*(9.d0-18.d0*rp1+2.d0*rp1**2)+4.d0*rp1*rp2**2)/rp3**2

        fun2z=-2.d0*distp(0,1)*(-729.d0+3798d0*rp1-2594d0*rp1**2+610d0*rp1**3&
       &-59.d0*rp1**4+2.d0*rp1**5+2.d0*rp2*(81.d0-783d0*rp1+600.d0*rp1**2&
       &-112.d0*rp1**3+6.d0*rp1**4)+4.d0*rp1*rp2**2*(45.d0-53.d0*rp1+6.d0*&
       &rp1**2)+16.d0*rp1**2*rp2**3)/rp3**3


           fun0=2.75d0*fun0/dd1+fun0z
           fun=2.75d0*fun/dd1+funz
           fun2=2.75d0*fun2/dd1+fun2z
                                                                        
               indorbp=indorb 
               do ic=1,9 
                  if(iocc(indshell+ic).eq.1) then 
                     indorbp=indorbp+1 
                     do i=1,3 
                z(indorbp,indt+i)=distp(0,1+ic)*rmu(i,0)*fun
                       if(ic.eq.1) then 
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost1g*fun0*(-60.d0*rmu(1,0)*rmu(3,0)**2+12.d0*rmu(1,0)*r(0)**2)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost1g*fun0*(-60.d0*rmu(2,0)*rmu(3,0)**2+12.d0*rmu(2,0)*r(0)**2)
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost1g*fun0*(80.d0*rmu(3,0)**3-48.d0*rmu(3,0)*r(0)**2)
                       endif
                       elseif(ic.eq.2) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-9.d0*rmu(1,0)**2*rmu(3,0)-3.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))                  
                       else 
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-3.d0*rmu(1,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))   
                       endif
                       elseif(ic.eq.3) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0))                  
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-3.d0*rmu(1,0)**2*rmu(3,0)-9.d0*rmu(2,0)**2*rmu(3,0)+4.d0*rmu(3,0)**3)
                       else 
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost2g*fun0*(-3.d0*rmu(2,0)*(rmu(1,0)**2+rmu(2,0)**2-4.d0*rmu(3,0)**2))
                       endif
                       elseif(ic.eq.4) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(-4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(3,0)**2))
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(4.d0*(rmu(2,0)**3-3.d0*rmu(2,0)*rmu(3,0)**2))
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(12.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0))
                       endif
                       elseif(ic.eq.5) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(-2.d0*rmu(2,0)*(3.d0*rmu(1,0)**2+rmu(2,0)**2-6.d0*rmu(3,0)**2))
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*(-2.d0*rmu(1,0)*(rmu(1,0)**2+3.d0*rmu(2,0)**2-6.d0*rmu(3,0)**2))
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost3g*fun0*24.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)                         
                       endif
                       elseif(ic.eq.6) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      -cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
                       endif
                       elseif(ic.eq.7) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*6.d0*rmu(1,0)*rmu(2,0)*rmu(3,0)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*3.d0*(rmu(1,0)**2-rmu(2,0)**2)*rmu(3,0)
                       else
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost4g*fun0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
                       endif                          
                       elseif(ic.eq.8) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(rmu(2,0)**3-3.d0*rmu(1,0)**2*rmu(2,0))
                       endif
                       elseif(ic.eq.9) then
                       if(i.eq.1) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(3.d0*rmu(1,0)**2*rmu(2,0)-rmu(2,0)**3)
                       elseif(i.eq.2) then
                          z(indorbp,indt+i)=z(indorbp,indt+i)           &
      +cost5g*fun0*4.d0*(rmu(1,0)**3-3.d0*rmu(1,0)*rmu(2,0)**2)
                       endif
                       endif
                           !enddo for i                                 
                    enddo 
      z(indorbp,indt+4)=distp(0,1+ic)*(10.d0*fun+fun2) 
                        !endif for iocc                                 
                 endif 
                       ! enddo fot ic                                   
             enddo 
                                                                        
             else

              indorbp=indorb 
              do ic=1,9 
               if(iocc(indshell+ic).eq.1) then 
                 indorbp=indorbp+1 
                 z(indorbp,indt+1:indt+4)=0.d0
               endif 
              enddo


             endif 
                                                                        
                !endif for indt                                         
         endif 
                                                                        
!        indpar=indpar+1 
!        indshell=indshell+9
!        indorb=indorbp 
                                                                        
                      
       case(199) 
! derivative of 200 LA COSTANTE                                         
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=0.d0 
               enddo 
            endif 
                                                                        
         if(typec.ne.1) z(indorbp,indt+1:indt+4)=0 
                                                                        
!        indshell=indshellp 
!        indorb=indorbp 
       case(200) 
!     THE  COSTANT                                                      
                                                                        
            indorbp=indorb+1 
            indshellp=indshell+1 
                                                                        
            if(iocc(indshellp).eq.1) then 
               do i=i0,indtm 
               z(indorbp,i)=1.d0 
               enddo 
            endif 
                                                                        
         if(typec.ne.1) then 
                  do i=1,3 
                     z(indorbp,indt+i)=0 
                  enddo 
                                                                        
        z(indorbp,indt+4)=0 
                !endif for indt                                         
         endif 
                                                                        
!        indshell=indshellp 
!        indorb=indorbp 
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW SECTION: determinant bump orbitals with user-defined cutoff. 
! For these orbitals the convolution function is a Fermi function.
! These orbitals are useful when computing determinant occupations in order
! to evaluate occupations of localized orbitals such as d-orbitals.
! They are used for the tool: dm_occ.x
! cutoff_p = input parameter which determines the cutoff in (a.u.)^2
! steep_fd = steepness of the Fermi function in (a.u.)^2
! NB: only simple gaussian orbitals s/p/d/f are inmplemented so far.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  case(801)  ! s orbitals, simple gaussian 

      
     indshellp=indshell+1 

     if(iocc(indshellp).eq.1) then 

        indorbp=indorb+1                                                                
        dd1=dd(indpar+1) 
        c=dd1**0.75d0*ratios
        
        do k=indtmin,indtm 
           cost=dd1*r(k)**2
           if( r(k)**2 .gt. cutoff_p + 4.d0*steep_fd) then ! avoid underflow
              distp(k,1) = 0.d0
           else
              cost_fd = 1.d0/(1.d0 + dexp(( r(k)**2-cutoff_p)/steep_fd ) )
              distp(k,1)=c*dexp(-cost)*cost_fd
              if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
           endif
        enddo

        do i=i0,indtm 
           z(indorbp,i)=distp(i,1) 
        enddo

     endif


  case(811)  !  p orbitals, simple gaussian

                                                               
     dd1=dd(indpar+1) 
     c=dd1**1.25d0*ratiop
     
     do k=indtmin,indtm 
        cost=dd1*r(k)**2
        if( r(k)**2 .gt. cutoff_p + 5.d0*steep_fd) then ! avoid underflow
           distp(k,1) = 0.d0
        else
           cost_fd = 1.d0/(1.d0 + dexp(( r(k)**2-cutoff_p)/steep_fd ) )
           distp(k,1)=c*dexp(-cost)*cost_fd
           if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
        endif
     enddo

     indorbp=indorb 

     do ic=1,3 
        if(iocc(indshell+ic).eq.1) then 
           indorbp=indorbp+1 
           do i=i0,indtm 
              z(indorbp,i)=rmu(ic,i)*distp(i,1) 
           enddo
        endif
     enddo
     

  case(821)  ! d orbitals, simple gaussian

      
     indorbp=indorb 
     indparp=indpar+1 
     dd1=dd(indparp) 
     c=dd1**1.75d0*ratiod

     do k=indtmin,indtm 
        cost=dd1*r(k)**2
        if( r(k)**2 .gt. cutoff_p + 5.d0*steep_fd) then ! avoid underflow
           distp(k,1) = 0.d0
        else
           cost_fd = 1.d0/(1.d0 + dexp(( r(k)**2-cutoff_p)/steep_fd ) )
           distp(k,1)=c*dexp(-cost)*cost_fd
           if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
        endif        
     enddo
     
     do i=indtmin,indtm 
        ! lz=0       
        distp(i,2)=(3.d0*rmu(3,i)**2-r(i)**2)*cost1d 
        ! lz=+/-2     
        distp(i,3)=(rmu(1,i)**2-rmu(2,i)**2)*cost2d 
        ! lz=+/-2            
        distp(i,4)=rmu(1,i)*rmu(2,i)*cost3d 
        ! lz=+/-1            
        distp(i,5)=rmu(2,i)*rmu(3,i)*cost3d 
        ! lz=+/-1            
        distp(i,6)=rmu(1,i)*rmu(3,i)*cost3d 
     enddo

     do ic=1,5 
        if(iocc(indshell+ic).eq.1) then 
           indorbp=indorbp+1 
           do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic) 
           enddo
        endif
     enddo

     
  case(831)  ! f orbitals, simple gaussian

      
     indorbp=indorb 
     indparp=indpar+1 

     dd1=dd(indparp)                                                                         
     c=dd1**2.25d0*ratiof

     do k=indtmin,indtm 
        cost=dd1*r(k)**2
        if( r(k)**2 .gt. cutoff_p + 5.d0*steep_fd) then ! avoid underflow
           distp(k,1) = 0.d0
        else
           cost_fd = 1.d0/(1.d0 + dexp(( r(k)**2-cutoff_p)/steep_fd ) )
           distp(k,1)=c*dexp(-cost)*cost_fd
           if(r(k).lt.0.d0) distp(k,1)=-distp(k,1)
        endif
     enddo

     do i=indtmin,indtm 
        distp(i,2)=cost1f*rmu(3,i)*(5.d0*rmu(3,i)**2-3.d0*r(i)**2)               
        ! lz=0        
        distp(i,3)=cost2f*rmu(1,i)*(5.d0*rmu(3,i)**2-r(i)**2)                    
        ! lz=+/-1         
        distp(i,4)=cost2f*rmu(2,i)*(5.d0*rmu(3,i)**2-r(i)**2)                    
        ! lz=+/-1         
        distp(i,5)=cost3f*rmu(3,i)*(rmu(1,i)**2-rmu(2,i)**2)                     
        ! lz=+/-2           
        distp(i,6)=cost3f*2.d0*rmu(3,i)*rmu(1,i)*rmu(2,i)                             
        ! lz=+/-2                   
        distp(i,7)=cost4f*rmu(1,i)*(rmu(1,i)**2-3.d0*rmu(2,i)**2)                
        ! lz=+/-3     
        distp(i,8)=cost4f*rmu(2,i)*(3.d0*rmu(1,i)**2-rmu(2,i)**2)                
        ! lz=+/-3     
     enddo


     do ic=1,7 
        if(iocc(indshell+ic).eq.1) then 
           indorbp=indorbp+1 
           do k=i0,indtm 
              z(indorbp,k)=distp(k,1)*distp(k,1+ic) 
           enddo
        endif
     enddo

      case default 
      write(6,*) 'WARNING makefun: orbital',iopt,'not found'                                                             
      iflagerr=1
                                                                        
      end select 
      
      return 
      
    END                                           
