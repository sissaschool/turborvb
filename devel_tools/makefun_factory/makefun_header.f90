!TL off
subroutine makefun(iopt,indt,i0,indtmin,indtm,typec,indpar               &
           &,indorb,indshell,nelskip,z,dd,zeta,r,rmu,distp               &
           &,iflagnorm_unused,cr)
    use constants
    implicit none
    integer iopt,indt,i,k,nelskip,indpar,indorbp                         &
           &,indorb,indshell,indshellp,ic,indtmin,i0                     &
           &,iflagnorm_unused,indparp,indtm,npower,typec                 &
           &,ii,jj,kk
    real*8 z(nelskip,i0:*),dd(*),zeta(*),rmu(3,0:indtm)                  &
           &,r(0:indtm)                                                  &
           &,distp(0:indtm,20),peff,fun,fun0,fun2                        &
           &,rp1,rp2,rp3,rp4,rp5,rp6,rp7,rp8                             &
           &,dd1,dd2,dd3,dd4,dd5,c,cr,funp,fun2p,funb                    &
           &,peff2,arg,c0,c1,cost,zv(6),yv(6),xv(6),r2,r4,r6 ! up to i

    integer :: count, multiplicity
    integer, parameter :: max_power = 20
    real*8 :: powers(3,-2:max_power,0:indtm)
  !
  ! indorb are the number of orbitals occupied before calling
  ! this subroutine
  !
  ! indpar is the number of variational parameters used
  ! before calling this subroutine
  !
  ! indshell is the index of the last  occupied orbital
  ! in the shell, characterized by occupation number iocc(indshell)
  !
  ! z(i,indt+4)   contains the laplacian of the orbital i
  ! z(i,indt+mu) contains the gradient of the orbital i (mu=1,2,3)
  ! In the following given a radial part of the orbital f(r)
  ! fun=1/r  d f(r)/d r

  !print *,__FILE__
  !print *,'makefun: iopt=',iopt
  !print *,'makefun: i=',i0,' a=',indtmin,' b=',indtm
  !print *,'makefun: indpar=',indpar,' indorb=',indorb,' indshell=',indshell
  !print *,'makefun: nelskip=',nelskip
