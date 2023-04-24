!TL off
subroutine ddot_b(n, a, ab, b, bb, ddotb)
    implicit none
    integer n,k
    real*8, intent(inout):: ddotb
    real*8, intent(in):: a(n),b(n)
    real*8, intent(inout)::  ab(n),bb(n)
    do k=1,n
    ab(k) = ab(k) + ddotb * b(k)
    bb(k) = bb(k) + ddotb * a(k)
    enddo
    ddotb=0.d0
    return
end
      
