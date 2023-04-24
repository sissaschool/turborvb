!TL off
subroutine derfc_b(arg1, arg1b, result1b)
    use constants, only : M_2_SQRTPI
    implicit none
    real*8 arg1, arg1b, result1b
    arg1b = arg1b - result1b * M_2_SQRTPI * dexp(-arg1**2)
    result1b = 0_8
    return
end

